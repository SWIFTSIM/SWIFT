/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timeline.h"
#include "timers.h"

#ifdef SWIFT_DEBUG_CHECKS
/* Temporary probe: a time_bin <= 0 in the limiter exchange reaches
 * part_is_active() on the receiver and asserts (slurm-66505880). Report
 * what was packed, what arrived, and which particle trips the check. */
#define LIMITER_TIMEBIN_PROBE
#endif

#ifdef LIMITER_TIMEBIN_PROBE
static integertime_t probe_task_ti_run(const struct task *t) {
  return (t == NULL) ? -1 : t->ti_run;
}
#endif

/**
 * @brief Recompute the hydro h_max_active of a (foreign) cell and its
 * progeny after their particles' time-bins have just been overwritten.
 *
 * cell_unpack_timebin() can lower a particle's time-bin (waking it up, e.g.
 * following a time-step limiter/sync propagation) without going through
 * recv_xv/recv_rho, which are the only tasks that normally refresh
 * h_max_active (see runner_do_recv_part()). Left stale, h_max_active can
 * under-count a particle that just became active, since it was excluded
 * from the max the last time h_max_active was computed. This trips the "h
 * larger than h_max_active" debug check in the limiter pair/self
 * interactions.
 *
 * @param c The #cell (and its progeny) whose time-bins were just unpacked.
 * @param e The #engine containing the current time information.
 */
static void cell_update_hydro_h_max_active(struct cell *c,
                                           const struct engine *e) {

  float h_max_active = 0.f;

  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_update_hydro_h_max_active(c->progeny[k], e);
        h_max_active = max(h_max_active, c->progeny[k]->hydro.h_max_active);
      }
    }
  } else {
    const struct part *parts = c->hydro.parts;
    for (int i = 0; i < c->hydro.count; ++i) {
      if (part_is_inhibited(&parts[i], e)) continue;
#ifdef LIMITER_TIMEBIN_PROBE
      if (parts[i].time_bin <= 0)
        message(
            "LIMITER_TIMEBIN_PROBE recv-part: cellID=%lld owner=%d depth=%d "
            "count=%d i=%d id=%lld time_bin=%d wakeup=%d "
            "to_be_synchronized=%d min_ngb_time_bin=%d ti_current=%lld "
            "max_active_bin=%d",
            c->cellID, c->nodeID, c->depth, c->hydro.count, i, parts[i].id,
            parts[i].time_bin, parts[i].limiter_data.wakeup,
            parts[i].limiter_data.to_be_synchronized,
            parts[i].limiter_data.min_ngb_time_bin, e->ti_current,
            e->max_active_bin);
#endif
      if (part_is_active(&parts[i], e))
        h_max_active = max(h_max_active, parts[i].h);
    }
  }

  c->hydro.h_max_active = h_max_active;
}

/**
 * @brief Pack the data needed by the time-step limiter loop prior to sending
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param buffer The array to allocate and fill.
 * @param timer Are we timing this ?
 */
void runner_do_pack_limiter(struct runner *r, struct cell *c, void **buffer,
                            const int timer) {

  const size_t count = c->hydro.count * sizeof(timebin_t);
  if (posix_memalign((void **)buffer, SWIFT_CACHE_ALIGNMENT, count) != 0)
    error("Error allocating timebin send buffer");

  cell_pack_timebin(c, (timebin_t *)*buffer);

#ifdef LIMITER_TIMEBIN_PROBE
  const struct engine *e = r->e;
  const timebin_t *packed = (const timebin_t *)*buffer;
  for (int i = 0; i < c->hydro.count; ++i) {
    if (packed[i] == time_bin_inhibited) continue;
    if (packed[i] <= 0) {
      const struct part *p = &c->hydro.parts[i];
      message(
          "LIMITER_TIMEBIN_PROBE pack: cellID=%lld depth=%d count=%d i=%d "
          "id=%lld packed_bin=%d live_bin=%d wakeup=%d to_be_synchronized=%d "
          "ti_current=%lld timestep_ti_run=%lld limiter_ti_run=%lld "
          "sync_ti_run=%lld",
          c->cellID, c->depth, c->hydro.count, i, p->id, packed[i], p->time_bin,
          p->limiter_data.wakeup, p->limiter_data.to_be_synchronized,
          e->ti_current, probe_task_ti_run(c->super->timestep),
          probe_task_ti_run(c->super->timestep_limiter),
          probe_task_ti_run(c->super->timestep_sync));
    }
  }
#endif
}

/**
 * @brief UnPack the data needed by the time-step limiter loop after receiving
 * it
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param buffer The array to read from and free.
 * @param timer Are we timing this ?
 */
void runner_do_unpack_limiter(struct runner *r, struct cell *c, void *buffer,
                              const int timer) {

  cell_unpack_timebin(c, (timebin_t *)buffer);

#ifdef LIMITER_TIMEBIN_PROBE
  {
    const timebin_t *received = (const timebin_t *)buffer;
    int n_bad = 0, first_bad = -1, last_bad = -1;
    for (int i = 0; i < c->hydro.count; ++i) {
      if (received[i] <= 0) {
        if (n_bad == 0) first_bad = i;
        last_bad = i;
        ++n_bad;
      }
    }
    if (n_bad > 0)
      message(
          "LIMITER_TIMEBIN_PROBE recv-cell: cellID=%lld owner=%d depth=%d "
          "split=%d count=%d n_bad=%d first_bad=%d last_bad=%d "
          "ti_current=%lld",
          c->cellID, c->nodeID, c->depth, c->split, c->hydro.count, n_bad,
          first_bad, last_bad, r->e->ti_current);
  }
#endif

  /* The time-bins we just overwrote may have woken particles up; make sure
     h_max_active reflects that before the limiter task reads it. */
  cell_update_hydro_h_max_active(c, r->e);

  free(buffer);
}

/**
 * @brief Pack the data needed by the gravity loop prior to sending
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param buffer The array to allocate and fill.
 * @param timer Are we timing this ?
 */
void runner_do_pack_gpart(struct runner *r, struct cell *c, void **buffer,
                          const int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Catch an undrifted gpart on the sender, before it becomes a mystery
   * stale slot on some other rank. c->grav.super/depth pin down whether
   * the packed range even matches what the drift task covered. */
  for (int i = 0; i < c->grav.count; ++i) {
    const struct gpart *gp = &c->grav.parts[i];
    if (gp->ti_drift != r->e->ti_current && !gpart_is_inhibited(gp, r->e))
      error(
          "Packing an undrifted gpart: i=%d count=%d c->cellID=%lld "
          "c->depth=%d c->split=%d c->grav.ti_old_part=%lld "
          "c->grav.ti_old_part_on_entry=%lld c->grav.count_drifted=%d "
          "c->grav.drift_force_on_entry=%d c->grav.super->cellID=%lld "
          "c->grav.super->depth=%d c->grav.super->grav.ti_old_part=%lld "
          "gp.ti_drift=%lld e->ti_current=%lld gp.type=%d",
          i, c->grav.count, c->cellID, c->depth, c->split, c->grav.ti_old_part,
          c->grav.ti_old_part_on_entry, c->grav.count_drifted,
          c->grav.drift_force_on_entry, c->grav.super->cellID,
          c->grav.super->depth, c->grav.super->grav.ti_old_part, gp->ti_drift,
          r->e->ti_current, gp->type);
  }
#endif

  const size_t count = c->grav.count * sizeof(struct gpart_foreign);
  if (posix_memalign((void **)buffer, SWIFT_CACHE_ALIGNMENT, count) != 0)
    error("Error allocating gpart send buffer");

  cell_pack_gpart(c, *buffer);
}

/**
 * @brief Pack the data needed by the fof loop prior to sending
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param buffer The array to allocate and fill.
 * @param timer Are we timing this ?
 */
void runner_do_pack_fof(struct runner *r, struct cell *c, void **buffer,
                        const int timer) {

  const size_t count = c->grav.count * sizeof(struct gpart_fof_foreign);
  if (posix_memalign((void **)buffer, SWIFT_CACHE_ALIGNMENT, count) != 0)
    error("Error allocating gpart send buffer");

  cell_pack_fof_gpart(c, *buffer);
}
