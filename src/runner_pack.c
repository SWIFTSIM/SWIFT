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
#include "timers.h"

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
 * @param max_active_bin The current maximal active time-bin.
 */
static void cell_update_hydro_h_max_active(struct cell *c,
                                           const timebin_t max_active_bin) {

  float h_max_active = 0.f;

  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_update_hydro_h_max_active(c->progeny[k], max_active_bin);
        h_max_active = max(h_max_active, c->progeny[k]->hydro.h_max_active);
      }
    }
  } else {
    const struct part *parts = c->hydro.parts;
    for (int i = 0; i < c->hydro.count; ++i) {
      if (parts[i].time_bin == time_bin_inhibited) continue;
      if (parts[i].time_bin <= max_active_bin)
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

  /* The time-bins we just overwrote may have woken particles up; make sure
     h_max_active reflects that before the limiter task reads it. */
  cell_update_hydro_h_max_active(c, r->e->max_active_bin);

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

  /* Snapshot count as count_packed for c and every descendant, before
   * star/sink formation (which runs after this pack) can change count.
   * Propagated to the foreign copy via task_subtype_grav_counts. */
  cell_snapshot_grav_count_packed(c);

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
