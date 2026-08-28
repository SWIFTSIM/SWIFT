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
