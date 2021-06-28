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
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "cell.h"
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
