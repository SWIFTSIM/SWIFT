/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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

/* This object's header. */
#include "engine.h"

/**
 * @brief Mapper function to init activate dm parts
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_init_dmpart_mapper(void *map_data, int num_elements,
                                     void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Initialize all the DM particles */
      cell_init_dmpart(c, e, /* force the init=*/1);
    }
  }
}

/**
 * @brief Initialize DM particles
 * forward to the current time.
 *
 * @param e The #engine.
 * @param drift_mpoles Do we want to drift all the multipoles as well?
 */
void engine_init_dmparts(struct engine *e) {

  const ticks tic = getticks();

  if (!e->restarting) {

    /* Normal case: We have a list of local cells with tasks to play with */

    if (e->s->nr_dmparts > 0) {
      threadpool_map(&e->threadpool, engine_do_init_dmpart_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }

  } else {

    /* When restarting, the list of local cells with tasks does not yet
       exist. We use the raw list of top-level cells instead */

    if (e->s->nr_dmparts > 0) {
      threadpool_map(&e->threadpool, engine_do_init_dmpart_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

