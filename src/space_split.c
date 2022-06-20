/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/**
 * @brief #threadpool mapper function to split cells if they contain
 *        too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void space_sort_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;
  const int nbits = 21;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];

    /* Sort this cell. */
    cell_split_sort(s, c,
                    c->hydro.parts - s->parts,
                    c->stars.parts - s->sparts,
                    c->black_holes.parts - s->bparts,
                    c->sinks.parts - s->sinks, nbits);

  }
}

#ifdef WITH_ZOOM_REGION

/**
 * @brief A wrapper for #threadpool mapper function to split background cells if
 * they contain too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void bkg_space_sort_mapper(void *map_data, int num_cells,
                                 void *extra_data) {
  space_sort_mapper(map_data, num_cells, extra_data);
}

/**
 * @brief A wrapper for #threadpool mapper function to split zoom cells if they
 * contain too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void zoom_space_sort_mapper(void *map_data, int num_cells,
                                  void *extra_data) {
  space_sort_mapper(map_data, num_cells, extra_data);
}

#endif


/**
 * @brief #threadpool mapper function to split cells if they contain
 *        too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void space_split_sort_mapper(void *map_data, int num_cells, void *extra_data,
                             int tid) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;
  const int nbits = 21;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];

    /* Finally, we can split this cell */
    cell_split_recursive(s, c, nbits - 1, c->loc, c->width, tid);

  }
}

#ifdef WITH_ZOOM_REGION

/**
 * @brief A wrapper for #threadpool mapper function to split background cells if
 * they contain too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void bkg_space_split_sort_mapper(void *map_data, int num_cells,
                                 void *extra_data, int tid) {
  space_split_sort_mapper(map_data, num_cells, extra_data, tid);
}

/**
 * @brief A wrapper for #threadpool mapper function to split zoom cells if they
 * contain too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void zoom_space_split_sort_mapper(void *map_data, int num_cells,
                                  void *extra_data, int tid) {
  space_split_sort_mapper(map_data, num_cells, extra_data, tid);
}

#endif

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * This is done in parallel using threads in the #threadpool.
 * Only do this for the local non-empty top-level cells.
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, int verbose) {

  const ticks tic = getticks();

  /* /\* Allocate the first level of progeny to avoid down time waiting on locks */
  /*  * (unused ones will be freed later) *\/ */
  /* for (int ind = 0; ind < s->nr_local_cells_with_particles; ind++) { */
  /*   struct cell *c = &s->cells_top[s->local_cells_with_particles_top[ind]]; */

  /*   /\* Only preallocate progeny for cells that will require splitting. *\/ */
  /*   if ((s->with_self_gravity && c->grav.count > space_splitsize) || */
  /*       (!s->with_self_gravity && */
  /*        (c->hydro.count > space_splitsize || */
  /*         c->stars.count > space_splitsize))){ */
      
  /*     space_getcells(s, 8, c->progeny); */
      
  /*   } */
  /* } */

  /* Sort particles ready for the cell splitting. */
#ifdef WITH_ZOOM_REGION
  if (s->with_zoom_region) {
    threadpool_map(&s->e->threadpool, bkg_space_sort_mapper,
                   s->zoom_props->local_bkg_cells_with_particles_top,
                   s->zoom_props->nr_local_bkg_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
    threadpool_map(&s->e->threadpool, zoom_space_sort_mapper,
                   s->zoom_props->local_zoom_cells_with_particles_top,
                   s->zoom_props->nr_local_zoom_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    threadpool_map(&s->e->threadpool, space_sort_mapper,
                   s->local_cells_with_particles_top,
                   s->nr_local_cells_with_particles, sizeof(int),
                   threadpool_auto_chunk_size, s);
  }
#else
  threadpool_map(&s->e->threadpool, space_sort_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, s);
#endif

    /* Sort particles ready for the cell splitting. */
#ifdef WITH_ZOOM_REGION
  if (s->with_zoom_region) {
    threadpool_map_with_tid(&s->e->threadpool, bkg_space_split_sort_mapper,
                   s->zoom_props->local_bkg_cells_with_particles_top,
                   s->zoom_props->nr_local_bkg_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
    threadpool_map_with_tid(&s->e->threadpool, zoom_space_split_sort_mapper,
                   s->zoom_props->local_zoom_cells_with_particles_top,
                   s->zoom_props->nr_local_zoom_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    threadpool_map_with_tid(&s->e->threadpool, space_split_sort_mapper,
                   s->local_cells_with_particles_top,
                   s->nr_local_cells_with_particles, sizeof(int),
                   threadpool_auto_chunk_size, s);
  }
#else
  threadpool_map_with_tid(&s->e->threadpool, space_split_sort_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, s);
#endif
  
  /* Calculate the properties of the cell heirarchy */
#ifdef WITH_ZOOM_REGION
  if (s->with_zoom_region) {
    threadpool_map(&s->e->threadpool, bkg_cell_props_mapper,
                   s->zoom_props->local_bkg_cells_with_particles_top,
                   s->zoom_props->nr_local_bkg_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
    threadpool_map(&s->e->threadpool, zoom_cell_props_mapper,
                   s->zoom_props->local_zoom_cells_with_particles_top,
                   s->zoom_props->nr_local_zoom_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    threadpool_map(&s->e->threadpool, cell_props_mapper,
                   s->local_cells_with_particles_top,
                   s->nr_local_cells_with_particles, sizeof(int),
                   threadpool_auto_chunk_size, s);
  }
#else
  threadpool_map(&s->e->threadpool, cell_props_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, s);
#endif

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
