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
void space_split_sort_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;
  const int nbits = 21;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
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
void bkg_space_split_sort_mapper(void *map_data, int num_cells,
                                 void *extra_data) {
  space_split_sort_mapper(map_data, num_cells, extra_data);
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
                                  void *extra_data) {
  space_split_sort_mapper(map_data, num_cells, extra_data);
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
void space_split_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;
  const int nbits = 21;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];

    /* Initialise the thread local copy of the cell */
    struct cell *temp_c = NULL;
    
    /* Allocare the tempoary cell so we can
       avoid memory movement overhead in sorts */
    if (swift_memalign("temp_cell", (void**)&temp_c,
                       cell_align,
                       sizeof(struct cell)) != 0)
      error("Error while allocating temporary memory for cell");

    /* Copy cell contents into temporary local cell */
    memcpy(temp_c, c, sizeof(struct cell));

    /* Split this cell */
    cell_split_recursive(s, c, nbits - 1, c->loc, c->width);
    
    /* Replace the cell with the local cell */
    cells_top[local_cells_with_particles[ind]] = *temp_c;

    /* Free up now unused cell,
     NOTE: causes an error so currently this is a memory leak. */
    //free(c);

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
void bkg_space_split_mapper(void *map_data, int num_cells, void *extra_data) {
  space_split_mapper(map_data, num_cells, extra_data);
}

/**
 * @brief A wrapper for #threadpool mapper function to split zoom cells if they
 * contain too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void zoom_space_split_mapper(void *map_data, int num_cells, void *extra_data) {
  space_split_mapper(map_data, num_cells, extra_data);
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

    /* Sort particles ready for the cell splitting. */
#ifdef WITH_ZOOM_REGION
  if (s->with_zoom_region) {
    threadpool_map(&s->e->threadpool, bkg_space_split_sort_mapper,
                   s->zoom_props->local_bkg_cells_with_particles_top,
                   s->zoom_props->nr_local_bkg_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
    threadpool_map(&s->e->threadpool, zoom_space_split_sort_mapper,
                   s->zoom_props->local_zoom_cells_with_particles_top,
                   s->zoom_props->nr_local_zoom_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    threadpool_map(&s->e->threadpool, space_split_sort_mapper,
                   s->local_cells_with_particles_top,
                   s->nr_local_cells_with_particles, sizeof(int),
                   threadpool_auto_chunk_size, s);
  }
#else
  threadpool_map(&s->e->threadpool, space_split_sort_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, s);
#endif

  /* Split the cell heirarchy. */
#ifdef WITH_ZOOM_REGION
  if (s->with_zoom_region) {
    threadpool_map(&s->e->threadpool, bkg_space_split_mapper,
                   s->zoom_props->local_bkg_cells_with_particles_top,
                   s->zoom_props->nr_local_bkg_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
    threadpool_map(&s->e->threadpool, zoom_space_split_mapper,
                   s->zoom_props->local_zoom_cells_with_particles_top,
                   s->zoom_props->nr_local_zoom_cells_with_particles,
                   sizeof(int), threadpool_auto_chunk_size, s);
  } else {
    threadpool_map(&s->e->threadpool, space_split_mapper,
                   s->local_cells_with_particles_top,
                   s->nr_local_cells_with_particles, sizeof(int),
                   threadpool_auto_chunk_size, s);
  }
#else
  threadpool_map(&s->e->threadpool, space_split_mapper,
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
