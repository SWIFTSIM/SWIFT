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
#include "active.h"
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "multipole.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/**
 * @brief Recursively split a cell.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param buff A buffer for particle sorting, should be of size at least
 *        c->hydro.count or @c NULL.
 * @param sbuff A buffer for particle sorting, should be of size at least
 *        c->stars.count or @c NULL.
 * @param bbuff A buffer for particle sorting, should be of size at least
 *        c->black_holes.count or @c NULL.
 * @param gbuff A buffer for particle sorting, should be of size at least
 *        c->grav.count or @c NULL.
 * @param sink_buff A buffer for particle sorting, should be of size at least
 *        c->sinks.count or @c NULL.
 */
void space_split_recursive(struct space *s, struct cell *c,
                           struct cell_buff *restrict buff,
                           struct cell_buff *restrict sbuff,
                           struct cell_buff *restrict bbuff,
                           struct cell_buff *restrict gbuff,
                           struct cell_buff *restrict sink_buff) {

  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  const int with_self_gravity = s->with_self_gravity;
  const int depth = c->depth;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct sink *sinks = c->sinks.parts;

  /* If the buff is NULL, allocate it, and remember to free it. */
  const int allocate_buffer = (buff == NULL && gbuff == NULL && sbuff == NULL &&
                               bbuff == NULL && sink_buff == NULL);
  if (allocate_buffer) {
    if (count > 0) {
      if (swift_memalign("tempbuff", (void **)&buff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * count) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (parts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
        buff[k].x[0] = parts[k].x[0];
        buff[k].x[1] = parts[k].x[1];
        buff[k].x[2] = parts[k].x[2];
      }
    }
    if (gcount > 0) {
      if (swift_memalign("tempgbuff", (void **)&gbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * gcount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
        if (gparts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (gparts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
        gbuff[k].x[0] = gparts[k].x[0];
        gbuff[k].x[1] = gparts[k].x[1];
        gbuff[k].x[2] = gparts[k].x[2];
      }
    }
    if (scount > 0) {
      if (swift_memalign("tempsbuff", (void **)&sbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * scount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
        if (sparts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (sparts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
        sbuff[k].x[0] = sparts[k].x[0];
        sbuff[k].x[1] = sparts[k].x[1];
        sbuff[k].x[2] = sparts[k].x[2];
      }
    }
    if (bcount > 0) {
      if (swift_memalign("tempbbuff", (void **)&bbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * bcount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < bcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
        if (bparts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (bparts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
        bbuff[k].x[0] = bparts[k].x[0];
        bbuff[k].x[1] = bparts[k].x[1];
        bbuff[k].x[2] = bparts[k].x[2];
      }
    }
    if (sink_count > 0) {
      if (swift_memalign("temp_sink_buff", (void **)&sink_buff,
                         SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * sink_count) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < sink_count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
        if (sinks[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (sinks[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
        sink_buff[k].x[0] = sinks[k].x[0];
        sink_buff[k].x[1] = sinks[k].x[1];
        sink_buff[k].x[2] = sinks[k].x[2];
      }
    }
  }

  /* If the depth is too large, we have a problem and should stop. */
  if (depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting cells, aborting. This is "
        "most likely due to having too many particles at the exact same "
        "position, making the construction of a tree impossible.",
        space_cell_maxdepth);
  }

  /* Split or let it be? */
  if ((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize))) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->hydro.count = 0;
      cp->grav.count = 0;
      cp->stars.count = 0;
      cp->sinks.count = 0;
      cp->black_holes.count = 0;
      cp->hydro.count_total = 0;
      cp->grav.count_total = 0;
      cp->sinks.count_total = 0;
      cp->stars.count_total = 0;
      cp->black_holes.count_total = 0;
      cp->hydro.ti_old_part = c->hydro.ti_old_part;
      cp->grav.ti_old_part = c->grav.ti_old_part;
      cp->grav.ti_old_multipole = c->grav.ti_old_multipole;
      cp->stars.ti_old_part = c->stars.ti_old_part;
      cp->sinks.ti_old_part = c->sinks.ti_old_part;
      cp->black_holes.ti_old_part = c->black_holes.ti_old_part;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->width[0] = c->width[0] / 2;
      cp->width[1] = c->width[1] / 2;
      cp->width[2] = c->width[2] / 2;
      cp->dmin = c->dmin / 2;
      if (k & 4) cp->loc[0] += cp->width[0];
      if (k & 2) cp->loc[1] += cp->width[1];
      if (k & 1) cp->loc[2] += cp->width[2];
      cp->depth = c->depth + 1;
      cp->split = 0;
      cp->hydro.h_max = 0.f;
      cp->hydro.h_max_active = 0.f;
      cp->hydro.dx_max_part = 0.f;
      cp->hydro.dx_max_sort = 0.f;
      cp->stars.h_max = 0.f;
      cp->stars.h_max_active = 0.f;
      cp->stars.dx_max_part = 0.f;
      cp->stars.dx_max_sort = 0.f;
      cp->sinks.r_cut_max = 0.f;
      cp->sinks.r_cut_max_active = 0.f;
      cp->sinks.dx_max_part = 0.f;
      cp->black_holes.h_max = 0.f;
      cp->black_holes.h_max_active = 0.f;
      cp->black_holes.dx_max_part = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->top = c->top;
      cp->super = NULL;
      cp->hydro.super = NULL;
      cp->grav.super = NULL;
      cp->flags = 0;
      star_formation_logger_init(&cp->stars.sfh);
#ifdef WITH_MPI
      cp->mpi.tag = -1;
#endif  // WITH_MPI
#ifdef WITH_ZOOM_REGION
      cp->tl_cell_type = c->tl_cell_type;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
      cell_assign_cell_index(cp, c);
#endif
    }

    /* Split the cell's particle data. */
    cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts,
               c->black_holes.parts - s->bparts, c->sinks.parts - s->sinks,
               buff, sbuff, bbuff, gbuff, sink_buff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff, *progeny_bbuff = bbuff,
                     *progeny_sink_buff = sink_buff;

    for (int k = 0; k < 8; k++) {

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->hydro.count == 0 && cp->grav.count == 0 && cp->stars.count == 0 &&
          cp->black_holes.count == 0 && cp->sinks.count == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {

        /* Recurse */
        space_split_recursive(s, cp, progeny_buff, progeny_sbuff, progeny_bbuff,
                              progeny_gbuff, progeny_sink_buff);

        /* Update the pointers in the buffers */
        progeny_buff += cp->hydro.count;
        progeny_gbuff += cp->grav.count;
        progeny_sbuff += cp->stars.count;
        progeny_bbuff += cp->black_holes.count;
        progeny_sink_buff += cp->sinks.count;
      }
    }
  }   /* Split or let it be? */

  /* Otherwise, clean up progeny. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
  }

  /* Clean up. */
  if (allocate_buffer) {
    if (buff != NULL) swift_free("tempbuff", buff);
    if (gbuff != NULL) swift_free("tempgbuff", gbuff);
    if (sbuff != NULL) swift_free("tempsbuff", sbuff);
    if (bbuff != NULL) swift_free("tempbbuff", bbuff);
    if (sink_buff != NULL) swift_free("temp_sink_buff", sink_buff);
  }
}

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

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    space_split_recursive(s, c, NULL, NULL, NULL, NULL, NULL);
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
