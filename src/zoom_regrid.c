/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "gravity_properties.h"
#include "scheduler.h"
#include "zoom_region.h"

#ifdef WITH_ZOOM_REGION
/**
 * @brief Re-build the top-level cell grid.
 *
 * @param s The #space.
 * @param verbose Print messages to stdout or not.
 */
void space_regrid_zoom(struct space *s,
                       struct gravity_props *gravity_properties, int verbose) {

  const size_t nr_parts = s->nr_parts;
  const size_t nr_sparts = s->nr_sparts;
  const size_t nr_bparts = s->nr_bparts;
  const size_t nr_sinks = s->nr_sinks;
  const ticks tic = getticks();
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the cells and get the current h_max, when using a zoom region
   * h_max needs to be set by the zoom cells. */
  // tic = getticks();
  const double zoom_cell_min = s->zoom_props->cell_min;
  float h_max = zoom_cell_min / kernel_gamma / space_stretch;
  if (nr_parts > 0) {

    /* Can we use the list of local non-empty top-level cells? */
    if (s->local_cells_with_particles_top != NULL) {
      for (int k = 0; k < s->nr_local_cells_with_particles; ++k) {
        const struct cell *c =
            &s->cells_top[s->local_cells_with_particles_top[k]];

        /* We only want to consider zoom cells to avoid over inflating the
         * smoothing length */
        if (c->type != zoom) continue;

        if (c->hydro.h_max > h_max) {
          h_max = c->hydro.h_max;
        }
        if (c->stars.h_max > h_max) {
          h_max = c->stars.h_max;
        }
        if (c->black_holes.h_max > h_max) {
          h_max = c->black_holes.h_max;
        }
        if (c->sinks.r_cut_max > h_max) {
          h_max = c->sinks.r_cut_max / kernel_gamma;
        }
      }

      /* Can we instead use all the top-level cells? */
    } else if (s->cells_top != NULL) {
      /* We only want to consider zoom cells to avoid over inflating the
       * smoothing length */
      for (int k = s->zoom_props->bkg_cell_offset; k < s->nr_cells; k++) {
        const struct cell *c = &s->cells_top[k];

        if (c->nodeID == engine_rank && c->hydro.h_max > h_max) {
          h_max = c->hydro.h_max;
        }
        if (c->nodeID == engine_rank && c->stars.h_max > h_max) {
          h_max = c->stars.h_max;
        }
        if (c->nodeID == engine_rank && c->black_holes.h_max > h_max) {
          h_max = c->black_holes.h_max;
        }
        if (c->nodeID == engine_rank && c->sinks.r_cut_max > h_max) {
          h_max = c->sinks.r_cut_max / kernel_gamma;
        }
      }

      /* Last option: run through the particles */
    } else {
      for (size_t k = 0; k < nr_parts; k++) {
        if (s->parts[k].h > h_max) h_max = s->parts[k].h;
      }
      for (size_t k = 0; k < nr_sparts; k++) {
        if (s->sparts[k].h > h_max) h_max = s->sparts[k].h;
      }
      for (size_t k = 0; k < nr_bparts; k++) {
        if (s->bparts[k].h > h_max) h_max = s->bparts[k].h;
      }
      for (size_t k = 0; k < nr_sinks; k++) {
        if (s->sinks[k].r_cut > h_max) h_max = s->sinks[k].r_cut / kernel_gamma;
      }
    }
  }

/* If we are running in parallel, make sure everybody agrees on
   how large the largest cell should be. */
#ifdef WITH_MPI
  {
    float buff;
    if (MPI_Allreduce(&h_max, &buff, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    h_max = buff;
  }
#endif
  if (verbose)
    message("h_max is %.3e (zoom_cell_min=%.3e).", h_max, zoom_cell_min);

  /* Get the new putative zoom cell dimensions. */
  const double dmax =
      max3(s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  const int zoom_cdim[3] = {
      (int)floor((dmax + 0.01 * zoom_cell_min) /
                 fmax(h_max * kernel_gamma * space_stretch, zoom_cell_min)),
      (int)floor((dmax + 0.01 * zoom_cell_min) /
                 fmax(h_max * kernel_gamma * space_stretch, zoom_cell_min)),
      (int)floor((dmax + 0.01 * zoom_cell_min) /
                 fmax(h_max * kernel_gamma * space_stretch, zoom_cell_min))};

/* In MPI-Land, changing the top-level cell size requires that the
 * global partition is recomputed and the particles redistributed.
 * Be prepared to do that. */
#ifdef WITH_MPI
  double oldwidth[3] = {0., 0., 0.};
  double oldcdim[3] = {0., 0., 0.};
  double oldzoomwidth[3] = {0., 0., 0.};
  double oldzoomcdim[3] = {0., 0., 0.};
  int *oldnodeIDs = NULL;
  if (zoom_cdim[0] < s->zoom_props->cdim[0] ||
      zoom_cdim[1] < s->zoom_props->cdim[1] ||
      zoom_cdim[2] < s->zoom_props->cdim[2]) {

    error("Regridding is currently not supported, decrease h_max!!");

    /* Capture state of current space. */
    oldcdim[0] = s->cdim[0];
    oldcdim[1] = s->cdim[1];
    oldcdim[2] = s->cdim[2];
    oldwidth[0] = s->width[0];
    oldwidth[1] = s->width[1];
    oldwidth[2] = s->width[2];
    oldzoomcdim[0] = s->zoom_props->cdim[0];
    oldzoomcdim[1] = s->zoom_props->cdim[1];
    oldzoomcdim[2] = s->zoom_props->cdim[2];
    oldzoomwidth[0] = s->zoom_props->width[0];
    oldzoomwidth[1] = s->zoom_props->width[1];
    oldzoomwidth[2] = s->zoom_props->width[2];

    if ((oldnodeIDs =
             (int *)swift_malloc("nodeIDs", sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate temporary nodeIDs.");

    int cid = 0;
    /* First loop over zoom cells */
    for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
      for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
        for (int k = 0; k < s->zoom_props->cdim[2]; k++) {
          cid = cell_getid(oldzoomcdim, i, j, k);
          oldnodeIDs[cid] = s->cells_top[cid].nodeID;
        }
      }
    }
    /* Now lets do the natural cells */
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {
          cid = cell_getid(oldcdim, i, j, k) + s->zoom_props->bkg_cell_offset;
          oldnodeIDs[cid] = s->cells_top[cid].nodeID;
        }
      }
    }
  }

  /* Are we about to allocate new top level cells without a regrid?
   * Can happen when restarting the application. */
  const int no_regrid = (s->cells_top == NULL && oldnodeIDs == NULL);
#endif /* WITH_MPI */

  /* Do we need to re-build the upper-level cells? */
  // tic = getticks();
  if (s->cells_top == NULL || zoom_cdim[0] < s->zoom_props->cdim[0] ||
      zoom_cdim[1] < s->zoom_props->cdim[1] ||
      zoom_cdim[2] < s->zoom_props->cdim[2]) {

    /* Free the old cells, if they were allocated. */
    if (s->cells_top != NULL) {
      space_free_cells(s);
      swift_free("local_cells_with_tasks_top", s->local_cells_with_tasks_top);
      swift_free("local_cells_top", s->local_cells_top);
      swift_free("local_zoom_cells_top", s->zoom_props->local_zoom_cells_top);
      swift_free("local_bkg_cells_top", s->zoom_props->local_bkg_cells_top);
      swift_free("local_zoom_cells_with_particles_top",
                 s->zoom_props->local_zoom_cells_with_particles_top);
      swift_free("local_bkg_cell_with_particless_top",
                 s->zoom_props->local_bkg_cells_with_particles_top);
      swift_free("local_buffer_cell_with_particless_top",
                 s->zoom_props->local_buffer_cells_with_particles_top);
      swift_free("cells_with_particles_top", s->cells_with_particles_top);
      swift_free("local_cells_with_particles_top",
                 s->local_cells_with_particles_top);
      swift_free("cells_top", s->cells_top);
      swift_free("multipoles_top", s->multipoles_top);
    }

    /* Also free the task arrays, these will be regenerated and we can use the
     * memory while copying the particle arrays. */
    if (s->e != NULL) scheduler_free_tasks(&s->e->sched);

    /* Setting the new zoom cdim. */
    for (int ijk = 0; ijk < 3; ijk++) {
      s->zoom_props->cdim[ijk] = zoom_cdim[ijk];
    }

    message("Constructing zoom region.");
    /* Compute the bounds of the zoom region from the DM particles
     * and define the background cells based on this. */
    construct_zoom_region(s, verbose);

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells =
      (s->cdim[0] * s->cdim[1] * s->cdim[2]) +
      (s->zoom_props->cdim[0] * s->zoom_props->cdim[1] *
       s->zoom_props->cdim[2]) +
      (s->zoom_props->buffer_cdim[0] * s->zoom_props->buffer_cdim[1] *
       s->zoom_props->buffer_cdim[2]);

    message("Total number of cells = %d", s->nr_cells);

    if (swift_memalign("cells_top", (void **)&s->cells_top, cell_align,
                       s->nr_cells * sizeof(struct cell)) != 0)
      error("Failed to allocate top-level cells.");
    bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

    /* Allocate the multipoles for the top-level cells. */
    if (s->with_self_gravity) {
      if (swift_memalign("multipoles_top", (void **)&s->multipoles_top,
                         multipole_align,
                         s->nr_cells * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate top-level multipoles.");
      bzero(s->multipoles_top, s->nr_cells * sizeof(struct gravity_tensors));
    }

    /* Allocate the indices of local cells */
    if (swift_memalign("local_cells_top", (void **)&s->local_cells_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells.");
    bzero(s->local_cells_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local cells with tasks */
    if (swift_memalign("local_cells_with_tasks_top",
                       (void **)&s->local_cells_with_tasks_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells with tasks.");
    bzero(s->local_cells_with_tasks_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of cells with particles */
    if (swift_memalign("cells_with_particles_top",
                       (void **)&s->cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of top-level cells with particles.");
    bzero(s->cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local cells with particles */
    if (swift_memalign("local_cells_with_particles_top",
                       (void **)&s->local_cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level cells with "
          "particles.");
    bzero(s->local_cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local zoom cells */
    if (swift_memalign("local_zoom_cells_top",
                       (void **)&s->zoom_props->local_zoom_cells_top,
                       SWIFT_STRUCT_ALIGNMENT,
                       s->zoom_props->nr_zoom_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level zoom cells.");
    bzero(s->zoom_props->local_zoom_cells_top,
          s->zoom_props->nr_zoom_cells * sizeof(int));

    /* Allocate the indices of local bkg cells */
    if (swift_memalign("local_bkg_cells_top",
                       (void **)&s->zoom_props->local_bkg_cells_top,
                       SWIFT_STRUCT_ALIGNMENT,
                       s->zoom_props->nr_bkg_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level background cells.");
    bzero(s->zoom_props->local_bkg_cells_top,
          s->zoom_props->nr_bkg_cells * sizeof(int));

    /* Allocate the indices of local buffer cells */
    if (s->zoom_props->with_buffer_cells) {
      if (swift_memalign("local_buffer_cells_top",
                         (void **)&s->zoom_props->local_buffer_cells_top,
                         SWIFT_STRUCT_ALIGNMENT,
                         s->zoom_props->nr_buffer_cells * sizeof(int)) != 0)
        error("Failed to allocate indices of local top-level background cells.");
      bzero(s->zoom_props->local_buffer_cells_top,
            s->zoom_props->nr_buffer_cells * sizeof(int));
    }

    /* Allocate the indices of local zoom cells with particles */
    if (swift_memalign(
            "local_zoom_cells_with_particles_top",
            (void **)&s->zoom_props->local_zoom_cells_with_particles_top,
            SWIFT_STRUCT_ALIGNMENT,
            s->zoom_props->nr_zoom_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level zoom cells with "
          "particles.");
    bzero(s->zoom_props->local_zoom_cells_with_particles_top,
          s->zoom_props->nr_zoom_cells * sizeof(int));

    /* Allocate the indices of local bkg cells with particles */
    if (swift_memalign(
            "local_bkg_cells_with_particles_top",
            (void **)&s->zoom_props->local_bkg_cells_with_particles_top,
            SWIFT_STRUCT_ALIGNMENT,
            s->zoom_props->nr_bkg_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level background cells with "
          "particles.");
    bzero(s->zoom_props->local_bkg_cells_with_particles_top,
          s->zoom_props->nr_bkg_cells * sizeof(int));

    /* Allocate the indices of local buffer cells with particles */
    if (s->zoom_props->with_buffer_cells) {
      if (swift_memalign(
                         "local_buffer_cells_with_particles_top",
                         (void **)&s->zoom_props->local_buffer_cells_with_particles_top,
                         SWIFT_STRUCT_ALIGNMENT,
                         s->zoom_props->nr_buffer_cells * sizeof(int)) != 0)
        error(
              "Failed to allocate indices of local top-level buffer cells with "
              "particles.");
      bzero(s->zoom_props->local_buffer_cells_with_particles_top,
            s->zoom_props->nr_buffer_cells * sizeof(int));
    }

    /* Set the cells' locks */
    for (int k = 0; k < s->nr_cells; k++) {
      if (lock_init(&s->cells_top[k].hydro.lock) != 0)
        error("Failed to init spinlock for hydro.");
      if (lock_init(&s->cells_top[k].grav.plock) != 0)
        error("Failed to init spinlock for gravity.");
      if (lock_init(&s->cells_top[k].grav.mlock) != 0)
        error("Failed to init spinlock for multipoles.");
      if (lock_init(&s->cells_top[k].grav.star_formation_lock) != 0)
        error("Failed to init spinlock for star formation (gpart).");
      if (lock_init(&s->cells_top[k].stars.lock) != 0)
        error("Failed to init spinlock for stars.");
      if (lock_init(&s->cells_top[k].sinks.lock) != 0)
        error("Failed to init spinlock for sinks.");
      if (lock_init(&s->cells_top[k].sinks.sink_formation_lock) != 0)
        error("Failed to init spinlock for sink formation.");
      if (lock_init(&s->cells_top[k].black_holes.lock) != 0)
        error("Failed to init spinlock for black holes.");
      if (lock_init(&s->cells_top[k].stars.star_formation_lock) != 0)
        error("Failed to init spinlock for star formation (spart).");
    }

    /* Construct both grids of cells */
    const double dmin = min3(s->width[0], s->width[1], s->width[2]);
    construct_tl_cells_with_zoom_region(s, s->cdim, dmin, ti_current,
                                        gravity_properties, verbose);

#ifdef WITH_MPI
    if (oldnodeIDs != NULL) {
      /* We have changed the top-level cell dimension, so need to redistribute
       * cells around the nodes. We repartition using the old space node
       * positions as a grid to resample. */
      if (s->e->nodeID == 0)
        message(
            "basic cell dimensions have increased - recalculating the "
            "global partition.");

      if (!partition_space_to_space_zoom(oldwidth, oldcdim, oldzoomwidth,
                                         oldzoomcdim, oldnodeIDs, s)) {

        /* Failed, try another technique that requires no settings. */
        message("Failed to get a new partition, trying less optimal method");
        struct partition initial_partition;
#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
        initial_partition.type = INITPART_METIS_NOWEIGHT;
#else
        initial_partition.type = INITPART_VECTORIZE;
#endif
        partition_initial_partition(&initial_partition, s->e->nodeID,
                                    s->e->nr_nodes, s);
      }

      /* Re-distribute the particles to their new nodes. */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);

      /* Finished with these. */
      swift_free("nodeIDs", oldnodeIDs);

    } else if (no_regrid && s->e != NULL) {
      /* If we have created the top-levels cells and not done an initial
       * partition (can happen when restarting), then the top-level cells
       * are not assigned to a node, we must do that and then associate the
       * particles with the cells. Note requires that
       * partition_store_celllist() was called once before, or just before
       * dumping the restart files.*/
      partition_restore_celllist(s, s->e->reparttype);

      /* Now re-distribute the particles, should just add to cells? */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);
    }
#endif /* WITH_MPI */

    // message( "rebuilding upper-level cells took %.3f %s." ,
    // clocks_from_ticks(double)(getticks() - tic), clocks_getunit());

  }      /* re-build upper-level cells? */
  else { /* Otherwise, just clean up the cells. */

    /* Free the old cells, if they were allocated. */
    space_free_cells(s);
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
#endif /* WITH_ZOOM_REGION */
