/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2024 Will J. Roper (w.roper@sussex.ac.uk)
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

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "scheduler.h"
#include "space.h"
#include "zoom_cell.h"
#include "zoom_init.h"
#include "zoom_regrid.h"

/**
 * @brief Re-build the top-level cell grid with a zoom region.
 *
 * @param s The #space.
 * @param gravity_properties The properties of gravity, used to calculate
 *                           neighbouring cells.
 * @param nr_nodes The number of MPI ranks.
 * @param verbose Print messages to stdout or not.
 */
void zoom_space_regrid(struct space *s, int verbose) {

  const ticks tic = getticks();
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* If this is our first regrid then we need to get the zoom region
   * geometry before moving on. */
  if (s->cells_top == NULL) {
    zoom_region_init(s, verbose);
  }

  /* Extract the zoom properties. */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* We need to check we have a big enough mesh for the zoom cells.
   * Sadly we have to wait to do this until the engine is attached to the space
   * before we can make this check without breaking the initialisation order,
   * or passing around extra information. Although annoying this just means
   * this error could be triggered after everything is set up instead of
   * during set up. */
  if (s->e != NULL &&
      s->dim[0] / s->e->grav_props->mesh_size > zoom_props->width[0]) {
    error(
        "Mesh too small given the size of top-level zoom cells (width= "
        "%.2f). Should be at least %d cells wide (Currently: %d).",
        zoom_props->width[0], (int)(dim[0] / zoom_props->width[0]),
        s->e->grav_props->mesh_size);
  }

  /* Get the current h_max. */
  double zoom_cell_min = zoom_props->cell_min;
  float h_max =
      get_current_hmax(s, zoom_props->local_zoom_cells_with_particles_top,
                       zoom_props->nr_local_zoom_cells_with_particles,
                       zoom_props->nr_zoom_cells, zoom_cell_min);

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

  /* check that we have at least 1 cell in each dimension */
  if (zoom_cdim[0] == 0 || zoom_cdim[1] == 0 || zoom_cdim[2] == 0) {
    error(
        "Zoom top level cell dimension of size 0 detected (cdim = [%i %i "
        "%i])!\nThis usually indicates a problem with the initial smoothing "
        "lengths of the particles, e.g. a smoothing length that is comparable "
        "in size to the zoom region.",
        zoom_cdim[0], zoom_cdim[1], zoom_cdim[2]);
  }

  /* Do we need to re-build the upper-level cells? */
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
      swift_free("void_cells_top", s->zoom_props->void_cells_top);
      swift_free("neighbour_cells_top", s->zoom_props->neighbour_cells_top);
      swift_free("cells_top", s->cells_top);
      swift_free("multipoles_top", s->multipoles_top);

      /* Setting the new zoom cdim. */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->zoom_props->cdim[ijk] = zoom_cdim[ijk];
      }

      /* Calculate the region geometry. */
      zoom_region_init(s, verbose);
    }

    /* Also free the task arrays, these will be regenerated and we can use the
     * memory while copying the particle arrays. */
    if (s->e != NULL) scheduler_free_tasks(&s->e->sched);

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells =
        (s->cdim[0] * s->cdim[1] * s->cdim[2]) +
        (s->zoom_props->cdim[0] * s->zoom_props->cdim[1] *
         s->zoom_props->cdim[2]) +
        (s->zoom_props->buffer_cdim[0] * s->zoom_props->buffer_cdim[1] *
         s->zoom_props->buffer_cdim[2]);

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
        error(
            "Failed to allocate indices of local top-level background cells.");
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
    zoom_construct_tl_cells(s, ti_current, verbose);

  } /* re-build upper-level cells? */
  else { /* Otherwise, just clean up the cells. */

    /* Free the old cells, if they were allocated. */
    space_free_cells(s);
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
