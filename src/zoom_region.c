/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
 *               2022 Will Roper (w.roper@sussex.ac.uk)
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
#include "zoom_region.h"

#include "../config.h"
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "multipole.h"
#include "proxy.h"
#include "space.h"

#include <float.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief 
 *
 * @param params Swift parameter structure.
 * @param s The space
 */
int get_cell_grids_with_buffer_cells(struct space *s,
                                     struct gravity_props *gravity_properties,
                                     double *max_dim, double ini_dim,
                                     int verbose) {

  /* The number of background cells needs to be odd. */
  for (int ijk = 0; ijk < 3; ijk++) {
    /* if (s->cdim[ijk] % 2 == 0) */
    /*   s->cdim[ijk] -= 1; */
    s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
    s->iwidth[ijk] = 1.0 / s->width[ijk];
  }

  /* Calculate how many background cells we need in the buffer region. The
   * goal is to have this as large as could be necessary, overshooting
   * isn't an issue. */
  const double max_distance =
    gravity_properties->r_s * gravity_properties->r_cut_max_ratio;

  /* Find the buffer region boundaries. The zoom region is already centred on
   * the middle of the box. */
  for (int ijk = 0; ijk < 3; ijk++) {

    /* Find the background cell containing lower and upper bounds of the zoom
     * regions "gravity reach". */
    int lower =
      (s->zoom_props->region_bounds[(ijk * 2)] - max_distance) * s->iwidth[ijk];
    int upper =
      (s->zoom_props->region_bounds[(ijk * 2) + 1] + max_distance) * s->iwidth[ijk];
    
    s->zoom_props->buffer_bounds[(ijk * 2)] = lower * s->width[ijk];
    s->zoom_props->buffer_bounds[(ijk * 2) + 1] = (upper + 1) * s->width[ijk];
  }

  /* Define the extent of the buffer region. */
  double buffer_dim =
    s->zoom_props->buffer_bounds[1] - s->zoom_props->buffer_bounds[0];

  /* Calculate the number of zoom regions covered by the buffer region. */
  int nr_zoom_regions = (int)(buffer_dim / ini_dim);

  /* If the region to buffer ratio is odd we need to have an odd number of
   * zoom regions. */
  if ((s->zoom_props->region_buffer_ratio % 2 == 1) &&
      (nr_zoom_regions % 2 == 0))
    nr_zoom_regions -= 1;

  /* Calculate the new zoom region dimension. */
  *max_dim = buffer_dim / nr_zoom_regions;

  /* Set the buffer cells properties. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->buffer_cdim[ijk] = nr_zoom_regions * s->zoom_props->region_buffer_ratio;
    s->zoom_props->buffer_width[ijk] = *max_dim / s->zoom_props->region_buffer_ratio;
    s->zoom_props->buffer_iwidth[ijk] =
      1.0 / s->zoom_props->buffer_width[ijk];
  }

  return ((*max_dim) < ini_dim);

}

/**
 * @brief Read parameter file for "ZoomRegion" properties, and initialize the
 * zoom_region_properties struct.
 *
 * @param params Swift parameter structure.
 * @param s The space
 */
void zoom_region_init(struct swift_params *params, struct space *s,
                      struct gravity_props *gravity_properties,
                      int verbose) {
#ifdef WITH_ZOOM_REGION
  /* Are we running with a zoom region? */
  s->with_zoom_region =
      parser_get_opt_param_int(params, "ZoomRegion:enable", 0);

  /* If so... */
  if (s->with_zoom_region) {

    /* Zoom region properties are stored in a structure. */
    s->zoom_props = (struct zoom_region_properties *)malloc(
        sizeof(struct zoom_region_properties));
    if (s->zoom_props == NULL)
      error("Error allocating memory for the zoom parameters.");

    /* Set the zoom cdim. */
    s->zoom_props->cdim[0] =
        parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                                 space_max_top_level_cells_default);
    s->zoom_props->cdim[1] =
        parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                                 space_max_top_level_cells_default);
    s->zoom_props->cdim[2] =
        parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                                 space_max_top_level_cells_default);

    /* Set the target background cdim, default is a negative value so that if no
     * value is given for a target then the zoom region defines the background
     * cell size. */
    s->zoom_props->bkg_cdim[0] =
        parser_get_opt_param_int(params,
                                 "ZoomRegion:bkg_top_level_cells",
                                 space_max_top_level_cells_default);
    s->zoom_props->bkg_cdim[1] =
        parser_get_opt_param_int(params,
                                 "ZoomRegion:bkg_top_level_cells",
                                 space_max_top_level_cells_default);
    s->zoom_props->bkg_cdim[2] =
        parser_get_opt_param_int(params,
                                 "ZoomRegion:bkg_top_level_cells",
                                 space_max_top_level_cells_default);
    
    /* Get the ratio between the zoom region size and buffer cell size.
     * Ignored if buffer cells aren't needed. */
    s->zoom_props->region_buffer_ratio =
        parser_get_opt_param_int(params,
                                 "ZoomRegion:region_dim_buffer_cell_ratio",
                                 1);

    /* Ensure we have been given a power of 2 times the region buffer ratio
     * for cdim. */
    if (!(((s->zoom_props->cdim[0] / s->zoom_props->region_buffer_ratio) &
          ((s->zoom_props->cdim[0] /s->zoom_props->region_buffer_ratio) - 1))
          == 0))
      error("Scheduler:max_top_level_cells must be a power "
            "of 2 times region_dim_buffer_cell_ratio (by default 1)"
            "when running with a zoom region!");

    /* Extract the zoom width boost factor (used to define the buffer around the
     * zoom region). */
    s->zoom_props->zoom_boost_factor =
        parser_get_opt_param_float(params,
                                   "ZoomRegion:buffer_region_ratio",
                                   1.5);
    message("Zoom region boost factor from params is %.2f", s->zoom_props->zoom_boost_factor);
    
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
    /* If we are doing a metis decomp are we using wedge in the background? */
    s->zoom_props->use_bkg_wedges =
        parser_get_opt_param_int(params,
                                 "DomainDecomposition:background_wedge_decomp",
                                 0);

    /* If we are doing a metis decomp are we doing each grid separately? */
    s->zoom_props->separate_decomps =
        parser_get_opt_param_int(params,
                                 "DomainDecomposition:separate_decomps",
                                 0);
#endif

    /* Define the background grid. NOTE: Can be updated later.*/
    for (int ijk = 0; ijk < 3; ijk++) {
      s->cdim[ijk] = s->zoom_props->bkg_cdim[ijk];
      s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
      s->iwidth[ijk] = 1.0 / s->width[ijk];
    }

    /* Initialise the number of wanders (unused if with_hydro == False)*/
    s->zoom_props->nr_wanderers = 0;

    /* Get an initial dimension for the zoom region and its geometric mid point.
     */
    double new_zoom_boundary[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
    double midpoint[3] = {0.0, 0.0, 0.0};
    const size_t nr_gparts = s->nr_gparts;
    double mtot = 0.0;
    double com[3] = {0.0, 0.0, 0.0};

    /* Get the shift from the ICs since this hasn't been applied yet. */
    double shift[3] = {0.0, 0.0, 0.0};
    parser_get_opt_param_double_array(params, "InitialConditions:shift", 3,
                                      shift);

    /* Find the min/max location in each dimension for each mask gravity
     * particle, and their COM. */
    for (size_t k = 0; k < nr_gparts; k++) {
      if (s->gparts[k].type != swift_type_dark_matter &&
          s->gparts[k].type != swift_type_gas) continue;

      /* Shift initial positions by IC shift. */
      const double x = s->gparts[k].x[0] + shift[0];
      const double y = s->gparts[k].x[1] + shift[1];
      const double z = s->gparts[k].x[2] + shift[2];

      /* Wrap if periodic. */
      if (s->periodic) {
        box_wrap(x, 0.0, s->dim[0]);
        box_wrap(y, 0.0, s->dim[1]);
        box_wrap(z, 0.0, s->dim[2]);
      }

      /* Ammend boundaries for this particle. */
      if (x < new_zoom_boundary[0]) new_zoom_boundary[0] = x;
      if (x > new_zoom_boundary[1]) new_zoom_boundary[1] = x;
      if (y < new_zoom_boundary[2]) new_zoom_boundary[2] = y;
      if (y > new_zoom_boundary[3]) new_zoom_boundary[3] = y;
      if (z < new_zoom_boundary[4]) new_zoom_boundary[4] = z;
      if (z > new_zoom_boundary[5]) new_zoom_boundary[5] = z;

      /* Total up mass and position for COM. */
      mtot += s->gparts[k].mass;
      com[0] += x * s->gparts[k].mass;
      com[1] += y * s->gparts[k].mass;
      com[2] += z * s->gparts[k].mass;
    }

#ifdef WITH_MPI
    /* Share answers amoungst nodes. */

    /* Boundary. */
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[0], 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[1], 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[2], 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[3], 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[4], 1, MPI_DOUBLE, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[5], 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    /* CoM. */
    MPI_Allreduce(MPI_IN_PLACE, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    /* Finalize CoM calcuation. */
    const double imass = 1.0 / mtot;
    com[0] *= imass;
    com[1] *= imass;
    com[2] *= imass;

    /* Store result. */
    for (int ijk = 0; ijk < 3; ijk++) s->zoom_props->com[ijk] = com[ijk];

    if (verbose) message("com: [%f %f %f]", com[0], com[1], com[2]);

    if (verbose)
      message(
          "initial_dim: [%f %f %f] initial_zoom_boundary: [%f-%f %f-%f %f-%f]",
          new_zoom_boundary[1] - new_zoom_boundary[0],
          new_zoom_boundary[3] - new_zoom_boundary[2],
          new_zoom_boundary[5] - new_zoom_boundary[4], new_zoom_boundary[0],
          new_zoom_boundary[1], new_zoom_boundary[2], new_zoom_boundary[3],
          new_zoom_boundary[4], new_zoom_boundary[5]);

    /* Get the initial dimensions and midpoint. */
    double ini_dims[3] = {0.0, 0.0, 0.0};
    for (int ijk = 0; ijk < 3; ijk++) {
      ini_dims[ijk] =
          (new_zoom_boundary[(ijk * 2) + 1] - new_zoom_boundary[ijk * 2]);
      midpoint[ijk] = new_zoom_boundary[(ijk * 2) + 1] - (ini_dims[ijk] / 2);
    }

    /* Throw an error if the zoom region extends over the box boundries.
     * TODO: This could be fixed automatically! */
    double shiftx = 0.;
    double shifty = 0.;
    double shiftz = 0.;
    if ((ini_dims[0] > s->dim[0] / 2) || (ini_dims[1] > s->dim[1] / 2) ||
        (ini_dims[2] > s->dim[2] / 2)) {
      if (ini_dims[0] > s->dim[0] / 2) shiftx = s->dim[0] / 2;
      if (ini_dims[1] > s->dim[1] / 2) shifty = s->dim[1] / 2;
      if (ini_dims[2] > s->dim[2] / 2) shiftz = s->dim[2] / 2;
      error(
          "Zoom region extends beyond the boundaries of the box. "
          "Shift the ICs by [%f, %f, %f]",
          shiftx, shifty, shiftz);
    }

    /* Calculate the shift needed to place the mid point of the high res
     * particles at the centre of the box. This shift is applied to the
     * particles in space_init in space.c */
    const double box_mid[3] = {s->dim[0] / 2, s->dim[1] / 2, s->dim[2] / 2};
    for (int ijk = 0; ijk < 3; ijk++) {
      s->zoom_props->zoom_shift[ijk] = box_mid[ijk] - midpoint[ijk];
    }
    if (verbose) {
      message("box_mid = [%f %f %f] midpoint = [%f %f %f]", box_mid[0],
              box_mid[1], box_mid[2], midpoint[0], midpoint[1], midpoint[2]);
      message("Need to shift the box by [%e, %e, %e] to centre the zoom region",
              s->zoom_props->zoom_shift[0], s->zoom_props->zoom_shift[1],
              s->zoom_props->zoom_shift[2]);
    }

    /* Let's shift the COM.
     * NOTE: boundaries are recalculated relative to box centre later. */
    for (int ijk = 0; ijk < 3; ijk++)
      s->zoom_props->com[ijk] += s->zoom_props->zoom_shift[ijk];

    /* Compute maximum side length of the zoom region, we need zoom dim to be
     * equal. */
    double ini_dim = max3(ini_dims[0], ini_dims[1], ini_dims[2]) *
                     s->zoom_props->zoom_boost_factor;

    /* If the zoom region is much smaller than a background cell we need to
     * construct the buffer cell region to limit the number of background
     * cells. */
    double max_dim = ini_dim;
    if (max_dim < s->width[0] / 2) {

      /* Set the initial zoom_region boundaries with boost factor.
       * The zoom region is already centred on the middle of the box */
      for (int ijk = 0; ijk < 3; ijk++) {
        /* Set the new boundaries. */
        s->zoom_props->region_bounds[(ijk * 2)] =
          (s->dim[ijk] / 2) - (ini_dim / 2);
        s->zoom_props->region_bounds[(ijk * 2) + 1] =
          (s->dim[ijk] / 2) + (ini_dim / 2);
      }

      /* Flag that we have buffer cells. */
      s->zoom_props->with_buffer_cells = 1;

      /* And initialise the count of empty background cells that house them. */
      s->zoom_props->nr_empty_cells = 0;

      /* Compute the cell grid properties. */
      if (get_cell_grids_with_buffer_cells(s, gravity_properties,
                                           &max_dim, ini_dim,
                                           verbose))
        error("Found a zoom region smaller than the high resolution particle "
              "distribution! Adjust the cell structure "
              "(ZoomRegion:bkg_top_level_cells, Scheduler:max_top_level_cells"
              " and ZoomRegion:region_dim_buffer_cell_ratio)");
      
    }

    /* If it is smaller but not drastically smaller we can simply tessalate
     * cells the size of the zoom region across the whole volume without
     * drastically effecting the size of the zoom region. */
    else if (max_dim < s->width[0]) {

      if (verbose)
        message("WARNING: bkg_top_level_cells (%d) resulted in a large "
                "increase in zoom region size. Falling back on "
                "tessalating zoom region across volume.", s->cdim[0]);

      /* Ensure an odd integer number of the zoom regions tessalate the box. */
      int nr_zoom_regions = (int)(s->dim[0] / max_dim);
      if (nr_zoom_regions % 2 == 0) nr_zoom_regions -= 1;
      max_dim = s->dim[0] / nr_zoom_regions;

      /* Redefine the background cells using this new width */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->cdim[ijk] = nr_zoom_regions;
        s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
        s->iwidth[ijk] = 1.0 / s->width[ijk];
      }

      /* Declare we have no buffer region. */
      s->zoom_props->with_buffer_cells = 0;

      /* Zero the buffer region. */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->zoom_props->buffer_bounds[(ijk * 2)] = 0;
        s->zoom_props->buffer_bounds[(ijk * 2) + 1] = 0;
        s->zoom_props->buffer_cdim[ijk] = 0;
        s->zoom_props->buffer_width[ijk] = 0;
        s->zoom_props->buffer_iwidth[ijk] = 0;
      } 
      
    }

    /* Otherwise, the zoom region is larger than a background cell and we need
     * to modify the cell structure to allow having multiple void cells
     * containing the zoom reigon and buffer cells the same size as the
     * background cells.
     * NOTE: with this the number of background cells is defined by geometry
     * and attempts to get as close as possible to the user defined cdim from
     * the parameter file. */
    else {

      /* First we need to define the zoom region width. */
      int nr_zoom_regions = (int)(s->dim[0] / max_dim);
      max_dim = s->dim[0] / nr_zoom_regions;

      /* Define the requested background cdim as the target. */
      int target_bkg_cdim = s->cdim[0];

      /* Now we can define the background grid. */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->cdim[ijk] =
          (int)floor((s->dim[ijk] + 0.1 * max_dim) / max_dim);
      }

      /* Compute the new number of a background cells. */
      int new_bkg_cdim = s->cdim[0];
      while (new_bkg_cdim < target_bkg_cdim) {
        new_bkg_cdim *= 2;
      }

      if (verbose)
        message("Modifying background cdim from %d to %d", s->cdim[0],
                new_bkg_cdim);

      /* Set the background cdim. */
      s->cdim[0] = new_bkg_cdim;
      s->cdim[1] = new_bkg_cdim;
      s->cdim[2] = new_bkg_cdim;

      /* Set the background cell width. */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
        s->iwidth[ijk] = 1.0 / s->width[ijk];
      }

      /* Declare we have no buffer region. */
      s->zoom_props->with_buffer_cells = 0;

      /* Zero the buffer region. */
      for (int ijk = 0; ijk < 3; ijk++) {
        s->zoom_props->buffer_bounds[(ijk * 2)] = 0;
        s->zoom_props->buffer_bounds[(ijk * 2) + 1] = 0;
        s->zoom_props->buffer_cdim[ijk] = 0;
        s->zoom_props->buffer_width[ijk] = 0;
        s->zoom_props->buffer_iwidth[ijk] = 0;
      } 
    }
    
    /* Find the new boundaries with this extra width and boost factor.
     * The zoom region is already centred on the middle of the box */
    for (int ijk = 0; ijk < 3; ijk++) {
      /* Set the new boundaries. */
      s->zoom_props->region_bounds[(ijk * 2)] =
          (s->dim[ijk] / 2) - (max_dim / 2);
      s->zoom_props->region_bounds[(ijk * 2) + 1] =
          (s->dim[ijk] / 2) + (max_dim / 2);

      /* Set the reigon dim. */
      s->zoom_props->dim[ijk] = max_dim;
    }
    
    if (verbose) {
      message("Initial buffer_region_size = [%.2f %.2f %.2f]",
              (max_dim - ini_dims[0]) / 2, (max_dim - ini_dims[1]) / 2,
              (max_dim - ini_dims[2]) / 2);
      message("Calculated buffer_region_ratio = %.2f",
              max_dim / max3(ini_dims[0], ini_dims[1], ini_dims[2]));
    }

    /* Let's be safe and error if we have drastically changed the size of the
     * buffer region. */
    if ((max_dim / ini_dim) >= 2)
      error("WARNING: The buffer region has to be 2x larger than requested."
            "Either increase ZoomRegion:buffer_region_ratio or increase the "
            "number of background cells.");
    
    /* Set the minimum allowed zoom cell width. */
    const double zoom_dmax = max3(s->zoom_props->dim[0], s->zoom_props->dim[1],
                                  s->zoom_props->dim[2]);
    s->zoom_props->cell_min = 0.99 * zoom_dmax / s->zoom_props->cdim[0];
  }
#endif
}

/**
 * @brief For a given particle location, what TL cell does it belong in nested
 *        grid?
 *
 * This is a simple helper function to reduce repition.
 *
 * @param cdim The cell grid dimensions.
 * @param bounds The edges of this nested region.
 * @param x, y, z Location of particle.
 * @param iwidth The width of a cell in this grid.
 * @param offset The offset of this cell type in cells_top.
 */
int cell_getid_with_bounds(const int *cdim, const double *bounds, 
                           const double x, const double y, const double z,
                           const double *iwidth, const int offset) {

  /* Get the cell ijk coordinates in this grid. */
  const int i = (x - bounds[0]) * iwidth[0];
  const int j = (y - bounds[2]) * iwidth[1];
  const int k = (z - bounds[4]) * iwidth[2];

  /* Which zoom TL cell are we in? */
  const int cell_id = cell_getid(cdim, i, j, k) + offset;

  return cell_id;
}


/**
 * @brief For a given particle location, what TL cell does it belong to?
 *
 * Slightly more complicated in the zoom case, as there are now two embedded TL
 * grids.
 *
 * First see if the particle is in the background grid, if it is a void cell
 * check the nested cell types.
 *
 * @param s The space.
 * @param x, y, z Location of particle.
 */
int cell_getid_zoom(const struct space *s, const double x, const double y,
                    const double z) {
#ifdef WITH_ZOOM_REGION
  int cell_id;

  /* Lets get some properties of the zoom region. */
  const struct zoom_region_properties *zoom_props = s->zoom_props;
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_region_bounds[6] = {
    zoom_props->region_bounds[0], zoom_props->region_bounds[1],
    zoom_props->region_bounds[2], zoom_props->region_bounds[3],
    zoom_props->region_bounds[4], zoom_props->region_bounds[5]};
  const int buffer_cell_offset = zoom_props->buffer_cell_offset;
  const double  buffer_bounds[6] = {
    zoom_props->buffer_bounds[0], zoom_props->buffer_bounds[1],
    zoom_props->buffer_bounds[2], zoom_props->buffer_bounds[3],
    zoom_props->buffer_bounds[4], zoom_props->buffer_bounds[5]};

  /* Here we go down the heirarchy to get the cell_id, it's marginally slower
   * but guarantees that void cells are handled properly. */
  
  /* Get the background cell ijk coordinates. */
  const int bkg_i = x * s->iwidth[0];
  const int bkg_j = y * s->iwidth[1];
  const int bkg_k = z * s->iwidth[2];

  /* Which background cell is this? */
  cell_id = cell_getid(s->cdim, bkg_i, bkg_j, bkg_k) + bkg_cell_offset;

  /* If this is a void cell we are in the zoom region. */
  if (s->cells_top[cell_id].subtype == void_cell) {

    /* Which zoom TL cell are we in? */
    cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_region_bounds,
                                     x, y, z, s->zoom_props->iwidth,
                                     /*offset*/0);
    
  }

  /* If this is a neighbour void cell we are in the buffer cells.
   * Otherwise, It's a legitimate background cell, and we'll return it. */
  else if (s->cells_top[cell_id].subtype == empty) {

    /* Which buffer TL cell are we in? */
    cell_id = cell_getid_with_bounds(s->zoom_props->buffer_cdim, buffer_bounds,
                                     x, y, z, s->zoom_props->buffer_iwidth,
                                     buffer_cell_offset);

    /* Here we need to check if this is the void buffer cell.
     * Otherwise, It's a legitimate buffer cell, and we'll return it. */
    if (s->cells_top[cell_id].subtype == void_cell) {
      
      /* Which zoom TL cell are we in? */
      cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_region_bounds,
                                       x, y, z, s->zoom_props->iwidth,
                                       /*offset*/0);
      
    } 
  }

#ifdef SWIFT_DEBUG_CHECKS
    if (cell_id < 0 || cell_id >= s->nr_cells)
      error("cell_id out of range: %i (%f %f %f)", cell_id, x, y, z);

    if (s->cells_top[cell_id].subtype == void_cell ||
        s->cells_top[cell_id].subtype == empty)
      error("void cell has been given a particle! (c->type=%d, "
            "x=%f, y=%f, z=%f)",
            s->cells_top[cell_id].type, x, y, z);
#endif

  return cell_id;
#else
  error("Using cell_getid_zoom but compiled without zoom regions enabled!");
#endif
}

#ifdef WITH_ZOOM_REGION
#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Run through all cells and ensure they have the correct cell type and
 * width for their position in s->cells_top.
 *
 * @param s The space.
 */
static void debug_cell_type(struct space *s) {

  /* Get the cells array and cell properties */
  struct cell *cells = s->cells_top;
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;
  const int buffer_offset = s->zoom_props->buffer_cell_offset;
  const double *zoom_width = s->zoom_props->width;
  const double *width = s->width;

  /* Loop over all cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Check cell type */
    if (cid < bkg_cell_offset && cells[cid].type != zoom)
      error(
          "Cell has the wrong cell type for it's array position (cid=%d, "
          "c->type=%d, "
          "s->zoom_props->bkg_cell_offset=%d)",
          cid, cells[cid].type, bkg_cell_offset);
    if (cid >= bkg_cell_offset && cells[cid].type == zoom)
      error(
          "Cell has the wrong cell type for it's array position (cid=%d, "
          "c->type=%d, "
          "s->zoom_props->bkg_cell_offset=%d)",
          cid, cells[cid].type, bkg_cell_offset);

    /* Check cell widths */
    for (int ijk = 0; ijk < 3; ijk++) {
      if (cid < bkg_cell_offset && cells[cid].width[ijk] != zoom_width[ijk])
        error(
            "Cell has the wrong cell width for it's array position (cid=%d, "
            "c->type=%d, "
            "s->zoom_props->bkg_cell_offset=%d, c->width=[%f %f %f], "
            "s->zoom_props->width=[%f %f %f])",
            cid, cells[cid].type, bkg_cell_offset, cells[cid].width[0],
            cells[cid].width[1], cells[cid].width[2], s->zoom_props->width[0],
            s->zoom_props->width[1], s->zoom_props->width[2]);
      if ((cid >= bkg_cell_offset && cid < buffer_offset) &&
          cells[cid].width[ijk] != width[ijk])
        error(
            "Cell has the wrong cell width for it's array position (cid=%d, "
            "c->type=%d, "
            "s->zoom_props->bkg_cell_offset=%d, c->width=[%f %f %f], "
            "s->zoom_props->width=[%f %f %f])",
            cid, cells[cid].type, bkg_cell_offset, cells[cid].width[0],
            cells[cid].width[1], cells[cid].width[2], s->zoom_props->width[0],
            s->zoom_props->width[1], s->zoom_props->width[2]);
    }
  }

  if (s->zoom_props->with_buffer_cells) {

    /* Loop over natural cells and ensure the cell boundaries and buffer
     * boundaries line up. */
    int found_i = 0;
    int found_j = 0;
    int found_k = 0;
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {
          const size_t cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;
          
          if (cells[cid].loc[0] == s->zoom_props->buffer_bounds[0])
            found_i = 1;
          
          if (cells[cid].loc[1] == s->zoom_props->buffer_bounds[2])
            found_j = 1;
          
          if (cells[cid].loc[2] == s->zoom_props->buffer_bounds[4])
            found_k = 1;
          
        }
      }
    }
    
    /* Report if we didn't find matching boundaries. */
    if (!found_i || !found_j || !found_k)
      error("The background cell and buffer region edges don't match!");
  }
}
#endif
#endif

/**
 * @brief Compute the extent/bounds of the zoom region using the high-res DM
 * particles.
 *
 * The min/max [x,y,z] for each particle is found, and the CoM of these
 * particles is computed.
 *
 * @param s The space.
 * @param nr_nodes The number of ranks.
 * @param verbose Are we talking?
 */
void construct_zoom_region(struct space *s, int nr_nodes, int verbose) {
#ifdef WITH_ZOOM_REGION
  /* Get the width of the zoom region, zoom dims are equal. */
  const double zoom_dim = s->zoom_props->dim[0];

  /* Let's set what we know about the zoom region. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->width[ijk] = zoom_dim / s->zoom_props->cdim[ijk];
    s->zoom_props->iwidth[ijk] = 1 / s->zoom_props->width[ijk];
    s->zoom_props->dim[ijk] = zoom_dim;
  }

  /* Overwrite the minimum allowed zoom cell width. */
  s->zoom_props->cell_min = 0.99 * zoom_dim / s->zoom_props->cdim[0];

  /* Resize the top level cells in the space. */
  const double dmax = max3(s->dim[0], s->dim[1], s->dim[2]);
  s->cell_min = 0.99 * dmax / s->cdim[0];

  /* Store cell number information. */
  s->zoom_props->bkg_cell_offset =
      s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2];
  s->zoom_props->nr_zoom_cells = s->zoom_props->bkg_cell_offset;
  s->zoom_props->nr_bkg_cells = s->cdim[0] * s->cdim[1] * s->cdim[2];
  s->zoom_props->buffer_cell_offset =
    s->zoom_props->bkg_cell_offset + s->zoom_props->nr_bkg_cells;
  s->zoom_props->nr_buffer_cells =
    s->zoom_props->buffer_cdim[0] * s->zoom_props->buffer_cdim[1] *
    s->zoom_props->buffer_cdim[2];

  /* Lets report what we have constructed. */
  if (verbose) {
    message(
        "set cell dimensions to zoom_cdim=[%d %d %d]"
        " background_cdim=[%d %d %d] buffer_cdim=[%d %d %d]",
        s->zoom_props->cdim[0], s->zoom_props->cdim[1], s->zoom_props->cdim[2],
        s->cdim[0], s->cdim[1], s->cdim[2], s->zoom_props->buffer_cdim[0],
        s->zoom_props->buffer_cdim[1], s->zoom_props->buffer_cdim[2]);
    message("nr_zoom_cells/bkg_cell_offset: %d nr_bkg_cells: %d "
            "nr_buffer_cells: %d",
            s->zoom_props->nr_zoom_cells, s->zoom_props->nr_bkg_cells,
            s->zoom_props->nr_buffer_cells);
    message("zoom_boundary: [%.2f-%.2f %.2f-%.2f %.2f-%.2f]",
            s->zoom_props->region_bounds[0], s->zoom_props->region_bounds[1],
            s->zoom_props->region_bounds[2], s->zoom_props->region_bounds[3],
            s->zoom_props->region_bounds[4], s->zoom_props->region_bounds[5]);
    message("buffer_boundary: [%.2f-%.2f %.2f-%.2f %.2f-%.2f]",
            s->zoom_props->buffer_bounds[0], s->zoom_props->buffer_bounds[1],
            s->zoom_props->buffer_bounds[2], s->zoom_props->buffer_bounds[3],
            s->zoom_props->buffer_bounds[4], s->zoom_props->buffer_bounds[5]);
    message(
        "zoom_region_dim: [%.2f %.2f %.2f] tl_cell_width: [%.2f %.2f %.2f] "
        "zoom_cell_width: [%.2f %.2f %.2f] buffer_cell_width: [%.2f %.2f %.2f]",
        s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2],
        s->width[0], s->width[1], s->width[2], s->zoom_props->width[0],
        s->zoom_props->width[1], s->zoom_props->width[2],
        s->zoom_props->buffer_width[0], s->zoom_props->buffer_width[1],
        s->zoom_props->buffer_width[2]);
  }

  /* Check we have enough cells for periodicity. */
  if (s->periodic && (s->cdim[0] < 3 || s->cdim[1] < 3 || s->cdim[2] < 3))
    error(
        "Must have at least 3 cells in each spatial dimension when periodicity "
        "is switched on (cdim=%d).\nThis error is often caused by any of the "
        "followings:\n"
        " - too few particles to generate a sensible grid,\n"
        " - the initial value of 'Scheduler:max_top_level_cells' is too "
        "small,\n"
        " - the (minimal) time-step is too large leading to particles with "
        "predicted smoothing lengths too large for the box size,\n"
        " - particles with velocities so large that they move by more than two "
        "box sizes per time-step.\n", s->cdim[0]);

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  /* What is the angular extent of a background cell? Here we use the,
   * number of nodes but double it so the number of wedges will be
   * (2 * nr_nodes) ** 2 ensuring the background is well divided without
   * being bogged down with the number of wedges. */
  double cell_angular_size = M_PI / 2 / (2 * nr_nodes);

  /* The number of slices in theta and phi. */
  s->zoom_props->theta_nslices = floor(2 * M_PI / cell_angular_size);
  s->zoom_props->phi_nslices = floor(M_PI / cell_angular_size);

  /* Calculate the size of a wedge in theta and phi. */
  s->zoom_props->theta_width = 2 * M_PI / s->zoom_props->theta_nslices;
  s->zoom_props->phi_width = M_PI / s->zoom_props->phi_nslices;

  /* How many wedges do we have in total? */
  s->zoom_props->nwedges =
    s->zoom_props->theta_nslices * s->zoom_props->phi_nslices;

  /* Allocate the wedge edge counts. */
  if (swift_memalign("wedge_edge_counts",
                     (void **)&s->zoom_props->nr_wedge_edges,
                     SWIFT_STRUCT_ALIGNMENT,
                     s->zoom_props->nwedges * sizeof(int)) != 0)
    error("Failed to allocate the number of wedge edges.");
  bzero(s->zoom_props->nr_wedge_edges, s->zoom_props->nwedges * sizeof(int));

  /* Allocate the wedge edge counts. */
  if (swift_memalign("wedge_edge_starts",
                     (void **)&s->zoom_props->wedge_edges_start,
                     SWIFT_STRUCT_ALIGNMENT,
                     s->zoom_props->nwedges * sizeof(int)) != 0)
    error("Failed to allocate the start pointer for wedge edges.");
  bzero(s->zoom_props->wedge_edges_start, s->zoom_props->nwedges * sizeof(int));
  
#endif

#endif
}

/**
 * @brief Build the TL cells, with a zoom region.
 *
 * This replaces the loop in space_regrid when running with a zoom region.
 *
 * Construct an additional set of TL "zoom" cells embedded within the TL cell
 * structure with the dimensions of each cell structure being the same (with
 * differing widths).
 *
 * Therefore the new TL cell structure is 2*cdim**3, with the "natural" TL cells
 * occupying the first half of the TL cell list, and the "zoom" TL cells
 * ocupying the second half.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void construct_tl_cells_with_zoom_region(
    struct space *s, const int *cdim, const float dmin,
    const integertime_t ti_current, struct gravity_props *gravity_properties,
    int verbose) {
#ifdef WITH_ZOOM_REGION

  /* Get some zoom region properties */
  const float dmin_zoom = min3(s->zoom_props->width[0], s->zoom_props->width[1],
                               s->zoom_props->width[2]);
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;
  const double zoom_region_bounds[6] = {
      s->zoom_props->region_bounds[0], s->zoom_props->region_bounds[1],
      s->zoom_props->region_bounds[2], s->zoom_props->region_bounds[3],
      s->zoom_props->region_bounds[4], s->zoom_props->region_bounds[5]};

  struct cell *restrict c;

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(s->zoom_props->cdim, i, j, k);

        /* Create the zoom cell and it's multipoles */
        c = &s->cells_top[cid];
        c->loc[0] = (i * s->zoom_props->width[0]) + zoom_region_bounds[0];
        c->loc[1] = (j * s->zoom_props->width[1]) + zoom_region_bounds[2];
        c->loc[2] = (k * s->zoom_props->width[2]) + zoom_region_bounds[4];
        c->width[0] = s->zoom_props->width[0];
        c->width[1] = s->zoom_props->width[1];
        c->width[2] = s->zoom_props->width[2];
        const size_t parent_cid =
          cell_getid(cdim,
                     (int)(c->loc[0] + (c->width[0] / 2)) / s->width[0],
                     (int)(c->loc[1] + (c->width[1] / 2)) / s->width[1],
                     (int)(c->loc[2] + (c->width[2] / 2)) / s->width[2]) +
          bkg_cell_offset;
        c->parent_bkg_cid = parent_cid;
        if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
        c->type = zoom;
        c->subtype = none;
        c->dmin = dmin_zoom;
        c->depth = 0;
        c->split = 0;
        c->hydro.count = 0;
        c->grav.count = 0;
        c->stars.count = 0;
        c->sinks.count = 0;
        c->top = c;
        c->super = c;
        c->hydro.super = c;
        c->grav.super = c;
        c->hydro.ti_old_part = ti_current;
        c->grav.ti_old_part = ti_current;
        c->stars.ti_old_part = ti_current;
        c->sinks.ti_old_part = ti_current;
        c->black_holes.ti_old_part = ti_current;
        c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
        c->mpi.tag = -1;
        c->mpi.recv = NULL;
        c->mpi.send = NULL;
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
          c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif
      }
    }
  }
  
  /* Loop over natural cells and set locations and initial values */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        const size_t cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;

        /* Natural top level cells. */
        c = &s->cells_top[cid];
        c->loc[0] = i * s->width[0];
        c->loc[1] = j * s->width[1];
        c->loc[2] = k * s->width[2];
        c->width[0] = s->width[0];
        c->width[1] = s->width[1];
        c->width[2] = s->width[2];
        c->dmin = dmin;
        c->parent_bkg_cid = cid;
        if (s->with_self_gravity)
          c->grav.multipole = &s->multipoles_top[cid];
        c->type = bkg;
        c->subtype = none;
        c->depth = 0;
        c->split = 0;
        c->hydro.count = 0;
        c->grav.count = 0;
        c->stars.count = 0;
        c->sinks.count = 0;
        c->top = c;
        c->super = c;
        c->hydro.super = c;
        c->grav.super = c;
        c->hydro.ti_old_part = ti_current;
        c->grav.ti_old_part = ti_current;
        c->stars.ti_old_part = ti_current;
        c->sinks.ti_old_part = ti_current;
        c->black_holes.ti_old_part = ti_current;
        c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
        c->mpi.tag = -1;
        c->mpi.recv = NULL;
        c->mpi.send = NULL;
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
          c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif

        /* Assign the cell type. */
        if (s->zoom_props->with_buffer_cells &&
            cell_inside_buffer_region(c, s)) {
          c->subtype = empty;
          s->zoom_props->nr_empty_cells++;
        } 
      }
    }
  }

  /* If we have a buffer region create buffer cells. */
  if (s->zoom_props->with_buffer_cells) {

    if (verbose)
      message("%i background cells are in the buffer region",
              s->zoom_props->nr_empty_cells);

    /* Get relevant information. */
    const float dmin_buffer = min3(s->zoom_props->buffer_width[0],
                                 s->zoom_props->buffer_width[1],
                                 s->zoom_props->buffer_width[2]);
    const int buffer_offset = s->zoom_props->buffer_cell_offset;
    const double buffer_bounds[6] = {
      s->zoom_props->buffer_bounds[0], s->zoom_props->buffer_bounds[1],
      s->zoom_props->buffer_bounds[2], s->zoom_props->buffer_bounds[3],
      s->zoom_props->buffer_bounds[4], s->zoom_props->buffer_bounds[5]};

    /* Loop over buffer cells and set locations and initial values */
    for (int i = 0; i < s->zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < s->zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < s->zoom_props->buffer_cdim[2]; k++) {
          const size_t cid =
            cell_getid(s->zoom_props->buffer_cdim, i, j, k) + buffer_offset;
          
          /* Buffer top level cells. */
          c = &s->cells_top[cid];
          c->loc[0] = (i * s->zoom_props->buffer_width[0]) + buffer_bounds[0];
          c->loc[1] = (j * s->zoom_props->buffer_width[1]) + buffer_bounds[2];
          c->loc[2] = (k * s->zoom_props->buffer_width[2]) + buffer_bounds[4];
          c->width[0] = s->zoom_props->buffer_width[0];
          c->width[1] = s->zoom_props->buffer_width[1];
          c->width[2] = s->zoom_props->buffer_width[2];
          c->dmin = dmin_buffer;
          const size_t parent_cid =
            cell_getid(s->zoom_props->buffer_cdim,
                       (int)((c->loc[0] + (c->width[0] / 2)) -
                             buffer_bounds[0]) * s->zoom_props->buffer_iwidth[0],
                       (int)((c->loc[1] + (c->width[1] / 2)) -
                             buffer_bounds[2]) * s->zoom_props->buffer_iwidth[1],
                       (int)((c->loc[2] + (c->width[2] / 2)) -
                             buffer_bounds[4]) *
                       s->zoom_props->buffer_iwidth[2]) + buffer_offset;
          c->parent_bkg_cid = parent_cid;
          if (s->with_self_gravity)
            c->grav.multipole = &s->multipoles_top[cid];
          c->type = buffer;
          c->subtype = none;
          c->depth = 0;
          c->split = 0;
          c->hydro.count = 0;
          c->grav.count = 0;
          c->stars.count = 0;
          c->sinks.count = 0;
          c->top = c;
          c->super = c;
          c->hydro.super = c;
          c->grav.super = c;
          c->hydro.ti_old_part = ti_current;
          c->grav.ti_old_part = ti_current;
          c->stars.ti_old_part = ti_current;
          c->sinks.ti_old_part = ti_current;
          c->black_holes.ti_old_part = ti_current;
          c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
          c->mpi.tag = -1;
          c->mpi.recv = NULL;
          c->mpi.send = NULL;
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
          c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
          cell_assign_top_level_cell_index(c, s);
#endif
        }
      }
    }
  }

  /* Now find what cells contain the zoom region. */
  find_void_cells(s, verbose);

  /* Now find what cells neighbour the zoom region. */
  find_neighbouring_cells(s, gravity_properties, verbose);

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  /* Find the number of edges we will need for the domain decomp. */
  find_vertex_edges(s, verbose);

  if (verbose)
    message("%i vertex 'edges' found in total", s->zoom_props->nr_edges);

#endif

#ifdef SWIFT_DEBUG_CHECKS
  
  /* Lets check all the cells are in the right place with the correct widths */
  debug_cell_type(s);
#endif

#endif
}

/**
 * @brief Find what TL cells contain the zoom region.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_void_cells(struct space *s, const int verbose) {
#ifdef WITH_ZOOM_REGION

  /* Get the right cell cdim. */
  int cdim[3];
  if (s->zoom_props->with_buffer_cells) {
    cdim[0] = s->zoom_props->buffer_cdim[0];
    cdim[1] = s->zoom_props->buffer_cdim[1];
    cdim[2] = s->zoom_props->buffer_cdim[2];
  } else {
    cdim[0] = s->cdim[0];
    cdim[1] = s->cdim[1];
    cdim[2] = s->cdim[2];
  }

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;
  
  /* Get the right offset and the number of cells we're dealing with. */
  int offset;
  int ncells;
  if (s->zoom_props->with_buffer_cells) {
    offset = s->zoom_props->buffer_cell_offset;
    ncells = s->zoom_props->nr_buffer_cells;
  } else {
    offset = s->zoom_props->bkg_cell_offset;
    ncells = s->zoom_props->nr_bkg_cells;
  }

  /* Allocate the indices of void cells */
  /* TODO: We can know ahead of time how many void cells there are. We don't
   * need to allocate this many.  */
  int void_count = 0;
  if (swift_memalign("void_cells_top",
                     (void **)&s->zoom_props->void_cells_top,
                     SWIFT_STRUCT_ALIGNMENT, ncells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->void_cells_top, ncells * sizeof(int));

  /* Loop over natural cells and find cells containing the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
        const size_t cid = cell_getid(cdim, i, j, k) + offset;

        /* Label this background cell. */
        if (cell_contains_zoom_region(&cells[cid], s)) {
          cells[cid].subtype = void_cell;
          s->zoom_props->void_cells_top[void_count++] = cid;
        }
      }
    }
  }

  /* Store the number of neighbour cells */
  s->zoom_props->nr_void_cells = void_count;

  if (void_count == 0)
    error("No void cells were found! (nr_buffer_cells=%d)",
          s->zoom_props->nr_buffer_cells);

  if (verbose)
    message("%i cells contain the zoom region", void_count);

#endif
}

/**
 * @brief Find what TL cells surround the zoom region.
 *
 * When interacting "natural" TL cells and "zoom" TL cells, it helps to know
 * what natural TL cells surround the zoom region. These cells then get tagged
 * as "neighbour" cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_neighbouring_cells(struct space *s,
                             struct gravity_props *gravity_properties,
                             const int verbose) {
#ifdef WITH_ZOOM_REGION
  
  /* Get the right cell cdim. */
  int cdim[3];
  double iwidth[3];
  if (s->zoom_props->with_buffer_cells) {
    cdim[0] = s->zoom_props->buffer_cdim[0];
    cdim[1] = s->zoom_props->buffer_cdim[1];
    cdim[2] = s->zoom_props->buffer_cdim[2];
    iwidth[0] = s->zoom_props->buffer_iwidth[0];
    iwidth[1] = s->zoom_props->buffer_iwidth[1];
    iwidth[2] = s->zoom_props->buffer_iwidth[2];
  } else {
    cdim[0] = s->cdim[0];
    cdim[1] = s->cdim[1];
    cdim[2] = s->cdim[2];
    iwidth[0] = s->iwidth[0];
    iwidth[1] = s->iwidth[1];
    iwidth[2] = s->iwidth[2];
  }

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;
  
  /* Get the right offset and the number of cells we're dealing with. */
  int offset;
  int ncells;
  if (s->zoom_props->with_buffer_cells) {
    offset = s->zoom_props->buffer_cell_offset;
    ncells = s->zoom_props->nr_buffer_cells;
  } else {
    offset = s->zoom_props->bkg_cell_offset;
    ncells = s->zoom_props->nr_bkg_cells;
  }

  /* Get gravity mesh distance. */
  const double max_distance = gravity_properties->r_s
    * gravity_properties->r_cut_max_ratio;

  /* At this point we can only define neighbour cells by cell properties,
   * leaving the fancy gravity distance criterion for task creation later.
   * Here we just make sure all possible neighbour cells are flagged
   * as such. */
  
  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
  const int delta_cells = ceil(max_distance * max3(iwidth[0],
                                                   iwidth[1],
                                                   iwidth[2])) + 1;

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Let's be verbose about this choice */
  if (verbose)
    message(
        "Looking for neighbouring natural cells up to %d natural top-level "
        "cells away from the zoom region (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Allocate the indices of neighbour background cells */
  if (swift_memalign("neighbour_cells_top",
                     (void **)&s->zoom_props->neighbour_cells_top,
                     SWIFT_STRUCT_ALIGNMENT, ncells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->neighbour_cells_top,  ncells * sizeof(int));

  int neighbour_count = 0;

  /* Loop over natural cells and find cells neighbouring the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

         /* Get this cell. */
        const size_t cid = cell_getid(cdim, i, j, k) + offset;
        struct cell *ci = &cells[cid];

        /* Skip non-void cells. */
        if (ci->subtype != void_cell) continue;

        /* Loop over every other cell within (Manhattan) range delta */
        for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

          /* Escape beyond range */
          if (ii < 0 || ii >= cdim[0]) continue;

          for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

            /* Escape beyond range */
            if (jj < 0 || jj >= cdim[1]) continue;

            for (int kk = k - delta_m; kk <= k + delta_p; kk++) {
              
              /* Escape beyond range */
              if (kk < 0 || kk >= cdim[2]) continue;

              /* Get this cell. */
              const int cjd = cell_getid(cdim, ii, jj, kk) + offset;

              /* Ensure this neighbour isn't a void cell or an already
               * counted neighbour. */
              if (cells[cjd].subtype != void_cell &&
                  cells[cjd].subtype != neighbour) {

                /* Record that we've found a neighbour. */
                cells[cjd].subtype = neighbour;
                s->zoom_props->neighbour_cells_top[neighbour_count++] = cjd;
              }
            } /* neighbour k loop */
          } /* neighbour j loop */
        } /* neighbour i loop */
      } /* natural k loop */
    } /* natural j loop */
  } /* natural i loop */

  /* Store the number of neighbour cells */
  s->zoom_props->nr_neighbour_cells = neighbour_count;

  if (verbose)
    message("%i cells neighbouring the zoom region", neighbour_count);
#endif
}

/**
 * @brief For METIS we need to work out how many edges each vertex (cell) has.
 *
 * @param verbose The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
void link_zoom_to_void(struct space *s, struct cell *c) {
  
  /* Set up some useful information. */
  double zoom_loc[3];

  /* We need to ensure this bottom level isn't treated like a
   * normal split cell since it's linked into top level "progeny". */
  /* c->split = 0; */

  /* Loop over the 8 progeny cells which are now the zoom cells. */
  for (int k = 0; k < 8; k++) {

    /* Establish the location of the fake progeny cell. */
    zoom_loc[0] = c->loc[0];
    zoom_loc[1] = c->loc[1];
    zoom_loc[2] = c->loc[2];
    if (k & 4) zoom_loc[0] += s->zoom_props->width[0];
    if (k & 2) zoom_loc[1] += s->zoom_props->width[1];
    if (k & 1) zoom_loc[2] += s->zoom_props->width[2];

    /* Which zoom cell are we in? */
    int cid = cell_getid_pos(s,
                             zoom_loc[0] + (s->zoom_props->width[0] / 2),
                             zoom_loc[1] + (s->zoom_props->width[0] / 2),
                             zoom_loc[2] + (s->zoom_props->width[0] / 2));

    /* Get the zoom cell. */
    struct cell *zoom_cell = &s->cells_top[cid];
    
    /* Link this zoom cell into the void cell hierarchy. */
    c->progeny[k] = zoom_cell;

    /* Flag this void cell "progeny" as the zoom cell's void cell parent. */
    zoom_cell->void_parent = c;
    
  }

  /* Interact the zoom cell multipoles with this cell. */
  if (s->with_self_gravity) {

    /* Reset everything */
    gravity_reset(c->grav.multipole);

    /* Compute CoM and bulk velocity from all progenies */
    double CoM[3] = {0., 0., 0.};
    double vel[3] = {0., 0., 0.};
    float max_delta_vel[3] = {0.f, 0.f, 0.f};
    float min_delta_vel[3] = {0.f, 0.f, 0.f};
    double mass = 0.;
    
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        const struct gravity_tensors *m = c->progeny[k]->grav.multipole;
        
        mass += m->m_pole.M_000;
        
        CoM[0] += m->CoM[0] * m->m_pole.M_000;
        CoM[1] += m->CoM[1] * m->m_pole.M_000;
        CoM[2] += m->CoM[2] * m->m_pole.M_000;

        vel[0] += m->m_pole.vel[0] * m->m_pole.M_000;
        vel[1] += m->m_pole.vel[1] * m->m_pole.M_000;
        vel[2] += m->m_pole.vel[2] * m->m_pole.M_000;

        max_delta_vel[0] = max(m->m_pole.max_delta_vel[0], max_delta_vel[0]);
        max_delta_vel[1] = max(m->m_pole.max_delta_vel[1], max_delta_vel[1]);
        max_delta_vel[2] = max(m->m_pole.max_delta_vel[2], max_delta_vel[2]);
        
        min_delta_vel[0] = min(m->m_pole.min_delta_vel[0], min_delta_vel[0]);
        min_delta_vel[1] = min(m->m_pole.min_delta_vel[1], min_delta_vel[1]);
        min_delta_vel[2] = min(m->m_pole.min_delta_vel[2], min_delta_vel[2]);
      }
    }
    
    /* Final operation on the CoM and bulk velocity */
    const double inv_mass = 1. / mass;
    c->grav.multipole->CoM[0] = CoM[0] * inv_mass;
    c->grav.multipole->CoM[1] = CoM[1] * inv_mass;
    c->grav.multipole->CoM[2] = CoM[2] * inv_mass;
    c->grav.multipole->m_pole.vel[0] = vel[0] * inv_mass;
    c->grav.multipole->m_pole.vel[1] = vel[1] * inv_mass;
    c->grav.multipole->m_pole.vel[2] = vel[2] * inv_mass;
    
    /* Min max velocity along each axis */
    c->grav.multipole->m_pole.max_delta_vel[0] = max_delta_vel[0];
    c->grav.multipole->m_pole.max_delta_vel[1] = max_delta_vel[1];
    c->grav.multipole->m_pole.max_delta_vel[2] = max_delta_vel[2];
    c->grav.multipole->m_pole.min_delta_vel[0] = min_delta_vel[0];
    c->grav.multipole->m_pole.min_delta_vel[1] = min_delta_vel[1];
    c->grav.multipole->m_pole.min_delta_vel[2] = min_delta_vel[2];
    
    /* Now shift progeny multipoles and add them up */
    struct multipole temp;
    double r_max = 0.;
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        const struct cell *cp = c->progeny[k];
        const struct multipole *m = &cp->grav.multipole->m_pole;
        
        /* Contribution to multipole */
        gravity_M2M(&temp, m, c->grav.multipole->CoM,
                    cp->grav.multipole->CoM);
        gravity_multipole_add(&c->grav.multipole->m_pole, &temp);
        
        /* Upper limit of max CoM<->gpart distance */
        const double dx =
          c->grav.multipole->CoM[0] - cp->grav.multipole->CoM[0];
        const double dy =
          c->grav.multipole->CoM[1] - cp->grav.multipole->CoM[1];
        const double dz =
          c->grav.multipole->CoM[2] - cp->grav.multipole->CoM[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        r_max = max(r_max, cp->grav.multipole->r_max + sqrt(r2));
      }
    }
    
    /* Alternative upper limit of max CoM<->gpart distance */
    const double dx =
      c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
      ? c->grav.multipole->CoM[0] - c->loc[0]
      : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
    const double dy =
      c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
      ? c->grav.multipole->CoM[1] - c->loc[1]
      : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
    const double dz =
      c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
      ? c->grav.multipole->CoM[2] - c->loc[2]
      : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];
    
    /* Take minimum of both limits */
    c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));
    
    /* Store the value at rebuild time */
    c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
    c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
    c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
    c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];
    
    /* Compute the multipole power */
    gravity_multipole_compute_power(&c->grav.multipole->m_pole);
    
  } /* Deal with gravity */
}

/**
 * @brief For METIS we need to work out how many edges each vertex (cell) has.
 * 
 *
 * @param verbose The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
void find_vertex_edges(struct space *s, const int verbose) {
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  /* Get some useful constants. */
  const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                            s->zoom_props->cdim[2]};

  /* Initialise edge count. */
  s->zoom_props->nr_edges = 0;
  int iedge;
  int all_iedge = 0;

  /* Find adjacency arrays for cells and wedges. */
  if (s->zoom_props->use_bkg_wedges) {
    iedge = 0;
    edge_loop(zoom_cdim, 0, s, /*adjncy*/ NULL, /*xadj*/ NULL,
              /*counts*/ NULL, /*edges*/ NULL, &iedge);
    all_iedge = iedge;
  } else if (s->zoom_props->separate_decomps){
    /* Otherwise, we need to find the edges in each individual level. */

    /* Zoom */
    iedge = 0;
    edge_loop(zoom_cdim, 0, s, /*adjncy*/ NULL, /*xadj*/ NULL,
              /*counts*/ NULL, /*edges*/ NULL, &iedge);
    all_iedge += iedge;

    /* Buffer, if we have them */
    if (s->zoom_props->with_buffer_cells) {
      iedge = 0;
      edge_loop(s->zoom_props->buffer_cdim, s->zoom_props->buffer_cell_offset,
                s, /*adjncy*/ NULL, /*xadj*/ NULL,
                /*counts*/ NULL, /*edges*/ NULL, &iedge);
      all_iedge += iedge;
    }
    
    /* Background */
    iedge = 0;
    edge_loop(s->cdim, s->zoom_props->bkg_cell_offset,
              s, /*adjncy*/ NULL, /*xadj*/ NULL,
              /*counts*/ NULL, /*edges*/ NULL, &iedge);
    all_iedge += iedge;
  } else {

    /* Otherwise, we only need zoom edges */
    iedge = 0;
    edge_loop(zoom_cdim, 0, s, /*adjncy*/ NULL, /*xadj*/ NULL,
              /*counts*/ NULL, /*edges*/ NULL, &iedge);
    all_iedge += iedge;
    
  }

  /* Set the total number of edges. */
  s->zoom_props->nr_edges = all_iedge;

#ifdef SWIFT_DEBUG_CHECKS
  for (int cid = 0; cid < s->zoom_props->nr_zoom_cells; cid++) {
    if (s->cells_top[cid].nr_vertex_edges == 0)
      error("Cell (%d) has no edges! (c->type=%d)", cid,
            s->cells_top[cid].type);
  }
#endif
#endif
}

/**
 * @brief Minimum distance between two TL cells with different sizes.
 *
 * @param ci, cj The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
double cell_min_dist2_diff_size(const struct cell *restrict ci,
                                const struct cell *restrict cj,
                                int periodic, const double dim[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->width[0] == cj->width[0]) error("x cells of same size!");
  if (ci->width[1] == cj->width[1]) error("y cells of same size!");
  if (ci->width[2] == cj->width[2]) error("z cells of same size!");
#endif

  /* We need to check if we need to consider periodicity since only
   * background cells are periodic. */
  if (ci->type == bkg || cj->type == bkg)
    periodic = periodic;
  else
    periodic = 0;

  const double cix = ci->loc[0] + ci->width[0] / 2.;
  const double ciy = ci->loc[1] + ci->width[1] / 2.;
  const double ciz = ci->loc[2] + ci->width[2] / 2.;

  const double cjx = cj->loc[0] + cj->width[0] / 2.;
  const double cjy = cj->loc[1] + cj->width[1] / 2.;
  const double cjz = cj->loc[2] + cj->width[2] / 2.;

  const double diag_ci2 = ci->width[0] * ci->width[0] +
                          ci->width[1] * ci->width[1] +
                          ci->width[2] * ci->width[2];
  const double diag_cj2 = cj->width[0] * cj->width[0] +
                          cj->width[1] * cj->width[1] +
                          cj->width[2] * cj->width[2];

  /* Get the distance between the cells */
  double dx = cix - cjx;
  double dy = ciy - cjy;
  double dz = ciz - cjz;

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  /* Minimal distance between any 2 particles in the two cells */
  const double dist2 = r2 - (diag_ci2 / 2. + diag_cj2 / 2.);

  return dist2;
}

/**
 * @brief Minimum distance between two TL cells.
 *
 * Generic wrapper, don't know if the TL cells are the same size or not at time
 * of calling.
 *
 * @param ci, cj The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
double cell_min_dist2(const struct cell *restrict ci,
                      const struct cell *restrict cj, int periodic,
                      const double dim[3]) {
#ifdef WITH_ZOOM_REGION
  double dist2;

  /* Two TL cells with the same size. */
  if (ci->width[0] == cj->width[0]) {
    dist2 = cell_min_dist2_same_size(ci, cj, periodic, dim);
  }
  /* Two cells with different sizes. */
  else {
    dist2 = cell_min_dist2_diff_size(ci, cj, periodic, dim);
  }

  return dist2;
#else
  return cell_min_dist2_same_size(ci, cj, periodic, dim);
#endif
}

#ifdef WITH_ZOOM_REGION

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for the natural/background cells.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 */
void engine_make_self_gravity_tasks_mapper_natural_cells(void *map_data,
                                                         int num_elements,
                                                         void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Some info about the zoom domain */
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[bkg_cell_offset].width[0],
      s->max_softening, s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  const int delta = max((int)(sqrt(3) * distance /
                          cells[bkg_cell_offset].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (periodic) {
    if (delta >= cdim[0] / 2) {
      if (cdim[0] % 2 == 0) {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2 - 1;
      } else {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2;
      }
    }
  } else {
    if (delta > cdim[0]) {
      delta_m = cdim[0];
      delta_p = cdim[0];
    }
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index, including background cell offset. */
    const int cid = (size_t)(map_data) + ind + bkg_cell_offset;

    /* Integer indices of the cell in the top-level grid */
    const int i = (cid - bkg_cell_offset) / (cdim[1] * cdim[2]);
    const int j = ((cid - bkg_cell_offset) / cdim[2]) % cdim[1];
    const int k = (cid - bkg_cell_offset) % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles and void cells. */
    if (ci->grav.count == 0 || ci->subtype == void_cell) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav_bkg, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Escape if non-periodic and beyond range */
      if (!periodic && (ii < 0 || ii >= cdim[0])) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

	/* Escape if non-periodic and beyond range */
        if (!periodic && (jj < 0 || jj >= cdim[1])) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Escape if non-periodic and beyond range */
          if (!periodic && (kk < 0 || kk >= cdim[2])) continue;

          /* Apply periodic BC (not harmful if not using periodic BC) */
	  const int iii = (ii + cdim[0]) % cdim[0];
	  const int jjj = (jj + cdim[1]) % cdim[1];
	  const int kkk = (kk + cdim[2]) % cdim[2];

          /* Get the second cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_cell_offset;
          struct cell *cj = &cells[cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 || cj->subtype == void_cell ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

#ifdef WITH_MPI
          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");
#endif

          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (periodic && min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0)) {

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav_bkg, 0, 0,
                              ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
            /* Ensure both cells are background cells */
            if (ci->type == zoom || cj->type == zoom) {
              error(
                  "Cell %d and cell %d are not background cells! "
                  "(ci->type=%d, cj->type=%d)",
                  cid, cjd, ci->type, cj->type);
            }
#ifdef WITH_MPI

            /* Let's cross-check that we had a proxy for that cell */
            if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[cj->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", cj->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == cj) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cjd);
            } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[ci->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", ci->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == ci) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cid);
            }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for the natural/background cells.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 */
void engine_make_self_gravity_tasks_mapper_buffer_cells(void *map_data,
                                                         int num_elements,
                                                         void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->zoom_props->buffer_cdim[0],
                       s->zoom_props->buffer_cdim[1],
                       s->zoom_props->buffer_cdim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Some info about the zoom domain */
  const int buffer_offset = s->zoom_props->buffer_cell_offset;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[buffer_offset].width[0],
      s->max_softening, s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  const int delta = max((int)(sqrt(3) * distance /
                          cells[buffer_offset].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index, including background cell offset. */
    const int cid = (size_t)(map_data) + ind + buffer_offset;

    /* Integer indices of the cell in the top-level grid */
    const int i = (cid - buffer_offset) / (cdim[1] * cdim[2]);
    const int j = ((cid - buffer_offset) / cdim[2]) % cdim[1];
    const int k = (cid - buffer_offset) % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles and void cells. */
    if (ci->grav.count == 0 || ci->subtype == void_cell) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav_bkg, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Buffer cells are never periodic. */
      if (ii < 0 || ii >= cdim[0]) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

        /* Buffer cells are never periodic. */
        if (jj < 0 || jj >= cdim[1]) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Buffer cells are never periodic. */
          if (kk < 0 || kk >= cdim[2]) continue;

          /* Get the second cell */
          const int cjd = cell_getid(cdim, ii, jj, kk) + buffer_offset;
          struct cell *cj = &cells[cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 || cj->subtype == void_cell ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

#ifdef WITH_MPI
          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");
#endif

          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (periodic && min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0)) {

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav_bkg, 0, 0,
                              ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
            /* Ensure both cells are background cells */
            if (ci->type == zoom || cj->type == zoom) {
              error(
                  "Cell %d and cell %d are not background cells! "
                  "(ci->type=%d, cj->type=%d)",
                  cid, cjd, ci->type, cj->type);
            }
#ifdef WITH_MPI

            /* Let's cross-check that we had a proxy for that cell */
            if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[cj->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", cj->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == cj) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cjd);
            } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[ci->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", ci->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == ci) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cid);
            }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 */
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                       s->zoom_props->cdim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[0].width[0], s->max_softening,
      s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  const int delta = max((int)(sqrt(3) * distance /
                          cells[0].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Zoom cells are never periodic, exit if beyond zoom region */
      if (ii < 0 || ii >= cdim[0]) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

        /* Zoom cells are never periodic, exit if beyond zoom region */
        if (jj < 0 || jj >= cdim[1]) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Zoom cells are never periodic, exit if beyond zoom region */
          if (kk < 0 || kk >= cdim[2]) continue;
          
          /* Get the second cell */
          const int cjd = cell_getid(cdim, ii, jj, kk);
          struct cell *cj = &cells[cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

#ifdef WITH_MPI
          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");
#endif
          
          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0? */
          if (periodic && min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0)) {

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
                              ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
            /* Ensure both cells are zoom cells */
            if (ci->type != zoom || cj->type != zoom) {
              error(
                  "Cell %d and cell %d are not zoom cells! "
                  "(ci->type=%d, cj->type=%d)",
                  cid, cjd, ci->type, cj->type);
            }
#ifdef WITH_MPI

            /* Let's cross-check that we had a proxy for that cell */
            if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[cj->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", cj->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == cj) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cjd);
            } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[ci->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", ci->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == ci) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cid);
            }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between natural level cells
 * and zoom level cells.
 *
 * This replaces the function in engine_maketasks when running with a zoom
 region.
 *
 * - All top-cells get a self task.
 * - All pairs of differing sized cells within range according to
 *   the multipole acceptance criterion get a pair task.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.a
 * @param extra_data The #engine.

 */
void engine_make_self_gravity_tasks_mapper_zoom_bkg(
    void *map_data, int num_elements, void *extra_data) {

  /* Useful local information */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int periodic = s->periodic;

  /* Get some info about the physics */
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Get the neighbouring background cells. */
  const int nr_neighbours = s->zoom_props->nr_neighbour_cells;
  const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

  /* Get the correct task label. */
  enum task_subtypes subtype;
  if (s->zoom_props->with_buffer_cells) {
    subtype = task_subtype_grav_zoombuff;
  } else {
    subtype = task_subtype_grav_zoombkg;
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;
    
    /* Loop over every neighbouring background cells */
    for (int k = 0; k < nr_neighbours; k++) {

      /* Get the cell */
      int cjd = neighbour_cells[k];
      struct cell *cj = &cells[cjd];

      /* Avoid duplicates, empty cells and completely foreign pairs */
      if (cid >= cjd || cj->grav.count == 0 || cj->subtype == void_cell ||
          (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

#ifdef SWIFT_DEBUG_CHECKS
      /* Ensure both cells are not in the same level */
      if (ci->type == cj->type) {
          error(
              "Cell %d and cell %d are the same cell type "
              "(One should be a neighbour)! "
              "(ci->type=%d, cj->type=%d)",
              cid, cjd, ci->type, cj->type);
        }
#endif

#ifdef WITH_MPI
          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");
#endif

      /* Minimal distance between any pair of particles */
      const double min_radius2 =
          cell_min_dist2_diff_size(ci, cj, periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0 ?*/
      if (periodic && min_radius2 > max_mesh_dist2) continue;

      /* Are the cells too close for a MM interaction ? */
      if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                /*is_tree_walk=*/0)) {
        
          /* Ok, we need to add a direct pair calculation */
          scheduler_addtask(sched, task_type_pair, subtype,
                            0, 0, ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI

        /* Let's cross-check that we had a proxy for that cell */
        if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[cj->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", cj->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == cj) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cjd);
        } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[ci->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", ci->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == ci) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cid);
        }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between natural level cells
 * and zoom level cells.
 *
 * This replaces the function in engine_maketasks when running with a zoom
 region.
 *
 * - All top-cells get a self task.
 * - All pairs of differing sized cells within range according to
 *   the multipole acceptance criterion get a pair task.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.a
 * @param extra_data The #engine.

 */
void engine_make_self_gravity_tasks_mapper_buffer_bkg(
    void *map_data, int num_elements, void *extra_data) {

  /* Useful local information */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;

  /* Some info about the zoom domain */
  const int buffer_offset = s->zoom_props->buffer_cell_offset;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int periodic = s->periodic;

  /* Get some info about the physics */
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind + buffer_offset;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles and void cells */
    if (ci->grav.count == 0 || ci->subtype == void_cell) continue;
    
    /* Loop over every neighbouring background cells */
    for (int cjd = bkg_offset; cjd < buffer_offset; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

      /* Avoid empty cells, void cells and completely foreign pairs */
      if (cj->grav.count == 0 || cj->subtype == void_cell ||
          (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

#ifdef WITH_MPI
      /* Recover the multipole information */
      const struct gravity_tensors *multi_i = ci->grav.multipole;
      const struct gravity_tensors *multi_j = cj->grav.multipole;
      
      if (multi_i == NULL && ci->nodeID != nodeID)
        error("Multipole of ci was not exchanged properly via the proxies");
      if (multi_j == NULL && cj->nodeID != nodeID)
        error("Multipole of cj was not exchanged properly via the proxies");
#endif

      /* Minimal distance between any pair of particles */
      const double min_radius2 =
          cell_min_dist2(ci, cj, periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0 ?*/
      if (periodic && min_radius2 > max_mesh_dist2) continue;

      /* Are the cells too close for a MM interaction ? */
      if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                /*is_tree_walk=*/0)) {
        
          /* Ok, we need to add a direct pair calculation */
          scheduler_addtask(sched, task_type_pair, task_subtype_grav_buffbkg,
                            0, 0, ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI

        /* Let's cross-check that we had a proxy for that cell */
        if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[cj->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", cj->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == cj) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cjd);
        } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[ci->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", ci->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == ci) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cid);
        }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
      }
    }
  }
}


/**
 * @brief Constructs the top-level pair tasks for the first hydro loop over
 * neighbours
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * This replaces the function in engine_maketasks but simply removes periodicity
 * (the zoom level never has periodicity), otherwise the function is identical.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_hydroloop_tasks_mapper_with_zoom(void *map_data,
                                                  int num_elements,
                                                  void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_sinks = (e->policy & engine_policy_sinks);
  const int with_black_holes = (e->policy & engine_policy_black_holes);

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->zoom_props->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without hydro or star particles */
    if ((ci->hydro.count == 0) && (!with_stars || ci->stars.count == 0) &&
        (!with_sinks || ci->sinks.count == 0) &&
        (!with_black_holes || ci->black_holes.count == 0))
      continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_density, 0, 0, ci,
                        NULL);
    }

    /* Now loop over all the neighbours of this cell */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (iii < 0 || iii >= cdim[0]) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (jjj < 0 || jjj >= cdim[1]) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (kkk < 0 || kkk >= cdim[2]) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Is that neighbour local and does it have gas or star particles ? */
          if ((cid >= cjd) ||
              ((cj->hydro.count == 0) &&
               (!with_feedback || cj->stars.count == 0) &&
               (!with_sinks || cj->sinks.count == 0) &&
               (!with_black_holes || cj->black_holes.count == 0)) ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Construct the pair task */
          const int sid = sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
          scheduler_addtask(sched, task_type_pair, task_subtype_density, sid, 0,
                            ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI

          /* Let's cross-check that we had a proxy for that cell */
          if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

            /* Find the proxy for this node */
            const int proxy_id = e->proxy_ind[cj->nodeID];
            if (proxy_id < 0)
              error("No proxy exists for that foreign node %d!", cj->nodeID);

            const struct proxy *p = &e->proxies[proxy_id];

            /* Check whether the cell exists in the proxy */
            int n = 0;
            for (n = 0; n < p->nr_cells_in; n++)
              if (p->cells_in[n] == cj) break;
            if (n == p->nr_cells_in)
              error(
                  "Cell %d not found in the proxy but trying to construct "
                  "hydro task!",
                  cjd);
          } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

            /* Find the proxy for this node */
            const int proxy_id = e->proxy_ind[ci->nodeID];
            if (proxy_id < 0)
              error("No proxy exists for that foreign node %d!", ci->nodeID);

            const struct proxy *p = &e->proxies[proxy_id];

            /* Check whether the cell exists in the proxy */
            int n = 0;
            for (n = 0; n < p->nr_cells_in; n++)
              if (p->cells_in[n] == ci) break;
            if (n == p->nr_cells_in)
              error(
                  "Cell %d not found in the proxy but trying to construct "
                  "hydro task!",
                  cid);
          }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level self + pair tasks for the FOF loop over
 * neighbours.
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * This replaces the function in engine_maketasks but simply removes periodicity
 * (the zoom level never has periodicity), otherwise the function is identical.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_fofloop_tasks_mapper_with_zoom(void *map_data,
                                                int num_elements,
                                                void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->zoom_props->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cells is local build a self-interaction */
    if (ci->nodeID == nodeID)
      scheduler_addtask(sched, task_type_fof_self, task_subtype_none, 0, 0, ci,
                        NULL);
    else
      continue;

    /* Now loop over all the neighbours of this cell */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (iii < 0 || iii >= cdim[0]) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (jjj < 0 || jjj >= cdim[1]) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (kkk < 0 || kkk >= cdim[2]) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Is that neighbour local and does it have particles ? */
          if (cid >= cjd || cj->grav.count == 0 || (ci->nodeID != cj->nodeID))
            continue;

          /* Construct the pair task */
          scheduler_addtask(sched, task_type_fof_pair, task_subtype_none, 0, 0,
                            ci, cj);
        }
      }
    }
  }
}

/**
 * @brief Add send tasks for the gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_grav The send_grav #task, if it has already been created.
 */
void engine_addtasks_send_zoom_gravity(struct engine *e, struct cell *ci,
                                       struct cell *cj,
                                       struct task *t_grav, int tag) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Early abort (are we below the level where tasks are)? */
  if (ci->type == zoom && !cell_get_flag(ci, cell_flag_has_tasks)) return;

  /* Create the tasks and their dependencies? */
  if (t_grav == NULL) {

    /* Make sure this cell is tagged. */
    cell_ensure_tagged(ci);

    t_grav = scheduler_addtask(s, task_type_send, task_subtype_gpart_void,
                               ci->mpi.tag, 0, ci, cj);

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.send, t_grav);

    /* We need to propagate the tag down the tree. */
    tag = ci->mpi.tag;

  }

  /* Check if any of the gravity tasks are for the target node. */
  for (l = ci->grav.grav; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    if (t_grav != NULL && ci->subtype == zoom && ci->grav.super == ci) {

      /* The sends should unlock the down pass. */
      scheduler_addunlock(s, t_grav, ci->grav.down);

      /* Drift before you send */
      scheduler_addunlock(s, ci->grav.drift, t_grav);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.send, t_grav);

    /* Assign the tag. */
    ci->mpi.tag = tag;
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_zoom_gravity(e, ci->progeny[k], cj,
                                          t_grav, tag);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The void #cell.
 * @param zoom_c The foreign zoom #cell.
 * @param t_grav The recv_gpart #task, if it has already been created.
 * @param tend The top-level time-step communication #task.
 */
void engine_addtasks_recv_zoom_gravity(struct engine *e, struct cell *c,
                                       struct cell *zoom_c,
                                       struct task *t_grav,
                                       struct task *tend) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;
  const int nodeID = zoom_c->nodeID;
  const int tag = zoom_c->mpi.tag;

  /* Early abort (are we below the level where tasks are)? */
  if (c->type == zoom && !cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Do we need to make a task? */
  if (t_grav == NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    /* Create the tasks. */
    t_grav = scheduler_addtask(s, task_type_recv, task_subtype_gpart_void,
                               tag, 0, c, zoom_c);
    engine_addlink(e, &c->mpi.recv, t_grav);
  }

  /* If we have tasks, link them. */
  if (t_grav != NULL && c->type == zoom && c->nodeID == nodeID) {
    engine_addlink(e, &c->mpi.recv, t_grav);

    /* Get the timestep exchange task if we don't have it yet. */
    if (tend == NULL) {
      for (struct link *ll = c->mpi.recv; ll != NULL; ll = ll->next) {
        if (ll->t->subtype == task_subtype_tend) {
          tend = ll->t;
          break;
        }
      }
    }

    for (struct link *l = c->grav.grav; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_grav, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_zoom_gravity(e, c->progeny[k], zoom_c,
                                          t_grav, tend);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Count the number gparts we need to send to each rank.
 *
 * @param The #cell.
 * @param e The #engine.
 * @param count How many particles to send?
 * @param nodeID Where are we sending them?
 */
int void_count_send_gparts(struct cell *c, struct engine *e, int count,
                            int nodeID) {

  /* Do we need to recurse? */
  if (c->type != zoom) {
    for (int k = 0; k < 8; k++) {
      count += void_count_send_gparts(c->progeny[k], e, count, nodeID);
    }
    return count;
  }

  /* We only want cells on the target node. */
  if (c->nodeID == nodeID) return count;

  /* Is this cell in the proxy? */
  struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
  for (int i = 0; i < p->nr_cells_out; i++) {
    if (p->cells_out[i] == c) {
      count += c->grav.count;
      break;
    }
  }

  return count;
}

/**
 * @brief Count the number gparts we need to recv from each rank.
 *
 * @param The #cell.
 * @param e The #engine.
 * @param counts The counts array to populate for each rank.
 * @param nodeID The ID of the sending node.
 */
int void_count_recv_gparts(struct cell *c, struct engine *e, int count,
                           int nodeID) {

  /* Do we need to recurse? */
  if (c->type != zoom) {
    for (int k = 0; k < 8; k++) {
      count += void_count_recv_gparts(c->progeny[k], e, count, nodeID);
    }
    return count;
  }

  /* Only need cells from the sending node. */
  if (c->nodeID == nodeID) return count;

  /* Is this cell in the proxy? */
  struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
  for (int i = 0; i < p->nr_cells_in; i++) {
    if (p->cells_in[i] == c) {
      count += c->grav.count;
      break;
    }
  }

  return count;
}

/**
 * @brief Count the number gparts we need to send to each rank.
 *
 * @param The #cell.
 * @param e The #engine.
 * @param count The number of particles currently in the buffer.
 * @param buff The buffer contain #gpart structs we need to send.
 * @param nodeID The foreign node we are sending to.
 */
int void_attach_send_gparts(struct cell *c, struct engine *e, int count,
                             struct gpart *buff, int nodeID) {

  /* Do we need to recurse? */
  if (c->type != zoom) {
    for (int k = 0; k < 8; k++) {
      count += void_attach_send_gparts(c->progeny[k], e, count, buff, nodeID);
    }
    return count;
  }

  /* Don't need cells not on the target node. */
  if (c->nodeID != nodeID) return count;

  /* Is this cell in the proxy? */
  struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
  for (int i = 0; i < p->nr_cells_out; i++) {
    if (p->cells_out[i] == c) {
      buff[count] = *c->grav.parts;
      count += c->grav.count;
      break;
    }
  }

  return count;
}

/**
 * @brief Are any of the zoom cell leaves on the target node active?
 *
 * @param The #cell.
 * @param e The #engine.
 * @param is_active The flag for whether the void cell has active.
 * @param nodeID The ID of the sending node.
 * @param send_or_recv which proxy do we search? 0 for send, 1 for recv.
 */
int void_is_active(struct cell *c, struct engine *e, int is_active,
                   int nodeID, int recv) {

  /* Do we need to recurse? */
  if (c->type != zoom) {
    for (int k = 0; k < 8; k++) {
      is_active |= void_is_active(c->progeny[k], e, is_active, nodeID,
                                  recv);
      if (is_active) break;
    }
    return is_active;
  }

  /* Don't need cells not on the target node. */
  if (c->nodeID != nodeID) return is_active;

  /* Is this cell in the proxy? */
  if (recv) {
    struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
    for (int i = 0; i < p->nr_cells_in; i++) {
      if (p->cells_in[i] == c) {
        /* Is this cell active? */
        return cell_is_active_gravity(c, e);
      }
    }
  } else {
    struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
    for (int i = 0; i < p->nr_cells_out; i++) {
      if (p->cells_out[i] == c) {
        /* Is this cell active? */
        return cell_is_active_gravity(c, e);
      }
    }
  }

  return is_active;
}

/**
 * @brief Get one of the zoom cell leaves in the target proxy.
 *
 * @param The #cell.
 * @param e The #engine.
 * @param is_active The flag for whether the void cell has active.
 * @param nodeID The ID of the sending node.
 * @param recv Which proxy do we search? 0 for send, 1 for recv.
 */
void void_get_zoom_on_node(struct cell *c, struct cell *zoom_c,
                           struct engine *e,
                           int nodeID, int recv) {

  /* Do we need to recurse? */
  if (c->type != zoom) {
    for (int k = 0; k < 8; k++) {
      void_get_zoom_on_node(c->progeny[k], zoom_c, e, nodeID, recv);
      if (zoom_c != NULL) break;
    }
    return;
  }

  /* Return the cell if on the target node. */
  if (c->nodeID == nodeID) {

    /* Is this cell in the proxy? */
    if (recv) {
      struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
      for (int i = 0; i < p->nr_cells_in; i++) {
        if (p->cells_in[i] == c) {
          zoom_c = c;
        }
      }
    } else {
      struct proxy *p = &e->proxies[e->proxy_ind[nodeID]];
      for (int i = 0; i < p->nr_cells_out; i++) {
        if (p->cells_out[i] == c) {
          zoom_c = c;
        }
      }
    }
  }
}

/**
 * @brief Construct send tasks for void cells.
 *
 * Each void cell gets a send for each node present in its zoom cell progeny
 * which is also in the proxy.
 *
 * @param e The #engine.
 */
void engine_addtasks_send_void(struct engine *e) {

  /* Get some things we will need. */
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int nr_voids = s->zoom_props->nr_void_cells;
  const int *void_cells = s->zoom_props->void_cells_top;
  const int nodeID = e->nodeID;
  const int nr_nodes = e->nr_nodes;

  /* Loop over ranks. */
  for (int inode = 0; inode < nr_nodes; inode++) {

    /* Skip this rank. */
    if (nodeID == inode) continue;

    /* Loop over void cells. */
    for (int n = 0; n < nr_voids; n++) {

      /* Get the void cell. */
      struct cell *void_c = &cells[void_cells[n]];

      message("Working on void cell in send %d", n);
      /* Get one of the zoom progeny for the target node. */
      struct cell *zoom_c = NULL;
      void_get_zoom_on_node(void_c, zoom_c, e, inode,
                            /*send_or_recv*/0);

      /* If there are no valid zoom progeny: skip. */
      if (zoom_c == NULL) continue;
      message("Got zoom cell in send");
      /* Make the send, link it and add unlocks. */
      engine_addtasks_send_zoom_gravity(e, void_c, zoom_c, /*tgrav*/NULL,
                                        /*tag*/-1);
    }
  }
}

/**
 * @brief Construct recv tasks for void cells.
 *
 * Each void cell gets a recv for each node present in its zoom cell progeny
 * which is also in the proxy.
 *
 * @param e The #engine.
 */
void engine_addtasks_recv_void(struct engine *e) {
#ifdef WITH_MPI
  /* Get some things we will need. */
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int nr_voids = s->zoom_props->nr_void_cells;
  const int *void_cells = s->zoom_props->void_cells_top;
  const int nodeID = e->nodeID;
  const int nr_nodes = e->nr_nodes;

  /* Loop over ranks. */
  for (int inode = 0; inode < nr_nodes; inode++) {

    /* Skip this rank. */
    if (nodeID == inode) continue;

    /* Loop over void cells. */
    for (int n = 0; n < nr_voids; n++) {

      /* Get the void cell. */
      struct cell *void_c = &cells[void_cells[n]];
      message("Working on void cell in recv %d", n);
      /* Get one of the zoom progeny for the target node. */
      struct cell *zoom_c = NULL;
      void_get_zoom_on_node(void_c, zoom_c, e, inode,
                            /*send_or_recv*/1);

      /* If there are no valid zoom progeny: skip. */
      if (zoom_c == NULL) continue;
      message("Got zoom cell in recv");
      /* Make the recv, link it and add unlocks.
       * Note that the tend is extracted at the right level of the heirarchy. */
      engine_addtasks_recv_zoom_gravity(e, void_c, zoom_c, /*tgrav*/NULL,
                                        /*tend*/NULL);
    }
  }
#endif
}

/**
 * @brief Activate the void cell send and receive tasks.
 *
 * @param e The #engine.
 */
void activate_void_tasks(struct engine *e) {
#ifdef WITH_MPI
  /* Get some things we will need. */
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int nr_voids = s->zoom_props->nr_void_cells;
  const int *void_cells = s->zoom_props->void_cells_top;
  const int nodeID = e->nodeID;
  const int nr_nodes = e->nr_nodes;

  /* Loop over ranks. */
  for (int inode = 0; inode < nr_nodes; inode++) {

    /* Skip this rank. */
    if (nodeID == inode) continue;

    /* Loop over void cells. */
    for (int n = 0; n < nr_voids; n++) {

      /* Get the void cell. */
      struct cell *void_c = &cells[void_cells[n]];

      /* Reset the void cell receive counter. */
      void_c->mpi.num_gparts_recvd = 0;

      /* Activate the receive if there is an active zoom cell to receive. */
      if (void_is_active(void_c, e, /*is_active*/0, inode, /*send_or_recv*/1)) {
        scheduler_activate_void_recv(&e->sched, void_c->mpi.recv,
                                     task_subtype_gpart_void,
                                     inode);
      }

      /* Activate the send if there is an active zoom cell to send. */
      if (void_is_active(void_c, e, /*is_active*/0, inode, /*send_or_recv*/0)) {
        scheduler_activate_send(&e->sched, void_c->mpi.send,
                                task_subtype_gpart_void, inode);
      }
    }
  }
#endif
}
#endif /* WITH_ZOOM_REGION */
