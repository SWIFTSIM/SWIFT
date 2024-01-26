/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
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

/* Includes */
#include <float.h>

/* Config */
#include <config.h>

/* Local includes */
#include "../cell.h"
#include "../engine.h"
#include "../gravity_properties.h"
#include "../space.h"
#include "zoom_init.h"

/* mpi headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Read parameter file for "ZoomRegion" properties.
 *
 * @param params Swift parameter structure.
 * @param props The zoom properties struct.
 */
void zoom_parse_params(struct swift_params *params,
                       struct zoom_region_properties *props) {
#ifdef WITH_ZOOM_REGION
  /* Set the zoom cdim. */
  int zoom_cdim = parser_get_opt_param_int(
      params, "ZoomRegion:zoom_top_level_cells", space_max_top_level_cells);
  for (int ijk = 0; ijk < 3; ijk++) {
    props->cdim[ijk] = zoom_cdim;
  }

  /* Set the target background cdim, default is a negative value so that if no
   * value is given for a target then the zoom region defines the background
   * cell size. */
  int bkg_cdim =
      parser_get_opt_param_int(params, "ZoomRegion:bkg_top_level_cells",
                               space_max_top_level_cells_default);
  for (int ijk = 0; ijk < 3; ijk++) {
    props->bkg_cdim[ijk] = bkg_cdim;
  }

  /* Get the ratio between the zoom region size and buffer cell size.
   * Ignored if buffer cells aren't needed.
   * NOTE: this has to be an integer to ensure the buffer and zoom cells align.
   * The buffer region is divided into region_buffer_ratio cells along each
   * axis. */
  props->region_buffer_ratio = parser_get_opt_param_int(
      params, "ZoomRegion:region_buffer_cell_ratio", 0);

  /* Ensure we have been given a power of 2 times the region buffer ratio
   * for cdim. */
  if (props->region_buffer_ratio) {
    if (!(((props->cdim[0] / props->region_buffer_ratio) &
           ((props->cdim[0] / props->region_buffer_ratio) - 1)) == 0))
      error(
          "Scheduler:max_top_level_cells must be a power "
          "of 2 times region_buffer_cell_ratio (by default 1)"
          "when running with a zoom region!");
  }

  /* Extract the zoom width boost factor (used to define the buffer around the
   * zoom region). */
  props->region_pad_factor =
      parser_get_opt_param_float(params, "ZoomRegion:region_pad_factor", 1.1);

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Compute the zoom region centre and boundaries.
 *
 * Finds the dimensions of the high resolution particle distribution and
 * computes the necessary shift to shift the zoom region to the centre of the
 * box. This shift is stored to be applied in space_init and for
 * transformation when writing out.
 *
 * @param s The space
 * @param params The SWIFT parameter structure.
 */
double zoom_get_region_dim_and_shift(struct space *s,
                                     struct swift_params *params) {

#ifdef WITH_ZOOM_REGION

  /* Initialise values we will need. */
  const size_t nr_gparts = s->nr_gparts;
  double min_bounds[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  double max_bounds[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  double midpoint[3] = {0.0, 0.0, 0.0};
  double com[3] = {0.0, 0.0, 0.0};
  double mtot = 0.0;
  double ini_dims[3] = {0.0, 0.0, 0.0};
  const double box_mid[3] = {s->dim[0] / 2, s->dim[1] / 2, s->dim[2] / 2};

  /* Get the shift from the ICs since this hasn't been applied yet. */
  double shift[3] = {0.0, 0.0, 0.0};
  parser_get_opt_param_double_array(params, "InitialConditions:shift", 3,
                                    shift);

  /* Find the min/max location in each dimension for each
   * high resolution gravity particle (non-background), and their COM. */
  for (size_t k = 0; k < nr_gparts; k++) {
    /* Skip background particles. */
    if (s->gparts[k].type == swift_type_dark_matter_background) {
      continue;
    }

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
    if (x > max_bounds[0]) max_bounds[0] = x;
    if (y > max_bounds[2]) max_bounds[2] = y;
    if (z > max_bounds[4]) max_bounds[4] = z;
    if (x < min_bounds[1]) min_bounds[1] = x;
    if (y < min_bounds[3]) min_bounds[3] = y;
    if (z < min_bounds[5]) min_bounds[5] = z;

    /* Total up mass and position for COM. */
    mtot += s->gparts[k].mass;
    com[0] += x * s->gparts[k].mass;
    com[1] += y * s->gparts[k].mass;
    com[2] += z * s->gparts[k].mass;
  }

#ifdef WITH_MPI
  /* Share answers amoungst nodes. */

  /* Boundary. */
  MPI_Allreduce(MPI_IN_PLACE, &min_bounds[0], 3, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max_bounds[1], 3, MPI_DOUBLE, MPI_MAX,
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

  /* Get the initial dimensions and midpoint. */
  for (int ijk = 0; ijk < 3; ijk++) {
    ini_dims[ijk] = max_bounds[ijk] - min_bounds[ijk];
    midpoint[ijk] = min_bounds[ijk] + (ini_dims[ijk] / 2);
  }

  /* Calculate the shift needed to place the mid point of the high res
   * particles at the centre of the box. This shift is applied to the
   * particles in space_init in space.c */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->zoom_shift[ijk] = box_mid[ijk] - midpoint[ijk];
  }

  /* Let's shift the COM.
   * NOTE: boundaries are recalculated relative to box centre later. */
  for (int ijk = 0; ijk < 3; ijk++)
    s->zoom_props->com[ijk] += s->zoom_props->zoom_shift[ijk];

  /* Compute maximum side length of the zoom region, we need zoom dim to be
   * equal. */
  double ini_dim = max3(ini_dims[0], ini_dims[1], ini_dims[2]);

  return ini_dim;

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Compute the cell properties for the large zoom region case.
 *
 * This case is used when the zoom region is larger than a background cell.
 * This will not respect the ZoomRegion:bkg_top_level_cells parameter but
 * instead treat it as a target.
 *
 * @param s The space
 * @param max_dim The dim of the zoom region to be modified.
 */
void zoom_get_cell_props_large_region(struct space *s, double *max_dim) {

#ifdef WITH_ZOOM_REGION

  /* First we need to define the zoom region width. */
  int nr_zoom_regions = (int)(s->dim[0] / *max_dim);
  (*max_dim) = s->dim[0] / nr_zoom_regions;

  /* Define the requested background cdim as the target. */
  int target_bkg_cdim = s->cdim[0];

  /* Now we can define the background grid. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->cdim[ijk] = (int)floor(s->dim[ijk] / *max_dim);
  }

  /* Compute the new number of a background cells. */
  int new_bkg_cdim = s->cdim[0];
  while (new_bkg_cdim < target_bkg_cdim) {
    new_bkg_cdim *= 2;
  }

  /* Set the background cdim. */
  s->cdim[0] = new_bkg_cdim;
  s->cdim[1] = new_bkg_cdim;
  s->cdim[2] = new_bkg_cdim;

  /* Set the background cell width. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
    s->iwidth[ijk] = 1.0 / s->width[ijk];
  }

  /* Zero the buffer region. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->buffer_lower_bounds[ijk] = 0;
    s->zoom_props->buffer_upper_bounds[ijk] = 0;
    s->zoom_props->buffer_cdim[ijk] = 0;
    s->zoom_props->buffer_width[ijk] = 0;
    s->zoom_props->buffer_iwidth[ijk] = 0;
  }

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Compute the zoom, background and buffer cell grid properties.
 *
 * @param s The space
 * @param grav_props Swift gravity properties.
 * @param max_dim The dim of the zoom region to be modified.
 * @param ini_dim The dim of the zoom region based on the particle distribution.
 * @param verbose Are we talking?
 */
int zoom_get_cell_props_with_buffer_cells(
    struct space *s, const struct gravity_props *grav_props, double *max_dim,
    const double ini_dim) {

#ifdef WITH_ZOOM_REGION

  /* Set the initial zoom_region boundaries with boost factor.
   * The zoom region is already centred on the middle of the box */
  for (int ijk = 0; ijk < 3; ijk++) {
    /* Set the new boundaries. */
    s->zoom_props->region_lower_bounds[ijk] =
        (s->dim[ijk] / 2) - (*max_dim / 2);
    s->zoom_props->region_upper_bounds[ijk] =
        (s->dim[ijk] / 2) + (*max_dim / 2);
  }

  /* Flag that we have buffer cells. */
  s->zoom_props->with_buffer_cells = 1;

  /* Calculate how many background cells we need in the buffer region. The
   * goal is to have this as large as could be necessary, overshooting
   * isn't an issue. */
  const double max_distance = grav_props->r_s * grav_props->r_cut_max_ratio;

  /* Find the buffer region boundaries. The zoom region is already centred on
   * the middle of the box. */
  for (int ijk = 0; ijk < 3; ijk++) {

    /* Find the background cell containing lower and upper bounds of the zoom
     * region's "gravity reach". */
    int lower = (s->zoom_props->region_lower_bounds[ijk] - max_distance) *
                s->iwidth[ijk];
    int upper = (s->zoom_props->region_upper_bounds[ijk] + max_distance) *
                s->iwidth[ijk];

    s->zoom_props->buffer_lower_bounds[ijk] = lower * s->width[ijk];
    s->zoom_props->buffer_upper_bounds[ijk] = (upper + 1) * s->width[ijk];
  }

  /* Define the extent of the buffer region. */
  double buffer_dim = s->zoom_props->buffer_upper_bounds[0] -
                      s->zoom_props->buffer_lower_bounds[0];

  /* Calculate the initial buffer region cdim accounting for how many buffer
   * cells we want in the zoom region. */
  int ini_buffer_cdim =
      (int)(floor(buffer_dim / *max_dim)) * s->zoom_props->region_buffer_ratio;

  /* Calculate the intial width of a buffer cell. */
  double ini_buffer_width = buffer_dim / ini_buffer_cdim;

  /* Now redefine the bounds of the zoom region based on the edges of the
   * buffer cells containing it. */
  for (int ijk = 0; ijk < 3; ijk++) {

    /* Find the background cell containing lower and upper bounds of the zoom
     * regions "gravity reach". */
    int lower = (s->zoom_props->region_lower_bounds[ijk] -
                 s->zoom_props->buffer_lower_bounds[ijk]) /
                ini_buffer_width;
    int upper = (s->zoom_props->region_upper_bounds[ijk] -
                 s->zoom_props->buffer_lower_bounds[ijk]) /
                ini_buffer_width;

    s->zoom_props->region_lower_bounds[ijk] = lower * ini_buffer_width;
    s->zoom_props->region_upper_bounds[ijk] = (upper + 1) * ini_buffer_width;
  }

  /* Calculate the new zoom region dimension. */
  *max_dim = s->zoom_props->region_upper_bounds[0] -
             s->zoom_props->region_lower_bounds[0];

  /* Set the buffer cells properties. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->buffer_cdim[ijk] = ini_buffer_cdim;
    s->zoom_props->buffer_width[ijk] = ini_buffer_width;
    s->zoom_props->buffer_iwidth[ijk] = 1.0 / s->zoom_props->buffer_width[ijk];
  }

  return ((*max_dim) < ini_dim);

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Compute cell properties without buffer cells.
 *
 * This will tesselate background cells the size of the zoom region across
 * the box.
 *
 * @param s The space
 * @param max_dim The dim of the zoom region to be modified.
 */
void zoom_get_cell_props_no_buffer_cells(struct space *s, double *max_dim) {

#ifdef WITH_ZOOM_REGION

  /* Ensure an odd integer number of the zoom regions tessalate the box. */
  int nr_zoom_regions = (int)(s->dim[0] / *max_dim);
  if (nr_zoom_regions % 2 == 0) nr_zoom_regions -= 1;
  *max_dim = s->dim[0] / nr_zoom_regions;

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
    s->zoom_props->buffer_lower_bounds[ijk] = 0;
    s->zoom_props->buffer_upper_bounds[ijk] = 0;
    s->zoom_props->buffer_cdim[ijk] = 0;
    s->zoom_props->buffer_width[ijk] = 0;
    s->zoom_props->buffer_iwidth[ijk] = 0;
  }

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Report Zoom Region Properties
 *
 * This function prints out a table containing the properties of the
 * zoom region, if it is enabled. The table includes information such as
 * dimensions, center, CDIM, background CDIM, buffer CDIM, region buffer ratio,
 * zoom boost factor, minimum zoom cell width, background cell width, buffer
 * width, and the number of wanderers.
 *
 * @param s The space
 */
void zoom_report_cell_properties(const struct space *s) {

#ifdef WITH_ZOOM_REGION

  struct zoom_region_properties *zoom_props = s->zoom_props;

  message("%25s = %f", "Zoom Region Pad Factor", zoom_props->region_pad_factor);
  message("%25s = [%f, %f, %f]", "Zoom Region Shift", zoom_props->zoom_shift[0],
          zoom_props->zoom_shift[1], zoom_props->zoom_shift[2]);
  message("%25s = [%f, %f, %f]", "Zoom Region Dimensions", zoom_props->dim[0],
          zoom_props->dim[1], zoom_props->dim[2]);
  message("%25s = [%f, %f, %f]", "Zoom Region Center",
          zoom_props->region_lower_bounds[0] + (zoom_props->dim[0] / 2),
          zoom_props->region_lower_bounds[1] + (zoom_props->dim[1] / 2),
          zoom_props->region_lower_bounds[2] + (zoom_props->dim[2] / 2));
  message(
      "%25s = [%f-%f, %f-%f, %f-%f]", "Zoom Region Bounds",
      zoom_props->region_lower_bounds[0], zoom_props->region_upper_bounds[0],
      zoom_props->region_lower_bounds[1], zoom_props->region_upper_bounds[1],
      zoom_props->region_lower_bounds[2], zoom_props->region_upper_bounds[2]);
  message("%25s = [%d, %d, %d]", "Zoom Region cdim", zoom_props->cdim[0],
          zoom_props->cdim[1], zoom_props->cdim[2]);
  message("%25s = [%f, %f, %f]", "Zoom Cell Width", zoom_props->width[0],
          zoom_props->width[1], zoom_props->width[2]);
  message("%25s = [%d, %d, %d]", "Background cdim", s->cdim[0], s->cdim[1],
          s->cdim[2]);
  message("%25s = [%f, %f, %f]", "Background Cell Width", s->width[0],
          s->width[1], s->width[2]);
  if (zoom_props->with_buffer_cells) {
    message("%25s = %d", "Region Buffer Ratio",
            zoom_props->region_buffer_ratio);
    message("%25s = [%d, %d, %d]", "Buffer cdim", zoom_props->buffer_cdim[0],
            zoom_props->buffer_cdim[1], zoom_props->buffer_cdim[2]);
    message("%25s = [%f, %f, %f]", "Buffer Width", zoom_props->buffer_width[0],
            zoom_props->buffer_width[1], zoom_props->buffer_width[2]);
    message("%25s = [%f, %f, %f]", "Buffer Dimensions",
            zoom_props->buffer_width[0] * zoom_props->buffer_cdim[0],
            zoom_props->buffer_width[1] * zoom_props->buffer_cdim[1],
            zoom_props->buffer_width[2] * zoom_props->buffer_cdim[2]);
  }

#endif /* WITH_ZOOM_REGION */
}

/**
 * @brief Initialise the zoom region.
 *
 * This will compute the cell grid properties ready for cell
 * cosntruction when space_regrid is called.
 *
 * @param params Swift parameter structure.
 * @param s The space
 */
void zoom_region_init(struct swift_params *params, struct space *s,
                      const struct gravity_props *grav_props,
                      const int verbose) {

#ifdef WITH_ZOOM_REGION
  /* Are we running with a zoom region? */
  s->with_zoom_region =
      parser_get_opt_param_int(params, "ZoomRegion:enable", 0);

  /* If not, we're done here */
  if (!s->with_zoom_region) {
    return;
  }

  /* Zoom region properties are stored in a structure. */
  s->zoom_props = (struct zoom_region_properties *)malloc(
      sizeof(struct zoom_region_properties));
  bzero(s->zoom_props, sizeof(struct zoom_region_properties));
  if (s->zoom_props == NULL)
    error("Error allocating memory for the zoom parameters.");

  /* Parse the parameter file and populate the properties struct. */
  zoom_parse_params(params, s->zoom_props);

  /* Compute the extent of the zoom region.
   * NOTE: this calculates the shift necessary to move the zoom region to
   * the centre of the box and stores it in s->zoom_props */
  double ini_dim = zoom_get_region_dim_and_shift(s, params);

  /* Include the requested padding around the high resolution particles. */
  double max_dim = ini_dim * s->zoom_props->region_pad_factor;

  /* Define the background grid.
   * NOTE: This can be updated below if max_dim > s->width[0]. In that event the
   * number of background cells is modified until an acceptable number is found.
   * See the note below. It can also be modified if max_dim < s->width[0] and
   * no buffer cells are being used. In that event the background cdim becomes
   * the number of zoom regions that tesselate the full box. */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->cdim[ijk] = s->zoom_props->bkg_cdim[ijk];
    s->width[ijk] = s->dim[ijk] / s->cdim[ijk];
    s->iwidth[ijk] = 1.0 / s->width[ijk];
  }

  /* Warn the user if they have turned off buffer cells with a small zoom
   * region. */
  if (max_dim < s->width[0] / 2 && s->zoom_props->region_buffer_ratio == 0) {
    error(
        "Running with a zoom region significantly smaller than a "
        "background cell (region_dim=%f, bkg_cell_width=%f) and no buffer "
        "cells, performance will be poor! Increase "
        "ZoomRegion:region_buffer_cell_ratio",
        max_dim, s->width[0]);
  }

  /* If we have a region larger than a background cell construct the zoom
   * region for that case regardless of buffer cell definition in the
   * parameter file. */
  if (max_dim > s->width[0]) {

    /* NOTE: for this case the number of background cells is defined by
     * the geometry but attempts to get as close as possible to the user
     * defined cdim from the parameter file. */
    zoom_get_cell_props_large_region(s, &max_dim);
  }

  /* If we have buffer cells: use them alongside the zoom and background
   * cells. */
  else if (s->zoom_props->region_buffer_ratio > 0) {

    /* Compute the cell grid properties. */
    if (zoom_get_cell_props_with_buffer_cells(s, grav_props, &max_dim, ini_dim))
      error(
          "Found a zoom region smaller than the high resolution particle "
          "distribution! Adjust the cell structure "
          "(ZoomRegion:bkg_top_level_cells, ZoomRegion:zoom_top_level_cells"
          " and ZoomRegion:region_buffer_cell_ratio)");

  }

  /* Otherwise we simply tessalate cells the size of the zoom region across
   * the whole volume without padding with buffer cells. */
  else {
    zoom_get_cell_props_no_buffer_cells(s, &max_dim);
  }

  /* Finally define the region boundaries in the centre of the box. */
  for (int ijk = 0; ijk < 3; ijk++) {
    /* Set the new boundaries. */
    s->zoom_props->region_lower_bounds[ijk] = (s->dim[ijk] / 2) - (max_dim / 2);
    s->zoom_props->region_upper_bounds[ijk] = (s->dim[ijk] / 2) + (max_dim / 2);

    /* Set the reigon dim. */
    s->zoom_props->dim[ijk] = max_dim;
  }

  /* Store what the true boost factor ended up being */
  s->zoom_props->region_pad_factor = max_dim / ini_dim;

  /* Let's be safe and error if we have drastically changed the size of the
  padding region. */
  if ((max_dim / ini_dim) >= 2)
    error(
        "WARNING: The pad region has to be 2x larger than requested."
        "Either increase ZoomRegion:region_pad_factor or increase the "
        "number of background cells.");

  /* Set zoom cell width */
  for (int ijk = 0; ijk < 3; ijk++) {
    s->zoom_props->width[ijk] =
        s->zoom_props->dim[ijk] / s->zoom_props->cdim[ijk];
  }

  /* Set the minimum allowed zoom cell width. */
  const double zoom_dmax =
      max3(s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  s->zoom_props->cell_min = 0.99 * zoom_dmax / s->zoom_props->cdim[0];

  /* Report what we have done */
  if (verbose) {
    zoom_report_cell_properties(s);
  }

  /* Make sure we have a compatible mesh size. */
  int min_mesh_size = max3((int)(s->dim[0] / s->zoom_props->cell_min),
                           (int)(s->dim[1] / s->zoom_props->cell_min),
                           (int)(s->dim[2] / s->zoom_props->cell_min));
  if (grav_props->mesh_size < min_mesh_size)
    error(
        "Mesh too small given the number of top-level cells. Should be at "
        "least %d cells wide.",
        min_mesh_size);

#endif
}
