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

/* Config */
#include <config.h>

/* Includes */
#include <float.h>

/* Local includes */
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "space.h"
#include "zoom.h"

/* mpi headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Declare the task diff grav constant. */
int zoom_bkg_subdepth_diff_grav = zoom_bkg_subdepth_diff_grav_default;

/**
 * @brief Read parameter file for "ZoomRegion" properties.
 *
 * @param params Swift parameter structure.
 * @param props The zoom properties struct.
 */
void zoom_parse_params(struct swift_params *params,
                       struct zoom_region_properties *props) {
  /* Set the zoom cdim. */
  props->zoom_cell_depth =
      parser_get_opt_param_int(params, "ZoomRegion:zoom_top_level_depth", 2);

  /* Set the target background cdim, default is a negative value so that if no
   * value is given for a target then the zoom region defines the background
   * cell size. */
  int bkg_cdim =
      parser_get_opt_param_int(params, "ZoomRegion:bkg_top_level_cells",
                               space_max_top_level_cells_default);
  for (int i = 0; i < 3; i++) {
    props->bkg_cdim[i] = bkg_cdim;
  }

  /* Get the ratio between the zoom region size and buffer cell size.
   * Ignored if buffer cells aren't needed. */
  props->buffer_cell_depth =
      parser_get_opt_param_int(params, "ZoomRegion:buffer_top_level_depth", 0);

  /* Ensure the buffer cell depth is less than the zoom cell depth. */
  if (props->buffer_cell_depth > props->zoom_cell_depth) {
    error("Buffer cell depth must be less than the zoom cell depth.");
  }

  /* Extract the zoom width boost factor (used to define the buffer around the
   * zoom region). */
  props->region_pad_factor =
      parser_get_opt_param_float(params, "ZoomRegion:region_pad_factor", 1.1);

  /* Extract the depth we'll split neighbour cells to. */
  props->neighbour_max_tree_depth = parser_get_opt_param_int(
      params, "ZoomRegion:neighbour_max_tree_depth", -1);

  /* Extract the minimum difference between the task level and the leaves
   * for background cells. */
  zoom_bkg_subdepth_diff_grav =
      parser_get_opt_param_int(params, "ZoomRegion:bkg_subdepth_diff_grav",
                               zoom_bkg_subdepth_diff_grav_default);

  /* Should we truncate the background? */
  props->truncate_background =
      parser_get_opt_param_int(params, "ZoomRegion:truncate_background", 0);

  /* If we are truncating the background we need the tidal factor and the
   * target tolerance. */
  if (props->truncate_background) {
    props->tidal_factor = parser_get_opt_param_float(
        params, "ZoomRegion:truncate_tidal_factor", 1.0);
    props->truncate_epsilon =
        parser_get_opt_param_float(params, "ZoomRegion:truncate_epsilon", 1e-3);
  } else {
    props->tidal_factor = 0.0;
    props->truncate_epsilon = 0.0;
  }
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
 */
double zoom_get_region_dim_and_shift(struct space *s) {

  /* Initialise values we will need. */
  const size_t nr_gparts = s->nr_gparts;
  double min_bounds[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  double max_bounds[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  double midpoint[3] = {0.0, 0.0, 0.0};
  double com[3] = {0.0, 0.0, 0.0};
  double mtot = 0.0;
  double ini_dims[3] = {0.0, 0.0, 0.0};
  const double box_mid[3] = {s->dim[0] / 2.0, s->dim[1] / 2.0, s->dim[2] / 2.0};

  /* Find the min/max location in each dimension for each
   * high resolution gravity particle (non-background), and their COM. */
  for (size_t k = 0; k < nr_gparts; k++) {
    /* Skip background particles. */
    if (s->gparts[k].type == swift_type_dark_matter_background) {
      continue;
    }

    /* Unpack the particle positions.
     * NOTE: these will have already been shifted by the user requested amount
     * in space_init if shift in the parameter file is non-zero. */
    const double x = s->gparts[k].x[0];
    const double y = s->gparts[k].x[1];
    const double z = s->gparts[k].x[2];

    /* Wrap if periodic. */
    if (s->periodic) {
      box_wrap(x, 0.0, s->dim[0]);
      box_wrap(y, 0.0, s->dim[1]);
      box_wrap(z, 0.0, s->dim[2]);
    }

    /* Ammend boundaries for this particle. */
    if (x > max_bounds[0]) max_bounds[0] = x;
    if (y > max_bounds[1]) max_bounds[1] = y;
    if (z > max_bounds[2]) max_bounds[2] = z;
    if (x < min_bounds[0]) min_bounds[0] = x;
    if (y < min_bounds[1]) min_bounds[1] = y;
    if (z < min_bounds[2]) min_bounds[2] = z;

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
  for (int i = 0; i < 3; i++) {
    ini_dims[i] = max_bounds[i] - min_bounds[i];
    midpoint[i] = min_bounds[i] + (ini_dims[i] / 2.0);
  }

  /* Calculate the shift needed to place the mid point of the high res
   * particles at the centre of the box. This shift is applied to the
   * particles in space_init in space.c */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->zoom_shift[i] = box_mid[i] - midpoint[i];
  }

  /* We shouldn't shift if the shift is incremental. */
  for (int i = 0; i < 3; i++) {
    if (fabs(s->zoom_props->zoom_shift[i]) < 0.01 * s->dim[i]) {
      s->zoom_props->zoom_shift[i] = 0.0;
    }
  }

  /* If the volume isn't periodic then we can't shift. */
  if (!s->periodic) {
    for (int i = 0; i < 3; i++) {
      if (fabs(s->zoom_props->zoom_shift[i]) > 0.01 * s->dim[i]) {
        error(
            "Cannot shift the zoom region to the centre of the box "
            "when the box is not periodic. Centre the CoM of the high "
            "resolution particles in the box. (shift=[%f, %f, %f], dim=[%f, "
            "%f, %f])",
            s->zoom_props->zoom_shift[0], s->zoom_props->zoom_shift[1],
            s->zoom_props->zoom_shift[2], s->dim[0], s->dim[1], s->dim[2]);
      }
    }
  }

  /* Let's shift the COM.
   * NOTE: boundaries are recalculated relative to box centre later. */
  for (int i = 0; i < 3; i++)
    s->zoom_props->com[i] += s->zoom_props->zoom_shift[i];

  /* Compute maximum side length of the zoom region, we need zoom dim to be
   * equal. */
  double ini_dim = max3(ini_dims[0], ini_dims[1], ini_dims[2]);

  return ini_dim;
}

/**
 * @brief Calculate the distance at which we can truncate the parent volume.
 *
 * This uses a simple geometric argument (based on the (L/R)^3 tidal criterion)
 * to find the distance from the zoom region at which the contributions from the
 * background drop below the desired accuracy.
 *
 * @param zoom_dim The zoom region dimensions.
 * @param tidal_factor The tidal factor accounting for anisotropies in the
 *     background (>1, higher means more background preserved, i.e. more
 *     accurate).
 * @param epsilon The desired accuracy.
 * @return The truncation distance.
 */
static double zoom_compute_bkg_truncate_dist(const double zoom_dim,
                                             const double tidal_factor,
                                             const double epsilon) {

  return tidal_factor * zoom_dim / pow(epsilon, 1.0 / 3.0);
}

/**
 * @brief Truncate the simulation volume to remove distant background.
 *
 * This removes all cells that are further away from the zoom region than the
 * truncation distance computed with zoom_compute_bkg_truncate_dist.
 *
 * @param s The #space.
 * @param verbose Whether to be verbose or not.
 */
void zoom_truncate_background(struct space *s, const double zoom_dim,
                              const int verbose) {

  /* Extract some useful pointers and information. */
  double tidal_factor = s->zoom_props->tidal_factor;
  double epsilon = s->zoom_props->truncate_epsilon;

  /* Compute the truncation distance. */
  const double r_trunc =
      zoom_compute_bkg_truncate_dist(zoom_dim, tidal_factor, epsilon);

  if (verbose)
    message(
        "Computed a truncation distance of %.2f internal units (with %.2f x "
        "%.2f * (%.1e)^(-1/3))",
        r_trunc, tidal_factor, zoom_dim, epsilon);

  /* If the truncation distance exceeds the box size we can't truncate. */
  if (r_trunc * 2.0 >= fmin(s->dim[0], fmin(s->dim[1], s->dim[2]))) {
    error(
        "Truncation distance (%.2e) exceeds box size (%.2e), cannot truncate. "
        "You probably don't need truncation in this case, turn off "
        "ZoomRegion:truncate_background.",
        r_trunc, fmin(s->dim[0], fmin(s->dim[1], s->dim[2])));
    return;
  }

  /* Define the new lower bounds of the box. */
  double new_min[3] = {s->zoom_props->com[0] - r_trunc,
                       s->zoom_props->com[1] - r_trunc,
                       s->zoom_props->com[2] - r_trunc};

  /* Shift all the particles so that they lie within the new truncated box. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    s->parts[k].x[0] -= new_min[0];
    s->parts[k].x[1] -= new_min[1];
    s->parts[k].x[2] -= new_min[2];
  }
  for (size_t k = 0; k < s->nr_gparts; k++) {
    s->gparts[k].x[0] -= new_min[0];
    s->gparts[k].x[1] -= new_min[1];
    s->gparts[k].x[2] -= new_min[2];
  }
  for (size_t k = 0; k < s->nr_sparts; k++) {
    s->sparts[k].x[0] -= new_min[0];
    s->sparts[k].x[1] -= new_min[1];
    s->sparts[k].x[2] -= new_min[2];
  }
  for (size_t k = 0; k < s->nr_bparts; k++) {
    s->bparts[k].x[0] -= new_min[0];
    s->bparts[k].x[1] -= new_min[1];
    s->bparts[k].x[2] -= new_min[2];
  }
  for (size_t k = 0; k < s->nr_sinks; k++) {
    s->sinks[k].x[0] -= new_min[0];
    s->sinks[k].x[1] -= new_min[1];
    s->sinks[k].x[2] -= new_min[2];
  }

  /* Set the new box dimensions. */
  for (int i = 0; i < 3; i++) {
    s->dim[i] = 2.0 * r_trunc;
    s->width[i] = s->dim[i] / s->zoom_props->bkg_cdim[i];
    s->iwidth[i] = 1.0 / s->width[i];
  }

  /* Include the new edge in the zoom shift. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->zoom_shift[i] -= r_trunc;
    s->zoom_props->com[i] = s->dim[i] / 2.0;
  }

  /* Loop over all the gparts and inhibit background particles that are
   * further away than the truncation distance. */
  int ntrunc = 0;
  for (size_t k = 0; k < s->nr_gparts; k++) {

    /* Skip non-background particles. */
    if (s->gparts[k].type != swift_type_dark_matter_background) {
      continue;
    }

    /* Inhibit background particles that are too far away. */
    if (s->gparts[k].x[0] < 0.0 || s->gparts[k].x[0] > s->dim[0] ||
        s->gparts[k].x[1] < 0.0 || s->gparts[k].x[1] > s->dim[1] ||
        s->gparts[k].x[2] < 0.0 || s->gparts[k].x[2] > s->dim[2]) {
      s->gparts[k].time_bin = time_bin_inhibited;
      ntrunc++;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that no particles are outside the new box. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    if (s->parts[k].x[0] < 0.0 || s->parts[k].x[0] > s->dim[0] ||
        s->parts[k].x[1] < 0.0 || s->parts[k].x[1] > s->dim[1] ||
        s->parts[k].x[2] < 0.0 || s->parts[k].x[2] > s->dim[2]) {
      error("Found a part (%zu) outside the truncated box. (%f, %f, %f)", k,
            s->parts[k].x[0], s->parts[k].x[1], s->parts[k].x[2]);
    }
  }
  for (size_t k = 0; k < s->nr_gparts; k++) {
    if (s->gparts[k].time_bin != time_bin_inhibited) {
      continue;
    }
    if (s->gparts[k].x[0] < 0.0 || s->gparts[k].x[0] > s->dim[0] ||
        s->gparts[k].x[1] < 0.0 || s->gparts[k].x[1] > s->dim[1] ||
        s->gparts[k].x[2] < 0.0 || s->gparts[k].x[2] > s->dim[2]) {
      error("Found a gpart (%zu) outside the truncated box. (%f, %f, %f)", k,
            s->gparts[k].x[0], s->gparts[k].x[1], s->gparts[k].x[2]);
    }
  }
  for (size_t k = 0; k < s->nr_sparts; k++) {
    if (s->sparts[k].x[0] < 0.0 || s->sparts[k].x[0] > s->dim[0] ||
        s->sparts[k].x[1] < 0.0 || s->sparts[k].x[1] > s->dim[1] ||
        s->sparts[k].x[2] < 0.0 || s->sparts[k].x[2] > s->dim[2]) {
      error("Found a spart (%zu) outside the truncated box. (%f, %f, %f)", k,
            s->sparts[k].x[0], s->sparts[k].x[1], s->sparts[k].x[2]);
    }
  }
  for (size_t k = 0; k < s->nr_bparts; k++) {
    if (s->bparts[k].x[0] < 0.0 || s->bparts[k].x[0] > s->dim[0] ||
        s->bparts[k].x[1] < 0.0 || s->bparts[k].x[1] > s->dim[1] ||
        s->bparts[k].x[2] < 0.0 || s->bparts[k].x[2] > s->dim[2]) {
      error("Found a bpart (%zu) outside the truncated box. (%f, %f, %f)", k,
            s->bparts[k].x[0], s->bparts[k].x[1], s->bparts[k].x[2]);
    }
  }
  for (size_t k = 0; k < s->nr_sinks; k++) {
    if (s->sinks[k].x[0] < 0.0 || s->sinks[k].x[0] > s->dim[0] ||
        s->sinks[k].x[1] < 0.0 || s->sinks[k].x[1] > s->dim[1] ||
        s->sinks[k].x[2] < 0.0 || s->sinks[k].x[2] > s->dim[2]) {
      error("Found a sink (%zu) outside the truncated box. (%f, %f, %f)", k,
            s->sinks[k].x[0], s->sinks[k].x[1], s->sinks[k].x[2]);
    }
  }
#endif

  if (verbose)
    message("Removing %d background particles out of %zu.", ntrunc,
            s->nr_gparts);
}

/**
 * @brief Compute the void region geometry.
 *
 * The void region is the region covered by background cells above the zoom
 * region. If the void region is sufficiently close to the zoom region size
 * then the two will be made equivalent later on. Otherwide the void region is
 * equivalent to the buffer region.
 *
 * This function will derive the bounds of the void region and return how many
 * zoom regions tesselate the void region.
 *
 * @param s The space
 * @param ini_max_dim The dim of the zoom region before tesselating the
 * volume.
 *
 * @return The number of zoom regions that tesselate the void region.
 */
int zoom_get_void_geometry(struct space *s, const double region_dim) {

  /* Get the lower and upper bounds of the zoom region based on the initial
   * dimensions and the padding factor. */
  double lower_bounds[3];
  double upper_bounds[3];
  for (int i = 0; i < 3; i++) {
    lower_bounds[i] = (s->dim[i] / 2) - (region_dim / 2.0);
    upper_bounds[i] = (s->dim[i] / 2) + (region_dim / 2.0);
  }

  /* Assign these temporary bounds to the zoom region bounds, we'll overwrite
   * these later but useful to have them for now. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->region_lower_bounds[i] = lower_bounds[i];
    s->zoom_props->region_upper_bounds[i] = upper_bounds[i];
  }

  /* Find the background cell edges that contain these bounds. */
  double void_lower_bounds[3];
  double void_upper_bounds[3];
  for (int i = 0; i < 3; i++) {
    int lower = (int)floor(lower_bounds[i] * s->iwidth[i]);
    int upper = (int)floor(upper_bounds[i] * s->iwidth[i]);
    void_lower_bounds[i] = lower * s->width[i];
    void_upper_bounds[i] = (upper + 1) * s->width[i];
  }

  /* Assign the void bounds. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->void_lower_bounds[i] = void_lower_bounds[i];
    s->zoom_props->void_upper_bounds[i] = void_upper_bounds[i];
  }

  /* Compute the void region dimensions. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->void_dim[i] = void_upper_bounds[i] - void_lower_bounds[i];
  }

  /* Compute the number of zoom regions that tesselate the void region. */
  int nr_zoom_regions = ceil(s->zoom_props->void_dim[0] / region_dim);

  return nr_zoom_regions;
}

/**
 * @brief Compute the number of child cells in a single parent cell at a given
 * depth.
 *
 * @param region_dim The dimension of the region.
 * @param parent_width The width of the parent cell.
 * @param child_depth The depth of the child cell within the parent.
 *
 * @return The number of child cells in a single parent cell.
 */
static int zoom_get_cdim_at_depth(double region_dim, double parent_width,
                                  int child_depth) {

  /* How many parent_widths are in the region? (ensure correct rounding) */
  int region_parent_cdim =
      floor((region_dim + (0.1 * parent_width)) / parent_width);

  /* We now know how many parent cells we have in the region, use this and the
   * depth of the zoom region to calculate the cdim (the number of parents
   * times the number of children in a parent. */
  return region_parent_cdim * pow(2, child_depth);
}

void zoom_get_geometry_no_buffer_cells(struct space *s) {

  /* If we have a buffer cell depth warn that we will ignore it. */
  if (s->zoom_props->buffer_cell_depth > 0) {
    message("No buffer cells are needed, ignoring buffer cell depth.");
    s->zoom_props->buffer_cell_depth = 0;
  }

  /* Zero the buffer region properties explictly. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->buffer_lower_bounds[i] = 0.0;
    s->zoom_props->buffer_upper_bounds[i] = 0.0;
    s->zoom_props->buffer_dim[i] = 0.0;
    s->zoom_props->buffer_cdim[i] = 0;
    s->zoom_props->buffer_width[i] = 0.0;
  }

  /* Match the zoom reigon bounds to the void region bounds. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->region_lower_bounds[i] = s->zoom_props->void_lower_bounds[i];
    s->zoom_props->region_upper_bounds[i] = s->zoom_props->void_upper_bounds[i];
  }

  /* Compute the zoom region dimensions. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->dim[i] = s->zoom_props->region_upper_bounds[i] -
                            s->zoom_props->region_lower_bounds[i];
  }

  /* Compute the number of zoom cells in the void region. */
  int cdim = zoom_get_cdim_at_depth(s->zoom_props->dim[0], s->width[0],
                                    s->zoom_props->zoom_cell_depth);

  /* Compute the zoom cdim and cell width. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->cdim[i] = cdim;
    s->zoom_props->width[i] = s->zoom_props->dim[i] / cdim;
    s->zoom_props->iwidth[i] = 1.0 / s->zoom_props->width[i];
  }
}

/**
 * @brief Compute the geometry of the zoom region with buffer cells.
 *
 * This function computes the geometry of the zoom region when buffer cells
 * are enabled. It calculates the bounds, dimensions, and cell widths for both
 * the buffer and zoom regions.
 *
 * Currently, buffer cells are not fully supported and this function will
 * simply through an error if called.
 *
 * @param s The space
 */
void zoom_get_geometry_with_buffer_cells(struct space *s) {

  error(
      "Buffer cells currently provide no performance benefit and carry a "
      "significant complexity cost. They are thus not currently fully "
      "supported. Set ZoomRegion:buffer_top_level_depth to 0 to disable "
      "buffer cells.");

  /* Ensure we have a buffer cell depth. */
  if (s->zoom_props->buffer_cell_depth == 0) {
    error(
        "Current cell structure requires buffer cells but not buffer cell "
        "has "
        "been given. ZoomRegion:buffer_top_level_depth must be greater than "
        "0.");
  }

  /* Match the buffer region bounds to the void region bounds. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->buffer_lower_bounds[i] = s->zoom_props->void_lower_bounds[i];
    s->zoom_props->buffer_upper_bounds[i] = s->zoom_props->void_upper_bounds[i];
  }

  /* Compute the buffer region dimensions. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->buffer_dim[i] = s->zoom_props->buffer_upper_bounds[i] -
                                   s->zoom_props->buffer_lower_bounds[i];
  }

  /* Compute the number of buffer cells in the void region. */
  int buffer_cdim =
      zoom_get_cdim_at_depth(s->zoom_props->buffer_dim[0], s->width[0],
                             s->zoom_props->buffer_cell_depth);

  /* Compute the buffer cdim and cell width. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->buffer_cdim[i] = buffer_cdim;
    s->zoom_props->buffer_width[i] = s->zoom_props->buffer_dim[i] / buffer_cdim;
    s->zoom_props->buffer_iwidth[i] = 1.0 / s->zoom_props->buffer_width[i];
  }

  /* Find the buffer cell edges that contain the zoom region bounds. */
  double region_lower_bounds[3];
  double region_upper_bounds[3];
  for (int i = 0; i < 3; i++) {
    int lower = (int)floor((s->zoom_props->region_lower_bounds[i] -
                            s->zoom_props->buffer_lower_bounds[i]) *
                           s->zoom_props->buffer_iwidth[i]);
    int upper = (int)floor((s->zoom_props->region_upper_bounds[i] -
                            s->zoom_props->buffer_lower_bounds[i]) *
                           s->zoom_props->buffer_iwidth[i]);
    region_lower_bounds[i] = lower * s->zoom_props->buffer_width[i] +
                             s->zoom_props->buffer_lower_bounds[i];
    region_upper_bounds[i] = (upper + 1) * s->zoom_props->buffer_width[i] +
                             s->zoom_props->buffer_lower_bounds[i];
  }

  /* Assign the new aligned zoom bounds. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->region_lower_bounds[i] = region_lower_bounds[i];
    s->zoom_props->region_upper_bounds[i] = region_upper_bounds[i];
  }

  /* Compute the zoom region dimensions. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->dim[i] = region_upper_bounds[i] - region_lower_bounds[i];
  }

  /* Compute the number of zoom cells in the zoom region. Here we need to
   * subtract the buffer depth from the user defined zoom depth, both are
   * defined from the background cells but the calculation below is from the
   * buffer cells down to the zoom level. */
  int cdim = zoom_get_cdim_at_depth(
      s->zoom_props->dim[0], s->zoom_props->buffer_width[0],
      s->zoom_props->zoom_cell_depth - s->zoom_props->buffer_cell_depth);

  /* Compute the zoom cdim and cell width. */
  for (int i = 0; i < 3; i++) {
    s->zoom_props->cdim[i] = cdim;
    s->zoom_props->width[i] = s->zoom_props->dim[i] / cdim;
    s->zoom_props->iwidth[i] = 1.0 / s->zoom_props->width[i];
  }
}

/**
 * @brief Report Zoom Region Properties
 *
 * This function prints out a table containing the properties of the
 * zoom region, if it is enabled. The table includes information such as
 * dimensions, center, CDIM, background CDIM, buffer CDIM, region buffer
 * ratio, zoom boost factor, minimum zoom cell width, background cell width,
 * buffer width, and the number of wanderers.
 *
 * @param s The space
 */
void zoom_report_cell_properties(const struct space *s) {

  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Compute the number of background cells along each side of the void
   * region. */
  for (int i = 0; i < 3; i++) {
    zoom_props->void_cdim[i] =
        (int)ceil(s->zoom_props->void_dim[i] * s->iwidth[i]);
  }

  /* Cdims */
  message("%28s = [%d, %d, %d]", "Background cdim", s->cdim[0], s->cdim[1],
          s->cdim[2]);
  if (zoom_props->with_buffer_cells)
    message("%28s = [%d, %d, %d]", "Buffer cdim", zoom_props->buffer_cdim[0],
            zoom_props->buffer_cdim[1], zoom_props->buffer_cdim[2]);
  message("%28s = [%d, %d, %d]", "Zoom cdim", zoom_props->cdim[0],
          zoom_props->cdim[1], zoom_props->cdim[2]);
  message("%28s = [%d, %d, %d]", "Void cdim", zoom_props->void_cdim[0],
          zoom_props->void_cdim[1], zoom_props->void_cdim[2]);

  /* Dimensions */
  message("%28s = [%f, %f, %f]", "Background Dimensions", s->dim[0], s->dim[1],
          s->dim[2]);
  if (zoom_props->with_buffer_cells)
    message("%28s = [%f, %f, %f]", "Buffer Region Dimensions",
            zoom_props->buffer_dim[0], zoom_props->buffer_dim[1],
            zoom_props->buffer_dim[2]);
  message("%28s = [%f, %f, %f]", "Zoom Region Dimensions", zoom_props->dim[0],
          zoom_props->dim[1], zoom_props->dim[2]);
  message("%28s = [%f, %f, %f]", "Void Region Dimensions",
          s->zoom_props->void_dim[0], s->zoom_props->void_dim[1],
          s->zoom_props->void_dim[2]);

  /* Cell Widths */
  message("%28s = [%f, %f, %f]", "Background Cell Width", s->width[0],
          s->width[1], s->width[2]);
  if (zoom_props->with_buffer_cells)
    message("%28s = [%f, %f, %f]", "Buffer Cell Width",
            zoom_props->buffer_width[0], zoom_props->buffer_width[1],
            zoom_props->buffer_width[2]);
  message("%28s = [%f, %f, %f]", "Zoom Cell Width", zoom_props->width[0],
          zoom_props->width[1], zoom_props->width[2]);

  /* Number of Cells */
  message("%28s = %d", "Number of Background Cells", zoom_props->nr_bkg_cells);
  if (zoom_props->with_buffer_cells)
    message("%28s = %d", "Number of Buffer Cells", zoom_props->nr_buffer_cells);
  message("%28s = %d", "Number of Zoom Cells", zoom_props->nr_zoom_cells);

  /* Bounds */
  if (zoom_props->with_buffer_cells)
    message(
        "%28s = [%f-%f, %f-%f, %f-%f]", "Buffer Bounds",
        zoom_props->buffer_lower_bounds[0], zoom_props->buffer_upper_bounds[0],
        zoom_props->buffer_lower_bounds[1], zoom_props->buffer_upper_bounds[1],
        zoom_props->buffer_lower_bounds[2], zoom_props->buffer_upper_bounds[2]);
  message(
      "%28s = [%f-%f, %f-%f, %f-%f]", "Zoom Region Bounds",
      zoom_props->region_lower_bounds[0], zoom_props->region_upper_bounds[0],
      zoom_props->region_lower_bounds[1], zoom_props->region_upper_bounds[1],
      zoom_props->region_lower_bounds[2], zoom_props->region_upper_bounds[2]);

  /* Depths */
  if (zoom_props->with_buffer_cells)
    message("%28s = %d", "Buffer Top Level Depth",
            zoom_props->buffer_cell_depth);
  message("%28s = %d", "Zoom Top Level Depth", zoom_props->zoom_cell_depth);
  message("%28s = %d", "Neighbour Max Tree Depth",
          zoom_props->neighbour_max_tree_depth);

  /* Assorted extra zoom properties */
  message("%28s = %f", "Zoom Region Pad Factor", zoom_props->region_pad_factor);
  message("%28s = [%f, %f, %f]", "Zoom Region Shift", zoom_props->zoom_shift[0],
          zoom_props->zoom_shift[1], zoom_props->zoom_shift[2]);
  message("%28s = [%f, %f, %f]", "Zoom Region Center",
          zoom_props->region_lower_bounds[0] + (zoom_props->dim[0] / 2.0),
          zoom_props->region_lower_bounds[1] + (zoom_props->dim[1] / 2.0),
          zoom_props->region_lower_bounds[2] + (zoom_props->dim[2] / 2.0));
}

/**
 * @brief Parse and set the zoom region properties.
 *
 * This function allocates the zoom region properties struct and populates it.
 *
 * If we're not running a zoom this function will do nothing.
 *
 * @param params Swift parameter structure.
 * @param s The space
 * @param verbose Are we talking?
 */
void zoom_props_init(struct swift_params *params, struct space *s,
                     const int verbose) {

  /* If not, we're done here */
  if (!s->with_zoom_region) {
    return;
  }

  /* Zoom region properties are stored in a structure. */
  s->zoom_props = (struct zoom_region_properties *)malloc(
      sizeof(struct zoom_region_properties));
  if (s->zoom_props == NULL)
    error("Error allocating memory for the zoom parameters.");
  bzero(s->zoom_props, sizeof(struct zoom_region_properties));

  /* Calculate the gravity mesh distance, we need this for buffer cells and
   * neighbour cell labbeling later on. */
  /* NOTE: when this is first called we don't have the gravity properties (and
   * the engine isn't attached to the space) yet so we need to read directly
   * from the params. */
  /* Get the mesh size */
  int mesh_size = parser_get_param_int(params, "Gravity:mesh_side_length");

  /* Calculate the maximum distance at which we have a gravity task based
   * on the . */
  float a_smooth = parser_get_opt_param_float(params, "Gravity:a_smooth", 1.25);
  float r_cut_max_ratio =
      parser_get_opt_param_float(params, "Gravity:r_cut_max", 4.5);
  float r_s = a_smooth * s->dim[0] / mesh_size;
  s->zoom_props->neighbour_distance = r_s * r_cut_max_ratio;

  /* Parse the parameter file and populate the properties struct. */
  zoom_parse_params(params, s->zoom_props);
}

/**
 * @brief Initialise the zoom region geometry.
 *
 * This will compute the cell grid properties ready for cell
 * cosntruction when zoom_construct_tl_cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void zoom_region_init(struct space *s, const int verbose) {

  /* Nothing to do if we are restarting, just report geometry and move on. */
  if (s->e != NULL && s->e->restarting) {
    if (verbose) zoom_report_cell_properties(s);
    return;
  }

  /* Update the neighbour distance in case the gravity props have changed. */
  if (s->e != NULL) {
    s->zoom_props->neighbour_distance =
        s->e->gravity_properties->r_s *
        s->e->gravity_properties->r_cut_max_ratio;
  }

  /* Compute the extent of the zoom region.
   * NOTE: this calculates the shift necessary to move the zoom region to
   * the centre of the box and stores it in s->zoom_props */
  double ini_dim = zoom_get_region_dim_and_shift(s);

  /* Apply the shift to the particles. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    s->parts[k].x[0] += s->zoom_props->zoom_shift[0];
    s->parts[k].x[1] += s->zoom_props->zoom_shift[1];
    s->parts[k].x[2] += s->zoom_props->zoom_shift[2];
  }
  for (size_t k = 0; k < s->nr_gparts; k++) {
    s->gparts[k].x[0] += s->zoom_props->zoom_shift[0];
    s->gparts[k].x[1] += s->zoom_props->zoom_shift[1];
    s->gparts[k].x[2] += s->zoom_props->zoom_shift[2];
  }
  for (size_t k = 0; k < s->nr_sparts; k++) {
    s->sparts[k].x[0] += s->zoom_props->zoom_shift[0];
    s->sparts[k].x[1] += s->zoom_props->zoom_shift[1];
    s->sparts[k].x[2] += s->zoom_props->zoom_shift[2];
  }
  for (size_t k = 0; k < s->nr_bparts; k++) {
    s->bparts[k].x[0] += s->zoom_props->zoom_shift[0];
    s->bparts[k].x[1] += s->zoom_props->zoom_shift[1];
    s->bparts[k].x[2] += s->zoom_props->zoom_shift[2];
  }
  for (size_t k = 0; k < s->nr_sinks; k++) {
    s->sinks[k].x[0] += s->zoom_props->zoom_shift[0];
    s->sinks[k].x[1] += s->zoom_props->zoom_shift[1];
    s->sinks[k].x[2] += s->zoom_props->zoom_shift[2];
  }

  /* Are we truncating? */
  if (s->zoom_props->truncate_background) {
    zoom_truncate_background(s, ini_dim * s->zoom_props->region_pad_factor,
                             verbose);
  }

  /* Include the requested padding around the high resolution particles. */
  double max_dim = ini_dim * s->zoom_props->region_pad_factor;

  /* Define the background grid (we'll treat this as gospel). */
  for (int i = 0; i < 3; i++) {
    s->cdim[i] = s->zoom_props->bkg_cdim[i];
    s->width[i] = s->dim[i] / s->cdim[i];
    s->iwidth[i] = 1.0 / s->width[i];
  }

  /* Compute the void region bounds and number of zoom regions that tesselate
   * it. */
  int nr_zoom_regions = zoom_get_void_geometry(s, max_dim);

  /* Check the user gave a sensible background cdim, if the number of zoom
   * regions is too high we will have to set up a unworkable number of buffer
   * cells. */
  if (nr_zoom_regions >= 64) {
    error(
        "Background cell size is too large relative to the zoom region! "
        "Increase ZoomRegion:bkg_top_level_cells (would have had %d zoom "
        "regions in the void region).",
        nr_zoom_regions);
  }

  /* If its alot but not silly just warn the user. */
  if (nr_zoom_regions >= 16) {
    warning(
        "Background cell size is large relative to the zoom region! "
        "(would have had %d zoom regions in the void region). ",
        nr_zoom_regions);
  }

  if (verbose) {
    message("Initial geometry gives %d zoom regions in the void region.",
            nr_zoom_regions);
  }

  /* Construct the zoom region geometry. */
  /* NOTE: here we entirely avoid any buffer cell considerations since they
   * provide no performance benefit and are for now vestigual. */
  zoom_get_geometry_no_buffer_cells(s);

  /* Store what the true boost factor ended up being */
  double input_pad_factor = s->zoom_props->region_pad_factor;
  s->zoom_props->region_pad_factor = s->zoom_props->dim[0] / ini_dim;

  /* Ensure we haven't got a zoom region smaller than the high resolution
   * particle distribution. */
  if (s->zoom_props->dim[0] < ini_dim) {
    error(
        "Found a zoom region smaller than the high resolution particle "
        "distribution! Adjust the cell structure "
        "(ZoomRegion:bkg_top_level_cells, ZoomRegion:zoom_top_level_cells)");
  }

  /* Let's be safe and warn if we have drastically changed the size of the
   * requested padding region. */
  if ((s->zoom_props->region_pad_factor / input_pad_factor) >= 2)
    warning(
        "The pad region has to be %d times larger than requested. "
        "Either increase ZoomRegion:region_pad_factor, increase the "
        "number of background cells, or increase the depths of the zoom "
        "cells.",
        (int)(s->zoom_props->region_pad_factor / input_pad_factor));

  /* If we didn't get an explicit neighbour cell depth we'll match the zoom
   * depth. */
  s->zoom_props->neighbour_max_tree_depth =
      (s->zoom_props->neighbour_max_tree_depth < 0)
          ? s->zoom_props->zoom_cell_depth
          : s->zoom_props->neighbour_max_tree_depth;

  /* The neighbour depth must be less the zoom depth or higher if given. */
  if (s->zoom_props->neighbour_max_tree_depth <
      s->zoom_props->zoom_cell_depth) {
    error(
        "Zoom neighbour cell depth (%d) must be greater than or equal to the "
        "zoom cell depth (%d).",
        s->zoom_props->neighbour_max_tree_depth,
        s->zoom_props->zoom_cell_depth);
  }

  /* Set the minimum allowed zoom cell width. */
  const double zoom_dmax =
      max3(s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  s->zoom_props->cell_min = 0.99 * zoom_dmax / s->zoom_props->cdim[0];

  /* Set the minimum background cell size. */
  const double dmax = max3(s->dim[0], s->dim[1], s->dim[2]);
  s->cell_min = 0.99 * dmax / s->cdim[0];

  /* Store cell numbers and offsets. */
  s->zoom_props->bkg_cell_offset =
      s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2];
  s->zoom_props->nr_zoom_cells = s->zoom_props->bkg_cell_offset;
  s->zoom_props->nr_bkg_cells = s->cdim[0] * s->cdim[1] * s->cdim[2];
  s->zoom_props->buffer_cell_offset =
      s->zoom_props->bkg_cell_offset + s->zoom_props->nr_bkg_cells;
  s->zoom_props->nr_buffer_cells = s->zoom_props->buffer_cdim[0] *
                                   s->zoom_props->buffer_cdim[1] *
                                   s->zoom_props->buffer_cdim[2];

  /* Report what we have done */
  if (verbose) {
    zoom_report_cell_properties(s);
  }
}
