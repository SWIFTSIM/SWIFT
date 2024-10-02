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
  /* Set the zoom cell depth in the background cells. */
  props->zoom_cell_depth =
      parser_get_param_int(params, "ZoomRegion:zoom_top_level_depth");

  /* Set the target background cdim. */
  int bkg_cdim =
      parser_get_opt_param_int(params, "ZoomRegion:bkg_top_level_cells",
                               space_max_top_level_cells_default);
  for (int i = 0; i < 3; i++) {
    props->bkg_cdim[i] = bkg_cdim;
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
 * @brief Compute cell properties.
 *
 * @param s The space
 * @param ini_max_dim The dim of the zoom region before tesselating the
 * volume.
 */
void zoom_get_cell_props(struct space *s, double ini_max_dim) {

  /* Define the background grid. */
  for (int i = 0; i < 3; i++) {
    s->cdim[i] = s->zoom_props->bkg_cdim[i];
    s->width[i] = s->dim[i] / s->cdim[i];
    s->iwidth[i] = 1.0 / s->width[i];
  }

  /* Find the edges of the background cells that contain the zoom region, i.e.
   * the lower and upper edges of the cells containing the edge of the zoom
   * region centred on the middle of the box. */
  for (int i = 0; i < 3; i++) {
    /* Get the edges of the zoom region. */
    const double region_lower_bound = (s->dim[i] / 2.0) - (ini_max_dim / 2.0);
    const double region_upper_bound = (s->dim[i] / 2.0) + (ini_max_dim / 2.0);

    /* Get the indices for the background cell the edges fall in. */
    const int lower_ind = (int)(region_lower_bound * s->iwidth[i]);
    const int upper_ind = (int)(region_upper_bound * s->iwidth[i]) + 1;

    /* Set the new boundaries. */
    s->zoom_props->region_lower_bounds[i] = lower_ind * s->width[i];
    s->zoom_props->region_upper_bounds[i] = upper_ind * s->width[i];

    /* Set the reigon dim. */
    s->zoom_props->dim[i] = s->zoom_props->region_upper_bounds[i] -
                            s->zoom_props->region_lower_bounds[i];

    /* Set zoom cell width */
    s->zoom_props->width[i] =
        s->width[i] / pow(2, s->zoom_props->zoom_cell_depth);
    s->zoom_props->iwidth[i] = 1.0 / s->zoom_props->width[i];

    /* Set the zoom cdim. */
    s->zoom_props->cdim[i] =
        (int)(s->zoom_props->dim[i] * s->zoom_props->iwidth[i] * 1.0001);
  }
}

/**
 * @brief Report Zoom Region Properties
 *
 * This function prints out a table containing the properties of the
 * zoom region, if it is enabled. The table includes information such as
 * dimensions, center, CDIM, background CDIM, zoom boost factor, minimum zoom
 * cell width, background cell width, and the number of wanderers.
 *
 * @param s The space
 */
void zoom_report_cell_properties(const struct space *s) {

  struct zoom_region_properties *zoom_props = s->zoom_props;

  message("%25s = %f", "Zoom Region Pad Factor", zoom_props->region_pad_factor);
  message("%25s = [%f, %f, %f]", "Zoom Region Shift", zoom_props->zoom_shift[0],
          zoom_props->zoom_shift[1], zoom_props->zoom_shift[2]);
  message("%25s = [%f, %f, %f]", "Zoom Region Dimensions", zoom_props->dim[0],
          zoom_props->dim[1], zoom_props->dim[2]);
  message("%25s = %d", "Zoom Depth in Void Tree", zoom_props->zoom_cell_depth);
  message("%25s = [%f, %f, %f]", "Zoom Region Center",
          zoom_props->region_lower_bounds[0] + (zoom_props->dim[0] / 2.0),
          zoom_props->region_lower_bounds[1] + (zoom_props->dim[1] / 2.0),
          zoom_props->region_lower_bounds[2] + (zoom_props->dim[2] / 2.0));
  message(
      "%25s = [%f-%f, %f-%f, %f-%f]", "Zoom Region Bounds",
      zoom_props->region_lower_bounds[0], zoom_props->region_upper_bounds[0],
      zoom_props->region_lower_bounds[1], zoom_props->region_upper_bounds[1],
      zoom_props->region_lower_bounds[2], zoom_props->region_upper_bounds[2]);
  message("%25s = [%d, %d, %d]", "Zoom Region cdim", zoom_props->cdim[0],
          zoom_props->cdim[1], zoom_props->cdim[2]);
  message("%25s = [%f, %f, %f]", "Zoom Cell Width", zoom_props->width[0],
          zoom_props->width[1], zoom_props->width[2]);
  message("%25s = %d", "Number of Zoom Cells", zoom_props->nr_zoom_cells);
  message("%25s = [%d, %d, %d]", "Background cdim", s->cdim[0], s->cdim[1],
          s->cdim[2]);
  message("%25s = [%f, %f, %f]", "Background Cell Width", s->width[0],
          s->width[1], s->width[2]);
  message("%25s = %d", "Number of Background Cells", zoom_props->nr_bkg_cells);
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

  /* We may now have the gravity properties and engine attached so update the
   * neighbour distance in case the gravity props have changed. */
  if (s->e != NULL) {
    s->zoom_props->neighbour_distance =
        s->e->gravity_properties->r_s *
        s->e->gravity_properties->r_cut_max_ratio;
  }

  /* Compute the extent of the zoom region.
   * NOTE: this calculates the shift necessary to move the zoom region to
   * the centre of the box and stores it in s->zoom_props */
  double ini_dim = zoom_get_region_dim_and_shift(s);

  /* Apply the shift we just calculated to the particles. */
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

  /* Include the requested padding around the high resolution particles. */
  double max_dim = ini_dim * s->zoom_props->region_pad_factor;

  /* Get the cell properties (cdim/widths). This will also set the zoom region
   * boundaries to be coincident on the edges of background cells. */
  zoom_get_cell_props(s, max_dim);

  /* Store what the true boost factor ended up being */
  double input_pad_factor = s->zoom_props->region_pad_factor;
  s->zoom_props->region_pad_factor = s->zoom_props->dim[0] / ini_dim;

  /* Ensure we haven't got a zoom region smaller than the high resolution
   * particle distribution. */
  if (s->zoom_props->dim[0] < ini_dim) {
    error(
        "Found a zoom region smaller than the high resolution particle "
        "distribution! Adjust the cell structure "
        "(ZoomRegion:bkg_top_level_cells and ZoomRegion:zoom_top_level_cells)");
  }

  /* Let's be safe and error if we have drastically changed the size of the
  padding region. */
  if ((s->zoom_props->region_pad_factor / input_pad_factor) >= 2)
    message(
        "WARNING: The pad region has to be 2x larger than requested (%f / %f = "
        "%f). Either increase ZoomRegion:region_pad_factor or increase the "
        "number of background cells.",
        s->zoom_props->region_pad_factor, input_pad_factor,
        s->zoom_props->region_pad_factor / input_pad_factor);

  /* If we didn't get an explicit neighbour cell depth we'll use the zoom
   * depth. */
  s->zoom_props->neighbour_max_tree_depth =
      (s->zoom_props->neighbour_max_tree_depth < 0)
          ? s->zoom_props->zoom_cell_depth
          : s->zoom_props->neighbour_max_tree_depth;

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

  /* Report what we have done */
  if (verbose) {
    zoom_report_cell_properties(s);
  }
}
