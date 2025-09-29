/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will J. Roper (w.roper@sussex.ac.uk)
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
#include <hdf5.h>

/* Local includes. */
#include "common_io.h"
#include "engine.h"
#include "error.h"
#include "space.h"
#include "zoom.h"

/**
 * @brief Write the zoom region metadata to the header of an HDF5 file.
 *
 * All the metadata written out here is in the internal frame, i.e. with
 * the zoom shift already applied to centre the zoom region in the box. The
 * coordinates on the other hand are shifted back to their original position
 * before writing out.
 *
 * @param root_grp The root HDF5 group of the file.
 * @param head_grp The header HDF5 group of the file.
 * @param e The #engine.
 */
void zoom_write_metadata(hid_t root_grp, hid_t head_grp,
                         const struct space *s) {

  /* Extract the zoom properties */
  const struct zoom_region_properties *zp = s->zoom_props;

  /* Write out the flag saying we have run a zoom simulation (or not if
   * the case may be) */
  io_write_attribute_i(head_grp, "ZoomIn", s->with_zoom_region);

  /* If we haven't run zoom we've written everything we need to know */
  if (!s->with_zoom_region) return;

  /* Create a group for the zoom region metadata */
  hid_t h_zoom =
      H5Gcreate(root_grp, "ZoomRegion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_zoom < 0) error("Failed to create ZoomRegion group in HDF5 file.");

  /* Remove the shift from the center of mass */
  double center[3] = {zp->com[0] - zp->zoom_shift[0],
                      zp->com[1] - zp->zoom_shift[1],
                      zp->com[2] - zp->zoom_shift[2]};

  /* Define the internal used centre (the centre of the box) */
  double internal_center[3] = {0.5 * s->dim[0], 0.5 * s->dim[1],
                               0.5 * s->dim[2]};

  /* Write out the rest of the data.*/
  io_write_attribute(h_zoom, "CentreOfMass", DOUBLE, center, 3);
  io_write_attribute(h_zoom, "Shift", DOUBLE, zp->applied_zoom_shift, 3);
  io_write_attribute(h_zoom, "Size", DOUBLE, zp->dim, 3);
  io_write_attribute(h_zoom, "CDim", INT, zp->cdim, 3);
  io_write_attribute_i(h_zoom, "NZoomCells", zp->nr_zoom_cells);
  io_write_attribute(h_zoom, "InternalLowerBounds", DOUBLE,
                     zp->region_lower_bounds, 3);
  io_write_attribute(h_zoom, "InternalUpperBounds", DOUBLE,
                     zp->region_upper_bounds, 3);
  io_write_attribute(h_zoom, "InternalCenter", DOUBLE, internal_center, 3);
  io_write_attribute_i(h_zoom, "Depth", zp->zoom_cell_depth);

  /* Write out the velocity shift dependant on cosmology */
  if (s->e->policy & engine_policy_cosmology) {
    io_write_attribute(h_zoom, "ComovingVelocityShift", FLOAT,
                       zp->applied_zoom_vel_shift, 3);
    io_write_attribute(h_zoom, "ScaleFactorLastShift", DOUBLE,
                       &zp->scale_factor_at_last_shift, 1);
  } else {
    io_write_attribute(h_zoom, "VelocityShift", FLOAT,
                       zp->applied_zoom_vel_shift, 3);
  }

  H5Gclose(h_zoom);
}

/**
 * @brief Undo the shift applied to particles to centre the zoom region.
 *
 * NOTE: In a periodic simulation this must have box_wrap called after the
 * unshift to ensure the particle is back in the box. We don't do this here
 * to remove possible duplication of effort.
 *
 * @param s The #space.
 * @param pos The position to unshift.
 */
void zoom_unshift_pos(const struct space *s, double pos[3]) {

  const struct zoom_region_properties *zp = s->zoom_props;

  /* Unshift the position. */
  pos[0] -= zp->zoom_shift[0];
  pos[1] -= zp->zoom_shift[1];
  pos[2] -= zp->zoom_shift[2];
}
