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
#include "io_properties.h"
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
  double center[3] = {zp->com[0] - zp->applied_zoom_shift[0],
                      zp->com[1] - zp->applied_zoom_shift[1],
                      zp->com[2] - zp->applied_zoom_shift[2]};

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
  if (zp->zoom_cell_min_cdim > 0)
    io_write_attribute_i(h_zoom, "RequestedMinimumCDim",
                         zp->zoom_cell_min_cdim);

  /* If we are running with cosmology write out the scale factor at the
   * last shift */
  if (s->e->policy & engine_policy_cosmology) {
    io_write_attribute(h_zoom, "ScaleFactorLastShift", DOUBLE,
                       &zp->scale_factor_at_last_shift, 1);
  }

  H5Gclose(h_zoom);
}

/**
 * @brief Write particle counts split between zoom and background cells.
 *
 * @param head_grp The snapshot Header group.
 * @param this_file Number of particles of each type in this file.
 * @param in_cells_this_file Number of particles of each type in zoom cells in
 * this file.
 * @param total Total number of particles of each type across all files.
 * @param in_cells_total Number of particles of each type in zoom cells across
 * all files.
 * @param num_fields Number of fields selected for each particle type.
 */
void zoom_write_particle_counts(
    hid_t head_grp, const long long this_file[swift_type_count],
    const long long in_cells_this_file[swift_type_count],
    const long long total[swift_type_count],
    const long long in_cells_total[swift_type_count],
    const int num_fields[swift_type_count]) {

  /* Derive per-file and global counts, omitting disabled particle types. */
  long long num_in_cells_this_file[swift_type_count] = {0};
  long long num_outside_cells_this_file[swift_type_count] = {0};
  long long num_in_cells_total[swift_type_count] = {0};
  long long num_outside_cells_total[swift_type_count] = {0};
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    if (num_fields[ptype] > 0) {
      num_in_cells_this_file[ptype] = in_cells_this_file[ptype];
      num_outside_cells_this_file[ptype] =
          this_file[ptype] - in_cells_this_file[ptype];
      num_in_cells_total[ptype] = in_cells_total[ptype];
      num_outside_cells_total[ptype] = total[ptype] - in_cells_total[ptype];
    }
  }

  /* Keep the original attributes and provide explicit file and total forms. */
  io_write_attribute(head_grp, "NumParticles_InCells", LONGLONG,
                     num_in_cells_this_file, swift_type_count);
  io_write_attribute(head_grp, "NumParticles_OutsideCells", LONGLONG,
                     num_outside_cells_this_file, swift_type_count);
  io_write_attribute(head_grp, "NumParticles_InCells_ThisFile", LONGLONG,
                     num_in_cells_this_file, swift_type_count);
  io_write_attribute(head_grp, "NumParticles_OutsideCells_ThisFile", LONGLONG,
                     num_outside_cells_this_file, swift_type_count);
  io_write_attribute(head_grp, "NumParticles_InCells_Total", LONGLONG,
                     num_in_cells_total, swift_type_count);
  io_write_attribute(head_grp, "NumParticles_OutsideCells_Total", LONGLONG,
                     num_outside_cells_total, swift_type_count);
}

/**
 * @brief Count local particles in zoom cells for each particle type.
 *
 * @param e The #engine.
 * @param subsample Whether each particle type is being subsampled.
 * @param subsample_fraction Fraction of each particle type to retain.
 * @param local Total local particle counts for each particle type.
 * @param local_in_cells Local counts in zoom cells for each particle type.
 */
void zoom_io_count_particles_in_cells(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    long long local_in_cells[swift_type_count]) {

  /* A non-zoom snapshot treats every local particle as being in-cell. */
  memcpy(local_in_cells, local, swift_type_count * sizeof(long long));
  if (!e->s->with_zoom_region) return;

  /* Count each particle type using the same snapshot subsampling rules. */
  const struct space *s = e->s;
  const int snap_num = e->snapshot_output_count;
  local_in_cells[swift_type_gas] = io_count_gas_in_zoom_to_write(
      s, subsample[swift_type_gas], subsample_fraction[swift_type_gas],
      snap_num);
  local_in_cells[swift_type_dark_matter] =
      io_count_dark_matter_in_zoom_to_write(
          s, subsample[swift_type_dark_matter],
          subsample_fraction[swift_type_dark_matter], snap_num);
  local_in_cells[swift_type_dark_matter_background] =
      io_count_background_dark_matter_in_zoom_to_write(
          s, subsample[swift_type_dark_matter_background],
          subsample_fraction[swift_type_dark_matter_background], snap_num);
  local_in_cells[swift_type_sink] = io_count_sinks_in_zoom_to_write(
      s, subsample[swift_type_sink], subsample_fraction[swift_type_sink],
      snap_num);
  local_in_cells[swift_type_stars] = io_count_stars_in_zoom_to_write(
      s, subsample[swift_type_stars], subsample_fraction[swift_type_stars],
      snap_num);
  local_in_cells[swift_type_black_hole] = io_count_black_holes_in_zoom_to_write(
      s, subsample[swift_type_black_hole],
      subsample_fraction[swift_type_black_hole], snap_num);
  local_in_cells[swift_type_neutrino] = io_count_neutrinos_in_zoom_to_write(
      s, subsample[swift_type_neutrino],
      subsample_fraction[swift_type_neutrino], snap_num);
}

#ifdef WITH_MPI
/**
 * @brief Prepare contiguous zoom and background regions for MPI snapshots.
 *
 * @param e The #engine.
 * @param subsample Whether each particle type is being subsampled.
 * @param subsample_fraction Fraction of each particle type to retain.
 * @param local Total local particle counts for each particle type.
 * @param total Total particle counts across all ranks.
 * @param offset Original rank offset for each particle type.
 * @param comm The MPI communicator.
 * @param layout The particle counts and offsets to populate.
 */
void zoom_io_prepare_particle_layout(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    const long long total[swift_type_count],
    const long long offset[swift_type_count], MPI_Comm comm,
    struct zoom_io_particle_layout *layout) {

  zoom_io_count_particles_in_cells(e, subsample, subsample_fraction, local,
                                   layout->local_in_cells);

  /* Preserve the conventional rank-ordered layout outside zoom runs. */
  if (!e->s->with_zoom_region) {
    memcpy(layout->total_in_cells, total, swift_type_count * sizeof(long long));
    memcpy(layout->offset_in_cells, offset,
           swift_type_count * sizeof(long long));
    if (e->nodeID == 0)
      bzero(layout->offset_in_cells, swift_type_count * sizeof(long long));
    memcpy(layout->offset_outside_cells, total,
           swift_type_count * sizeof(long long));
    return;
  }

  /* Compute each rank's prefix and the global size of the zoom region. */
  bzero(layout->offset_in_cells, swift_type_count * sizeof(long long));
  MPI_Exscan(layout->local_in_cells, layout->offset_in_cells, swift_type_count,
             MPI_LONG_LONG_INT, MPI_SUM, comm);
  if (e->nodeID == 0)
    bzero(layout->offset_in_cells, swift_type_count * sizeof(long long));
  MPI_Allreduce(layout->local_in_cells, layout->total_in_cells,
                swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);

  /* Particles in background cells follow all particles in zoom cells. Their
   * rank prefix is the total rank prefix minus the zoom-cell rank prefix. */
  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    layout->offset_outside_cells[ptype] = layout->total_in_cells[ptype] +
                                          (e->nodeID == 0 ? 0 : offset[ptype]) -
                                          layout->offset_in_cells[ptype];
}
#endif

/**
 * @brief Advance all particle pointers to the background segment.
 *
 * @param props The I/O properties to update.
 * @param offset Number of local particles in zoom cells preceding particles in
 * background cells.
 */
void zoom_io_offset_io_props(struct io_props *props, size_t offset) {
  /* Advance every possible backing array represented by the I/O property. */
  if (props->field != NULL) props->field += offset * props->partSize;
  if (props->parts != NULL) props->parts += offset;
  if (props->xparts != NULL) props->xparts += offset;
  if (props->gparts != NULL) props->gparts += offset;
  if (props->sparts != NULL) props->sparts += offset;
  if (props->bparts != NULL) props->bparts += offset;
  if (props->sinks != NULL) props->sinks += offset;
}

/**
 * @brief Write local particle data from zoom and background cells into their
 * global regions.
 *
 * @param h_data The HDF5 dataset to write.
 * @param h_memspace The HDF5 memory dataspace.
 * @param h_filespace The HDF5 file dataspace.
 * @param h_type The HDF5 datatype.
 * @param buffer The local particle data.
 * @param rank Rank of the HDF5 dataspace.
 * @param dimension Number of values per particle.
 * @param count Total number of local particles.
 * @param count_in_cells Number of local particles in zoom cells.
 * @param offset_in_cells Global offset of this rank's particles in zoom cells.
 * @param offset_outside_cells Global offset of this rank's particles in
 * background cells.
 * @param field_name Name of the particle field being written.
 */
void zoom_io_write_serial_particle_regions(
    hid_t h_data, hid_t h_memspace, hid_t h_filespace, hid_t h_type,
    const void *buffer, int rank, int dimension, size_t count,
    size_t count_in_cells, long long offset_in_cells,
    long long offset_outside_cells, const char *field_name) {

  const hsize_t region_sizes[2] = {count_in_cells, count - count_in_cells};
  const hsize_t memory_offsets[2] = {0, count_in_cells};
  const hsize_t file_offsets[2] = {offset_in_cells, offset_outside_cells};
  hsize_t shape[2] = {0, dimension};
  hsize_t offset[2] = {0, 0};
  if (rank == 1) shape[1] = 0;

  /* Write the zoom segment first, followed by the background segment. */
  for (int region = 0; region < 2; ++region) {
    shape[0] = region_sizes[region];
    offset[0] = memory_offsets[region];

    /* Select matching slabs in memory and in the global file layout. */
    if (shape[0] > 0) {
      herr_t err = H5Sselect_hyperslab(h_memspace, H5S_SELECT_SET, offset, NULL,
                                       shape, NULL);
      if (err < 0) error("Error selecting memory space for '%s'.", field_name);
      offset[0] = file_offsets[region];
      err = H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offset, NULL,
                                shape, NULL);
      if (err < 0) error("Error selecting file space for '%s'.", field_name);
    } else {
      H5Sselect_none(h_memspace);
      H5Sselect_none(h_filespace);
    }

    /* Participate in both collective writes, including empty regions. */
    const herr_t err =
        H5Dwrite(h_data, h_type, h_memspace, h_filespace, H5P_DEFAULT, buffer);
    if (err < 0) error("Error while writing data array '%s'.", field_name);
  }
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
  pos[0] -= zp->applied_zoom_shift[0];
  pos[1] -= zp->applied_zoom_shift[1];
  pos[2] -= zp->applied_zoom_shift[2];
}
