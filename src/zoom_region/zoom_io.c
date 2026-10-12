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
 * @brief Write particle counts split by zoom and background cells.
 *
 * This function writes the number of particles in zoom cells and outside zoom
 * cells (in the background) for each particle type to the snapshot header. This
 * information can be used with the cell lookup table to read specific
 * particles based on their location in the zoom region.
 *
 * Background particles are not included in the cell look up tree since they
 * should be sufficiently few in number to read without issue.
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

  /* Write out the counts as attributes in the header group. */
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
 * The count uses the same subsampling decision and snapshot number as the
 * particle writer. This ensures the calculated counts match the particles from
 * zoom cells that will actually be placed in the output buffer.
 *
 * @param e The #engine.
 * @param subsample Whether each particle type is subsampled in this snapshot.
 * @param subsample_fraction Fraction of each particle type retained when
 * subsampling.
 * @param local Total local particle counts for each particle type.
 * @param local_in_cells (return) Number of particles in local zoom cells for
 * each particle type.
 */
void zoom_io_count_particles_in_cells(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    long long local_in_cells[swift_type_count]) {

  const struct space *s = e->s;
  const int snap_num = e->snapshot_output_count;

  /* Count each type using the snapshot's subsampling selection. */
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
 * The normal MPI layout concatenates each rank's complete particle array. In a
 * zoom each rank's particle array contains particles from zoom cells first,
 * followed by particles from background cells. We therefore cannot simply
 * concatenate each rank's particle array. Instead, we compute the number of
 * particles from zoom and background cells on each rank and their offsets in
 * the output array.
 *
 * @param e The #engine.
 * @param subsample Whether each particle type is subsampled in this snapshot.
 * @param subsample_fraction Fraction of each particle type retained when
 * subsampling.
 * @param local Number of local particles of each type written to the snapshot.
 * @param total Number of particles of each type across all ranks.
 * @param offset Offset of this rank in the normal rank-concatenated layout.
 * @param comm MPI communicator used by the snapshot writers.
 * @param local_in_cells (return) Number of particles in local zoom cells.
 * @param total_in_cells (return) Number of particles in zoom cells across all
 * ranks.
 * @param offset_in_cells (return) Offset in the output array of this rank's
 * particles from zoom cells.
 * @param offset_outside_cells (optional return) Offset in the output array of
 * this rank's particles from background cells. May be NULL.
 */
void zoom_io_prepare_particle_layout(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    const long long total[swift_type_count],
    const long long offset[swift_type_count], MPI_Comm comm,
    long long local_in_cells[swift_type_count],
    long long total_in_cells[swift_type_count],
    long long offset_in_cells[swift_type_count],
    long long offset_outside_cells[swift_type_count]) {

  zoom_io_count_particles_in_cells(e, subsample, subsample_fraction, local,
                                   local_in_cells);

  /* Compute each rank's offset and the total number from zoom cells. */
  int rank;
  MPI_Comm_rank(comm, &rank);
  bzero(offset_in_cells, swift_type_count * sizeof(long long));
  MPI_Exscan(local_in_cells, offset_in_cells, swift_type_count,
             MPI_LONG_LONG_INT, MPI_SUM, comm);
  if (rank == 0) {
    bzero(offset_in_cells, swift_type_count * sizeof(long long));
  }
  MPI_Allreduce(local_in_cells, total_in_cells, swift_type_count,
                MPI_LONG_LONG_INT, MPI_SUM, comm);

  /* Particles from background cells follow all particles from zoom cells. Their
   * rank offset is the total rank offset minus the zoom-cell rank offset. */
  if (offset_outside_cells != NULL) {
    for (int ptype = 0; ptype < swift_type_count; ++ptype) {
      offset_outside_cells[ptype] = total_in_cells[ptype] +
                                    (rank == 0 ? 0 : offset[ptype]) -
                                    offset_in_cells[ptype];
    }
  }
}
#endif

/**
 * @brief Advance all particle pointers in an #io_props by a given offset.
 *
 * In zoom simulations we need to write particles from zoom and background cells
 * to separate locations in the output array. This function jumps along the
 * particle arrays by the given offset so that the next write starts with the
 * correct particle.
 *
 * @param props I/O properties whose backing pointers are updated in place.
 * @param offset Number of local particles from zoom cells preceding the first
 * particle from a background cell.
 */
void zoom_io_advance_particle_pointers(struct io_props *props, size_t offset) {
  /* Advance whichever backing arrays are used by this field. */
  if (props->field != NULL) {
    props->field += offset * props->partSize;
  }
  if (props->parts != NULL) {
    props->parts += offset;
  }
  if (props->xparts != NULL) {
    props->xparts += offset;
  }
  if (props->gparts != NULL) {
    props->gparts += offset;
  }
  if (props->sparts != NULL) {
    props->sparts += offset;
  }
  if (props->bparts != NULL) {
    props->bparts += offset;
  }
  if (props->sinks != NULL) {
    props->sinks += offset;
  }
}

/**
 * @brief Add virtual mappings for one rank's zoom and background particles.
 *
 * A distributed snapshot stores each rank's zoom particles before its
 * background particles. The virtual snapshot instead stores ALL zoom particles
 * before ALL background particles, so each source dataset requires two
 * mappings.
 *
 * @param h_prop The virtual dataset creation property list.
 * @param h_space The virtual dataset dataspace.
 * @param h_source_space The source dataset dataspace.
 * @param file_name The source file name.
 * @param source_dataset_name The source dataset name.
 * @param dimension Number of values stored per particle.
 * @param particle_count Total number of particles in the source dataset.
 * @param count_in_cells Number of source particles in zoom cells.
 * @param start_in_cells (in/out) Current offset for particles from zoom cells
 * in the virtual dataset, advanced by @p count_in_cells.
 * @param start_outside_cells (in/out) Current offset for particles from
 * background cells in the virtual dataset, advanced by their number.
 */
void zoom_io_map_virtual_particle_regions(
    hid_t h_prop, hid_t h_space, hid_t h_source_space, const char *file_name,
    const char *source_dataset_name, int dimension, hsize_t particle_count,
    hsize_t count_in_cells, hsize_t start_in_cells[2],
    hsize_t start_outside_cells[2]) {

  /* Ensure the two portions form a valid partition of the source dataset. */
  if (count_in_cells > particle_count) {
    error("Invalid zoom particle count for '%s' (%llu > %llu).",
          source_dataset_name, (unsigned long long)count_in_cells,
          (unsigned long long)particle_count);
  }

  /* Describe both source portions and their virtual dataset offsets. */
  const hsize_t region_sizes[2] = {count_in_cells,
                                   particle_count - count_in_cells};
  const hsize_t source_offsets[2] = {0, count_in_cells};
  hsize_t *destination_offsets[2] = {start_in_cells, start_outside_cells};
  hsize_t count[2] = {0, dimension > 1 ? dimension : 0};
  hsize_t source_offset[2] = {0, 0};

  /* Map both source portions to their positions in the virtual dataset. */
  for (int region = 0; region < 2; ++region) {
    count[0] = region_sizes[region];
    source_offset[0] = source_offsets[region];

    if (count[0] > 0) {
      /* Select this region's destination in the virtual dataset. */
      herr_t err =
          H5Sselect_hyperslab(h_space, H5S_SELECT_SET,
                              destination_offsets[region], NULL, count, NULL);
      if (err < 0) {
        error("Error selecting virtual space for '%s'.", source_dataset_name);
      }

      /* Select the matching contiguous region in the rank's source file. */
      err = H5Sselect_hyperslab(h_source_space, H5S_SELECT_SET, source_offset,
                                NULL, count, NULL);
      if (err < 0) {
        error("Error selecting source space for '%s'.", source_dataset_name);
      }

      /* Link the selected source region into the virtual dataset. */
      err = H5Pset_virtual(h_prop, h_space, file_name, source_dataset_name,
                           h_source_space);
      if (err < 0) {
        error("Error mapping virtual dataset '%s'.", source_dataset_name);
      }
    }

    /* Advance this destination offset ready for the next rank. */
    destination_offsets[region][0] += count[0];
  }
}

/**
 * @brief Write local particle data from zoom and background cells to the output
 * file.
 *
 * The local memory buffer contains particles from zoom cells followed by
 * particles from background cells. The output file instead contains particles
 * from zoom cells on all ranks followed by particles from background cells on
 * all ranks. This function selects the corresponding memory and file
 * hyperslabs and writes each portion separately.
 *
 * Both writes are issued even when one region is empty. In that case an empty
 * HDF5 selection is used.
 *
 * @param h_data The HDF5 dataset to write.
 * @param h_memspace The HDF5 memory dataspace.
 * @param h_filespace The HDF5 file dataspace.
 * @param h_type The HDF5 datatype.
 * @param buffer Local particle data, ordered by zoom then background cells.
 * @param rank Rank of the HDF5 dataspace.
 * @param dimension Number of values stored per particle.
 * @param count Total number of local particles.
 * @param count_in_cells Number of local particles in zoom cells.
 * @param offset_in_cells Offset in the output file of the local particles from
 * zoom cells.
 * @param offset_outside_cells Offset in the output file of the local particles
 * from background cells.
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
  if (rank == 1) {
    shape[1] = 0;
  }

  /* Write particles from zoom cells first, then those from background cells. */
  for (int region = 0; region < 2; ++region) {
    shape[0] = region_sizes[region];
    offset[0] = memory_offsets[region];
    if (shape[0] > 0) {
      herr_t err = H5Sselect_hyperslab(h_memspace, H5S_SELECT_SET, offset, NULL,
                                       shape, NULL);
      if (err < 0) {
        error("Error selecting memory space for '%s'.", field_name);
      }
      offset[0] = file_offsets[region];
      err = H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offset, NULL,
                                shape, NULL);
      if (err < 0) {
        error("Error selecting file space for '%s'.", field_name);
      }
    } else {
      H5Sselect_none(h_memspace);
      H5Sselect_none(h_filespace);
    }

    /* Issue both writes, using empty selections for empty regions. */
    const herr_t err =
        H5Dwrite(h_data, h_type, h_memspace, h_filespace, H5P_DEFAULT, buffer);
    if (err < 0) {
      error("Error while writing data array '%s'.", field_name);
    }
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
