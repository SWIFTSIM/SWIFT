/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "common_io.h"

/* Local includes. */
#include "cell.h"
#include "timeline.h"
#include "units.h"

#if defined(HAVE_HDF5)

#include <hdf5.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

static long long cell_count_non_inhibited_gas(const struct cell* c) {
  const int total_count = c->hydro.count;
  struct part* parts = c->hydro.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((parts[i].time_bin != time_bin_inhibited) &&
        (parts[i].time_bin != time_bin_not_created)) {
      ++count;
    }
  }
  return count;
}

static long long cell_count_non_inhibited_dark_matter(const struct cell* c) {
  const int total_count = c->grav.count;
  struct gpart* gparts = c->grav.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((gparts[i].time_bin != time_bin_inhibited) &&
        (gparts[i].time_bin != time_bin_not_created) &&
        (gparts[i].type == swift_type_dark_matter)) {
      ++count;
    }
  }
  return count;
}

static long long cell_count_non_inhibited_background_dark_matter(
    const struct cell* c) {
  const int total_count = c->grav.count;
  struct gpart* gparts = c->grav.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((gparts[i].time_bin != time_bin_inhibited) &&
        (gparts[i].time_bin != time_bin_not_created) &&
        (gparts[i].type == swift_type_dark_matter_background)) {
      ++count;
    }
  }
  return count;
}

static long long cell_count_non_inhibited_stars(const struct cell* c) {
  const int total_count = c->stars.count;
  struct spart* sparts = c->stars.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((sparts[i].time_bin != time_bin_inhibited) &&
        (sparts[i].time_bin != time_bin_not_created)) {
      ++count;
    }
  }
  return count;
}

static long long cell_count_non_inhibited_black_holes(const struct cell* c) {
  const int total_count = c->black_holes.count;
  struct bpart* bparts = c->black_holes.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((bparts[i].time_bin != time_bin_inhibited) &&
        (bparts[i].time_bin != time_bin_not_created)) {
      ++count;
    }
  }
  return count;
}

static long long cell_count_non_inhibited_sinks(const struct cell* c) {
  const int total_count = c->sinks.count;
  struct sink* sinks = c->sinks.parts;
  long long count = 0;
  for (int i = 0; i < total_count; ++i) {
    if ((sinks[i].time_bin != time_bin_inhibited) &&
        (sinks[i].time_bin != time_bin_not_created)) {
      ++count;
    }
  }
  return count;
}

/**
 * @brief Write a single 1D array to a hdf5 group.
 *
 * This creates a simple Nx1 array with a chunk size of 1024x1.
 * The Fletcher-32 filter is applied to the array.
 *
 * @param h_grp The open hdf5 group.
 * @param n The number of elements in the array.
 * @param array The data to write.
 * @param type The type of the data to write.
 * @param name The name of the array.
 * @param array_content The name of the parent group (only used for error
 * messages).
 */
void io_write_array(hid_t h_grp, const int n, const void* array,
                    const enum IO_DATA_TYPE type, const char* name,
                    const char* array_content) {

  /* Create memory space */
  const hsize_t shape[2] = {(hsize_t)n, 1};
  hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for %s %s", name, array_content);
  hid_t h_err = H5Sset_extent_simple(h_space, 1, shape, shape);
  if (h_err < 0)
    error("Error while changing shape of %s %s data space.", name,
          array_content);

  /* Dataset type */
  hid_t h_type = H5Tcopy(io_hdf5_type(type));

  const hsize_t chunk[2] = {(1024 > n ? n : 1024), 1};
  hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);
  h_err = H5Pset_chunk(h_prop, 1, chunk);
  if (h_err < 0)
    error("Error while setting chunk shapes of %s %s data space.", name,
          array_content);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting check-sum filter on %s %s data space.", name,
          array_content);

  /* Write */
  hid_t h_data =
      H5Dcreate(h_grp, name, h_type, h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0)
    error("Error while creating dataspace for %s %s.", name, array_content);
  h_err = H5Dwrite(h_data, io_hdf5_type(type), h_space, H5S_ALL, H5P_DEFAULT,
                   array);
  if (h_err < 0) error("Error while writing %s %s.", name, array_content);
  H5Tclose(h_type);
  H5Dclose(h_data);
  H5Pclose(h_prop);
  H5Sclose(h_space);
}

void io_write_cell_offsets(hid_t h_grp, const int cdim[3], const double dim[3],
                           const struct cell* cells_top, const int nr_cells,
                           const double width[3], const int nodeID,
                           const int distributed,
                           const long long global_counts[swift_type_count],
                           const long long global_offsets[swift_type_count],
                           const int num_fields[swift_type_count],
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units) {

#ifdef SWIFT_DEBUG_CHECKS
  if (distributed) {
    if (global_offsets[0] != 0 || global_offsets[1] != 0 ||
        global_offsets[2] != 0 || global_offsets[3] != 0 ||
        global_offsets[4] != 0 || global_offsets[5] != 0)
      error("Global offset non-zero in the distributed io case!");
  }
#endif

  /* Abort if we don't have any cells yet (i.e. haven't constructed the space)
   */
  if (nr_cells == 0) return;

  double cell_width[3] = {width[0], width[1], width[2]};

  /* Temporary memory for the cell-by-cell information */
  double* centres = NULL;
  centres = (double*)malloc(3 * nr_cells * sizeof(double));

  /* Temporary memory for the cell files ID */
  int* files = NULL;
  files = (int*)malloc(nr_cells * sizeof(int));

  /* Count of particles in each cell */
  long long *count_part = NULL, *count_gpart = NULL,
            *count_background_gpart = NULL, *count_spart = NULL,
            *count_bpart = NULL, *count_sink = NULL;
  count_part = (long long*)malloc(nr_cells * sizeof(long long));
  count_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_background_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_spart = (long long*)malloc(nr_cells * sizeof(long long));
  count_bpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_sink = (long long*)malloc(nr_cells * sizeof(long long));

  /* Global offsets of particles in each cell */
  long long *offset_part = NULL, *offset_gpart = NULL,
            *offset_background_gpart = NULL, *offset_spart = NULL,
            *offset_bpart = NULL, *offset_sink = NULL;
  offset_part = (long long*)malloc(nr_cells * sizeof(long long));
  offset_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_background_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_spart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_bpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_sink = (long long*)malloc(nr_cells * sizeof(long long));

  /* Offsets of the 0^th element */
  offset_part[0] = 0;
  offset_gpart[0] = 0;
  offset_background_gpart[0] = 0;
  offset_spart[0] = 0;
  offset_bpart[0] = 0;
  offset_sink[0] = 0;

  /* Collect the cell information of *local* cells */
  long long local_offset_part = 0;
  long long local_offset_gpart = 0;
  long long local_offset_background_gpart = 0;
  long long local_offset_spart = 0;
  long long local_offset_bpart = 0;
  long long local_offset_sink = 0;
  for (int i = 0; i < nr_cells; ++i) {

    /* Store in which file this cell will be found */
    if (distributed) {
      files[i] = cells_top[i].nodeID;
    } else {
      files[i] = 0;
    }

    /* Is the cell on this node (i.e. we have full information */
    if (cells_top[i].nodeID == nodeID) {

      /* Centre of each cell */
      centres[i * 3 + 0] = cells_top[i].loc[0] + cell_width[0] * 0.5;
      centres[i * 3 + 1] = cells_top[i].loc[1] + cell_width[1] * 0.5;
      centres[i * 3 + 2] = cells_top[i].loc[2] + cell_width[2] * 0.5;

      /* Finish by box wrapping to match what is done to the particles */
      centres[i * 3 + 0] = box_wrap(centres[i * 3 + 0], 0.0, dim[0]);
      centres[i * 3 + 1] = box_wrap(centres[i * 3 + 1], 0.0, dim[1]);
      centres[i * 3 + 2] = box_wrap(centres[i * 3 + 2], 0.0, dim[2]);

      /* Count real particles that will be written */
      count_part[i] = cell_count_non_inhibited_gas(&cells_top[i]);
      count_gpart[i] = cell_count_non_inhibited_dark_matter(&cells_top[i]);
      count_background_gpart[i] =
          cell_count_non_inhibited_background_dark_matter(&cells_top[i]);
      count_spart[i] = cell_count_non_inhibited_stars(&cells_top[i]);
      count_bpart[i] = cell_count_non_inhibited_black_holes(&cells_top[i]);
      count_sink[i] = cell_count_non_inhibited_sinks(&cells_top[i]);

      /* Offsets including the global offset of all particles on this MPI rank
       * Note that in the distributed case, the global offsets are 0 such that
       * we actually compute the offset in the file written by this rank. */
      offset_part[i] = local_offset_part + global_offsets[swift_type_gas];
      offset_gpart[i] =
          local_offset_gpart + global_offsets[swift_type_dark_matter];
      offset_background_gpart[i] =
          local_offset_background_gpart +
          global_offsets[swift_type_dark_matter_background];
      offset_spart[i] = local_offset_spart + global_offsets[swift_type_stars];
      offset_bpart[i] =
          local_offset_bpart + global_offsets[swift_type_black_hole];
      offset_sink[i] = local_offset_sink + global_offsets[swift_type_sink];

      local_offset_part += count_part[i];
      local_offset_gpart += count_gpart[i];
      local_offset_background_gpart += count_background_gpart[i];
      local_offset_spart += count_spart[i];
      local_offset_bpart += count_bpart[i];
      local_offset_sink += count_sink[i];

    } else {

      /* Just zero everything for the foregin cells */

      centres[i * 3 + 0] = 0.;
      centres[i * 3 + 1] = 0.;
      centres[i * 3 + 2] = 0.;

      count_part[i] = 0;
      count_gpart[i] = 0;
      count_background_gpart[i] = 0;
      count_spart[i] = 0;
      count_bpart[i] = 0;
      count_sink[i] = 0;

      offset_part[i] = 0;
      offset_gpart[i] = 0;
      offset_background_gpart[i] = 0;
      offset_spart[i] = 0;
      offset_bpart[i] = 0;
      offset_sink[i] = 0;
    }
  }

#ifdef WITH_MPI
  /* Now, reduce all the arrays. Note that we use a bit-wise OR here. This
     is safe as we made sure only local cells have non-zero values. */
  MPI_Allreduce(MPI_IN_PLACE, count_part, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count_gpart, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count_background_gpart, nr_cells,
                MPI_LONG_LONG_INT, MPI_BOR, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count_sink, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count_spart, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, count_bpart, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, offset_part, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, offset_gpart, nr_cells, MPI_LONG_LONG_INT,
                MPI_BOR, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, offset_background_gpart, nr_cells,
                MPI_LONG_LONG_INT, MPI_BOR, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, offset_sink, nr_cells, MPI_LONG_LONG_INT, MPI_BOR,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, offset_spart, nr_cells, MPI_LONG_LONG_INT,
                MPI_BOR, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, offset_bpart, nr_cells, MPI_LONG_LONG_INT,
                MPI_BOR, MPI_COMM_WORLD);

  /* For the centres we use a sum as MPI does not like bit-wise operations
     on floating point numbers */
  MPI_Allreduce(MPI_IN_PLACE, centres, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  /* When writing a single file, only rank 0 writes the meta-data */
  if ((distributed) || (!distributed && nodeID == 0)) {

    /* Unit conversion if necessary */
    const double factor = units_conversion_factor(
        internal_units, snapshot_units, UNIT_CONV_LENGTH);
    if (factor != 1.) {

      /* Convert the cell centres */
      for (int i = 0; i < nr_cells; ++i) {
        centres[i * 3 + 0] *= factor;
        centres[i * 3 + 1] *= factor;
        centres[i * 3 + 2] *= factor;
      }

      /* Convert the cell widths */
      cell_width[0] *= factor;
      cell_width[1] *= factor;
      cell_width[2] *= factor;
    }

    /* Write some meta-information first */
    hid_t h_subgrp =
        H5Gcreate(h_grp, "Meta-data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_subgrp < 0) error("Error while creating meta-data sub-group");
    io_write_attribute(h_subgrp, "nr_cells", INT, &nr_cells, 1);
    io_write_attribute(h_subgrp, "size", DOUBLE, cell_width, 3);
    io_write_attribute(h_subgrp, "dimension", INT, cdim, 3);
    H5Gclose(h_subgrp);

    /* Write the centres to the group */
    hsize_t shape[2] = {(hsize_t)nr_cells, 3};
    hid_t h_space = H5Screate(H5S_SIMPLE);
    if (h_space < 0) error("Error while creating data space for cell centres");
    hid_t h_err = H5Sset_extent_simple(h_space, 2, shape, shape);
    if (h_err < 0)
      error("Error while changing shape of gas offsets data space.");
    hid_t h_data = H5Dcreate(h_grp, "Centres", io_hdf5_type(DOUBLE), h_space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_data < 0) error("Error while creating dataspace for gas offsets.");
    h_err = H5Dwrite(h_data, io_hdf5_type(DOUBLE), h_space, H5S_ALL,
                     H5P_DEFAULT, centres);
    if (h_err < 0) error("Error while writing centres.");
    H5Dclose(h_data);
    H5Sclose(h_space);

    /* Group containing the offsets and counts for each particle type */
    hid_t h_grp_offsets = H5Gcreate(h_grp, "OffsetsInFile", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_offsets < 0) error("Error while creating offsets sub-group");
    hid_t h_grp_files =
        H5Gcreate(h_grp, "Files", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_files < 0) error("Error while creating filess sub-group");
    hid_t h_grp_counts =
        H5Gcreate(h_grp, "Counts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_counts < 0) error("Error while creating counts sub-group");

    if (global_counts[swift_type_gas] > 0 && num_fields[swift_type_gas] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType0", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_part, LONGLONG,
                     "PartType0", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_part, LONGLONG, "PartType0",
                     "counts");
    }

    if (global_counts[swift_type_dark_matter] > 0 &&
        num_fields[swift_type_dark_matter] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType1", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_gpart, LONGLONG,
                     "PartType1", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_gpart, LONGLONG, "PartType1",
                     "counts");
    }

    if (global_counts[swift_type_dark_matter_background] > 0 &&
        num_fields[swift_type_dark_matter_background] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType2", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_background_gpart, LONGLONG,
                     "PartType2", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_background_gpart, LONGLONG,
                     "PartType2", "counts");
    }

    if (global_counts[swift_type_sink] > 0 && num_fields[swift_type_sink] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType3", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_sink, LONGLONG,
                     "PartType3", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_sink, LONGLONG, "PartType3",
                     "counts");
    }

    if (global_counts[swift_type_stars] > 0 &&
        num_fields[swift_type_stars] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType4", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_spart, LONGLONG,
                     "PartType4", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_spart, LONGLONG, "PartType4",
                     "counts");
    }

    if (global_counts[swift_type_black_hole] > 0 &&
        num_fields[swift_type_black_hole] > 0) {
      io_write_array(h_grp_files, nr_cells, files, INT, "PartType5", "files");
      io_write_array(h_grp_offsets, nr_cells, offset_bpart, LONGLONG,
                     "PartType5", "offsets");
      io_write_array(h_grp_counts, nr_cells, count_bpart, LONGLONG, "PartType5",
                     "counts");
    }

    H5Gclose(h_grp_offsets);
    H5Gclose(h_grp_files);
    H5Gclose(h_grp_counts);
  }

  /* Free everything we allocated */
  free(centres);
  free(files);
  free(count_part);
  free(count_gpart);
  free(count_background_gpart);
  free(count_spart);
  free(count_bpart);
  free(count_sink);
  free(offset_part);
  free(offset_gpart);
  free(offset_background_gpart);
  free(offset_spart);
  free(offset_bpart);
  free(offset_sink);
}

#endif /* HAVE_HDF5 */
