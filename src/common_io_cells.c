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
#include "random.h"
#include "timeline.h"
#include "units.h"

/* Standard includes */
#include <float.h>

/**
 * @brief Count the non-inhibted particles in the cell and return the
 * min/max positions
 *
 * @param c The #cell.
 * @param subsample Are we subsampling the output?
 * @param susample_ratio Fraction of particles to write when sub-sampling.
 * @param snap_num The snapshot number (used for the sampling random draws).
 * @param min_pos (return) The min position of all particles to write.
 * @param max_pos (return) The max position of all particles to write.
 */
#define CELL_COUNT_NON_INHIBITED_PARTICLES(TYPE, CELL_TYPE)                   \
  cell_count_non_inhibited_##TYPE(                                            \
      const struct cell* c, const int subsample, const float subsample_ratio, \
      const int snap_num, double min_pos[3], double max_pos[3]) {             \
                                                                              \
    const int total_count = c->CELL_TYPE.count;                               \
    const struct TYPE* parts = c->CELL_TYPE.parts;                            \
    long long count = 0;                                                      \
    min_pos[0] = min_pos[1] = min_pos[2] = DBL_MAX;                           \
    max_pos[0] = max_pos[1] = max_pos[2] = -DBL_MAX;                          \
                                                                              \
    for (int i = 0; i < total_count; ++i) {                                   \
      if ((parts[i].time_bin != time_bin_inhibited) &&                        \
          (parts[i].time_bin != time_bin_not_created)) {                      \
                                                                              \
        /* When subsampling, select particles at random */                    \
        if (subsample) {                                                      \
          const float r = random_unit_interval(                               \
              parts[i].id, snap_num, random_number_snapshot_sampling);        \
          if (r > subsample_ratio) continue;                                  \
        }                                                                     \
                                                                              \
        ++count;                                                              \
                                                                              \
        min_pos[0] = min(parts[i].x[0], min_pos[0]);                          \
        min_pos[1] = min(parts[i].x[1], min_pos[1]);                          \
        min_pos[2] = min(parts[i].x[2], min_pos[2]);                          \
                                                                              \
        max_pos[0] = max(parts[i].x[0], max_pos[0]);                          \
        max_pos[1] = max(parts[i].x[1], max_pos[1]);                          \
        max_pos[2] = max(parts[i].x[2], max_pos[2]);                          \
      }                                                                       \
    }                                                                         \
    return count;                                                             \
  }

static long long CELL_COUNT_NON_INHIBITED_PARTICLES(part, hydro);
static long long CELL_COUNT_NON_INHIBITED_PARTICLES(spart, stars);
static long long CELL_COUNT_NON_INHIBITED_PARTICLES(bpart, black_holes);
static long long CELL_COUNT_NON_INHIBITED_PARTICLES(sink, sinks);

/**
 * @brief Count the non-inhibted g-particles in the cell and return the
 * min/max positions
 *
 * @param c The #cell.
 * @param subsample Are we subsampling the output?
 * @param susample_ratio Fraction of particles to write when sub-sampling.
 * @param snap_num The snapshot number (used for the sampling random draws).
 * @param min_pos (return) The min position of all particles to write.
 * @param max_pos (return) The max position of all particles to write.
 */
#define CELL_COUNT_NON_INHIBITED_GPARTICLES(TYPE, PART_TYPE)                  \
  cell_count_non_inhibited_##TYPE(                                            \
      const struct cell* c, const int subsample, const float subsample_ratio, \
      const int snap_num, double min_pos[3], double max_pos[3]) {             \
                                                                              \
    const int total_count = c->grav.count;                                    \
    const struct gpart* gparts = c->grav.parts;                               \
    long long count = 0;                                                      \
    min_pos[0] = min_pos[1] = min_pos[2] = DBL_MAX;                           \
    max_pos[0] = max_pos[1] = max_pos[2] = -DBL_MAX;                          \
                                                                              \
    for (int i = 0; i < total_count; ++i) {                                   \
      if ((gparts[i].time_bin != time_bin_inhibited) &&                       \
          (gparts[i].time_bin != time_bin_not_created) &&                     \
          (gparts[i].type == PART_TYPE)) {                                    \
                                                                              \
        /* When subsampling, select particles at random */                    \
        if (subsample) {                                                      \
          const float r =                                                     \
              random_unit_interval(gparts[i].id_or_neg_offset, snap_num,      \
                                   random_number_snapshot_sampling);          \
          if (r > subsample_ratio) continue;                                  \
        }                                                                     \
                                                                              \
        ++count;                                                              \
                                                                              \
        min_pos[0] = min(gparts[i].x[0], min_pos[0]);                         \
        min_pos[1] = min(gparts[i].x[1], min_pos[1]);                         \
        min_pos[2] = min(gparts[i].x[2], min_pos[2]);                         \
                                                                              \
        max_pos[0] = max(gparts[i].x[0], max_pos[0]);                         \
        max_pos[1] = max(gparts[i].x[1], max_pos[1]);                         \
        max_pos[2] = max(gparts[i].x[2], max_pos[2]);                         \
      }                                                                       \
    }                                                                         \
    return count;                                                             \
  }

static long long CELL_COUNT_NON_INHIBITED_GPARTICLES(dark_matter,
                                                     swift_type_dark_matter);
static long long CELL_COUNT_NON_INHIBITED_GPARTICLES(
    background_dark_matter, swift_type_dark_matter_background);
static long long CELL_COUNT_NON_INHIBITED_GPARTICLES(neutrinos,
                                                     swift_type_neutrino);

/**
 * @brief Count the number of local non-inhibited particles to write.
 *
 * Takes into account downsampling.
 *
 * @param s The #space.
 * @param subsample Are we subsampling?
 * @param subsample_ratio The fraction of particle to keep when subsampling.
 * @param snap_num The snapshot number to use as random seed.
 */
#define IO_COUNT_PARTICLES_TO_WRITE(NAME, TYPE)                               \
  io_count_##NAME##_to_write(const struct space* s, const int subsample,      \
                             const float subsample_ratio,                     \
                             const int snap_num) {                            \
    long long count = 0;                                                      \
    for (int i = 0; i < s->nr_local_cells; ++i) {                             \
      double dummy1[3], dummy2[3];                                            \
      const struct cell* c = &s->cells_top[s->local_cells_top[i]];            \
      count += cell_count_non_inhibited_##TYPE(c, subsample, subsample_ratio, \
                                               snap_num, dummy1, dummy2);     \
    }                                                                         \
    return count;                                                             \
  }

long long IO_COUNT_PARTICLES_TO_WRITE(gas, part);
long long IO_COUNT_PARTICLES_TO_WRITE(dark_matter, dark_matter);
long long IO_COUNT_PARTICLES_TO_WRITE(background_dark_matter,
                                      background_dark_matter);
long long IO_COUNT_PARTICLES_TO_WRITE(stars, spart);
long long IO_COUNT_PARTICLES_TO_WRITE(sinks, sink);
long long IO_COUNT_PARTICLES_TO_WRITE(black_holes, bpart);
long long IO_COUNT_PARTICLES_TO_WRITE(neutrinos, neutrinos);

#if defined(HAVE_HDF5)

#include <hdf5.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Write a single rank 1 or rank 2 array to a hdf5 group.
 *
 * This creates a simple Nxdim array with a chunk size of 1024xdim.
 * The Fletcher-32 filter is applied to the array.
 *
 * @param h_grp The open hdf5 group.
 * @param n The number of elements in the array.
 * @param dim The dimensionality of each element.
 * @param array The data to write.
 * @param type The type of the data to write.
 * @param name The name of the array.
 * @param array_content The name of the parent group (only used for error
 * messages).
 */
void io_write_array(hid_t h_grp, const int n, const int dim, const void* array,
                    const enum IO_DATA_TYPE type, const char* name,
                    const char* array_content) {

  /* Create memory space */
  const hsize_t shape[2] = {(hsize_t)n, dim};
  hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for %s %s", name, array_content);
  hid_t h_err = H5Sset_extent_simple(h_space, dim > 1 ? 2 : 1, shape, shape);
  if (h_err < 0)
    error("Error while changing shape of %s %s data space.", name,
          array_content);

  /* Dataset type */
  hid_t h_type = H5Tcopy(io_hdf5_type(type));

  const hsize_t chunk[2] = {(1024 > n ? n : 1024), dim};
  hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);
  h_err = H5Pset_chunk(h_prop, dim > 1 ? 2 : 1, chunk);
  if (h_err < 0)
    error("Error while setting chunk shapes of %s %s data space.", name,
          array_content);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting check-sum filter on %s %s data space.", name,
          array_content);

  /* Impose SHUFFLE compression */
  h_err = H5Pset_shuffle(h_prop);
  if (h_err < 0)
    error("Error while setting shuffling options for field '%s'.", name);

  /* Impose GZIP compression */
  h_err = H5Pset_deflate(h_prop, 4);
  if (h_err < 0)
    error("Error while setting compression options for field '%s'.", name);

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

/**
 * @brief Compute and write the top-level cell counts and offsets meta-data.
 *
 * @param h_grp the hdf5 group to write to.
 * @param cdim The number of top-level cells along each axis.
 * @param dim The box size.
 * @param cells_top The top-level cells.
 * @param nr_cells The number of top-level cells.
 * @param distributed Is this a distributed snapshot?
 * @param subsample Are we subsampling the different particle types?
 * @param subsample_fraction The fraction of particles to keep when subsampling.
 * @param snap_num The snapshot number used as subsampling random seed.
 * @param global_counts The total number of particles across all nodes.
 * @param global_offsets The offsets of this node into the global list of
 * particles.
 * @param to_write Whether a given particle type should be written to the cell
 * info.
 * @param numFields The number of fields to write for each particle type.
 * @param internal_units The internal unit system.
 * @param snapshot_units The snapshot unit system.
 */
void io_write_cell_offsets(hid_t h_grp, const int cdim[3], const double dim[3],
                           const struct cell* cells_top, const int nr_cells,
                           const double width[3], const int nodeID,
                           const int distributed,
                           const int subsample[swift_type_count],
                           const float subsample_fraction[swift_type_count],
                           const int snap_num,
                           const long long global_counts[swift_type_count],
                           const long long global_offsets[swift_type_count],
                           const int to_write[swift_type_count],
                           const int num_fields[swift_type_count],
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units) {

#ifdef SWIFT_DEBUG_CHECKS
  if (distributed) {
    if (global_offsets[0] != 0 || global_offsets[1] != 0 ||
        global_offsets[2] != 0 || global_offsets[3] != 0 ||
        global_offsets[4] != 0 || global_offsets[5] != 0 ||
        global_offsets[6] != 0)
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

  /* Temporary memory for the min position of particles in the cells */
  double *min_part_pos = NULL, *min_gpart_pos = NULL,
         *min_gpart_background_pos = NULL, *min_spart_pos = NULL,
         *min_bpart_pos = NULL, *min_sink_pos = NULL, *min_nupart_pos = NULL;
  min_part_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_gpart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_gpart_background_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_spart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_bpart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_sink_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  min_nupart_pos = (double*)calloc(3 * nr_cells, sizeof(double));

  /* Temporary memory for the max position of particles in the cells */
  double *max_part_pos = NULL, *max_gpart_pos = NULL,
         *max_gpart_background_pos = NULL, *max_spart_pos = NULL,
         *max_bpart_pos = NULL, *max_sink_pos = NULL, *max_nupart_pos = NULL;
  max_part_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_gpart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_gpart_background_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_spart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_bpart_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_sink_pos = (double*)calloc(3 * nr_cells, sizeof(double));
  max_nupart_pos = (double*)calloc(3 * nr_cells, sizeof(double));

  /* Count of particles in each cell */
  long long *count_part = NULL, *count_gpart = NULL,
            *count_background_gpart = NULL, *count_spart = NULL,
            *count_bpart = NULL, *count_sink = NULL, *count_nupart = NULL;
  count_part = (long long*)malloc(nr_cells * sizeof(long long));
  count_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_background_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_spart = (long long*)malloc(nr_cells * sizeof(long long));
  count_bpart = (long long*)malloc(nr_cells * sizeof(long long));
  count_sink = (long long*)malloc(nr_cells * sizeof(long long));
  count_nupart = (long long*)malloc(nr_cells * sizeof(long long));

  /* Global offsets of particles in each cell */
  long long *offset_part = NULL, *offset_gpart = NULL,
            *offset_background_gpart = NULL, *offset_spart = NULL,
            *offset_bpart = NULL, *offset_sink = NULL, *offset_nupart = NULL;
  offset_part = (long long*)malloc(nr_cells * sizeof(long long));
  offset_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_background_gpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_spart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_bpart = (long long*)malloc(nr_cells * sizeof(long long));
  offset_sink = (long long*)malloc(nr_cells * sizeof(long long));
  offset_nupart = (long long*)malloc(nr_cells * sizeof(long long));

  /* Offsets of the 0^th element */
  offset_part[0] = 0;
  offset_gpart[0] = 0;
  offset_background_gpart[0] = 0;
  offset_spart[0] = 0;
  offset_bpart[0] = 0;
  offset_sink[0] = 0;
  offset_nupart[0] = 0;

  /* Collect the cell information of *local* cells */
  long long local_offset_part = 0;
  long long local_offset_gpart = 0;
  long long local_offset_background_gpart = 0;
  long long local_offset_spart = 0;
  long long local_offset_bpart = 0;
  long long local_offset_sink = 0;
  long long local_offset_nupart = 0;
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

      /* Count real particles that will be written and collect the min/max
       * positions */
      count_part[i] = cell_count_non_inhibited_part(
          &cells_top[i], subsample[swift_type_gas],
          subsample_fraction[swift_type_gas], snap_num, &min_part_pos[i * 3],
          &max_part_pos[i * 3]);

      count_gpart[i] = cell_count_non_inhibited_dark_matter(
          &cells_top[i], subsample[swift_type_dark_matter],
          subsample_fraction[swift_type_dark_matter], snap_num,
          &min_gpart_pos[i * 3], &max_gpart_pos[i * 3]);

      count_background_gpart[i] =
          cell_count_non_inhibited_background_dark_matter(
              &cells_top[i], subsample[swift_type_dark_matter_background],
              subsample_fraction[swift_type_dark_matter_background], snap_num,
              &min_gpart_background_pos[i * 3],
              &max_gpart_background_pos[i * 3]);

      count_spart[i] = cell_count_non_inhibited_spart(
          &cells_top[i], subsample[swift_type_stars],
          subsample_fraction[swift_type_stars], snap_num, &min_spart_pos[i * 3],
          &max_spart_pos[i * 3]);

      count_bpart[i] = cell_count_non_inhibited_bpart(
          &cells_top[i], subsample[swift_type_black_hole],
          subsample_fraction[swift_type_black_hole], snap_num,
          &min_bpart_pos[i * 3], &max_bpart_pos[i * 3]);

      count_sink[i] = cell_count_non_inhibited_sink(
          &cells_top[i], subsample[swift_type_sink],
          subsample_fraction[swift_type_sink], snap_num, &min_sink_pos[i * 3],
          &max_sink_pos[i * 3]);

      count_nupart[i] = cell_count_non_inhibited_neutrinos(
          &cells_top[i], subsample[swift_type_neutrino],
          subsample_fraction[swift_type_neutrino], snap_num,
          &min_nupart_pos[i * 3], &max_nupart_pos[i * 3]);

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
      offset_nupart[i] =
          local_offset_nupart + global_offsets[swift_type_neutrino];

      local_offset_part += count_part[i];
      local_offset_gpart += count_gpart[i];
      local_offset_background_gpart += count_background_gpart[i];
      local_offset_spart += count_spart[i];
      local_offset_bpart += count_bpart[i];
      local_offset_sink += count_sink[i];
      local_offset_nupart += count_nupart[i];

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
      count_nupart[i] = 0;

      offset_part[i] = 0;
      offset_gpart[i] = 0;
      offset_background_gpart[i] = 0;
      offset_spart[i] = 0;
      offset_bpart[i] = 0;
      offset_sink[i] = 0;
      offset_nupart[i] = 0;
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
  MPI_Allreduce(MPI_IN_PLACE, count_nupart, nr_cells, MPI_LONG_LONG_INT,
                MPI_BOR, MPI_COMM_WORLD);

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
  MPI_Allreduce(MPI_IN_PLACE, offset_nupart, nr_cells, MPI_LONG_LONG_INT,
                MPI_BOR, MPI_COMM_WORLD);

  /* For the centres we use a sum as MPI does not like bit-wise operations
     on floating point numbers */
  MPI_Allreduce(MPI_IN_PLACE, centres, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, min_part_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_gpart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_gpart_background_pos, 3 * nr_cells,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_spart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_bpart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_sink_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, min_nupart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(MPI_IN_PLACE, max_part_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_gpart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_gpart_background_pos, 3 * nr_cells,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_spart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_bpart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_sink_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, max_nupart_pos, 3 * nr_cells, MPI_DOUBLE, MPI_SUM,
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

      /* Convert the particle envelopes */
      for (int i = 0; i < nr_cells; ++i) {
        for (int k = 0; k < 3; ++k) {
          min_part_pos[i * 3 + k] *= factor;
          max_part_pos[i * 3 + k] *= factor;
          min_gpart_pos[i * 3 + k] *= factor;
          max_gpart_pos[i * 3 + k] *= factor;
          min_gpart_background_pos[i * 3 + k] *= factor;
          max_gpart_background_pos[i * 3 + k] *= factor;
          min_sink_pos[i * 3 + k] *= factor;
          max_sink_pos[i * 3 + k] *= factor;
          min_spart_pos[i * 3 + k] *= factor;
          max_spart_pos[i * 3 + k] *= factor;
          min_bpart_pos[i * 3 + k] *= factor;
          max_bpart_pos[i * 3 + k] *= factor;
          min_nupart_pos[i * 3 + k] *= factor;
          max_nupart_pos[i * 3 + k] *= factor;
        }
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
    io_write_array(h_grp, nr_cells, /*dim=*/3, centres, DOUBLE, "Centres",
                   "Cell centres");

    /* Group containing the offsets and counts for each particle type */
    hid_t h_grp_offsets = H5Gcreate(h_grp, "OffsetsInFile", H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_offsets < 0) error("Error while creating offsets sub-group");
    hid_t h_grp_files =
        H5Gcreate(h_grp, "Files", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_files < 0) error("Error while creating filess sub-group");
    hid_t h_grp_counts =
        H5Gcreate(h_grp, "Counts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t h_grp_min_pos =
        H5Gcreate(h_grp, "MinPositions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t h_grp_max_pos =
        H5Gcreate(h_grp, "MaxPositions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (h_grp_counts < 0) error("Error while creating counts sub-group");

    if (to_write[swift_type_gas] > 0 && num_fields[swift_type_gas] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType0",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_part, LONGLONG,
                     "PartType0", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_part, LONGLONG,
                     "PartType0", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_part_pos, DOUBLE,
                     "PartType0", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_part_pos, DOUBLE,
                     "PartType0", "max_pos");
    }

    if (to_write[swift_type_dark_matter] > 0 &&
        num_fields[swift_type_dark_matter] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType1",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_gpart, LONGLONG,
                     "PartType1", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_gpart, LONGLONG,
                     "PartType1", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_gpart_pos, DOUBLE,
                     "PartType1", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_gpart_pos, DOUBLE,
                     "PartType1", "max_pos");
    }

    if (to_write[swift_type_dark_matter_background] > 0 &&
        num_fields[swift_type_dark_matter_background] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType2",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1,
                     offset_background_gpart, LONGLONG, "PartType2", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_background_gpart,
                     LONGLONG, "PartType2", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3,
                     min_gpart_background_pos, DOUBLE, "PartType2", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3,
                     max_gpart_background_pos, DOUBLE, "PartType2", "max_pos");
    }

    if (to_write[swift_type_sink] > 0 && num_fields[swift_type_sink] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType3",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_sink, LONGLONG,
                     "PartType3", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_sink, LONGLONG,
                     "PartType3", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_sink_pos, DOUBLE,
                     "PartType3", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_sink_pos, DOUBLE,
                     "PartType3", "max_pos");
    }

    if (to_write[swift_type_stars] > 0 && num_fields[swift_type_stars] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType4",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_spart, LONGLONG,
                     "PartType4", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_spart, LONGLONG,
                     "PartType4", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_spart_pos, DOUBLE,
                     "PartType4", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_spart_pos, DOUBLE,
                     "PartType4", "max_pos");
    }

    if (to_write[swift_type_black_hole] > 0 &&
        num_fields[swift_type_black_hole] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType5",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_bpart, LONGLONG,
                     "PartType5", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_bpart, LONGLONG,
                     "PartType5", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_bpart_pos, DOUBLE,
                     "PartType5", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_bpart_pos, DOUBLE,
                     "PartType5", "max_pos");
    }

    if (to_write[swift_type_neutrino] > 0 &&
        num_fields[swift_type_neutrino] > 0) {
      io_write_array(h_grp_files, nr_cells, /*dim=*/1, files, INT, "PartType6",
                     "files");
      io_write_array(h_grp_offsets, nr_cells, /*dim=*/1, offset_nupart,
                     LONGLONG, "PartType6", "offsets");
      io_write_array(h_grp_counts, nr_cells, /*dim=*/1, count_nupart, LONGLONG,
                     "PartType6", "counts");
      io_write_array(h_grp_min_pos, nr_cells, /*dim=*/3, min_nupart_pos, DOUBLE,
                     "PartType6", "min_pos");
      io_write_array(h_grp_max_pos, nr_cells, /*dim=*/3, max_nupart_pos, DOUBLE,
                     "PartType6", "max_pos");
    }

    H5Gclose(h_grp_offsets);
    H5Gclose(h_grp_files);
    H5Gclose(h_grp_counts);
    H5Gclose(h_grp_min_pos);
    H5Gclose(h_grp_max_pos);
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
  free(count_nupart);
  free(offset_part);
  free(offset_gpart);
  free(offset_background_gpart);
  free(offset_spart);
  free(offset_bpart);
  free(offset_sink);
  free(offset_nupart);
  free(min_part_pos);
  free(min_gpart_pos);
  free(min_gpart_background_pos);
  free(min_spart_pos);
  free(min_bpart_pos);
  free(min_sink_pos);
  free(min_nupart_pos);
  free(max_part_pos);
  free(max_gpart_pos);
  free(max_gpart_background_pos);
  free(max_spart_pos);
  free(max_bpart_pos);
  free(max_sink_pos);
  free(max_nupart_pos);
}

#endif /* HAVE_HDF5 */
