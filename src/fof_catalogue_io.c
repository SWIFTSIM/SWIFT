/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifdef HAVE_HDF5

/* Some standard headers. */
#include <libgen.h>

/* Local headers */
#include "engine.h"
#include "fof.h"
#include "hydro_io.h"
#include "tools.h"
#include "version.h"

void write_fof_hdf5_header(hid_t h_file, const struct engine* e,
                           const long long num_groups_total,
                           const long long num_groups_this_file,
                           const struct fof_props* props,
                           const int virtual_file) {

  /* Open header to write simulation properties */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_TIME);
  const double factor_length = units_conversion_factor(
      e->internal_units, e->snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute_d(h_grp, "Time", dblTime);
  io_write_attribute_d(h_grp, "Dimension", (int)hydro_dimension);
  io_write_attribute_d(h_grp, "Redshift", e->cosmology->z);
  io_write_attribute_d(h_grp, "Scale-factor", e->cosmology->a);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* We write rank 0's hostname so that it is uniform across all files. */
  char systemname[256] = {0};
  if (e->nodeID == 0) sprintf(systemname, "%s", hostname());
#ifdef WITH_MPI
  if (!virtual_file) MPI_Bcast(systemname, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
  io_write_attribute_s(h_grp, "System", systemname);
  io_write_attribute(h_grp, "Shift", DOUBLE, e->s->initial_shift, 3);

  /* Write out the particle types */
  io_write_part_type_names(h_grp);

  /* Write out the time-base */
  if (e->policy & engine_policy_cosmology) {
    io_write_attribute_d(h_grp, "TimeBase_dloga", e->time_base);
    const double delta_t = cosmology_get_timebase(e->cosmology, e->ti_current);
    io_write_attribute_d(h_grp, "TimeBase_dt", delta_t);
  } else {
    io_write_attribute_d(h_grp, "TimeBase_dloga", 0);
    io_write_attribute_d(h_grp, "TimeBase_dt", e->time_base);
  }

  /* Store the time at which the snapshot was written */
  time_t tm = time(NULL);
  struct tm* timeinfo = localtime(&tm);
  char snapshot_date[64];
  strftime(snapshot_date, 64, "%T %F %Z", timeinfo);
  io_write_attribute_s(h_grp, "SnapshotDate", snapshot_date);

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
  long long N_total[swift_type_count] = {0};
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
  }
  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total", UINT, numParticles,
                     swift_type_count);
  io_write_attribute(h_grp, "NumPart_Total_HighWord", UINT,
                     numParticlesHighWord, swift_type_count);
  io_write_attribute(h_grp, "TotalNumberOfParticles", LONGLONG, N_total,
                     swift_type_count);
  double MassTable[swift_type_count] = {0};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  io_write_attribute(h_grp, "InitialMassTable", DOUBLE,
                     e->s->initial_mean_mass_particles, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute_i(h_grp, "NumFilesPerSnapshot", e->nr_nodes);
  io_write_attribute_i(h_grp, "ThisFile", e->nodeID);
  io_write_attribute_s(h_grp, "SelectOutput", "Default");
  io_write_attribute_i(h_grp, "Virtual", virtual_file);
  const int to_write[swift_type_count] = {0};
  io_write_attribute(h_grp, "CanHaveTypes", INT, to_write, swift_type_count);
  io_write_attribute_s(h_grp, "OutputType", "FOF");

  /* FOF-specific counters */
  io_write_attribute_ll(h_grp, "NumGroups_Total", num_groups_total);
  io_write_attribute_ll(h_grp, "NumGroups_ThisFile", num_groups_this_file);

  /* Close group */
  H5Gclose(h_grp);

  /* Copy metadata from ICs to the file */
  ic_info_write_hdf5(e->ics_metadata, h_file);

  /* Write all the meta-data */
  io_write_meta_data(h_file, e, e->internal_units, e->snapshot_units,
                     /*fof=*/1);
}

void write_virtual_fof_hdf5_array(
    const struct engine* e, hid_t grp, const char* fileName_base,
    const char* partTypeGroupName, const struct io_props props,
    const size_t N_total, const long long* N_counts,
    const enum lossy_compression_schemes lossy_compression,
    const struct unit_system* internal_units,
    const struct unit_system* snapshot_units) {

#if H5_VERSION_GE(1, 10, 0)

  /* Create data space */
  hid_t h_space;
  if (N_total > 0)
    h_space = H5Screate(H5S_SIMPLE);
  else
    h_space = H5Screate(H5S_NULL);

  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  int rank = 0;
  hsize_t shape[2];
  hsize_t source_shape[2];
  hsize_t start[2] = {0, 0};
  hsize_t count[2];
  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N_total;
    shape[1] = props.dimension;
    source_shape[0] = 0;
    source_shape[1] = props.dimension;
    count[0] = 0;
    count[1] = props.dimension;

  } else {
    rank = 1;
    shape[0] = N_total;
    shape[1] = 0;
    source_shape[0] = 0;
    source_shape[1] = 0;
    count[0] = 0;
    count[1] = 0;
  }

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, NULL);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset type */
  hid_t h_type = H5Tcopy(io_hdf5_type(props.type));

  /* Dataset properties */
  hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Create filters and set compression level if we have something to write */
  char comp_buffer[32] = "None";

  /* The name of the dataset to map to in the other files */
  char source_dataset_name[256];
  sprintf(source_dataset_name, "Groups/%s", props.name);

  /* Construct a relative base name */
  char fileName_relative_base[256];
  int pos_last_slash = strlen(fileName_base) - 1;
  for (/* */; pos_last_slash >= 0; --pos_last_slash)
    if (fileName_base[pos_last_slash] == '/') break;

  sprintf(fileName_relative_base, "%s", &fileName_base[pos_last_slash + 1]);

  /* Create all the virtual mappings */
  for (int i = 0; i < e->nr_nodes; ++i) {

    /* Get the number of particles of this type written on this rank */
    count[0] = N_counts[i];

    /* Select the space in the virtual file */
    h_err = H5Sselect_hyperslab(h_space, H5S_SELECT_SET, start, /*stride=*/NULL,
                                count, /*block=*/NULL);
    if (h_err < 0) error("Error selecting hyper-slab in the virtual file");

    /* Select the space in the (already existing) source file */
    source_shape[0] = count[0];
    hid_t h_source_space = H5Screate_simple(rank, source_shape, NULL);
    if (h_source_space < 0) error("Error creating space in the source file");

    char fileName[1024];
    sprintf(fileName, "%s.%d.hdf5", fileName_relative_base, i);

    /* Make the virtual link */
    h_err = H5Pset_virtual(h_prop, h_space, fileName, source_dataset_name,
                           h_source_space);
    if (h_err < 0) error("Error setting the virtual properties");

    H5Sclose(h_source_space);

    /* Move to the next slab (i.e. next file) */
    start[0] += count[0];
  }

  /* Create virtual dataset */
  const hid_t h_data = H5Dcreate(grp, props.name, h_type, h_space, H5P_DEFAULT,
                                 h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, props.units,
                              props.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, props.units);
  io_write_attribute_f(h_data, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(h_data, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(h_data, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(h_data, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(h_data, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(h_data, "h-scale exponent", 0.f);
  io_write_attribute_f(h_data, "a-scale exponent", props.scale_factor_exponent);
  io_write_attribute_s(h_data, "Expression for physical CGS units", buffer);
  io_write_attribute_s(h_data, "Lossy compression filter", comp_buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(snapshot_units, props.units);
  io_write_attribute_d(
      h_data,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      h_data,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, props.scale_factor_exponent));

#ifdef SWIFT_DEBUG_CHECKS
  if (strlen(props.description) == 0)
    error("Invalid (empty) description of the field '%s'", props.name);
#endif

  /* Write the full description */
  io_write_attribute_s(h_data, "Description", props.description);

  /* Close everything */
  H5Tclose(h_type);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);

#endif
}

void write_fof_virtual_file(const struct fof_props* props,
                            const size_t num_groups_total,
                            const long long* N_counts, const struct engine* e) {
#if H5_VERSION_GE(1, 10, 0)

  /* Create the filename */
  char file_name[512];
  char file_name_base[512];
  char subdir_name[265];
  sprintf(subdir_name, "%s_%04i", props->base_name, e->snapshot_output_count);
  const char* base = basename(subdir_name);
  sprintf(file_name_base, "%s/%s", subdir_name, base);
  sprintf(file_name, "%s/%s.hdf5", subdir_name, base);

  /* Open HDF5 file with the chosen parameters */
  hid_t h_file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", file_name);

  /* Start by writing the header */
  write_fof_hdf5_header(h_file, e, num_groups_total, num_groups_total, props,
                        /*virtual_file=*/1);
  hid_t h_grp =
      H5Gcreate(h_file, "/Groups", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating groups group.\n");

  struct io_props output_prop;
  output_prop = io_make_output_field_("Masses", DOUBLE, 1, UNIT_CONV_MASS, 0.f,
                                      (char*)props->group_mass, sizeof(double),
                                      "FOF group masses");
  write_virtual_fof_hdf5_array(e, h_grp, file_name_base, "Groups", output_prop,
                               num_groups_total, N_counts,
                               compression_write_lossless, e->internal_units,
                               e->snapshot_units);
  output_prop =
      io_make_output_field_("Centres", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f,
                            (char*)props->group_centre_of_mass,
                            3 * sizeof(double), "FOF group centres of mass");
  write_virtual_fof_hdf5_array(e, h_grp, file_name_base, "Groups", output_prop,
                               num_groups_total, N_counts,
                               compression_write_lossless, e->internal_units,
                               e->snapshot_units);
  output_prop = io_make_output_field_(
      "GroupIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
      (char*)props->group_index, sizeof(size_t), "FOF group IDs");
  write_virtual_fof_hdf5_array(e, h_grp, file_name_base, "Groups", output_prop,
                               num_groups_total, N_counts,
                               compression_write_lossless, e->internal_units,
                               e->snapshot_units);
  output_prop =
      io_make_output_field_("Sizes", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                            (char*)props->final_group_size, sizeof(long long),
                            "FOF group length (number of particles)");
  write_virtual_fof_hdf5_array(e, h_grp, file_name_base, "Groups", output_prop,
                               num_groups_total, N_counts,
                               compression_write_lossless, e->internal_units,
                               e->snapshot_units);

  /* Close everything */
  H5Gclose(h_grp);
  H5Fclose(h_file);
#endif
}

void write_fof_hdf5_array(
    const struct engine* e, hid_t grp, const char* fileName,
    const char* partTypeGroupName, const struct io_props props, const size_t N,
    const enum lossy_compression_schemes lossy_compression,
    const struct unit_system* internal_units,
    const struct unit_system* snapshot_units) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t num_elements = N * props.dimension;

  /* message("Writing '%s' array...", props.name); */

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");
  /* Copy the particle data to the temporary buffer */
  io_copy_temp_buffer(temp, e, props, N, internal_units, snapshot_units);

  /* Create data space */
  hid_t h_space;
  if (N > 0)
    h_space = H5Screate(H5S_SIMPLE);
  else
    h_space = H5Screate(H5S_NULL);

  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  /* Decide what chunk size to use based on compression */
  int log2_chunk_size = 20;

  int rank;
  hsize_t shape[2];
  hsize_t chunk_shape[2];

  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << log2_chunk_size;
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << log2_chunk_size;
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset type */
  hid_t h_type = H5Tcopy(io_hdf5_type(props.type));

  /* Dataset properties */
  hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Create filters and set compression level if we have something to write */
  char comp_buffer[32] = "None";
  if (N > 0) {

    /* Set chunk size */
    h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
    if (h_err < 0)
      error("Error while setting chunk size (%llu, %llu) for field '%s'.",
            (unsigned long long)chunk_shape[0],
            (unsigned long long)chunk_shape[1], props.name);

    /* Are we imposing some form of lossy compression filter? */
    if (lossy_compression != compression_write_lossless)
      set_hdf5_lossy_compression(&h_prop, &h_type, lossy_compression,
                                 props.name, comp_buffer);

    /* Impose GZIP data compression */
    if (e->snapshot_compression > 0) {
      h_err = H5Pset_shuffle(h_prop);
      if (h_err < 0)
        error("Error while setting shuffling options for field '%s'.",
              props.name);

      h_err = H5Pset_deflate(h_prop, e->snapshot_compression);
      if (h_err < 0)
        error("Error while setting compression options for field '%s'.",
              props.name);
    }

    /* Impose check-sum to verify data corruption */
    h_err = H5Pset_fletcher32(h_prop);
    if (h_err < 0)
      error("Error while setting checksum options for field '%s'.", props.name);
  }

  /* Create dataset */
  const hid_t h_data = H5Dcreate(grp, props.name, h_type, h_space, H5P_DEFAULT,
                                 h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_space, H5S_ALL,
                   H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while writing data array '%s'.", props.name);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, props.units,
                              props.scale_factor_exponent);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, props.units);
  io_write_attribute_f(h_data, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(h_data, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(h_data, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(h_data, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(h_data, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(h_data, "h-scale exponent", 0.f);
  io_write_attribute_f(h_data, "a-scale exponent", props.scale_factor_exponent);
  io_write_attribute_s(h_data, "Expression for physical CGS units", buffer);
  io_write_attribute_s(h_data, "Lossy compression filter", comp_buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double factor =
      units_cgs_conversion_factor(snapshot_units, props.units);
  io_write_attribute_d(
      h_data,
      "Conversion factor to CGS (not including cosmological corrections)",
      factor);
  io_write_attribute_d(
      h_data,
      "Conversion factor to physical CGS (including cosmological corrections)",
      factor * pow(e->cosmology->a, props.scale_factor_exponent));

#ifdef SWIFT_DEBUG_CHECKS
  if (strlen(props.description) == 0)
    error("Invalid (empty) description of the field '%s'", props.name);
#endif

  /* Write the full description */
  io_write_attribute_s(h_data, "Description", props.description);

  /* Free and close everything */
  swift_free("writebuff", temp);
  H5Tclose(h_type);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
}

void write_fof_hdf5_catalogue(const struct fof_props* props,
                              const size_t num_groups, const struct engine* e) {

  char file_name[512];
#ifdef WITH_MPI
  char subdir_name[265];
  sprintf(subdir_name, "%s_%04i", props->base_name, e->snapshot_output_count);
  if (e->nodeID == 0) safe_checkdir(subdir_name, /*create=*/1);
  MPI_Barrier(MPI_COMM_WORLD);

  const char* base = basename(subdir_name);
  sprintf(file_name, "%s/%s.%d.hdf5", subdir_name, base, e->nodeID);
#else
  sprintf(file_name, "%s_%04i.hdf5", props->base_name,
          e->snapshot_output_count);
#endif

  hid_t h_file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", file_name);

  /* Compute the number of groups */
  long long num_groups_local = num_groups;
  long long num_groups_total = num_groups;
#ifdef WITH_MPI
  MPI_Allreduce(&num_groups, &num_groups_total, 1, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);

  /* Rank 0 collects the number of groups written by each rank */
  long long* N_counts = (long long*)malloc(e->nr_nodes * sizeof(long long));
  MPI_Gather(&num_groups_local, 1, MPI_LONG_LONG_INT, N_counts, 1,
             MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
#endif

  /* Start by writing the header */
  write_fof_hdf5_header(h_file, e, num_groups_total, num_groups_local, props,
                        /*virtual_file=*/0);

  hid_t h_grp =
      H5Gcreate(h_file, "/Groups", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating groups group.\n");

  struct io_props output_prop;
  output_prop = io_make_output_field_("Masses", DOUBLE, 1, UNIT_CONV_MASS, 0.f,
                                      (char*)props->group_mass, sizeof(double),
                                      "FOF group masses");
  write_fof_hdf5_array(e, h_grp, file_name, "Groups", output_prop,
                       num_groups_local, compression_write_lossless,
                       e->internal_units, e->snapshot_units);
  output_prop =
      io_make_output_field_("Centres", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f,
                            (char*)props->group_centre_of_mass,
                            3 * sizeof(double), "FOF group centres of mass");
  write_fof_hdf5_array(e, h_grp, file_name, "Groups", output_prop,
                       num_groups_local, compression_write_lossless,
                       e->internal_units, e->snapshot_units);
  output_prop = io_make_output_field_(
      "GroupIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
      (char*)props->group_index, sizeof(size_t), "FOF group IDs");
  write_fof_hdf5_array(e, h_grp, file_name, "Groups", output_prop,
                       num_groups_local, compression_write_lossless,
                       e->internal_units, e->snapshot_units);
  output_prop =
      io_make_output_field_("Sizes", LONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                            (char*)props->final_group_size, sizeof(long long),
                            "FOF group length (number of particles)");
  write_fof_hdf5_array(e, h_grp, file_name, "Groups", output_prop,
                       num_groups_local, compression_write_lossless,
                       e->internal_units, e->snapshot_units);

  /* Close everything */
  H5Gclose(h_grp);
  H5Fclose(h_file);

#ifdef WITH_MPI
#if H5_VERSION_GE(1, 10, 0)

  /* Write the virtual meta-file */
  if (e->nodeID == 0)
    write_fof_virtual_file(props, num_groups_total, N_counts, e);
#endif

  /* Free the counts-per-rank array */
  free(N_counts);

  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

#endif /* HAVE_HDF5 */
