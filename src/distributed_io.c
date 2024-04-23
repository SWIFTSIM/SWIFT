/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#if defined(HAVE_HDF5) && defined(WITH_MPI)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

/* This object's header. */
#include "distributed_io.h"

/* Local includes. */
#include "black_holes_io.h"
#include "chemistry_io.h"
#include "common_io.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_compression.h"
#include "io_properties.h"
#include "memuse.h"
#include "output_list.h"
#include "output_options.h"
#include "part.h"
#include "part_type.h"
#include "sink_io.h"
#include "star_formation_io.h"
#include "stars_io.h"
#include "tools.h"
#include "units.h"
#include "version.h"
#include "xmf.h"

/* Are we timing the i/o? */
//#define IO_SPEED_MEASUREMENT

/* Max number of entries that can be written for a given particle type */
static const int io_max_size_output_list = 100;

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param e The #engine we are writing from.
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param props The #io_props of the field to read
 * @param N The number of particles to write.
 * @param lossy_compression Level of lossy compression to use for this field.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 * the part array will be written once the structures have been stabilized.
 */
void write_distributed_array(
    const struct engine* e, hid_t grp, const char* fileName,
    const char* partTypeGroupName, const struct io_props props, const size_t N,
    const enum lossy_compression_schemes lossy_compression,
    const struct unit_system* internal_units,
    const struct unit_system* snapshot_units) {

#ifdef IO_SPEED_MEASUREMENT
  const ticks tic_total = getticks();
#endif

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t num_elements = N * props.dimension;

  /* message("Writing '%s' array...", props.name); */

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (swift_memalign("writebuff", (void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

#ifdef IO_SPEED_MEASUREMENT
  ticks tic = getticks();
#endif

  /* Copy the particle data to the temporary buffer */
  io_copy_temp_buffer(temp, e, props, N, internal_units, snapshot_units);

#ifdef IO_SPEED_MEASUREMENT
  if (engine_rank == IO_SPEED_MEASUREMENT || IO_SPEED_MEASUREMENT == -1)
    message("Copying for '%s' took %.3f %s.", props.name,
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#endif

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

#ifdef IO_SPEED_MEASUREMENT
  tic = getticks();
#endif

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_space, H5S_ALL,
                   H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while writing data array '%s'.", props.name);

#ifdef IO_SPEED_MEASUREMENT
  ticks toc = getticks();
  float ms = clocks_from_ticks(toc - tic);
  int megaBytes = N * props.dimension * typeSize / (1024 * 1024);
  if (engine_rank == IO_SPEED_MEASUREMENT || IO_SPEED_MEASUREMENT == -1)
    message(
        "H5Dwrite for '%s' (%d MB) on rank %d took %.3f %s (speed = %f MB/s).",
        props.name, megaBytes, engine_rank, ms, clocks_getunit(),
        megaBytes / (ms / 1000.));
#endif

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

#ifdef IO_SPEED_MEASUREMENT
  if (engine_rank == IO_SPEED_MEASUREMENT || IO_SPEED_MEASUREMENT == -1)
    message("'%s' took %.3f %s.", props.name,
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#endif
}

/**
 * @brief Prepares an array in the snapshot.
 *
 * @param e The #engine we are writing from.
 * @param grp The HDF5 grp to write to.
 * @param fileName The name of the file we are writing to.
 * @param xmfFile The (opened) XMF file we are appending to.
 * @param partTypeGroupName The name of the group we are writing to.
 * @param props The #io_props of the field to write.
 * @param N_total The total number of particles to write in this array.
 * @param snapshot_units The units used for the data in this snapshot.
 */
void write_array_virtual(struct engine* e, hid_t grp, const char* fileName_base,
                         FILE* xmfFile, char* partTypeGroupName,
                         struct io_props props, long long N_total,
                         const long long* N_counts, const int num_ranks,
                         const int ptype,
                         const enum lossy_compression_schemes lossy_compression,
                         const struct unit_system* snapshot_units) {

#if H5_VERSION_GE(1, 10, 0)

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
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

  /* Are we imposing some form of lossy compression filter? */
  char comp_buffer[32] = "None";
  if (lossy_compression != compression_write_lossless)
    sprintf(comp_buffer, "%s",
            lossy_compression_schemes_names[lossy_compression]);

  /* The name of the dataset to map to in the other files */
  char source_dataset_name[256];
  sprintf(source_dataset_name, "PartType%d/%s", ptype, props.name);

  /* Construct a relative base name */
  char fileName_relative_base[256];
  int pos_last_slash = strlen(fileName_base) - 1;
  for (/* */; pos_last_slash >= 0; --pos_last_slash)
    if (fileName_base[pos_last_slash] == '/') break;

  sprintf(fileName_relative_base, "%s", &fileName_base[pos_last_slash + 1]);

  /* Create all the virtual mappings */
  for (int i = 0; i < num_ranks; ++i) {

    /* Get the number of particles of this type written on this rank */
    count[0] = N_counts[i * swift_type_count + ptype];

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

  /* Add a line to the XMF */
  if (xmfFile != NULL) {
    char fileName[1024];
    sprintf(fileName, "%s.hdf5", fileName_base);
    xmf_write_line(xmfFile, fileName, /*distributed=*/1, partTypeGroupName,
                   props.name, N_total, props.dimension, props.type);
  }

  /* Close everything */
  H5Tclose(h_type);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);

#else
  error(
      "Function cannot be called when the code is compiled with hdf5 older "
      "than 1.10.0");
#endif
}

/**
 * @brief Prepares a file for a parallel write.
 *
 * @param e The #engine.
 * @param fileName The file name to write to.
 * @param N_total The total number of particles of each type to write.
 * @param numFields The number of fields to write for each particle type.
 * @param internal_units The #unit_system used internally.
 * @param snapshot_units The #unit_system used in the snapshots.
 * @param fof Is this a snapshot related to a stand-alone FOF call?
 * @param subsample_any Are any fields being subsampled?
 * @param subsample_fraction The subsampling fraction of each particle type.
 */
void write_virtual_file(struct engine* e, const char* fileName_base,
                        const char* xmfFileName,
                        const long long N_total[swift_type_count],
                        const long long* N_counts, const int num_ranks,
                        const int to_write[swift_type_count],
                        const int numFields[swift_type_count],
                        char current_selection_name[FIELD_BUFFER_SIZE],
                        const struct unit_system* internal_units,
                        const struct unit_system* snapshot_units, const int fof,
                        const int subsample_any,
                        const float subsample_fraction[swift_type_count]) {

#if H5_VERSION_GE(1, 10, 0)

  struct output_options* output_options = e->output_options;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif
  const int with_rt = e->policy & engine_policy_rt;

  FILE* xmfFile = 0;
  int numFiles = 1;

  /* First time, we need to create the XMF file */
  if (e->snapshot_output_count == 0) xmf_create_file(xmfFileName);

  /* Prepare the XMF file for the new entry */
  xmfFile = xmf_prepare_file(xmfFileName);

  char fileName[1024];
  sprintf(fileName, "%s.hdf5", fileName_base);

  /* Write the part of the XMF file corresponding to this
   * specific output */
  xmf_write_outputheader(xmfFile, fileName, e->time);

  /* Open HDF5 file with the chosen parameters */
  hid_t h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  hid_t h_grp =
      H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_TIME);
  const double factor_length =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  const int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
  io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
  io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);
  io_write_attribute_s(h_grp, "System", hostname());
  io_write_attribute(h_grp, "Shift", DOUBLE, e->s->initial_shift, 3);

  /* Write out the particle types */
  io_write_part_type_names(h_grp);

  /* Write out the time-base */
  if (with_cosmology) {
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
  long long numParticlesThisFile[swift_type_count] = {0};
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};

  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);

    if (numFields[ptype] == 0) {
      numParticlesThisFile[ptype] = 0;
    } else {
      numParticlesThisFile[ptype] = N_total[ptype];
    }
  }

  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, numParticlesThisFile,
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
  io_write_attribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);
  io_write_attribute_i(h_grp, "ThisFile", 0);
  io_write_attribute_s(h_grp, "SelectOutput", current_selection_name);
  io_write_attribute_i(h_grp, "Virtual", 1);
  io_write_attribute(h_grp, "CanHaveTypes", INT, to_write, swift_type_count);

  if (subsample_any) {
    io_write_attribute_s(h_grp, "OutputType", "SubSampled");
    io_write_attribute(h_grp, "SubSampleFractions", FLOAT, subsample_fraction,
                       swift_type_count);
  } else {
    io_write_attribute_s(h_grp, "OutputType", "FullVolume");
  }

  /* Close header */
  H5Gclose(h_grp);

  /* Copy metadata from ICs to the file */
  ic_info_write_hdf5(e->ics_metadata, h_file);

  /* Write all the meta-data */
  io_write_meta_data(h_file, e, internal_units, snapshot_units, fof);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if there are
     * (a) no particles of this kind in this run, or
     * (b) if we have disabled every field of this particle type. */
    if (!to_write[ptype] || numFields[ptype] == 0) continue;

    /* Add the global information for that particle type to
     * the XMF meta-file */
    xmf_write_groupheader(xmfFile, fileName, /*distributed=*/1, N_total[ptype],
                          (enum part_type)ptype);

    /* Create the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0)
      error("Error while creating particle group %s.", partTypeGroupName);

    /* Add an alias name for convenience */
    char aliasName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(aliasName, PARTICLE_GROUP_BUFFER_SIZE, "/%sParticles",
             part_type_names[ptype]);
    hid_t h_err = H5Lcreate_soft(partTypeGroupName, h_grp, aliasName,
                                 H5P_DEFAULT, H5P_DEFAULT);
    if (h_err < 0) error("Error while creating alias for particle group.\n");

    /* Write the number of particles as an attribute */
    io_write_attribute_ll(h_grp, "NumberOfParticles", N_total[ptype]);
    io_write_attribute_ll(h_grp, "TotalNumberOfParticles", N_total[ptype]);

    int num_fields = 0;
    struct io_props list[io_max_size_output_list];
    bzero(list, io_max_size_output_list * sizeof(struct io_props));

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        io_select_hydro_fields(NULL, NULL, with_cosmology, with_cooling,
                               with_temperature, with_fof, with_stf, with_rt, e,
                               &num_fields, list);
        break;

      case swift_type_dark_matter:
      case swift_type_dark_matter_background:
        io_select_dm_fields(NULL, NULL, with_fof, with_stf, e, &num_fields,
                            list);
        break;

      case swift_type_neutrino:
        io_select_neutrino_fields(NULL, NULL, with_fof, with_stf, e,
                                  &num_fields, list);
        break;

      case swift_type_sink:
        io_select_sink_fields(NULL, with_cosmology, with_fof, with_stf, e,
                              &num_fields, list);
        break;

      case swift_type_stars:
        io_select_star_fields(NULL, with_cosmology, with_fof, with_stf, with_rt,
                              e, &num_fields, list);
        break;

      case swift_type_black_hole:
        io_select_bh_fields(NULL, with_cosmology, with_fof, with_stf, e,
                            &num_fields, list);
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Verify we are not going to crash when writing below */
    if (num_fields >= io_max_size_output_list)
      error("Too many fields to write for particle type %d", ptype);
    for (int i = 0; i < num_fields; ++i) {
      if (!list[i].is_used) error("List of field contains an empty entry!");
      if (!list[i].dimension)
        error("Dimension of field '%s' is <= 1!", list[i].name);
    }

    /* Did the user specify a non-standard default for the entire particle
     * type? */
    const enum lossy_compression_schemes compression_level_current_default =
        output_options_get_ptype_default_compression(
            output_options->select_output, current_selection_name,
            (enum part_type)ptype, e->verbose);

    /* Prepare everything that is not cancelled */
    int num_fields_written = 0;
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      const enum lossy_compression_schemes compression_level =
          output_options_get_field_compression(
              output_options, current_selection_name, list[i].name,
              (enum part_type)ptype, compression_level_current_default,
              e->verbose);

      if (compression_level != compression_do_not_write) {
        write_array_virtual(e, h_grp, fileName_base, xmfFile, partTypeGroupName,
                            list[i], N_total[ptype], N_counts, num_ranks, ptype,
                            compression_level, snapshot_units);
        num_fields_written++;
      }
    }

    /* Only write this now that we know exactly how many fields there are. */
    io_write_attribute_i(h_grp, "NumberOfFields", num_fields_written);

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
  }

  /* Write LXMF file descriptor */
  xmf_write_outputfooter(xmfFile, e->snapshot_output_count, e->time);

  /* Close the file for now */
  H5Fclose(h_file);

#else
  error(
      "Function cannot be called when the code is compiled with hdf5 older "
      "than 1.10.0");
#endif
}

/**
 * @brief Writes a snapshot distributed into multiple files.
 *
 * @param e The engine containing all the system.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param fof Is this a snapshot related to a stand-alone FOF call?
 * @param mpi_rank The rank number of the calling MPI rank.
 * @param mpi_size the number of MPI ranks.
 * @param comm The communicator used by the MPI ranks.
 * @param info The MPI information object.
 *
 * Creates a series of HDF5 output files (1 per MPI node) as a snapshot.
 * Writes the particles contained in the engine.
 * If such files already exist, it is erased and replaced by the new one.
 * The companion XMF file is also updated accordingly.
 */
void write_output_distributed(struct engine* e,
                              const struct unit_system* internal_units,
                              const struct unit_system* snapshot_units,
                              const int fof, const int mpi_rank,
                              const int mpi_size, MPI_Comm comm,
                              MPI_Info info) {

  hid_t h_file = 0, h_grp = 0;
  int numFiles = mpi_size;
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct sink* sinks = e->s->sinks;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  struct output_options* output_options = e->output_options;
  struct output_list* output_list = e->output_list_snapshots;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
  const int with_DM_background = e->s->with_DM_background;
  const int with_DM = e->s->with_DM;
  const int with_neutrinos = e->s->with_neutrinos;
  const int with_rt = e->policy & engine_policy_rt;
  const int with_hydro = (e->policy & engine_policy_hydro) ? 1 : 0;
  const int with_stars = (e->policy & engine_policy_stars) ? 1 : 0;
  const int with_black_hole = (e->policy & engine_policy_black_holes) ? 1 : 0;
  const int with_sink = (e->policy & engine_policy_sinks) ? 1 : 0;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nsinks = e->s->nr_sinks;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nblackholes = e->s->nr_bparts;

  /* Determine if we are writing a reduced snapshot, and if so which
   * output selection type to use */
  char current_selection_name[FIELD_BUFFER_SIZE] =
      select_output_header_default_name;
  if (output_list) {
    /* Users could have specified a different Select Output scheme for each
     * snapshot. */
    output_list_get_current_select_output(output_list, current_selection_name);
  }

  int snap_count = -1;
  int number_digits = -1;
  if (output_list && output_list->alternative_labels_on) {
    snap_count = output_list->snapshot_labels[snap_count];
    number_digits = 0;
  } else if (e->snapshot_invoke_stf) {
    snap_count = e->stf_output_count;
    number_digits = 4;
  } else {
    snap_count = e->snapshot_output_count;
    number_digits = 4;
  }

  /* Directory and file name */
  char dirName[1024];
  char fileName[1024];
  char fileName_base[1024];
  char xmfFileName[FILENAME_BUFFER_SIZE];
  char snapshot_subdir_name[FILENAME_BUFFER_SIZE];
  char snapshot_base_name[FILENAME_BUFFER_SIZE];

  output_options_get_basename(output_options, current_selection_name,
                              e->snapshot_subdir, e->snapshot_base_name,
                              snapshot_subdir_name, snapshot_base_name);

  io_get_snapshot_filename(
      fileName, xmfFileName, output_list, e->snapshot_invoke_stf,
      e->stf_output_count, e->snapshot_output_count, e->snapshot_subdir,
      snapshot_subdir_name, e->snapshot_base_name, snapshot_base_name);

  /* Are we using a sub-dir? */
  if (strnlen(e->snapshot_subdir, PARSER_MAX_LINE_SIZE) > 0) {
    sprintf(dirName, "%s/%s_%0*d", snapshot_subdir_name, snapshot_base_name,
            number_digits, snap_count);

    sprintf(fileName, "%s/%s_%0*d/%s_%0*d.%d.hdf5", snapshot_subdir_name,
            snapshot_base_name, number_digits, snap_count, snapshot_base_name,
            number_digits, snap_count, mpi_rank);

    sprintf(fileName_base, "%s/%s_%0*d/%s_%0*d", snapshot_subdir_name,
            snapshot_base_name, number_digits, snap_count, snapshot_base_name,
            number_digits, snap_count);

  } else {
    sprintf(dirName, "%s_%0*d", snapshot_base_name, number_digits, snap_count);

    sprintf(fileName, "%s_%0*d/%s_%0*d.%d.hdf5", snapshot_base_name,
            number_digits, snap_count, snapshot_base_name, number_digits,
            snap_count, mpi_rank);

    sprintf(fileName_base, "%s_%0*d/%s_%0*d", snapshot_base_name, number_digits,
            snap_count, snapshot_base_name, number_digits, snap_count);
  }

  /* Create the directory */
  if (mpi_rank == 0) safe_checkdir(snapshot_subdir_name, /*create=*/1);
  if (mpi_rank == 0) safe_checkdir(dirName, /*create=*/1);
  MPI_Barrier(comm);

  /* Do we want to sub-sample any of the arrays */
  int subsample[swift_type_count];
  float subsample_fraction[swift_type_count];
  for (int i = 0; i < swift_type_count; ++i) {
    subsample[i] = 0;
    subsample_fraction[i] = 1.f;
  }

  output_options_get_subsampling(
      output_options, current_selection_name, e->snapshot_subsample,
      e->snapshot_subsample_fraction, subsample, subsample_fraction);

  /* Is any particle type being subsampled? */
  int subsample_any = 0;
  for (int i = 0; i < swift_type_count; ++i) {
    subsample_any += subsample[i];
    if (!subsample[i]) subsample_fraction[i] = 1.f;
  }

  /* Number of particles that we will write */
  size_t Ngas_written, Ndm_written, Ndm_background, Ndm_neutrino,
      Nsinks_written, Nstars_written, Nblackholes_written;

  if (subsample[swift_type_gas]) {
    Ngas_written = io_count_gas_to_write(e->s, /*subsample=*/1,
                                         subsample_fraction[swift_type_gas],
                                         e->snapshot_output_count);
  } else {
    Ngas_written =
        e->s->nr_parts - e->s->nr_inhibited_parts - e->s->nr_extra_parts;
  }

  if (subsample[swift_type_stars]) {
    Nstars_written = io_count_stars_to_write(
        e->s, /*subsample=*/1, subsample_fraction[swift_type_stars],
        e->snapshot_output_count);
  } else {
    Nstars_written =
        e->s->nr_sparts - e->s->nr_inhibited_sparts - e->s->nr_extra_sparts;
  }

  if (subsample[swift_type_black_hole]) {
    Nblackholes_written = io_count_black_holes_to_write(
        e->s, /*subsample=*/1, subsample_fraction[swift_type_black_hole],
        e->snapshot_output_count);
  } else {
    Nblackholes_written =
        e->s->nr_bparts - e->s->nr_inhibited_bparts - e->s->nr_extra_bparts;
  }

  if (subsample[swift_type_sink]) {
    Nsinks_written = io_count_sinks_to_write(
        e->s, /*subsample=*/1, subsample_fraction[swift_type_sink],
        e->snapshot_output_count);
  } else {
    Nsinks_written =
        e->s->nr_sinks - e->s->nr_inhibited_sinks - e->s->nr_extra_sinks;
  }

  Ndm_written = io_count_dark_matter_to_write(
      e->s, subsample[swift_type_dark_matter],
      subsample_fraction[swift_type_dark_matter], e->snapshot_output_count);

  if (with_DM_background) {
    Ndm_background = io_count_background_dark_matter_to_write(
        e->s, subsample[swift_type_dark_matter_background],
        subsample_fraction[swift_type_dark_matter_background],
        e->snapshot_output_count);
  } else {
    Ndm_background = 0;
  }

  if (with_neutrinos) {
    Ndm_neutrino = io_count_neutrinos_to_write(
        e->s, subsample[swift_type_neutrino],
        subsample_fraction[swift_type_neutrino], e->snapshot_output_count);
  } else {
    Ndm_neutrino = 0;
  }

  /* Compute offset in the file and total number of particles */
  const long long N[swift_type_count] = {
      Ngas_written,   Ndm_written,         Ndm_background, Nsinks_written,
      Nstars_written, Nblackholes_written, Ndm_neutrino};

  /* Gather the total number of particles to write */
  long long N_total[swift_type_count] = {0};
  MPI_Allreduce(N, N_total, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);

  /* Collect the number of particles written by each rank */
  long long* N_counts =
      (long long*)malloc(mpi_size * swift_type_count * sizeof(long long));
  MPI_Gather(N, swift_type_count, MPI_LONG_LONG_INT, N_counts, swift_type_count,
             MPI_LONG_LONG_INT, 0, comm);

  /* List what fields to write.
   * Note that we want to want to write a 0-size dataset for some species
   * in case future snapshots will contain them (e.g. star formation) */
  const int to_write[swift_type_count] = {
      with_hydro, with_DM,         with_DM_background, with_sink,
      with_stars, with_black_hole, with_neutrinos

  };

  /* Use a single Lustre stripe with a rank-based OST offset? */
  if (e->snapshot_lustre_OST_count != 0) {

    /* Use a random offset to avoid placing things in the same OSTs. We do
     * this to keep the use of OSTs balanced, much like using -1 for the
     * stripe. */
    int offset = rand() % e->snapshot_lustre_OST_count;
    MPI_Bcast(&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);

    char string[1200];
    sprintf(string, "lfs setstripe -c 1 -i %d %s",
            ((e->nodeID + offset) % e->snapshot_lustre_OST_count), fileName);
    const int result = system(string);
    if (result != 0) {
      message("lfs setstripe command returned error code %d", result);
    }
  }

  /* Open file */
  /* message("Opening file '%s'.", fileName); */
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Convert basic output information to snapshot units */
  const double factor_time =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_TIME);
  const double factor_length =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);
  const double dblTime = e->time * factor_time;
  const double dim[3] = {e->s->dim[0] * factor_length,
                         e->s->dim[1] * factor_length,
                         e->s->dim[2] * factor_length};

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, dim, 3);
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  const int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);
  io_write_attribute(h_grp, "Redshift", DOUBLE, &e->cosmology->z, 1);
  io_write_attribute(h_grp, "Scale-factor", DOUBLE, &e->cosmology->a, 1);
  io_write_attribute_s(h_grp, "Code", "SWIFT");
  io_write_attribute_s(h_grp, "RunName", e->run_name);

  /* We write rank 0's hostname so that it is uniform across all files. */
  char systemname[256] = {0};
  if (mpi_rank == 0) sprintf(systemname, "%s", hostname());
  MPI_Bcast(systemname, 256, MPI_CHAR, 0, comm);
  io_write_attribute_s(h_grp, "System", systemname);
  io_write_attribute(h_grp, "Shift", DOUBLE, e->s->initial_shift, 3);

  /* Write out the particle types */
  io_write_part_type_names(h_grp);

  /* Write out the time-base */
  if (with_cosmology) {
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

  /* GADGET-2 legacy values:  Number of particles of each type */
  long long numParticlesThisFile[swift_type_count] = {0};
  unsigned int numParticles[swift_type_count] = {0};
  unsigned int numParticlesHighWord[swift_type_count] = {0};

  /* Total number of fields to write per ptype */
  int numFields[swift_type_count] = {0};

  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);

    numFields[ptype] = output_options_get_num_fields_to_write(
        output_options, current_selection_name, ptype);

    if (numFields[ptype] == 0) {
      numParticlesThisFile[ptype] = 0;
    } else {
      numParticlesThisFile[ptype] = N[ptype];
    }
  }

  io_write_attribute(h_grp, "NumPart_ThisFile", LONGLONG, numParticlesThisFile,
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
  io_write_attribute_i(h_grp, "NumFilesPerSnapshot", numFiles);
  io_write_attribute_i(h_grp, "ThisFile", mpi_rank);
  io_write_attribute_s(h_grp, "SelectOutput", current_selection_name);
  io_write_attribute_i(h_grp, "Virtual", 0);
  io_write_attribute(h_grp, "CanHaveTypes", INT, to_write, swift_type_count);

  if (subsample_any) {
    io_write_attribute_s(h_grp, "OutputType", "SubSampled");
    io_write_attribute(h_grp, "SubSampleFractions", FLOAT, subsample_fraction,
                       swift_type_count);
  } else {
    io_write_attribute_s(h_grp, "OutputType", "FullVolume");
  }

  /* Close header */
  H5Gclose(h_grp);

  /* Copy metadata from ICs to the file */
  ic_info_write_hdf5(e->ics_metadata, h_file);

  /* Write all the meta-data */
  io_write_meta_data(h_file, e, internal_units, snapshot_units, fof);

  /* Now write the top-level cell structure
   * We use a global offset of 0 here. This means that the cells will write
   * their offset with respect to the start of the file they belong to and
   * not a global offset */
  long long global_offsets[swift_type_count] = {0};
  h_grp = H5Gcreate(h_file, "/Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cells group");

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp, e->s->cdim, e->s->dim, e->s->cells_top,
                        e->s->nr_cells, e->s->width, mpi_rank,
                        /*distributed=*/1, subsample, subsample_fraction,
                        e->snapshot_output_count, N_total, global_offsets,
                        to_write, numFields, internal_units, snapshot_units);
  H5Gclose(h_grp);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if there are
     * (a) no particles of this kind in this run, or
     * (b) if we have disabled every field of this particle type. */
    if (!to_write[ptype] || numFields[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating particle group.\n");

    /* Add an alias name for convenience */
    char aliasName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(aliasName, PARTICLE_GROUP_BUFFER_SIZE, "/%sParticles",
             part_type_names[ptype]);
    hid_t h_err = H5Lcreate_soft(partTypeGroupName, h_grp, aliasName,
                                 H5P_DEFAULT, H5P_DEFAULT);
    if (h_err < 0) error("Error while creating alias for particle group.\n");

    /* Write the number of particles as an attribute */
    io_write_attribute_ll(h_grp, "NumberOfParticles", N[ptype]);
    io_write_attribute_ll(h_grp, "TotalNumberOfParticles", N_total[ptype]);

    int num_fields = 0;
    struct io_props list[io_max_size_output_list];
    bzero(list, io_max_size_output_list * sizeof(struct io_props));
    size_t Nparticles = 0;

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct sink* sinks_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas: {
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          Nparticles = Ngas;

          /* Select the fields to write */
          io_select_hydro_fields(parts, xparts, with_cosmology, with_cooling,
                                 with_temperature, with_fof, with_stf, with_rt,
                                 e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Ngas_written;

          /* Allocate temporary arrays */
          if (swift_memalign("parts_written", (void**)&parts_written,
                             part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporary memory for parts");
          if (swift_memalign("xparts_written", (void**)&xparts_written,
                             xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporary memory for xparts");

          /* Collect the particles we want to write */
          io_collect_parts_to_write(
              parts, xparts, parts_written, xparts_written,
              subsample[swift_type_gas], subsample_fraction[swift_type_gas],
              e->snapshot_output_count, Ngas, Ngas_written);

          /* Select the fields to write */
          io_select_hydro_fields(parts_written, xparts_written, with_cosmology,
                                 with_cooling, with_temperature, with_fof,
                                 with_stf, with_rt, e, &num_fields, list);
        }
      } break;

      case swift_type_dark_matter: {
        if (Ntot == Ndm_written) {

          /* This is a DM-only run without background or inhibited particles
           * or neutrinos */
          Nparticles = Ntot;

          /* Select the fields to write */
          io_select_dm_fields(gparts, e->s->gpart_group_data, with_fof,
                              with_stf, e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Ndm_written;

          /* Allocate temporary array */
          if (swift_memalign("gparts_written", (void**)&gparts_written,
                             gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporary memory for gparts");

          if (with_stf) {
            if (swift_memalign(
                    "gpart_group_written", (void**)&gpart_group_data_written,
                    gpart_align,
                    Ndm_written * sizeof(struct velociraptor_gpart_data)) != 0)
              error(
                  "Error while allocating temporary memory for gparts STF "
                  "data");
          }

          /* Collect the non-inhibited DM particles from gpart */
          io_collect_gparts_to_write(
              gparts, e->s->gpart_group_data, gparts_written,
              gpart_group_data_written, subsample[swift_type_dark_matter],
              subsample_fraction[swift_type_dark_matter],
              e->snapshot_output_count, Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          io_select_dm_fields(gparts_written, gpart_group_data_written,
                              with_fof, with_stf, e, &num_fields, list);
        }
      } break;

      case swift_type_dark_matter_background: {

        /* Ok, we need to fish out the particles we want */
        Nparticles = Ndm_background;

        /* Allocate temporary array */
        if (swift_memalign("gparts_written", (void**)&gparts_written,
                           gpart_align,
                           Ndm_background * sizeof(struct gpart)) != 0)
          error("Error while allocating temporart memory for gparts");

        if (with_stf) {
          if (swift_memalign(
                  "gpart_group_written", (void**)&gpart_group_data_written,
                  gpart_align,
                  Ndm_background * sizeof(struct velociraptor_gpart_data)) != 0)
            error(
                "Error while allocating temporart memory for gparts STF "
                "data");
        }

        /* Collect the non-inhibited DM particles from gpart */
        io_collect_gparts_background_to_write(
            gparts, e->s->gpart_group_data, gparts_written,
            gpart_group_data_written,
            subsample[swift_type_dark_matter_background],
            subsample_fraction[swift_type_dark_matter_background],
            e->snapshot_output_count, Ntot, Ndm_background, with_stf);

        /* Select the fields to write */
        io_select_dm_fields(gparts_written, gpart_group_data_written, with_fof,
                            with_stf, e, &num_fields, list);
      } break;

      case swift_type_neutrino: {

        /* Ok, we need to fish out the particles we want */
        Nparticles = Ndm_neutrino;

        /* Allocate temporary array */
        if (swift_memalign("gparts_written", (void**)&gparts_written,
                           gpart_align,
                           Ndm_neutrino * sizeof(struct gpart)) != 0)
          error("Error while allocating temporart memory for gparts");

        if (with_stf) {
          if (swift_memalign(
                  "gpart_group_written", (void**)&gpart_group_data_written,
                  gpart_align,
                  Ndm_neutrino * sizeof(struct velociraptor_gpart_data)) != 0)
            error(
                "Error while allocating temporart memory for gparts STF "
                "data");
        }

        /* Collect the non-inhibited DM particles from gpart */
        io_collect_gparts_neutrino_to_write(
            gparts, e->s->gpart_group_data, gparts_written,
            gpart_group_data_written, subsample[swift_type_neutrino],
            subsample_fraction[swift_type_neutrino], e->snapshot_output_count,
            Ntot, Ndm_neutrino, with_stf);

        /* Select the fields to write */
        io_select_neutrino_fields(gparts_written, gpart_group_data_written,
                                  with_fof, with_stf, e, &num_fields, list);
      } break;

      case swift_type_sink: {
        if (Nsinks == Nsinks_written) {

          /* No inhibted particles: easy case */
          Nparticles = Nsinks;

          /* Select the fields to write */
          io_select_sink_fields(sinks, with_cosmology, with_fof, with_stf, e,
                                &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nsinks_written;

          /* Allocate temporary arrays */
          if (swift_memalign("sinks_written", (void**)&sinks_written,
                             sink_align,
                             Nsinks_written * sizeof(struct sink)) != 0)
            error("Error while allocating temporary memory for sinks");

          /* Collect the particles we want to write */
          io_collect_sinks_to_write(
              sinks, sinks_written, subsample[swift_type_sink],
              subsample_fraction[swift_type_sink], e->snapshot_output_count,
              Nsinks, Nsinks_written);

          /* Select the fields to write */
          io_select_sink_fields(sinks_written, with_cosmology, with_fof,
                                with_stf, e, &num_fields, list);
        }
      } break;

      case swift_type_stars: {
        if (Nstars == Nstars_written) {

          /* No inhibted particles: easy case */
          Nparticles = Nstars;

          /* Select the fields to write */
          io_select_star_fields(sparts, with_cosmology, with_fof, with_stf,
                                with_rt, e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nstars_written;

          /* Allocate temporary arrays */
          if (swift_memalign("sparts_written", (void**)&sparts_written,
                             spart_align,
                             Nstars_written * sizeof(struct spart)) != 0)
            error("Error while allocating temporary memory for sparts");

          /* Collect the particles we want to write */
          io_collect_sparts_to_write(
              sparts, sparts_written, subsample[swift_type_stars],
              subsample_fraction[swift_type_stars], e->snapshot_output_count,
              Nstars, Nstars_written);

          /* Select the fields to write */
          io_select_star_fields(sparts_written, with_cosmology, with_fof,
                                with_stf, with_rt, e, &num_fields, list);
        }
      } break;

      case swift_type_black_hole: {
        if (Nblackholes == Nblackholes_written) {

          /* No inhibted particles: easy case */
          Nparticles = Nblackholes;

          /* Select the fields to write */
          io_select_bh_fields(bparts, with_cosmology, with_fof, with_stf, e,
                              &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          Nparticles = Nblackholes_written;

          /* Allocate temporary arrays */
          if (swift_memalign("bparts_written", (void**)&bparts_written,
                             bpart_align,
                             Nblackholes_written * sizeof(struct bpart)) != 0)
            error("Error while allocating temporary memory for bparts");

          /* Collect the particles we want to write */
          io_collect_bparts_to_write(
              bparts, bparts_written, subsample[swift_type_black_hole],
              subsample_fraction[swift_type_black_hole],
              e->snapshot_output_count, Nblackholes, Nblackholes_written);

          /* Select the fields to write */
          io_select_bh_fields(bparts_written, with_cosmology, with_fof,
                              with_stf, e, &num_fields, list);
        }
      } break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Verify we are not going to crash when writing below */
    if (num_fields >= io_max_size_output_list)
      error("Too many fields to write for particle type %d", ptype);
    for (int i = 0; i < num_fields; ++i) {
      if (!list[i].is_used) error("List of field contains an empty entry!");
      if (!list[i].dimension)
        error("Dimension of field '%s' is <= 1!", list[i].name);
    }

    /* Did the user specify a non-standard default for the entire particle
     * type? */
    const enum lossy_compression_schemes compression_level_current_default =
        output_options_get_ptype_default_compression(
            output_options->select_output, current_selection_name,
            (enum part_type)ptype, e->verbose);

    /* Write everything that is not cancelled */
    int num_fields_written = 0;
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      const enum lossy_compression_schemes compression_level =
          output_options_get_field_compression(
              output_options, current_selection_name, list[i].name,
              (enum part_type)ptype, compression_level_current_default,
              e->verbose);

      if (compression_level != compression_do_not_write) {
        write_distributed_array(e, h_grp, fileName, partTypeGroupName, list[i],
                                Nparticles, compression_level, internal_units,
                                snapshot_units);
        num_fields_written++;
      }
    }

    /* Only write this now that we know exactly how many fields there are. */
    io_write_attribute_i(h_grp, "NumberOfFields", num_fields_written);

    /* Free temporary arrays */
    if (parts_written) swift_free("parts_written", parts_written);
    if (xparts_written) swift_free("xparts_written", xparts_written);
    if (gparts_written) swift_free("gparts_written", gparts_written);
    if (gpart_group_data_written)
      swift_free("gpart_group_written", gpart_group_data_written);
    if (sinks_written) swift_free("sinks_written", sinks_written);
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

#if H5_VERSION_GE(1, 10, 0)

  /* Write the virtual meta-file */
  if (mpi_rank == 0)
    write_virtual_file(e, fileName_base, xmfFileName, N_total, N_counts,
                       mpi_size, to_write, numFields, current_selection_name,
                       internal_units, snapshot_units, fof, subsample_any,
                       subsample_fraction);

  /* Make sure nobody is allowed to progress until rank 0 is done. */
  MPI_Barrier(comm);

  /* Now write the top-level cell structure in the virtual file
   * but this time, it is *not* distributed. i.e. all the offsets are
   * in the virtual file */
  hid_t h_file_cells = 0, h_grp_cells = 0;
  if (mpi_rank == 0) {

    char fileName_virtual[1030];
    sprintf(fileName_virtual, "%s.hdf5", fileName_base);

    /* Open the snapshot on rank 0 */
    h_file_cells = H5Fopen(fileName_virtual, H5F_ACC_RDWR, H5P_DEFAULT);
    if (h_file_cells < 0)
      error("Error while opening file '%s' on rank %d.", fileName_virtual,
            mpi_rank);

    /* Create the group we want in the file */
    h_grp_cells = H5Gcreate(h_file_cells, "/Cells", H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (h_grp_cells < 0) error("Error while creating cells group");
  }

  /* We need to recompute the offsets since they are now with respect
   * to a single file. */
  for (int i = 0; i < swift_type_count; ++i) global_offsets[i] = 0;
  MPI_Exscan(N, global_offsets, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM,
             comm);

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp_cells, e->s->cdim, e->s->dim, e->s->cells_top,
                        e->s->nr_cells, e->s->width, mpi_rank,
                        /*distributed=*/0, subsample, subsample_fraction,
                        e->snapshot_output_count, N_total, global_offsets,
                        to_write, numFields, internal_units, snapshot_units);

  /* Close everything */
  if (mpi_rank == 0) {
    H5Gclose(h_grp_cells);
    H5Fclose(h_file_cells);
  }

#endif

  /* Free the counts-per-rank array */
  free(N_counts);

  /* Make sure nobody is allowed to progress until everyone is done. */
  MPI_Barrier(comm);

  e->snapshot_output_count++;
  if (e->snapshot_invoke_stf) e->stf_output_count++;
}

#endif /* HAVE_HDF5 && WITH_MPI */
