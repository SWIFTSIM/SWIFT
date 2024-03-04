/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

#if defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* This object's header. */
#include "single_io.h"

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
#include "mhd_io.h"
#include "output_list.h"
#include "output_options.h"
#include "part.h"
#include "part_type.h"
#include "rt_io.h"
#include "sink_io.h"
#include "star_formation_io.h"
#include "stars_io.h"
#include "tools.h"
#include "units.h"
#include "version.h"
#include "xmf.h"

/* Max number of entries that can be written for a given particle type */
static const int io_max_size_output_list = 100;

/**
 * @brief Reads a data array from a given HDF5 group.
 *
 * @param h_grp The group from which to read.
 * @param prop The #io_props of the field to read
 * @param N The number of particles.
 * @param internal_units The #unit_system used internally
 * @param ic_units The #unit_system used in the ICs
 * @param cleanup_h Are we removing h-factors from the ICs?
 * @param cleanup_sqrt_a Are we cleaning-up the sqrt(a) factors in the Gadget
 * IC velocities?
 * @param h The value of the reduced Hubble constant.
 * @param a The current value of the scale-factor.
 *
 * @todo A better version using HDF5 hyper-slabs to read the file directly into
 * the part array will be written once the structures have been stabilized.
 */
void read_array_single(hid_t h_grp, const struct io_props props, size_t N,
                       const struct unit_system* internal_units,
                       const struct unit_system* ic_units, int cleanup_h,
                       int cleanup_sqrt_a, double h, double a) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;
  const size_t num_elements = N * props.dimension;

  /* Check whether the dataspace exists or not */
  const htri_t exist = H5Lexists(h_grp, props.name, 0);
  if (exist < 0) {
    error("Error while checking the existence of data set '%s'.", props.name);
  } else if (exist == 0) {
    if (props.importance == COMPULSORY) {
      error("Compulsory data set '%s' not present in the file.", props.name);
    } else {

      /* Create a single instance of the default value */
      float* temp = (float*)malloc(copySize);
      for (int i = 0; i < props.dimension; ++i) temp[i] = props.default_value;

      /* Copy it everywhere in the particle array */
      for (size_t i = 0; i < N; ++i)
        memcpy(props.field + i * props.partSize, temp, copySize);

      free(temp);
      return;
    }
  }

  /* message("Reading %s '%s' array...", */
  /*         props.importance == COMPULSORY ? "compulsory" : "optional  ", */
  /*         props.name); */

  /* Open data space */
  const hid_t h_data = H5Dopen(h_grp, props.name, H5P_DEFAULT);
  if (h_data < 0) error("Error while opening data space '%s'.", props.name);

  /* Allocate temporary buffer */
  void* temp = malloc(num_elements * typeSize);
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Read HDF5 dataspace in temporary buffer */
  /* Dirty version that happens to work for vectors but should be improved */
  /* Using HDF5 dataspaces would be better */
  const hid_t h_err = H5Dread(h_data, io_hdf5_type(props.type), H5S_ALL,
                              H5S_ALL, H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while reading data array '%s'.", props.name);

  /* Unit conversion if necessary */
  const double unit_factor =
      units_conversion_factor(ic_units, internal_units, props.units);
  if (unit_factor != 1. && exist != 0) {

    /* message("Converting ! factor=%e", factor); */

    if (io_is_double_precision(props.type)) {
      double* temp_d = (double*)temp;
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= unit_factor;

    } else {
      float* temp_f = (float*)temp;

#ifdef SWIFT_DEBUG_CHECKS
      float maximum = 0.f;
      float minimum = FLT_MAX;
#endif

      /* Loop that converts the Units */
      for (size_t i = 0; i < num_elements; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Find the absolute minimum and maximum values */
        const float abstemp_f = fabsf(temp_f[i]);
        if (abstemp_f != 0.f) {
          maximum = max(maximum, abstemp_f);
          minimum = min(minimum, abstemp_f);
        }
#endif

        /* Convert the float units */
        temp_f[i] *= unit_factor;
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* The two possible errors: larger than float or smaller
       * than float precision. */
      if (unit_factor * maximum > FLT_MAX) {
        error("Unit conversion results in numbers larger than floats");
      } else if (unit_factor * minimum < FLT_MIN) {
        error("Numbers smaller than float precision");
      }
#endif
    }
  }

  /* Clean-up h if necessary */
  const float h_factor_exp = units_h_factor(internal_units, props.units);
  if (cleanup_h && h_factor_exp != 0.f && exist != 0) {

    /* message("Multipltying '%s' by h^%f=%f", props.name, h_factor_exp,
     * h_factor); */

    if (io_is_double_precision(props.type)) {
      double* temp_d = (double*)temp;
      const double h_factor = pow(h, h_factor_exp);
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= h_factor;
    } else {
      float* temp_f = (float*)temp;
      const float h_factor = pow(h, h_factor_exp);
      for (size_t i = 0; i < num_elements; ++i) temp_f[i] *= h_factor;
    }
  }

  /* Clean-up a if necessary */
  if (cleanup_sqrt_a && a != 1. && (strcmp(props.name, "Velocities") == 0)) {

    if (io_is_double_precision(props.type)) {
      double* temp_d = (double*)temp;
      const double vel_factor = sqrt(a);
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= vel_factor;
    } else {
      float* temp_f = (float*)temp;
      const float vel_factor = sqrt(a);
      for (size_t i = 0; i < num_elements; ++i) temp_f[i] *= vel_factor;
    }
  }

  /* Copy temporary buffer to particle data */
  char* temp_c = (char*)temp;
  for (size_t i = 0; i < N; ++i)
    memcpy(props.field + i * props.partSize, &temp_c[i * copySize], copySize);

  /* Free and close everything */
  free(temp);
  H5Dclose(h_data);
}

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param e The #engine we are writing from.
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param partTypeGroupName The name of the group containing the particles in
 * the HDF5 file.
 * @param props The #io_props of the field to read
 * @param N The number of particles to write.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 * the part array will be written once the structures have been stabilized.
 */
void write_array_single(const struct engine* e, hid_t grp, const char* fileName,
                        FILE* xmfFile, const char* partTypeGroupName,
                        const struct io_props props, const size_t N,
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
  const hid_t h_space = H5Screate(H5S_SIMPLE);
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

    /* Impose GZIP and shuffle data compression */
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

  /* Write XMF description for this data set */
  if (xmfFile != NULL)
    xmf_write_line(xmfFile, fileName, /*distributed=*/0, partTypeGroupName,
                   props.name, N, props.dimension, props.type);

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

/**
 * @brief Reads an HDF5 initial condition file (GADGET-3 type)
 *
 * @param fileName The file to read.
 * @param internal_units The system units used internally
 * @param dim (output) The dimension of the volume.
 * @param parts (output) Array of #part particles.
 * @param gparts (output) Array of #gpart particles.
 * @param sinks (output) Array of #sink particles.
 * @param sparts (output) Array of #spart particles.
 * @param bparts (output) Array of #bpart particles.
 * @param Ngas (output) number of Gas particles read.
 * @param Ngparts (output) The number of #gpart read.
 * @param Ngparts_background (output) The number of background #gpart (type 2)
 * @param Nnuparts (output) The number of neutrino #gpart (type 6)
 * @param Nsinks (output) The number of #sink read.
 * @param Nstars (output) The number of #spart read.
 * @param Nblackholes (output) The number of #bpart read.
 * @param flag_entropy (output) 1 if the ICs contained Entropy in the
 * InternalEnergy field
 * @param with_hydro Are we reading gas particles ?
 * @param with_gravity Are we reading/creating #gpart arrays ?
 * @param with_sink Are we reading sink particles ?
 * @param with_stars Are we reading star particles ?
 * @param with_black_hole Are we reading black hole particles ?
 * @param with_cosmology Are we running with cosmology ?
 * @param cleanup_h Are we cleaning-up h-factors from the quantities we read?
 * @param cleanup_sqrt_a Are we cleaning-up the sqrt(a) factors in the Gadget
 * IC velocities?
 * @param h The value of the reduced Hubble constant to use for correction.
 * @param a The current value of the scale-factor.
 * @prarm n_threads The number of threads to use for the temporary threadpool.
 * @param dry_run If 1, don't read the particle. Only allocates the arrays.
 * @param remap_ids Are we ignoring the ICs' IDs and remapping them to [1, N[ ?
 * @param ics_metadata Will store metadata group copied from the ICs file
 *
 * Opens the HDF5 file fileName and reads the particles contained
 * in the parts array. N is the returned number of particles found
 * in the file.
 *
 * @warning Can not read snapshot distributed over more than 1 file !!!
 * @todo Read snapshots distributed in more than one file.
 */
void read_ic_single(
    const char* fileName, const struct unit_system* internal_units,
    double dim[3], struct part** parts, struct gpart** gparts,
    struct sink** sinks, struct spart** sparts, struct bpart** bparts,
    size_t* Ngas, size_t* Ngparts, size_t* Ngparts_background, size_t* Nnuparts,
    size_t* Nsinks, size_t* Nstars, size_t* Nblackholes, int* flag_entropy,
    const int with_hydro, const int with_gravity, const int with_sink,
    const int with_stars, const int with_black_holes, const int with_cosmology,
    const int cleanup_h, const int cleanup_sqrt_a, const double h,
    const double a, const int n_threads, const int dry_run, const int remap_ids,
    struct ic_info* ics_metadata) {

  hid_t h_file = 0, h_grp = 0;
  /* GADGET has only cubic boxes (in cosmological mode) */
  double boxSize[3] = {0.0, -1.0, -1.0};
  long long numParticles[swift_type_count] = {0};
  long long numParticles_highWord[swift_type_count] = {0};
  size_t N[swift_type_count] = {0};
  int dimension = 3; /* Assume 3D if nothing is specified */
  size_t Ndm = 0;
  size_t Ndm_background = 0;
  size_t Ndm_neutrino = 0;

  /* Initialise counters */
  *Ngas = 0, *Ngparts = 0, *Ngparts_background = 0, *Nstars = 0,
  *Nblackholes = 0, *Nsinks = 0, *Nnuparts = 0;

  /* Open file */
  /* message("Opening file '%s' as IC.", fileName); */
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) error("Error while opening file '%s'.", fileName);

  /* Open header to read simulation properties */
  /* message("Reading file header..."); */
  h_grp = H5Gopen(h_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening file header\n");

  /* Check the dimensionality of the ICs (if the info exists) */
  const hid_t hid_dim = H5Aexists(h_grp, "Dimension");
  if (hid_dim < 0)
    error("Error while testing existance of 'Dimension' attribute");
  if (hid_dim > 0) io_read_attribute(h_grp, "Dimension", INT, &dimension);
  if (dimension != hydro_dimension)
    error("ICs dimensionality (%dD) does not match code dimensionality (%dD)",
          dimension, (int)hydro_dimension);

  /* Check whether the number of files is specified (if the info exists) */
  const hid_t hid_files = H5Aexists(h_grp, "NumFilesPerSnapshot");
  int num_files = 1;
  if (hid_files < 0)
    error(
        "Error while testing the existance of 'NumFilesPerSnapshot' attribute");
  if (hid_files > 0)
    io_read_attribute(h_grp, "NumFilesPerSnapshot", INT, &num_files);
  if (num_files != 1)
    error(
        "ICs are split over multiples files (%d). SWIFT cannot handle this "
        "case. The script /tools/combine_ics.py is availalbe in the repository "
        "to combine files into a valid input file.",
        num_files);

  /* Read the relevant information and print status */
  int flag_entropy_temp[swift_type_count];
  io_read_attribute(h_grp, "Flag_Entropy_ICs", INT, flag_entropy_temp);
  *flag_entropy = flag_entropy_temp[0];
  io_read_attribute(h_grp, "BoxSize", DOUBLE, boxSize);
  io_read_attribute(h_grp, "NumPart_Total", LONGLONG, numParticles);
  io_read_attribute(h_grp, "NumPart_Total_HighWord", LONGLONG,
                    numParticles_highWord);

  /* Check that the user is not doing something silly when they e.g. restart
   * from a snapshot by asserting that the current scale-factor (from
   * parameter file) and the redshift in the header are consistent */
  if (with_cosmology) {
    io_assert_valid_header_cosmology(h_grp, a);
  }

  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    N[ptype] = (numParticles[ptype]) + (numParticles_highWord[ptype] << 32);

  /* Get the box size if not cubic */
  dim[0] = boxSize[0];
  dim[1] = (boxSize[1] < 0) ? boxSize[0] : boxSize[1];
  dim[2] = (boxSize[2] < 0) ? boxSize[0] : boxSize[2];

  /* Change box size in the 1D and 2D case */
  if (hydro_dimension == 2)
    dim[2] = min(dim[0], dim[1]);
  else if (hydro_dimension == 1)
    dim[2] = dim[1] = dim[0];

  /* Convert the box size if we want to clean-up h-factors */
  if (cleanup_h) {
    dim[0] /= h;
    dim[1] /= h;
    dim[2] /= h;
  }

  /* message("Found %d particles in a %speriodic box of size [%f %f %f].",  */
  /* 	  *N, (periodic ? "": "non-"), dim[0], dim[1], dim[2]);  */

  /* Close header */
  H5Gclose(h_grp);

  /* Read the unit system used in the ICs */
  struct unit_system* ic_units =
      (struct unit_system*)malloc(sizeof(struct unit_system));
  if (ic_units == NULL) error("Unable to allocate memory for IC unit system");
  io_read_unit_system(h_file, ic_units, internal_units, 0);

  /* Tell the user if a conversion will be needed */
  if (units_are_equal(ic_units, internal_units)) {

    message("IC and internal units match. No conversion needed.");

  } else {

    message("Conversion needed from:");
    message("(ICs) Unit system: U_M =      %e g.", ic_units->UnitMass_in_cgs);
    message("(ICs) Unit system: U_L =      %e cm.",
            ic_units->UnitLength_in_cgs);
    message("(ICs) Unit system: U_t =      %e s.", ic_units->UnitTime_in_cgs);
    message("(ICs) Unit system: U_I =      %e A.",
            ic_units->UnitCurrent_in_cgs);
    message("(ICs) Unit system: U_T =      %e K.",
            ic_units->UnitTemperature_in_cgs);
    message("to:");
    message("(internal) Unit system: U_M = %e g.",
            internal_units->UnitMass_in_cgs);
    message("(internal) Unit system: U_L = %e cm.",
            internal_units->UnitLength_in_cgs);
    message("(internal) Unit system: U_t = %e s.",
            internal_units->UnitTime_in_cgs);
    message("(internal) Unit system: U_I = %e A.",
            internal_units->UnitCurrent_in_cgs);
    message("(internal) Unit system: U_T = %e K.",
            internal_units->UnitTemperature_in_cgs);
  }

  /* Read metadata from ICs file */
  ic_info_read_hdf5(ics_metadata, h_file);

  /* Convert the dimensions of the box */
  for (int j = 0; j < 3; j++)
    dim[j] *=
        units_conversion_factor(ic_units, internal_units, UNIT_CONV_LENGTH);

  /* Allocate memory to store SPH particles */
  if (with_hydro) {
    *Ngas = N[swift_type_gas];
    if (swift_memalign("parts", (void**)parts, part_align,
                       *Ngas * sizeof(struct part)) != 0)
      error("Error while allocating memory for SPH particles");
    bzero(*parts, *Ngas * sizeof(struct part));
  }

  /* Allocate memory to store star particles */
  if (with_sink) {
    *Nsinks = N[swift_type_sink];
    if (swift_memalign("sinks", (void**)sinks, sink_align,
                       *Nsinks * sizeof(struct sink)) != 0)
      error("Error while allocating memory for sink particles");
    bzero(*sinks, *Nsinks * sizeof(struct sink));
  }

  /* Allocate memory to store star particles */
  if (with_stars) {
    *Nstars = N[swift_type_stars];
    if (swift_memalign("sparts", (void**)sparts, spart_align,
                       *Nstars * sizeof(struct spart)) != 0)
      error("Error while allocating memory for stars particles");
    bzero(*sparts, *Nstars * sizeof(struct spart));
  }

  /* Allocate memory to store black hole particles */
  if (with_black_holes) {
    *Nblackholes = N[swift_type_black_hole];
    if (swift_memalign("bparts", (void**)bparts, bpart_align,
                       *Nblackholes * sizeof(struct bpart)) != 0)
      error("Error while allocating memory for black hole particles");
    bzero(*bparts, *Nblackholes * sizeof(struct bpart));
  }

  /* Allocate memory to store all gravity particles */
  if (with_gravity) {
    Ndm = N[swift_type_dark_matter];
    Ndm_background = N[swift_type_dark_matter_background];
    Ndm_neutrino = N[swift_type_neutrino];
    *Ngparts = (with_hydro ? N[swift_type_gas] : 0) +
               N[swift_type_dark_matter] +
               N[swift_type_dark_matter_background] + N[swift_type_neutrino] +
               (with_sink ? N[swift_type_sink] : 0) +
               (with_stars ? N[swift_type_stars] : 0) +
               (with_black_holes ? N[swift_type_black_hole] : 0);
    *Ngparts_background = Ndm_background;
    *Nnuparts = Ndm_neutrino;
    if (swift_memalign("gparts", (void**)gparts, gpart_align,
                       *Ngparts * sizeof(struct gpart)) != 0)
      error("Error while allocating memory for gravity particles");
    bzero(*gparts, *Ngparts * sizeof(struct gpart));
  }

  /* message("Allocated %8.2f MB for particles.", *N * sizeof(struct part) /
   * (1024.*1024.)); */

  /* message("BoxSize = %lf", dim[0]); */
  /* message("NumPart = [%zd, %zd] Total = %zd", *Ngas, Ndm, *Ngparts); */

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gopen(h_file, partTypeGroupName, H5P_DEFAULT);
    if (h_grp < 0)
      error("Error while opening particle group %s.", partTypeGroupName);

    int num_fields = 0;
    struct io_props list[io_max_size_output_list];
    bzero(list, io_max_size_output_list * sizeof(struct io_props));
    size_t Nparticles = 0;

    /* Read particle fields into the structure */
    switch (ptype) {

      case swift_type_gas:
        if (with_hydro) {
          Nparticles = *Ngas;
          hydro_read_particles(*parts, list, &num_fields);
          num_fields += mhd_read_particles(*parts, list + num_fields);
          num_fields += chemistry_read_particles(*parts, list + num_fields);
          num_fields += rt_read_particles(*parts, list + num_fields);
        }
        break;

      case swift_type_dark_matter:
        if (with_gravity) {
          Nparticles = Ndm;
          darkmatter_read_particles(*gparts, list, &num_fields);
        }
        break;

      case swift_type_dark_matter_background:
        if (with_gravity) {
          Nparticles = Ndm_background;
          darkmatter_read_particles(*gparts + Ndm, list, &num_fields);
        }
        break;

      case swift_type_neutrino:
        if (with_gravity) {
          Nparticles = Ndm_neutrino;
          darkmatter_read_particles(*gparts + Ndm + Ndm_background, list,
                                    &num_fields);
        }
        break;

      case swift_type_sink:
        if (with_sink) {
          Nparticles = *Nsinks;
          sink_read_particles(*sinks, list, &num_fields);
        }
        break;

      case swift_type_stars:
        if (with_stars) {
          Nparticles = *Nstars;
          stars_read_particles(*sparts, list, &num_fields);
          num_fields +=
              star_formation_read_particles(*sparts, list + num_fields);
          num_fields += rt_read_stars(*sparts, list + num_fields);
        }
        break;

      case swift_type_black_hole:
        if (with_black_holes) {
          Nparticles = *Nblackholes;
          black_holes_read_particles(*bparts, list, &num_fields);
        }
        break;

      default:
        message("Particle Type %d not yet supported. Particles ignored", ptype);
    }

    /* Read everything */
    if (!dry_run)
      for (int i = 0; i < num_fields; ++i) {
        /* If we are remapping ParticleIDs later, don't need to read them. */
        if (remap_ids && strcmp(list[i].name, "ParticleIDs") == 0) continue;

        /* Read array. */
        read_array_single(h_grp, list[i], Nparticles, internal_units, ic_units,
                          cleanup_h, cleanup_sqrt_a, h, a);
      }

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* If we are remapping ParticleIDs later, start by setting them to 1. */
  if (remap_ids) io_set_ids_to_one(*gparts, *Ngparts);

  /* Duplicate the parts for gravity */
  if (!dry_run && with_gravity) {

    /* Let's initialise a bit of thread parallelism here */
    struct threadpool tp;
    threadpool_init(&tp, n_threads);

    /* Prepare the DM particles */
    io_prepare_dm_gparts(&tp, *gparts, Ndm);

    /* Prepare the DM background particles */
    io_prepare_dm_background_gparts(&tp, *gparts + Ndm, Ndm_background);

    /* Prepare the DM neutrino particles */
    io_prepare_dm_neutrino_gparts(&tp, *gparts + Ndm + Ndm_background,
                                  Ndm_neutrino);

    /* Duplicate the hydro particles into gparts */
    if (with_hydro)
      io_duplicate_hydro_gparts(&tp, *parts, *gparts, *Ngas,
                                Ndm + Ndm_background + Ndm_neutrino);

    /* Duplicate the star particles into gparts */
    if (with_sink)
      io_duplicate_sinks_gparts(&tp, *sinks, *gparts, *Nsinks,
                                Ndm + Ndm_background + Ndm_neutrino + *Ngas);

    /* Duplicate the star particles into gparts */
    if (with_stars)
      io_duplicate_stars_gparts(
          &tp, *sparts, *gparts, *Nstars,
          Ndm + Ndm_background + Ndm_neutrino + *Ngas + *Nsinks);

    /* Duplicate the black hole particles into gparts */
    if (with_black_holes)
      io_duplicate_black_holes_gparts(
          &tp, *bparts, *gparts, *Nblackholes,
          Ndm + Ndm_background + Ndm_neutrino + *Ngas + *Nsinks + *Nstars);

    threadpool_clean(&tp);
  }

  /* message("Done Reading particles..."); */

  /* Clean up */
  free(ic_units);

  /* Close file */
  H5Fclose(h_file);
}

/**
 * @brief Writes an HDF5 output file (GADGET-3 type) with its XMF descriptor
 *
 * @param e The engine containing all the system.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param fof Is this a snapshot related to a stand-alone FOF call?
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_single(struct engine* e,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units,
                         const int fof) {

  hid_t h_file = 0, h_grp = 0;
  int numFiles = 1;
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  const struct bpart* bparts = e->s->bparts;
  const struct sink* sinks = e->s->sinks;
  struct output_options* output_options = e->output_options;
  struct output_list* output_list = e->output_list_snapshots;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
  const int with_fof = e->policy & engine_policy_fof;
  const int with_DM_background = e->s->with_DM_background;
  const int with_DM = e->s->with_DM;
  const int with_neutrinos = e->s->with_neutrinos;
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
  const int with_rt = e->policy & engine_policy_rt;

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Nsinks = e->s->nr_sinks;
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

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
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

  /* Create the directory */
  safe_checkdir(snapshot_subdir_name, /*create=*/1);

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

  /* First time, we need to create the XMF file */
  if (e->snapshot_output_count == 0) xmf_create_file(xmfFileName);

  /* Prepare the XMF file for the new entry */
  FILE* xmfFile = 0;
  xmfFile = xmf_prepare_file(xmfFileName);

  /* Write the part corresponding to this specific output */
  xmf_write_outputheader(xmfFile, fileName, e->time);

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

  /* Format things in a Gadget-friendly array */
  const long long N_total[swift_type_count] = {
      (long long)Ngas_written,   (long long)Ndm_written,
      (long long)Ndm_background, (long long)Nsinks_written,
      (long long)Nstars_written, (long long)Nblackholes_written,
      (long long)Ndm_neutrino};

  /* List what fields to write.
   * Note that we want to want to write a 0-size dataset for some species
   * in case future snapshots will contain them (e.g. star formation) */
  const int to_write[swift_type_count] = {
      with_hydro, with_DM,         with_DM_background, with_sink,
      with_stars, with_black_hole, with_neutrinos

  };

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

  /* GADGET-2 legacy values: number of particles of each type */
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

  /* Now write the top-level cell structure */
  long long global_offsets[swift_type_count] = {0};
  h_grp = H5Gcreate(h_file, "/Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cells group");

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp, e->s->cdim, e->s->dim, e->s->cells_top,
                        e->s->nr_cells, e->s->width, e->nodeID,
                        /*distributed=*/0, subsample, subsample_fraction,
                        e->snapshot_output_count, N_total, global_offsets,
                        to_write, numFields, internal_units, snapshot_units);
  H5Gclose(h_grp);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if there are
     * (a) no particles of this kind in this run, or
     * (b) if we have disabled every field of this particle type. */
    if (!to_write[ptype] || numFields[ptype] == 0) continue;

    /* Add the global information for that particle type to the XMF meta-file */
    xmf_write_groupheader(xmfFile, fileName, /*distributed=*/0,
                          numParticles[ptype], (enum part_type)ptype);

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

    /* Write the number of particles and fields as an attribute */
    io_write_attribute_ll(h_grp, "NumberOfParticles", N_total[ptype]);
    io_write_attribute_ll(h_grp, "TotalNumberOfParticles", N_total[ptype]);

    int num_fields = 0;
    struct io_props list[io_max_size_output_list];
    bzero(list, io_max_size_output_list * sizeof(struct io_props));
    size_t N = 0;

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;
    struct sink* sinks_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas: {
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          N = Ngas;

          /* Select the fields to write */
          io_select_hydro_fields(parts, xparts, with_cosmology, with_cooling,
                                 with_temperature, with_fof, with_stf, with_rt,
                                 e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ngas_written;

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

          /* This is a DM-only run without background or inhibited particles or
           * neutrinos */
          N = Ntot;

          /* Select the fields to write */
          io_select_dm_fields(gparts, e->s->gpart_group_data, with_fof,
                              with_stf, e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ndm_written;

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
        N = Ndm_background;

        /* Allocate temporary array */
        if (swift_memalign("gparts_written", (void**)&gparts_written,
                           gpart_align,
                           Ndm_background * sizeof(struct gpart)) != 0)
          error("Error while allocating temporary memory for gparts");

        if (with_stf) {
          if (swift_memalign(
                  "gpart_group_written", (void**)&gpart_group_data_written,
                  gpart_align,
                  Ndm_background * sizeof(struct velociraptor_gpart_data)) != 0)
            error(
                "Error while allocating temporary memory for gparts STF "
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
        N = Ndm_neutrino;

        /* Allocate temporary array */
        if (swift_memalign("gparts_written", (void**)&gparts_written,
                           gpart_align,
                           Ndm_neutrino * sizeof(struct gpart)) != 0)
          error("Error while allocating temporary memory for gparts");

        if (with_stf) {
          if (swift_memalign(
                  "gpart_group_written", (void**)&gpart_group_data_written,
                  gpart_align,
                  Ndm_neutrino * sizeof(struct velociraptor_gpart_data)) != 0)
            error(
                "Error while allocating temporary memory for gparts STF "
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
          N = Nsinks;

          /* Select the fields to write */
          io_select_sink_fields(sinks, with_cosmology, with_fof, with_stf, e,
                                &num_fields, list);
        } else {

          /* Ok, we need to fish out the particles we want */
          N = Nsinks_written;

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

          /* No inhibited particles: easy case */
          N = Nstars;

          /* Select the fields to write */
          io_select_star_fields(sparts, with_cosmology, with_fof, with_stf,
                                with_rt, e, &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Nstars_written;

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

          /* No inhibited particles: easy case */
          N = Nblackholes;

          /* Select the fields to write */
          io_select_bh_fields(bparts, with_cosmology, with_fof, with_stf, e,
                              &num_fields, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Nblackholes_written;

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
        write_array_single(e, h_grp, fileName, xmfFile, partTypeGroupName,
                           list[i], N, compression_level, internal_units,
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
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);
    if (sinks_written) swift_free("sinks_written", sinks_written);

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
  } /* ends loop over particle types */

  /* Write LXMF file descriptor */
  xmf_write_outputfooter(xmfFile, e->snapshot_output_count, e->time);

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  e->snapshot_output_count++;
  if (e->snapshot_invoke_stf) e->stf_output_count++;
}

#endif /* HAVE_HDF5 && !WITH_MPI */
