/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#include "chemistry_io.h"
#include "common_io.h"
#include "cooling_io.h"
#include "dimension.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_properties.h"
#include "kernel_hydro.h"
#include "part.h"
#include "part_type.h"
#include "star_formation_io.h"
#include "stars_io.h"
#include "tracers_io.h"
#include "units.h"
#include "velociraptor_io.h"
#include "xmf.h"

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
void readArray(hid_t h_grp, const struct io_props props, size_t N,
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
      /* message("Optional data set '%s' not present. Zeroing this particle
       * props...", name);	   */

      for (size_t i = 0; i < N; ++i)
        memset(props.field + i * props.partSize, 0, copySize);

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
void writeArray(const struct engine* e, hid_t grp, char* fileName,
                FILE* xmfFile, char* partTypeGroupName,
                const struct io_props props, size_t N,
                const struct unit_system* internal_units,
                const struct unit_system* snapshot_units) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t num_elements = N * props.dimension;

  /* message("Writing '%s' array...", props.name); */

  /* Allocate temporary buffer */
  void* temp = NULL;
  if (posix_memalign((void**)&temp, IO_BUFFER_ALIGNMENT,
                     num_elements * typeSize) != 0)
    error("Unable to allocate temporary i/o buffer");

  /* Copy the particle data to the temporary buffer */
  io_copy_temp_buffer(temp, e, props, N, internal_units, snapshot_units);

  /* Create data space */
  const hid_t h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0)
    error("Error while creating data space for field '%s'.", props.name);

  int rank;
  hsize_t shape[2];
  hsize_t chunk_shape[2];

  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = props.dimension;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 20; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  hid_t h_err = H5Sset_extent_simple(h_space, rank, shape, shape);
  if (h_err < 0)
    error("Error while changing data space shape for field '%s'.", props.name);

  /* Dataset properties */
  const hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Set chunk size */
  h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
  if (h_err < 0)
    error("Error while setting chunk size (%llu, %llu) for field '%s'.",
          chunk_shape[0], chunk_shape[1], props.name);

  /* Impose check-sum to verify data corruption */
  h_err = H5Pset_fletcher32(h_prop);
  if (h_err < 0)
    error("Error while setting checksum options for field '%s'.", props.name);

  /* Impose data compression */
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

  /* Create dataset */
  const hid_t h_data = H5Dcreate(grp, props.name, io_hdf5_type(props.type),
                                 h_space, H5P_DEFAULT, h_prop, H5P_DEFAULT);
  if (h_data < 0) error("Error while creating dataspace '%s'.", props.name);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_space, H5S_ALL,
                   H5P_DEFAULT, temp);
  if (h_err < 0) error("Error while writing data array '%s'.", props.name);

  /* Write XMF description for this data set */
  if (xmfFile != NULL)
    xmf_write_line(xmfFile, fileName, partTypeGroupName, props.name, N,
                   props.dimension, props.type);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE];
  units_cgs_conversion_string(buffer, snapshot_units, props.units);
  io_write_attribute_d(
      h_data, "CGS conversion factor",
      units_cgs_conversion_factor(snapshot_units, props.units));
  io_write_attribute_f(h_data, "h-scale exponent", 0);
  io_write_attribute_f(h_data, "a-scale exponent",
                       units_a_factor(snapshot_units, props.units));
  io_write_attribute_s(h_data, "Conversion factor", buffer);

  /* Free and close everything */
  free(temp);
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
 * @param sparts (output) Array of #spart particles.
 * @param Ngas (output) number of Gas particles read.
 * @param Ngparts (output) The number of #gpart read.
 * @param Nstars (output) The number of #spart read.
 * @param flag_entropy (output) 1 if the ICs contained Entropy in the
 * InternalEnergy field
 * @param with_hydro Are we reading gas particles ?
 * @param with_gravity Are we reading/creating #gpart arrays ?
 * @param with_stars Are we reading star particles ?
 * @param cleanup_h Are we cleaning-up h-factors from the quantities we read?
 * @param cleanup_sqrt_a Are we cleaning-up the sqrt(a) factors in the Gadget
 * IC velocities?
 * @param h The value of the reduced Hubble constant to use for correction.
 * @param a The current value of the scale-factor.
 * @prarm n_threads The number of threads to use for the temporary threadpool.
 * @param dry_run If 1, don't read the particle. Only allocates the arrays.
 *
 * Opens the HDF5 file fileName and reads the particles contained
 * in the parts array. N is the returned number of particles found
 * in the file.
 *
 * @warning Can not read snapshot distributed over more than 1 file !!!
 * @todo Read snapshots distributed in more than one file.
 */
void read_ic_single(const char* fileName,
                    const struct unit_system* internal_units, double dim[3],
                    struct part** parts, struct gpart** gparts,
                    struct spart** sparts, size_t* Ngas, size_t* Ngparts,
                    size_t* Nstars, int* flag_entropy, int with_hydro,
                    int with_gravity, int with_stars, int cleanup_h,
                    int cleanup_sqrt_a, double h, double a, int n_threads,
                    int dry_run) {

  hid_t h_file = 0, h_grp = 0;
  /* GADGET has only cubic boxes (in cosmological mode) */
  double boxSize[3] = {0.0, -1.0, -1.0};
  /* GADGET has 6 particle types. We only keep the type 0 & 1 for now...*/
  long long numParticles[swift_type_count] = {0};
  long long numParticles_highWord[swift_type_count] = {0};
  size_t N[swift_type_count] = {0};
  int dimension = 3; /* Assume 3D if nothing is specified */
  size_t Ndm = 0;

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
  int flag_entropy_temp[6];
  io_read_attribute(h_grp, "Flag_Entropy_ICs", INT, flag_entropy_temp);
  *flag_entropy = flag_entropy_temp[0];
  io_read_attribute(h_grp, "BoxSize", DOUBLE, boxSize);
  io_read_attribute(h_grp, "NumPart_Total", LONGLONG, numParticles);
  io_read_attribute(h_grp, "NumPart_Total_HighWord", LONGLONG,
                    numParticles_highWord);

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

  /* Convert the dimensions of the box */
  for (int j = 0; j < 3; j++)
    dim[j] *=
        units_conversion_factor(ic_units, internal_units, UNIT_CONV_LENGTH);

  /* Allocate memory to store SPH particles */
  if (with_hydro) {
    *Ngas = N[swift_type_gas];
    if (posix_memalign((void**)parts, part_align,
                       *Ngas * sizeof(struct part)) != 0)
      error("Error while allocating memory for SPH particles");
    bzero(*parts, *Ngas * sizeof(struct part));
  }

  /* Allocate memory to store star particles */
  if (with_stars) {
    *Nstars = N[swift_type_stars];
    if (posix_memalign((void**)sparts, spart_align,
                       *Nstars * sizeof(struct spart)) != 0)
      error("Error while allocating memory for stars particles");
    bzero(*sparts, *Nstars * sizeof(struct spart));
  }

  /* Allocate memory to store all gravity particles */
  if (with_gravity) {
    Ndm = N[swift_type_dark_matter];
    *Ngparts = (with_hydro ? N[swift_type_gas] : 0) +
               N[swift_type_dark_matter] +
               (with_stars ? N[swift_type_stars] : 0);
    if (posix_memalign((void**)gparts, gpart_align,
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
    struct io_props list[100];
    size_t Nparticles = 0;

    /* Read particle fields into the structure */
    switch (ptype) {

      case swift_type_gas:
        if (with_hydro) {
          Nparticles = *Ngas;
          hydro_read_particles(*parts, list, &num_fields);
          num_fields += chemistry_read_particles(*parts, list + num_fields);
        }
        break;

      case swift_type_dark_matter:
        if (with_gravity) {
          Nparticles = Ndm;
          darkmatter_read_particles(*gparts, list, &num_fields);
        }
        break;

      case swift_type_stars:
        if (with_stars) {
          Nparticles = *Nstars;
          stars_read_particles(*sparts, list, &num_fields);
        }
        break;

      default:
        message("Particle Type %d not yet supported. Particles ignored", ptype);
    }

    /* Read everything */
    if (!dry_run)
      for (int i = 0; i < num_fields; ++i)
        readArray(h_grp, list[i], Nparticles, internal_units, ic_units,
                  cleanup_h, cleanup_sqrt_a, h, a);

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* Duplicate the parts for gravity */
  if (!dry_run && with_gravity) {

    /* Let's initialise a bit of thread parallelism here */
    struct threadpool tp;
    threadpool_init(&tp, n_threads);

    /* Prepare the DM particles */
    io_prepare_dm_gparts(&tp, *gparts, Ndm);

    /* Duplicate the hydro particles into gparts */
    if (with_hydro) io_duplicate_hydro_gparts(&tp, *parts, *gparts, *Ngas, Ndm);

    /* Duplicate the star particles into gparts */
    if (with_stars)
      io_duplicate_stars_gparts(&tp, *sparts, *gparts, *Nstars, Ndm + *Ngas);

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
 * @param baseName The common part of the snapshot file name.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_single(struct engine* e, const char* baseName,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units) {

  hid_t h_file = 0, h_grp = 0;
  int numFiles = 1;
  const struct part* parts = e->s->parts;
  const struct xpart* xparts = e->s->xparts;
  const struct gpart* gparts = e->s->gparts;
  const struct spart* sparts = e->s->sparts;
  struct swift_params* params = e->parameter_file;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int with_cooling = e->policy & engine_policy_cooling;
  const int with_temperature = e->policy & engine_policy_temperature;
#ifdef HAVE_VELOCIRAPTOR
  const int with_stf = (e->policy & engine_policy_structure_finding) &&
                       (e->s->gpart_group_data != NULL);
#else
  const int with_stf = 0;
#endif

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  // const size_t Nbaryons = Ngas + Nstars;
  // const size_t Ndm = Ntot > 0 ? Ntot - Nbaryons : 0;

  /* Number of particles that we will write */
  const size_t Ntot_written =
      e->s->nr_gparts - e->s->nr_inhibited_gparts - e->s->nr_extra_gparts;
  const size_t Ngas_written =
      e->s->nr_parts - e->s->nr_inhibited_parts - e->s->nr_extra_parts;
  const size_t Nstars_written =
      e->s->nr_sparts - e->s->nr_inhibited_sparts - e->s->nr_extra_sparts;
  const size_t Nbaryons_written = Ngas_written + Nstars_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written : 0;

  /* Format things in a Gadget-friendly array */
  long long N_total[swift_type_count] = {(long long)Ngas_written,
                                         (long long)Ndm_written,
                                         0,
                                         0,
                                         (long long)Nstars_written,
                                         0};

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  if (e->snapshot_int_time_label_on)
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%06i.hdf5", baseName,
             (int)round(e->time));
  else
    snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.hdf5", baseName,
             e->snapshot_output_count);

  /* First time, we need to create the XMF file */
  if (e->snapshot_output_count == 0) xmf_create_file(baseName);

  /* Prepare the XMF file for the new entry */
  FILE* xmfFile = 0;
  xmfFile = xmf_prepare_file(baseName);

  /* Write the part corresponding to this specific output */
  xmf_write_outputheader(xmfFile, fileName, e->time);

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
  time_t tm = time(NULL);
  io_write_attribute_s(h_grp, "Snapshot date", ctime(&tm));

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
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
  double MassTable[swift_type_count] = {0};
  io_write_attribute(h_grp, "MassTable", DOUBLE, MassTable, swift_type_count);
  unsigned int flagEntropy[swift_type_count] = {0};
  flagEntropy[0] = writeEntropyFlag();
  io_write_attribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                     swift_type_count);
  io_write_attribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

  /* Close header */
  H5Gclose(h_grp);

  /* Print the code version */
  io_write_code_description(h_file);

  /* Print the run's policy */
  io_write_engine_policy(h_file, e);

  /* Print the SPH parameters */
  if (e->policy & engine_policy_hydro) {
    h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating SPH group");
    hydro_props_print_snapshot(h_grp, e->hydro_properties);
    hydro_write_flavour(h_grp);
    H5Gclose(h_grp);
  }

  /* Print the subgrid parameters */
  h_grp = H5Gcreate(h_file, "/SubgridScheme", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating subgrid group");
  entropy_floor_write_flavour(h_grp);
  cooling_write_flavour(h_grp, e->cooling_func);
  chemistry_write_flavour(h_grp);
  tracers_write_flavour(h_grp);
  H5Gclose(h_grp);

  /* Print the gravity parameters */
  if (e->policy & engine_policy_self_gravity) {
    h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating gravity group");
    gravity_props_print_snapshot(h_grp, e->gravity_properties);
    H5Gclose(h_grp);
  }

  /* Print the stellar parameters */
  if (e->policy & engine_policy_stars) {
    h_grp = H5Gcreate(h_file, "/StarsScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating stars group");
    stars_props_print_snapshot(h_grp, e->stars_properties);
    H5Gclose(h_grp);
  }

  /* Print the cosmological model  */
  h_grp =
      H5Gcreate(h_file, "/Cosmology", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cosmology group");
  if (e->policy & engine_policy_cosmology)
    io_write_attribute_i(h_grp, "Cosmological run", 1);
  else
    io_write_attribute_i(h_grp, "Cosmological run", 0);
  cosmology_write_model(h_grp, e->cosmology);
  H5Gclose(h_grp);

  /* Print the runtime parameters */
  h_grp =
      H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 1);
  H5Gclose(h_grp);

  /* Print the runtime unused parameters */
  h_grp = H5Gcreate(h_file, "/UnusedParameters", H5P_DEFAULT, H5P_DEFAULT,
                    H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp, 0);
  H5Gclose(h_grp);

  /* Print the system of Units used in the spashot */
  io_write_unit_system(h_file, snapshot_units, "Units");

  /* Print the system of Units used internally */
  io_write_unit_system(h_file, internal_units, "InternalCodeUnits");

  /* Now write the top-level cell structure */
  long long global_offsets[swift_type_count] = {0};
  h_grp = H5Gcreate(h_file, "/Cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating cells group");

  /* Write the location of the particles in the arrays */
  io_write_cell_offsets(h_grp, e->s->cdim, e->s->cells_top, e->s->nr_cells,
                        e->s->width, e->nodeID, N_total, global_offsets,
                        internal_units, snapshot_units);
  H5Gclose(h_grp);

  /* Tell the user if a conversion will be needed */
  if (e->verbose) {
    if (units_are_equal(snapshot_units, internal_units)) {

      message("Snapshot and internal units match. No conversion needed.");

    } else {

      message("Conversion needed from:");
      message("(Snapshot) Unit system: U_M =      %e g.",
              snapshot_units->UnitMass_in_cgs);
      message("(Snapshot) Unit system: U_L =      %e cm.",
              snapshot_units->UnitLength_in_cgs);
      message("(Snapshot) Unit system: U_t =      %e s.",
              snapshot_units->UnitTime_in_cgs);
      message("(Snapshot) Unit system: U_I =      %e A.",
              snapshot_units->UnitCurrent_in_cgs);
      message("(Snapshot) Unit system: U_T =      %e K.",
              snapshot_units->UnitTemperature_in_cgs);
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
  }

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (numParticles[ptype] == 0) continue;

    /* Add the global information for that particle type to the XMF meta-file */
    xmf_write_groupheader(xmfFile, fileName, numParticles[ptype],
                          (enum part_type)ptype);

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating particle group.\n");

    int num_fields = 0;
    struct io_props list[100];
    size_t N = 0;

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct spart* sparts_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas: {
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          N = Ngas;
          hydro_write_particles(parts, xparts, list, &num_fields);
          num_fields += chemistry_write_particles(parts, list + num_fields);
          if (with_cooling || with_temperature) {
            num_fields += cooling_write_particles(
                parts, xparts, list + num_fields, e->cooling_func);
          }
          if (with_stf) {
            num_fields +=
                velociraptor_write_parts(parts, xparts, list + num_fields);
          }
          num_fields += tracers_write_particles(
              parts, xparts, list + num_fields, with_cosmology);
          num_fields +=
              star_formation_write_particles(parts, xparts, list + num_fields);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ngas_written;

          /* Allocate temporary arrays */
          if (posix_memalign((void**)&parts_written, part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporart memory for parts");
          if (posix_memalign((void**)&xparts_written, xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporart memory for xparts");

          /* Collect the particles we want to write */
          io_collect_parts_to_write(parts, xparts, parts_written,
                                    xparts_written, Ngas, Ngas_written);

          /* Select the fields to write */
          hydro_write_particles(parts_written, xparts_written, list,
                                &num_fields);
          num_fields +=
              chemistry_write_particles(parts_written, list + num_fields);
          if (with_cooling || with_temperature) {
            num_fields +=
                cooling_write_particles(parts_written, xparts_written,
                                        list + num_fields, e->cooling_func);
          }
          if (with_stf) {
            num_fields += velociraptor_write_parts(
                parts_written, xparts_written, list + num_fields);
          }
          num_fields += tracers_write_particles(
              parts_written, xparts_written, list + num_fields, with_cosmology);
          num_fields += star_formation_write_particles(
              parts_written, xparts_written, list + num_fields);
        }
      } break;

      case swift_type_dark_matter: {
        if (Ntot == Ndm_written) {

          /* This is a DM-only run without inhibited particles */
          N = Ntot;
          darkmatter_write_particles(gparts, list, &num_fields);
          if (with_stf) {
            num_fields += velociraptor_write_gparts(e->s->gpart_group_data,
                                                    list + num_fields);
          }
        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ndm_written;

          /* Allocate temporary array */
          if (posix_memalign((void**)&gparts_written, gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporart memory for gparts");

          if (with_stf) {
            if (posix_memalign(
                    (void**)&gpart_group_data_written, gpart_align,
                    Ndm_written * sizeof(struct velociraptor_gpart_data)) != 0)
              error(
                  "Error while allocating temporart memory for gparts STF "
                  "data");
          }

          /* Collect the non-inhibited DM particles from gpart */
          io_collect_gparts_to_write(gparts, e->s->gpart_group_data,
                                     gparts_written, gpart_group_data_written,
                                     Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          darkmatter_write_particles(gparts_written, list, &num_fields);
          if (with_stf) {
            num_fields += velociraptor_write_gparts(gpart_group_data_written,
                                                    list + num_fields);
          }
        }
      } break;

      case swift_type_stars: {
        if (Nstars == Nstars_written) {

          /* No inhibted particles: easy case */
          N = Nstars;
          stars_write_particles(sparts, list, &num_fields);
          num_fields += chemistry_write_sparticles(sparts, list + num_fields);
          num_fields += tracers_write_sparticles(sparts, list + num_fields,
                                                 with_cosmology);
          if (with_stf) {
            num_fields += velociraptor_write_sparts(sparts, list + num_fields);
          }
        } else {

          /* Ok, we need to fish out the particles we want */
          N = Nstars_written;

          /* Allocate temporary arrays */
          if (posix_memalign((void**)&sparts_written, spart_align,
                             Nstars_written * sizeof(struct spart)) != 0)
            error("Error while allocating temporart memory for sparts");

          /* Collect the particles we want to write */
          io_collect_sparts_to_write(sparts, sparts_written, Nstars,
                                     Nstars_written);

          /* Select the fields to write */
          stars_write_particles(sparts_written, list, &num_fields);
          num_fields +=
              chemistry_write_sparticles(sparts_written, list + num_fields);
          num_fields += tracers_write_sparticles(
              sparts_written, list + num_fields, with_cosmology);
          if (with_stf) {
            num_fields +=
                velociraptor_write_sparts(sparts_written, list + num_fields);
          }
        }
      } break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Write everything that is not cancelled */
    for (int i = 0; i < num_fields; ++i) {

      /* Did the user cancel this field? */
      char field[PARSER_MAX_LINE_SIZE];
      sprintf(field, "SelectOutput:%s_%s", list[i].name,
              part_type_names[ptype]);
      int should_write = parser_get_opt_param_int(params, field, 1);

      if (should_write)
        writeArray(e, h_grp, fileName, xmfFile, partTypeGroupName, list[i], N,
                   internal_units, snapshot_units);
    }

    /* Free temporary arrays */
    if (parts_written) free(parts_written);
    if (xparts_written) free(xparts_written);
    if (gparts_written) free(gparts_written);
    if (gpart_group_data_written) free(gpart_group_data_written);
    if (sparts_written) free(sparts_written);

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
  }

  /* Write LXMF file descriptor */
  xmf_write_outputfooter(xmfFile, e->snapshot_output_count, e->time);

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  e->snapshot_output_count++;
}

#endif /* HAVE_HDF5 */
