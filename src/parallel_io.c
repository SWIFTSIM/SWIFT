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

#if defined(HAVE_HDF5) && defined(WITH_MPI) && defined(HAVE_PARALLEL_HDF5)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "parallel_io.h"

/* Local includes. */
#include "common_io.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_properties.h"
#include "kernel_hydro.h"
#include "part.h"
#include "stars_io.h"
#include "units.h"
#include "xmf.h"

/**
 * @brief Reads a data array from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the array to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param N The number of particles.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurrence of the field of
 *interest in the parts array
 * @param importance If COMPULSORY, the data must be present in the IC file. If
 *OPTIONAL, the array will be zeroed when the data is not present.
 *
 * @todo A better version using HDF5 hyper-slabs to read the file directly into
 *the part array
 * will be written once the structures have been stabilized.
 *
 * Calls #error() if an error occurs.
 */
void readArray(hid_t grp, const struct io_props prop, size_t N,
               long long N_total, long long offset,
               const struct unit_system* internal_units,
               const struct unit_system* ic_units) {

  const size_t typeSize = io_sizeof_type(prop.type);
  const size_t copySize = typeSize * prop.dimension;
  const size_t num_elements = N * prop.dimension;

  /* Check whether the dataspace exists or not */
  const htri_t exist = H5Lexists(grp, prop.name, 0);
  if (exist < 0) {
    error("Error while checking the existence of data set '%s'.", prop.name);
  } else if (exist == 0) {
    if (prop.importance == COMPULSORY) {
      error("Compulsory data set '%s' not present in the file.", prop.name);
    } else {
      for (size_t i = 0; i < N; ++i)
        memset(prop.field + i * prop.partSize, 0, copySize);
      return;
    }
  }

  /* message("Reading %s '%s' array...", */
  /*         prop.importance == COMPULSORY ? "compulsory" : "optional  ", */
  /*         prop.name); */

  /* Open data space in file */
  const hid_t h_data = H5Dopen2(grp, prop.name, H5P_DEFAULT);
  if (h_data < 0) error("Error while opening data space '%s'.", prop.name);

  /* Check data type */
  const hid_t h_type = H5Dget_type(h_data);
  if (h_type < 0) error("Unable to retrieve data type from the file");
  /* if (!H5Tequal(h_type, hdf5_type(type))) */
  /*   error("Non-matching types between the code and the file"); */

  /* Allocate temporary buffer */
  void* temp = malloc(num_elements * typeSize);
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Prepare information for hyper-slab */
  hsize_t shape[2], offsets[2];
  int rank;
  if (prop.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = prop.dimension;
    offsets[0] = offset;
    offsets[1] = 0;
  } else {
    rank = 2;
    shape[0] = N;
    shape[1] = 1;
    offsets[0] = offset;
    offsets[1] = 0;
  }

  /* Create data space in memory */
  const hid_t h_memspace = H5Screate_simple(rank, shape, NULL);

  /* Select hyper-slab in file */
  const hid_t h_filespace = H5Dget_space(h_data);
  H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);

  /* Set collective reading properties */
  const hid_t h_plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE);

  /* Read HDF5 dataspace in temporary buffer */
  /* Dirty version that happens to work for vectors but should be improved */
  /* Using HDF5 dataspaces would be better */
  const hid_t h_err = H5Dread(h_data, io_hdf5_type(prop.type), h_memspace,
                              h_filespace, h_plist_id, temp);
  if (h_err < 0) {
    error("Error while reading data array '%s'.", prop.name);
  }

  /* Unit conversion if necessary */
  const double factor =
      units_conversion_factor(ic_units, internal_units, prop.units);
  if (factor != 1. && exist != 0) {

    /* message("Converting ! factor=%e", factor); */

    if (io_is_double_precision(prop.type)) {
      double* temp_d = temp;
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= factor;
    } else {
      float* temp_f = temp;
      for (size_t i = 0; i < num_elements; ++i) temp_f[i] *= factor;
    }
  }

  /* Copy temporary buffer to particle data */
  char* temp_c = temp;
  for (size_t i = 0; i < N; ++i)
    memcpy(prop.field + i * prop.partSize, &temp_c[i * copySize], copySize);

  /* Free and close everything */
  free(temp);
  H5Pclose(h_plist_id);
  H5Sclose(h_filespace);
  H5Sclose(h_memspace);
  H5Tclose(h_type);
  H5Dclose(h_data);
}

/*-----------------------------------------------------------------------------
 * Routines writing an output file
 *-----------------------------------------------------------------------------*/

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param e The #engine we are writing from.
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param N The number of particles to write.
 * @param N_total Total number of particles across all cores
 * @param offset Offset in the array where this mpi task starts writing
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 * the part array will be written once the structures have been stabilized.
 *
 */
void writeArray(struct engine* e, hid_t grp, char* fileName, FILE* xmfFile,
                char* partTypeGroupName, const struct io_props props, size_t N,
                long long N_total, int mpi_rank, long long offset,
                const struct unit_system* internal_units,
                const struct unit_system* snapshot_units) {

  const size_t typeSize = io_sizeof_type(props.type);
  const size_t copySize = typeSize * props.dimension;
  const size_t num_elements = N * props.dimension;

  /* message("Writing '%s' array...", props.name); */

  /* Allocate temporary buffer */
  void* temp = malloc(num_elements * io_sizeof_type(props.type));
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Copy particle data to temporary buffer */
  if (props.convert_part == NULL &&
      props.convert_gpart == NULL) { /* No conversion */

    char* temp_c = temp;
    for (size_t i = 0; i < N; ++i)
      memcpy(&temp_c[i * copySize], props.field + i * props.partSize, copySize);

  } else if (props.convert_part != NULL) { /* conversion (for parts)*/

    float* temp_f = temp;
    for (size_t i = 0; i < N; ++i)
      temp_f[i] = props.convert_part(e, &props.parts[i]);

  } else if (props.convert_gpart != NULL) { /* conversion (for gparts)*/

    float* temp_f = temp;
    for (size_t i = 0; i < N; ++i)
      temp_f[i] = props.convert_gpart(e, &props.gparts[i]);
  }

  /* Unit conversion if necessary */
  const double factor =
      units_conversion_factor(internal_units, snapshot_units, props.units);
  if (factor != 1.) {

    /* message("Converting ! factor=%e", factor); */

    if (io_is_double_precision(props.type)) {
      double* temp_d = temp;
      for (size_t i = 0; i < num_elements; ++i) temp_d[i] *= factor;
    } else {
      float* temp_f = temp;
      for (size_t i = 0; i < num_elements; ++i) temp_f[i] *= factor;
    }
  }

  /* Create data space */
  const hid_t h_memspace = H5Screate(H5S_SIMPLE);
  if (h_memspace < 0) {
    error("Error while creating data space (memory) for field '%s'.",
          props.name);
  }

  hid_t h_filespace = H5Screate(H5S_SIMPLE);
  if (h_filespace < 0) {
    error("Error while creating data space (file) for field '%s'.", props.name);
  }

  int rank;
  hsize_t shape[2];
  hsize_t shape_total[2];
  hsize_t offsets[2];
  if (props.dimension > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = props.dimension;
    shape_total[0] = N_total;
    shape_total[1] = props.dimension;
    offsets[0] = offset;
    offsets[1] = 0;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    shape_total[0] = N_total;
    shape_total[1] = 0;
    offsets[0] = offset;
    offsets[1] = 0;
  }

  /* Change shape of memory data space */
  hid_t h_err = H5Sset_extent_simple(h_memspace, rank, shape, NULL);
  if (h_err < 0) {
    error("Error while changing data space (memory) shape for field '%s'.",
          props.name);
  }

  /* Change shape of file data space */
  h_err = H5Sset_extent_simple(h_filespace, rank, shape_total, NULL);
  if (h_err < 0) {
    error("Error while changing data space (file) shape for field '%s'.",
          props.name);
  }

  /* Create dataset */
  const hid_t h_data =
      H5Dcreate(grp, props.name, io_hdf5_type(props.type), h_filespace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_data < 0) {
    error("Error while creating dataset '%s'.", props.name);
  }

  H5Sclose(h_filespace);
  h_filespace = H5Dget_space(h_data);
  H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);

  /* Create property list for collective dataset write.    */
  const hid_t h_plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, io_hdf5_type(props.type), h_memspace, h_filespace,
                   h_plist_id, temp);
  if (h_err < 0) {
    error("Error while writing data array '%s'.", props.name);
  }

  /* Write XMF description for this data set */
  if (mpi_rank == 0)
    xmf_write_line(xmfFile, fileName, partTypeGroupName, props.name, N_total,
                   props.dimension, props.type);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE];
  units_cgs_conversion_string(buffer, snapshot_units, props.units);
  io_write_attribute_d(
      h_data, "CGS conversion factor",
      units_cgs_conversion_factor(snapshot_units, props.units));
  io_write_attribute_f(h_data, "h-scale exponent",
                       units_h_factor(snapshot_units, props.units));
  io_write_attribute_f(h_data, "a-scale exponent",
                       units_a_factor(snapshot_units, props.units));
  io_write_attribute_s(h_data, "Conversion factor", buffer);

  /* Free and close everything */
  free(temp);
  H5Dclose(h_data);
  H5Pclose(h_plist_id);
  H5Sclose(h_memspace);
  H5Sclose(h_filespace);
}

/**
 * @brief Reads an HDF5 initial condition file (GADGET-3 type) in parallel
 *
 * @param fileName The file to read.
 * @param internal_units The system units used internally
 * @param dim (output) The dimension of the volume read from the file.
 * @param parts (output) The array of #part read from the file.
 * @param N (output) The number of particles read from the file.
 * @param periodic (output) 1 if the volume is periodic, 0 if not.
 * @param flag_entropy (output) 1 if the ICs contained Entropy in the
 * InternalEnergy field
 * @param mpi_rank The MPI rank of this node
 * @param mpi_size The number of MPI ranks
 * @param comm The MPI communicator
 * @param info The MPI information object
 * @param dry_run If 1, don't read the particle. Only allocates the arrays.
 *
 * Opens the HDF5 file fileName and reads the particles contained
 * in the parts array. N is the returned number of particles found
 * in the file.
 *
 * @warning Can not read snapshot distributed over more than 1 file !!!
 * @todo Read snapshots distributed in more than one file.
 *
 * Calls #error() if an error occurs.
 *
 */
void read_ic_parallel(char* fileName, const struct unit_system* internal_units,
                      double dim[3], struct part** parts, struct gpart** gparts,
                      struct spart** sparts, size_t* Ngas, size_t* Ngparts,
                      size_t* Nstars, int* periodic, int* flag_entropy,
                      int with_hydro, int with_gravity, int with_stars,
                      int mpi_rank, int mpi_size, MPI_Comm comm, MPI_Info info,
                      int dry_run) {

  hid_t h_file = 0, h_grp = 0;
  /* GADGET has only cubic boxes (in cosmological mode) */
  double boxSize[3] = {0.0, -1.0, -1.0};
  long long numParticles[swift_type_count] = {0};
  long long numParticles_highWord[swift_type_count] = {0};
  size_t N[swift_type_count] = {0};
  long long N_total[swift_type_count] = {0};
  long long offset[swift_type_count] = {0};
  int dimension = 3; /* Assume 3D if nothing is specified */
  size_t Ndm = 0;

  /* Open file */
  /* message("Opening file '%s' as IC.", fileName); */
  hid_t h_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(h_plist_id, comm, info);
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, h_plist_id);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Open header to read simulation properties */
  /* message("Reading runtime parameters..."); */
  h_grp = H5Gopen(h_file, "/RuntimePars", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening runtime parameters\n");

  /* Read the relevant information */
  io_read_attribute(h_grp, "PeriodicBoundariesOn", INT, periodic);

  /* Close runtime parameters */
  H5Gclose(h_grp);

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

  /* Read the relevant information and print status */
  int flag_entropy_temp[6];
  io_read_attribute(h_grp, "Flag_Entropy_ICs", INT, flag_entropy_temp);
  *flag_entropy = flag_entropy_temp[0];
  io_read_attribute(h_grp, "BoxSize", DOUBLE, boxSize);
  io_read_attribute(h_grp, "NumPart_Total", LONGLONG, numParticles);
  io_read_attribute(h_grp, "NumPart_Total_HighWord", LONGLONG,
                    numParticles_highWord);

  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    N_total[ptype] =
        (numParticles[ptype]) + (numParticles_highWord[ptype] << 32);

  /* Get the box size if not cubic */
  dim[0] = boxSize[0];
  dim[1] = (boxSize[1] < 0) ? boxSize[0] : boxSize[1];
  dim[2] = (boxSize[2] < 0) ? boxSize[0] : boxSize[2];

  /* Change box size in the 1D and 2D case */
  if (hydro_dimension == 2)
    dim[2] = min(dim[0], dim[1]);
  else if (hydro_dimension == 1)
    dim[2] = dim[1] = dim[0];

  /* message("Found %lld particles in a %speriodic box of size [%f %f %f].", */
  /* 	  N_total[0], (periodic ? "": "non-"), dim[0], */
  /* 	  dim[1], dim[2]); */

  /* Divide the particles among the tasks. */
  for (int ptype = 0; ptype < swift_type_count; ++ptype) {
    offset[ptype] = mpi_rank * N_total[ptype] / mpi_size;
    N[ptype] = (mpi_rank + 1) * N_total[ptype] / mpi_size - offset[ptype];
  }

  /* Close header */
  H5Gclose(h_grp);

  /* Read the unit system used in the ICs */
  struct unit_system* ic_units = malloc(sizeof(struct unit_system));
  if (ic_units == NULL) error("Unable to allocate memory for IC unit system");
  io_read_unit_system(h_file, ic_units);

  /* Tell the user if a conversion will be needed */
  if (mpi_rank == 0) {
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
  }

  /* Convert the dimensions of the box */
  for (int j = 0; j < 3; j++)
    dim[j] *=
        units_conversion_factor(ic_units, internal_units, UNIT_CONV_LENGTH);

  /* Allocate memory to store SPH particles */
  if (with_hydro) {
    *Ngas = N[0];
    if (posix_memalign((void*)parts, part_align,
                       (*Ngas) * sizeof(struct part)) != 0)
      error("Error while allocating memory for particles");
    bzero(*parts, *Ngas * sizeof(struct part));
  }

  /* Allocate memory to store star particles */
  if (with_stars) {
    *Nstars = N[swift_type_star];
    if (posix_memalign((void*)sparts, spart_align,
                       *Nstars * sizeof(struct spart)) != 0)
      error("Error while allocating memory for star particles");
    bzero(*sparts, *Nstars * sizeof(struct spart));
  }

  /* Allocate memory to store gravity particles */
  if (with_gravity) {
    Ndm = N[1];
    *Ngparts = (with_hydro ? N[swift_type_gas] : 0) +
               N[swift_type_dark_matter] +
               (with_stars ? N[swift_type_star] : 0);
    if (posix_memalign((void*)gparts, gpart_align,
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
    if (N_total[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gopen(h_file, partTypeGroupName, H5P_DEFAULT);
    if (h_grp < 0) {
      error("Error while opening particle group %s.", partTypeGroupName);
    }

    int num_fields = 0;
    struct io_props list[100];
    size_t Nparticles = 0;

    /* Read particle fields into the particle structure */
    switch (ptype) {

      case swift_type_gas:
        if (with_hydro) {
          Nparticles = *Ngas;
          hydro_read_particles(*parts, list, &num_fields);
        }
        break;

      case swift_type_dark_matter:
        if (with_gravity) {
          Nparticles = Ndm;
          darkmatter_read_particles(*gparts, list, &num_fields);
        }
        break;

      case swift_type_star:
        if (with_stars) {
          Nparticles = *Nstars;
          star_read_particles(*sparts, list, &num_fields);
        }
        break;

      default:
        message("Particle Type %d not yet supported. Particles ignored", ptype);
    }

    /* Read everything */
    if (!dry_run)
      for (int i = 0; i < num_fields; ++i)
        readArray(h_grp, list[i], Nparticles, N_total[ptype], offset[ptype],
                  internal_units, ic_units);

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* Prepare the DM particles */
  if (!dry_run && with_gravity) io_prepare_dm_gparts(*gparts, Ndm);

  /* Duplicate the hydro particles into gparts */
  if (!dry_run && with_gravity && with_hydro)
    io_duplicate_hydro_gparts(*parts, *gparts, *Ngas, Ndm);

  /* Duplicate the star particles into gparts */
  if (!dry_run && with_gravity && with_stars)
    io_duplicate_star_gparts(*sparts, *gparts, *Nstars, Ndm + *Ngas);

  /* message("Done Reading particles..."); */

  /* Clean up */
  free(ic_units);

  /* Close property handler */
  H5Pclose(h_plist_id);

  /* Close file */
  H5Fclose(h_file);
}

/**
 * @brief Writes an HDF5 output file (GADGET-3 type) with
 *its XMF descriptor
 *
 * @param e The engine containing all the system.
 * @param baseName The common part of the snapshot file name.
 * @param internal_units The #unit_system used internally
 * @param snapshot_units The #unit_system used in the snapshots
 * @param mpi_rank The MPI rank of this node.
 * @param mpi_size The number of MPI ranks.
 * @param comm The MPI communicator.
 * @param info The MPI information object
 *
 * Creates an HDF5 output file and writes the particles
 * contained in the engine. If such a file already exists, it is
 * erased and replaced by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_parallel(struct engine* e, const char* baseName,
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units,
                           int mpi_rank, int mpi_size, MPI_Comm comm,
                           MPI_Info info) {

  hid_t h_file = 0, h_grp = 0;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  const size_t Ntot = e->s->nr_gparts;
  int periodic = e->s->periodic;
  int numFiles = 1;
  struct part* parts = e->s->parts;
  struct gpart* gparts = e->s->gparts;
  struct gpart* dmparts = NULL;
  struct spart* sparts = e->s->sparts;
  static int outputCount = 0;
  FILE* xmfFile = 0;

  /* Number of unassociated gparts */
  const size_t Ndm = Ntot > 0 ? Ntot - (Ngas + Nstars) : 0;

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%s_%04i.hdf5", baseName,
           outputCount);

  /* First time, we need to create the XMF file */
  if (outputCount == 0 && mpi_rank == 0) xmf_create_file(baseName);

  /* Prepare the XMF file for the new entry */
  if (mpi_rank == 0) xmfFile = xmf_prepare_file(baseName);

  /* Open HDF5 file */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Compute offset in the file and total number of
   * particles */
  size_t N[swift_type_count] = {Ngas, Ndm, 0, 0, Nstars, 0};
  long long N_total[swift_type_count] = {0};
  long long offset[swift_type_count] = {0};
  MPI_Exscan(&N, &offset, swift_type_count, MPI_LONG_LONG_INT, MPI_SUM, comm);
  for (int ptype = 0; ptype < swift_type_count; ++ptype)
    N_total[ptype] = offset[ptype] + N[ptype];

  /* The last rank now has the correct N_total. Let's
   * broadcast from there */
  MPI_Bcast(&N_total, 6, MPI_LONG_LONG_INT, mpi_size - 1, comm);

  /* Now everybody konws its offset and the total number of
   * particles of each
   * type */

  /* Write the part of the XMF file corresponding to this
   * specific output */
  if (mpi_rank == 0) xmf_write_outputheader(xmfFile, fileName, e->time);

  /* Open header to write simulation properties */
  /* message("Writing runtime parameters..."); */
  h_grp =
      H5Gcreate(h_file, "/RuntimePars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating runtime parameters group\n");

  /* Write the relevant information */
  io_write_attribute(h_grp, "PeriodicBoundariesOn", INT, &periodic, 1);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Print the relevant information and print status */
  io_write_attribute(h_grp, "BoxSize", DOUBLE, e->s->dim, 3);
  double dblTime = e->time;
  io_write_attribute(h_grp, "Time", DOUBLE, &dblTime, 1);
  int dimension = (int)hydro_dimension;
  io_write_attribute(h_grp, "Dimension", INT, &dimension, 1);

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
  double MassTable[6] = {0., 0., 0., 0., 0., 0.};
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

  /* Print the SPH parameters */
  if (e->policy & engine_policy_hydro) {
    h_grp = H5Gcreate(h_file, "/HydroScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating SPH group");
    hydro_props_print_snapshot(h_grp, e->hydro_properties);
    writeSPHflavour(h_grp);
    H5Gclose(h_grp);
  }

  /* Print the gravity parameters */
  if (e->policy & engine_policy_self_gravity) {
    h_grp = H5Gcreate(h_file, "/GravityScheme", H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) error("Error while creating gravity group");
    gravity_props_print_snapshot(h_grp, e->gravity_properties);
    H5Gclose(h_grp);
  }

  /* Print the runtime parameters */
  h_grp =
      H5Gcreate(h_file, "/Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating parameters group");
  parser_write_params_to_hdf5(e->parameter_file, h_grp);
  H5Gclose(h_grp);

  /* Print the system of Units used in the spashot */
  io_write_unit_system(h_file, snapshot_units, "Units");

  /* Print the system of Units used internally */
  io_write_unit_system(h_file, internal_units, "InternalCodeUnits");

  /* Tell the user if a conversion will be needed */
  if (e->verbose && mpi_rank == 0) {
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
    if (N_total[ptype] == 0) continue;

    /* Add the global information for that particle type to
     * the XMF meta-file */
    if (mpi_rank == 0)
      xmf_write_groupheader(xmfFile, fileName, N_total[ptype],
                            (enum part_type)ptype);

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) {
      error("Error while opening particle group %s.", partTypeGroupName);
    }

    int num_fields = 0;
    struct io_props list[100];
    size_t Nparticles = 0;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        Nparticles = Ngas;
        hydro_write_particles(parts, list, &num_fields);
        break;

      case swift_type_dark_matter:
        /* Allocate temporary array */
        if (posix_memalign((void*)&dmparts, gpart_align,
                           Ndm * sizeof(struct gpart)) != 0)
          error(
              "Error while allocating temporart memory for "
              "DM particles");
        bzero(dmparts, Ndm * sizeof(struct gpart));

        /* Collect the DM particles from gpart */
        io_collect_dm_gparts(gparts, Ntot, dmparts, Ndm);

        /* Write DM particles */
        Nparticles = Ndm;
        darkmatter_write_particles(dmparts, list, &num_fields);
        break;

      case swift_type_star:
        Nparticles = Nstars;
        star_write_particles(sparts, list, &num_fields);
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Write everything */
    for (int i = 0; i < num_fields; ++i)
      writeArray(e, h_grp, fileName, xmfFile, partTypeGroupName, list[i],
                 Nparticles, N_total[ptype], mpi_rank, offset[ptype],
                 internal_units, snapshot_units);

    /* Free temporary array */
    if (dmparts) {
      free(dmparts);
      dmparts = 0;
    }

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    if (mpi_rank == 0) xmf_write_groupfooter(xmfFile, (enum part_type)ptype);
  }

  /* Write LXMF file descriptor */
  if (mpi_rank == 0) xmf_write_outputfooter(xmfFile, outputCount, e->time);

  /* message("Done writing particles..."); */

  /* Close property descriptor */
  H5Pclose(plist_id);

  /* Close file */
  H5Fclose(h_file);

  ++outputCount;
}

#endif /* HAVE_HDF5 */
