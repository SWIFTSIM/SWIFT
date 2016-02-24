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

/* This object's header. */
#include "single_io.h"

/* Local includes. */
#include "const.h"
#include "common_io.h"
#include "error.h"

/*-----------------------------------------------------------------------------
 * Routines reading an IC file
 *-----------------------------------------------------------------------------*/

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
void readArrayBackEnd(hid_t grp, char* name, enum DATA_TYPE type, int N,
                      int dim, char* part_c, enum DATA_IMPORTANCE importance) {
  hid_t h_data = 0, h_err = 0, h_type = 0;
  htri_t exist = 0;
  void* temp;
  int i = 0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  const size_t partSize = sizeof(struct part);
  char* temp_c = 0;

  /* Check whether the dataspace exists or not */
  exist = H5Lexists(grp, name, 0);
  if (exist < 0) {
    error("Error while checking the existence of data set '%s'.", name);
  } else if (exist == 0) {
    if (importance == COMPULSORY) {
      error("Compulsory data set '%s' not present in the file.", name);
    } else {
      /* message("Optional data set '%s' not present. Zeroing this particle
       * field...", name);	   */

      for (i = 0; i < N; ++i) memset(part_c + i * partSize, 0, copySize);

      return;
    }
  }

  /* message( "Reading %s '%s' array...", importance == COMPULSORY ?
   * "compulsory": "optional  ", name); */

  /* Open data space */
  h_data = H5Dopen(grp, name, H5P_DEFAULT);
  if (h_data < 0) {
    error("Error while opening data space '%s'.", name);
  }

  /* Check data type */
  h_type = H5Dget_type(h_data);
  if (h_type < 0) error("Unable to retrieve data type from the file");
  // if (!H5Tequal(h_type, hdf5Type(type)))
  //  error("Non-matching types between the code and the file");

  /* Allocate temporary buffer */
  temp = malloc(N * dim * typeSize);
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Read HDF5 dataspace in temporary buffer */
  /* Dirty version that happens to work for vectors but should be improved */
  /* Using HDF5 dataspaces would be better */
  h_err = H5Dread(h_data, hdf5Type(type), H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if (h_err < 0) {
    error("Error while reading data array '%s'.", name);
  }

  /* Copy temporary buffer to particle data */
  temp_c = temp;
  for (i = 0; i < N; ++i)
    memcpy(part_c + i * partSize, &temp_c[i * copySize], copySize);

  /* Free and close everything */
  free(temp);
  H5Tclose(h_type);
  H5Dclose(h_data);
}

/*-----------------------------------------------------------------------------
 * Routines writing an output file
 *-----------------------------------------------------------------------------*/

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurrence of the field of
 *interest in the parts array
 * @param us The UnitSystem currently in use
 * @param convFactor The UnitConversionFactor for this array
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 *the part array
 * will be written once the structures have been stabilized.
 *
 * Calls #error() if an error occurs.
 */
void writeArrayBackEnd(hid_t grp, char* fileName, FILE* xmfFile, char* name,
                       enum DATA_TYPE type, int N, int dim, char* part_c,
                       struct UnitSystem* us,
                       enum UnitConversionFactor convFactor) {
  hid_t h_data = 0, h_err = 0, h_space = 0, h_prop = 0;
  void* temp = 0;
  int i = 0, rank = 0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  const size_t partSize = sizeof(struct part);
  char* temp_c = 0;
  hsize_t shape[2];
  hsize_t chunk_shape[2];
  char buffer[150];

  /* message("Writing '%s' array...", name); */

  /* Allocate temporary buffer */
  temp = malloc(N * dim * sizeOfType(type));
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Copy particle data to temporary buffer */
  temp_c = temp;
  for (i = 0; i < N; ++i)
    memcpy(&temp_c[i * copySize], part_c + i * partSize, copySize);

  /* Create data space */
  h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0) {
    error("Error while creating data space for field '%s'.", name);
  }

  if (dim > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = dim;
    chunk_shape[0] = 1 << 16; /* Just a guess...*/
    chunk_shape[1] = dim;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 16; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Change shape of data space */
  h_err = H5Sset_extent_simple(h_space, rank, shape, NULL);
  if (h_err < 0) {
    error("Error while changing data space shape for field '%s'.", name);
  }

  /* Dataset properties */
  h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Set chunk size */
  h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
  if (h_err < 0) {
    error("Error while setting chunk size (%lld, %lld) for field '%s'.",
          chunk_shape[0], chunk_shape[1], name);
  }

  /* Impose data compression */
  h_err = H5Pset_deflate(h_prop, 4);
  if (h_err < 0) {
    error("Error while setting compression options for field '%s'.", name);
  }

  /* Create dataset */
  h_data = H5Dcreate(grp, name, hdf5Type(type), h_space, H5P_DEFAULT, h_prop,
                     H5P_DEFAULT);
  if (h_data < 0) {
    error("Error while creating dataspace '%s'.", name);
  }

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, hdf5Type(type), h_space, H5S_ALL, H5P_DEFAULT, temp);
  if (h_err < 0) {
    error("Error while writing data array '%s'.", name);
  }

  /* Write XMF description for this data set */
  writeXMFline(xmfFile, fileName, name, N, dim, type);

  /* Write unit conversion factors for this data set */
  conversionString(buffer, us, convFactor);
  writeAttribute_d(h_data, "CGS conversion factor",
                   conversionFactor(us, convFactor));
  writeAttribute_f(h_data, "h-scale exponent", hFactor(us, convFactor));
  writeAttribute_f(h_data, "a-scale exponent", aFactor(us, convFactor));
  writeAttribute_s(h_data, "Conversion factor", buffer);

  /* Free and close everything */
  free(temp);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
}

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group from which to read.
 * @param name The name of the array to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param N The number of particles.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part The array of particles to fill
 * @param N_total Unused parameter in non-MPI mode
 * @param offset Unused parameter in non-MPI mode
 * @param field The name of the field (C code name as defined in part.h) to fill
 * @param importance Is the data compulsory or not
 *
 */
#define readArray(grp, name, type, N, dim, part, N_total, offset, field, \
                  importance)                                            \
  readArrayBackEnd(grp, name, type, N, dim, (char*)(&(part[0]).field),   \
                   importance)

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part A (char*) pointer on the first occurrence of the field of
 * interest in the parts array
 * @param N_total Unused parameter in non-MPI mode
 * @param mpi_rank Unused parameter in non-MPI mode
 * @param offset Unused parameter in non-MPI mode
 * @param field The name (code name) of the field to read from.
 * @param us The UnitSystem currently in use
 * @param convFactor The UnitConversionFactor for this array
 *
 */
#define writeArray(grp, fileName, xmfFile, name, type, N, dim, part, N_total, \
                   mpi_rank, offset, field, us, convFactor)                   \
  writeArrayBackEnd(grp, fileName, xmfFile, name, type, N, dim,               \
                    (char*)(&(part[0]).field), us, convFactor)

/* Import the right hydro definition */
#include "hydro_io.h"

/**
 * @brief Reads an HDF5 initial condition file (GADGET-3 type)
 *
 * @param fileName The file to read.
 * @param dim (output) The dimension of the volume read from the file.
 * @param parts (output) The array of #part read from the file.
 * @param N (output) The number of particles read from the file.
 * @param periodic (output) 1 if the volume is periodic, 0 if not.
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
void read_ic_single(char* fileName, double dim[3], struct part** parts, int* N,
                    int* periodic) {
  hid_t h_file = 0, h_grp = 0;
  double boxSize[3] = {0.0, -1.0, -1.0};
  /* GADGET has only cubic boxes (in cosmological mode) */
  int numParticles[6] = {0};
  /* GADGET has 6 particle types. We only keep the type 0*/

  /* Open file */
  /* message("Opening file '%s' as IC.", fileName); */
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Open header to read simulation properties */
  /* message("Reading runtime parameters..."); */
  h_grp = H5Gopen(h_file, "/RuntimePars", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening runtime parameters\n");

  /* Read the relevant information */
  readAttribute(h_grp, "PeriodicBoundariesOn", INT, periodic);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to read simulation properties */
  /* message("Reading file header..."); */
  h_grp = H5Gopen(h_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening file header\n");

  /* Read the relevant information and print status */
  readAttribute(h_grp, "BoxSize", DOUBLE, boxSize);
  readAttribute(h_grp, "NumPart_Total", UINT, numParticles);

  *N = numParticles[0];
  dim[0] = boxSize[0];
  dim[1] = (boxSize[1] < 0) ? boxSize[0] : boxSize[1];
  dim[2] = (boxSize[2] < 0) ? boxSize[0] : boxSize[2];

  /* message("Found %d particles in a %speriodic box of size [%f %f %f].",  */
  /* 	 *N, (periodic ? "": "non-"), dim[0], dim[1], dim[2]); */

  /* Close header */
  H5Gclose(h_grp);

  /* Allocate memory to store particles */
  if (posix_memalign((void*)parts, part_align, *N * sizeof(struct part)) != 0)
    error("Error while allocating memory for particles");
  bzero(*parts, *N * sizeof(struct part));

  /* message("Allocated %8.2f MB for particles.", *N * sizeof(struct part) /
   * (1024.*1024.)); */

  /* Open SPH particles group */
  /* message("Reading particle arrays..."); */
  h_grp = H5Gopen(h_file, "/PartType0", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening particle group.\n");

  /* Read particle fields into the particle structure */
  hydro_read_particles(h_grp, *N, *N, 0, *parts);

  /* Close particle group */
  H5Gclose(h_grp);

  /* message("Done Reading particles..."); */

  /* Close file */
  H5Fclose(h_file);
}

/**
 * @brief Writes an HDF5 output file (GADGET-3 type) with its XMF descriptor
 *
 * @param e The engine containing all the system.
 * @param us The UnitSystem used for the conversion of units in the output
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_single(struct engine* e, struct UnitSystem* us) {

  hid_t h_file = 0, h_grp = 0, h_grpsph = 0;
  int N = e->s->nr_parts;
  int periodic = e->s->periodic;
  int numParticles[6] = {N, 0};
  int numParticlesHighWord[6] = {0};
  int numFiles = 1;
  struct part* parts = e->s->parts;
  FILE* xmfFile = 0;
  static int outputCount = 0;

  /* File name */
  char fileName[200];
  sprintf(fileName, "output_%03i.hdf5", outputCount);

  /* First time, we need to create the XMF file */
  if (outputCount == 0) createXMFfile();

  /* Prepare the XMF file for the new entry */
  xmfFile = prepareXMFfile();

  /* Write the part corresponding to this specific output */
  writeXMFheader(xmfFile, N, fileName, e->time);

  /* Open file */
  /* message("Opening file '%s'.", fileName); */
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Open header to write simulation properties */
  /* message("Writing runtime parameters..."); */
  h_grp =
      H5Gcreate(h_file, "/RuntimePars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating runtime parameters group\n");

  /* Write the relevant information */
  writeAttribute(h_grp, "PeriodicBoundariesOn", INT, &periodic, 1);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Print the relevant information and print status */
  writeAttribute(h_grp, "BoxSize", DOUBLE, e->s->dim, 3);
  writeAttribute(h_grp, "NumPart_ThisFile", UINT, numParticles, 6);
  double dblTime = e->time;
  writeAttribute(h_grp, "Time", DOUBLE, &dblTime, 1);

  /* GADGET-2 legacy values */
  writeAttribute(h_grp, "NumPart_Total", UINT, numParticles, 6);
  writeAttribute(h_grp, "NumPart_Total_HighWord", UINT, numParticlesHighWord,
                 6);
  double MassTable[6] = {0., 0., 0., 0., 0., 0.};
  writeAttribute(h_grp, "MassTable", DOUBLE, MassTable, 6);
  writeAttribute(h_grp, "Flag_Entropy_ICs", UINT, numParticlesHighWord, 6);
  writeAttribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

  /* Close header */
  H5Gclose(h_grp);

  /* Print the code version */
  writeCodeDescription(h_file);

  /* Print the SPH parameters */
  h_grpsph = H5Gcreate(h_file, "/SPH", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grpsph < 0) error("Error while creating SPH group");
  writeSPHflavour(h_grpsph);
  H5Gclose(h_grpsph);

  /* Print the system of Units */
  writeUnitSystem(h_file, us);

  /* Create SPH particles group */
  /* message("Writing particle arrays..."); */
  h_grp =
      H5Gcreate(h_file, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating particle group.\n");

  /* Write particle fields from the particle structure */
  hydro_write_particles(h_grp, fileName, xmfFile, N, N, 0, 0, parts, us);

  /* Close particle group */
  H5Gclose(h_grp);

  /* Write LXMF file descriptor */
  writeXMFfooter(xmfFile);

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  ++outputCount;
}

#endif /* HAVE_HDF5 */
