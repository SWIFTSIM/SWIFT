/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

#if defined(HAVE_HDF5) && defined(WITH_MPI)


/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <hdf5.h>
#include <math.h>
#include <mpi.h>

#include "const.h"
#include "cycle.h"
#include "lock.h"
#include "task.h"
#include "part.h"
#include "space.h"
#include "scheduler.h"
#include "engine.h"
#include "error.h"
#include "kernel.h"

/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum DATA_TYPE{INT, LONG, LONGLONG, UINT, ULONG, ULONGLONG, FLOAT, DOUBLE, CHAR};

/**
 * @brief The two sorts of data present in the GADGET IC files: compulsory to start a run or optional.
 *
 */
enum DATA_IMPORTANCE{COMPULSORY=1, OPTIONAL=0};

/**
 * @brief Converts a C data type to the HDF5 equivalent. 
 *
 * This function is a trivial wrapper around the HDF5 types but allows
 * to change the exact storage types matching the code types in a transparent way.
 */
static hid_t hdf5Type(enum DATA_TYPE type)
{
  switch(type)
    {
    case INT: return H5T_NATIVE_INT;
    case UINT: return H5T_NATIVE_UINT;
    case LONG: return H5T_NATIVE_LONG;
    case ULONG: return H5T_NATIVE_ULONG;
    case LONGLONG: return H5T_NATIVE_LLONG;
    case ULONGLONG: return H5T_NATIVE_ULLONG;
    case FLOAT: return H5T_NATIVE_FLOAT;
    case DOUBLE: return H5T_NATIVE_DOUBLE;
    case CHAR: return H5T_C_S1;
    default: error("Unknown type"); return 0;
    }
}

/**
 * @brief Returns the memory size of the data type
 */
static size_t sizeOfType(enum DATA_TYPE type)
{
  switch(type)
    {
    case INT: return sizeof(int);
    case UINT: return sizeof(unsigned int);
    case LONG: return sizeof(long);
    case ULONG: return sizeof(unsigned long);
    case LONGLONG: return sizeof(long long);
    case ULONGLONG: return sizeof(unsigned long long);
    case FLOAT: return sizeof(float);
    case DOUBLE: return sizeof(double);
    case CHAR: return sizeof(char);
    default: error("Unknown type"); return 0;
    }
}

/*-----------------------------------------------------------------------------
 * Routines reading an IC file
 *-----------------------------------------------------------------------------*/


/**
 * @brief Reads an attribute from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the attribute to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param data (output) The attribute read from the HDF5 group.
 *
 * Calls #error() if an error occurs.
 */
void readAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data)
{
  hid_t h_attr=0, h_err=0;

  h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if(h_attr < 0)
    {
      error( "Error while opening attribute '%s'" , name );
    }

  h_err = H5Aread(h_attr, hdf5Type(type), data);
  if(h_err < 0)
    {
      error( "Error while reading attribute '%s'" , name );
    }

  H5Aclose(h_attr);
}


/**
 * @brief Reads a data array from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the array to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param N The number of particles.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurence of the field of interest in the parts array
 * @param importance If COMPULSORY, the data must be present in the IC file. If OPTIONAL, the array will be zeroed when the data is not present.
 *
 * @todo A better version using HDF5 hyperslabs to read the file directly into the part array 
 * will be written once the strucutres have been stabilized.
 *  
 * Calls #error() if an error occurs.
 */
void readArrayBackEnd(hid_t grp, char* name, enum DATA_TYPE type, int N, int dim, char* part_c, enum DATA_IMPORTANCE importance)
{
  hid_t h_data=0, h_err=0, h_type=0;
  htri_t exist=0;
  void* temp;
  int i=0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  const size_t partSize = sizeof(struct part);
  char* temp_c = 0;

  /* Check whether the dataspace exists or not */
  exist = H5Lexists(grp, name, 0);
  if(exist < 0)
    {
      error( "Error while checking the existence of data set '%s'." , name );
    }
  else if(exist == 0)
    {
      if(importance == COMPULSORY)
	{
	  error( "Compulsory data set '%s' not present in the file." , name );
	}
      else
	{
	  /* message("Optional data set '%s' not present. Zeroing this particle field...", name);	   */
	  
	  for(i=0; i<N; ++i)
	    memset(part_c+i*partSize, 0, copySize);
	  
	  return;
	}
   }

  /* message( "Reading %s '%s' array...", importance == COMPULSORY ? "compulsory": "optional  ", name); */

  /* Open data space */
  h_data = H5Dopen1(grp, name);
  if(h_data < 0)
    {
      error( "Error while opening data space '%s'." , name );
    }

  /* Check data type */
  h_type = H5Dget_type(h_data);
  if(h_type < 0)
    error("Unable to retrieve data type from the file");
  if(!H5Tequal(h_type, hdf5Type(type)))
    error("Non-matching types between the code and the file");
  
  /* Allocate temporary buffer */
  temp = malloc(N * dim * sizeOfType(type));
  if(temp == NULL)
    error("Unable to allocate memory for temporary buffer");

  /* Read HDF5 dataspace in temporary buffer */
  /* Dirty version that happens to work for vectors but should be improved */
  /* Using HDF5 dataspaces would be better */
  h_err = H5Dread(h_data, hdf5Type(type), H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if(h_err < 0)
    {
      error( "Error while reading data array '%s'." , name );
    }

  /* Copy temporary buffer to particle data */
  temp_c = temp;
  for(i=0; i<N; ++i)
    memcpy(part_c+i*partSize, &temp_c[i*copySize], copySize);
  
  /* Free and close everything */
  free(temp);
  H5Tclose(h_type);
  H5Dclose(h_data);
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
 * @param field The name of the field (C code name as defined in part.h) to fill
 * @param importance Is the data compulsory or not
 *
 */
#define readArray(grp, name, type, N, dim, part, field, importance) readArrayBackEnd(grp, name, type, N, dim, (char*)(&(part[0]).field), importance)

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
 * @todo Read snaphsots distributed in more than one file.
 *
 * Calls #error() if an error occurs.
 *
 */
void read_ic_parallel ( char* fileName, double dim[3], struct part **parts,  int* N, int* periodic)
{
  hid_t h_file=0, h_grp=0;
  double boxSize[3]={0.0,-1.0,-1.0};         /* GADGET has only cubic boxes (in cosmological mode) */
  int numParticles[6]={0};   /* GADGET has 6 particle types. We only keep the type 0*/

  /* Open file */
  /* message("Opening file '%s' as IC.", fileName); */
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(h_file < 0)
    {
      error( "Error while opening file '%s'." , fileName );
    }

  /* Open header to read simulation properties */
  /* message("Reading runtime parameters..."); */
  h_grp = H5Gopen1(h_file, "/RuntimePars");
  if(h_grp < 0)
    error("Error while opening runtime parameters\n");

  /* Read the relevant information */
  readAttribute(h_grp, "PeriodicBoundariesOn", INT, periodic);

  /* Close runtime parameters */
  H5Gclose(h_grp);
  
  /* Open header to read simulation properties */
  /* message("Reading file header..."); */
  h_grp = H5Gopen1(h_file, "/Header");
  if(h_grp < 0)
    error("Error while opening file header\n");
    
  /* Read the relevant information and print status */
  readAttribute(h_grp, "BoxSize", DOUBLE, boxSize);
  readAttribute(h_grp, "NumPart_Total", UINT, numParticles);

  *N = numParticles[0];
  dim[0] = boxSize[0];
  dim[1] = ( boxSize[1] < 0 ) ? boxSize[0] : boxSize[1];
  dim[2] = ( boxSize[2] < 0 ) ? boxSize[0] : boxSize[2];

  /* message("Found %d particles in a %speriodic box of size [%f %f %f].",  */
  /* 	 *N, (periodic ? "": "non-"), dim[0], dim[1], dim[2]); */

  /* Close header */
  H5Gclose(h_grp);

  /* Allocate memory to store particles */
  if(posix_memalign( (void*)parts , part_align , *N * sizeof(struct part)) != 0)
    error("Error while allocating memory for particles");
  bzero( *parts , *N * sizeof(struct part) );

  /* message("Allocated %8.2f MB for particles.", *N * sizeof(struct part) / (1024.*1024.)); */
		  
  /* Open SPH particles group */
  /* message("Reading particle arrays..."); */
  h_grp = H5Gopen1(h_file, "/PartType0");
  if(h_grp < 0)
    error( "Error while opening particle group.\n");

  /* Read arrays */
  readArray(h_grp, "Coordinates", DOUBLE, *N, 3, *parts, x, COMPULSORY);
  readArray(h_grp, "Velocities", FLOAT, *N, 3, *parts, v, COMPULSORY);
  readArray(h_grp, "Masses", FLOAT, *N, 1, *parts, mass, COMPULSORY);
  readArray(h_grp, "SmoothingLength", FLOAT, *N, 1, *parts, h, COMPULSORY);
  readArray(h_grp, "InternalEnergy", FLOAT, *N, 1, *parts, u, COMPULSORY);
  readArray(h_grp, "ParticleIDs", ULONGLONG, *N, 1, *parts, id, COMPULSORY);
  readArray(h_grp, "TimeStep", FLOAT, *N, 1, *parts, dt, OPTIONAL);
  readArray(h_grp, "Acceleration", FLOAT, *N, 3, *parts, a, OPTIONAL);
  readArray(h_grp, "Density", FLOAT, *N, 1, *parts, rho, OPTIONAL );

  /* Close particle group */
  H5Gclose(h_grp);

  /* message("Done Reading particles..."); */

  /* Close file */
  H5Fclose(h_file);
}


/*-----------------------------------------------------------------------------
 * Routines writing an output file
 *-----------------------------------------------------------------------------*/


/* File descriptor output routine */
void createXMFfile();
FILE* prepareXMFfile();
void writeXMFfooter(FILE* xmfFile);
void writeXMFheader(FILE* xmfFile, long long N, char* hdfFileName,  float time);
void writeXMFline(FILE* xmfFile, char* fileName, char* name, long long N, int dim, enum DATA_TYPE type);


/**
 * @brief Write an attribute to a given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param name The name of the attribute to write.
 * @param type The #DATA_TYPE of the attribute.
 * @param data The attribute to write.
 * @param num The number of elements to write
 *
 * Calls #error() if an error occurs.
 */
void writeAttribute(hid_t grp, char* name, enum DATA_TYPE type, void* data, int num)
{
  hid_t h_space=0, h_attr=0, h_err=0;
  hsize_t dim[1]={num};

  h_space = H5Screate(H5S_SIMPLE);
  if(h_space < 0)
    {
      error( "Error while creating dataspace for attribute '%s'." , name );
    }

  h_err = H5Sset_extent_simple(h_space, 1, dim, NULL);
  if(h_err < 0)
    {
      error( "Error while changing dataspace shape for attribute '%s'." , name );
    }

  h_attr = H5Acreate1(grp, name, hdf5Type(type), h_space, H5P_DEFAULT);
  if(h_attr < 0)
    {
      error( "Error while creating attribute '%s'.", name );
    }

  h_err = H5Awrite(h_attr, hdf5Type(type), data);
  if(h_err < 0)
    {
      error( "Error while reading attribute '%s'." , name );
    }

  H5Sclose(h_space);
  H5Aclose(h_attr);
}

/**
 * @brief Write a string as an attribute to a given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param name The name of the attribute to write.
 * @param str The string to write.
 * @param length The length of the string
 *
 * Calls #error() if an error occurs.
 */
void writeStringAttribute(hid_t grp, char* name, char* str, int length)
{
  hid_t h_space=0, h_attr=0, h_err=0, h_type=0;

  h_space = H5Screate(H5S_SCALAR);
  if(h_space < 0)
    {
      error( "Error while creating dataspace for attribute '%s'." , name );
    }

  h_type = H5Tcopy(H5T_C_S1);
  if(h_type < 0)
    {
      error( "Error while copying datatype 'H5T_C_S1'." );
    }

  h_err = H5Tset_size(h_type, length);
  if(h_err < 0)
    {
      error( "Error while resizing attribute tyep to '%i'." , length );
    }

  h_attr = H5Acreate1(grp, name, h_type, h_space, H5P_DEFAULT);
  if(h_attr < 0)
    {
      error( "Error while creating attribute '%s'." , name );
    }

  h_err = H5Awrite(h_attr, h_type, str );
  if(h_err < 0)
    {
      error( "Error while reading attribute '%s'." , name );
    }

  H5Tclose(h_type);
  H5Sclose(h_space);
  H5Aclose(h_attr);
}

/**
 * @brief Writes a double value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_d(hid_t grp, char* name, double data)
{
  writeAttribute(grp, name, DOUBLE, &data, 1);
}

/**
 * @brief Writes a float value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_f(hid_t grp, char* name, float data)
{
  writeAttribute(grp, name, FLOAT, &data, 1);
}

/**
 * @brief Writes an int value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */

void writeAttribute_i(hid_t grp, char* name, int data)
{
  writeAttribute(grp, name, INT, &data, 1);
}

/**
 * @brief Writes a long value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param data The value to write
 */
void writeAttribute_l(hid_t grp, char* name, long data)
{
  writeAttribute(grp, name, LONG, &data, 1);
}

/**
 * @brief Writes a string value as an attribute
 * @param grp The groupm in which to write
 * @param name The name of the attribute
 * @param str The string to write
 */
void writeAttribute_s(hid_t grp, char* name, char* str)
{
  writeStringAttribute(grp, name, str, strlen(str));
}


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
 * @param N_total Total number of particles across all cores
 * @param offset Offset in the array where this mpi task starts writing
 * @param part_c A (char*) pointer on the first occurence of the field of interest in the parts array
 *
 * @todo A better version using HDF5 hyperslabs to write the file directly from the part array
 * will be written once the strucutres have been stabilized.
 *
 * Calls #error() if an error occurs.
 */
void writeArrayBackEnd(hid_t grp, char* fileName, FILE* xmfFile, char* name, enum DATA_TYPE type, int N, int dim, int N_total, int mpi_rank, int offset, char* part_c)
{
  hid_t h_data=0, h_err=0, h_memspace=0, h_filespace=0, h_plist_id=0;
  void* temp = 0;
  int i=0, rank=0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  const size_t partSize = sizeof(struct part);
  char* temp_c = 0;
  hsize_t shape[2], shape_total[2], offsets[2];

  /* message("Writing '%s' array...", name); */

  /* Allocate temporary buffer */
  temp = malloc(N * dim * sizeOfType(type));
  if(temp == NULL)
    error("Unable to allocate memory for temporary buffer");

  /* Copy particle data to temporary buffer */
  temp_c = temp;
  for(i=0; i<N; ++i)
    memcpy(&temp_c[i*copySize], part_c+i*partSize, copySize);

  /* Create data space */
  h_memspace = H5Screate(H5S_SIMPLE);
  if(h_memspace < 0)
    {
      error( "Error while creating data space (memory) for field '%s'." , name );
    }

  h_filespace = H5Screate(H5S_SIMPLE);
  if(h_filespace < 0)
    {
      error( "Error while creating data space (file) for field '%s'." , name );
    }
  
  if(dim > 1)
    {
      rank = 2;
      shape[0] = N; shape[1] = dim;
      shape_total[0] = N_total; shape_total[1] = dim;
      offsets[0] = offset; offsets[1] = 0;
    }
  else
    {
      rank = 1;
      shape[0] = N; shape[1] = 0;
      shape_total[0] = N_total; shape_total[1] = 0;
      offsets[0] = offset; offsets[1] = 0;
    }
  
  /* Change shape of memory data space */
  h_err = H5Sset_extent_simple(h_memspace, rank, shape, NULL);
  if(h_err < 0)
    {
      error( "Error while changing data space (memory) shape for field '%s'." , name );
    }

  /* Change shape of file data space */
  h_err = H5Sset_extent_simple(h_filespace, rank, shape_total, NULL);
  if(h_err < 0)
    {
      error( "Error while changing data space (file) shape for field '%s'." , name );
    }
  
  /* Create dataset */
  h_data = H5Dcreate1(grp, name, hdf5Type(type), h_filespace, H5P_DEFAULT);
  if(h_data < 0)
    {
      error( "Error while creating dataspace '%s'." , name );
    }
  
  H5Sclose(h_filespace);
  h_filespace = H5Dget_space(h_data);
  H5Sselect_hyperslab(h_filespace, H5S_SELECT_SET, offsets, NULL, shape, NULL);

  /* Create property list for collective dataset write.    */
  h_plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE);

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, hdf5Type(type), h_memspace, h_filespace, h_plist_id, temp);
  if(h_err < 0)
    {
      error( "Error while writing data array '%s'." , name );
    }

  /* Write XMF description for this data set */
  if(mpi_rank == 0)
    writeXMFline(xmfFile, fileName, name, N_total, dim, type);

  
  /* Free and close everything */
  free(temp);
  H5Dclose(h_data);
  H5Pclose(h_plist_id);
  H5Sclose(h_memspace);
}

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write from this core.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param N_total Total number of particles across all cores
 * @param mpi_rank The MPI task rank calling the function
 * @param offset Offset in the array where this mpi task starts writing
 * @param part A (char*) pointer on the first occurence of the field of interest in the parts array
 * @param field The name (code name) of the field to read from.
 *
 */
#define writeArray(grp, fileName, xmfFile, name, type, N, dim, part, N_total, mpi_rank, offset, field) writeArrayBackEnd(grp, fileName, xmfFile, name, type, N, dim, N_total, mpi_rank, offset, (char*)(&(part[0]).field))

/**
 * @brief Writes an HDF5 output file (GADGET-3 type) with its XMF descriptor
 *
 * @param e The engine containing all the system.
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one. 
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_parallel (struct engine *e, int mpi_rank, int mpi_size, MPI_Comm comm, MPI_Info info)
{
  
  hid_t h_file=0, h_grp=0, h_grpsph=0;
  int N = e->s->nr_parts;
  int periodic = e->s->periodic;
  unsigned int numParticles[6]={N,0};
  unsigned int numParticlesHighWord[6]={0};
  unsigned int flagEntropy[6]={0};
  long long N_total = 0, offset = 0;
  double offset_d = 0., N_d = 0., N_total_d = 0.;
  int numFiles = 1;
  struct part* parts = e->s->parts;
  FILE* xmfFile = 0;
  static int outputCount = 0;
  
  /* File name */
  char fileName[200];
  sprintf(fileName, "output_%03i.hdf5", outputCount);

  /* First time, we need to create the XMF file */
  if(outputCount == 0 && mpi_rank == 0)
    createXMFfile();
  
  /* Prepare the XMF file for the new entry */
  if(mpi_rank == 0)
    xmfFile = prepareXMFfile();

  /* Open HDF5 file */
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  if(h_file < 0)
    {
      error( "Error while opening file '%s'." , fileName );
    }

  /* Compute offset in the file and total number of particles */
  /* Done using double to allow for up to 2^50=10^15 particles */
  N_d = (double)N;
  MPI_Exscan(&N_d, &offset_d, 1, MPI_DOUBLE, MPI_SUM, comm);
  N_total_d = offset_d + N_d;
  MPI_Bcast(&N_total_d, 1, MPI_DOUBLE, mpi_size-1, comm);
  if(N_total_d > 1.e15)
    error("Error while computing the offest for parallel output: Simulation has more than 10^15 particles.\n");
  N_total = (long long) N_total_d;
  offset = (long long) offset_d;

  /* Write the part of the XMF file corresponding to this specific output */
  if(mpi_rank == 0)
    writeXMFheader(xmfFile, N_total, fileName, e->time);

  /* Open header to write simulation properties */
  /* message("Writing runtime parameters..."); */
  h_grp = H5Gcreate1(h_file, "/RuntimePars", 0);
  if(h_grp < 0)
    error("Error while creating runtime parameters group\n");

  /* Write the relevant information */
  writeAttribute(h_grp, "PeriodicBoundariesOn", INT, &periodic, 1);

  /* Close runtime parameters */
  H5Gclose(h_grp);
  
  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate1(h_file, "/Header", 0);
  if(h_grp < 0)
    error("Error while creating file header\n");
    
  /* Print the relevant information and print status */
  writeAttribute(h_grp, "BoxSize", DOUBLE, e->s->dim, 1);
  writeAttribute(h_grp, "NumPart_ThisFile", UINT, numParticles, 6);
  writeAttribute(h_grp, "Time", DOUBLE, &e->time, 1);

  /* GADGET-2 legacy values */
  numParticles[0] = (unsigned int) N_total ;
  writeAttribute(h_grp, "NumPart_Total", UINT, numParticles, 6);
  numParticlesHighWord[0] = (unsigned int) (N_total >> 32);
  writeAttribute(h_grp, "NumPart_Total_HighWord", UINT, numParticlesHighWord, 6);
  double MassTable[6] = {0., 0., 0., 0., 0., 0.};
  writeAttribute(h_grp, "MassTable", DOUBLE, MassTable, 6);
  writeAttribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy, 6);
  writeAttribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

  /* Print the SPH parameters */
  h_grpsph = H5Gcreate1(h_file, "/Header/SPH", 0);
  if(h_grpsph < 0)
    error("Error while creating SPH header\n");

  writeAttribute_f(h_grpsph, "Kernel eta", const_eta_kernel);
  writeAttribute_f(h_grpsph, "Weighted N_ngb", kernel_nwneigh);
  writeAttribute_f(h_grpsph, "Delta N_ngb", const_delta_nwneigh);
  writeAttribute_f(h_grpsph, "Hydro gamma", const_hydro_gamma);

#ifdef LEGACY_GADGET2_SPH
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model", "(No treatment) Legacy Gadget-2 as in Springel (2005)");  
  writeAttribute_s(h_grpsph, "Viscosity Model", "Legacy Gadget-2 as in Springel (2005)");  
  writeAttribute_f(h_grpsph, "Viscosity alpha", const_viscosity_alpha);  
  writeAttribute_f(h_grpsph, "Viscosity beta", 3.f);  
#else
  writeAttribute_s(h_grpsph, "Thermal Conductivity Model", "Price (2008) without switch");  
  writeAttribute_f(h_grpsph, "Thermal Conductivity alpha", const_conductivity_alpha);  
  writeAttribute_s(h_grpsph, "Viscosity Model", "Morris & Monaghan (1997), Rosswog, Davies, Thielemann & Piran (2000) with additional Balsara (1995) switch");  
  writeAttribute_f(h_grpsph, "Viscosity alpha_min", const_viscosity_alpha_min);  
  writeAttribute_f(h_grpsph, "Viscosity alpha_max", const_viscosity_alpha_max);  
  writeAttribute_f(h_grpsph, "Viscosity beta", 2.f);  
  writeAttribute_f(h_grpsph, "Viscosity decay length", const_viscosity_length);  
#endif

  writeAttribute_f(h_grpsph, "CFL parameter", const_cfl);  
  writeAttribute_f(h_grpsph, "Maximal ln(Delta h) change over dt", const_ln_max_h_change);  
  writeAttribute_f(h_grpsph, "Maximal Delta h change over dt", exp(const_ln_max_h_change));  
  writeAttribute_f(h_grpsph, "Maximal Delta u change over dt", const_max_u_change);  
  writeAttribute_s(h_grpsph, "Kernel", kernel_name);

  /* Close headers */
  H5Gclose(h_grpsph);
  H5Gclose(h_grp);
		  
  /* Create SPH particles group */
  /* message("Writing particle arrays..."); */
  h_grp = H5Gcreate1(h_file, "/PartType0", 0);
  if(h_grp < 0)
    error( "Error while creating particle group.\n");

  /* Write arrays */
  writeArray(h_grp, fileName, xmfFile, "Coordinates", DOUBLE, N, 3, parts, N_total, mpi_rank, offset, x);
  writeArray(h_grp, fileName, xmfFile, "Velocities", FLOAT, N, 3, parts, N_total, mpi_rank, offset, v);
  writeArray(h_grp, fileName, xmfFile, "Masses", FLOAT, N, 1, parts, N_total, mpi_rank, offset, mass);
  writeArray(h_grp, fileName, xmfFile, "SmoothingLength", FLOAT, N, 1, parts, N_total, mpi_rank, offset, h);
  writeArray(h_grp, fileName, xmfFile, "InternalEnergy", FLOAT, N, 1, parts, N_total, mpi_rank, offset, u);
  writeArray(h_grp, fileName, xmfFile, "ParticleIDs", ULONGLONG, N, 1, parts, N_total, mpi_rank, offset, id);
  writeArray(h_grp, fileName, xmfFile, "TimeStep", FLOAT, N, 1, parts, N_total, mpi_rank, offset, dt);
  writeArray(h_grp, fileName, xmfFile, "Acceleration", FLOAT, N, 3, parts, N_total, mpi_rank, offset, a);
  writeArray(h_grp, fileName, xmfFile, "Density", FLOAT, N, 1, parts, N_total, mpi_rank, offset, rho);

  /* Close particle group */
  H5Gclose(h_grp);

  /* Write LXMF file descriptor */
  if(mpi_rank == 0)
    writeXMFfooter(xmfFile);

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  ++outputCount;
}



/* ------------------------------------------------------------------------------------------------ 
 * This part writes the XMF file descriptor enabling a visualisation through ParaView
 * ------------------------------------------------------------------------------------------------ */

/**
 * @brief Writes the lines corresponding to an array of the HDF5 output
 *
 * @param xmfFile The file in which to write
 * @param fileName The name of the HDF5 file associated to this XMF descriptor.
 * @param name The name of the array in the HDF5 file.
 * @param N The number of particles.
 * @param dim The dimension of the quantity (1 for scalars, 3 for vectors).
 * @param type The type of the data to write.
 *
 * @todo Treat the types in a better way.
 */
void writeXMFline(FILE* xmfFile, char* fileName, char* name, long long N, int dim, enum DATA_TYPE type )
{
  fprintf(xmfFile, "<Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Node\">\n", name, dim == 1 ? "Scalar": "Vector");
  if(dim == 1)
    fprintf(xmfFile, "<DataItem Dimensions=\"%lld\" NumberType=\"Double\" Precision=\"%d\" Format=\"HDF\">%s:/PartType0/%s</DataItem>\n", N, type==FLOAT ? 4:8, fileName, name);
  else
    fprintf(xmfFile, "<DataItem Dimensions=\"%lld %d\" NumberType=\"Double\" Precision=\"%d\" Format=\"HDF\">%s:/PartType0/%s</DataItem>\n", N, dim, type==FLOAT ? 4:8, fileName, name);
  fprintf(xmfFile, "</Attribute>\n");
}


/**
 * @brief Prepares the XMF file for the new entry
 *
 * Creates a temporary file on the disk in order to copy the right lines.
 *
 * @todo Use a proper XML library to avoid stupid copies.
 */
FILE* prepareXMFfile()
{
  char buffer[1024];

  FILE* xmfFile = fopen("output.xmf", "r");
  FILE* tempFile = fopen("output_temp.xmf", "w");

  if(xmfFile == NULL)
    error("Unable to open current XMF file.");

  if(tempFile == NULL)
    error("Unable to open temporary file.");


  /* First we make a temporary copy of the XMF file and count the lines */
  int counter = 0;
  while (fgets(buffer, 1024, xmfFile) != NULL)
    {
      counter++;
      fprintf(tempFile, "%s", buffer);
    }
  fclose(tempFile);
  fclose(xmfFile);
  
  /* We then copy the XMF file back up to the closing lines */
  xmfFile = fopen("output.xmf", "w");
  tempFile = fopen("output_temp.xmf", "r");

  if(xmfFile == NULL)
    error("Unable to open current XMF file.");

  if(tempFile == NULL)
    error("Unable to open temporary file.");

  int i = 0;
  while (fgets(buffer, 1024, tempFile) != NULL && i < counter - 3)
    {
      i++;
      fprintf(xmfFile, "%s", buffer);
    }
  fprintf(xmfFile, "\n");
  fclose(tempFile);
  remove("output_temp.xmf");
 
  return xmfFile;
}

/**
 * @brief Writes the begin of the XMF file
 *
 * @todo Exploit the XML nature of the XMF format to write a proper XML writer and simplify all the XMF-related stuff.
 */
void createXMFfile()
{
  FILE* xmfFile = fopen("output.xmf", "w");

  fprintf(xmfFile, "<?xml version=\"1.0\" ?> \n");
  fprintf(xmfFile, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
  fprintf(xmfFile, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.1\">\n");
  fprintf(xmfFile, "<Domain>\n");
  fprintf(xmfFile, "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n\n");

  fprintf(xmfFile, "</Grid>\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");

  fclose(xmfFile);
}


/**
 * @brief Writes the part of the XMF entry presenting the geometry of the snapshot
 *
 * @param xmfFile The file to write in.
 * @param Nparts The number of particles.
 * @param hdfFileName The name of the HDF5 file corresponding to this output.
 * @param time The current simulation time.
 */
void writeXMFheader(FILE* xmfFile, long long Nparts, char* hdfFileName, float time)
{
  /* Write end of file */
  
  fprintf(xmfFile, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");
  fprintf(xmfFile, "<Time Type=\"Single\" Value=\"%f\"/>\n", time);
  fprintf(xmfFile, "<Grid Name=\"Gas\" GridType=\"Uniform\">\n");
  fprintf(xmfFile, "<Topology TopologyType=\"Polyvertex\" Dimensions=\"%lld\"/>\n", Nparts);
  fprintf(xmfFile, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf(xmfFile, "<DataItem Dimensions=\"%lld 3\" NumberType=\"Double\" Precision=\"8\" Format=\"HDF\">%s:/PartType0/Coordinates</DataItem>\n", Nparts, hdfFileName);
  fprintf(xmfFile, "</Geometry>");
}


/**
 * @brief Writes the end of the XMF file (closes all open markups)
 *
 * @param xmfFile The file to write in.
 */
void writeXMFfooter(FILE* xmfFile)
{
  /* Write end of the section of this time step */
  
  fprintf(xmfFile, "\n</Grid>\n");
  fprintf(xmfFile, "</Grid>\n");
  fprintf(xmfFile, "\n</Grid>\n");
  fprintf(xmfFile, "</Domain>\n");
  fprintf(xmfFile, "</Xdmf>\n");
  
  fclose(xmfFile);
}

#endif  /* HAVE_HDF5 */


