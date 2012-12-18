/*******************************************************************************
 * This file is part of GadgetSMP.
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

#ifdef HAVE_HDF5


/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <hdf5.h>


#include "lock.h"
#include "task.h"
#include "part.h"
#include "space.h"

/* Error macro. */
#define error(s) { fprintf( stderr , "%s:%s():%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }


/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum DATA_TYPE{INT, LONG, LONGLONG, UINT, ULONG, ULONGLONG, FLOAT, DOUBLE};

/**
 * @brief Converts a C data type to the HDF5 equivalent. 
 *
 * This function is a trivial wrapper around the HDF5 types but allows
 * to change the exact storage types matching the code types in a transparent way.
 */
hid_t hdf5Type(enum DATA_TYPE type)
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
    default: error("Unknown type");
    }
}


/**
 * @brief Returns the memory size of the data type
 */
size_t sizeOfType(enum DATA_TYPE type)
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
    default: error("Unknown type");
    }
}


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
      char buf[100];
      sprintf(buf, "Error while opening attribute '%s'\n", name);
      error(buf);
    }

  h_err = H5Aread(h_attr, hdf5Type(type), data);
  if(h_err < 0)
    {
      char buf[100];
      sprintf(buf, "Error while reading attribute '%s'\n", name);
      error(buf);
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
 *
 * @todo A better version using HDF5 hyperslabs to read the file directly into the part array 
 * will be written once the strucutres have been stabilized.
 *  
 * Calls #error() if an error occurs.
 */
void readArrayBackEnd(hid_t grp, char* name, enum DATA_TYPE type, int N, int dim, char* part_c)
{
  hid_t h_data=0, h_err=0, h_type=0;
  void* temp;
  int i=0;
  size_t typeSize = sizeOfType(type);
  size_t copySize = typeSize * dim;
  size_t partSize = sizeof(struct part);
  char* temp_c = 0;

  printf("readArray: Reading '%s' array...\n", name);

  /* Open data space */
  h_data = H5Dopen(grp, name);
  if(h_data < 0)
    {
      char buf[100];
      sprintf(buf, "Error while opening data space '%s'\n", name);
      error(buf);
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
      char buf[100];
      sprintf(buf, "Error while reading data array '%s'\n", name);
      error(buf);
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
 *
 */
#define readArray(grp, name, type, N, dim, part, field) readArrayBackEnd(grp, name, type, N, dim, (char*)(&(part[0]).field))

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
void read_ic ( char* fileName, double dim[3], struct part **parts,  int* N, int* periodic)
{
  hid_t h_file=0, h_grp=0;
  double boxSize=0.;         /* GADGET has only cubic boxes (in cosmological mode) */
  int numParticles[6]={0};   /* GADGET has 6 particle types. We only keep the type 0*/

  /* Open file */
  printf("read_ic: Opening file '%s' as IC.\n", fileName);
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(h_file < 0)
    {
      char buf[200];
      sprintf(buf, "Error while opening file '%s'", fileName);
      error(buf);
    }

  /* Open header to read simulation properties */
  printf("read_ic: Reading runtime parameters...\n");
  h_grp = H5Gopen(h_file, "/RuntimePars");
  if(h_grp < 0)
    error("Error while opening runtime parameters\n");

  /* Read the relevant information */
  readAttribute(h_grp, "PeriodicBoundariesOn", INT, periodic);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to read simulation properties */
  printf("read_ic: Reading file header...\n");
  h_grp = H5Gopen(h_file, "/Header");
  if(h_grp < 0)
    error("Error while opening file header\n");
    
  /* Read the relevant information and print status */
  readAttribute(h_grp, "BoxSize", DOUBLE, &boxSize);
  readAttribute(h_grp, "NumPart_Total", UINT, numParticles);

  *N = numParticles[0];
  dim[0] = dim[1] = dim[2] = boxSize;

  printf("read_ic: Found %d particles in a %speriodic box of size [%f %f %f]\n", 
  	 *N, (periodic ? "": "non-"), dim[0], dim[1], dim[2]);

  /* Close header */
  H5Gclose(h_grp);

  /* Allocate memory to store particles */
  if(posix_memalign( (void*)parts , 32 , *N * sizeof(struct part)) != 0)
    error("Error while allocating memory for particles");
  bzero( *parts , *N * sizeof(struct part) );

  printf("read_ic: Allocated %8.2f MB for particles.\n", *N * sizeof(struct part) / (1024.*1024.));
		  
  /* Open SPH particles group */
  printf("read_ic: Reading particle arrays...\n");
  h_grp = H5Gopen(h_file, "/PartType0");
  if(h_grp < 0)
    error( "Error while opening particle group.\n");

  /* Read arrays */
  readArray(h_grp, "Coordinates", DOUBLE, *N, 3, *parts, x);
  readArray(h_grp, "Velocity", FLOAT, *N, 3, *parts, v);
  readArray(h_grp, "Mass", FLOAT, *N, 1, *parts, mass);
  readArray(h_grp, "SmoothingLength", FLOAT, *N, 1, *parts, h);
  readArray(h_grp, "InternalEnergy", FLOAT, *N, 1, *parts, u);
  readArray(h_grp, "ParticleIDs", ULONGLONG, *N, 1, *parts, id);

  /* Close particle group */
  H5Gclose(h_grp);

  printf("read_ic: Done Reading particles...\n");

  /* Close file */
  H5Fclose(h_file);
}

#endif  /* HAVE_HDF5 */


