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
#include "io_types.h"


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
      char buf[100];
      sprintf(buf, "Error while checking the existence of data set '%s'\n", name);
      error(buf);
    }
  else if(exist == 0)
    {
      if(importance == COMPULSORY)
	{
	  char buf[100];
	  sprintf(buf, "Compulsory data set '%s' not present in the file.\n", name);
	  error(buf);
	}
      else
	{
	  printf("readArray: Optional data set '%s' not present. Zeroing this particle field...\n", name);	  
	  
	  for(i=0; i<N; ++i)
	    memset(part_c+i*partSize, 0, copySize);
	  
	  return;
	}
   }

  printf("readArray: Reading %s '%s' array...\n", importance == COMPULSORY ? "compulsory": "optional  ", name);

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
  readArray(h_grp, "Coordinates", DOUBLE, *N, 3, *parts, x, COMPULSORY);
  readArray(h_grp, "Velocity", FLOAT, *N, 3, *parts, v, COMPULSORY);
  readArray(h_grp, "Mass", FLOAT, *N, 1, *parts, mass, COMPULSORY);
  readArray(h_grp, "SmoothingLength", FLOAT, *N, 1, *parts, h, COMPULSORY);
  readArray(h_grp, "InternalEnergy", FLOAT, *N, 1, *parts, u, COMPULSORY);
  readArray(h_grp, "ParticleIDs", ULONGLONG, *N, 1, *parts, id, COMPULSORY);
  readArray(h_grp, "TimeStep", FLOAT, *N, 1, *parts, dt, OPTIONAL);
  readArray(h_grp, "Acceleration", FLOAT, *N, 3, *parts, a, OPTIONAL);
  readArray(h_grp, "Density", FLOAT, *N, 1, *parts, rho, OPTIONAL );
  readArray(h_grp, "Acceleration", FLOAT, *N, 3, *parts, a, OPTIONAL );

  /* Close particle group */
  H5Gclose(h_grp);

  printf("read_ic: Done Reading particles...\n");

  /* Close file */
  H5Fclose(h_file);
}

#endif  /* HAVE_HDF5 */


