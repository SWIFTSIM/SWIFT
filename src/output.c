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
      char buf[100];
      sprintf(buf, "Error while creating dataspace for attribute '%s'\n", name);
      error(buf);
    }

  h_err = H5Sset_extent_simple(h_space, 1, dim, NULL);
  if(h_err < 0)
    {
      char buf[100];
      sprintf(buf, "Error while changing dataspace shape for attribute '%s'\n", name);
      error(buf);
    }

  h_attr = H5Acreate(grp, name, hdf5Type(type), h_space, H5P_DEFAULT);
  if(h_attr < 0)
    {
      char buf[100];
      sprintf(buf, "Error while creating attribute '%s'\n", name);
      error(buf);
    }

  h_err = H5Awrite(h_attr, hdf5Type(type), data);
  if(h_err < 0)
    {
      char buf[100];
      sprintf(buf, "Error while reading attribute '%s'\n", name);
      error(buf);
    }

  H5Sclose(h_space);
  H5Aclose(h_attr);
}



/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurence of the field of interest in the parts array
 *
 * @todo A better version using HDF5 hyperslabs to write the file directly from the part array
 * will be written once the strucutres have been stabilized.
 *
 * Calls #error() if an error occurs.
 */
void writeArrayBackEnd(hid_t grp, char* name, enum DATA_TYPE type, int N, int dim, char* part_c)
{
  hid_t h_data=0, h_err=0, h_space=0;
  void* temp = 0;
  int i=0, rank=0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  const size_t partSize = sizeof(struct part);
  char* temp_c = 0;
  hsize_t shape[2];

  printf("writeArray: Writing '%s' array...\n", name);

  /* Allocate temporary buffer */
  temp = malloc(N * dim * sizeOfType(type));
  if(temp == NULL)
    error("Unable to allocate memory for temporary buffer");

  /* Copy temporary buffer to particle data */
  temp_c = temp;
  for(i=0; i<N; ++i)
    memcpy(&temp_c[i*copySize], part_c+i*partSize, copySize);

  /* Create data space */
  h_space = H5Screate(H5S_SIMPLE);
  if(h_space < 0)
    {
      char buf[100];
      sprintf(buf, "Error while creating data space for field '%s'\n", name);
      error(buf);      
    }
  
  if(dim > 1)
    {
      rank = 2;
      shape[0] = N; shape[1] = dim;
    }
  else
    {
      rank = 1;
      shape[0] = N; shape[1] = 0;
    }
  
  /* Change shape of data space */
  h_err = H5Sset_extent_simple(h_space, rank, shape, NULL);
  if(h_err < 0)
    {
      char buf[100];
      sprintf(buf, "Error while changing data space shape for field '%s'\n", name);
      error(buf);      
    }
  
  /* Create dataset */
  h_data = H5Dcreate(grp, name, hdf5Type(type), h_space, H5P_DEFAULT);
  if(h_data < 0)
    {
      char buf[100];
      sprintf(buf, "Error while creating dataspace '%s'\n", name);
      error(buf);
    }
  
  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, hdf5Type(type), h_space, H5S_ALL, H5P_DEFAULT, temp);
  if(h_err < 0)
    {
      char buf[100];
      sprintf(buf, "Error while reading data array '%s'\n", name);
      error(buf);
    }
  
  /* Free and close everything */
  free(temp);
  H5Dclose(h_data);
  H5Sclose(h_space);
}

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group in which to write.
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurence of the field of interest in the parts array
 *
 */
#define writeArray(grp, name, type, N, dim, part, field) writeArrayBackEnd(grp, name, type, N, dim, (char*)(&(part[0]).field))

/**
 * @brief Writes an HDF5 output file (GADGET-3 type)
 *
 * @param fileName The file to write.
 * @param dim The dimension of the volume written to the file.
 * @param parts The array of #part to write in the file.
 * @param N The number of particles to write.
 * @param periodic 1 if the volume is periodic, 0 if not.
 *
 * Creates the HDF5 file fileName and writess the particles contained
 * in the parts array. If such a file already exists, it is erased and replaced
 * by the new one.
 *
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output ( char* fileName, double dim[3], struct part *parts,  int N, int periodic)
{
  hid_t h_file=0, h_grp=0;
  int numParticles[6]={N,0};

  /* Open file */
  printf("write_output: Opening file '%s'.\n", fileName);
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  if(h_file < 0)
    {
      char buf[200];
      sprintf(buf, "Error while opening file '%s'", fileName);
      error(buf);
    }

  /* Open header to read simulation properties */
  printf("write_output: Writing runtime parameters...\n");
  h_grp = H5Gcreate(h_file, "/RuntimePars", 0);
  if(h_grp < 0)
    error("Error while creating runtime parameters group\n");

  /* Write the relevant information */
  writeAttribute(h_grp, "PeriodicBoundariesOn", INT, &periodic, 1);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  (void) sizeOfType(INT);
  
  /* Open header to write simulation properties */
  printf("write_output: Writing file header...\n");
  h_grp = H5Gcreate(h_file, "/Header", 0);
  if(h_grp < 0)
    error("Error while creating file header\n");
    
  /* Read the relevant information and print status */
  writeAttribute(h_grp, "BoxSize", DOUBLE, &dim[0], 1);
  writeAttribute(h_grp, "NumPart_Total", UINT, numParticles, 6);

  /* Close header */
  H5Gclose(h_grp);
		  
  /* Create SPH particles group */
  printf("write_output: Writing particle arrays...\n");
  h_grp = H5Gcreate(h_file, "/PartType0", 0);
  if(h_grp < 0)
    error( "Error while creating particle group.\n");

  /* Write arrays */
  writeArray(h_grp, "Coordinates", DOUBLE, N, 3, parts, x);
  writeArray(h_grp, "Velocity", FLOAT, N, 3, parts, v);
  writeArray(h_grp, "Mass", FLOAT, N, 1, parts, mass);
  writeArray(h_grp, "SmoothingLength", FLOAT, N, 1, parts, h);
  writeArray(h_grp, "InternalEnergy", FLOAT, N, 1, parts, u);
  writeArray(h_grp, "ParticleIDs", ULONGLONG, N, 1, parts, id);
  writeArray(h_grp, "TimeStep", FLOAT, N, 1, parts, dt);
  writeArray(h_grp, "Acceleration", FLOAT, N, 3, parts, a);
  writeArray(h_grp, "Density", FLOAT, N, 1, parts, rho);

  /* Close particle group */
  H5Gclose(h_grp);

  printf("write_output: Done writing particles...\n");

  /* Close file */
  H5Fclose(h_file);
}

#endif  /* HAVE_HDF5 */


