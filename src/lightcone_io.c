/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Some standard headers. */
#include <hdf5.h>

/* This object's header. */
#include "lightcone_io.h"

/* Local headers */
#include "error.h"
#include "gravity.h"
#include "lightcone.h"
#include "particle_buffer.h"


/**
 * @brief Write data to a HDF5 dataset, appending along first axis if it already exists
 */
void append_dataset(hid_t loc_id, const char *name, hid_t mem_type_id, const int rank,
                    const hsize_t *dims, const hsize_t num_written, const void *data) {
  
  const int max_rank = 2;
  if(rank > max_rank)error("HDF5 dataset has too may dimensions. Increase max_rank.");

  /* If we have zero elements to append, there's nothing to do */
  if(dims[0] == 0)return;

  /* Determine size of the dataset after we append our data */
  hsize_t full_dims[max_rank];
  for(int i=0; i<rank; i+=1)
    full_dims[i] = dims[i];
  full_dims[0] += num_written;

  /* Determine maximum size in each dimension */
  hsize_t max_dims[max_rank];
  for(int i=1; i<rank; i+=1)
    max_dims[i] = full_dims[i];
  max_dims[0] = H5S_UNLIMITED;

  /* Determine chunk size in each dimension */
  const hsize_t chunk_size = 102400;
  hsize_t chunk_dims[max_rank];
  for(int i=1; i<rank; i+=1)
    chunk_dims[i] = full_dims[i];
  chunk_dims[0] = chunk_size;

  /* Find offset to region to write in each dimension */
  hsize_t offset[max_rank];
  for(int i=1; i<rank; i+=1)
    offset[i] = 0;
  offset[0] = num_written;

  hid_t dataset_id;
  hid_t file_space_id;
  if(num_written==0) {

    /* We need to create a new dataset */
    file_space_id = H5Screate_simple(rank, full_dims, max_dims);    
    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(prop_id, rank, chunk_dims);
    dataset_id = H5Dcreate(loc_id, name, mem_type_id, file_space_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);
    if(dataset_id < 0)error("Failed to create new dataset: %s", name);
    H5Pclose(prop_id);

  } else {

    /* We're appending to an existing dataset */
    dataset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    if(dataset_id < 0)error("Failed to open existing dataset: %s", name);
    if(H5Dset_extent(dataset_id, full_dims) < 0)error("Unable to extend dataset: %s", name);
    file_space_id = H5Dget_space(dataset_id);

  }

  /* Create memory dataspace */
  hid_t mem_space_id = H5Screate_simple(rank, dims, NULL);

  /* Select region to write in the file */
  if(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset, NULL, dims, NULL) < 0)
    error("Failed to select region in dataset: %s", name);

  /* Write the data */
  if(H5Dwrite(dataset_id, mem_type_id, mem_space_id, file_space_id, H5P_DEFAULT, data) < 0)
    error("Failed to write dataset: %s", name);

  /* Clean up*/
  H5Sclose(file_space_id);
  H5Sclose(mem_space_id);
  H5Dclose(dataset_id);
}


hid_t init_write(struct lightcone_props *props, hid_t file_id, int ptype,
                 size_t *num_written, size_t *num_to_write) {
  
  /* Number of particles already written to the file */
  *num_written = props->num_particles_written_to_file[ptype];
  
  /* Number of buffered particles */
  *num_to_write = props->buffer[ptype].total_num_elements;

  /* Create or open the HDF5 group for this particle type */
  const char *name = part_type_names[ptype];
  hid_t group_id;
  if(*num_written > 0) {
    group_id = H5Gopen(file_id, name, H5P_DEFAULT);
    if(group_id < 0)error("Failed to open existing group: %s", name);
  } else {
    group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(group_id < 0)error("Failed to create new group: %s", name);
  }
  return group_id;
}

/**
 * @brief Store gas properties to write to the lightcone
 */
void lightcone_store_gas(const struct gpart *gp, const struct part *p,
                         const struct xpart *xp, const double x_cross[3],
                         struct lightcone_gas_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
}


/**
 * @brief Append buffered gas particles to the output file.
 */
void lightcone_write_gas(struct lightcone_props *props, hid_t file_id,
                         int ptype) {

  /* Open group and get number and offset of particles to write */
  size_t num_written, num_to_write;
  hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);

  /* Allocate output arrays */
  double *pos = malloc(3*num_to_write*sizeof(double));
  long long *id = malloc(num_to_write*sizeof(long long));

  /* Loop over blocks of buffered particles and copy to output arrays */
  size_t num_elements;
  size_t offset = 0;
  struct particle_buffer_block *block = NULL;
  struct lightcone_gas_data *data;
  do {
    particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &data);
    for(size_t i=0; i<num_elements; i+=1) {
      id[offset]      = data[i].id;
      pos[3*offset+0] = data[i].x[0];
      pos[3*offset+1] = data[i].x[1];
      pos[3*offset+2] = data[i].x[2];
      offset += 1;
    }
  } while(block);

  /* Write the data */
  hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) 3};
  append_dataset(group_id, "Coordinates", H5T_NATIVE_DOUBLE, 2, dims, num_written, pos);
  append_dataset(group_id, "ParticleIDs", H5T_NATIVE_LLONG, 1, dims, num_written, id);

  /* Clean up */
  free(pos);
  free(id);
  H5Gclose(group_id);
}


/**
 * @brief Store dark matter properties to write to the lightcone
 */
void lightcone_store_dark_matter(const struct gpart *gp, const double x_cross[3],
                                 struct lightcone_dark_matter_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
}


/**
 * @brief Append buffered dark matter particles to the output file.
 */
void lightcone_write_dark_matter(struct lightcone_props *props, hid_t file_id,
                                 int ptype) {

  /* Open group and get number and offset of particles to write */
  size_t num_written, num_to_write;
  hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);

  /* Allocate output arrays */
  double *pos = malloc(3*num_to_write*sizeof(double));
  long long *id = malloc(num_to_write*sizeof(long long));

  /* Loop over blocks of buffered particles and copy to output arrays */
  size_t num_elements;
  size_t offset = 0;
  struct particle_buffer_block *block = NULL;
  struct lightcone_dark_matter_data *data;
  do {
    particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &data);
    for(size_t i=0; i<num_elements; i+=1) {
      id[offset]      = data[i].id;
      pos[3*offset+0] = data[i].x[0];
      pos[3*offset+1] = data[i].x[1];
      pos[3*offset+2] = data[i].x[2];
      offset += 1;
    }
  } while(block);

  /* Write the data */
  hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) 3};
  append_dataset(group_id, "Coordinates", H5T_NATIVE_DOUBLE, 2, dims, num_written, pos);
  append_dataset(group_id, "ParticleIDs", H5T_NATIVE_LLONG, 1, dims, num_written, id);

  /* Clean up */
  free(pos);
  free(id);
  H5Gclose(group_id);
}


/**
 * @brief Store star properties to write to the lightcone
 */
void lightcone_store_stars(const struct gpart *gp, const struct spart *sp,
                           const double x_cross[3],
                           struct lightcone_stars_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
}


/**
 * @brief Append buffered star particles to the output file.
 */
void lightcone_write_stars(struct lightcone_props *props, hid_t file_id,
                           int ptype) {

  /* Open group and get number and offset of particles to write */
  size_t num_written, num_to_write;
  hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);

  /* Allocate output arrays */
  double *pos = malloc(3*num_to_write*sizeof(double));
  long long *id = malloc(num_to_write*sizeof(long long));

  /* Loop over blocks of buffered particles and copy to output arrays */
  size_t num_elements;
  size_t offset = 0;
  struct particle_buffer_block *block = NULL;
  struct lightcone_stars_data *data;
  do {
    particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &data);
    for(size_t i=0; i<num_elements; i+=1) {
      id[offset]      = data[i].id;
      pos[3*offset+0] = data[i].x[0];
      pos[3*offset+1] = data[i].x[1];
      pos[3*offset+2] = data[i].x[2];
      offset += 1;
    }
  } while(block);

  /* Write the data */
  hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) 3};
  append_dataset(group_id, "Coordinates", H5T_NATIVE_DOUBLE, 2, dims, num_written, pos);
  append_dataset(group_id, "ParticleIDs", H5T_NATIVE_LLONG, 1, dims, num_written, id);

  /* Clean up */
  free(pos);
  free(id);
  H5Gclose(group_id);
}


/**
 * @brief Store black hole properties to write to the lightcone
 */
void lightcone_store_black_hole(const struct gpart *gp, const struct bpart *bp,
                                const double x_cross[3],
                                struct lightcone_black_hole_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
}


/**
 * @brief Append buffered black hole particles to the output file.
 */
void lightcone_write_black_hole(struct lightcone_props *props, hid_t file_id,
                           int ptype) {

  /* Open group and get number and offset of particles to write */
  size_t num_written, num_to_write;
  hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);

  /* Allocate output arrays */
  double *pos = malloc(3*num_to_write*sizeof(double));
  long long *id = malloc(num_to_write*sizeof(long long));

  /* Loop over blocks of buffered particles and copy to output arrays */
  size_t num_elements;
  size_t offset = 0;
  struct particle_buffer_block *block = NULL;
  struct lightcone_black_hole_data *data;
  do {
    particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &data);
    for(size_t i=0; i<num_elements; i+=1) {
      id[offset]      = data[i].id;
      pos[3*offset+0] = data[i].x[0];
      pos[3*offset+1] = data[i].x[1];
      pos[3*offset+2] = data[i].x[2];
      offset += 1;
    }
  } while(block);

  /* Write the data */
  hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) 3};
  append_dataset(group_id, "Coordinates", H5T_NATIVE_DOUBLE, 2, dims, num_written, pos);
  append_dataset(group_id, "ParticleIDs", H5T_NATIVE_LLONG, 1, dims, num_written, id);

  /* Clean up */
  free(pos);
  free(id);
  H5Gclose(group_id);
}


/**
 * @brief Store neutrino properties to write to the lightcone
 */
void lightcone_store_neutrino(const struct gpart *gp, const double x_cross[3],
                              struct lightcone_neutrino_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
}


/**
 * @brief Append buffered black hole particles to the output file.
 */
void lightcone_write_neutrino(struct lightcone_props *props, hid_t file_id,
                           int ptype) {

  /* Open group and get number and offset of particles to write */
  size_t num_written, num_to_write;
  hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);

  /* Allocate output arrays */
  double *pos = malloc(3*num_to_write*sizeof(double));
  long long *id = malloc(num_to_write*sizeof(long long));

  /* Loop over blocks of buffered particles and copy to output arrays */
  size_t num_elements;
  size_t offset = 0;
  struct particle_buffer_block *block = NULL;
  struct lightcone_neutrino_data *data;
  do {
    particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &data);
    for(size_t i=0; i<num_elements; i+=1) {
      id[offset]      = data[i].id;
      pos[3*offset+0] = data[i].x[0];
      pos[3*offset+1] = data[i].x[1];
      pos[3*offset+2] = data[i].x[2];
      offset += 1;
    }
  } while(block);

  /* Write the data */
  hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) 3};
  append_dataset(group_id, "Coordinates", H5T_NATIVE_DOUBLE, 2, dims, num_written, pos);
  append_dataset(group_id, "ParticleIDs", H5T_NATIVE_LLONG, 1, dims, num_written, id);

  /* Clean up */
  free(pos);
  free(id);
  H5Gclose(group_id);
}
