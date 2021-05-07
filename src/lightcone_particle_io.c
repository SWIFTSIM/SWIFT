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
#include "lightcone_particle_io.h"

/* Local headers */
#include "error.h"
#include "gravity.h"
#include "lightcone.h"
#include "particle_buffer.h"

/* Array of output fields */
#define max_fields 20
static int num_fields[swift_type_count];
static struct lightcone_io_props field[swift_type_count][max_fields];


/**
 * @brief Make an array of output fields for each particle type
 */
void lightcone_io_make_output_fields(void) {
  
  for(int i=0; i<swift_type_count; i+=1)
    num_fields[i] = 0;

  int n;

  /* Gas */
#define OFFSET(x) offsetof(struct lightcone_gas_data, x)
  n = 0;
  field[swift_type_gas][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_gas][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_gas][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_gas] = n;
#undef OFFSET

  /* DM */
#define OFFSET(x) offsetof(struct lightcone_dark_matter_data, x)
  n = 0;
  field[swift_type_dark_matter][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_dark_matter][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_dark_matter][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_dark_matter] = n;
#undef OFFSET

  /* DM background (uses same struct as dark matter) */
#define OFFSET(x) offsetof(struct lightcone_dark_matter_data, x)
  n = 0;
  field[swift_type_dark_matter_background][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_dark_matter_background][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_dark_matter_background][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_dark_matter_background] = n;
#undef OFFSET

  /* Stars */
#define OFFSET(x) offsetof(struct lightcone_stars_data, x)
  n = 0;
  field[swift_type_stars][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_stars][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_stars][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_stars] = n;
#undef OFFSET

  /* Black holes */
#define OFFSET(x) offsetof(struct lightcone_black_hole_data, x)
  n = 0;
  field[swift_type_black_hole][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_black_hole][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_black_hole][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_black_hole] = n;
#undef OFFSET

  /* Neutrinos */
#define OFFSET(x) offsetof(struct lightcone_neutrino_data, x)
  n = 0;
  field[swift_type_neutrino][n++] = lightcone_io_make_output_field("ParticleIDs", LONGLONG, 1, OFFSET(id),   UNIT_CONV_NO_UNITS, 0.0);
  field[swift_type_neutrino][n++] = lightcone_io_make_output_field("Coordinates", DOUBLE,   3, OFFSET(x),    UNIT_CONV_LENGTH, 1.0);
  field[swift_type_neutrino][n++] = lightcone_io_make_output_field("Masses",      DOUBLE,   1, OFFSET(mass), UNIT_CONV_MASS, 0.0);
  num_fields[swift_type_neutrino] = n;
#undef OFFSET

}

/*
  Functions to store particle properties in the lightcone_*_data structs.

  These should determine whether the particle should be included in the
  lightcone and, if so, copy the needed quantities into the struct and
  return 1. If the particle should be discarded the function should 
  return 0.
  
 */

/**
 * @brief Store gas properties to write to the lightcone
 */
int lightcone_store_gas(const struct gpart *gp, const struct part *p,
                        const struct xpart *xp, const double a_cross,
                        const double x_cross[3], struct lightcone_gas_data *data) {
  data->id = p->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->mass = p->mass;

  return 1;
}


/**
 * @brief Store dark matter properties to write to the lightcone
 */
int lightcone_store_dark_matter(const struct gpart *gp, const double a_cross, 
                                const double x_cross[3],
                                struct lightcone_dark_matter_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->mass = gp->mass;
  
  return 1;
}


/**
 * @brief Store star properties to write to the lightcone
 */
int lightcone_store_stars(const struct gpart *gp, const struct spart *sp,
                          const double a_cross, const double x_cross[3],
                          struct lightcone_stars_data *data) {
  data->id = sp->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->mass = sp->mass;

  return 1;
}


/**
 * @brief Store black hole properties to write to the lightcone
 */
int lightcone_store_black_hole(const struct gpart *gp, const struct bpart *bp,
                               const double a_cross, const double x_cross[3],
                               struct lightcone_black_hole_data *data) {
  data->id = bp->id;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->mass = bp->mass;

  return 1;
}


/**
 * @brief Store neutrino properties to write to the lightcone
 */
int lightcone_store_neutrino(const struct gpart *gp, const double a_cross,
                             const double x_cross[3], struct lightcone_neutrino_data *data) {
  data->id = gp->id_or_neg_offset;
  data->x[0] = x_cross[0];
  data->x[1] = x_cross[1];
  data->x[2] = x_cross[2];
  data->mass = gp->mass;
  
  return 1;
}


/**
 * @brief Write data to a HDF5 dataset, appending along first axis if it already exists
 */
void append_dataset(const struct unit_system *snapshot_units,
                    enum unit_conversion_factor units, float scale_factor_exponent,
                    hid_t loc_id, const char *name, hid_t mem_type_id, hsize_t chunk_size,
                    const int rank, const hsize_t *dims, const hsize_t num_written,
                    const void *data) {
  
  const int max_rank = 2;
  if(rank > max_rank)error("HDF5 dataset has too may dimensions. Increase max_rank.");
  if(rank < 1)error("HDF5 dataset must be at least one dimensional");

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
  hsize_t chunk_dims[max_rank];
  for(int i=1; i<rank; i+=1)
    chunk_dims[i] = full_dims[i];
  chunk_dims[0] = (hsize_t) chunk_size;

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

    /* Write unit conversion factors for this data set */
    char buffer[FIELD_BUFFER_SIZE] = {0};
    units_cgs_conversion_string(buffer, snapshot_units, units,
                                scale_factor_exponent);
    float baseUnitsExp[5];
    units_get_base_unit_exponents_array(baseUnitsExp, units);
    io_write_attribute_f(dataset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
    io_write_attribute_f(dataset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
    io_write_attribute_f(dataset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
    io_write_attribute_f(dataset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
    io_write_attribute_f(dataset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
    io_write_attribute_f(dataset_id, "h-scale exponent", 0.f);
    io_write_attribute_f(dataset_id, "a-scale exponent", scale_factor_exponent);
    io_write_attribute_s(dataset_id, "Expression for physical CGS units", buffer);

    /* Write the actual number this conversion factor corresponds to */
    const double factor = units_cgs_conversion_factor(snapshot_units, units);
    io_write_attribute_d(dataset_id,
                         "Conversion factor to CGS (not including cosmological corrections)",
                         factor);

    /* Note that we can't write the conversion factor including cosmological corrections
       as an attribute because it will be different for each particle. */

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
  *num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);

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
 * @brief Append buffered particles to the output file.
 */
void lightcone_write_particles(struct lightcone_props *props,
                               const struct unit_system *internal_units,
                               const struct unit_system *snapshot_units,
                               int ptype, hid_t file_id) {
  
  if(num_fields[ptype] > 0) {

    /* Open group and get number and offset of particles to write */
    size_t num_written, num_to_write;
    hid_t group_id = init_write(props, file_id, ptype, &num_written, &num_to_write);
    
    /* Get size of the data struct for this type */
    const size_t data_struct_size = lightcone_io_struct_size(ptype);
      
    /* Loop over output fields */
    for(int field_nr=0; field_nr<num_fields[ptype]; field_nr +=1) {
        
      /* Find output field info */
      struct lightcone_io_props *f = &field[ptype][field_nr];
      hid_t dtype_id = io_hdf5_type(f->type);           /* HDF5 data type */
      size_t type_size = io_sizeof_type(f->type);       /* Bytes per value */
      const size_t field_size = f->dimension*type_size; /* Bytes per particle */

      /* Find unit conversion factor for this quantity */
      const double conversion_factor =
        units_conversion_factor(internal_units, snapshot_units, f->units);

      /* Allocate output buffer */
      char *outbuf = malloc(num_to_write*field_size);
      if(!outbuf)error("Unable to allocate lightcone output buffer");
      char *outptr = outbuf;

      /* Loop over blocks of buffered particles and copy to output array */
      size_t num_elements;
      struct particle_buffer_block *block = NULL;
      char *block_data;
      do {
        particle_buffer_iterate(&props->buffer[ptype], &block, &num_elements, (void **) &block_data);
        for(size_t i=0; i<num_elements; i+=1) {
          char *src = block_data+i*data_struct_size+f->offset;
          char *dest = outptr;
          memcpy(dest, src, field_size);
          outptr += field_size;
        }
      } while(block);

      /* Convert units if necessary */
      if(conversion_factor != 1.0) {
        const size_t nr_values = num_to_write*f->dimension;
        switch(f->type) {
        case INT: {
          int *values = (int *) outbuf;
          for(size_t i=0; i<nr_values;i+=1)
            values[i] *= conversion_factor;
        } break;
        case LONGLONG: {
          long long *values = (long long *) outbuf;
          for(size_t i=0; i<nr_values;i+=1)
            values[i] *= conversion_factor;
        } break;
        case FLOAT: {
          float *values = (float *) outbuf;
          for(size_t i=0; i<nr_values;i+=1)
            values[i] *= conversion_factor;
        } break;
        case DOUBLE: {
          double *values = (double *) outbuf;
          for(size_t i=0; i<nr_values;i+=1)
            values[i] *= conversion_factor;
        } break;
        default:
          error("Unhandled data type");
        }
      }        

      /* Write the data */
      const hsize_t chunk_size = props->hdf5_chunk_size;
      hsize_t dims[] = {(hsize_t) num_to_write, (hsize_t) f->dimension};
      int rank = 1;
      if(f->dimension > 1)rank = 2;
      append_dataset(snapshot_units, f->units, f->scale_factor_exponent,
                     group_id, f->name, dtype_id, chunk_size,
                     rank, dims, num_written, outbuf);
      
      /* Free the output buffer */
      free(outbuf);
      
    } /* Next field */
    
    /* If all fields are done, we can close the particle type group */
    H5Gclose(group_id);
  }
}
