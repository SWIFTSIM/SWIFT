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
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "lightcone/lightcone_map.h"

/* Local headers */
#include "align.h"
#include "common_io.h"
#include "error.h"
#include "memuse.h"
#include "restart.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

void lightcone_map_init(struct lightcone_map *map, const int nside,
                        const pixel_index_t total_nr_pix,
                        const pixel_index_t pix_per_rank,
                        const pixel_index_t local_nr_pix,
                        const pixel_index_t local_pix_offset,
                        const double r_min, const double r_max,
                        struct lightcone_map_type type) {

  /*Store number of pixels in the map etc */
  map->nside = nside;
  map->total_nr_pix = total_nr_pix;
  map->pix_per_rank = pix_per_rank;
  map->local_nr_pix = local_nr_pix;
  map->local_pix_offset = local_pix_offset;

  /* Pixel data is initially not allocated */
  map->data = NULL;

  /* Store resolution parameter, shell size, units etc */
  map->r_min = r_min;
  map->r_max = r_max;
  map->type = type;

  /* Store factor used to retrieve values from the update buffer */
  map->buffer_scale_factor_inv = 1.0 / (type.buffer_scale_factor);

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  /* Initialize total for consistency check */
  map->total = 0.0;
#endif
}

/**
 * @brief Deallocate the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 */
void lightcone_map_clean(struct lightcone_map *map) {

  if (map->data) lightcone_map_free_pixels(map);
}

/**
 * @brief Allocate (and maybe initialize) the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 * @param zero_pixels if true, set allocated pixels to zero
 */
void lightcone_map_allocate_pixels(struct lightcone_map *map,
                                   const int zero_pixels) {

  if (swift_memalign("lightcone_map_pixels", (void **)&map->data,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(double) * map->local_nr_pix) != 0)
    error("Failed to allocate lightcone map pixel data");

  if (zero_pixels) {
    for (pixel_index_t i = 0; i < map->local_nr_pix; i += 1) map->data[i] = 0.0;
  }
}

void lightcone_map_free_pixels(struct lightcone_map *map) {

  swift_free("lightcone_map_pixels", (void *)map->data);
  map->data = NULL;
}

/**
 * @brief Dump lightcone_map struct to the output stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to write to.
 */
void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream) {

  /* Write the struct */
  restart_write_blocks((void *)map, sizeof(struct lightcone_map), 1, stream,
                       "lightcone_map", "lightcone_map");

  /* Write the pixel data if it is allocated */
  if (map->data)
    restart_write_blocks((void *)map->data, sizeof(double), map->local_nr_pix,
                         stream, "lightcone_map_data", "lightcone_map_data");
}

/**
 * @brief Restore lightcone_map struct from the input stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to read from.
 */
void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream) {

  /* Read the struct */
  restart_read_blocks((void *)map, sizeof(struct lightcone_map), 1, stream,
                      NULL, "lightcone_map");

  /* Read the pixel data if it was allocated.
     map->data from the restart file is not a valid pointer now but we can
     check if it is not null to see if the pixel data block was written out. */
  if (map->data) {
    lightcone_map_allocate_pixels(map, /* zero_pixels = */ 0);
    restart_read_blocks((void *)map->data, sizeof(double), map->local_nr_pix,
                        stream, NULL, "lightcone_map");
  }
}

#ifdef HAVE_HDF5
/**
 * @brief Write a lightcone map to a HDF5 file
 *
 * @param map the #lightcone_map structure
 * @param loc a HDF5 file or group identifier to write to
 * @param name the name of the dataset to create
 */
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id,
                         const char *name,
                         const struct unit_system *internal_units,
                         const struct unit_system *snapshot_units,
                         const int collective, const int gzip_level,
                         const int chunk_size,
                         enum lossy_compression_schemes compression) {

#ifdef WITH_MPI
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

  /* Find unit conversion factor for this quantity */
  const double map_conversion_factor =
      units_conversion_factor(internal_units, snapshot_units, map->type.units);

  /* Convert units of pixel data if necessary */
  if (map_conversion_factor != 1.0) {
    for (pixel_index_t i = 0; i < map->local_nr_pix; i += 1)
      map->data[i] *= map_conversion_factor;
  }

  /* Get conversion factor for shell radii */
  const double length_conversion_factor =
      units_conversion_factor(internal_units, snapshot_units, UNIT_CONV_LENGTH);

  /* Create dataspace in memory corresponding to local pixels */
  const hsize_t mem_dims[1] = {(hsize_t)map->local_nr_pix};
  hid_t mem_space_id = H5Screate_simple(1, mem_dims, NULL);
  if (mem_space_id < 0) error("Unable to create memory dataspace");

  /* Create dataspace corresponding to the part of the map in this file. */
  hsize_t file_dims[1];
  if (collective) {
    /* For collective writes the file contains the full map */
    file_dims[0] = (hsize_t)map->total_nr_pix;
  } else {
    /* For distributed writes the file contains the local pixels only */
    file_dims[0] = (hsize_t)map->local_nr_pix;
  }
  hid_t file_space_id = H5Screate_simple(1, file_dims, NULL);
  if (file_space_id < 0) error("Unable to create file dataspace");

    /* Select the part of the dataset in the file to write to */
#ifdef WITH_MPI
#ifdef HAVE_PARALLEL_HDF5
  if (collective) {
    const pixel_index_t pixel_offset = map->local_pix_offset;
    const hsize_t start[1] = {(hsize_t)pixel_offset};
    const hsize_t count[1] = {(hsize_t)map->local_nr_pix};
    if (H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, NULL, count,
                            NULL) < 0)
      error("Unable to select part of file dataspace to write to");
  }
#else
  if (collective)
    error("Writing lightcone maps with collective I/O requires parallel HDF5");
#endif
#endif

  /* Property list for creating the dataset */
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);

  /* Data type to write in the file */
  hid_t dtype_id = H5Tcopy(H5T_NATIVE_DOUBLE);

  /* In non-collective mode we may want to enable compression, which requires
     a chunked dataset layout. If no compression filters are in use we use
     a contiguous layout. */
  if (!collective) {

    /* Set lossy compression, if requested. This might change the
       output data type and add to the property list. */
    char filter_name[32];
    set_hdf5_lossy_compression(&prop_id, &dtype_id, compression, name,
                               filter_name);

    /* Set lossless compression */
    if (gzip_level > 0) {
      H5Pset_shuffle(prop_id);
      if (H5Pset_deflate(prop_id, gzip_level) < 0)
        error("Unable to set HDF5 deflate filter for healpix map");
    }

    /* Set the chunk size, but only if we're using filters */
    if (H5Pget_nfilters(prop_id) > 0) {
      hsize_t dim[1];
      if (map->local_nr_pix > chunk_size) {
        dim[0] = (hsize_t)chunk_size;
      } else {
        dim[0] = (hsize_t)map->local_nr_pix;
      }
      if (H5Pset_chunk(prop_id, 1, dim) < 0)
        error("Unable to set HDF5 chunk size for healpix map");
    }
  }

  /* Create the dataset */
  hid_t dset_id = H5Dcreate(loc_id, name, dtype_id, file_space_id, H5P_DEFAULT,
                            prop_id, H5P_DEFAULT);
  H5Pclose(prop_id);
  H5Tclose(dtype_id);
  if (dset_id < 0) error("Unable to create dataset %s", name);

  /* Write attributes */
  io_write_attribute_i(dset_id, "nside", map->nside);
  io_write_attribute_ll(dset_id, "number_of_pixels",
                        (long long)map->total_nr_pix);
  io_write_attribute_s(dset_id, "pixel_ordering_scheme", "ring");
  io_write_attribute_d(dset_id, "comoving_inner_radius",
                       map->r_min * length_conversion_factor);
  io_write_attribute_d(dset_id, "comoving_outer_radius",
                       map->r_max * length_conversion_factor);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, map->type.units, 0.f);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, map->type.units);
  io_write_attribute_f(dset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset_id, "h-scale exponent", 0.f);
  io_write_attribute_f(dset_id, "a-scale exponent", 0.f);
  io_write_attribute_s(dset_id, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double cgs_factor =
      units_cgs_conversion_factor(snapshot_units, map->type.units);
  io_write_attribute_d(
      dset_id,
      "Conversion factor to CGS (not including cosmological corrections)",
      cgs_factor);

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
  /* Consistency check: will write out expected sum over pixels */
  double total = map->total;
#ifdef WITH_MPI
  MPI_Allreduce(&map->total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  total *= map_conversion_factor;
  io_write_attribute_d(dset_id, "expected_sum", total);
#endif

  /* Set up property list for the write */
  hid_t h_plist_id = H5Pcreate(H5P_DATASET_XFER);
#if defined(WITH_MPI)
#ifdef HAVE_PARALLEL_HDF5
  if (collective) {
    if (H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE) < 0)
      error("Unable to set collective transfer mode");
  }
#else
  if (collective)
    error("Writing lightcone maps with MPI requires parallel HDF5");
#endif
#endif

  /* Write the data */
  if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id,
               h_plist_id, map->data) < 0)
    error("Unable to write dataset %s", name);

  /* Tidy up */
  H5Dclose(dset_id);
  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);
  H5Pclose(h_plist_id);
}
#endif /* HAVE_HDF5*/
