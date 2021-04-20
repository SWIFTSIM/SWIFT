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

#ifndef SWIFT_LIGHTCONE_MAP_H
#define SWIFT_LIGHTCONE_MAP_H

/* Config parameters. */
#include "../config.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* Local headers */
#include "particle_buffer.h"


/**
 * @brief Struct to store a contribution to a healpix map pixel
 */
struct lightcone_map_contribution {

  /*! Healpix pixel identifier */
  size_t pixel;
  
  /*! Amount to contribute to the pixel */
  double value;
};


/**
 * @brief Struct to store a single lightcone healpix map
 */
struct lightcone_map {

  /*! Buffer to store contributions to the healpix map */
  struct particle_buffer buffer;

  /*! Block size in the particle buffer */
  size_t elements_per_block;

  /*! Total pixels in the map */
  size_t total_nr_pix;

  /*! Number of pixels stored on this node */
  size_t local_nr_pix;
  
  /*! Number of pixels per rank (last node has any extra) */
  size_t pix_per_rank;

  /*! Local healpix map data */
  double *data;

};


/**
 * @brief Add a value to the buffer for a healpix map
 *
 * @param map the #lightcone_map to update
 * @param pixel the pixel index to update
 * @param value the value to add
 *
 */
__attribute__((always_inline)) INLINE static void lightcone_map_buffer_update(struct lightcone_map *map,
                                                                              const size_t pixel,
                                                                              const double value) {
  struct lightcone_map_contribution contr;
  contr.pixel = pixel;
  contr.value = value;
  particle_buffer_append(&map->buffer, &contr);
}


void lightcone_map_init(struct lightcone_map *map, const int nside,
                        const size_t elements_per_block);

void lightcone_map_clean(struct lightcone_map *map);

void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream);

void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream);

void lightcone_map_allocate_pixels(struct lightcone_map *map, const int zero_pixels);

void lightcone_map_free_pixels(struct lightcone_map *map);

void lightcone_map_update_from_buffer(struct lightcone_map *map);

#ifdef HAVE_HDF5
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id, const char *name);
#endif

#endif /* #ifndef SWIFT_LIGHTCONE_MAP_H */
