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

#ifndef SWIFT_LIGHTCONE_SHELL_H
#define SWIFT_LIGHTCONE_SHELL_H

/* Standard headers */
#include <stdio.h>

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "cosmology.h"
#include "error.h"
#include "lightcone/lightcone_map.h"
#include "lightcone/lightcone_map_types.h"
#include "lightcone/pixel_index.h"
#include "lightcone/projected_kernel.h"
#include "particle_buffer.h"

enum lightcone_shell_state {
  shell_uninitialized,
  shell_current,
  shell_complete,
};

union lightcone_map_buffer_entry {
  int i;
  float f;
};

/**
 * @brief Encode an angle in the range 0 to 2pi as an int
 *
 * @param angle the angle to encode
 */
__attribute__((always_inline)) INLINE static int angle_to_int(
    const double angle) {

  if (angle < 0.0 || angle > 2 * M_PI) error("angle is out of range!");
  const double fac = ((1 << 30) - 1) / M_PI;
  return (int)(angle * fac);
}

/**
 * @brief Convert an encoded angle back to a double
 *
 * @param i the int containing the angle
 */
__attribute__((always_inline)) INLINE static double int_to_angle(const int i) {

  const double fac = M_PI / ((1 << 30) - 1);
  return i * fac;
}

/**
 * @brief Information about a particle type contributing to the lightcone
 *
 * For each Swift particle type we store how many lightcone maps that type
 * contributes to and their indexes in the array of lightcone_maps structs
 * associated with each lightcone_shell.
 *
 * We also record the number of bytes needed to store one update: updates
 * consist of the angular coordinates of the particle, its angular smoothing
 * radius, and the quantities contributed to the lightcone maps.
 *
 */
struct lightcone_particle_type {

  /*! Number of lightcone maps this particle type contributes to */
  int nr_maps;

  /*! Number of smoothed this particle type contributes to */
  int nr_smoothed_maps;

  /*! Number of un-smoothed this particle type contributes to */
  int nr_unsmoothed_maps;

  /*! Indices of the lightcone maps this particle type contributes to.
    Smoothed maps will be stored first in the array. */
  int *map_index;

  /*! Amount of data to store per particle: theta, phi, radius and the value to
   * add to each healpix map */
  size_t buffer_element_size;
};

/**
 * @brief Information about each lightcone shell
 *
 * Each shell contains one lightcone_map for each healpix map
 * we're making. This is where the pixel data is stored while
 * the current simulation timestep overlaps the shell's redshift
 * range.
 *
 * Each shell also contains one particle_buffer per particle type,
 * which stores the updates to be applied to the pixel data.
 * Updates are accumulated in the buffers during each time step
 * and applied at the end of the step.
 *
 */
struct lightcone_shell {

  /*! State of this shell */
  enum lightcone_shell_state state;

  /*! Inner radius of shell */
  double rmin;

  /*! Outer radius of shell */
  double rmax;

  /*! Minimum expansion factor for this shell */
  double amin;

  /*! Maximum expansion factor for this shell */
  double amax;

  /*! Number of maps associated with this shell */
  int nr_maps;

  /*! Array of lightcone maps for this shell */
  struct lightcone_map *map;

  /*! Buffers to store the map updates for each particle type */
  struct particle_buffer buffer[swift_type_count];

  /*! Healpix nside parameter */
  int nside;

  /*! Total pixels in the maps */
  pixel_index_t total_nr_pix;

  /*! Number of pixels per map stored on this node */
  pixel_index_t local_nr_pix;

  /*! Offset of the first pixel stored on this rank */
  pixel_index_t local_pix_offset;

  /*! Number of pixels per rank (last node has any extra) */
  pixel_index_t pix_per_rank;
};

struct lightcone_shell *lightcone_shell_array_init(
    const struct cosmology *cosmo, const char *radius_file, int nr_maps,
    struct lightcone_map_type *map_type, int nside, pixel_index_t total_nr_pix,
    struct lightcone_particle_type *part_type, size_t elements_per_block,
    int *nr_shells_out);

void lightcone_shell_array_free(struct lightcone_shell *shell, int nr_shells);

void lightcone_shell_array_dump(const struct lightcone_shell *shell,
                                int nr_shells, FILE *stream);

struct lightcone_shell *lightcone_shell_array_restore(
    FILE *stream, int nr_shells, struct lightcone_particle_type *part_type,
    size_t elements_per_block);

void lightcone_shell_flush_map_updates(
    struct lightcone_shell *shell, struct threadpool *tp,
    struct lightcone_particle_type *part_type,
    const double max_map_update_send_size_mb,
    struct projected_kernel_table *kernel_table, int verbose);

#endif /* SWIFT_LIGHTCONE_SHELL_H */
