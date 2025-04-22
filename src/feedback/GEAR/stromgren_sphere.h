/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_STROMGREN_SPHERE_GEAR_H
#define SWIFT_STROMGREN_SPHERE_GEAR_H

#include <float.h>

#include "inline.h"

#define GEAR_STROMGREN_NUMBER_NEIGHBOURS 50

/**
 * @brief Fields used in gas HII ionization feedback
 */
struct stromgren_shell_data {

  /* Distance */
  float distance;

  /* Ionizing photon rate to fully ionize this particle. This is *huge* number
     --> float is not enough */
  double Delta_N_dot;
};


/**
 * @brief For a given feedback scheme, sets all fields in the stromgren_shell
 * array at their default values
 *
 * @param stromgren An array of ray structs
 * @param max_number_of_rays Maximum number of rays in the feedback scheme
 */
__attribute__((always_inline)) INLINE static void stromgren_shell_init(
    struct stromgren_shell_data* stromgren, const int max_number_of_part) {

  /* Set all fields in the ray struct at their default values */
  for (int i = 0; i < max_number_of_part; i++) {
    stromgren[i].distance = FLT_MAX;
    stromgren[i].Delta_N_dot = 0;
  }
}

/**
 * @brief Minimises the distance between the star and its neighbours. In the
 * resulting ray array, the gas neighbours will be sorted based on their
 * separation from the star. The left-most element will store the information
 * of the particle with the smallest distance to the star, and the right-most
 * with the largest.
 *
 * @param r Comoving distance between the two particles
 * @param ray Ray data
 * @param N_shell_arr Size of stromgren_shell array
 * @param gas_part_id ID of the gas particle
 * @param m Gas particle mass
 */
__attribute__((always_inline)) INLINE static void stromgren_sort_distance(
    const float r, struct stromgren_shell_data* stromgren, const int N_shell_arr,
    const double Delta_N_dot) {

  int insert_index = -1;

  /* Find the first element with a distance greater than r.
     Loop from left to right. Find the left-most element whose  current min
     length is larger than r. Note that N_ray_arr is equal to min(how many
     times this function has been called at this time-step, the maximum number
     of rays per particle) */
  for (int i = 0; i < N_shell_arr; i++) {
    if (r < stromgren[i].distance) {
      insert_index = i;
      break;
    }
  }

  /* If found something to update */
  if (insert_index != -1) {

    /* Shift all elements from the insertion point rightward by one. (Loop form
       right to left). */
    for (int i = N_shell_arr - 2; i >= insert_index; i--) {
      stromgren[i + 1] = stromgren[i];
    }

    /* Insert the new particle's data */
    stromgren[insert_index].distance = r;
    stromgren[insert_index].Delta_N_dot = Delta_N_dot;
  }
}

#endif /* SWIFT_STROMGREN_SPHERE_GEAR_H */
