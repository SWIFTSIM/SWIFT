/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2019 Folkert Nobels    (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/
#ifndef SWIFT_RANDOM_H
#define SWIFT_RANDOM_H

/* COde configuration */
#include "../config.h"

/* Random 123 headers */
#include "Random123/threefry.h"

/* Standard header */
#include <stdlib.h>
#include <limits.h>

/**
 * @brief The categories of random number generated.
 *
 * The values of the fields are carefully chose numbers
 * but it is not super important what the values are, 
 * as long as for the different types the numbers are 
 * integers and not the same
 */
enum random_number_type {
  random_number_star_formation = 0,
  random_number_stellar_feedback = 1, 
  random_number_stellar_enrichment = 2,
  random_number_BH_feedback = 3,
  random_number_BH_swallow = 4,
  random_number_stellar_winds = 5,
  random_number_SNIa_feedback = 6,
  random_number_HII_regions = 7
};

/**
 * @brief Returns a pseudo-random number in the range (0, 1).
 *
 * We generate numbers that are always reproducible for a given particle ID and
 * simulation time (on the integer time-line). If more than one number per
 * time-step per particle is needed, additional randomness can be obtained by
 * using the type argument. This random number generator has a very small spacing
 * of 2^-64, and therefore can be used for processes with a very low propability
 *
 * @param id The ID of the particle for which to generate a number.
 * @param ti_current The time (on the time-line) for which to generate a number.
 * @param type The #random_number_type to generate.
 * @return a random number in the open interval (0, 1.).
 */
INLINE static double random_unit_interval(const long long int id,
                                          const integertime_t ti_current,
                                          const enum random_number_type type) {

  /* Default used limits */
  static const unsigned long long int max_long = ULLONG_MAX;
  static const double max_long_inv = 1. / max_long;
  static const int number_of_rounds = 12;
  
  /* Prepare the state, the initial state depends on the type */
  threefry4x64_ctr_t ctr = {{type,0}};
  threefry4x64_key_t key = {{id, ti_current}};
    
  /* Generate the random number */
  threefry4x64_ctr_t rand = threefry4x64_R(number_of_rounds,ctr, key);

  return rand.v[0] * max_long_inv;
  
}


#endif /* SWIFT_RANDOM_H */
