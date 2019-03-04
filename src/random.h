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

/* Standard header */
#include <stdlib.h>

/**
 * @brief The categories of random number generated.
 *
 * The values of the fields are carefully chose numbers
 * the numbers are very large primes such that the IDs
 * will not have a prime factorization with this coefficient
 * this results in a very high period for the random number
 * generator.
 * Only change when you know what you are doing, changing
 * the numbers to bad values will break the random number
 * generator.
 * In case new numbers need to be added other possible
 * numbers could be:
 * 4947009007, 5947309451, 6977309513
 */
enum random_number_type {
  random_number_star_formation = 0LL,
  random_number_stellar_feedback = 3947008991LL,
  random_number_stellar_enrichment = 2936881973LL,
  random_number_BH_feedback = 1640531371LL
};

/**
 * @brief Returns a pseudo-random number in the range [0, 1[.
 *
 * We generate numbers that are always reproducible for a given particle ID and
 * simulation time (on the integer time-line). If more than one number per
 * time-step per particle is needed, additional randomness can be obtained by
 * using the type argument.
 *
 * @param id The ID of the particle for which to generate a number.
 * @param ti_current The time (on the time-line) for which to generate a number.
 * @param type The #random_number_type to generate.
 * @return a random number in the interval [0, 1.[.
 */
INLINE static double random_unit_interval(const long long int id,
                                          const integertime_t ti_current,
                                          const enum random_number_type type) {

  /* Range used for the seeds. Best if prime */
  static const long long seed_range = RAND_MAX;
  static const double RAND_MAX_inv = 1. / ((double)RAND_MAX);
  static const long long mwc_number = (1LL << 32) - 1LL;

  /* Calculate the seed */
  /* WARNING: Only change the math if you really know what you are doing!
   * The numbers are carefully chosen prime numbers that prevent correlation
   * with either the current integer time or the particle IDs. The current
   * method also prevents any correlation between different random number
   * types.
   * The calculation overflows on purpose.
   * 1. The first step is calculating the seed by using a multiply with carry
   * (MWC) method, this method depends on the type of random number and
   * this therefore also prevents that there is any correlation between
   * the different types of random numbers.
   * 2. After this we use the 64 bit Xorshift method to randomize the seeds
   * even more.
   * 3. We calculate a prime multiplication for the id with a quadratic
   * term.
   * 4. We calculate the seed by using a Quadratic congruential generator,
   * in which we use the id part and the current time step bin.
   */
  unsigned long long number = ti_current;

  /* Multiply with carry (MWC), (adviced variables by NR) */
  number = 4294957665LL * (number & (mwc_number)) + (number >> 32);

  /* 64-bit Xorshift (adviced variables by NR) */
  number ^= number << 21;
  number ^= number >> 35;
  number ^= number << 4;

  /* Add constant to ID */
  const unsigned long long idt = id + type;

  /* Nonlinear congruential generator */
  const unsigned long long idpart =
      3457LL * idt + 593LL * idt * ti_current + 5417LL * idt * idt;
  unsigned int seed =
      (937LL * number + 5171LL * number * number + idpart + 1109LL) %
      9996361LL % seed_range;

  /* Generate a random number between 0 and 1. */
  return rand_r(&seed) * RAND_MAX_inv;
}

#endif /* SWIFT_RANDOM_H */
