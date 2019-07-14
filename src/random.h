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
#include <stdint.h>
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
 * 5947309451, 6977309513
 */
enum random_number_type {
  random_number_star_formation = 0LL,
  random_number_stellar_feedback = 3947008991LL,
  random_number_stellar_enrichment = 2936881973LL,
  random_number_BH_feedback = 1640531371LL,
  random_number_BH_swallow = 4947009007LL
};

#include <errno.h>
#include <ieee754.h>
#include <limits.h>
#include <sys/types.h>

/* Inline the default RNG functions to avoid costly function calls. */

INLINE static int inl_rand_r(unsigned int *seed) {
  unsigned int next = *seed;
  int result;
  next *= 1103515245;
  next += 12345;
  result = (unsigned int)(next / 65536) % 2048;
  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int)(next / 65536) % 1024;
  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int)(next / 65536) % 1024;
  *seed = next;
  return result;
}

INLINE void inl_drand48_iterate(unsigned short int xsubi[3]) {
  uint64_t X;
  uint64_t result;
  const uint64_t __a = 0x5deece66dull;
  const uint16_t __c = 0xb;

  /* Do the real work.  We choose a data type which contains at least
     48 bits.  Because we compute the modulus it does not care how
     many bits really are computed.  */
  X = (uint64_t)xsubi[2] << 32 | (uint32_t)xsubi[1] << 16 | xsubi[0];
  result = X * __a + __c;
  xsubi[0] = result & 0xffff;
  xsubi[1] = (result >> 16) & 0xffff;
  xsubi[2] = (result >> 32) & 0xffff;
}

INLINE double inl_erand48(unsigned short int xsubi[3]) {
  union ieee754_double temp;

  /* Compute next state.  */
  inl_drand48_iterate(xsubi);

  /* Construct a positive double with the 48 random bits distributed over
     its fractional part so the resulting FP number is [0.0,1.0).  */
  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
  temp.ieee.mantissa1 = ((xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);

  /* Please note the lower 4 bits of mantissa1 are always 0.  */
  return temp.d - 1.0;
}

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
INLINE static double random_unit_interval(long long int id,
                                          const integertime_t ti_current,
                                          const enum random_number_type type) {
  /* Start by packing the state into a sequence of 16-bit seeds for rand_r. */
  uint16_t buff[9];
  id += type;
  memcpy(&buff[0], &id, 8);
  memcpy(&buff[4], &ti_current, 8);
  buff[8] = 6178;

  /* Shuffle the buffer values, this will be our source of entropy for
     the erand48 generator. */
  uint32_t seed16 = 0;
  for (int k = 0; k < 9; k++) {
    seed16 ^= buff[k];
    inl_rand_r(&seed16);
  }
  for (int k = 0; k < 9; k++) buff[k] ^= inl_rand_r(&seed16) & 0xffff;

  /* Do three steps of erand48() over the state generated previously. */
  uint16_t seed48[3] = {0, 0, 0};
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) seed48[j] ^= buff[3 * k + j];
    inl_erand48(seed48);
  }

  /* Generate one final value, this is our output. */
  return inl_erand48(seed48);
}

INLINE static double random_unit_interval_previous(
    const long long int id, const integertime_t ti_current,
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
