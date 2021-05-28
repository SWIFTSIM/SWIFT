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

/* Code configuration */
#include "../config.h"

/* Standard header */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

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
 * 65610001
 * 126247697
 * 193877777
 * 303595777
 * 384160001
 * 406586897
 * 562448657
 */
enum random_number_type {
  random_number_star_formation = 0LL,
  random_number_sink_formation = 5947309451LL,
  random_number_stellar_feedback_1 = 3947008991LL,
  random_number_stellar_feedback_2 = 6977309513LL,
  random_number_stellar_feedback_3 = 9762399103LL,
  random_number_isotropic_SNII_feedback_ray_theta = 3298327511LL,
  random_number_isotropic_SNII_feedback_ray_phi = 6311114273LL,
  random_number_isotropic_AGN_feedback_ray_theta = 8899891613LL,
  random_number_isotropic_AGN_feedback_ray_phi = 10594523341LL,
  random_number_stellar_enrichment = 2936881973LL,
  random_number_BH_feedback = 1640531371LL,
  random_number_BH_swallow = 4947009007LL,
  random_number_BH_reposition = 59969537LL,
  random_number_SIDM = 2189989727LL,
  random_number_SIDM_theta = 1640531371LL,
  random_number_SIDM_phi = 4947009007LL
};

#ifndef __APPLE__

#include <errno.h>
#include <ieee754.h>
#include <limits.h>

/* Inline the default RNG functions to avoid costly function calls. These
   functions are minor modifications, but functional equivalents, of their glibc
   counterparts. */

INLINE static int inl_rand_r(uint32_t *seed) {
  uint32_t next = *seed;
  int result;
  next *= 1103515245;
  next += 12345;
  result = (uint32_t)(next / 65536) % 2048;
  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (uint32_t)(next / 65536) % 1024;
  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (uint32_t)(next / 65536) % 1024;
  *seed = next;
  return result;
}

INLINE static void inl_drand48_iterate(uint16_t xsubi[3]) {
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

INLINE static double inl_erand48(uint16_t xsubi[3]) {
  union ieee754_double temp;

  /* Compute next state.  */
  inl_drand48_iterate(xsubi);

  /* Construct a positive double with the 48 random bits distributed over
     its fractional part so the resulting FP number is [0.0,1.0).  */
  temp.ieee.negative = 0;
  temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
  temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
  temp.ieee.mantissa1 = (((uint32_t)xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);

  /* Please note the lower 4 bits of mantissa1 are always 0.  */
  return temp.d - 1.0;
}

#else

/* In the case of OSX, we default to the platform's
   default implementation. */

INLINE static int inl_rand_r(uint32_t *seed) { return rand_r(seed); }

INLINE static double inl_erand48(uint16_t xsubi[3]) { return erand48(xsubi); }

#endif

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
INLINE static double random_unit_interval(int64_t id,
                                          const integertime_t ti_current,
                                          const enum random_number_type type) {

  /* Start by packing the state into a sequence of 16-bit seeds for rand_r. */
  uint16_t buff[9];
  id += type;
  memcpy(&buff[0], &id, 8);
  memcpy(&buff[4], &ti_current, 8);

  /* The inputs give 16 bytes of state, but we need a multiple of 6 for the
     calls to erand48(), so we add an additional aribrary constant two-byte
     value to get 18 bytes of state. */
  buff[8] = 6178;

  /* Use the random seed to generate a new random number */
  buff[0] = buff[0] ^ (uint16_t)SWIFT_RANDOM_SEED_XOR;

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

INLINE static double random_unit_interval_two_IDs(
    int64_t id_star, int64_t id_gas, const integertime_t ti_current,
    const enum random_number_type type) {

  /* We need to combine the gas and star IDs such that we do not get correlation
   * for same id_star + id_gas pairs, because of this we combine everything
   * nonlinearly */
  int64_t input_id = (id_star * id_gas + id_star * ti_current +
                      id_gas * ti_current * ti_current) %
                     INT64_MAX;

  return random_unit_interval(input_id, ti_current, type);
}

INLINE static double random_unit_interval_part_ID_and_ray_idx(
    int64_t id_star, const int ray_idx, const integertime_t ti_current,
    const enum random_number_type type) {

  /* For better mixing, we apply a non-linear transformation y=x^3 */
  const long long ray_idx_3 = ray_idx * ray_idx * ray_idx;
  return random_unit_interval_two_IDs(id_star, ray_idx_3, ti_current, type);
}

#endif /* SWIFT_RANDOM_H */
