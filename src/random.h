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
#include <config.h>

/* Standard header */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

/* Local headers */
#include "sincos.h"

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
 */
enum random_number_type {
  random_number_star_formation = 0LL,
  random_number_sink_formation = 5947309451LL,
  random_number_stellar_feedback_1 = 3947008991LL,
  random_number_stellar_feedback_2 = 6977309513LL,
  random_number_stellar_feedback_3 = 9762399103LL,
  random_number_isotropic_SNII_feedback_ray_theta = 3298327511LL,
  random_number_isotropic_SNII_feedback_ray_phi = 6311114273LL,
  random_number_isotropic_SNIa_feedback_ray_theta = 11134675471LL,
  random_number_isotropic_SNIa_feedback_ray_phi = 5165786851LL,
  random_number_isotropic_AGN_feedback_ray_theta = 8899891613LL,
  random_number_isotropic_AGN_feedback_ray_phi = 10594523341LL,
  random_number_stellar_enrichment = 2936881973LL,
  random_number_BH_feedback = 1640531371LL,
  random_number_BH_swallow = 4947009007LL,
  random_number_BH_reposition = 59969537LL,
  random_number_BH_spin = 193877777LL,
  random_number_BH_kick = 303595777LL,
  random_number_snapshot_sampling = 6561001LL,
  random_number_stellar_winds = 5947309451LL,
  random_number_HII_regions = 8134165677LL,
  random_number_enrichment_1 = 7245742351LL,
  random_number_enrichment_2 = 1156895281LL,
  random_number_enrichment_3 = 2189989727LL,
  random_number_mosaic_powerlaw = 406586897LL,
  random_number_mosaic_schechter = 562448657LL,
  random_number_mosaic_poisson = 384160001LL,
  random_number_powerspectrum_split = 126247697LL,
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

/**
 * @brief Returns a pseudo-random number in the range [0, 1[.
 *
 * We generate numbers that are always reproducible for a given pair of particle
 * IDs and simulation time (on the integer time-line). If more than one number
 * per time-step per particle is needed, additional randomness can be obtained
 * by using the type argument.
 *
 * @param id_star The ID of the first particle for which to generate a number.
 * @param id_gas The ID of the second particle for which to generate a number.
 * @param ti_current The time (on the time-line) for which to generate a number.
 * @param type The #random_number_type to generate.
 * @return a random number in the interval [0, 1.[.
 */
INLINE static double random_unit_interval_two_IDs(
    const int64_t id_star, const int64_t id_gas, const integertime_t ti_current,
    const enum random_number_type type) {

  /* We need to combine the gas and star IDs such that we do not get correlation
   * for same id_star + id_gas pairs, because of this we combine everything
   * nonlinearly */
  int64_t input_id = (id_star * id_gas + id_star * ti_current +
                      id_gas * ti_current * ti_current) %
                     INT64_MAX;

  return random_unit_interval(input_id, ti_current, type);
}

/**
 * @brief Returns a pseudo-random number in the range [0, 1[.
 *
 * We generate numbers that are always reproducible for a given particle
 * ID, integer index and simulation time (on the integer time-line). If more
 * than one number per time-step per particle is needed, additional randomness
 * can be obtained by using the type argument.
 *
 * @param id the particle ID.
 * @param index A integer added to the randomness.
 * @param ti_current The current integer time.
 * @param type The type of random number (aka.  #random_number_type seed).
 */
INLINE static double random_unit_interval_part_ID_and_index(
    const int64_t id, const int index, const integertime_t ti_current,
    const enum random_number_type type) {

  /* For better mixing, we apply a non-linear transformation y=1+x^3 */
  const long long index_3 = index * index * index;
  const long long index_3_one = index_3 + 1LL;

  return random_unit_interval_two_IDs(id, index_3_one, ti_current, type);
}

/**
 * @brief Return a random integer following a Poisson distribution.
 *
 * Uses the Knuth-Junhao method to avoid underflow when lambda is large.
 *
 * @param id The ID of the particle for which to generate a number.
 * @param lambda The parameter of the Poisson distribution.
 * @param ti_current The time (on the time-line) for which to generate a number.
 * @param type The #random_number_type to generate.
 */
INLINE static int random_poisson(const int64_t id, const double lambda,
                                 const integertime_t ti_current,
                                 const enum random_number_type type) {

  const double step = 500.;
  const double exp_step = exp(step);

  double lambda_left = lambda;
  int k = 0;
  double p = 1.;

  do {

    k++;
    const double r =
        random_unit_interval_part_ID_and_index(id, k, ti_current, type);
    p *= r;

    while (p < 1. && lambda_left > 0.) {

      if (lambda_left > step) {
        p *= exp_step;
        lambda_left -= step;
      } else {
        p *= exp(lambda_left);
        lambda_left = 0.;
      }
    };
  } while (p > 1.);

  return k - 1;
}
/**
 * @brief Generates a unit vector within a cone.
 *
 * We draw a random unit vector around the unit vector a = (a0, a1, a2),
 * such that it is uniformly distributed in solid angle from 0 opening
 * angle (perfectly aligned with the vector a) to a maximum of
 * opening_angle radians.
 *
 * @param id_bh The ID of the (BH) particle which is doing the jet feedback.
 * @param ti_current The time (on the time-line) for which to generate the
 *                    random kick direction.
 * @param type The #random_number_type to generate.
 * @param opening_angle Opening angle of the cone.
 * @param a Reference direction that defines the cone.
 * @param rand_cone_direction Return value.
 */
INLINE static void random_direction_in_cone(const int64_t id_bh,
                                            const integertime_t ti_current,
                                            const enum random_number_type type,
                                            const float opening_angle,
                                            const float a[3],
                                            float rand_cone_direction[3]) {

  /* We want to draw a random unit vector from a cone around the unit
   * vector a = (a0, a1, a2). We do this in a frame x'y'z', where the z'
   * axis is aligned with the a vector, i.e. the a vector is one of its
   * basis vectors. The choice of the other two axes, x' and y', is
   * arbitrary. We can use any two unit vectors that are orthogonal to the
   * a vector. These two vectors can be obtained using cross products.
   * However, we need to first start from an initial vector that is not
   * perfectly aligned with a. The choice of this vector is also arbitrary;
   * we choose a vector that is perfectly aligned with the smallest
   * component of the spin vector a.  */
  float init_unit[3] = {0.f, 0.f, 1.f};

  /* Find which of the x, y or z is the smallest components of a. We also
   * need to take into account the possibility that two of the components
   * are equal. In this case it doesn't matter which of the two we
   * choose. */
  const float a0_abs = fabsf(a[0]);
  const float a1_abs = fabsf(a[1]);
  const float a2_abs = fabsf(a[2]);
  if (((a0_abs < a1_abs) && (a0_abs < a2_abs)) ||
      ((a0_abs == a1_abs) && (a0_abs < a2_abs))) {
    init_unit[0] = 1.f;
  } else if (((a1_abs < a2_abs) && (a1_abs < a0_abs)) ||
             ((a1_abs == a2_abs) && (a1_abs < a0_abs))) {
    init_unit[1] = 1.f;
  } else if (((a2_abs < a0_abs) && (a2_abs < a1_abs)) ||
             ((a2_abs == a0_abs) && (a2_abs < a1_abs))) {
    init_unit[2] = 1.f;
  }

  /* Using this vector and the a vector, we can find the first basis
   * vector x (alongside a) that will also be orthogonal to a, by using a
   * cross product, such that x = init_unit cross a. This vector needs to be
   * normalized to get an actual unit vector. */
  float basis_vec_x[3];
  basis_vec_x[0] = init_unit[1] * a[2] - init_unit[2] * a[1];
  basis_vec_x[1] = init_unit[2] * a[0] - init_unit[0] * a[2];
  basis_vec_x[2] = init_unit[0] * a[1] - init_unit[1] * a[0];
  const float basis_vec_x_magn =
      sqrtf(basis_vec_x[0] * basis_vec_x[0] + basis_vec_x[1] * basis_vec_x[1] +
            basis_vec_x[2] * basis_vec_x[2]);
  basis_vec_x[0] = basis_vec_x[0] / basis_vec_x_magn;
  basis_vec_x[1] = basis_vec_x[1] / basis_vec_x_magn;
  basis_vec_x[2] = basis_vec_x[2] / basis_vec_x_magn;

  /* The other basis vector, y, follows as x cross a. It is already
   * normalized since x and a are orthogonal. */
  const float basis_vec_y[3] = {basis_vec_x[1] * a[2] - basis_vec_x[2] * a[1],
                                basis_vec_x[2] * a[0] - basis_vec_x[0] * a[2],
                                basis_vec_x[0] * a[1] - basis_vec_x[1] * a[0]};

  /* Draw a random cosine confined to the range [cos(opening_angle), 1] */
  const float rand_cos_theta =
      1.f - (1.f - cosf(opening_angle)) *
                random_unit_interval(id_bh, ti_current, type);

  /* Get the corresponding sine */
  const float rand_sin_theta =
      sqrtf(max(0.f, (1.f - rand_cos_theta) * (1.f + rand_cos_theta)));

  /* Get a random equitorial angle from [0, 180] deg */
  const float rand_phi = ((float)(2. * M_PI)) *
                         random_unit_interval(id_bh * id_bh, ti_current, type);

  float rand_cos_phi, rand_sin_phi;
  sincosf(rand_phi, &rand_sin_phi, &rand_cos_phi);

  /* We now calculate the direction of the final vector by picking a random
   * direction within a cone around a, with basis_vec_x and basis_vec_y
   * playing the role of x and y unit vectors, and a playing the role of z
   * unit vector. In other words, in vector notation rand_cone_direction =
   * sin(theta)cos(phi) * basis_vec_x + sin(theta)sin(phi) * basis_vec_y
   * + cos(theta) * a . */
  rand_cone_direction[0] = rand_sin_theta * rand_cos_phi * basis_vec_x[0] +
                           rand_sin_theta * rand_sin_phi * basis_vec_y[0] +
                           rand_cos_theta * a[0],
  rand_cone_direction[1] = rand_sin_theta * rand_cos_phi * basis_vec_x[1] +
                           rand_sin_theta * rand_sin_phi * basis_vec_y[1] +
                           rand_cos_theta * a[1],
  rand_cone_direction[2] = rand_sin_theta * rand_cos_phi * basis_vec_x[2] +
                           rand_sin_theta * rand_sin_phi * basis_vec_y[2] +
                           rand_cos_theta * a[2];
}

#endif /* SWIFT_RANDOM_H */
