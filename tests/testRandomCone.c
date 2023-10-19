/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2022 Filip Husko (filip.husko@durham.ac.uk).
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
#include <config.h>

/* System includes. */
#include <fenv.h>

/* Local headers. */
#include "swift.h"

/* Number of different maximal opening angles to test between 0 and pi/2. */
const int N_cos = 20;

/* Number of random cones to draw for each opening angle. */
const int N_cos_cone = 30;

/* Cubical grid size when checking the cone function along axes of the Cartesian
 * grid. Should be at least 2. */
const int N_cube = 5;

/**
 * @brief Test to check whether the function that generates random directions
 * within a cone actually generates vectors only within that cone. It also tests
 * whether they are uniformly distributed in solid angle.
 *
 * @param id_bh The ID of a black hole particle around whose spin vector a given
 * cone is drawn.
 * @param ti_current Current time of the simulation.
 * @param type Random number type used.
 * @param opening_angle The opening angle of the cone (in radians).
 * @param unit_vector The vector that defines where the cone is pointing.
 * @param N_test How many random directions to draw within the cone.
 * @param N_bins How many bins to distribute these directions into when testing
 * the uniformity of the distribution.
 * @param tolerance The tolerance of each bin relative to the expected value.
 */
float test_cone(int64_t id_bh, const integertime_t ti_current,
                const enum random_number_type type, double opening_angle,
                float unit_vector[3]) {

  /* Compute cosine that corresponds to the maximum opening angle */
  const double cos_theta_max = cos(opening_angle);

  /* Initialize an array that will hold a random vector every step */
  float rand_vector[3];

  /* Generate a random unit vector within a cone around unit_vector  */
  random_direction_in_cone(id_bh, ti_current, type, opening_angle, unit_vector,
                           rand_vector);

  /* Check that this vector is actually within the cone we want  */
  const double cos_rand_unit = rand_vector[0] * unit_vector[0] +
                               rand_vector[1] * unit_vector[1] +
                               rand_vector[2] * unit_vector[2];
  if (cos_rand_unit < 0.99999 * cos_theta_max) {
    printf("Cos_opening_angle is: %f, Random cos is: %f\n", cos_theta_max,
           cos_rand_unit);
    error("Generated random unit vector is outside cone.");
  }

  return cos_rand_unit;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Log the swift random seed */
  message("SWIFT random seed = %d", SWIFT_RANDOM_SEED_XOR);

  /* Test the random-vector-in-cone function, for different values of opening
   * angle from 0 to pi/2 (in radians). For each of these opening angles we draw
   * some number of cones, and test whether each of those has a uniform
   * distribution of randomly drawn vectors within it.*/
  for (int i = 1; i < N_cos; ++i) {

    /* Opening angle to use */
    const double opening_angle = 0.5 * M_PI * (double)i / (double)N_cos;

    /* Do the test for N_cos_cone cones with this opening angle */
    for (int l = 0; l < N_cos_cone; ++l) {

      /* Generate an id for the bh and a time. We do this for every opening
       * angle and every cone. */
      const long long id_bh = rand() * (1LL << 31) + rand();
      const integertime_t ti_current = rand() * (1LL << 31) + rand();

      /* Generate a random unit vector that defines a cone, along with the
       * opening angle. */
      float unit_vector[3];
      const double cos_unit =
          random_unit_interval(id_bh, ti_current, random_number_BH_kick);
      const double sin_unit = sqrtf(max(0., (1. - cos_unit) * (1. + cos_unit)));
      const double phi_unit =
          (2. * M_PI) * random_unit_interval(id_bh * id_bh, ti_current,
                                             random_number_BH_kick);
      unit_vector[0] = sin_unit * cos(phi_unit);
      unit_vector[1] = sin_unit * sin(phi_unit);
      unit_vector[2] = cos_unit;

      /* Do the test. */
      test_cone(id_bh, ti_current, random_number_BH_kick, opening_angle,
                unit_vector);
    }
  }

  /* Repeat the same test but with a larger number of random directions and
   * bins, but for just one opening angle and one randomly generated cone */
  const double opening_angle_0 = 0.2;

  /* Compute cosine that corresponds to the maximum opening angle */
  const double cos_theta_max = cos(opening_angle_0);

  /* Generate a random unit vector that defines a cone, along with the
   * opening angle. */
  const long long id_bh_0 = rand() * (1LL << 31) + rand();
  const integertime_t ti_current_0 = rand() * (1LL << 31) + rand();

  float unit_vector_0[3];
  const double cos_unit =
      random_unit_interval(id_bh_0, ti_current_0, random_number_BH_kick);
  const double sin_unit = sqrtf(max(0., (1. - cos_unit) * (1. + cos_unit)));
  const double phi_unit =
      (2. * M_PI) * random_unit_interval(id_bh_0 * id_bh_0, ti_current_0,
                                         random_number_BH_kick);
  unit_vector_0[0] = sin_unit * cos(phi_unit);
  unit_vector_0[1] = sin_unit * sin(phi_unit);
  unit_vector_0[2] = cos_unit;

  /* Some parameters to test the uniformity of drawn vectors */
  int N_test = 10000000;
  int N_bins = 100;
  float tolerance = 0.05;

  /* Initialize an array that will hold the binned number of drawn cosines,
     i.e. this is the probability density function that we wish to test. */
  double binned_cosines[N_bins];
  for (int j = 0; j < N_bins; ++j) {
    binned_cosines[j] = 0.;
  }

  /* Draw N_test vectors and bin them to test uniformity */
  for (int k = 0; k < N_test; ++k) {

    const long long id_bh = rand() * (1LL << 31) + rand();
    const integertime_t ti_current = rand() * (1LL << 31) + rand();

    /* Do the test, with a newly generated BH id and time */
    const float cos_rand_unit =
        test_cone(id_bh, ti_current, random_number_BH_kick, opening_angle_0,
                  unit_vector_0);

    /* Add the unit vector to the probability density function array. The solid
     * angle subtended by some angle theta grows as (1-cos(theta)). Furthermore,
     * we are limited to the spherical cap defined by the angles [0, theta_max].
     * Therefore the variable which we expect to be uniformly distributed is (1
     * - cos(theta)) / (1 - cos(theta_max)). */
    double uniform_variable = (1. - cos_rand_unit) / (1 - cos_theta_max);
    for (int j = 0; j < N_bins; ++j) {
      if ((uniform_variable > (double)j / (double)N_bins) &&
          (uniform_variable < (double)(j + 1) / (double)N_bins)) {
        binned_cosines[j] = binned_cosines[j] + 1. / (double)N_test;
      }
    }
  }

  /* Check whether the binned quantity is really uniformly distributed. If it
   * is, the density (value) of each bin should be 1/N_bin. */
  for (int j = 0; j < N_bins; ++j) {
    if ((binned_cosines[j] < (1. - tolerance) / (double)N_bins) ||
        (binned_cosines[j] > (1. + tolerance) / (double)N_bins)) {
      error(
          "Generated distribution of random unit vectors within a cone exceeds "
          "the limit imposed by the tolerance.");
    }
  }

  /* We now repeat the same process, but we do not generate random unit vectors
   * to define the cones. Instead, we sample unit vectors along the grid
   * [-N_cube, -N_cube + 1, ..., N_cube -1, N_cube] ^ 3. This can be, e.g. [-2,
   * -1, 0, 1, 2] ^ 3 (N_cube should be at least 2). This makes sure that the
   * function that generates random unit vectors is well-defined if the unit
   * vectors that define the cones point along any of the Cartesian axes, or if
   * any of their components are equal. Here we use a fixed opening angle of
   * 0.1, since we assume that the earlier test passing means that the function
   * correctly does what it should for all opening angles. */
  const double opening_angle = 0.1;
  for (int x = -N_cube; x < N_cube + 1; ++x) {
    for (int y = -N_cube; y < N_cube + 1; ++y) {
      for (int z = -N_cube; z < N_cube + 1; ++z) {

        /* Create our unit vector on this point of the grid */
        float unit_vector[3] = {(float)x, (float)y, (float)z};
        float unit_vector_norm =
            sqrtf((float)(x * x) + (float)(y * y) + (float)(z * z));

        /* Only do the test if the norm is >0, i.e. if we are not at the origin
         * of the coordinate frame. */
        if (unit_vector_norm > 0) {

          /* Generate an id for the bh and a time. We do this for every opening
           * angle and every cone. */
          const long long id_bh = rand() * (1LL << 31) + rand();
          const integertime_t ti_current = rand() * (1LL << 31) + rand();

          /* Normalize the unit vector. */
          unit_vector[0] = unit_vector[0] / unit_vector_norm;
          unit_vector[1] = unit_vector[1] / unit_vector_norm;
          unit_vector[2] = unit_vector[2] / unit_vector_norm;

          /* Do the test. */
          test_cone(id_bh, ti_current, random_number_BH_kick, opening_angle,
                    unit_vector);
        }
      }
    }
  }

  return 0;
}
