/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#include "../config.h"

#include <fenv.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Test to check that the pseodo-random numbers in SWIFT are random
 * enough for our purpose.
 *
 * @param argc Unused
 * @param argv Unused
 * @return 0 if everything is fine, 1 if random numbers are not random enough.
 */
int main(int argc, char* argv[]) {

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

  /* Time-step size */
  const int time_bin = 30;

  const double boundary[6] = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
  int count[6] = {0};
  unsigned long long int total = 0;

  /* Try a few different values for the ID */
  for (int i = 0; i < 10; ++i) {

    const long long id = rand() * (1LL << 31) + rand();
    const integertime_t increment = (1LL << time_bin);

    message("Testing id=%lld time_bin=%d", id, time_bin);

    /* Check that the numbers are uniform over the full-range of useful
     * time-steps */
    for (integertime_t ti_current = 0LL; ti_current < max_nr_timesteps;
         ti_current += increment) {

      total += 1;
      ti_current += increment;

      const double r =
          random_unit_interval(id, ti_current, random_number_star_formation);

      /* Count the number of random numbers below the boundaries */
      if (r < boundary[0]) count[0] += 1;
      if (r < boundary[1]) count[1] += 1;
      if (r < boundary[2]) count[2] += 1;
      if (r < boundary[3]) count[3] += 1;
      if (r < boundary[4]) count[4] += 1;
      if (r < boundary[5]) count[5] += 1;
    }

    /* Print counted number of random numbers below the boundaries */
    message("Categories           | %6.0e %6.0e %6.0e %6.0e %6.0e %6.0e",
            boundary[0], boundary[1], boundary[2], boundary[3], boundary[4],
            boundary[5]);
    message("total = %12lld | %6d %6d %6d %6d %6d %6d", total, count[0],
            count[1], count[2], count[3], count[4], count[5]);
    /* Calculate the expected amount of random numbers in this range */
    double expected_result[6];
    int expected_result_int[6];
    expected_result[0] = total * boundary[0];
    expected_result[1] = total * boundary[1];
    expected_result[2] = total * boundary[2];
    expected_result[3] = total * boundary[3];
    expected_result[4] = total * boundary[4];
    expected_result[5] = total * boundary[5];

    expected_result_int[0] = (int)expected_result[0];
    expected_result_int[1] = (int)expected_result[1];
    expected_result_int[2] = (int)expected_result[2];
    expected_result_int[3] = (int)expected_result[3];
    expected_result_int[4] = (int)expected_result[4];
    expected_result_int[5] = (int)expected_result[5];

    /* Print the expected numbers */
    message("expected             | %6d %6d %6d %6d %6d %6d",
            expected_result_int[0], expected_result_int[1],
            expected_result_int[2], expected_result_int[3],
            expected_result_int[4], expected_result_int[5]);

    int std_expected_result[6];

    /* Calculate the allowed standard error deviation the maximum of:
     * 1. the standard error of the expected number doing sqrt(N_expected)
     * 2. The standard error of the counted number doing sqrt(N_count)
     * 3. 1 to prevent low number statistics to crash for 1 while expected
     *    close to zero.
     *
     * 1 and 2 are for large numbers essentially the same but for small numbers
     * it becomes imporatant (e.g. count=6 expected=.9, allowed 5+.9 so 6
     * fails, but sqrt(6) ~ 2.5 so it should be fine) */
    std_expected_result[0] =
        (int)max3(sqrt(expected_result[0]), 1, sqrt(count[0]));
    std_expected_result[1] =
        (int)max3(sqrt(expected_result[1]), 1, sqrt(count[1]));
    std_expected_result[2] =
        (int)max3(sqrt(expected_result[2]), 1, sqrt(count[2]));
    std_expected_result[3] =
        (int)max3(sqrt(expected_result[3]), 1, sqrt(count[3]));
    std_expected_result[4] =
        (int)max3(sqrt(expected_result[4]), 1, sqrt(count[4]));
    std_expected_result[5] =
        (int)max3(sqrt(expected_result[5]), 1, sqrt(count[5]));

    /* We want 5 sigma (can be changed if necessary) */
    const int numb_sigma = 5;

    /* Print the differences and the 5 sigma differences */
    message("Difference           | %6d %6d %6d %6d %6d %6d",
            abs(expected_result_int[0] - count[0]),
            abs(expected_result_int[1] - count[1]),
            abs(expected_result_int[2] - count[2]),
            abs(expected_result_int[3] - count[3]),
            abs(expected_result_int[4] - count[4]),
            abs(expected_result_int[5] - count[5]));

    message("5 sigma difference   | %6d %6d %6d %6d %6d %6d",
            numb_sigma * std_expected_result[0],
            numb_sigma * std_expected_result[1],
            numb_sigma * std_expected_result[2],
            numb_sigma * std_expected_result[3],
            numb_sigma * std_expected_result[4],
            numb_sigma * std_expected_result[5]);

    /* Fail if it is not within numb_sigma (5) of the expected difference. */
    if (count[0] >
            expected_result_int[0] + numb_sigma * std_expected_result[0] ||
        count[0] <
            expected_result_int[0] - numb_sigma * std_expected_result[0] ||
        count[1] >
            expected_result_int[1] + numb_sigma * std_expected_result[1] ||
        count[1] <
            expected_result_int[1] - numb_sigma * std_expected_result[1] ||
        count[2] >
            expected_result_int[2] + numb_sigma * std_expected_result[2] ||
        count[2] <
            expected_result_int[2] - numb_sigma * std_expected_result[2] ||
        count[3] >
            expected_result_int[3] + numb_sigma * std_expected_result[3] ||
        count[3] <
            expected_result_int[3] - numb_sigma * std_expected_result[3] ||
        count[4] >
            expected_result_int[4] + numb_sigma * std_expected_result[4] ||
        count[4] <
            expected_result_int[4] - numb_sigma * std_expected_result[4] ||
        count[5] >
            expected_result_int[5] + numb_sigma * std_expected_result[5] ||
        count[5] <
            expected_result_int[5] - numb_sigma * std_expected_result[5]) {
      message("Not all criteria satisfied!");
      return 1;
    }
  }

  message("All good!");
  return 0;
}
