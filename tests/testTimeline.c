/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#include "../config.h"
#include "timeline.h"
#include "tools.h"

/* How many times to repeat randomly drawn tests. */
#define NREPEAT 1048576

/**
 * @brief test get_integer_time_end() function.
 * For each particle time bin, pick a random valid time_end for
 * the dt given by the particle bin; then set a random ti_current
 * by subtracting some time interval < dt from the expected time_end
 * and see whether the recovered time_end matches up.
 *
 * @param bin_min lowest bin to start test from
 * @param bin_max highest bin to run test with
 * @param max_nr_timesteps_test maximal number of timesteps for a sim
 */
void test_get_integer_time_end(timebin_t bin_min, timebin_t bin_max,
                               integertime_t max_nr_timesteps_test) {

  integertime_t dti, step, max_step, set_time_end, ti_current, displacement;
  integertime_t time_end_recovered;

  for (timebin_t bin = bin_min; bin < bin_max; bin++) {

    dti = get_integer_timestep(bin);

    /* Check random state in timeline */
    /* ------------------------------ */

    for (int r = 0; r < NREPEAT; r++) {
      /* First pick a place to set this time_end on the timeline. */

      /* we can't have more than this many steps of this size */
      max_step = max_nr_timesteps_test / dti;

      if (max_step == 0)
        error("max step == 0? bin %d max_nr_steps %lld", bin,
              max_nr_timesteps_test);

      /* Set the time_end at any step in between there */
      step = (integertime_t)(random_uniform(1., max_step));
      set_time_end = step * dti;

      /* Do some safety checks */
      if (set_time_end % dti != 0) error("time_end not divisible by dti?");
      if (set_time_end > max_nr_timesteps_test)
        error("Time end > max_nr_timesteps?");
      if (set_time_end < (integertime_t)0) error("Time end < 0?");

      /* Now mimick a "current time" by removing a fraction of dti from
       * the step, and see whether we recover the correct time_end */
      displacement = (integertime_t)(random_uniform(0., 1. - 1e-12) * dti);
      ti_current = set_time_end - displacement;

      /* Another round of safety checks */
      if (ti_current > set_time_end)
        error(
            "current>time_end? current=%lld time_end=%lld dti=%lld "
            "displacement=%lld bin=%d\n",
            ti_current, set_time_end, dti, displacement, bin);
      if (ti_current < 0)
        error(
            "current<0? current=%lld time_end=%lld dti=%lld "
            "displacement=%lld bin=%d\n",
            ti_current, set_time_end, dti, displacement, bin);

      /* Now the actual check. */
      time_end_recovered = get_integer_time_end(ti_current, bin);

      if (time_end_recovered != set_time_end) {
        error(
            "time_end incorrect: expect=%lld got=%lld diff=%lld; current=%lld "
            "displacement=%lld, dti=%lld, bin=%d",
            set_time_end, time_end_recovered, set_time_end - time_end_recovered,
            ti_current, displacement, dti, bin);
      }
    }

    /* Check time_end = 0 */
    /* ------------------ */
    set_time_end = 0;
    ti_current = 0;
    time_end_recovered = get_integer_time_end(ti_current, bin);
    if (time_end_recovered != set_time_end) {
      error(
          "time_end incorrect: expect=%lld got=%lld diff=%lld; current=%lld "
          "bin=%d",
          set_time_end, time_end_recovered, set_time_end - time_end_recovered,
          ti_current, bin);
    }

    /* Check time_end = final_step */
    /* --------------------------- */
    set_time_end = max_nr_timesteps_test;

    if (dti < NREPEAT) {
      /* Check all possible time states before the end */
      for (integertime_t delta = 1; delta < dti; delta++) {
        ti_current = max_nr_timesteps_test - delta;
        time_end_recovered = get_integer_time_end(ti_current, bin);
        if (time_end_recovered != set_time_end) {
          error(
              "time_end incorrect: expect=%lld got=%lld diff=%lld; "
              "current=%lld bin=%d",
              set_time_end, time_end_recovered,
              set_time_end - time_end_recovered, ti_current, bin);
        }
      }
    } else {
      /* Draw random states again */
      for (int r = 0; r < NREPEAT; r++) {
        /* Do some safety checks */
        if (set_time_end % dti != 0) error("time_end not divisible by dti?");

        /* Now mimick a "current time" by removing a fraction of dti from
         * the step, and see whether we recover the correct time_end */
        displacement = (integertime_t)(random_uniform(0., 1. - 1e-12) * dti);
        ti_current = set_time_end - displacement;

        /* Another round of safety checks */
        if (ti_current > set_time_end)
          error(
              "current>time_end? current=%lld time_end=%lld dti=%lld "
              "displacement=%lld bin=%d\n",
              ti_current, set_time_end, dti, displacement, bin);
        if (ti_current < 0)
          error(
              "current<0? current=%lld time_end=%lld dti=%lld "
              "displacement=%lld bin=%d\n",
              ti_current, set_time_end, dti, displacement, bin);

        /* Now the actual check. */
        time_end_recovered = get_integer_time_end(ti_current, bin);

        if (time_end_recovered != set_time_end) {
          error(
              "time_end incorrect: expect=%lld got=%lld diff=%lld; "
              "current=%lld "
              "displacement=%lld, dti=%lld, bin=%d",
              set_time_end, time_end_recovered,
              set_time_end - time_end_recovered, ti_current, displacement, dti,
              bin);
        }
      }
    }
  }

  printf("Passed %s\n", __func__);
}

/**
 * @brief test get_time_bin by converting all possible bins
 * into timesteps and trying to recover the initial bin
 * by calling get_time_bin()
 * @param bin_min lowest bin to start test from
 * @param bin_max highest bin to run test with
 *
 **/
void test_get_time_bin(timebin_t bin_min, timebin_t bin_max) {

  for (timebin_t bin = bin_min; bin < bin_max; bin++) {
    integertime_t dti = get_integer_timestep(bin);
    timebin_t bin_recovered = get_time_bin(dti);
    if (bin != bin_recovered) {
      error("Expected bin=%d, got=%d", bin, bin_recovered);
    }
  }

  printf("Passed %s\n", __func__);
}

/**
 * @brief test get_integer_time_begin() function.
 * For each particle time bin, pick a random valid time_begin for
 * the dt given by the particle bin; then set a random ti_current
 * by adding some time interval < dt to the expected time_begin
 * and see whether the recovered time_begin matches up.
 *
 * @param bin_min lowest bin to start test from
 * @param bin_max highest bin to run test with
 * @param max_nr_timesteps_test maximal number of timesteps for a sim
 */
void test_get_integer_time_begin(timebin_t bin_min, timebin_t bin_max,
                                 integertime_t max_nr_timesteps_test) {

  integertime_t dti, step, max_step, set_time_begin, ti_current, displacement;
  integertime_t time_begin_recovered;

  for (timebin_t bin = bin_min; bin < bin_max; bin++) {

    dti = get_integer_timestep(bin);

    /* Check random state in timeline */
    /* ------------------------------ */

    for (int r = 0; r < NREPEAT; r++) {
      /* First pick a place to set this time_end on the timeline. */

      /* we can't have more than this many steps of this size */
      max_step = max_nr_timesteps_test / dti;
      if (max_step == 0)
        error("max step == 0? bin %d max_nr_steps %lld", bin,
              max_nr_timesteps_test);

      /* Set the time_end at any step in between there */
      step = (integertime_t)(random_uniform(0., max_step));
      set_time_begin = step * dti;

      /* Do some safety checks */
      if (set_time_begin % dti != 0)
        error("set time_begin not divisible by dti?");
      if (set_time_begin >= max_nr_timesteps_test)
        error("Time begin %lld >= max_nr_timesteps %lld?", set_time_begin,
              max_nr_timesteps);
      if (set_time_begin < (integertime_t)0)
        error("Time begin < 0? set_time_begin=%lld", set_time_begin);

      /* Now mimick a "current time" by adding a fraction to dti from
       * the step, and see whether we recover the correct time_begin */
      displacement = (integertime_t)(random_uniform(0., 1.) * dti);
      ti_current = set_time_begin + displacement;

      /* Another round of safety checks */
      if (ti_current < set_time_begin)
        printf(
            "current<time_begin? current=%lld time_end=%lld dti=%lld "
            "displacement=%lld bin=%d\n",
            ti_current, set_time_begin, dti, displacement, bin);

      /* Now the actual check. */
      time_begin_recovered = get_integer_time_begin(ti_current, bin);

      if (ti_current == set_time_begin) {
        /* If the displacement was zero, then the function shall return
         * the beginning of the timestep that ends at ti_current */
        if (time_begin_recovered + dti != set_time_begin)
          error(
              "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
              "current=%lld "
              "displacement=%lld, dti=%lld",
              /*expect=*/set_time_begin - dti,
              /*got=*/time_begin_recovered,
              /*diff=*/set_time_begin - dti - time_begin_recovered, ti_current,
              displacement, dti);
      } else {
        if (time_begin_recovered != set_time_begin)
          error(
              "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
              "current=%lld displacement=%lld, dti=%lld",
              set_time_begin, time_begin_recovered,
              set_time_begin - time_begin_recovered, ti_current, displacement,
              dti);
      }
    }

    /* Check time_begin = 0 */
    /* ------------------ */
    set_time_begin = 0;
    ti_current = 0;
    time_begin_recovered = get_integer_time_begin(ti_current, bin);
    if (time_begin_recovered != set_time_begin)
      error(
          "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
          "current=%lld dti=%lld",
          set_time_begin, time_begin_recovered,
          set_time_begin - time_begin_recovered, ti_current, dti);

    /* Check time_begin = final_step */
    /* ----------------------------- */
    /* This is a tad nonsensial since in this scenario, we're
     * at an integer time > max_nr_timesteps */
    set_time_begin = max_nr_timesteps_test;

    if (dti < NREPEAT) {
      /* Check all possible time states before the end */
      for (integertime_t delta = 1; delta < dti; delta++) {
        ti_current = max_nr_timesteps_test + delta;
        time_begin_recovered = get_integer_time_begin(ti_current, bin);
        if (time_begin_recovered != set_time_begin)
          error(
              "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
              "current=%lld displacement=%lld, dti=%lld",
              set_time_begin, time_begin_recovered,
              set_time_begin - time_begin_recovered, ti_current, displacement,
              dti);
      }
    } else {
      /* Draw random states again */
      for (int r = 0; r < NREPEAT; r++) {

        /* Now mimick a "current time" by removing a fraction of dti from
         * the step, and see whether we recover the correct time_end */
        displacement = (integertime_t)(random_uniform(0., 1.) * dti);
        ti_current = set_time_begin + displacement;

        /* Another round of safety checks */
        if (ti_current < set_time_begin)
          error(
              "current>time_begin? current=%lld time_begin=%lld dti=%lld "
              "displacement=%lld bin=%d\n",
              ti_current, set_time_begin, dti, displacement, bin);
        if (ti_current < 0)
          error(
              "current<0? current=%lld time_begin=%lld dti=%lld "
              "displacement=%lld bin=%d\n",
              ti_current, set_time_begin, dti, displacement, bin);

        /* Now the actual check. */
        time_begin_recovered = get_integer_time_begin(ti_current, bin);

        if (ti_current == set_time_begin) {
          /* If the displacement was zero, then the function shall return
           * the beginning of the timestep that ends at ti_current */
          if (time_begin_recovered + dti != set_time_begin)
            error(
                "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
                "current=%lld "
                "displacement=%lld, dti=%lld",
                /*expect=*/set_time_begin - dti,
                /*got=*/time_begin_recovered,
                /*diff=*/set_time_begin - dti - time_begin_recovered,
                ti_current, displacement, dti);
        } else {
          if (time_begin_recovered != set_time_begin)
            error(
                "time_begin incorrect: expect=%lld got=%lld diff=%lld; "
                "current=%lld displacement=%lld, dti=%lld",
                set_time_begin, time_begin_recovered,
                set_time_begin - time_begin_recovered, ti_current, displacement,
                dti);
        }
      }
    }
  }

  printf("Passed %s\n", __func__);
}

/**
 * @brief test get_max_active_bin by randomly choosing current times
 * and manually checking whether a higher bin could be active
 *
 * @param bin_max highest bin to run test with
 * @param max_nr_timesteps_test maximal number of timesteps for a sim
 **/
void test_get_max_active_bin(timebin_t bin_max,
                             integertime_t max_nr_timesteps_test) {

  integertime_t ti_current, dti;
  timebin_t max_active_bin, testbin;

  /* Test random ti_currents */
  /* ----------------------- */
  for (int r = 0; r < NREPEAT; r++) {
    /* test time 0 later, so use time >= 1 */
    ti_current = (integertime_t)(random_uniform(1., max_nr_timesteps_test));

    max_active_bin = get_max_active_bin(ti_current);
    testbin = 0;
    dti = get_integer_timestep(max_active_bin);
    for (timebin_t bin = max_active_bin; bin < bin_max && dti <= ti_current;
         bin++) {
      if (ti_current % dti == 0) testbin = bin;
      /* update dti here for exit condition */
      dti = get_integer_timestep(bin + 1);
    }
    if (testbin > max_active_bin)
      error("Found higher active bin: time=%lld max_active_bin=%d found=%d",
            ti_current, max_active_bin, testbin);
  }

  /* Test first 2^bin_max integertimes */
  /* --------------------------------- */
  for (timebin_t bin = 1; bin < bin_max; bin++) {
    ti_current = get_integer_timestep(bin);
    max_active_bin = get_max_active_bin(ti_current);

    if (max_active_bin != bin)
      error("Got wrong max_active_bin: Expect=%d got=%d time=%lld", bin,
            max_active_bin, ti_current);
  }

  /* Test final 2^bin_max integertimes */
  /* --------------------------------- */
  for (timebin_t bin = 1; bin < bin_max; bin++) {
    dti = get_integer_timestep(bin);
    ti_current = max_nr_timesteps_test - dti;
    if (ti_current == 0)
      error("Testing bin which only fits in once into entire timeline");
    max_active_bin = get_max_active_bin(ti_current);

    if (max_active_bin != bin) {
      error("Got wrong max_active_bin: Expect=%d got=%d time=%lld", bin,
            max_active_bin, ti_current);
    }
  }

  /* Make sure everything is active at time zero */
  /* ------------------------------------------- */
  max_active_bin = get_max_active_bin(0);
  /* Note: this should be num_time_bins, not bin_max! */
  if (max_active_bin != num_time_bins)
    error("Didn't get max_active_bin = num_time_bins at t=0");

  printf("Passed %s\n", __func__);
}

/**
 * @brief test get_min_active_bin by randomly choosing current times
 * for a pair of small and big bins
 *
 * @param bin_min smallest bin to run test with
 * @param bin_max highest bin to run test with
 * @param max_nr_timesteps_test maximal number of timesteps for a sim
 **/
void test_get_min_active_bin(timebin_t bin_min, timebin_t bin_max,
                             integertime_t max_nr_timesteps_test) {

  integertime_t dti_lo, dti_hi;
  integertime_t step, max_step;
  integertime_t ti_old, ti_current;
  timebin_t min_active_bin;

  /* Test random ti_currents */
  /* ----------------------- */

  for (timebin_t bin_lo = bin_min; bin_lo < bin_max - 1; bin_lo++) {
    dti_lo = get_integer_timestep(bin_lo);

    for (timebin_t bin_hi = bin_lo + 1; bin_hi < bin_max; bin_hi++) {
      dti_hi = get_integer_timestep(bin_hi);
      max_step = (max_nr_timesteps_test / dti_hi);

      if (NREPEAT / bin_max < max_step) {
        /* Do random draws */
        for (int r = 0; r < NREPEAT / bin_max; r++) {

          step = (integertime_t)(random_uniform(0., max_step));
          ti_current = step * dti_hi;
          ti_old = ti_current - dti_lo;

          if (ti_current % dti_hi != 0)
            error("Time current not divisible by dti_hi");
          if (ti_current % dti_lo != 0)
            error("Time current not divisible by dti_lo");
          if (ti_current <= ti_old)
            error("ti_current=%lld <= ti_old=%lld | bins %d %d; ", ti_current,
                  ti_old, bin_lo, bin_hi);

          min_active_bin = get_min_active_bin(ti_current, ti_old);
          if (min_active_bin != bin_lo)
            error("Got wrong min active bin. Expect=%d got=%d", bin_lo,
                  min_active_bin);
        }
      } else {
        /* systematically check every possibility */
        for (step = 0; step <= max_step; step++) {

          ti_current = step * dti_hi;
          ti_old = ti_current - dti_lo;

          if (ti_current % dti_hi != 0)
            error("Time current not divisible by dti_hi");
          if (ti_current % dti_lo != 0)
            error("Time current not divisible by dti_lo");
          if (ti_current <= ti_old)
            error("ti_current=%lld <= ti_old=%lld | bins %d %d; ", ti_current,
                  ti_old, bin_lo, bin_hi);

          min_active_bin = get_min_active_bin(ti_current, ti_old);
          if (min_active_bin != bin_lo)
            error("Got wrong min active bin. Expect=%d got=%d", bin_lo,
                  min_active_bin);
        }
      }
    }
  }

  printf("Passed %s\n", __func__);
}

/**
 * @brief Check the timeline functions.
 */
int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* run test with the limits we use in SWIFT */
  /* ---------------------------------------- */
  test_get_time_bin(1, num_time_bins);
  test_get_integer_time_begin(1, num_time_bins, max_nr_timesteps);
  test_get_integer_time_end(1, num_time_bins, max_nr_timesteps);
  test_get_max_active_bin(num_time_bins, max_nr_timesteps);
  test_get_min_active_bin(1, num_time_bins, max_nr_timesteps);

  /* run test beyond the limits we use in SWIFT */
  /* ------------------------------------------ */
  /* Get maximal bin and number of timesteps based on type size */
  size_t timebin_bits = sizeof(timebin_t) * 8;
  timebin_t max_num_time_bins = 0;
  for (size_t i = 0; i < timebin_bits - 1; i++) {
    max_num_time_bins |= (1 << i);
  }
  size_t timestep_bits = sizeof(integertime_t) * 8;

  /* timestep_bits -1 because of how << works,
   * additional -1 because integertime_t is signed
   * additional -1 so we can do any timestep at least twice */
  max_num_time_bins = min((size_t)max_num_time_bins, timestep_bits - 3);

  if (num_time_bins < max_num_time_bins) {
    /* Use analogous definition as in timeline.h here */
    integertime_t max_nr_timesteps_test = (1LL << (max_num_time_bins + 1));

    test_get_time_bin(num_time_bins, max_num_time_bins);
    test_get_integer_time_begin(num_time_bins, max_num_time_bins,
                                max_nr_timesteps_test);
    test_get_integer_time_end(num_time_bins, max_num_time_bins,
                              max_nr_timesteps_test);
    test_get_max_active_bin(max_num_time_bins, max_nr_timesteps_test);
    test_get_min_active_bin(1, max_num_time_bins, max_nr_timesteps_test);
  }

  return 0;
}
