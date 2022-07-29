/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_TIMELINE_H
#define SWIFT_TIMELINE_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "inline.h"
#include "intrinsics.h"

#include <math.h>
#include <stdint.h>

typedef long long integertime_t;
typedef int8_t timebin_t;

/*! The number of time bins */
#define num_time_bins 56

/*! The maximal number of timesteps in a simulation */
#define max_nr_timesteps (1LL << (num_time_bins + 1))

/*! Fictious time-bin to hold inhibited particles */
#define time_bin_inhibited (num_time_bins + 2)

/*! Fictious time-bin to hold particles not yet created */
#define time_bin_not_created (num_time_bins + 3)

/*! Fictitious time-bin for particles not awaken */
#define time_bin_not_awake (-num_time_bins)

/* Maximal difference in time-bins between neighbouring particles */
#define time_bin_neighbour_max_delta_bin 2

/**
 * @brief Returns the integer time interval corresponding to a time bin
 *
 * @param bin The time bin of interest.
 */
__attribute__((const)) static INLINE integertime_t
get_integer_timestep(const timebin_t bin) {

  if (bin <= 0) return 0;
  return 1LL << (bin + 1);
}

/**
 * @brief Returns the time bin corresponding to a given time_step size.
 *
 * Assumes that integertime_t maps to an unsigned long long.
 * Given our definitions, this is log_2 of the time_step rounded down minus one.
 *
 * We use a fast (but exact for any non-zero value) logarithm in base 2
 * calculation based on the bit representation of the number:
 * log_2(x) = (number of bits in the type) - (number of leading 0-bits in x) - 1
 *
 * @param time_step An integer time-step length.
 */
__attribute__((const)) static INLINE timebin_t
get_time_bin(const integertime_t time_step) {

  /* ((int) log_2(time_step)) - 1 */
  return (timebin_t)((8 * sizeof(integertime_t) - 2) -
                     intrinsics_clzll((unsigned long long)time_step));
}

/**
 * @brief Returns the physical time interval corresponding to a time bin.
 *
 * @param bin The time bin of interest.
 * @param time_base the minimal time-step size of the simulation.
 */
__attribute__((const)) static INLINE double get_timestep(
    const timebin_t bin, const double time_base) {

  return get_integer_timestep(bin) * time_base;
}

/**
 * @brief Returns the integer time corresponding to the start of the time-step
 * given by a time-bin.
 * If the current time is a possible beginning for the given time-bin, return
 * the current time minus the time-step size.
 *
 * @param ti_current The current time on the integer time line.
 * @param bin The time bin of interest.
 */
__attribute__((const)) static INLINE integertime_t
get_integer_time_begin(const integertime_t ti_current, const timebin_t bin) {

  const integertime_t dti = get_integer_timestep(bin);
  if (dti == 0)
    return 0;
  else
    return dti * ((ti_current - 1) / dti);
}

/**
 * @brief Returns the integer time corresponding to the end of the time-step
 * given by a time-bin.
 * If the current time is a possible end for the given time-bin, return the
 * current time.
 *
 * @param ti_current The current time on the integer time line.
 * @param bin The time bin of interest.
 */
__attribute__((const)) static INLINE integertime_t
get_integer_time_end(const integertime_t ti_current, const timebin_t bin) {

  const integertime_t dti = get_integer_timestep(bin);
  if (dti == 0)
    return 0;
  else {
    const integertime_t mod = ti_current % dti;
    if (mod == 0)
      return ti_current;
    else
      return ti_current - mod + dti;
  }
}

/**
 * @brief Returns the highest active time bin at a given point on the time line.
 *
 * @param time The current point on the time line.
 */
__attribute__((const)) static INLINE timebin_t
get_max_active_bin(const integertime_t time) {

  if (time == 0) return num_time_bins;

  timebin_t bin = 1;
  while (!((1LL << (bin + 1)) & time)) ++bin;

  return bin;
}

/**
 * @brief Returns the lowest active time bin at a given point on the time line.
 *
 * @param ti_current The current point on the time line.
 * @param ti_old The last synchronisation point on the time line.
 */
__attribute__((const)) static INLINE timebin_t
get_min_active_bin(const integertime_t ti_current, const integertime_t ti_old) {

  const timebin_t min_bin = get_max_active_bin(ti_current - ti_old);
  return min_bin;
}

#endif /* SWIFT_TIMELINE_H */
