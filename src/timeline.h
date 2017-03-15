/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Local headers. */
#include "inline.h"
#include "intrinsics.h"

#include <math.h>

typedef long long integertime_t;
typedef char timebin_t;

/*! The number of time bins */
#define num_time_bins 56

/*! The maximal number of timesteps in a simulation */
#define max_nr_timesteps (1LL << (num_time_bins + 1))

/*! Fictious time-bin to hold inhibited particles */
#define time_bin_inhibited (num_time_bins + 2)

/**
 * @brief Returns the integer time interval corresponding to a time bin
 *
 * @param bin The time bin of interest.
 */
static INLINE integertime_t get_integer_timestep(timebin_t bin) {

  if (bin <= 0) return 0;
  return 1LL << (bin + 1);
}

/**
 * @brief Returns the time bin corresponding to a given time_step size.
 *
 * Assumes that integertime_t maps to an unsigned long long.
 */
static INLINE timebin_t get_time_bin(integertime_t time_step) {

  /* ((int) log_2(time_step)) - 1 */
  return (timebin_t)(62 - intrinsics_clzll(time_step));
}

/**
 * @brief Returns the physical time interval corresponding to a time bin.
 *
 * @param bin The time bin of interest.
 * @param timeBase the minimal time-step size of the simulation.
 */
static INLINE double get_timestep(timebin_t bin, double timeBase) {

  return get_integer_timestep(bin) * timeBase;
}

/**
 * @brief Returns the integer time corresponding to the start of the time-step
 * given by a time-bin.
 *
 * @param ti_current The current time on the integer time line.
 * @param bin The time bin of interest.
 */
static INLINE integertime_t get_integer_time_begin(integertime_t ti_current,
                                                   timebin_t bin) {

  const integertime_t dti = get_integer_timestep(bin);
  if (dti == 0)
    return 0;
  else
    return dti * ((ti_current - 1) / dti);
}

/**
 * @brief Returns the integer time corresponding to the start of the time-step
 * given by a time-bin.
 *
 * @param ti_current The current time on the integer time line.
 * @param bin The time bin of interest.
 */
static INLINE integertime_t get_integer_time_end(integertime_t ti_current,
                                                 timebin_t bin) {

  const integertime_t dti = get_integer_timestep(bin);
  if (dti == 0)
    return 0;
  else
    return dti * ceil((double)ti_current / (double)dti);
}

/**
 * @brief Returns the highest active time bin at a given point on the time line.
 *
 * @param time The current point on the time line.
 */
static INLINE timebin_t get_max_active_bin(integertime_t time) {

  if (time == 0) return num_time_bins;

  timebin_t bin = 1;
  while (!((1LL << (bin + 1)) & time)) ++bin;

  return bin;
}

#endif /* SWIFT_TIMELINE_H */
