/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H
#define SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H

/**
 * @file src/rt/GEAR/rt_stellar_emission_rate.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * stellar radiation emission rates related functions.
 */

/**
 * @brief Main function for getting the stellar emission rate and updating it
 * in the spart.
 *
 * @param sp Star particle to work on.
 * @param age_beg Age of the stars at the beginning of the step
 * @param age_end Age of the stars at the end of the step
 * @param rt_props RT properties struct
 * @param phys_const struct holding physical constants
 * @param internal_units units struct containing internal units
 */

__attribute__((always_inline)) INLINE static void rt_set_stellar_emission_rate(
    struct spart* restrict sp, double age_beg, double age_end,
    const struct rt_props* rt_props, const struct phys_const* phys_const,
    const struct unit_system* internal_units) {

  if (rt_props->use_const_emission_rates) {
    /* The read-in constant stellar emisison rates are in units of L_sol.
     * But they have been read in assuming they are in cgs. Convert this
     * only now to proper internal units to avoid float overflows. We only
     * store the energy that is to be distributed from this spart to its
     * neighbours in this step in internal units.*/
    const double solar_luminosity = 3.826e33; /* erg/s */
    const double dt = (age_end - age_beg);
    for (int g = 0; g < RT_NGROUPS; g++) {
      double emission_rate_internal_units =
          rt_props->stellar_const_emission_rates[g] * solar_luminosity;
      sp->rt_data.emission_this_step[g] = emission_rate_internal_units * dt;
    }
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  sp->rt_data.debug_emission_rate_set += 1;
#endif
}

#endif /* SWIFT_RT_STELLAR_EMISSION_RATE_GEAR_H */
