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
#ifndef SWIFT_RT_STELLAR_EMISSION_RATE_DEBUG_H
#define SWIFT_RT_STELLAR_EMISSION_RATE_DEBUG_H

/**
 * @file src/rt/debug/rt_stellar_emission_rate.h
 * @brief Main header file for the debug radiative transfer scheme
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

  sp->rt_data.debug_emission_rate_set += 1;
}

#endif /* SWIFT_RT_STELLAR_EMISSION_RATE_DEBUG_H */
