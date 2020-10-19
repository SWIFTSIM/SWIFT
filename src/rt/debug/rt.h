/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_DEBUG_H
#define SWIFT_RT_DEBUG_H

/**
 * @file src/rt/debug/rt.h
 * @brief Main header file for the debug radiative transfer scheme.
 */

/**
 * @brief Update the photon number of a particle, i.e. compute
 *        E^{n+1} = E^n + dt * dE_* / dt
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p) {

  p->rt_data.photon_number_updated = 1;
  p->rt_data.calls_tot += 1;
  p->rt_data.calls_per_step += 1;
}

/**
 * @brief Initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {

  p->rt_data.calls_per_step = 0;
  p->rt_data.iact_stars_inject = 0;
  p->rt_data.calls_self_inject = 0;
  p->rt_data.calls_pair_inject = 0;
  p->rt_data.photon_number_updated = 0;
}

/**
 * @brief First initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {

  p->rt_data.calls_tot = 0;
  rt_init_part(p);
}

/**
 * @brief Initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {

  sp->rt_data.calls_per_step = 0;
  sp->rt_data.iact_hydro_inject = 0;
  sp->rt_data.calls_self_inject = 0;
  sp->rt_data.calls_pair_inject = 0;
}

/**
 * @brief First initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {

  sp->rt_data.calls_tot = 0;
  rt_init_spart(sp);
}

#endif /* SWIFT_RT_DEBUG_H */
