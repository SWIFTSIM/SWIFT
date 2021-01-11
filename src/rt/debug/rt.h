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
 * @brief Initialisation of the RT density loop related particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief Reset of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {

  p->rt_data.calls_per_step = 0;
  p->rt_data.iact_stars_inject = 0;
  p->rt_data.calls_iact_gradient = 0;
  p->rt_data.calls_iact_transport = 0;
  p->rt_data.photon_number_updated = 0;
}

/**
 * @brief First initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {

  p->rt_data.calls_tot = 0;
  rt_init_part(p);
  rt_reset_part(p);
}

/**
 * @brief Initialisation of the RT density loop related particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {

  /* reset everything */
  sp->rt_data.calls_per_step = 0;
  sp->rt_data.iact_hydro_inject = 0;
  sp->rt_data.emission_rate_set = 0;
}

/**
 * @brief First initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {

  sp->rt_data.calls_tot = 0;
  rt_init_spart(sp);
  rt_reset_spart(sp);
}

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
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is initialized
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param time current system time
 * @param star_age age of the star *at the end of the step*
 * @param dt star time step
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp, double time,
                                 double star_age, double dt) {

  /* first reset old values */
  rt_reset_spart(sp);
  sp->rt_data.calls_tot += 1;
  sp->rt_data.calls_per_step += 1;

  if (time == 0.) {
    /* if this is the zeroth step, time is still at zero.
     * Do some bogus stuff for now. */
    /* TODO: check that this covers every possible case */
    star_age += 2 * dt;
  }
  if (star_age - dt >= 0.) {
    sp->rt_data.emission_rate_set += 1;
  } else {
    error(
        "Got negative time when setting emission rates?"
        " %10.3g %10.3g",
        star_age, dt);
  }
}

#endif /* SWIFT_RT_DEBUG_H */
