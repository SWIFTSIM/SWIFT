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

// #include "cosmology.h"

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
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is initialized
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param star_age the star particle's age
 * @param dt current time step size
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(
    struct spart* restrict sp, 
    const struct cosmology *cosmo,
    int with_cosmology,
    const integertime_t ti_current,
    double time,
    double time_base
    ) {

  /* get star's age and time step for stellar emission rates */
  const integertime_t ti_begin = get_integer_time_begin(ti_current - 1, sp->time_bin);
  const integertime_t ti_step = get_integer_timestep(sp->time_bin);

  /* Get particle time-step */
  double dt_star;
  if (with_cosmology) { 
    dt_star = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
  } else {
    dt_star = get_timestep(sp->time_bin, time_base);
  }

  /* Calculate age of the star at current time */
  double star_age_end_of_step;
  if (with_cosmology) {
    star_age_end_of_step = cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->birth_scale_factor, cosmo->a);
  } else {
    star_age_end_of_step = time - (double)sp->birth_time;
  }


  if (ti_current == 0){
    /* if this is the zeroth step, time is still at zero.
     * Do some bogus stuff for now. */
    star_age_end_of_step += 2*dt_star;
  } 
  if (star_age_end_of_step - dt_star >= 0.){
    sp->rt_data.emission_rate_set += 1;
  } else {
    printf("%lld %lld %lld\n", ti_begin, ti_step, ti_current);
    error("Got negative time when setting emission rates? %10.3g %10.3g %10.3g", star_age_end_of_step, dt_star, star_age_end_of_step - dt_star);
  }




  if (sp->id == 137500){
    printf("--- Set Emission Rate 137500 in cell.c\n");
  }

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
 *
 * @param sp star particle
 * @param reset_emission_rate whether to reset the stellar emission
 *        rate.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp, int reset_emission_rate) {

  /* reset everything */
  sp->rt_data.calls_per_step = 0;
  sp->rt_data.iact_hydro_inject = 0;
  sp->rt_data.calls_self_inject = 0;
  sp->rt_data.calls_pair_inject = 0;
  if (reset_emission_rate) sp->rt_data.emission_rate_set = 0;
}

/**
 * @brief First initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {

  sp->rt_data.calls_tot = 0;
  rt_init_spart(sp, /*reset_emission_rate =*/1);
}

#endif /* SWIFT_RT_DEBUG_H */
