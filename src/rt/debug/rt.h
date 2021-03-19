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

#include "rt_debugging.h"
#include "rt_stellar_emission_rate.h"
#include "rt_thermochemistry.h"

/**
 * @file src/rt/debug/rt.h
 * @brief Main header file for the debug radiative transfer scheme.
 */

/**
 * @brief Initialisation of the RT density loop related particle data.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief Reset of the RT hydro particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually. Also, if debugging checks are active, an
 * extra call to rt_reset_part is made in space_convert_rt_quantities() after
 * the zeroth time step is finished.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  p->rt_data.iact_stars_inject = 0;

  p->rt_data.calls_iact_gradient = 0;
  p->rt_data.calls_iact_transport = 0;
  p->rt_data.injection_check = 0;

  p->rt_data.injection_done = 0;
  p->rt_data.gradients_done = 0;
  p->rt_data.transport_done = 0;
  p->rt_data.thermochem_done = 0;
}

/**
 * @brief First initialisation of the RT hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {

  rt_init_part(p);
  rt_reset_part(p);
  p->rt_data.radiation_absorbed_tot = 0ULL;
}

/**
 * @brief Initialisation of the RT density loop related star particle data.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually. Also, if debugging checks are active, an
 * extra call to rt_reset_spart is made in space_convert_rt_quantities() after
 * the zeroth time step is finished.

 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {

  /* reset everything */

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  sp->rt_data.iact_hydro_inject = 0;

  sp->rt_data.emission_rate_set = 0;
  sp->rt_data.injection_check = 0;
}

/**
 * @brief First initialisation of the RT star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {

  rt_init_spart(sp);
  rt_reset_spart(sp);
  sp->rt_data.radiation_emitted_tot = 0ULL;
}

/**
 * @brief Split the RT data of a particle into n pieces
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void rt_split_part(struct part* p,
                                                                double n) {}

/**
 * @brief Exception handle a hydro part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_part_has_no_neighbours(
    struct part* p){};

/**
 * @brief Exception handle a star part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_spart_has_no_neighbours(
    struct spart* sp){};

/**
 * @brief Update the photon number of a particle, i.e. compute
 *  E^{n+1} = E^n + dt * dE_* / dt. This function finalises
 *  the injection step.
 *
 * @param p particle to work on
 * @param props struct #rt_props that contains global RT properties
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p,
                                   struct rt_props* props) {

  if (props->do_all_parts_have_stars_checks && p->rt_data.injection_check != 1)
    error("called ghost1 when injection check count is %d; ID=%lld",
          p->rt_data.injection_check, p->id);
  p->rt_data.injection_done += 1;
}

/**
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is being reset
 *        (during start-up and during stars ghost if spart is active)
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

  /* Skip initial fake time-step */
  if (dt == 0.0l) return;

  if (time == 0.l) {
    /* if function is called before the first actual step, time is still
     * at zero unless specified otherwise in parameter file.*/
    star_age = dt;
  }

  /* now get the emission rates */
  double star_age_begin_of_step = star_age - dt;
  star_age_begin_of_step = max(0.l, star_age_begin_of_step);
  rt_set_stellar_emission_rate(sp, star_age_begin_of_step, star_age);
}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_finalise_gradient(
    struct part* restrict p) {

  if (p->rt_data.injection_done != 1)
    error(
        "Called finalise gradient on particle "
        "where injection count = %d",
        p->rt_data.injection_done);

  if (p->rt_data.calls_iact_gradient == 0)
    error(
        "Called finalise gradient on particle "
        "with iact gradient count = %d",
        p->rt_data.calls_iact_gradient);

  p->rt_data.gradients_done += 1;
}

/**
 * @brief finishes up the transport step
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p) {

  if (p->rt_data.injection_done != 1)
    error(
        "Trying to do finalise_transport when "
        "injection count is %d",
        p->rt_data.injection_done);

  if (p->rt_data.gradients_done != 1)
    error(
        "Trying to do finalise_transport when "
        "rt_finalise_gradient count is %d",
        p->rt_data.gradients_done);

  if (p->rt_data.calls_iact_gradient == 0)
    error(
        "Called finalise transport on particle "
        "with iact gradient count = %d",
        p->rt_data.calls_iact_gradient);

  if (p->rt_data.calls_iact_transport == 0)
    error(
        "Called finalise transport on particle "
        "with iact transport count = %d",
        p->rt_data.calls_iact_transport);

  p->rt_data.transport_done += 1;
}

/**
 * @brief Do the thermochemistry on a particle.
 *
 * This function wraps around rt_do_thermochemistry function.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_tchem(
    struct part* restrict p) {

  rt_do_thermochemistry(p);
}

/**
 * @brief Clean the allocated memory inside the RT properties struct.
 *
 * @param props the #rt_props.
 */
__attribute__((always_inline)) INLINE static void rt_clean(
    struct rt_props* props) {}

#endif /* SWIFT_RT_DEBUG_H */
