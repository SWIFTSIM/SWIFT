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
 *
 * @param p the particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  p->rt_data.debug_iact_stars_inject = 0;
  p->rt_data.debug_iact_stars_inject_prep = 0;

  p->rt_data.debug_injection_check = 0;
  p->rt_data.debug_calls_iact_gradient_interaction = 0;
  p->rt_data.debug_calls_iact_transport_interaction = 0;

  p->rt_data.debug_injection_done = 0;
  p->rt_data.debug_gradients_done = 0;
  p->rt_data.debug_transport_done = 0;
  p->rt_data.debug_thermochem_done = 0;
}

/**
 * @brief First initialisation of the RT hydro particle data.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {

  rt_init_part(p);
  rt_reset_part(p);
  p->rt_data.debug_radiation_absorbed_tot = 0ULL;
  p->rt_data.debug_iact_stars_inject_prep_tot = 0ULL;
}

/**
 * @brief Initialisation of the RT density loop related star particle data.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually. Also, if debugging checks are active, an
 * extra call to rt_reset_spart is made in space_convert_rt_quantities() after
 * the zeroth time step is finished.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {

  /* reset everything */

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  sp->rt_data.debug_iact_hydro_inject = 0;
  sp->rt_data.debug_iact_hydro_inject_prep = 0;

  sp->rt_data.debug_emission_rate_set = 0;
  sp->rt_data.debug_injection_check = 0;
}

/**
 * @brief First initialisation of the RT star particle data.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {

  rt_init_spart(sp);
  rt_reset_spart(sp);
  sp->rt_data.debug_radiation_emitted_tot = 0ULL;
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
    struct part* p) {
  message("WARNING: found particle without neighbours");
};

/**
 * @brief Exception handle a star part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_spart_has_no_neighbours(
    struct spart* sp){};

/**
 * @brief Do checks/conversions on particles on startup.
 *
 * @param p The particle to work on
 * @param rtp The RT properties struct
 */
__attribute__((always_inline)) INLINE static void rt_convert_quantities(
    struct part* p, const struct rt_props* rtp){};

/**
 * @brief Computes the next radiative transfer time step size
 * of a given particle (during timestep tasks)
 *
 * @param p particle to work on
 * @param rt_props the RT properties struct
 * @param cosmo the cosmology
 */
__attribute__((always_inline)) INLINE static float rt_compute_timestep(
    const struct part* restrict p, const struct rt_props* restrict rt_props,
    const struct cosmology* restrict cosmo) {

  return FLT_MAX;
}

/**
 * @brief Compute the time-step length for an RT step of a particle.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the rt integration. (internal units).
 */
__attribute__((always_inline)) INLINE static double rt_part_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology* cosmo) {
  return 0.0;
}

/**
 * @brief This function finalises the injection step.
 *
 * @param p particle to work on
 * @param props struct #rt_props that contains global RT properties
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p,
                                   struct rt_props* props) {

  if (props->debug_do_all_parts_have_stars_checks &&
      p->rt_data.debug_injection_check != 1)
    error("called ghost1 when injection check count is %d; ID=%lld",
          p->rt_data.debug_injection_check, p->id);
  p->rt_data.debug_injection_done += 1;
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
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param internal_units struct holding internal units
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp, double time,
                                 double star_age, double dt,
                                 const struct rt_props* rt_props,
                                 const struct phys_const* phys_const,
                                 const struct unit_system* internal_units) {

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
  rt_set_stellar_emission_rate(sp, star_age_begin_of_step, star_age, rt_props,
                               phys_const, internal_units);
}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_end_gradient(
    struct part* restrict p) {

  if (p->rt_data.debug_injection_done != 1)
    error(
        "Called finalise gradient on particle "
        "where injection count = %d",
        p->rt_data.debug_injection_done);

  if (p->rt_data.debug_calls_iact_gradient_interaction == 0)
    message(
        "WARNING: Called finalise gradient on particle "
        "with iact gradient count from rt_iact = %d",
        p->rt_data.debug_calls_iact_gradient_interaction);

  p->rt_data.debug_gradients_done += 1;
}

/**
 * @brief finishes up the transport step
 *
 * @param p particle to work on
 * @param dt the current time step of the particle
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p, const double dt) {

  if (p->rt_data.debug_injection_done != 1)
    error(
        "Trying to do finalise_transport when "
        "injection count is %d",
        p->rt_data.debug_injection_done);

  if (p->rt_data.debug_gradients_done != 1)
    error(
        "Trying to do finalise_transport when "
        "rt_finalise_gradient count is %d",
        p->rt_data.debug_gradients_done);

  if (p->rt_data.debug_calls_iact_transport_interaction == 0)
    message(
        "WARNING: Called finalise transport on particle "
        "with iact transport count from rt_iact = %d",
        p->rt_data.debug_calls_iact_transport_interaction);

  p->rt_data.debug_transport_done += 1;
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
