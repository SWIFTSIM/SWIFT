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

/**
 * @file src/rt/debug/rt.h
 * @brief Main header file for the debug radiative transfer scheme.
 */

/**
 * @brief Compute the photon emission rates for this stellar particle.
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

  sp->rt_data.debug_emission_rate_set += 1;

  /* rt_set_stellar_emission_rate(sp, star_age_begin_of_step, star_age,
   *                              rt_props, phys_const, internal_units); */
}

/**
 * @brief Initialisation of the RT density loop related particle data.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually.
 *
 * @param p Particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief Reset the RT hydro particle data not related to the hydro density.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually. To reset RT data needed in each RT sub-cycle,
 * use rt_reset_part_each_subcycle().
 *
 * @param p the particle to work on
 * @param cosmo Cosmology.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p, const struct cosmology* cosmo) {

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  p->rt_data.debug_iact_stars_inject = 0;
  p->rt_data.debug_nsubcycles = 0;
  p->rt_data.debug_kicked = 0;
}

/**
 * @brief Reset RT particle data which needs to be reset each sub-cycle.
 *
 * @param p the particle to work on
 * @param cosmo Cosmology.
 * @param dt the current particle RT time step
 */
__attribute__((always_inline)) INLINE static void rt_reset_part_each_subcycle(
    struct part* restrict p, const struct cosmology* cosmo, double dt) {

  rt_debugging_reset_each_subcycle(p);
}

/**
 * @brief First initialisation of the RT hydro particle data.
 *
 * @param p particle to work on
 * @param cosmo #cosmology data structure.
 * @param rt_props RT properties struct
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p, const struct cosmology* cosmo,
    const struct rt_props* restrict rt_props) {

  rt_init_part(p);
  rt_reset_part(p, cosmo);
  rt_reset_part_each_subcycle(p, cosmo, 0.);
  p->rt_data.debug_radiation_absorbed_tot = 0ULL;
}

/**
 * @brief Initialisation of the RT density loop related star particle data.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {

  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  sp->rt_data.debug_iact_hydro_inject_prep = 0;
  sp->rt_data.debug_iact_hydro_inject = 0;
  sp->rt_data.debug_emission_rate_set = 0;
}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {}

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
                                                                double n) {
  error("RT can't run with split particles for now.");
}

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
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void rt_spart_has_no_neighbours(
    struct spart* sp){};

/**
 * @brief Do checks/conversions on particles on startup.
 *
 * @param p The particle to work on
 * @param rtp The RT properties struct
 * @param hydro_props The hydro properties struct
 * @param phys_const physical constants struct
 * @param us unit_system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_convert_quantities(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {}

/**
 * @brief Computes the next radiative transfer time step size
 * of a given particle (during timestep tasks)
 *
 * @param p Particle to work on.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static float rt_compute_timestep(
    const struct part* restrict p, const struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us) {

  return FLT_MAX;
}

/**
 * @brief Computes the next radiative transfer time step size
 * of a given star particle (during timestep tasks).
 *
 * @param sp spart to work on
 * @param rt_props the RT properties struct
 * @param cosmo the cosmology
 */
__attribute__((always_inline)) INLINE static float rt_compute_spart_timestep(
    const struct spart* restrict sp, const struct rt_props* restrict rt_props,
    const struct cosmology* restrict cosmo) {

  return FLT_MAX;
}

/**
 * @brief Compute the time-step length for an RT step of a particle from given
 * integer times ti_beg and ti_end. This time-step length is then used to
 * compute the actual time integration of the transport/force step and the
 * thermochemistry. This is not used to determine the time-step length during
 * the time-step tasks.
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
__attribute__((always_inline)) INLINE static void rt_finalise_injection(
    struct part* restrict p, struct rt_props* props) {

  rt_debug_sequence_check(p, 1, "rt_ghost1/rt_finalise_injection");

  p->rt_data.debug_injection_done += 1;
}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 * @param cosmo #cosmology data structure.
 */
__attribute__((always_inline)) INLINE static void rt_end_gradient(
    struct part* restrict p, const struct cosmology* cosmo) {

  rt_debug_sequence_check(p, 2, __func__);

  if (p->rt_data.debug_calls_iact_gradient_interaction == 0)
    message(
        "WARNING: Called finalise gradient on particle %lld"
        "with iact gradient count from rt_iact = %d",
        p->id, p->rt_data.debug_calls_iact_gradient_interaction);

  p->rt_data.debug_gradients_done += 1;
}

/**
 * @brief finishes up the transport step
 *
 * @param p particle to work on
 * @param dt the current time step of the particle
 * @param cosmo #cosmology data structure.
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p, struct rt_props* rtp, const double dt,
    const struct cosmology* restrict cosmo) {
  rt_debug_sequence_check(p, 3, __func__);

  if (p->rt_data.debug_calls_iact_transport_interaction == 0)
    message(
        "WARNING: Called finalise transport on particle %lld"
        "with iact transport count from rt_iact = %d",
        p->id, p->rt_data.debug_calls_iact_transport_interaction);

  p->rt_data.debug_transport_done += 1;
}

/**
 * @brief Do the thermochemistry on a particle.
 *
 * This function wraps around rt_do_thermochemistry function.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
__attribute__((always_inline)) INLINE static void rt_tchem(
    struct part* restrict p, struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const double dt) {

  rt_debug_sequence_check(p, 4, __func__);
  p->rt_data.debug_thermochem_done += 1;

  /* rt_do_thermochemistry(p); */
}

/**
 * @brief Extra operations done during the kick. This needs to be
 * done before the particle mass is updated in the hydro_kick_extra.
 *
 * @param p Particle to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 */
__attribute__((always_inline)) INLINE static void rt_kick_extra(
    struct part* p, float dt_therm, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {

  /* Don't account for timestep_sync backward kicks */
  if (dt_therm >= 0.f && dt_grav >= 0.f && dt_hydro >= 0.f &&
      dt_kick_corr >= 0.f) {

    rt_debug_sequence_check(p, 0, __func__);
    p->rt_data.debug_kicked += 1;
  }
}

/**
 * @brief Prepare a particle for the !HYDRO! force calculation.
 * E.g. for the meshless schemes, we need to take into account the
 * mass fluxes of the constituent species between particles.
 * NOTE: don't call this during rt_init_part or rt_reset_part,
 * follow the hydro_prepare_force logic.
 *
 * @param p particle to work on
 **/
__attribute__((always_inline)) INLINE static void rt_prepare_force(
    struct part* p) {}

/**
 * @brief Extra operations to be done during the drift
 *
 * @param p Particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void rt_predict_extra(
    struct part* p, struct xpart* xp, float dt_drift) {}

/**
 * @brief Clean the allocated memory inside the RT properties struct.
 *
 * @param props the #rt_props.
 * @param restart did we restart?
 */
__attribute__((always_inline)) INLINE static void rt_clean(
    struct rt_props* props, int restart) {}

#endif /* SWIFT_RT_DEBUG_H */
