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
#ifndef SWIFT_RT_GEAR_H
#define SWIFT_RT_GEAR_H

#include "rt_debugging.h"
#include "rt_flux.h"
#include "rt_gradients.h"
#include "rt_properties.h"
/* #include "rt_slope_limiters_cell.h" [> skipped for now <] */
#include "rt_stellar_emission_model.h"
#include "rt_thermochemistry.h"

#include <float.h>

/**
 * @file src/rt/GEAR/rt.h
 * @brief Main header file for the GEAR M1 Closure radiative transfer scheme.
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

#ifdef SWIFT_RT_DEBUG_CHECKS
  sp->rt_data.debug_emission_rate_set += 1;
#endif

  /* Skip initial fake time-step */
  if (dt == 0.0l) return;

  if (time == 0.l) {
    /* if function is called before the first actual step, time is still
     * at zero unless specified otherwise in parameter file.*/
    /* We're going to need the star age later for more sophistiscated models,
     * but for now the compiler won't let me get away with keeping this here,
     * so keep it as a comment. */
    /* star_age = dt; */
  }

  /* TODO: this is for later, when we use more sophisticated models. */
  /* now get the emission rates */
  /* double star_age_begin_of_step = star_age - dt; */
  /* star_age_begin_of_step = max(0.l, star_age_begin_of_step); */

  double emission_this_step[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) emission_this_step[g] = 0.;

  if (rt_props->stellar_emission_model == rt_stellar_emission_model_const) {
    rt_get_emission_this_step_const(emission_this_step,
                                    rt_props->stellar_const_emission_rates, dt);
  } else if (rt_props->stellar_emission_model ==
             rt_stellar_emission_model_IlievTest) {
    rt_get_emission_this_step_IlievTest(
        emission_this_step, sp->mass, dt, rt_props->photon_number_integral,
        rt_props->average_photon_energy, phys_const, internal_units);
  } else {
    error("Unknown stellar emission rate model %d",
          rt_props->stellar_emission_model);
  }

  for (int g = 0; g < RT_NGROUPS; g++)
    sp->rt_data.emission_this_step[g] = emission_this_step[g];
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
 * @param p particle to work on
 * @param cosmo Cosmology.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p, const struct cosmology* cosmo) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  p->rt_data.debug_iact_stars_inject = 0;
  p->rt_data.debug_nsubcycles = 0;
  p->rt_data.debug_kicked = 0;
#endif
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

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_debugging_reset_each_subcycle(p);
#endif

  rt_gradients_init(p);
  /* the Gizmo-style slope limiting doesn't help for RT as is,
   * so we're skipping it for now. */
  /* rt_slope_limit_cell_init(p); */

  p->rt_data.flux_dt = dt;
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

  /* Don't reset conserved quantities here! ICs will be overwritten */
  rt_init_part(p);
  rt_reset_part(p, cosmo);
  rt_part_reset_mass_fluxes(p);
  rt_reset_part_each_subcycle(p, cosmo, 0.);
  rt_part_reset_fluxes(p);

#ifdef SWIFT_RT_DEBUG_CHECKS
  p->rt_data.debug_radiation_absorbed_tot = 0ULL;
#endif
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

  for (int i = 0; i < 8; i++) {
    sp->rt_data.octant_weights[i] = 0.f;
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* reset this here as well as in the rt_debugging_checks_end_of_step()
   * routine to test task dependencies are done right */
  sp->rt_data.debug_iact_hydro_inject_prep = 0;
  sp->rt_data.debug_iact_hydro_inject = 0;
  sp->rt_data.debug_emission_rate_set = 0;

  for (int g = 0; g < RT_NGROUPS; g++) {
    sp->rt_data.debug_injected_energy[g] = 0.f;
  }
  for (int g = 0; g < RT_NGROUPS; g++) {
    sp->rt_data.emission_this_step[g] = 0.f;
  }
  sp->rt_data.debug_psi_sum = 0.f;
#endif
}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 *
 * @param sp star particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    sp->rt_data.emission_this_step[g] = 0.f;
  }
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
#ifdef SWIFT_RT_DEBUG_CHECKS
  sp->rt_data.debug_radiation_emitted_tot = 0ULL;
  for (int g = 0; g < RT_NGROUPS; g++) {
    sp->rt_data.debug_injected_energy_tot[g] = 0.f;
  }
#endif
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
    struct spart* sp) {

  /* Reset energy to be injected so that global statistics
   * checks still work */
  for (int g = 0; g < RT_NGROUPS; g++) {
    sp->rt_data.emission_this_step[g] = 0.f;
  }
  message("WARNING: found star without neighbours");
};

/**
 * @brief Do checks/conversions on particles on startup.
 *
 * @param p The particle to work on
 * @param rt_props The RT properties struct
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
    const struct cosmology* restrict cosmo) {

  /* If we're reducing the speed of light, then we may encounter
   * photon fluxes which are way too high than the physically
   * allowable limit. This can lead to catastrophic problems for
   * the propagation of photons, as the pressure tensor assumes
   * the upper limit to be respected. So check this and correct
   * it if necessary.
   * We only read in conserved quantities, so only check those. */

  struct rt_part_data* rtd = &p->rt_data;
  const float Vinv = 1.f / p->geometry.volume;

  /* If we read in radiation energy, we read in
   * total energy and store it as energy density.
   * Same for fluxes.
   * Correct that now. */
  for (int g = 0; g < RT_NGROUPS; g++) {
    rtd->radiation[g].energy_density *= Vinv;
    rtd->radiation[g].flux[0] *= Vinv;
    rtd->radiation[g].flux[1] *= Vinv;
    rtd->radiation[g].flux[2] *= Vinv;

    /* Additional check with possible exit for ICs */
    rt_check_unphysical_state_ICs(p, g, &rtd->radiation[g].energy_density,
                                  rtd->radiation[g].flux,
                                  phys_const->const_speed_light_c);
    /* Check for too high fluxes */
    rt_check_unphysical_state(&rtd->radiation[g].energy_density,
                              rtd->radiation[g].flux, /*e_old =*/0.f,
                              /*callloc=*/0);
  }

  /* If we're setting up ionising equilibrium initial conditions,
   * then the particles need to have their densities known first.
   * So we can call the mass fractions initialization now. */
  rt_tchem_first_init_part(p, rt_props, hydro_props, phys_const, us, cosmo);
}

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

  /* just mimic the gizmo particle "size" for now */
  const float psize = cosmo->a * cosmo->a *
                      powf(p->geometry.volume / hydro_dimension_unit_sphere,
                           hydro_dimension_inv);
  float dt = psize * rt_params.reduced_speed_of_light_inverse *
             rt_props->CFL_condition;

  if (rt_props->skip_thermochemistry) return dt;

  float dt_cool = FLT_MAX;
  if (rt_props->f_limit_cooling_time > 0.f)
    /* Note: cooling time may be negative if the gas is being heated */
    dt_cool = rt_props->f_limit_cooling_time *
              rt_tchem_get_tchem_time(p, xp, rt_props, cosmo, hydro_props,
                                      phys_const, us);

  return min(dt, fabsf(dt_cool));
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

  /* For now, the only thing we care about is the upper threshold for stars. */
  return rt_props->stars_max_timestep;
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
  if (with_cosmology) {
    return cosmology_get_delta_time(cosmo, ti_beg, ti_end);
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief This function finalises the injection step.
 *
 * @param p particle to work on
 * @param props struct #rt_props that contains global RT properties
 */
__attribute__((always_inline)) INLINE static void rt_finalise_injection(
    struct part* restrict p, struct rt_props* props) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_debug_sequence_check(p, 1, "rt_ghost1/rt_finalise_injection");

  p->rt_data.debug_injection_done += 1;
#endif

  for (int g = 0; g < RT_NGROUPS; g++) {
    rt_check_unphysical_state(&p->rt_data.radiation[g].energy_density,
                              p->rt_data.radiation[g].flux, /*e_old=*/0.f,
                              /*callloc=*/3);
  }
}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 * @param cosmo #cosmology data structure.
 */
__attribute__((always_inline)) INLINE static void rt_end_gradient(
    struct part* restrict p, const struct cosmology* cosmo) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_debug_sequence_check(p, 2, __func__);

  if (p->rt_data.debug_calls_iact_gradient_interaction == 0)
    message(
        "WARNING: Called finalise gradient on particle %lld"
        "with iact gradient count from rt_iact = %d",
        p->id, p->rt_data.debug_calls_iact_gradient_interaction);

  p->rt_data.debug_gradients_done += 1;
#endif

  rt_finalise_gradient_part(p);
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

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_debug_sequence_check(p, 3, __func__);

  if (p->rt_data.debug_calls_iact_transport_interaction == 0)
    message(
        "WARNING: Called finalise transport on particle %lld"
        "with iact transport count from rt_iact = %d",
        p->id, p->rt_data.debug_calls_iact_transport_interaction);

  p->rt_data.debug_transport_done += 1;
#endif

  struct rt_part_data* restrict rtd = &p->rt_data;
  const float Vinv = 1.f / p->geometry.volume;

  /* Do not redshift if we have a constant spectrum (type == 0) */
  const float redshift_factor =
      (rtp->stellar_spectrum_type == 0) ? 0. : cosmo->H * dt;

  for (int g = 0; g < RT_NGROUPS; g++) {
    const float e_old = rtd->radiation[g].energy_density;

    /* Note: in this scheme, we're updating d/dt (U * V) + sum F * A * dt = 0.
     * So we'll need the division by the volume here. */

    rtd->radiation[g].energy_density += rtd->flux[g].energy * Vinv;
    rtd->radiation[g].energy_density -=
        rtd->radiation[g].energy_density *
        redshift_factor;  // Energy lost due to redshift

    rtd->radiation[g].flux[0] += rtd->flux[g].flux[0] * Vinv;
    rtd->radiation[g].flux[0] -=
        rtd->radiation[g].flux[0] *
        redshift_factor;  // Energy lost due to redshift

    rtd->radiation[g].flux[1] += rtd->flux[g].flux[1] * Vinv;
    rtd->radiation[g].flux[1] -=
        rtd->radiation[g].flux[1] *
        redshift_factor;  // Energy lost due to redshift

    rtd->radiation[g].flux[2] += rtd->flux[g].flux[2] * Vinv;
    rtd->radiation[g].flux[2] -=
        rtd->radiation[g].flux[2] *
        redshift_factor;  // Energy lost due to redshift

    rt_check_unphysical_state(&rtd->radiation[g].energy_density,
                              rtd->radiation[g].flux, e_old, /*callloc=*/4);
  }

  /* Reset the fluxes now that they have been applied. */
  rt_part_reset_fluxes(p);
  /* Mark the particle as inactive for now. */
  rtd->flux_dt = -1.f;
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

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_debug_sequence_check(p, 4, __func__);
  p->rt_data.debug_thermochem_done += 1;
#endif

  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */
  rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us, dt,
                        0);
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

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Don't account for timestep_sync backward kicks */
  if (dt_therm >= 0.f && dt_grav >= 0.f && dt_hydro >= 0.f &&
      dt_kick_corr >= 0.f) {

    rt_debug_sequence_check(p, 0, __func__);
    p->rt_data.debug_kicked += 1;
  }
#endif

  /* Note: We need to mimick here what Gizmo does for the mass fluxes.
   * The relevant time scale is the hydro time step for the mass fluxes,
   * not the RT times. We also need to prevent the kick to apply the mass
   * fluxes twice, so exit if the particle time step < 0 */

  if (p->flux.dt > 0.0f) {
    /* Update the mass fraction changes due to interparticle fluxes */
    const float current_mass = p->conserved.mass;

    if ((current_mass <= 0.f) || (p->rho <= 0.f)) {
      /* Deal with vacuum. Let hydro deal with actuall mass < 0, just do your
       * mass fractions thing. */
      p->rt_data.tchem.mass_fraction_HI = 0.f;
      p->rt_data.tchem.mass_fraction_HII = 0.f;
      p->rt_data.tchem.mass_fraction_HeI = 0.f;
      p->rt_data.tchem.mass_fraction_HeII = 0.f;
      p->rt_data.tchem.mass_fraction_HeIII = 0.f;
      rt_part_reset_mass_fluxes(p);
      return;
    }

    const float current_mass_HI =
        current_mass * p->rt_data.tchem.mass_fraction_HI;
    const float current_mass_HII =
        current_mass * p->rt_data.tchem.mass_fraction_HII;
    const float current_mass_HeI =
        current_mass * p->rt_data.tchem.mass_fraction_HeI;
    const float current_mass_HeII =
        current_mass * p->rt_data.tchem.mass_fraction_HeII;
    const float current_mass_HeIII =
        current_mass * p->rt_data.tchem.mass_fraction_HeIII;

    /* At this point, we're exchanging (time integrated) mass fluxes,
     * which in rare cases can lead to unphysical results, i.e. negative
     * masses. Make sure we prevent unphysical solutions propagating by
     * enforcing a minumum of zero. */
    const float new_mass_HI =
        max(current_mass_HI + p->rt_data.mass_flux.HI, 0.f);
    const float new_mass_HII =
        max(current_mass_HII + p->rt_data.mass_flux.HII, 0.f);
    const float new_mass_HeI =
        max(current_mass_HeI + p->rt_data.mass_flux.HeI, 0.f);
    const float new_mass_HeII =
        max(current_mass_HeII + p->rt_data.mass_flux.HeII, 0.f);
    const float new_mass_HeIII =
        max(current_mass_HeIII + p->rt_data.mass_flux.HeIII, 0.f);

    const float new_mass_tot = new_mass_HI + new_mass_HII + new_mass_HeI +
                               new_mass_HeII + new_mass_HeIII;

    /* During the initial fake time step, if the mass fractions haven't been set
     * up yet, we could encounter divisions by zero here, so skip that. If it
     * isn't the initial time step, the error will be caught down the line by
     * another call to rt_check_unphysical_mass_fractions() (not the one 10
     * lines below this) */
    if (new_mass_tot == 0.f) return;

    const float new_mass_tot_inv = 1.f / new_mass_tot;

    p->rt_data.tchem.mass_fraction_HI = new_mass_HI * new_mass_tot_inv;
    p->rt_data.tchem.mass_fraction_HII = new_mass_HII * new_mass_tot_inv;
    p->rt_data.tchem.mass_fraction_HeI = new_mass_HeI * new_mass_tot_inv;
    p->rt_data.tchem.mass_fraction_HeII = new_mass_HeII * new_mass_tot_inv;
    p->rt_data.tchem.mass_fraction_HeIII = new_mass_HeIII * new_mass_tot_inv;

    /* Reset fluxes after they have been applied, so they can be collected
     * again even when particle is inactive. */
    rt_part_reset_mass_fluxes(p);

    /* Don't update actual particle mass, that'll be done in the
     * hydro_kick_extra calls */
  }

  rt_check_unphysical_mass_fractions(p);
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
    struct part* p, struct xpart* xp, float dt_drift) {

  float dx[3] = {xp->v_full[0] * dt_drift, xp->v_full[1] * dt_drift,
                 xp->v_full[2] * dt_drift};

  for (int g = 0; g < RT_NGROUPS; g++) {
    float Unew[4];
    rt_gradients_predict_drift(p, Unew, g, dx);
    p->rt_data.radiation[g].energy_density = Unew[0];
    p->rt_data.radiation[g].flux[0] = Unew[1];
    p->rt_data.radiation[g].flux[1] = Unew[2];
    p->rt_data.radiation[g].flux[2] = Unew[3];
  }
}

/**
 * @brief Clean the allocated memory inside the RT properties struct.
 *
 * @param props the #rt_props.
 * @param restart did we restart?
 */
__attribute__((always_inline)) INLINE static void rt_clean(
    struct rt_props* props, int restart) {

  /* If we were restarting, free-ing manually will lead to
   * segfaults since we didn't malloc the stuff */
  if (!restart) {
    /* Clean up grackle data. This is a call to a grackle function */
    local_free_chemistry_data(&props->grackle_chemistry_data,
                              &props->grackle_chemistry_rates);

    for (int g = 0; g < RT_NGROUPS; g++) {
      free(props->energy_weighted_cross_sections[g]);
      free(props->number_weighted_cross_sections[g]);
    }
    free(props->energy_weighted_cross_sections);
    free(props->number_weighted_cross_sections);
  }
}

#endif /* SWIFT_RT_GEAR_H */
