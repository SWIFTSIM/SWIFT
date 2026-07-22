/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_GEAR_COOLING_SUBGRID_H
#define SWIFT_GEAR_COOLING_SUBGRID_H

#include "../../feedback/GEAR/radiation.h"
#include "chemistry.h"
#include "cooling.h"
#include "cooling_utils.h"
#include "hydro.h"
#include "minmax.h"

/**
 * @file src/cooling/grackle/cooling_gear_subgrid.h
 * @brief Subgrid model for GEAR cooling, independent from grackle.
 */

/**
 * @brief Ionize gas particles due to stellar feedback.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current simulation time.
 * @return 1 if the particle was held at the subgrid-ionized floor this
 *         call, 0 otherwise. The caller uses this to skip Grackle's own
 *         chemistry/cooling solve for this step -- otherwise Grackle's
 *         ODE integrator, unaware of the external forcing, pulls the
 *         energy back down within the same step it was just floored.
 */
INLINE static int cooling_ionize_part_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time) {
#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
  if (!radiation_is_part_tagged_as_ionized(p, xp)) {
    return 0;
  }

#if COOLING_GRACKLE_MODE > 0
  if (cooling->HII_couple_ionization_rate) {
    /* Rate-coupled path: cooling_copy_to_grackle() injects this particle's
       Gamma_HI/RT_heating_rate and Grackle's own solver computes the
       resulting state, so the energy is not forced here. The tag is not
       reset here either: that happens later, in cooling_new_energy after
       cooling_copy_to_grackle has consumed it, since resetting first would
       erase it before Grackle ever sees it. Return value means "did we
       force the energy," not "is this particle ionized" -- cooling_new_energy
       still runs the solve. */
    return 0;
  }
#endif

  /* Specific internal energy this particle is held at while ionized
     (shared with radiation_get_part_rate_to_fully_ionize, which uses the
     same value to evaluate the temperature-dependent case-B
     recombination coefficient -- keeping the two consistent). */
  const float u_new = radiation_get_part_ionized_internal_energy(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

  /* Now update the gas internal energy state */
  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

  /* Deactivate the cooling as we are heating the particle */
  hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);

#if COOLING_GRACKLE_MODE > 0
  /* If we have the non-equilibrium cooling, we can also set these values
and
     let grackle solve the network. */
  xp->cooling_data.HI_frac = 0.0;
  xp->cooling_data.HII_frac = 1.0;
#endif

  /* Keep the particle flagged (and re-floored above, every step) until
     the ionizing star's next HII rebuild -- reset only once that window
     has elapsed. */
  if (time >= radiation_get_part_ionized_end_time(p, xp)) {
    radiation_reset_part_ionized_tag(p, xp);
  }
  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Update the gas properties with subgrid physics.
 *
 * GEAR ionize gas particles due to stellar feedback.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the particle extra data
 * @param dt The time-step of this particle.
 * @param dt_therm The time-step operator used for thermal quantities.
 * @param time The current simulation time.
 * @return 1 if the particle was held at the subgrid-ionized floor this
 *         call, 0 otherwise. See #cooling_ionize_part_subgrid.
 */
INLINE static int cooling_update_part_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time) {

  /* Apply ionization */
  const int ionized_this_step = cooling_ionize_part_subgrid(
      phys_const, us, cosmo, hydro_props, pressure_floor, cooling, p, xp, dt,
      dt_therm, time);

  /* TODO (future plan): Apply space-time varying UV background */

  return ionized_this_step;
}

/**
 * @brief Debug-only: hold every non-ionized particle fixed at
 * IONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K, bypassing
 * Grackle's solve entirely. Combined with
 * IONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K, this reproduces
 * STARBENCH's idealized two-fixed-state assumption, isolating the
 * ionization scheme's geometry/timing from real thermal physics. Recomputed
 * every step, not cached, mirroring #cooling_ionize_part_subgrid's own
 * re-flooring.
 *
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param cosmo The #cosmology.
 * @param hydro_props The #hydro_props.
 * @param pressure_floor Properties of the pressure floor.
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param u_out (return) The forced internal energy, only set if this
 *        returns 1.
 * @return 1 if the particle's energy was forced (caller should skip
 *         Grackle's own solve and return u_out directly), 0 otherwise
 *         (i.e. always 0 unless compiled with
 *         IONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K).
 */
INLINE static int cooling_debug_fix_neutral_temperature_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, float *u_out) {
#ifdef IONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);
  const double T_o_internal =
      (double)IONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K *
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  const float u_o = cooling_internal_energy_from_T(
      T_o_internal, mu, phys_const->const_boltzmann_k,
      phys_const->const_proton_mass);
  hydro_set_physical_internal_energy(p, xp, cosmo, u_o);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_o);
  hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);
  *u_out = u_o;
  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Compute Grackle's RT_heating_rate/RT_HI_ionization_rate for a
 * particle tagged by GEAR's rate-coupled HII feedback
 * (GEARFeedback:HII_couple_ionization_rate).
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param time_units Internal-time-to-cgs-time conversion factor.
 * @param heating_rate_cgs (return) Heating rate, in raw cgs.
 * @param HI_ionization_rate (return) Photoionization rate coefficient, in
 *        internal 1/time (Grackle's own expected unit for this field).
 * @return 1 if this particle is rate-coupled and both rates were set, 0
 *         otherwise -- the caller should then fall back to its own generic
 *         RT fields.
 */
INLINE static int cooling_get_rate_coupled_RT_fields_subgrid(
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp, double time_units, double *heating_rate_cgs,
    double *HI_ionization_rate) {
#if COOLING_GRACKLE_MODE > 0
  if (!cooling->HII_couple_ionization_rate ||
      !radiation_is_part_tagged_as_ionized(p, xp)) {
    return 0;
  }

  const double Gamma_HI =
      radiation_get_part_photoionization_rate_coefficient(p, xp);
  *HI_ionization_rate = Gamma_HI;

  /* excess_photon_energy_HI is cached in cgs already (see
     radiation_get_individual_star_mean_excess_photon_energy_HI). */
  const double Gamma_HI_cgs = Gamma_HI / time_units;
  const double E_excess_cgs =
      (double)radiation_get_part_excess_photon_energy_HI(p, xp);
  *heating_rate_cgs = Gamma_HI_cgs * E_excess_cgs;
  return 1;
#else
  return 0;
#endif
}

/**
 * @brief Expire a GEAR rate-coupled HII tag once Grackle has consumed it for
 * this step's solve (see #cooling_get_rate_coupled_RT_fields_subgrid,
 * called earlier via cooling_copy_to_grackle -- resetting the tag before
 * that call would erase it before Grackle ever sees the rate).
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param time The current simulation time.
 */
INLINE static void cooling_expire_rate_coupled_tag_subgrid(
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double time) {
#if COOLING_GRACKLE_MODE > 0
  if (cooling->HII_couple_ionization_rate &&
      radiation_is_part_tagged_as_ionized(p, xp) &&
      time >= radiation_get_part_ionized_end_time(p, xp)) {
    radiation_reset_part_ionized_tag(p, xp);
  }
#endif
}

#endif /* SWIFT_GEAR_COOLING_SUBGRID_H */
