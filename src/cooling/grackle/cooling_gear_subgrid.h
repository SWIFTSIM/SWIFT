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
 */
INLINE static void cooling_ionize_part_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time) {
#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
  if (!radiation_is_part_tagged_as_ionized(p, xp)) {
    return;
  }

  const double m_p = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  /* Get the internal energy increase to ionization */
  const double N_H = radiation_get_part_number_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, p, xp);

  const double E_ion =
      2.17872e-11 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  const double Delta_u_ionized = N_H * E_ion / hydro_get_mass(p);

  /* Get internal energy due to collisions */
  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(p);
  const double Z_sun = 0.02;
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

  /* Here we need to treat the cases Z << Z_sun otherwise we have T < 0
   */
  const double ten_to_four_K =
      1e4 * units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  double T_collisional;
  if (Z >= Z_sun * 1e-3) {
    const double tmp = 0.86 / (1 + 0.22 * log(Z / Z_sun));
    T_collisional = ten_to_four_K * min(6.62, tmp);
  } else {
    T_collisional = 6.62 * ten_to_four_K;  // High-temp asymptote
  }

  /* Convert T to internal energy */
  const double u_collisional =
      cooling_internal_energy_from_T(T_collisional, mu, k_B, m_p);

  /* The internal engergy is the min of the energy required to fully
ionize
     and the equilibrium temperature in HII regions */
  const float u_new = min(Delta_u_ionized, u_collisional);

  /* message("u = %e, Delta_u_ionized = %e, u_coll = %e, T_col = %e", */
  /*         u, Delta_u_ionized, u_collisional, T_collisional); */

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

  /* Should we provide the rates? */
#endif

  /* Keep the particle flagged (and re-floored above, every step) until
     the ionizing star's next HII rebuild -- reset only once that window
     has elapsed. */
  if (time >= radiation_get_part_ionized_end_time(p, xp)) {
    radiation_reset_part_ionized_tag(p, xp);
  }
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
 */
INLINE static void cooling_update_part_subgrid(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp, double dt, double dt_therm, double time) {

  /* Apply ionization */
  cooling_ionize_part_subgrid(phys_const, us, cosmo, hydro_props,
                              pressure_floor, cooling, p, xp, dt, dt_therm,
                              time);

  /* TODO (future plan): Apply space-time varying UV background */
}

#endif /* SWIFT_GEAR_COOLING_SUBGRID_H */
