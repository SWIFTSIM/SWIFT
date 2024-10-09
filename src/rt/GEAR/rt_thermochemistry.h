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
#ifndef SWIFT_RT_GEAR_THERMOCHEMISTRY_H
#define SWIFT_RT_GEAR_THERMOCHEMISTRY_H

#include "rt_grackle_utils.h"
#include "rt_interaction_cross_sections.h"
#include "rt_interaction_rates.h"
#include "rt_ionization_equilibrium.h"

/**
 * @file src/rt/GEAR/rt_thermochemistry.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * thermochemistry related functions.
 */

/**
 * @brief initialize particle quantities relevant for the thermochemistry.
 *
 * @param p part to work with
 * @param rt_props rt_properties struct
 * @param hydro_props hydro properties struct
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_tchem_first_init_part(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {

  if (rt_props->set_equilibrium_initial_ionization_mass_fractions) {
    float XHI, XHII, XHeI, XHeII, XHeIII;
    rt_ion_equil_get_mass_fractions(&XHI, &XHII, &XHeI, &XHeII, &XHeIII, p,
                                    rt_props, hydro_props, phys_const, us,
                                    cosmo);
    p->rt_data.tchem.mass_fraction_HI = XHI;
    p->rt_data.tchem.mass_fraction_HII = XHII;
    p->rt_data.tchem.mass_fraction_HeI = XHeI;
    p->rt_data.tchem.mass_fraction_HeII = XHeII;
    p->rt_data.tchem.mass_fraction_HeIII = XHeIII;
  } else if (rt_props->set_initial_ionization_mass_fractions) {
    p->rt_data.tchem.mass_fraction_HI = rt_props->mass_fraction_HI_init;
    p->rt_data.tchem.mass_fraction_HII = rt_props->mass_fraction_HII_init;
    p->rt_data.tchem.mass_fraction_HeI = rt_props->mass_fraction_HeI_init;
    p->rt_data.tchem.mass_fraction_HeII = rt_props->mass_fraction_HeII_init;
    p->rt_data.tchem.mass_fraction_HeIII = rt_props->mass_fraction_HeIII_init;
  }

  /* pretend you have nonzero density so the check doesn't reset the mass
   * fractions */
  p->rho = 1.f;
  /* Check that we didn't do something stupid */
  rt_check_unphysical_mass_fractions(p);
  p->rho = 0.f;

  /* Check that the Hydrogen and Helium mass fractions correspond to those
   * provided by the user in the parameter file. This mass fraction is also
   * passed down to grackle internally, so it is error-prone if left
   * unchecked. */
  const float mH =
      p->rt_data.tchem.mass_fraction_HI + p->rt_data.tchem.mass_fraction_HII;
  if (fabsf(mH - rt_props->hydrogen_mass_fraction) > 1e-4)
    error("Got wrong Hydrogen mass fraction: Got =%.6f provided in yml =%.6f",
          mH, rt_props->hydrogen_mass_fraction);
  const float mHe = p->rt_data.tchem.mass_fraction_HeI +
                    p->rt_data.tchem.mass_fraction_HeII +
                    p->rt_data.tchem.mass_fraction_HeIII;
  if (fabsf(mHe - rt_props->helium_mass_fraction) > 1e-4)
    error("Got wrong Helium mass fraction: Got =%.6f provided in yml =%.6f",
          mHe, rt_props->helium_mass_fraction);
}

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 * @param depth recursion depth
 */
INLINE static void rt_do_thermochemistry(
    struct part* restrict p, struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us, const double dt, int depth) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* Nothing to do here? */
  if (rt_props->skip_thermochemistry) return;
  if (dt == 0.) return;

  /* This is where the fun begins */
  /* ---------------------------- */

  /* initialize data so it'll be in scope */
  grackle_field_data particle_grackle_data;

  gr_float density = hydro_get_physical_density(p, cosmo);
  /* In rare cases, unphysical solutions can arise with negative densities
   * which won't be fixed in the hydro part until further down the dependency
   * graph. Also, we can have vacuum, in which case we have nothing to do here.
   * So exit early if that is the case. */
  if (density <= 0.) return;

  const float u_minimal = hydro_props->minimal_internal_energy;

  /* Physical internal energy */
  gr_float internal_energy_phys =
      hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  const float u_old = internal_energy;

  gr_float species_densities[6];
  rt_tchem_get_species_densities(p, density, species_densities);

  float radiation_energy_density[RT_NGROUPS];
  rt_part_get_physical_radiation_energy_density(p, radiation_energy_density,
                                                cosmo);

  gr_float iact_rates[5];
  rt_get_interaction_rates_for_grackle(
      iact_rates, radiation_energy_density, species_densities,
      rt_props->average_photon_energy, rt_props->energy_weighted_cross_sections,
      rt_props->number_weighted_cross_sections, phys_const, us);

  /* Put all the data into a grackle field struct */
  rt_get_grackle_particle_fields(&particle_grackle_data, density,
                                 internal_energy, species_densities,
                                 iact_rates);

  /* solve chemistry */
  if (local_solve_chemistry(
          &rt_props->grackle_chemistry_data, &rt_props->grackle_chemistry_rates,
          &rt_props->grackle_units, &particle_grackle_data, dt) == 0)
    error("Error in solve_chemistry.");

  /* copy updated grackle data to particle */
  /* update particle internal energy. Grackle had access by reference
   * to internal_energy */
  internal_energy_phys = particle_grackle_data.internal_energy[0];

  const float u_new = max(internal_energy_phys, u_minimal);

  /* Re-do thermochemistry? */
  if ((rt_props->max_tchem_recursion > depth) &&
      (fabsf(u_old - u_new) > 0.1 * u_old)) {
    /* Note that grackle already has internal "10% rules". But sometimes, they
     * may not suffice. */
    rt_clean_grackle_fields(&particle_grackle_data);
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us,
                          0.5 * dt, depth + 1);
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us,
                          0.5 * dt, depth + 1);
    return;
  }

  /* If we're good, update the particle data from grackle results */
  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);

  /* Update mass fractions */
  const gr_float one_over_rho = 1. / density;
  p->rt_data.tchem.mass_fraction_HI =
      particle_grackle_data.HI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HII =
      particle_grackle_data.HII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeI =
      particle_grackle_data.HeI_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeII =
      particle_grackle_data.HeII_density[0] * one_over_rho;
  p->rt_data.tchem.mass_fraction_HeIII =
      particle_grackle_data.HeIII_density[0] * one_over_rho;

  rt_check_unphysical_mass_fractions(p);

  /* Update radiation fields */
  /* First get absorption rates at the start and the end of the step */
  double absorption_rates[RT_NGROUPS];
  rt_get_absorption_rates(
      absorption_rates, species_densities, rt_props->average_photon_energy,
      rt_props->number_weighted_cross_sections, phys_const, us);

  gr_float species_densities_new[6];
  species_densities_new[0] = particle_grackle_data.HI_density[0];
  species_densities_new[1] = particle_grackle_data.HII_density[0];
  species_densities_new[2] = particle_grackle_data.HeI_density[0];
  species_densities_new[3] = particle_grackle_data.HeII_density[0];
  species_densities_new[4] = particle_grackle_data.HeIII_density[0];
  species_densities_new[5] = particle_grackle_data.e_density[0];
  double absorption_rates_new[RT_NGROUPS];
  rt_get_absorption_rates(absorption_rates_new, species_densities_new,
                          rt_props->average_photon_energy,
                          rt_props->number_weighted_cross_sections, phys_const,
                          us);

  /* Now remove absorbed radiation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float E_old = p->rt_data.radiation[g].energy_density;
    double f = dt * 0.5 * (absorption_rates[g] + absorption_rates_new[g]);
    f = min(1., f);
    f = max(0., f);
    p->rt_data.radiation[g].energy_density *= (1. - f);
    for (int i = 0; i < 3; i++) {
      p->rt_data.radiation[g].flux[i] *= (1. - f);
    }

    rt_check_unphysical_state(&p->rt_data.radiation[g].energy_density,
                              p->rt_data.radiation[g].flux, E_old,
                              /*callloc=*/2);
  }
  /* Clean up after yourself. */
  rt_clean_grackle_fields(&particle_grackle_data);
}

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 */
__attribute__((always_inline)) INLINE static float rt_tchem_get_tchem_time(
    const struct part* restrict p, const struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* initialize data so it'll be in scope */
  grackle_field_data particle_grackle_data;

  gr_float density = hydro_get_physical_density(p, cosmo);
  const float u_minimal = hydro_props->minimal_internal_energy;

  /* Physical internal energy */
  gr_float internal_energy_phys =
      hydro_get_physical_internal_energy(p, xp, cosmo);
  gr_float internal_energy = max(internal_energy_phys, u_minimal);

  gr_float species_densities[6];
  rt_tchem_get_species_densities(p, density, species_densities);

  float radiation_energy_density[RT_NGROUPS];
  rt_part_get_physical_radiation_energy_density(p, radiation_energy_density,
                                                cosmo);

  gr_float iact_rates[5];
  rt_get_interaction_rates_for_grackle(
      iact_rates, radiation_energy_density, species_densities,
      rt_props->average_photon_energy, rt_props->energy_weighted_cross_sections,
      rt_props->number_weighted_cross_sections, phys_const, us);

  rt_get_grackle_particle_fields(&particle_grackle_data, density,
                                 internal_energy, species_densities,
                                 iact_rates);

  /* Compute 'cooling' time */
  gr_float tchem_time;
  if (local_calculate_cooling_time(
          &rt_props->grackle_chemistry_data, &rt_props->grackle_chemistry_rates,
          &rt_props->grackle_units, &particle_grackle_data, &tchem_time) == 0)
    error("Error in calculate_cooling_time.");
  /* Clean up after yourself. */
  rt_clean_grackle_fields(&particle_grackle_data);

  return (float)tchem_time;
}

#endif /* SWIFT_RT_GEAR_THERMOCHEMISTRY_H */
