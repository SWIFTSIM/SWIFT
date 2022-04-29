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

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3

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
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_tchem_first_init_part(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {

  if (rt_props->set_equilibrium_initial_ionization_mass_fractions) {
    float XHI, XHII, XHeI, XHeII, XHeIII;
    rt_ion_equil_get_mass_fractions(&XHI, &XHII, &XHeI, &XHeII, &XHeIII, p,
                                    rt_props, phys_const, us, cosmo);
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

  /* Check that we didn't do something stupid */
  rt_check_unphysical_mass_fractions(p);
}

__attribute__((always_inline)) INLINE static void rt_tchem_copy_data_to_grackle(
    grackle_field_data* grackle_field, int grid_dimension[GRACKLE_RANK],
    int grid_start[GRACKLE_RANK], int grid_end[GRACKLE_RANK], gr_float* density,
    gr_float* internal_energy, gr_float species_densities[6],
    gr_float iact_rates[5]) {

  grackle_field->grid_dx = 0.;
  grackle_field->grid_rank = GRACKLE_RANK;
  grackle_field->grid_dimension = grid_dimension;
  grackle_field->grid_start = grid_start;
  grackle_field->grid_end = grid_end;

  /* initialize density */
  grackle_field->density = density;
  grackle_field->internal_energy = internal_energy;
  /* grackle 3.0 doc: "Currently not used" */
  grackle_field->x_velocity = NULL;
  grackle_field->y_velocity = NULL;
  grackle_field->z_velocity = NULL;

  grackle_field->HI_density = &species_densities[0];
  grackle_field->HII_density = &species_densities[1];
  grackle_field->HeI_density = &species_densities[2];
  grackle_field->HeII_density = &species_densities[3];
  grackle_field->HeIII_density = &species_densities[4];
  grackle_field->e_density = &species_densities[5];

  /* general particle data */
  grackle_field->volumetric_heating_rate = NULL;
  grackle_field->specific_heating_rate = NULL;

  grackle_field->RT_heating_rate = &iact_rates[0];
  grackle_field->RT_HI_ionization_rate = &iact_rates[1];
  grackle_field->RT_HeI_ionization_rate = &iact_rates[2];
  grackle_field->RT_HeII_ionization_rate = &iact_rates[3];
  grackle_field->RT_H2_dissociation_rate = &iact_rates[4];

  grackle_field->metal_density = NULL;
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
 */
static void rt_do_thermochemistry(struct part* restrict p,
                                  struct xpart* restrict xp,
                                  struct rt_props* rt_props,
                                  const struct cosmology* restrict cosmo,
                                  const struct hydro_props* hydro_props,
                                  const struct phys_const* restrict phys_const,
                                  const struct unit_system* restrict us,
                                  const double dt) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* Nothing to do here? */
  if (rt_props->skip_thermochemistry) return;
  if (dt == 0.) return;

  /* This is where the fun begins */
  /* ---------------------------- */

  /* initialize data so it'll be in scope */
  grackle_field_data particle_grackle_data;

  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  gr_float density = hydro_get_physical_density(p, cosmo);
  const float u_minimal = hydro_props->minimal_internal_energy;
  gr_float internal_energy =
      max(hydro_get_physical_internal_energy(p, xp, cosmo), u_minimal);

  const float u_old = internal_energy;
  gr_float species_densities[6];
  rt_tchem_get_species_densities(p, density, species_densities);

  gr_float iact_rates[5];
  float iact_rates_by_frequency_bin[RT_NGROUPS];
  rt_tchem_get_interaction_rates(iact_rates, iact_rates_by_frequency_bin, p,
                                 species_densities, rt_props, phys_const, us,
                                 cosmo);

  rt_tchem_copy_data_to_grackle(
      &particle_grackle_data, grid_dimension, grid_start, grid_end, &density,
      &internal_energy, species_densities, iact_rates);

  /* solve chemistry */
  /* Note: grackle_rates is a global variable defined by grackle itself.
   * Using a manually allocd and initialized variable here fails with MPI
   * for some reason. */
  if (local_solve_chemistry(&rt_props->grackle_chemistry_data, &grackle_rates,
                            &rt_props->grackle_units, &particle_grackle_data,
                            dt) == 0)
    error("Error in solve_chemistry.");

  /* update particle internal energy. Grackle had access by reference
   * to internal_energy */
  const float u_new = max(internal_energy, u_minimal);

  /* Redo if we changed by more than 10% ? */
  if (fabsf(u_new - u_old) > 0.1f * u_old) {
    /* Redo because we changed by more than 10% ! */
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us,
                          0.5 * dt);
    rt_do_thermochemistry(p, xp, rt_props, cosmo, hydro_props, phys_const, us,
                          0.5 * dt);
    return;
  }

  /* If we're good, update the particle data from grackle results */
  hydro_set_internal_energy(p, u_new);

  /* update radiation fields */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float e_old = p->rt_data.radiation[g].energy_density;
    const float factor_new = (1.f - dt * iact_rates_by_frequency_bin[g]);
    p->rt_data.radiation[g].energy_density *= factor_new;
    for (int i = 0; i < 3; i++) {
      p->rt_data.radiation[g].flux[i] *= factor_new;
    }
    rt_check_unphysical_state(&p->rt_data.radiation[g].energy_density,
                              p->rt_data.radiation[g].flux, e_old,
                              /*callloc=*/2);
  }

  /* copy updated grackle data to particle */
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
static float rt_tchem_get_tchem_time(
    const struct part* restrict p, const struct xpart* restrict xp,
    struct rt_props* rt_props, const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us) {
  /* Note: Can't pass rt_props as const struct because of grackle
   * accessinging its properties there */

  /* initialize data so it'll be in scope */
  grackle_field_data particle_grackle_data;

  int grid_dimension[GRACKLE_RANK] = {GRACKLE_NPART, 1, 1};
  int grid_start[GRACKLE_RANK] = {0, 0, 0};
  int grid_end[GRACKLE_RANK] = {GRACKLE_NPART - 1, 0, 0};

  gr_float density = hydro_get_physical_density(p, cosmo);
  const float u_minimal = hydro_props->minimal_internal_energy;
  gr_float internal_energy =
      max(hydro_get_physical_internal_energy(p, xp, cosmo), u_minimal);

  gr_float species_densities[6];
  rt_tchem_get_species_densities(p, density, species_densities);

  gr_float iact_rates[5];
  float iact_rates_by_frequency_bin[RT_NGROUPS];
  rt_tchem_get_interaction_rates(iact_rates, iact_rates_by_frequency_bin, p,
                                 species_densities, rt_props, phys_const, us,
                                 cosmo);

  rt_tchem_copy_data_to_grackle(
      &particle_grackle_data, grid_dimension, grid_start, grid_end, &density,
      &internal_energy, species_densities, iact_rates);

  /* Compute 'cooling' time */
  /* Note: grackle_rates is a global variable defined by grackle itself.
   * Using a manually allocd and initialized variable here fails with MPI
   * for some reason. */
  gr_float tchem_time;
  if (local_calculate_cooling_time(&rt_props->grackle_chemistry_data,
                                   &grackle_rates, &rt_props->grackle_units,
                                   &particle_grackle_data, &tchem_time) == 0)
    error("Error in calculate_cooling_time.");

  return (float)tchem_time;
}

#endif /* SWIFT_RT_GEAR_THERMOCHEMISTRY_H */
