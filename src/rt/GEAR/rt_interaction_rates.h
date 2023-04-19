/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GEAR_INTERACTION_RATES_H
#define SWIFT_RT_GEAR_INTERACTION_RATES_H

#include "rt_parameters.h"
#include "rt_properties.h"
#include "rt_species.h"
#include "rt_thermochemistry_utils.h"

/**
 * @file src/rt/GEAR/rt_interaction_rates.h
 * @brief header file concerning photoionization and photoheating rates
 **/

/**
 * @brief compute the heating, ionization, and dissassociation rates
 * for the particle radiation field as needed by grackle.
 *
 * @param rates (return) Interaction rates for grackle. [0]: heating rate.
 * [1]: HI ionization. [2]: HeI ionization. [3]: HeII ionization.
 * [4]: H2 dissociation.
 * @param energy_density energy densities of each photon group [internal units]
 * @param species_densities physical densities of all species [internal units]
 * @param average_photon_energy mean photon energy in group, in erg
 * @param cse energy weighted photon interaction cross sections, in cm^2
 * @param csn number weighted photon interaction cross sections, in cm^2
 * @param phys_const physical constants struct
 * @param us internal units struct
 **/
__attribute__((always_inline)) INLINE static void
rt_get_interaction_rates_for_grackle(
    gr_float rates[5], float energy_density[RT_NGROUPS],
    gr_float species_densities[6],
    const double average_photon_energy[RT_NGROUPS], double **cse, double **csn,
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us) {

  rates[0] = 0.; /* Needs to be in [erg / s / cm^3 / nHI] for grackle. */
  rates[1] = 0.; /* [1 / time_units] */
  rates[2] = 0.; /* [1 / time_units] */
  rates[3] = 0.; /* [1 / time_units] */
  rates[4] = 0.; /* [1 / time_units] */

  double E_ion_cgs[rt_ionizing_species_count];
  rt_species_get_ionizing_energy(E_ion_cgs);

  /* Get some conversions and constants first. */
  const double c_cgs = rt_params.reduced_speed_of_light *
                       units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
  const double to_energy_density_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_DENSITY);
  const double inv_time_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);

  /* get species number densities in cgs */
  double ns_cgs[rt_ionizing_species_count]; /* in cm^-3 */
  rt_tchem_get_ionizing_species_number_densities(ns_cgs, species_densities,
                                                 phys_const, us);

  /* store photoionization rate for each species here */
  double ionization_rates[rt_ionizing_species_count];
  for (int s = 0; s < rt_ionizing_species_count; s++) {
    ionization_rates[s] = 0.;
  }

  for (int g = 0; g < RT_NGROUPS; g++) {

    /* Sum results for this group over all species */
    double heating_rate_group_cgs = 0.;
    const double Eg = energy_density[g] * to_energy_density_cgs;
    const double Emean_g = average_photon_energy[g];
    const double Ng = (Emean_g > 0.) ? Eg / Emean_g : 0.;

    for (int s = 0; s < rt_ionizing_species_count; s++) {
      /* All quantities here are in cgs. */
      heating_rate_group_cgs +=
          (cse[g][s] * Emean_g - E_ion_cgs[s] * csn[g][s]) * ns_cgs[s];
      ionization_rates[s] += csn[g][s] * Ng * c_cgs;
    }
    rates[0] += heating_rate_group_cgs * Ng * c_cgs;
  }

  /* Convert into correct units. */
  const double nHI = ns_cgs[rt_ionizing_species_HI];
  if (nHI > 0.)
    rates[0] /= nHI;
  else
    rates[0] = 0.;

  for (int s = 0; s < rt_ionizing_species_count; s++) {
    ionization_rates[s] /= inv_time_cgs; /* internal units T^-1 */
#ifdef SWIFT_RT_DEBUG_CHECKS
    if (ionization_rates[s] < 0.)
      error("unphysical ion rate spec %d - %.4g", s, ionization_rates[s]);
#endif
  }

  /* We're done. Write the results in correct place */
  rates[1] = ionization_rates[0];
  rates[2] = ionization_rates[1];
  rates[3] = ionization_rates[2];
  /* rates[4] = skipped for now */

#ifdef SWIFT_RT_DEBUG_CHECKS
  for (int i = 0; i < 5; i++) {
    if (rates[i] < 0.) error("unphysical rate %d %.4g", i, rates[i]);
  }
#endif
}

/**
 * @brief compute the rates at which the photons get absorbed/destroyed
 * during interactions with gas.
 *
 * @param absorption_rates (return) the energy absorption rates in
 * internal units for each photon group.
 * @param species_densities the physical densities of all traced species
 *[internal units]
 * @param average_photon_energy mean photon energy in group, in erg
 * @param csn number weighted photon interaction cross sections, in cm^2
 * @param phys_const physical constants struct
 * @param us internal units struct
 **/
__attribute__((always_inline)) INLINE static void rt_get_absorption_rates(
    double absorption_rates[RT_NGROUPS], gr_float species_densities[6],
    const double average_photon_energy[RT_NGROUPS], double **csn,
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us) {

  for (int g = 0; g < RT_NGROUPS; g++) absorption_rates[g] = 0.;

  double E_ion_cgs[rt_ionizing_species_count];
  rt_species_get_ionizing_energy(E_ion_cgs);

  /* Get some conversions and constants first. */
  const double c_cgs = rt_params.reduced_speed_of_light *
                       units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
  const double inv_time_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);

  double ns_cgs[rt_ionizing_species_count]; /* in cm^-3 */
  rt_tchem_get_ionizing_species_number_densities(ns_cgs, species_densities,
                                                 phys_const, us);

  for (int g = 0; g < RT_NGROUPS; g++) {
    for (int s = 0; s < rt_ionizing_species_count; s++) {
      absorption_rates[g] += c_cgs * csn[g][s] * ns_cgs[s];
    }
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    absorption_rates[g] /= inv_time_cgs;
#ifdef SWIFT_RT_DEBUG_CHECKS
    if (absorption_rates[g] < 0.)
      error("unphysical rate %d - %.4g", g, absorption_rates[g]);
#endif
  }
}

#endif /* SWIFT_RT_GEAR_INTERACION_RATES_H */
