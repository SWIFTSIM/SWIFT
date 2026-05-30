/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_COOLING_GRACKLE_COOLING_UTILS_H
#define SWIFT_COOLING_GRACKLE_COOLING_UTILS_H
/**
 * @file src/cooling/grackle/cooling_utils.h
 * @brief Cooling utilities functions for grackle.
 */

#include "chemistry.h"
#include "cooling_properties.h"
#include "hydro.h"
#include "units.h"

/**
 * Compute gas mean molecular weight.
 *
 * @param u Internal energy in physical units
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param hydro_properties The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @return Mean molecular weight.
 */
__attribute__((always_inline)) INLINE static double
cooling_get_equilibrium_mean_molecular_weight(const float u, const struct phys_const *phys_const,
                                  const struct hydro_props *hydro_props,
                                  const struct cooling_function_data *cooling) {

  const double m_H = phys_const->const_proton_mass;

  /* Grackle mode 0: Use temperature-based molecular weight calculation */
  const double k_B = phys_const->const_boltzmann_k;
  const double H_frac = cooling->HydrogenFractionByMass;

  /* Internal energy and temperature-to-mean molecular weight calculation for
   * mode 0 */
  const double T_over_mu =
      (hydro_gamma_minus_one * u * m_H) / k_B;

  const double T_transition = hydro_props->hydrogen_ionization_temperature;
  const double mu_neutral = hydro_props->mu_neutral;
  const double mu_ionised = hydro_props->mu_ionised;
  const double mu_transition = 4.0 / (8.0 - 5.0 * (1.0 - H_frac));

  double mu = 0;

  /* Are we above or below the HII -> HI transition? */
  if (T_over_mu > (T_transition + 1.0) / mu_ionised) {
    mu = mu_ionised;
  } else if (T_over_mu < (T_transition - 1.) / mu_neutral) {
    mu = mu_neutral;
  } else {
    mu = mu_transition;
  }
  return mu;
}

/**
 * @brief compute the (physical) specific internal energy of an ideal gas for
 * given temperature and mean molecular weight.
 *
 * @param T Temperature of the gas.
 * @param mu Mean molecular weight of the gas.
 * @param kB Boltzmann constant.
 * @param mp Proton mass.
 * */
__attribute__((always_inline)) INLINE static float
cooling_internal_energy_from_T(const double T, const double mu, const double kB,
                               const double mp) {
  return kB * T * hydro_one_over_gamma_minus_one / (mu * mp);
}

/**
 * @brief compute the temperature of an ideal gas for a given specific internal
 * energy and mean molecular weight.
 *
 * @param u Specific internal energy of the gas.
 * @param mu Mean molecular weight of the gas.
 * @param kB Boltzmann constant.
 * @param mp Proton mass.
 * */
__attribute__((always_inline)) INLINE static float
cooling_temperature_from_internal_energy(const double u, const double mu,
                                         const double kB, const double mp) {
  return u * hydro_gamma_minus_one * mu * mp / kB;
}
#endif /* SWIFT_COOLING_GRACKLE_COOLING_UTILS_H */
