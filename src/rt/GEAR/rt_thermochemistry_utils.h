/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_THERMOCHEMISTRY_UTILS_H
#define SWIFT_RT_THERMOCHEMISTRY_UTILS_H

/**
 * @file src/rt/GEAR/rt_thermochemistry_utils.h
 * @brief thermochemistry utilities and misc functions.
 * */

/**
 * @brief compute the mean molecular weight mu for given
 * hydrogen and helium mass fractions
 * @param XHI mass fraction of HI
 * @param XHII mass fraction of HII
 * @param XHeI mass fraction of HeI
 * @param XHeII mass fraction of HeII
 * @param XHeIII mass fraction of HeII
 */
__attribute__((always_inline)) INLINE static double
rt_tchem_get_mean_molecular_weight(float XHI, float XHII, float XHeI,
                                   float XHeII, float XHeIII) {

  /* 1/mu = sum_j X_j / A_j * (1 + E_j)
   * A_H    = 1, E_H    = 0
   * A_Hp   = 1, E_Hp   = 1
   * A_He   = 4, E_He   = 0
   * A_Hep  = 4, E_Hep  = 1
   * A_Hepp = 4, E_Hepp = 2             */
  const double one_over_mu =
      XHI + 2.0 * XHII + 0.25 * XHeI + 0.5 * XHeII + 0.75 * XHeIII;

  return (1.0 / one_over_mu);
}

/**
 * @brief compute the temperature of an ideal gas for a given
 * specific internal energy and mean molecular weight
 *
 * @param u specific internal energy of the gas
 * @param mu mean molecular weight of the gas
 * @param kB Boltzmann constant
 * @param mp proton mass
 * */
__attribute__((always_inline)) INLINE static float
rt_tchem_temperature_from_internal_energy(double u, double mu, const double kB,
                                          const double mp) {

  return u * hydro_gamma_minus_one * mu * mp / kB;
}

/**
 * @brief compute the (physical) specific internal energy
 * of an ideal gas for given temperature and mean molecular
 * weight.
 *
 * @param u specific internal energy of the gas
 * @param mu mean molecular weight of the gas
 * @param kB Boltzmann constant
 * @param mp proton mass
 * */
__attribute__((always_inline)) INLINE static float
rt_tchem_internal_energy_from_T(const double T, const double mu,
                                const double kB, const double mp) {

  return kB * T * hydro_one_over_gamma_minus_one / mu / mp;
}

/**
 * @brief compute the derivative w.r.t. temperature of the
 * (physical) specific internal energy of an ideal gas for
 * given temperature and mean molecular weight
 *
 * @param u specific internal energy of the gas
 * @param mu mean molecular weight of the gas
 * @param kB Boltzmann constant
 * @param mp proton mass
 * */
__attribute__((always_inline)) INLINE static float rt_tchem_internal_energy_dT(
    double mu, const double kB, const double mp) {

  const double dudT = kB * hydro_one_over_gamma_minus_one / mu / mp;

  return dudT;
}

#endif /* SWIFT_RT_THERMOCHEMISTRY_UTILS_H */
