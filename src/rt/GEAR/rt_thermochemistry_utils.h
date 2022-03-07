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
#ifndef SWIFT_RT_GEAR_THERMOCHEMISTRY_UTILS_H
#define SWIFT_RT_GEAR_THERMOCHEMISTRY_UTILS_H

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

  return kB * T * hydro_one_over_gamma_minus_one / (mu * mp);
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

  const double dudT = kB * hydro_one_over_gamma_minus_one / (mu * mp);

  return dudT;
}

/**
 * @brief get the densities of all species and electrons.
 *
 * @param p particle to use
 * @param rho particle physical density
 * @param species_densities array to write densities in
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_get_species_densities(const struct part* restrict p, gr_float rho,
                               gr_float species_densities[6]) {

  /* TODO: check units here */
  species_densities[0] = p->rt_data.tchem.mass_fraction_HI * rho;
  species_densities[1] = p->rt_data.tchem.mass_fraction_HII * rho;
  species_densities[2] = p->rt_data.tchem.mass_fraction_HeI * rho;
  species_densities[3] = p->rt_data.tchem.mass_fraction_HeII * rho;
  species_densities[4] = p->rt_data.tchem.mass_fraction_HeIII * rho;

  /* nHII = rho_HII / m_p
   * nHeII = rho_HeII / 4 m_p
   * nHeIII = rho_HeIII / 4 m_p
   * ne = nHII + nHeII + 2 * nHeIII
   * But: it is grackle convention to use rho_e = n_e * m_p */
  const gr_float rho_e = species_densities[1] + 0.25 * species_densities[3] +
                         0.5 * species_densities[4];
  species_densities[5] = rho_e;
}

/**
 * @brief get the (physical) temperature of the gas
 *
 * @param p particle to use
 * @param phys_const physical constants struct
 * @param cosmo cosmology struct
 **/
__attribute__((always_inline)) INLINE static double
rt_tchem_get_gas_temperature(const struct part* restrict p,
                             const struct phys_const* restrict phys_const,
                             const struct cosmology* restrict cosmo) {

  const double kB = phys_const->const_boltzmann_k;
  const double mp = phys_const->const_proton_mass;

  const float XHI = p->rt_data.tchem.mass_fraction_HI;
  const float XHII = p->rt_data.tchem.mass_fraction_HII;
  const float XHeI = p->rt_data.tchem.mass_fraction_HeI;
  const float XHeII = p->rt_data.tchem.mass_fraction_HeII;
  const float XHeIII = p->rt_data.tchem.mass_fraction_HeIII;

  const double mu =
      rt_tchem_get_mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII);
  const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);

  double T = rt_tchem_temperature_from_internal_energy(u, mu, kB, mp);

  return T;
}
#endif /* SWIFT_RT_GEAR_THERMOCHEMISTRY_UTILS_H */
