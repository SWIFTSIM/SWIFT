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
 * @file src/cooling/grackle/cooling.c
 * @brief Cooling using the GRACKLE 3.1.1 library.
 */

#include "hydro.h"
#include "cooling_properties.h"
#include "chemistry.h"

/**
 * Compute the gas hydrogen mass fraction.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static float cooling_get_hydrogen_mass_fraction(
										      const struct cooling_function_data* cooling,
										      const struct part* p, const struct xpart* xp) {

#if COOLING_GRACKLE_MODE == 0
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Approximation that takes into account metals */
  return cooling->HydrogenFractionByMass - Z;
#else

  const struct cooling_xpart_data* cool_data = &xp->cooling_data;

  /* Mode 1-3 have at least thes constributions */
  float X_H = cool_data->HI_frac + cool_data->HII_frac;

#if COOLING_GRACKLE_MODE > 1
  X_H += cool_data->HM_frac + cool_data->H2I_frac + cool_data->H2II_frac;
#endif

  return X_H;
#endif
}

/**
 * Compute gas mean molecular weight.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_properties The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Mean molecular weight.
 */
__attribute__((always_inline)) INLINE static double cooling_get_mean_molecular_weight(const struct phys_const* phys_const, const struct unit_system* us,   const struct cosmology* cosmo,  const struct hydro_props* hydro_props,
const struct cooling_function_data* cooling, const struct part* p, const struct xpart *xp) {

  const double m_H = phys_const->const_proton_mass;
  
  /* Grackle mode 0: Use temperature-based molecular weight calculation */
#if COOLING_GRACKLE_MODE == 0
  const double k_B = phys_const->const_boltzmann_k;
  const double H_frac = cooling->HydrogenFractionByMass;  // Hydrogen fraction from metadata
    
  /* Internal energy and temperature-to-mean molecular weight calculation for mode 0 */
  const double u =  hydro_get_drifted_physical_internal_energy(p, cosmo);
  const double T_over_mu = (hydro_gamma_minus_one * u * m_H) / k_B;  // Adjust this as needed

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

  /* Grackle mode 1: Only HI, HII, HeI, HeII, and HeIII species */
#elif COOLING_GRACKLE_MODE == 1
  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);

  /* Extract mass fractions for various species from the cooling data */
  const double XHI    = cool_data->HI_frac;
  const double XHII   = cool_data->HII_frac;
  const double XHeI   = cool_data->HeI_frac;
  const double XHeII  = cool_data->HeII_frac;
  const double XHeIII = cool_data->HeIII_frac;
  
  const double nHI   = XHI * rho / m_H;
  const double nHII  = XHII * rho / m_H;
  const double nHeI  = XHeI * rho / (4 * m_H);  // He is ~4 times heavier than H
  const double nHeII = XHeII * rho / (4 * m_H);
  const double nHeIII = XHeIII * rho / (4 * m_H);

  /* Calculate total number of electrons (nel) */
  const double nel = nHII + nHeII + 2 * nHeIII;

  /* Compute the molecular weight using the species mass fractions */
  const double total_density = nHI + nHII + nHeI + nHeII + nHeIII + nel;
  const double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4) / total_density;

  return mu;

  /* Grackle mode 2: Add species like H2I and H2II to the molecular weight calculation */
#elif COOLING_GRACKLE_MODE == 2
  const double XHeI   = cool_data->HeI_frac;
  const double XHeII  = cool_data->HeII_frac;
  const double XHeIII = cool_data->HeIII_frac;
  const double XH2I   = cool_data->H2I_frac;
  const double XH2II  = cool_data->H2II_frac;
  const double XHI    = cool_data->HI_frac;
  const double XHII   = cool_data->HII_frac;

  const double nHI   = XHI * rho / m_H;
  const double nHII  = XHII * rho / m_H;
  const double nHeI  = XHeI * rho / (4 * m_H);  // He is ~4 times heavier than H
  const double nHeII = XHeII * rho / (4 * m_H);
  const double nHeIII = XHeIII * rho / (4 * m_H);

  const double nH2I  = XH2I * rho / (2 * m_H);  // H2 is 2 times the mass of H
  const double nH2II = XH2II * rho / (2 * m_H);

  /* Calculate total number of electrons (nel) */
  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;

  /* Compute the molecular weight using the species mass fractions */
  const double total_density = nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nel;
  double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2) / total_density;

  return mu;

  /* Grackle mode 3: Add species like H2I, H2II, HDI, and electrons to the molecular weight calculation */
#elif COOLING_GRACKLE_MODE == 3
  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);
  const double XHeI   = cool_data->HeI_frac;
  const double XHeII  = cool_data->HeII_frac;
  const double XHeIII = cool_data->HeIII_frac;
  const double XH2I   = cool_data->H2I_frac;
  const double XH2II  = cool_data->H2II_frac;
  const double XHI    = cool_data->HI_frac;
  const double XHII   = cool_data->HII_frac;
  const double XHDI   = cool_data->HDI_frac;
  
  const double nHI   = XHI * rho / m_H;
  const double nHII  = XHII * rho / m_H;
  const double nHeI  = XHeI * rho / (4 * m_H);  // He is ~4 times heavier than H
  const double nHeII = XHeII * rho / (4 * m_H);
  const double nHeIII = XHeIII * rho / (4 * m_H);

  const double nH2I  = XH2I * rho / (2 * m_H);  // H2 is 2 times the mass of H
  const double nH2II = XH2II * rho / (2 * m_H);

  const double nHDI  = XHDI * rho / (3 * m_H);  // HD is 3 times the mass of H

  /* Calculate total number of electrons (nel) */
  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;

  /* Compute the molecular weight using the species mass fractions */
  const double total_density = nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nHDI + nel;
  double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHDI * 3) / total_density;

  return mu;
#else
  /* Default case: Error if cooling model is not valid */
#error "Invalid COOLING_GRACKLE_MODE"
#endif
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
__attribute__((always_inline)) INLINE static float cooling_internal_energy_from_T(
    const double T, const double mu,
    const double kB, const double mp) {
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
__attribute__((always_inline)) INLINE static float cooling_temperature_from_internal_energy(
    const double u, const double mu,
    const double kB, const double mp) {
  return u * hydro_gamma_minus_one * mu * mp / kB;
}
#endif /* SWIFT_COOLING_GRACKLE_COOLING_UTILS_H */
