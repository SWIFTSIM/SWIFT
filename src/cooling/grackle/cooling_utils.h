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

#include <math.h>

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
cooling_get_equilibrium_mean_molecular_weight(
    const float u, const struct phys_const *phys_const,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling) {

  const double m_H = phys_const->const_proton_mass;

  /* Grackle mode 0: Use temperature-based molecular weight calculation */
  const double k_B = phys_const->const_boltzmann_k;
  const double H_frac = cooling->HydrogenFractionByMass;

  /* Internal energy and temperature-to-mean molecular weight calculation for
   * mode 0 */
  const double T_over_mu = (hydro_gamma_minus_one * u * m_H) / k_B;

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
 * Compute gas mean molecular weight.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Mean molecular weight.
 */
__attribute__((always_inline)) INLINE static double
cooling_get_mean_molecular_weight(const struct phys_const *phys_const,
                                  const struct unit_system *us,
                                  const struct cosmology *cosmo,
                                  const struct hydro_props *hydro_props,
                                  const struct cooling_function_data *cooling,
                                  const struct part *p,
                                  const struct xpart *xp) {

  /* Grackle mode 0: Use temperature-based molecular weight calculation */
#if COOLING_GRACKLE_MODE == 0
  const double u = hydro_get_drifted_physical_internal_energy(p, cosmo);
  const double mu = cooling_get_equilibrium_mean_molecular_weight(
      u, phys_const, hydro_props, cooling);
  return mu;

#elif COOLING_GRACKLE_MODE >= 1
  /* HI, HII, HeI, HeII, HeIII are tracked in every mode >= 1 -- shared by
     modes 1-3 instead of re-declared (and, for modes 2-3, previously
     forgotten) in each one. */
  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);
  const double m_H = phys_const->const_proton_mass;

  const double XHI = cool_data->HI_frac;
  const double XHII = cool_data->HII_frac;
  const double XHeI = cool_data->HeI_frac;
  const double XHeII = cool_data->HeII_frac;
  const double XHeIII = cool_data->HeIII_frac;

  const double nHI = XHI * rho / m_H;
  const double nHII = XHII * rho / m_H;
  const double nHeI = XHeI * rho / (4 * m_H);  // He is ~4 times heavier than H
  const double nHeII = XHeII * rho / (4 * m_H);
  const double nHeIII = XHeIII * rho / (4 * m_H);

  /* Match Grackle's own mu (calculate_temperature.c: MU_METAL = 16). */
  const double MU_METAL = 16.0;
  const double n_metal =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) * rho /
      (MU_METAL * m_H);

#if COOLING_GRACKLE_MODE == 1
  const double nel = nHII + nHeII + 2 * nHeIII;
  const double total_density =
      nHI + nHII + nHeI + nHeII + nHeIII + nel + n_metal;
  const double mu =
      ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4) / total_density;
  return mu;

#else /* COOLING_GRACKLE_MODE >= 2: also track H2I, H2II, HM */
  const double XH2I = cool_data->H2I_frac;
  const double XH2II = cool_data->H2II_frac;
  const double nH2I = XH2I * rho / (2 * m_H);  // H2 is 2 times the mass of H
  const double nH2II = XH2II * rho / (2 * m_H);

  const double XHM = cool_data->HM_frac;
  const double nHM = XHM * rho / m_H;  // HM (H-) has ~the mass of H

#if COOLING_GRACKLE_MODE == 2
  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;
  const double total_density =
      nHI + nHII + nHeI + nHeII + nHeIII + nH2I + nH2II + nHM + nel + n_metal;
  const double mu =
      ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHM) /
      total_density;
  return mu;

#else  /* COOLING_GRACKLE_MODE == 3: also track HDI */
  const double XHDI = cool_data->HDI_frac;
  const double nHDI = XHDI * rho / (3 * m_H);  // HD is 3 times the mass of H

  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;
  const double total_density = nHI + nHII + nHeI + nHeII + nHeIII + nH2I +
                               nH2II + nHM + nHDI + nel + n_metal;
  const double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 +
                     (nH2I + nH2II) * 2 + nHM + nHDI * 3) /
                    total_density;
  return mu;
#endif /* COOLING_GRACKLE_MODE == 3 */
#endif /* COOLING_GRACKLE_MODE >= 2 */
#endif /* COOLING_GRACKLE_MODE >= 1 */

#if COOLING_GRACKLE_MODE < 0 || COOLING_GRACKLE_MODE > 3
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

#if COOLING_GRACKLE_MODE >= 2
/**
 * @brief Grackle's H2 rotational/vibrational effective adiabatic index
 * (calculate_pressure.c, "Correct for Gamma from H2"), given the H2 and
 * non-H2 number densities and a seed temperature to evaluate the
 * vibrational term at.
 *
 * Shared by the forward (internal energy -> temperature) and inverse
 * (temperature -> internal energy) conversions, so the two stay exact
 * inverses of each other instead of drifting apart if edited separately.
 *
 * @param nH2 Number density of H2 (H2I + H2II).
 * @param number_density_noH2 Number density of everything except H2 that
 *   enters Grackle's own Gamma correction (see calculate_pressure.c).
 * @param T_seed Temperature to evaluate the vibrational term at -- the
 *   forward direction seeds this with its own fixed-gamma estimate, the
 *   inverse direction already knows the exact target temperature.
 * @return Grackle's effective Gamma1.
 */
__attribute__((always_inline)) INLINE static double cooling_h2_effective_gamma(
    const double nH2, const double number_density_noH2, const double T_seed) {

  const double T_safe = T_seed > 1.0 ? T_seed : 1.0;

  /* Rotational-only default; only refine with the vibrational term if H2
     is non-trace and the gas isn't so cold the mode is frozen out
     (Grackle's own x < 10 guard -- avoids exp() overflow at low T). */
  double GammaH2Inverse = 0.5 * 5.0;
  if (nH2 > 0.0 && number_density_noH2 > 0.0 &&
      nH2 / number_density_noH2 > 1e-3) {
    const double x = 6100.0 / T_safe;
    if (x < 10.0) {
      const double ex = exp(x);
      GammaH2Inverse =
          0.5 * (5.0 + 2.0 * x * x * ex / ((ex - 1.0) * (ex - 1.0)));
    }
  }

  const double GammaInverse = hydro_one_over_gamma_minus_one;
  return 1.0 + (nH2 + number_density_noH2) /
                   (nH2 * GammaH2Inverse + number_density_noH2 * GammaInverse);
}

/**
 * @brief compute the gas temperature with Grackle's H2 rotational/
 * vibrational effective-gamma correction (calculate_pressure.c,
 * "Correct for Gamma from H2"), which
 * cooling_temperature_from_internal_energy() does not apply.
 *
 * Self-contained: computes mu and the Gamma1 correction from the same
 * number densities in one pass, rather than calling
 * cooling_get_mean_molecular_weight() and recomputing them a second time.
 * The fixed-gamma temperature from that mu is Grackle's own "default
 * Gamma" estimate (verified to agree with Grackle to ~1e-4% whenever H2 is
 * trace), so it doubles as the seed for Grackle's own iteration. The
 * corrected temperature is then that seed rescaled by (Gamma1-1)/(gamma-1),
 * exactly as Grackle rescales its pressure array.
 *
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param u The particle's (already-drifted, physical) specific internal
 *   energy, so callers that already have it don't pay for it twice.
 * @return Temperature, H2-Gamma corrected.
 */
__attribute__((always_inline)) INLINE static float
cooling_get_temperature_h2_gamma_corrected(const struct phys_const *phys_const,
                                           const struct cosmology *cosmo,
                                           const struct part *p,
                                           const struct xpart *xp,
                                           const double u) {

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  const double nHI = cool_data->HI_frac * rho / m_H;
  const double nHII = cool_data->HII_frac * rho / m_H;
  const double nHeI = cool_data->HeI_frac * rho / (4 * m_H);
  const double nHeII = cool_data->HeII_frac * rho / (4 * m_H);
  const double nHeIII = cool_data->HeIII_frac * rho / (4 * m_H);
  const double nHM = cool_data->HM_frac * rho / m_H;
  const double nH2I = cool_data->H2I_frac * rho / (2 * m_H);
  const double nH2II = cool_data->H2II_frac * rho / (2 * m_H);
  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;
  const double nH2 = nH2I + nH2II;

  const double MU_METAL = 16.0;
  const double n_metal =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) * rho /
      (MU_METAL * m_H);

  /* Grackle's "number_density" local to the Gamma correction -- everything
     except H2 (tracked separately as nH2). */
  const double number_density_noH2 =
      nHI + nHII + nHeI + nHeII + nHeIII + nHM + nel;

#if COOLING_GRACKLE_MODE == 2
  const double total_density = number_density_noH2 + nH2 + n_metal;
  const double mu =
      ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHM) /
      total_density;
#else /* COOLING_GRACKLE_MODE == 3: also track HDI */
  const double nHDI = cool_data->HDI_frac * rho / (3 * m_H);
  const double total_density = number_density_noH2 + nH2 + nHDI + n_metal;
  const double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 +
                     (nH2I + nH2II) * 2 + nHM + nHDI * 3) /
                    total_density;
#endif

  const double T_default =
      cooling_temperature_from_internal_energy(u, mu, k_B, m_H);
  const double Gamma1 =
      cooling_h2_effective_gamma(nH2, number_density_noH2, T_default);

  return T_default * (Gamma1 - 1.0) / hydro_gamma_minus_one;
}

/**
 * @brief invert cooling_get_temperature_h2_gamma_corrected(): the specific
 * internal energy a particle with this composition must have for
 * cooling_get_temperature() to report exactly T_target.
 *
 * Needed because cooling_agora_cmb_floor_internal_energy() sets a target
 * temperature, not the other way around -- using the plain (uncorrected)
 * inverse for MODE >= 2 would under-shoot the floor for H2-rich gas, since
 * cooling_get_temperature() would then report a lower, Gamma1-corrected
 * value for that same energy instead of T_target.
 *
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param T_target Desired temperature, e.g. the AGORA CMB floor.
 * @return Specific internal energy that reproduces T_target under
 *   cooling_get_temperature_h2_gamma_corrected().
 */
__attribute__((always_inline)) INLINE static double
cooling_get_internal_energy_h2_gamma_corrected(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct part *p, const struct xpart *xp, const double T_target) {

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  const double nHI = cool_data->HI_frac * rho / m_H;
  const double nHII = cool_data->HII_frac * rho / m_H;
  const double nHeI = cool_data->HeI_frac * rho / (4 * m_H);
  const double nHeII = cool_data->HeII_frac * rho / (4 * m_H);
  const double nHeIII = cool_data->HeIII_frac * rho / (4 * m_H);
  const double nHM = cool_data->HM_frac * rho / m_H;
  const double nH2I = cool_data->H2I_frac * rho / (2 * m_H);
  const double nH2II = cool_data->H2II_frac * rho / (2 * m_H);
  const double nel = nHII + nHeII + 2 * nHeIII + nH2II;
  const double nH2 = nH2I + nH2II;

  const double MU_METAL = 16.0;
  const double n_metal =
      chemistry_get_total_metal_mass_fraction_for_cooling(p) * rho /
      (MU_METAL * m_H);

  const double number_density_noH2 =
      nHI + nHII + nHeI + nHeII + nHeIII + nHM + nel;

#if COOLING_GRACKLE_MODE == 2
  const double total_density = number_density_noH2 + nH2 + n_metal;
  const double mu =
      ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 + (nH2I + nH2II) * 2 + nHM) /
      total_density;
#else /* COOLING_GRACKLE_MODE == 3: also track HDI */
  const double nHDI = cool_data->HDI_frac * rho / (3 * m_H);
  const double total_density = number_density_noH2 + nH2 + nHDI + n_metal;
  const double mu = ((nHI + nHII) + (nHeI + nHeII + nHeIII) * 4 +
                     (nH2I + nH2II) * 2 + nHM + nHDI * 3) /
                    total_density;
#endif

  /* T_target is already exact, so it is a better Gamma1 seed than the
     forward direction's own fixed-gamma approximation gets to use. */
  const double Gamma1 =
      cooling_h2_effective_gamma(nH2, number_density_noH2, T_target);
  const double u_default =
      cooling_internal_energy_from_T(T_target, mu, k_B, m_H);

  return u_default * hydro_gamma_minus_one / (Gamma1 - 1.0);
}
#endif /* COOLING_GRACKLE_MODE >= 2 */

/**
 * @brief compute the AGORA redshift-dependent CMB-floor specific internal
 * energy for a particle, from its current composition.
 *
 * Call again after the composition changes (e.g. a Grackle chemistry solve)
 * to keep it consistent with cooling_get_temperature()'s later mu.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Specific internal energy corresponding to T_CMB,0 * (1+z).
 */
__attribute__((always_inline)) INLINE static double
cooling_agora_cmb_floor_internal_energy(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp) {

  const double z = (cooling->redshift == -1) ? cosmo->z : cooling->redshift;
  const double T_CMB_agora =
      CMB_TEMPERATURE_AT_REDSHIFT_0_IN_KELVIN * (z + 1.0);

#if COOLING_GRACKLE_MODE >= 2
  return cooling_get_internal_energy_h2_gamma_corrected(phys_const, cosmo, p,
                                                        xp, T_CMB_agora);
#else
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

  return cooling_internal_energy_from_T(T_CMB_agora, mu, k_B, m_H);
#endif
}
#endif /* SWIFT_COOLING_GRACKLE_COOLING_UTILS_H */
