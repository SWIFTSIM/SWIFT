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
#ifndef SWIFT_COOLING_GEAJAY_COOLING_UTILS_H
#define SWIFT_COOLING_GEAJAY_COOLING_UTILS_H
/**
 * @file src/cooling/geajay/cooling_utils.h
 * @brief Cooling utilities functions for geajay.
 */

#include "chemistry.h"
#include "cooling_properties.h"
#include "error.h"
#include "hydro.h"
#include "units.h"

#include <math.h>

/**
 * Compute the Hydrogen mass fraction.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static float
cooling_get_hydrogen_mass_fraction(const struct cooling_function_data *cooling,
                                   const struct part *p,
                                   const struct xpart *xp) {

#if COOLING_GEAJAY_MODE == 0
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Approximation that takes into account metals */
  return cooling->HydrogenFractionByMass - Z;
#else

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const float X_H = cool_data->HI_frac + cool_data->HII_frac;
  return X_H;
#endif
}

/**
 * Compute the Helium mass fraction.
 *
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Hydrogen mass fraction.
 */
__attribute__((always_inline)) INLINE static float
cooling_get_helium_mass_fraction(const struct cooling_function_data *cooling,
                                 const struct part *p, const struct xpart *xp) {

#if COOLING_GEAJAY_MODE == 0
  const float Z = chemistry_get_total_metal_mass_fraction_for_cooling(p);

  /* Approximation that takes into account metals */
  const float X_H = cooling->HydrogenFractionByMass - Z;
  return 1 - X_H;
#else

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const float X_He =
      cool_data->HeI_frac + cool_data->HeII_frac + cool_data->HeIII_frac;
  return X_He;
#endif
}

/**
 * Compute gas mean molecular weight.
 *
 * @TODO: This function may be moved into a future improved Grackle version.
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
 * @brief Compute the mean molecular weight of a particle, as geajay uses it.
 *
 * Calls geajay_particle_mean_molecular_weight() at the drifted state. Stops the
 * run if geajay refuses the particle.
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
double cooling_get_mean_molecular_weight(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double cooling_get_mean_molecular_weight_at_energy(
    const struct cosmology *cosmo, const struct cooling_function_data *cooling,
    const struct part *p, const struct xpart *xp, const double u);

/**
 * @brief compute the (physical) specific internal energy of an ideal gas for
 * given temperature and mean molecular weight.
 *
 * mu == 0 is a legitimate output for a particle with no tracked species at
 * all, and would turn this division into +-Inf, so floor it here once
 * rather than make every caller know that.
 *
 * @param T Temperature of the gas.
 * @param mu Mean molecular weight of the gas.
 * @param kB Boltzmann constant.
 * @param mp Proton mass.
 * */
__attribute__((always_inline)) INLINE static float
cooling_internal_energy_from_T(const double T, const double mu, const double kB,
                               const double mp) {
  const double mu_safe = mu > 0.0 ? mu : 1.e-20;
  return kB * T * hydro_one_over_gamma_minus_one / (mu_safe * mp);
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

#if COOLING_GEAJAY_MODE >= 2
/**
 * @brief The H2 and non-H2 number densities of the H2 Gamma correction.
 */
struct cooling_h2_species_densities {
  double nH2;
  double number_density_noH2;
};

/**
 * @brief compute the H2 and non-H2 number densities that the H2 Gamma
 * correction of geajay uses.
 *
 * The non-H2 species are those of the geajay solver: H, He and the electron
 * field of the particle, without H-. Floors the denominator to a tiny number
 * so that a particle with every species fraction at zero gives no 0/0.
 *
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return nH2 and number_density_noH2 for this particle.
 */
__attribute__((always_inline)) INLINE static struct cooling_h2_species_densities
cooling_get_h2_species_densities(const struct phys_const *phys_const,
                                 const struct cosmology *cosmo,
                                 const struct part *p, const struct xpart *xp) {

  const struct cooling_xpart_data *cool_data = &xp->cooling_data;
  const double rho = hydro_get_physical_density(p, cosmo);
  const double m_H = phys_const->const_proton_mass;

  const double nHI = cool_data->HI_frac * rho / m_H;
  const double nHII = cool_data->HII_frac * rho / m_H;
  const double nHeI = cool_data->HeI_frac * rho / (4 * m_H);
  const double nHeII = cool_data->HeII_frac * rho / (4 * m_H);
  const double nHeIII = cool_data->HeIII_frac * rho / (4 * m_H);
  const double nel = cool_data->e_frac * rho / m_H;
  const double nH2 =
      (cool_data->H2I_frac + cool_data->H2II_frac) * rho / (2 * m_H);

  const double tiny_number = 1.e-20;
  double number_density_noH2 = nHI + nHII + nHeI + nHeII + nHeIII + nel;
  if (number_density_noH2 == 0.0) number_density_noH2 = tiny_number;

  return (struct cooling_h2_species_densities){
      .nH2 = nH2, .number_density_noH2 = number_density_noH2};
}

/**
 * @brief Grackle's H2 rotational/vibrational effective adiabatic index
 * (calculate_pressure.c, "Correct for Gamma from H2").
 *
 * Used by the inverse conversion only: the temperature itself comes from
 * geajay_particle_temperature(), which applies the same correction.
 *
 * @param nH2 Number density of H2 (H2I + H2II).
 * @param number_density_noH2 Number density of everything except H2 that
 *   enters Grackle's own Gamma correction (see calculate_pressure.c).
 *   Already floored to a tiny positive value by
 *   cooling_get_h2_species_densities(), so the ratio below never divides
 *   by zero.
 * @param T_seed Temperature to evaluate the vibrational term at.
 * @return Grackle's effective Gamma1.
 */
__attribute__((always_inline)) INLINE static double cooling_h2_effective_gamma(
    const double nH2, const double number_density_noH2, const double T_seed) {

  const double T_safe = T_seed > 1.0 ? T_seed : 1.0;

  /* Rotational-only default. Refine with the vibrational term only if H2
     is non-trace and the gas isn't cold enough to freeze the mode out
     (Grackle's own x < 10 guard avoids exp() overflow at low T). */
  double GammaH2Inverse = 0.5 * 5.0;
  if (nH2 / number_density_noH2 > 1e-3) {
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
 * @brief invert the H2 Gamma correction of geajay_particle_temperature(): the
 * specific internal energy a particle with this composition must have for
 * cooling_get_temperature() to report T_target.
 *
 * Needed because cooling_agora_cmb_floor_internal_energy() sets a target
 * temperature, not the other way around, and geajay has no inverse. The plain
 * (uncorrected) inverse would under-shoot the floor for H2-rich gas. The
 * vibrational term is evaluated at T_target, not at the uncorrected
 * temperature of the forward direction.
 *
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param mu The mean molecular weight of the particle.
 * @param T_target Desired temperature, e.g. the AGORA CMB floor.
 * @return Specific internal energy that reproduces T_target.
 */
__attribute__((always_inline)) INLINE static double
cooling_get_internal_energy_h2_gamma_corrected(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct part *p, const struct xpart *xp, const double mu,
    const double T_target) {

  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  const struct cooling_h2_species_densities sd =
      cooling_get_h2_species_densities(phys_const, cosmo, p, xp);

  const double Gamma1 =
      cooling_h2_effective_gamma(sd.nH2, sd.number_density_noH2, T_target);
  const double u_default =
      cooling_internal_energy_from_T(T_target, mu, k_B, m_H);

  return u_default * hydro_gamma_minus_one / (Gamma1 - 1.0);
}
#endif /* COOLING_GEAJAY_MODE >= 2 */

/**
 * @brief compute the specific internal energy a particle with this
 * composition must have for cooling_get_temperature() to report T_target: the
 * T-to-u counterpart of cooling_get_temperature() itself.
 *
 * In mode 0 the mean molecular weight depends on the energy through the Cloudy
 * table, so the energy is found by fixed-point iteration. The result is exact
 * where mu is continuous. Where mu jumps, no energy may give T_target. The
 * iteration then alternates between two energies and the larger one is
 * returned, so T is up to about 0.5% above T_target and never below it. The
 * run stops only if the energy is not positive.
 *
 * Every caller that needs to turn a target temperature into the u that
 * actually produces it (currently just the AGORA CMB floor, potentially
 * others later) should go through this one function rather than repeat
 * the COOLING_GEAJAY_MODE dispatch itself.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param T_target Desired temperature.
 * @return Specific internal energy that reproduces T_target under
 *   cooling_get_temperature().
 */
__attribute__((always_inline)) INLINE static double
cooling_get_internal_energy_from_temperature(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp, const double T_target) {

#if COOLING_GEAJAY_MODE == 0
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  double u = hydro_get_drifted_physical_internal_energy(p, cosmo);
  double u_prev = u;
  for (int i = 0; i < 20; i++) {
    const double mu =
        cooling_get_mean_molecular_weight_at_energy(cosmo, cooling, p, xp, u);
    const double u_next =
        cooling_internal_energy_from_T(T_target, mu, k_B, m_H);
    if (!(u_next > 0.))
      error("Particle %lld: invalid floor energy %e for T = %e K.", p->id,
            u_next, T_target);
    if (fabs(u_next - u) <= 1.e-6 * u_next) return u_next;
    /* mu(u) has jumps: the iteration can alternate between two energies. */
    const int cycle = (i > 0 && fabs(u_next - u_prev) <= 1.e-6 * u_next);
    u_prev = u;
    u = u_next;
    if (cycle) break;
  }
  /* Keep the energy that gives T at or above T_target. */
  return fmax(u, u_prev);
#else
  const double mu = cooling_get_mean_molecular_weight(
      phys_const, us, cosmo, hydro_props, cooling, p, xp);

#if COOLING_GEAJAY_MODE >= 2
  return cooling_get_internal_energy_h2_gamma_corrected(phys_const, cosmo, p,
                                                        xp, mu, T_target);
#else
  const double m_H = phys_const->const_proton_mass;
  const double k_B = phys_const->const_boltzmann_k;

  return cooling_internal_energy_from_T(T_target, mu, k_B, m_H);
#endif
#endif
}

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
  const double T_CMB_agora = phys_const->const_T_CMB_0 * (z + 1.0);

  return cooling_get_internal_energy_from_temperature(
      phys_const, us, cosmo, hydro_props, cooling, p, xp, T_CMB_agora);
}
#endif /* SWIFT_COOLING_GEAJAY_COOLING_UTILS_H */
