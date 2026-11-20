/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_ISRF_EXTINCTION_GEAR_H
#define SWIFT_RADIATION_ISRF_EXTINCTION_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_isrf_extinction.h
 * @brief Receiver-side LW/PE dust extinction of the ISRF for GEAR: extinction
 * path, dust opacity and extinction factors. Inline, because the star-gas
 * pair loop calls them for every pair.
 */

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "part.h"
#include "physical_constants.h"
#include "radiation.h"
#include "units.h"

#include <math.h>

/**
 * @brief Comoving path length of the receiver-side LW/PE dust extinction
 * column, chosen by #isrf_extinction_path_mechanism.
 *
 * The arguments are the union over all mechanisms; each case uses a subset.
 * The call is inside the star-gas pair loop. The temperature_capped_jeans
 * mechanism is costly there (one temperature and a square root per pair).
 *
 * @param fb_props Properties of the feedback scheme.
 * @param p The receiving #part.
 * @param xp The receiving #xpart (tracked species, for the temperature).
 * @param r Comoving separation of the illuminating star and the receiver.
 * @param cosmo The current cosmological model.
 * @param phys_const The physical constants.
 * @param hydro_props The hydro scheme properties.
 * @param us Unit system.
 * @param cooling The cooling function properties.
 * @return Comoving extinction path length.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_comoving_extinction_path(
    const struct feedback_props *fb_props, const struct part *p,
    const struct xpart *xp, const float r, const struct cosmology *cosmo,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cooling_function_data *cooling) {

  /* One support radius: the longest admissible path, and the fallback for a
   * degenerate gas state. */
  const float h_gas = p->h * kernel_gamma;

  switch ((enum isrf_extinction_path_mechanism)
              fb_props->ISRF_extinction_path_mechanism) {

    case isrf_extinction_path_constant_kernel_path:
      return fb_props->ISRF_extinction_path_in_kernel_radii * h_gas;

    case isrf_extinction_path_pair_separation:
      return r;

    case isrf_extinction_path_temperature_capped_jeans: {
      const float T = cooling_get_temperature(phys_const, hydro_props, us,
                                              cosmo, cooling, p, xp);
      const float rho_phys = hydro_get_physical_density(p, cosmo);

      /* Return the cap: a clamped division is unsafe under -ffast-math. */
      if (T <= 0.f || rho_phys <= 0.f) return h_gas;

      /* c_s^2 scales with T, so the cap multiplies it by T_cap/T. */
      const float cs = hydro_get_physical_soundspeed(p, cosmo);
      const float T_ratio =
          min(1.f, fb_props->ISRF_extinction_jeans_temperature_cap_K / T);
      const float cs2_capped = cs * cs * T_ratio;

      /* lambda_J = sqrt(pi c_s^2 / (G rho)), physical, then comoving. */
      const float lambda_J_phys = sqrtf(
          (float)(M_PI * cs2_capped / (phys_const->const_newton_G * rho_phys)));
      return min(lambda_J_phys * cosmo->a_inv, h_gas);
    }
  }

  error("Unknown GEARFeedback:ISRF_extinction_path mechanism %d.",
        (int)fb_props->ISRF_extinction_path_mechanism);
  return 0.f;
}

/**
 * @brief Comoving gas column density at a gas particle: its own density times
 * the extinction path.
 *
 * @param p The #part.
 * @param extinction_path Comoving path length, from
 * #radiation_get_comoving_extinction_path.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_comoving_gas_column_density_at_part(const struct part *p,
                                                  const float extinction_path) {
  return extinction_path * p->rho;
}

/**
 * @brief Band dust mass opacity in internal units (area/mass):
 * `sigma_d_band * D(Z) / (mu_H * m_H)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants (for the proton mass).
 * @param D_relative Dust-to-gas ratio relative to the Milky Way, from
 * #cooling_get_dust_to_gas_ratio_relative_to_MW.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @return Dust mass opacity, internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_mass_opacity(const struct unit_system *us,
                                const struct phys_const *phys_const,
                                float D_relative, float sigma_d_band_cgs) {
  return sigma_d_band_cgs * D_relative /
         (RADIATION_MU_H * phys_const->const_proton_mass *
          units_cgs_conversion_factor(us, UNIT_CONV_AREA));
}

/**
 * @brief Band dust extinction factor `exp(-kappa_eff * Sigma_gas_p)`.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param D_relative Dust-to-gas ratio relative to the Milky Way.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us,
                                     const struct phys_const *phys_const,
                                     float D_relative, float sigma_d_band_cgs,
                                     float Sigma_gas_p) {

  const float kappa_eff = radiation_get_dust_mass_opacity(
      us, phys_const, D_relative, sigma_d_band_cgs);
  return expf(-kappa_eff * Sigma_gas_p);
}

/**
 * @brief Local linear dust absorption rate `kappa_eff(Z) * rho` in 1/length.
 * It is local: no column density is involved.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param D_relative Dust-to-gas ratio relative to the Milky Way.
 * @param rho_p Physical gas density, internal units.
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @return Local linear dust absorption rate, internal units (1/length).
 */
__attribute__((always_inline)) INLINE static float
radiation_get_part_linear_absorption_rate(const struct unit_system *us,
                                          const struct phys_const *phys_const,
                                          float D_relative, float rho_p,
                                          float sigma_d_band_cgs) {
  return radiation_get_dust_mass_opacity(us, phys_const, D_relative,
                                         sigma_d_band_cgs) *
         rho_p;
}

/**
 * @brief Receiver-side LW/PE dust extinction factors for a gas particle.
 *
 * @param us Unit system.
 * @param phys_const Physical constants.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param cooling The cooling function properties (for the dust-to-gas
 * ratio, see #cooling_get_dust_to_gas_ratio_relative_to_MW).
 * @param extinction_path Comoving extinction path length, from
 * #radiation_get_comoving_extinction_path.
 * @param extinction (return) Extinction factor of each operator, indexed by
 * #radiation_isrf_operator.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_part_ISRF_extinction_factors(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct part *p, float Z,
    const struct cooling_function_data *cooling, const float extinction_path,
    float extinction[ISRF_OPERATOR_COUNT]) {

  const float Sigma_gas_p =
      radiation_get_comoving_gas_column_density_at_part(p, extinction_path) *
      cosmo->a2_inv;
  const float D_relative = cooling_get_dust_to_gas_ratio_relative_to_MW(
      cooling, Z, hydro_get_physical_density(p, cosmo));

  extinction[ISRF_OPERATOR_PE] = radiation_get_dust_extinction_factor(
      us, phys_const, D_relative, RADIATION_SIGMA_D_PE_CGS, Sigma_gas_p);
  extinction[ISRF_OPERATOR_LW] = radiation_get_dust_extinction_factor(
      us, phys_const, D_relative, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p);
}

#endif /* SWIFT_RADIATION_ISRF_EXTINCTION_GEAR_H */
