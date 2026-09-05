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
/**
 * @file src/feedback/GEAR/radiation_pressure.c
 * @brief Star-side infrared radiation pressure for GEAR: the local gas
 * column density around a star (Sobolev approximation), the resulting IR
 * dust opacity and optical depth, and the radiation pressure they imply.
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "error.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "radiation.h"

#include <math.h>

/**
 * Compute the gas comoving column density at the star's location using the
 * Sobolev approximation.
 *
 * @param sp The #spart.
 * @return Comoving gas column density at the star's location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_star(const struct spart *sp) {
  /* enrichment_weight is the star's SPH-averaged local gas density. */
  const float rho_gas = sp->feedback_data.enrichment_weight;
  const float grad_rho[3] = {sp->feedback_data.grad_rho_star[0],
                             sp->feedback_data.grad_rho_star[1],
                             sp->feedback_data.grad_rho_star[2]};
  const float norm_grad_rho =
      sqrtf(grad_rho[0] * grad_rho[0] + grad_rho[1] * grad_rho[1] +
            grad_rho[2] * grad_rho[2]);

  /* Cap the Sobolev length rho/|grad rho| at the kernel support radius
     rather than letting it blow up towards infinity for a locally uniform
     density field (zero or near-zero gradient, e.g. an unperturbed
     glass/grid IC, where the raw ratio is dominated by SPH summation
     noise, not a resolved trend). Capping at h_gas, rather than switching
     to a Jeans-length estimate, is the resolution-robust choice across
     this model's wide production mass range. norm_grad_rho == 0 returns
     h_gas directly rather than dividing by 0: the final return below
     already zeroes the whole expression out for a star with no gas
     neighbours (rho_gas == 0 there), so no separate rho_gas guard is
     needed. */
  const float h_gas = sp->h * kernel_gamma;
  const float sobolev_length =
      norm_grad_rho > 0.0f ? fminf(rho_gas / norm_grad_rho, h_gas) : h_gas;
  const float length_gas = h_gas + sobolev_length;
  return length_gas * rho_gas;
}

/**
 * Compute a metallicity-scaled dust opacity around a star, in physical
 * internal units, from its cgs value at solar metallicity.
 *
 * @param sp The #spart.
 * @param us Unit system.
 * @param kappa0_cgs Opacity at solar metallicity, in cm^2/g.
 * @return Gas opacity around the star, in internal units.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_physical_opacity(const struct spart *sp,
                               const struct unit_system *us,
                               const float kappa0_cgs) {
  const float Z_gas = sp->feedback_data.Z_star;
  const float Z_sun = 0.02;
  const float value = kappa0_cgs *
                      units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
                      units_cgs_conversion_factor(us, UNIT_CONV_AREA);
  return value * Z_gas / Z_sun;
}

/**
 * Compute the physical optical depth around a star for a given opacity.
 *
 * @param sp The #spart.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param kappa0_cgs Opacity at solar metallicity, in cm^2/g.
 * @return Gas optical depth around the star.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_physical_optical_depth(const struct spart *sp,
                                     const struct unit_system *us,
                                     const struct cosmology *cosmo,
                                     const float kappa0_cgs) {
  const float Sigma_gas_c =
      radiation_get_comoving_gas_column_density_at_star(sp);
  const float Sigma_gas_p = Sigma_gas_c * cosmo->a2_inv;
  const float kappa = radiation_get_physical_opacity(sp, us, kappa0_cgs);
  return kappa * Sigma_gas_p;
}

/**
 * Compute the physical radiation pressure emitted by the star.
 *
 * LEBRON momentum coupling (Hopkins, Quataert & Murray 2012, MNRAS 421,
 * 3488, Sec 2.1; Hopkins et al. 2014, MNRAS 445, 581, App A): dot_p =
 * (1-exp(-tau_NUV)) * (1+tau_IR) * L_bol/c -- fraction of the non-ionizing
 * continuum absorbed before dust reprocessing (was assumed always 1), times
 * the IR-trapping boost, sharing one Sobolev column. kappa_NUV=1800
 * cm^2/g*(Z/Zsun) is one flux-mean opacity standing in for the whole
 * non-ionizing continuum (912A-3um) -- no band-by-band transport here, so
 * no per-band split either; 1800 is where the 2012 (single flux-mean over
 * that range) and 2020 (NUV sub-band) papers agree. Both are population-
 * (STARBURST99), not single-star-, calibrated, same caveat kappa_IR
 * already carries.
 *
 * @param sp The #spart.
 * @param Delta_t The current #spart timestep.
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @return Radiation pressure emittied by the star.
 */
__attribute__((always_inline)) INLINE float
radiation_get_star_physical_radiation_pressure(
    const struct spart *sp, const float Delta_t,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo) {

  const float tau_IR =
      radiation_get_physical_optical_depth(sp, us, cosmo, 10.0f);
  const float tau_NUV =
      radiation_get_physical_optical_depth(sp, us, cosmo, 1800.0f);
  const double L_bol = sp->feedback_data.radiation.L_bol;
  const double c = phys_const->const_speed_light_c;

  /* Double, not float: under -ffast-math the compiler may reassociate this
     product (e.g. L_bol*(1+tau_IR) before dividing by c), which overflows
     float32 at reachable extreme inputs (L_bol~1e38, tau_IR~19) even
     though the true result does not; negligible cost at one call per
     star feedback event. */
  const double p_rad = (double)Delta_t * L_bol / c *
                       (1.0 - exp(-(double)tau_NUV)) * (1.0 + (double)tau_IR);
  return (float)p_rad;
}

/**
 * Comoving gas column density AT a receiving gas particle's own location:
 * the RECEIVER-side generalization of
 * #radiation_get_comoving_gas_column_density_at_star, for Imladris (Smith
 * 2026, arXiv 2604.00100) Eq. 39's dust extinction of the LW/FUV bands.
 *
 * Phase-1 simplification, not yet the full generalization: uses the
 * kernel-radius cap unconditionally instead of a resolved density
 * gradient (`p->rho / |grad_rho|`, capped at `h_gas`, as the star-side
 * function does). #radiation_get_comoving_gas_column_density_at_star's
 * own norm_grad_rho == 0 branch already falls back to exactly this
 * (`length_gas = 2*h_gas`) when no gradient is resolved; this function
 * takes that fallback unconditionally rather than adding a new gas-gas
 * density-loop pass (a real, ~30-call-site touch of
 * src/runner_doiact_functions_hydro.h, per src/chemistry/GEAR/chemistry_
 * iact.h's own runner_iact_chemistry precedent) just to compute a
 * gradient that phase 2's Yukawa propagation will need to add anyway
 * (for its own harmonic-mean kappa_ij, .claude/dev/design-lw-fuv-
 * injection.md). Deferred there rather than duplicated here.
 *
 * @param p The #part.
 * @return Comoving gas column density at the particle's own location.
 */
__attribute__((always_inline)) INLINE float
radiation_get_comoving_gas_column_density_at_part(const struct part *p) {
  const float h_gas = p->h * kernel_gamma;
  return 2.0f * h_gas * p->rho;
}

/**
 * Dust-to-gas ratio relative to the Milky Way, exactly matching Grackle's
 * own internal convention (cool1d_multi_g.F: dust2gas(i) = fgr *
 * metallicity(i), metallicity(i) = metal(i,j,k)/d(i,j,k)/z_solar, when
 * chemistry_data.use_dust_density_field=0, the default) rather than an
 * independently-computed fit (e.g. Remy-Ruyer et al. 2014): our own
 * extinction's assumed dust abundance must track whatever Grackle's own
 * dust_chemistry=1-coupled channels (PE heating, H2-formation-on-dust,
 * dust recombination cooling) assume for the SAME gas, or the two would
 * disagree about how much dust is actually present. D(Z)/D(Z_sun) =
 * (fgr*Z/z_solar)/(fgr*1) = Z/z_solar: fgr (local_dust_to_gas_ratio)
 * cancels out of this ratio regardless of its configured value, so it is
 * not read here; only Grackle's own z_solar
 * (#RADIATION_GRACKLE_SOLAR_METAL_FRACTION) matters. Pure linear scaling
 * (Grackle applies no broken power law).
 *
 * @param Z Gas metal mass fraction.
 * @return Dust-to-gas ratio relative to the Milky Way.
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_to_gas_ratio_relative_to_MW(float Z) {
  return max(Z, 0.f) / RADIATION_GRACKLE_SOLAR_METAL_FRACTION;
}

/**
 * Band-specific dust extinction factor for the receiver-side LW/FUV
 * attenuation (Imladris Eq. 39): exp(-kappa_eff * Sigma_gas_p), with
 * kappa_eff = sigma_d_band * D(Z) / (mu_H * m_H) (see
 * .claude/dev/design-lw-fuv-injection.md's "Reconciling Imladris's
 * extinction with our own Yukawa propagation" for the unit-fix
 * derivation: sigma_d_band is a per-hydrogen-nucleon cross-section,
 * cm^2, not a mass opacity, so dividing by mu_H*m_H is required, not
 * optional). mu_H = 1.4 (mean mass per H nucleon, He folded in): the
 * design doc's own sanity check (sigma_d/2.3e-24 ~ 390-650 cm^2/g)
 * implicitly commits to this value, since 1.4*RADIATION_HYDROGEN_MASS_CGS
 * = 2.34e-24 g.
 *
 * @param us Unit system.
 * @param Z Gas metal mass fraction.
 * @param sigma_d_band_cgs Band-specific dust cross-section per hydrogen
 * nucleon, cm^2 (#RADIATION_SIGMA_D_FUV_CGS or #RADIATION_SIGMA_D_LW_CGS).
 * @param Sigma_gas_p Physical gas column density, internal units.
 * @return Dust extinction factor, in (0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_get_dust_extinction_factor(const struct unit_system *us, float Z,
                                     float sigma_d_band_cgs,
                                     float Sigma_gas_p) {

  const float D_relative = radiation_get_dust_to_gas_ratio_relative_to_MW(Z);
  const float kappa_eff_cgs = sigma_d_band_cgs * D_relative /
                              (RADIATION_MU_H * RADIATION_HYDROGEN_MASS_CGS);
  const float kappa_eff = kappa_eff_cgs *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS) /
                          units_cgs_conversion_factor(us, UNIT_CONV_AREA);
  const float tau = kappa_eff * Sigma_gas_p;
  return expf(-tau);
}

/**
 * @brief Receiver-side LW/FUV dust extinction factors for a gas particle
 * (Imladris Eq. 39), one per band.
 *
 * @param us Unit system.
 * @param cosmo The current cosmological model.
 * @param p The receiving #part.
 * @param Z The receiving particle's own metal mass fraction.
 * @param extinction_FUV (return) FUV-band (6-11.2 eV) extinction factor.
 * @param extinction_LW (return) Lyman-Werner-band (11.2-13.6 eV)
 * extinction factor.
 */
__attribute__((always_inline)) INLINE void
radiation_get_part_LW_FUV_extinction_factors(const struct unit_system *us,
                                             const struct cosmology *cosmo,
                                             const struct part *p, float Z,
                                             float *extinction_FUV,
                                             float *extinction_LW) {

  const float Sigma_gas_c =
      radiation_get_comoving_gas_column_density_at_part(p);
  const float Sigma_gas_p = Sigma_gas_c * cosmo->a2_inv;

  *extinction_FUV = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_FUV_CGS, Sigma_gas_p);
  *extinction_LW = radiation_get_dust_extinction_factor(
      us, Z, RADIATION_SIGMA_D_LW_CGS, Sigma_gas_p);
}
