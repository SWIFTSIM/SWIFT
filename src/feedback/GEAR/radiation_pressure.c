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
#include "radiation.h"

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

  /* A locally uniform density field (zero or near-zero gradient, e.g.
     unperturbed glass/grid ICs) makes the Sobolev length rho/|grad rho|
     undefined, or dominated by SPH summation noise rather than a resolved
     trend: fall back to just the kernel support radius in that case (see
     RADIATION_MIN_RELATIVE_DENSITY_GRADIENT). */
  const float h_gas = sp->h * kernel_gamma;
  const float grad_rho_floor =
      RADIATION_MIN_RELATIVE_DENSITY_GRADIENT * rho_gas / h_gas;
  const float sobolev_length =
      norm_grad_rho > grad_rho_floor ? rho_gas / norm_grad_rho : 0.0f;
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
  const float L_bol = sp->feedback_data.radiation.L_bol;
  const float c = phys_const->const_speed_light_c;

  return Delta_t * L_bol / c * (1.0f - expf(-tau_NUV)) * (1.0f + tau_IR);
}
