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
#ifndef SWIFT_RT_UNPHYSICAL_GEAR_H
#define SWIFT_RT_UNPHYSICAL_GEAR_H

#include "minmax.h"

/**
 * @file src/rt/GEAR/rt_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

/**
 * @brief Check that the float sum of squares of a vector is exact to float
 * rounding: no square overflows (components up to 2^63), and no square that
 * matters underflows (below 2^-63 only for components under 2^-13 of the
 * largest). Radiation fluxes leave that range in cgs and in galactic units.
 *
 * @param v the vector
 */
__attribute__((always_inline)) INLINE static int rt_float_squares_in_range(
    const float v[3]) {

  const float a[3] = {fabsf(v[0]), fabsf(v[1]), fabsf(v[2])};
  const float amax = max3(a[0], a[1], a[2]);
  if (amax == 0.f) return 1;
  if (amax > 0x1p63f || amax < 0x1p-63f) return 0;
  const float negligible = amax * 0x1p-13f;
  for (int k = 0; k < 3; k++)
    if (a[k] < 0x1p-63f && a[k] > negligible) return 0;
  return 1;
}

/**
 * @brief Squared norm of a flux vector. Computed in double when its float
 * squares would leave the float range, in float otherwise.
 *
 * @param F the flux vector
 */
__attribute__((always_inline)) INLINE static double rt_flux_norm2(
    const float F[3]) {

  if (rt_float_squares_in_range(F))
    return F[0] * F[0] + F[1] * F[1] + F[2] * F[2];
  return (double)F[0] * F[0] + (double)F[1] * F[1] + (double)F[2] * F[2];
}

/**
 * @brief Norm of a flux vector. Computed in double when its float squares
 * would leave the float range, in float otherwise. A norm below FLT_MIN is
 * returned as 0, since its inverse would overflow.
 *
 * @param F the flux vector
 */
__attribute__((always_inline)) INLINE static float rt_flux_norm(
    const float F[3]) {

  if (rt_float_squares_in_range(F))
    return sqrtf(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
  const double norm =
      sqrt((double)F[0] * F[0] + (double)F[1] * F[1] + (double)F[2] * F[2]);
  return (norm < FLT_MIN) ? 0.f : (float)norm;
}

/**
 * @brief check for and correct if needed unphysical
 * values for a radiation state.
 *
 * @param energy_density pointer to the radiation energy density
 * @param flux pointer to radiation flux (3 dimensional)
 * @param e_old energy density before change to check. Set = 0 if not available
 * @param callloc integer indentifier where this function was called from
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_state(
    float *energy_density, float *flux, const float e_old, int callloc) {

  /* Check for negative energies */
  /* Note to self for printouts: Maximal allowable F = E * c.
   * In some cases, e.g. while cooling, we don't modify the fluxes,
   * so you can get an estimate of what the photon energy used to be
   * by dividing the printed out fluxes by the speed of light in
   * code units */
#ifdef SWIFT_DEBUG_CHECKS
  float ratio = -1.f;
  char print = 0;
  if (e_old == 0.) {
    if (*energy_density < -1.e-4f) print = 1;
  } else {
    if (fabsf(*energy_density) > 1.e-30) {
      if (*energy_density < -1.e-4f * fabsf(e_old)) print = 1;
      ratio = fabsf(*energy_density / e_old);
    }
  }
  /* callloc = 1 is gradient extrapolation. Don't print out those. */
  if (callloc == 1) print = 0;
  if (print)
    message("Fixing unphysical energy case %d | %.6e | %.6e %.6e %.6e | %.6e",
            callloc, *energy_density, flux[0], flux[1], flux[2], ratio);
#endif
  if (isinf(*energy_density) || isnan(*energy_density))
    error("Got inf/nan radiation energy case %d | %.6e | %.6e %.6e %.6e",
          callloc, *energy_density, flux[0], flux[1], flux[2]);

  if (*energy_density <= 0.f) {
    *energy_density = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
    return;
  }

  /* Check for too high fluxes */
  const double flux_norm = sqrt(rt_flux_norm2(flux));
  const double flux_max = rt_params.reduced_speed_of_light * *energy_density;
  if (flux_norm > flux_max) {
    const double correct = flux_max / flux_norm;
    flux[0] *= correct;
    flux[1] *= correct;
    flux[2] *= correct;
  }
}

/**
 * @brief Do additional checks after reading in initial conditions, and exit on
 * error.
 *
 * @param p particle we're checking
 * @param group current photon group we're checking
 * @param energy_density pointer to the radiation energy density
 * @param flux pointer to radiation flux (3 dimensional)
 * @param c the speed of light (in internal units). NOT the reduced speed of
 * light.
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_state_ICs(
    const struct part *restrict p, int group, float *energy_density,
    float *flux, const double c) {

  /* Nothing to do here. The other unphysical check will catch other problems.
   */
  if (*energy_density == 0.f) return;

  /* Check for negative energies */
  if (*energy_density < 0.f)
    error(
        "Found particle with negative energy density after reading in ICs: "
        "pid= %lld group=%d E=%.6g",
        p->id, group, *energy_density);
  if (*energy_density > FLT_MAX || isnan(*energy_density))
    error("Got inf/nan energy_density: %g", *energy_density);

  /* Check for too high fluxes */
  const double flux_norm = sqrt(rt_flux_norm2(flux));
  const float flux_max = c * *energy_density;
  if (flux_max > FLT_MAX || isnan(flux_max))
    error("Got inf/nan flux_max: %g", flux_max);
  if (flux_norm > FLT_MAX || isnan(flux_norm))
    error("Got inf/nan flux_norm: %g", flux_norm);
  if (flux_norm > flux_max * 1.0001) {
    error(
        "Found too high radiation flux for a particle: pid=%lld, group=%d, "
        "have=%.6g, max=%.6g",
        p->id, group, flux_norm, flux_max);
  }
}

/**
 * @brief check for and correct if needed unphysical
 * values for a flux in the sense of hyperbolic conservation laws
 *
 * @param flux hyperbolic flux: 4 components (photon energy +
 *        photon flux) x 3 dimensions each
 */
__attribute__((always_inline)) INLINE static void
rt_check_unphysical_hyperbolic_flux(float flux[4][3]) {

#ifdef SWIFT_DEBUG_CHECKS
  int nans = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (isnan(flux[i][j])) {
        nans += 1;
        break;
      }
    }
  }

  if (nans) {
    message(
        "Fixing unphysical hyperbolic flux:"
        " %.3e %.3e %.3e | %.3e %.3e %.3e |"
        " %.3e %.3e %.3e | %.3e %.3e %.3e",
        flux[0][0], flux[0][1], flux[0][2], flux[1][0], flux[1][1], flux[1][2],
        flux[2][0], flux[2][1], flux[2][2], flux[3][0], flux[3][1], flux[3][2]);
  }
#endif

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (isnan(flux[i][j])) {
        flux[i][j] = 0.f;
      }
    }
  }
}

/**
 * @brief check whether gas species mass fractions have physical
 * values and correct small errors if necessary.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void
rt_check_unphysical_mass_fractions(struct part *restrict p) {

  /* For now, catch either mass or rho being zero. At the moment, they are not
   * necessarily both zero. For example, an unphysical check may zero out both
   * mass and density when it becomes negative in a hydro step. Once that
   * happens, particles may gain mass through flux exchanges with other active
   * particles, while they themselves remain inactive. The density of such
   * inactive particles however remains zero until the particle is active
   * again. See issue #833. */

  if (hydro_get_mass(p) <= 0.f || p->rho <= 0.f) {
    /* Deal with unphysical situations and vacuum. */
    p->rt_data.tchem.mass_fraction_HI = RT_GEAR_TINY_MASS_FRACTION;
    p->rt_data.tchem.mass_fraction_HII = RT_GEAR_TINY_MASS_FRACTION;
    p->rt_data.tchem.mass_fraction_HeI = RT_GEAR_TINY_MASS_FRACTION;
    p->rt_data.tchem.mass_fraction_HeII = RT_GEAR_TINY_MASS_FRACTION;
    p->rt_data.tchem.mass_fraction_HeIII = RT_GEAR_TINY_MASS_FRACTION;
    return;
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (p->rt_data.tchem.mass_fraction_HI < -1e4)
    message("WARNING: Got negative HI mass fraction?");
  if (p->rt_data.tchem.mass_fraction_HII < -1e4)
    message("WARNING: Got negative HII mass fraction?");
  if (p->rt_data.tchem.mass_fraction_HeI < -1e4)
    message("WARNING: Got negative HeI mass fraction?");
  if (p->rt_data.tchem.mass_fraction_HeII < -1e4)
    message("WARNING: Got negative HeII mass fraction?");
  if (p->rt_data.tchem.mass_fraction_HeIII < -1e4)
    message("WARNING: Got negative HeIII mass fraction?");
#endif

  /* TODO: this should be a for loop with mass fractions being enums. */
  p->rt_data.tchem.mass_fraction_HI =
      max(p->rt_data.tchem.mass_fraction_HI, RT_GEAR_TINY_MASS_FRACTION);
  p->rt_data.tchem.mass_fraction_HII =
      max(p->rt_data.tchem.mass_fraction_HII, RT_GEAR_TINY_MASS_FRACTION);
  p->rt_data.tchem.mass_fraction_HeI =
      max(p->rt_data.tchem.mass_fraction_HeI, RT_GEAR_TINY_MASS_FRACTION);
  p->rt_data.tchem.mass_fraction_HeII =
      max(p->rt_data.tchem.mass_fraction_HeII, RT_GEAR_TINY_MASS_FRACTION);
  p->rt_data.tchem.mass_fraction_HeIII =
      max(p->rt_data.tchem.mass_fraction_HeIII, RT_GEAR_TINY_MASS_FRACTION);

  const float XHI = p->rt_data.tchem.mass_fraction_HI;
  const float XHII = p->rt_data.tchem.mass_fraction_HII;
  const float XHeI = p->rt_data.tchem.mass_fraction_HeI;
  const float XHeII = p->rt_data.tchem.mass_fraction_HeII;
  const float XHeIII = p->rt_data.tchem.mass_fraction_HeIII;

  const float Xtot = XHI + XHII + XHeI + XHeII + XHeIII;

  /* Make sure we sum up to 1. TODO: Assuming we have no metals. */
  if (fabsf(Xtot - 1.f) > 1e-3)
    error("Got total mass fraction of gas = %.6g", Xtot);
}

#endif /* SWIFT_RT_UNPHYSICAL_GEAR_H */
