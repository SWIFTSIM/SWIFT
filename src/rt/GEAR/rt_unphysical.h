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

/**
 * @file src/rt/GEAR/rt_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

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
    float* energy_density, float* flux, const float e_old, int callloc) {

  /* Check for negative energies */
  /* Note to self for printouts: Maximal allowable F = E * c.
   * In some cases, e.g. while cooling, we don't modify the fluxes,
   * so you can get an estimate of what the photon energy used to be
   * by dividing the printed out fluxes by the speed of light in
   * code units */
#ifdef SWIFT_DEBUG_CHECKS
  float ratio = 2.;
  if (e_old != 0.f) ratio = fabsf(*energy_density / e_old);
  /* callloc = 1 is gradient extrapolation. Don't print out those. */
  if (*energy_density < -1e-2f && fabsf(ratio - 1.f) > 1.e-3f && callloc != 1)
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
  const float flux2 = flux[0] * flux[0] + flux[1] * flux[1] + flux[2] * flux[2];
  const float flux_norm = sqrtf(flux2);
  const float flux_max = rt_params.reduced_speed_of_light * *energy_density;
  if (flux_norm > flux_max) {
    const float correct = flux_max / flux_norm;
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
    const struct part* restrict p, int group, float* energy_density,
    float* flux, const double c) {

  /* Nothing to do here. The other unphysical check will catch other problems.
   */
  if (*energy_density == 0.f) return;

  /* Check for negative energies */
  if (*energy_density < 0.f)
    error(
        "Found particle with negative energy density after reading in ICs: "
        "pid= %lld group=%d E=%.6g",
        p->id, group, *energy_density);

  /* Check for too high fluxes */
  const float flux2 = flux[0] * flux[0] + flux[1] * flux[1] + flux[2] * flux[2];
  const float flux_norm = sqrtf(flux2);
  const float flux_max = c * *energy_density;
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
rt_check_unphysical_mass_fractions(struct part* restrict p) {

/* GRACKLE doesn't really like exact zeroes, so use something
 * comparatively small instead. */
#define RT_GEAR_TINY_MASS_FRACTION 1.e-6

  if (p->rt_data.tchem.mass_fraction_HI < 0.f) {
    if (p->rt_data.tchem.mass_fraction_HI < -1e4)
      message("WARNING: Got negative HI mass fraction?");
    p->rt_data.tchem.mass_fraction_HI = RT_GEAR_TINY_MASS_FRACTION;
  }
  if (p->rt_data.tchem.mass_fraction_HII < 0.f) {
    if (p->rt_data.tchem.mass_fraction_HII < -1e4)
      message("WARNING: Got negative HII mass fraction?");
    p->rt_data.tchem.mass_fraction_HII = RT_GEAR_TINY_MASS_FRACTION;
  }
  if (p->rt_data.tchem.mass_fraction_HeI < 0.f) {
    if (p->rt_data.tchem.mass_fraction_HeI < -1e4)
      message("WARNING: Got negative HeI mass fraction?");
    p->rt_data.tchem.mass_fraction_HeI = RT_GEAR_TINY_MASS_FRACTION;
  }
  if (p->rt_data.tchem.mass_fraction_HeII < 0.f) {
    if (p->rt_data.tchem.mass_fraction_HeII < -1e4)
      message("WARNING: Got negative HeII mass fraction?");
    p->rt_data.tchem.mass_fraction_HeII = RT_GEAR_TINY_MASS_FRACTION;
  }
  if (p->rt_data.tchem.mass_fraction_HeIII < 0.f) {
    if (p->rt_data.tchem.mass_fraction_HeIII < -1e4)
      message("WARNING: Got negative HeIII mass fraction?");
    p->rt_data.tchem.mass_fraction_HeIII = RT_GEAR_TINY_MASS_FRACTION;
  }

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
