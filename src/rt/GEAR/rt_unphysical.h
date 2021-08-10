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
 * values for a photon density state
 *
 * @param energy pointer to the photon energy density
 * @param flux pointer to the photon flux density
 * @param c integer identifier where this function was called from
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_density(
    float* energy, float* flux, int c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (*energy < 0.f && fabs(*energy) > 1.e-1)
    message("Fixing unphysical energy case%d %.3g | %.3g %.3g %.3g", c, *energy,
            flux[0], flux[1], flux[2]);
#endif
  if (*energy <= 0.f) {
    *energy = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
    return;
  }

  /* const float flux2 = flux[0] * flux[0] + flux[1] * flux[1] + flux[2] *
   * flux[2]; */
  /* const float flux_norm = sqrtf(flux2); */
  /* const float flux_max = rt_params.reduced_speed_of_light * *energy; */
  /* if (flux_norm > flux_max) { */
  /*   const float correct = flux_max / flux_norm; */
  /*   flux[0] *= correct; */
  /*   flux[1] *= correct; */
  /*   flux[2] *= correct; */
  /* } */
}

/**
 * @brief check for and correct if needed unphysical
 * values for a photon conserved state
 *
 * @param energy pointer to the photon energy
 * @param flux pointer to photon fluxes (3 dimensional)
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_conserved(
    float* energy, float* flux) {

#ifdef SWIFT_DEBUG_CHECKS
  if (*energy < 0.f && fabs(*energy) > 1.e-1)
    message("Fixing unphysical energy %.3g | %.3g %.3g %.3g", *energy, flux[0],
            flux[1], flux[2]);
#endif
  if (*energy <= 0.f) {
    *energy = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
    return;
  }

  /* const float flux2 = flux[0] * flux[0] + flux[1] * flux[1] + flux[2] *
   * flux[2]; */
  /* const float flux_norm = sqrtf(flux2); */
  /* const float flux_max = rt_params.reduced_speed_of_light * *energy; */
  /* if (flux_norm > flux_max) { */
  /*   const float correct = flux_max / flux_norm; */
  /*   flux[0] *= correct; */
  /*   flux[1] *= correct; */
  /*   flux[2] *= correct; */
  /* } */
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

  int nans = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (flux[i][j] != flux[i][j]) {
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

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (flux[i][j] != flux[i][j]) {
        flux[i][j] = 0.f;
      }
    }
  }
}
#endif /* SWIFT_RT_UNPHYSICAL_GEAR_H */
