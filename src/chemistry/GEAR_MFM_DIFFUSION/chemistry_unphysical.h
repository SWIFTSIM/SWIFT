/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_UNPHYSICAL_H
#define SWIFT_CHEMiSTRY_GEAR_MFM_DIFFUSION_UNPHYSICAL_H

#include "error.h"
#include "inline.h"

/**
 * @file src/chemistry/GEAR/chemistry_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

/**
 * @brief check for and correct if needed unphysical
 * values for a diffusion state.
 *
 * @param metal_density pointer to the radiation energy density
 * @param n_old metal density before change to check. Set = 0 if not available
 * @param callloc integer indentifier where this function was called from
 */
__attribute__((always_inline)) INLINE static void chemistry_check_unphysical_state(
    double* metal_density, const double n_old, int callloc) {

  /* Check for negative metal densities */
  /* Note to self for printouts: Maximal allowable F = E * c.
   * In some cases, e.g. while cooling, we don't modify the fluxes,
   * so you can get an estimate of what the photon energy used to be
   * by dividing the printed out fluxes by the speed of light in
   * code units */
#ifdef SWIFT_DEBUG_CHECKS
  float ratio = -1.f;
  char print = 0;
  if (n_old == 0.0) {
    if (*metal_density < -1.e-20) print = 1;
  } else {
    /* TODO: understadn that and transpose to chemistry */
    /* if (fabs(*metal_density) > 1.e-30) { */
    /*   if (*energy_density < -1.e-20 * fabs(e_old)) print = 1; */
    /*   ratio = fabsf(*energy_density / e_old); */
    /* } */
  }
  /* callloc = 1 is gradient extrapolation. Don't print out those. */
  if (callloc == 1) print = 0;
  if (print)
    message("Fixing unphysical metal density case %d | %.6e | %.6e",
            callloc, *metal_density, ratio);
#endif
  if (isinf(*metal_density) || isnan(*metal_density))
    error("Got inf/nan metal density diffusion case %d | %.6e ",
          callloc, *metal_density);

  if (*metal_density < 0.0) {
    *metal_density = 0.0;
    return;
  }
}

/**
 * @brief check for and correct if needed unphysical
 * values for a flux in the sense of parabolic conservation laws
 *
 * @param flux  flux: 3 components
 */
__attribute__((always_inline)) INLINE static void
chemistry_check_unphysical_diffusion_flux(double flux[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  int nans = 0;
  for (int i = 0; i < 3; i++) {
    if (isnan(flux[i])) {
      nans += 1;
      break;
    }
  }

  if (nans) {
    message(
        "Fixing unphysical diffusion flux:"
        " %.3e %.3e %.3e",
        flux[0], flux[1], flux[2]);
  }
#endif

  for (int i = 0; i < 3; i++) {
    if (isnan(flux[i])) {
      flux[i] = 0.f;
    }
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_UNPHYSICAL_H */
