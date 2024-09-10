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
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_UNPHYSICAL_H

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
 * Check that the metal density does not exceed the particle mass. To do so, we
 * whould check all metallicities at the same time. We should also ensure that
 * the sum of the metallicities does not exceed the particle mass.
 *
 * @param metal_density pointer to the radiation energy density
 * @param n_old metal density before change to check. Set = 0 if not available
 * @param callloc integer indentifier where this function was called from
 */
__attribute__((always_inline)) INLINE static void
chemistry_check_unphysical_state(double* metal_mass, const double mZ_old, const double gas_mass,
                                 int callloc) {

  /* Check for negative metal densities/masses */
#ifdef SWIFT_DEBUG_CHECKS
  /* float ratio = -1.f; */
  char print = 0;
  if (mZ_old == 0.0) {
    if (*metal_mass < -1.e-20) print = 1;
  }
  /* callloc = 1 is gradient extrapolation. Don't print out those. */
  if (callloc == 1) print = 0;
  if (print)
    error("Fixing unphysical metal density/mass case %d | %.6e | %.6e", callloc,
            *metal_mass, mZ_old);
#endif
  if (isinf(*metal_mass) || isnan(*metal_mass))
    error("Got inf/nan metal density/mass diffusion case %d | %.6e ", callloc,
          *metal_mass);

  if (*metal_mass < 0.0) {
    *metal_mass = 0.0;
    return;
  }

  if (*metal_mass > gas_mass) {
    error("Metal mass bigger than gas mass ! case %d | %e | %e", callloc,
	  *metal_mass, mZ_old);
    if (mZ_old <= gas_mass) {
      *metal_mass = mZ_old;
    } else {
      *metal_mass = 0.0;
    }
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
