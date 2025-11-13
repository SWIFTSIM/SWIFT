/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UNPHYSICAL_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UNPHYSICAL_H

#include "error.h"
#include "inline.h"
#include "part.h"

/**
 * @file src/chemistry/GEAR/chemistry_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

/**
 * @brief Check for, and correct, if needed, unphysical values for a diffusion
 * state (metal mass).
 *
 * @param metal_mass Pointer to the metal mass.
 * @param mZ_old Metal density before change to check. Set = 0 if not
 * available.
 * @param gas_mass #part mass.
 * @param callloc Integer indentifier where this function was called from
 *
 * 0: after density compuations (chemistry_end_density()),
 * 1: during gradients extrapolation (chemistry_gradients_predict()),
 * 2: after diffusion computations (chemistry_end_force()),
 * 3: during drift (chemistry_predict_extra()).
 *
 * @param element Integer identifier for the metal specie to correct.
 * @param id Particle id.
 */
__attribute__((always_inline)) INLINE static void
chemistry_check_unphysical_state(double* metal_mass, const double mZ_old,
                                 const double gas_mass, int callloc,
                                 const int element, const long long id) {

  if (isinf(*metal_mass) || isnan(*metal_mass))
    error("[%lld, %d] Got inf/nan metal density/mass diffusion case %d | %.6e ",
          id, element, callloc, *metal_mass);

  /* Fix negative masses */
  const double metal_mass_fraction = *metal_mass / gas_mass;
  const double Z_old = mZ_old / gas_mass;
  if (metal_mass_fraction < GEAR_NEGATIVE_METAL_MASS_FRACTION_TOLERANCE) {
    if (callloc == 1) {
      /* Do not extrapolate, use 0th order reconstruction. */
      *metal_mass =
          (Z_old >= GEAR_NEGATIVE_METAL_MASS_FRACTION_TOLERANCE) ? mZ_old : 0.0;
    } else {
      /* Note: Correcting metal masses afterwards can artificially create metal
         mass out of nothing. This mass creation might is never compensated and
         can lead to huge metal mass creation, bigger than gas mass. */
      error("[%lld, %d] Negative metal density/mass case %d | %e %e | %e %e",
            id, element, callloc, *metal_mass, metal_mass_fraction, mZ_old,
            Z_old);
      /* metal_mass = 0.0; */
    }
  }

  if (*metal_mass > gas_mass) {
    if (callloc == 1 && mZ_old <= gas_mass) {
      /* Do not extrapolate, use 0th order reconstruction. */
      *metal_mass = mZ_old;
    } else {
      *metal_mass /= 1.1 * mZ_old / gas_mass;
      warning(
          "[%lld, %d] Metal mass bigger than gas mass ! case %d | %e | %e | %e",
          id, element, callloc, *metal_mass, mZ_old, gas_mass);
    }
  }
}

/**
 * @brief Check for and correct if needed unphysical values for a diffusion 3D
 * flux.
 *
 * @param flux flux: 3 components diffusion flux
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

/**
 * @brief Check for, and correct if needed, unphysical total metal mass.
 *
 * @param p The #part to check.
 * @param callloc integer indentifier where this function was called from
 *
 * 0: after density compuations (chemistry_end_density()).
 * 2: after diffusion computations (chemistry_end_force()).
 * 3: during drift (chemistry_predict_extra()).
 *
 */
__attribute__((always_inline)) INLINE static void
chemistry_check_unphysical_total_metal_mass(struct part* restrict p,
                                            int callloc) {
  struct chemistry_part_data* chd = &p->chemistry_data;

  /* Verify that the total metal mass does not exceed the part's mass */
  const float gas_mass = hydro_get_mass(p);
  const double m_Z_tot = chemistry_get_total_metal_mass_fraction(p) * gas_mass;

  if (m_Z_tot > gas_mass) {
    /* Rescale the elements */
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      chd->metal_mass[i] /= 1e1 * m_Z_tot / gas_mass;
    }
    warning(
        "[%lld, %i] Total metal mass grew larger than the particle mass! "
        "Rescaling the element masses. m_Z_tot = %e, m = %e"
        " m_z_0 = %e, m_z_1 = %e, m_z_2 = %e, m_z_3 = %e, m_z_4 = %e, "
        "m_z_5 = %e, m_z_6 = %e, m_z_7 = %e, m_z_8 = %e, m_z_9 = %e",
        p->id, callloc, m_Z_tot, gas_mass, chd->metal_mass[0],
        chd->metal_mass[1], chd->metal_mass[2], chd->metal_mass[3],
        chd->metal_mass[4], chd->metal_mass[5], chd->metal_mass[6],
        chd->metal_mass[7], chd->metal_mass[8], chd->metal_mass[9]);
  }
}

/**
 * @brief Check for and correct if needed unphysical values for a flux in the
 * sense of hyperbolic conservation laws.
 *
 * @param flux Hyperbolic flux: 4 components (metal density +
 *        metal flux) x 3 dimensions each.
 */
__attribute__((always_inline)) INLINE static void
chemistry_check_unphysical_hyperbolic_flux(double  flux[4][3]) {

/* #ifdef SWIFT_DEBUG_CHECKS */
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
/* #endif */

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      if (isnan(flux[i][j])) {
        flux[i][j] = 0.f;
      }
    }
  }
}


#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UNPHYSICAL_H */
