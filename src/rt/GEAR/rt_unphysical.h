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
 * @param U the density state
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_density(float U[4]){

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (U[0] < 0.f) message("Got unphysical energy %.3g | %.3g %.3g %.3g", U[0], U[1], U[2], U[3]);
#endif
  if (U[0] <= 0.f){
    U[0] = 0.f;
    U[1] = 0.f;
    U[2] = 0.f;
    U[3] = 0.f;
  }
}

/**
 * @brief check for and correct if needed unphysical 
 * values for a photon conserved state
 *
 * @param energy pointer to the photon energy
 * @param flux pointer to photon fluxes (3 dimensional)
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_conserved(float* energy, float* flux){

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (*energy < 0.f) message("Got unphysical energy %.3g | %.3g %.3g %.3g", *energy, flux[0], flux[1], flux[2]);
#endif
  if (*energy <= 0.f) {
    *energy = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
  }
}


#endif /* SWIFT_RT_UNPHYSICAL_GEAR_H */
