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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_RIEMANN_CHECKS_H
#define SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_RIEMANN_CHECKS_H

#include "error.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/parabolic/chemistry_riemann_checks.h
 * @brief File containing functions to check that in/output of the Riemann
 * solver are meaningful, e.g. no NaNs or negative densities.
 * */

/**
 * @brief Check if the given state vector is physically valid.
 *
 * A valid state vector has 5 finite elements (not NaN), and a positive density
 * and pressure (first and fifth element).
 *
 * @param W Hydro state vector.
 * @param U Diffusion state vector.
 * @return 0 if the state vector is valid, 1 otherwise.
 */
__attribute__((always_inline)) INLINE static int chemistry_riemann_check_state(
    const float W[5], const double U) {

  int errorFlag = 0;

  /* check the density: should be finite and positive */
  if (W[0] != W[0]) {
    message("NaN density!");
    errorFlag = 1;
  }
  if (W[0] < 0.f) {
    message("Negative density!");
    errorFlag = 1;
  }

  /* check the velocities: should be finite */
  if (W[1] != W[1]) {
    message("NaN x velocity!");
    errorFlag = 1;
  }
  if (W[2] != W[2]) {
    message("NaN y velocity!");
    errorFlag = 1;
  }
  if (W[3] != W[3]) {
    message("NaN z velocity!");
    errorFlag = 1;
  }

  /* check the pressure: should be positive and finite */
  if (W[4] != W[4]) {
    message("NaN pressure!");
    errorFlag = 1;
  }
  if (W[4] < 0.f) {
    message("Negative pressure!");
    errorFlag = 1;
  }

  /* check the metal density: should be finite and positive */
  if (U != U) {
    message("NaN metal density!");
    errorFlag = 1;
  }
  if (U < 0.0) {
    message("Negative metal density!");
    errorFlag = 1;
  }

  return errorFlag;
}

/**
 * @brief Check that the given vector is physically valid.
 *
 * A valid vector has 3 finite elements (not NaN).
 *
 * @param x Vector to check.
 * @return 0 if the vector is valid, 1 otherwise.
 */
__attribute__((always_inline)) INLINE static int chemistry_riemann_check_vector(
    const float x[3]) {

  int errorFlag = 0;

  /* check that all components are finite */
  if (x[0] != x[0]) {
    message("NaN x component!");
    errorFlag = 1;
  }
  if (x[1] != x[1]) {
    message("NaN y component!");
    errorFlag = 1;
  }
  if (x[2] != x[2]) {
    message("NaN z component!");
    errorFlag = 1;
  }

  return errorFlag;
}

/**
 * @brief Check that the given input for a Riemann solver is physically valid.
 *
 * Physically valid input consists of a valid left and right state, and valid
 * surface normal and surface velocity vectors.
 * If no valid input is provided, an error is thrown.
 *
 * @param WL Left hydro state vector.
 * @param WR Right hydro state vector.
 * @param UL Left diffusion state vector (metal density).
 * @param UR Right diffusion state vector (metal density).
 * @param n Surface normal vector.
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_check_input(
    const float WL[5], const float WR[5], const double UL, const double UR,
    const float n[3]) {

  int errorFlag = 0;

  if (chemistry_riemann_check_state(WL, UL)) {
    message("Invalid left state!");
    errorFlag = 1;
  }

  if (chemistry_riemann_check_state(WR, UR)) {
    message("Invalid right state!");
    errorFlag = 1;
  }

  if (chemistry_riemann_check_vector(n)) {
    message("Invalid surface normal vector!");
    errorFlag = 1;
  }

  if (errorFlag) {
    message("WL: %g %g %g %g %g, UL = %e", WL[0], WL[1], WL[2], WL[3], WL[4],
            UL);
    message("WR: %g %g %g %g %g, UR = %e", WR[0], WR[1], WR[2], WR[3], WR[4],
            UR);
    message("n: %g %g %g", n[0], n[1], n[2]);
    error("Invalid Riemann solver input!");
  }
}

/**
 * @brief Check that the given output from a Riemann solver is physically valid.
 *
 * Physically valid output consists of 5 finite (not NaN) flux components.
 * If no valid output is provided, an error is thrown.
 *
 * @param WL Left hydro state vector.
 * @param WR Right hydro state vector.
 * @param UL Left diffusion state vector (metal density).
 * @param UR Right diffusion state vector (metal density).
 * @param n Surface normal vector.
 * @param totflux Riemann solver flux result.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_check_output(const float WL[5], const float WR[5],
                               const double UL, const double UR,
                               const float n[3], const double* totflux) {

  int errorFlag = 0;

  /* Check that metal mass flux is finite */
  if (totflux[0] != totflux[0]) {
    message("NaN mass flux!");
    errorFlag = 1;
  }

  if (errorFlag) {
    message("WL: %g %g %g %g %g, UL = %e", WL[0], WL[1], WL[2], WL[3], WL[4],
            UL);
    message("WR: %g %g %g %g %g, UR = %e", WR[0], WR[1], WR[2], WR[3], WR[4],
            UR);
    message("n: %g %g %g", n[0], n[1], n[2]);
    /* message("vij: %g %g %g", vij[0], vij[1], vij[2]); */
    message("totflux: %g %g %g %g %g", totflux[0], totflux[1], totflux[2],
            totflux[3], totflux[4]);
    error("Invalid Riemann solver output!");
  }
}
#endif /* SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_RIEMANN_CHECKS_H */
