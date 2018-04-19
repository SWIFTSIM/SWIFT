/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

#ifndef SWIFT_RIEMANN_CHECKS_H
#define SWIFT_RIEMANN_CHECKS_H

#include "error.h"

#ifdef SWIFT_DEBUG_CHECKS

/**
 * @brief Check if the given state vector is physically valid.
 *
 * A valid state vector has 5 finite elements (not NaN), and a positive density
 * and pressure (first and fifth element).
 *
 * @param W State vector.
 * @return 0 if the state vector is valid, 1 otherwise.
 */
__attribute__((always_inline)) INLINE static int riemann_check_state(
    const float *W) {

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
__attribute__((always_inline)) INLINE static int riemann_check_vector(
    const float *x) {

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
 * @param WL Left state vector.
 * @param WR Right state vector.
 * @param n Surface normal vector.
 * @param vij Surface velocity vector.
 */
__attribute__((always_inline)) INLINE static void riemann_check_input(
    const float *WL, const float *WR, const float *n, const float *vij) {

  int errorFlag = 0;

  if (riemann_check_state(WL)) {
    message("Invalid left state!");
    errorFlag = 1;
  }

  if (riemann_check_state(WR)) {
    message("Invalid right state!");
    errorFlag = 1;
  }

  if (riemann_check_vector(n)) {
    message("Invalid surface normal vector!");
    errorFlag = 1;
  }

  if (riemann_check_vector(vij)) {
    message("Invalid face velocity vector!");
    errorFlag = 1;
  }

  if (errorFlag) {
    message("WL: %g %g %g %g %g", WL[0], WL[1], WL[2], WL[3], WL[4]);
    message("WR: %g %g %g %g %g", WR[0], WR[1], WR[2], WR[3], WR[4]);
    message("n: %g %g %g", n[0], n[1], n[2]);
    message("vij: %g %g %g", vij[0], vij[1], vij[2]);
    error("Invalid Riemann solver input!");
  }
}

/**
 * @brief Check that the given output from a Riemann solver is physically valid.
 *
 * Physically valid output consists of 5 finite (not NaN) flux components.
 * If no valid output is provided, an error is thrown.
 *
 * @param WL Left state vector.
 * @param WR Right state vector.
 * @param n Surface normal vector.
 * @param vij Surface velocity vector.
 * @param totflux Riemann solver flux result.
 */
__attribute__((always_inline)) INLINE static void riemann_check_output(
    const float *WL, const float *WR, const float *n, const float *vij,
    const float *totflux) {

  int errorFlag = 0;

  /* check that all the fluxes are finite */
  if (totflux[0] != totflux[0]) {
    message("NaN mass flux!");
    errorFlag = 1;
  }
  if (totflux[1] != totflux[1]) {
    message("NaN x momentum flux!");
    errorFlag = 1;
  }
  if (totflux[2] != totflux[2]) {
    message("NaN y momentum flux!");
    errorFlag = 1;
  }
  if (totflux[3] != totflux[3]) {
    message("NaN z momentum flux!");
    errorFlag = 1;
  }
  if (totflux[4] != totflux[4]) {
    message("NaN energy flux!");
    errorFlag = 1;
  }

  if (errorFlag) {
    message("WL: %g %g %g %g %g", WL[0], WL[1], WL[2], WL[3], WL[4]);
    message("WR: %g %g %g %g %g", WR[0], WR[1], WR[2], WR[3], WR[4]);
    message("n: %g %g %g", n[0], n[1], n[2]);
    message("vij: %g %g %g", vij[0], vij[1], vij[2]);
    message("totflux: %g %g %g %g %g", totflux[0], totflux[1], totflux[2],
            totflux[3], totflux[4]);
    error("Invalid Riemann solver output!");
  }
}

#endif

#endif /* SWIFT_RIEMANN_CHECKS_H */
