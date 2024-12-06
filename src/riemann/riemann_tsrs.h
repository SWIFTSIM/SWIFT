/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_RIEMANN_TSRS_H
#define SWIFT_RIEMANN_TSRS_H

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "error.h"
#include "minmax.h"
#include "riemann_checks.h"
#include "riemann_common.h"
#include "riemann_vacuum.h"

#ifndef EOS_IDEAL_GAS
#error \
    "The TSRS Riemann solver currently only supports an ideal gas equation of state. Either select this equation of state, or try using another Riemann solver!"
#endif

__attribute__((always_inline)) INLINE static float riemann_g(float rho_state,
                                                             float P_state,
                                                             float P_guess) {
  const float A = hydro_two_over_gamma_plus_one / rho_state;
  const float B = hydro_gamma_minus_one_over_gamma_plus_one * P_state;
  return sqrtf(A / (P_guess + B));
}

/**
 * @brief Solve the Riemann problem using the Two Shock Riemann Solver
 *
 * By assuming 2 shock waves, we can get a direct (but still not exact) estimate
 * for the pressure in the star region, eliminating the iterative procedure.
 *
 * Like the TRRS, this implementations only assumes the two-shock approximation
 * to determine the estimate for P_star. The other variables in the star region
 * and the sampling of the solution happens identically to the exact Riemann
 * solver. This makes the TSRS more accurate according to Ivkovic (2023).
 *
 * We still need an initial estimate for P_star, from which we then compute a
 * more accurate estimate assuming two shock waves. Many choices are possible,
 * and here we opt for the value provided by the Primitive Value Riemann Solver
 * (Toro 2009).
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param Whalf Empty state vector in which the result will be stored
 * @param n_unit Normal vector of the interface
 * @return The pressure in the intermediate region (P_star)
 */
__attribute__((always_inline)) INLINE static float riemann_solver_solve(
    const float* WL, const float* WR, float* Whalf, const float* n_unit) {
  float aL, aR;
  float P0, gL, gR;
  float vL, vR;
  float P_star;

  /* calculate the velocities along the interface normal */
  vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* calculate the sound speeds */
  aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  if (riemann_is_vacuum(WL, WR, vL, vR, aL, aR)) {
    return riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
  }

  /* calculate the initial pressure */
  P0 = fmaxf(0.f, 0.5f * (WL[4] + WR[4]) -
                      0.125f * (vR - vL) * (WL[0] + WR[0]) * (aL + aR));

  /* Calculate P_star using the initial estimate P0 and the two shock
   * approximation */
  gL = riemann_g(WL[0], WL[4], P0);
  gR = riemann_g(WR[0], WR[4], P0);
  P_star = (gL * WL[4] + gR * WR[4] + vL - vR) / (gL + gR);

  /* sample the solution */
  riemann_sample(WL, WR, n_unit, vL, vR, aL, aR, P_star, Whalf);

  return P_star;
}

__attribute__((always_inline)) INLINE static float riemann_solve_for_flux(
    const float* Wi, const float* Wj, const float* n_unit, const float* vij,
    float* totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(Wi, Wj, n_unit, vij);
#endif

  float Whalf[5];
  const float P_star = riemann_solver_solve(Wi, Wj, Whalf, n_unit);
  riemann_flux_from_half_state(Whalf, vij, n_unit, totflux);

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(Wi, Wj, n_unit, vij, totflux);
#endif

  return P_star;
}

__attribute__((always_inline)) INLINE static void
riemann_solve_for_middle_state_flux(const float* Wi, const float* Wj,
                                    const float* n_unit, const float* vij,
                                    float* totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(Wi, Wj, n_unit, vij);
#endif

  if (Wi[0] == 0.0f || Wj[0] == 0.0f) {
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  /* calculate the velocities along the interface normal */
  const float vL = Wi[1] * n_unit[0] + Wi[2] * n_unit[1] + Wi[3] * n_unit[2];
  const float vR = Wj[1] * n_unit[0] + Wj[2] * n_unit[1] + Wj[3] * n_unit[2];

  /* calculate the sound speeds */
  const float aL = sqrtf(hydro_gamma * Wi[4] / Wi[0]);
  const float aR = sqrtf(hydro_gamma * Wj[4] / Wj[0]);

  if (riemann_is_vacuum(Wi, Wj, vL, vR, aL, aR)) {
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  /* calculate the initial pressure */
  float P0 = fmaxf(0.f, 0.5f * (Wi[4] + Wj[4]) -
                            0.125f * (vR - vL) * (Wi[0] + Wj[0]) * (aL + aR));

  /* Calculate P_star using the initial estimate P0 and the two shock
   * approximation */
  float gL = riemann_g(Wi[0], Wi[4], P0);
  float gR = riemann_g(Wj[0], Wj[4], P0);
  float P_star = (gL * Wi[4] + gR * Wj[4] + vL - vR) / (gL + gR);
  float v_star =
      0.5 * (vL + vR + (P_star - Wj[4]) * gR - (P_star - Wi[4]) * gL);

  totflux[0] = 0.0f;
  totflux[1] = P_star * n_unit[0];
  totflux[2] = P_star * n_unit[1];
  totflux[3] = P_star * n_unit[2];
  const float vface =
      vij[0] * n_unit[0] + vij[1] * n_unit[1] + vij[2] * n_unit[2];
  totflux[4] = P_star * (v_star + vface);

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(Wi, Wj, n_unit, vij, totflux);
#endif
}

#endif /* SWIFT_RIEMANN_TSRS_H */
