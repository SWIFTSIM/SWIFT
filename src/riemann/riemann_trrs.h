/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

#ifndef SWIFT_RIEMANN_TRRS_H
#define SWIFT_RIEMANN_TRRS_H

#include "adiabatic_index.h"

/* frequently used combinations of hydro_gamma */
#define const_riemann_gp1d2g (0.5f * (hydro_gamma + 1.0f) / hydro_gamma)
#define const_riemann_gm1d2g (0.5f * hydro_gamma_minus_one / hydro_gamma)
#define const_riemann_gm1dgp1 (hydro_gamma_minus_one / (hydro_gamma + 1.0f))
#define const_riemann_tdgp1 (2.0f / (hydro_gamma + 1.0f))
#define const_riemann_tdgm1 (2.0f / hydro_gamma_minus_one)
#define const_riemann_gm1d2 (0.5f * hydro_gamma_minus_one)
#define const_riemann_tgdgm1 (2.0f * hydro_gamma / hydro_gamma_minus_one)
#define const_riemann_ginv (1.0f / hydro_gamma)

/**
 * @brief Solve the Riemann problem using the Two Rarefaction Riemann Solver
 *
 * By assuming 2 rarefaction waves, we can analytically solve for the pressure
 * and velocity in the intermediate region, eliminating the iterative procedure.
 *
 * According to Toro: 'The two-rarefaction approximation is generally quite
 * robust; (...) The TRRS is in fact exact when both non-linear waves are
 * actually rarefaction waves.'
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param Whalf Empty state vector in which the result will be stored
 * @param n_unit Normal vector of the interface
 */
__attribute__((always_inline)) INLINE static void riemann_solver_solve(
    float* WL, float* WR, float* Whalf, float* n_unit) {
  float aL, aR;
  float PLR;
  float vL, vR;
  float ustar, pstar;
  float vhalf;
  float pdpR, SHR, STR;
  float pdpL, SHL, STL;

  /* calculate the velocities along the interface normal */
  vL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  vR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* calculate the sound speeds */
  aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* calculate the velocity and pressure in the intermediate state */
  PLR = pow(WL[4] / WR[4], const_riemann_gm1d2g);
  ustar = (PLR * vL / aL + vR / aR + const_riemann_tdgm1 * (PLR - 1.0f)) /
          (PLR / aL + 1.0f / aR);
  pstar = 0.5f * (WL[4] * pow(1.0f + const_riemann_gm1d2 / aL * (vL - ustar),
                              const_riemann_tgdgm1) +
                  WR[4] * pow(1.0f + const_riemann_gm1d2 / aR * (ustar - vR),
                              const_riemann_tgdgm1));

  /* sample the solution */
  if (ustar < 0.0f) {
    /* right state */
    Whalf[1] = WR[1];
    Whalf[2] = WR[2];
    Whalf[3] = WR[3];
    pdpR = pstar / WR[4];
    /* always a rarefaction wave, that's the approximation */
    SHR = vR + aR;
    if (SHR > 0.0f) {
      STR = ustar + aR * pow(pdpR, const_riemann_gm1d2g);
      if (STR <= 0.0f) {
        Whalf[0] =
            WR[0] * pow(const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR,
                        const_riemann_tdgm1);
        vhalf = const_riemann_tdgp1 * (-aR + const_riemann_gm1d2 * vR) - vR;
        Whalf[4] =
            WR[4] * pow(const_riemann_tdgp1 - const_riemann_gm1dgp1 / aR * vR,
                        const_riemann_tgdgm1);
      } else {
        Whalf[0] = WR[0] * pow(pdpR, const_riemann_ginv);
        vhalf = ustar - vR;
        Whalf[4] = pstar;
      }
    } else {
      Whalf[0] = WR[0];
      vhalf = 0.0f;
      Whalf[4] = WR[4];
    }
  } else {
    /* left state */
    Whalf[1] = WL[1];
    Whalf[2] = WL[2];
    Whalf[3] = WL[3];
    pdpL = pstar / WL[4];
    /* rarefaction wave */
    SHL = vL - aL;
    if (SHL < 0.0f) {
      STL = ustar - aL * pow(pdpL, const_riemann_gm1d2g);
      if (STL > 0.0f) {
        Whalf[0] =
            WL[0] * pow(const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL,
                        const_riemann_tdgm1);
        vhalf = const_riemann_tdgp1 * (aL + const_riemann_gm1d2 * vL) - vL;
        Whalf[4] =
            WL[4] * pow(const_riemann_tdgp1 + const_riemann_gm1dgp1 / aL * vL,
                        const_riemann_tgdgm1);
      } else {
        Whalf[0] = WL[0] * pow(pdpL, const_riemann_ginv);
        vhalf = ustar - vL;
        Whalf[4] = pstar;
      }
    } else {
      Whalf[0] = WL[0];
      vhalf = 0.0f;
      Whalf[4] = WL[4];
    }
  }

  /* add the velocity solution along the interface normal to the velocities */
  Whalf[1] += vhalf * n_unit[0];
  Whalf[2] += vhalf * n_unit[1];
  Whalf[3] += vhalf * n_unit[2];
}

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    float* Wi, float* Wj, float* n_unit, float* vij, float* totflux) {

  float Whalf[5];
  float flux[5][3];
  float vtot[3];
  float rhoe;

  riemann_solver_solve(Wi, Wj, Whalf, n_unit);

  flux[0][0] = Whalf[0] * Whalf[1];
  flux[0][1] = Whalf[0] * Whalf[2];
  flux[0][2] = Whalf[0] * Whalf[3];

  vtot[0] = Whalf[1] + vij[0];
  vtot[1] = Whalf[2] + vij[1];
  vtot[2] = Whalf[3] + vij[2];
  flux[1][0] = Whalf[0] * vtot[0] * Whalf[1] + Whalf[4];
  flux[1][1] = Whalf[0] * vtot[0] * Whalf[2];
  flux[1][2] = Whalf[0] * vtot[0] * Whalf[3];
  flux[2][0] = Whalf[0] * vtot[1] * Whalf[1];
  flux[2][1] = Whalf[0] * vtot[1] * Whalf[2] + Whalf[4];
  flux[2][2] = Whalf[0] * vtot[1] * Whalf[3];
  flux[3][0] = Whalf[0] * vtot[2] * Whalf[1];
  flux[3][1] = Whalf[0] * vtot[2] * Whalf[2];
  flux[3][2] = Whalf[0] * vtot[2] * Whalf[3] + Whalf[4];

  /* eqn. (15) */
  /* F_P = \rho e ( \vec{v} - \vec{v_{ij}} ) + P \vec{v} */
  /* \rho e = P / (\gamma-1) + 1/2 \rho \vec{v}^2 */
  rhoe = Whalf[4] / hydro_gamma_minus_one +
         0.5f * Whalf[0] *
             (vtot[0] * vtot[0] + vtot[1] * vtot[1] + vtot[2] * vtot[2]);
  flux[4][0] = rhoe * Whalf[1] + Whalf[4] * vtot[0];
  flux[4][1] = rhoe * Whalf[2] + Whalf[4] * vtot[1];
  flux[4][2] = rhoe * Whalf[3] + Whalf[4] * vtot[2];

  totflux[0] =
      flux[0][0] * n_unit[0] + flux[0][1] * n_unit[1] + flux[0][2] * n_unit[2];
  totflux[1] =
      flux[1][0] * n_unit[0] + flux[1][1] * n_unit[1] + flux[1][2] * n_unit[2];
  totflux[2] =
      flux[2][0] * n_unit[0] + flux[2][1] * n_unit[1] + flux[2][2] * n_unit[2];
  totflux[3] =
      flux[3][0] * n_unit[0] + flux[3][1] * n_unit[1] + flux[3][2] * n_unit[2];
  totflux[4] =
      flux[4][0] * n_unit[0] + flux[4][1] * n_unit[1] + flux[4][2] * n_unit[2];
}

#endif /* SWIFT_RIEMANN_TRRS_H */
