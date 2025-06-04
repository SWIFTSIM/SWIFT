/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_RIEMANN_VACUUM_H
#define SWIFT_RIEMANN_VACUUM_H

#include "riemann_common.h"
/**
 * @brief Check if the given input states are vacuum or will generate vacuum
 */
__attribute__((always_inline)) INLINE static int riemann_is_vacuum(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR) {

  /* vacuum */
  if (!WL[0] || !WR[0]) return 1;

  /* vacuum generation */
  else if (hydro_two_over_gamma_minus_one * aL +
               hydro_two_over_gamma_minus_one * aR <=
           vR - vL)
    return 1;

  /* no vacuum */
  else
    return 0;
}

/**
 * @brief Vacuum Riemann solver, based on section 4.6 in Toro
 *
 * @param WL The left state vector
 * @param WR The right state vector
 * @param vL The left velocity along the interface normal
 * @param vR The right velocity along the interface normal
 * @param aL The left sound speed
 * @param aR The right sound speed
 * @param Whalf Empty state vector to store the solution in
 * @param n_unit Normal vector of the interface
 * @return The pressure in the middle state (P_star)
 */
__attribute__((always_inline)) INLINE static float riemann_solve_vacuum(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR,
    float* Whalf, const float* n_unit) {

  float SL, SR;
  float vhalf = 0.f;
  float P_star = 0.f;

  if (!WR[0] && !WL[0]) {
    /* if both states are vacuum, the solution is also vacuum */
    Whalf[0] = 0.0f;
    Whalf[1] = 0.0f;
    Whalf[2] = 0.0f;
    Whalf[3] = 0.0f;
    Whalf[4] = 0.0f;
  } else if (!WR[0]) {
    /* vacuum right state */
    Whalf[1] = WL[1];
    Whalf[2] = WL[2];
    Whalf[3] = WL[3];
    P_star = WL[4] * pow_two_gamma_over_gamma_minus_one(
                         hydro_two_over_gamma_plus_one +
                         hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
    if (vL < aL) {
      SL = vL + hydro_two_over_gamma_minus_one * aL;
      if (SL > 0.0f) {
        Whalf[0] =
            WL[0] * pow_two_over_gamma_minus_one(
                        hydro_two_over_gamma_plus_one +
                        hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
        vhalf = hydro_two_over_gamma_plus_one *
                    (aL + hydro_gamma_minus_one_over_two * vL) -
                vL;
        Whalf[4] = P_star;
      } else {
        Whalf[0] = 0.0f;
        Whalf[1] = 0.0f;
        Whalf[2] = 0.0f;
        Whalf[3] = 0.0f;
        Whalf[4] = 0.0f;
      }
    } else {
      Whalf[0] = WL[0];
      Whalf[4] = WL[4];
    }
  } else if (!WL[0]) {
    /* vacuum left state */
    Whalf[1] = WR[1];
    Whalf[2] = WR[2];
    Whalf[3] = WR[3];
    P_star = WR[4] * pow_two_gamma_over_gamma_minus_one(
                         hydro_two_over_gamma_plus_one -
                         hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
    if (-aR < vR) {
      SR = vR - hydro_two_over_gamma_minus_one * aR;
      if (SR >= 0.0f) {
        Whalf[0] = 0.0f;
        Whalf[1] = 0.0f;
        Whalf[2] = 0.0f;
        Whalf[3] = 0.0f;
        Whalf[4] = 0.0f;
      } else {
        Whalf[0] =
            WR[0] * pow_two_over_gamma_minus_one(
                        hydro_two_over_gamma_plus_one -
                        hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
        vhalf = hydro_two_over_gamma_plus_one *
                    (-aR + hydro_gamma_minus_one_over_two * vR) -
                vR;
        Whalf[4] = P_star;
      }
    } else {
      Whalf[0] = WR[0];
      Whalf[4] = WR[4];
    }
  } else {
    /* vacuum generation */
    SR = vR - hydro_two_over_gamma_minus_one * aR;
    SL = vL + hydro_two_over_gamma_minus_one * aL;
    if (SR > 0.0f && SL < 0.0f) {
      Whalf[0] = 0.0f;
      Whalf[1] = 0.0f;
      Whalf[2] = 0.0f;
      Whalf[3] = 0.0f;
      Whalf[4] = 0.0f;
    } else {
      if (SL >= 0.0f) {
        Whalf[1] = WL[1];
        Whalf[2] = WL[2];
        Whalf[3] = WL[3];
        P_star =
            WL[4] * pow_two_gamma_over_gamma_minus_one(
                        hydro_two_over_gamma_plus_one +
                        hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
        if (aL > vL) {
          Whalf[0] =
              WL[0] * pow_two_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one +
                          hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
          vhalf = hydro_two_over_gamma_plus_one *
                      (aL + hydro_gamma_minus_one_over_two * vL) -
                  vL;
          Whalf[4] = P_star;
        } else {
          Whalf[0] = WL[0];
          Whalf[4] = WL[4];
        }
      } else {
        Whalf[1] = WR[1];
        Whalf[2] = WR[2];
        Whalf[3] = WR[3];
        P_star =
            WR[4] * pow_two_gamma_over_gamma_minus_one(
                        hydro_two_over_gamma_plus_one -
                        hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
        if (-aR < vR) {
          Whalf[0] =
              WR[0] * pow_two_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one -
                          hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
          vhalf = hydro_two_over_gamma_plus_one *
                      (-aR + hydro_gamma_minus_one_over_two * vR) -
                  vR;
          Whalf[4] = P_star;
        } else {
          Whalf[0] = WR[0];
          Whalf[4] = WR[4];
        }
      }
    }
  }

  /* Add the velocity solution along the interface normal to the velocities */
  Whalf[1] += vhalf * n_unit[0];
  Whalf[2] += vhalf * n_unit[1];
  Whalf[3] += vhalf * n_unit[2];

  return P_star;
}

/**
 * @brief Solve the vacuum Riemann problem and return the fluxes
 * @return The pressure of the middle state (P_star).
 */
__attribute__((always_inline)) INLINE static float riemann_solve_vacuum_flux(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR,
    const float* n_unit, const float* vij, float* totflux) {

  float Whalf[5];
  float P_star = riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);
  riemann_flux_from_half_state(Whalf, vij, n_unit, totflux);
  return P_star;
}

#endif /* SWIFT_RIEMANN_VACUUM_H */
