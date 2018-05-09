/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 */
__attribute__((always_inline)) INLINE static void riemann_solve_vacuum(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR,
    float* Whalf, const float* n_unit) {

  float SL, SR;
  float vhalf;

  if (!WR[0] && !WL[0]) {
    /* if both states are vacuum, the solution is also vacuum */
    Whalf[0] = 0.0f;
    Whalf[1] = 0.0f;
    Whalf[2] = 0.0f;
    Whalf[3] = 0.0f;
    Whalf[4] = 0.0f;
    return;
  }
  if (!WR[0]) {
    Whalf[1] = WL[1];
    Whalf[2] = WL[2];
    Whalf[3] = WL[3];
    /* vacuum right state */
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
        Whalf[4] =
            WL[4] * pow_two_gamma_over_gamma_minus_one(
                        hydro_two_over_gamma_plus_one +
                        hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
      } else {
        Whalf[0] = 0.0f;
        Whalf[1] = 0.0f;
        Whalf[2] = 0.0f;
        Whalf[3] = 0.0f;
        Whalf[4] = 0.0f;
        return;
      }
    } else {
      Whalf[0] = WL[0];
      vhalf = 0.0f;
      Whalf[4] = WL[4];
    }
  } else {
    if (!WL[0]) {
      Whalf[1] = WR[1];
      Whalf[2] = WR[2];
      Whalf[3] = WR[3];
      /* vacuum left state */
      if (-aR < vR) {
        SR = vR - hydro_two_over_gamma_minus_one * aR;
        if (SR >= 0.0f) {
          Whalf[0] = 0.0f;
          Whalf[1] = 0.0f;
          Whalf[2] = 0.0f;
          Whalf[3] = 0.0f;
          Whalf[4] = 0.0f;
          return;
        } else {
          Whalf[0] =
              WR[0] * pow_two_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one -
                          hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
          vhalf = hydro_two_over_gamma_plus_one *
                      (-aR + hydro_gamma_minus_one_over_two * vR) -
                  vR;
          Whalf[4] =
              WR[4] * pow_two_gamma_over_gamma_minus_one(
                          hydro_two_over_gamma_plus_one -
                          hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
        }
      } else {
        Whalf[0] = WR[0];
        vhalf = 0.0f;
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
        return;
      } else {
        if (SL >= 0.0f) {
          Whalf[1] = WL[1];
          Whalf[2] = WL[2];
          Whalf[3] = WL[3];
          if (aL > vL) {
            Whalf[0] = WL[0] *
                       pow_two_over_gamma_minus_one(
                           hydro_two_over_gamma_plus_one +
                           hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
            vhalf = hydro_two_over_gamma_plus_one *
                        (aL + hydro_gamma_minus_one_over_two * vL) -
                    vL;
            Whalf[4] = WL[4] *
                       pow_two_gamma_over_gamma_minus_one(
                           hydro_two_over_gamma_plus_one +
                           hydro_gamma_minus_one_over_gamma_plus_one / aL * vL);
          } else {
            Whalf[0] = WL[0];
            vhalf = 0.0f;
            Whalf[4] = WL[4];
          }
        } else {
          Whalf[1] = WR[1];
          Whalf[2] = WR[2];
          Whalf[3] = WR[3];
          if (-aR < vR) {
            Whalf[0] = WR[0] *
                       pow_two_over_gamma_minus_one(
                           hydro_two_over_gamma_plus_one -
                           hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
            vhalf = hydro_two_over_gamma_plus_one *
                        (-aR + hydro_gamma_minus_one_over_two * vR) -
                    vR;
            Whalf[4] = WR[4] *
                       pow_two_gamma_over_gamma_minus_one(
                           hydro_two_over_gamma_plus_one -
                           hydro_gamma_minus_one_over_gamma_plus_one / aR * vR);
          } else {
            Whalf[0] = WR[0];
            vhalf = 0.0f;
            Whalf[4] = WR[4];
          }
        }
      }
    }
  }

  /* Add the velocity solution along the interface normal to the velocities */
  Whalf[1] += vhalf * n_unit[0];
  Whalf[2] += vhalf * n_unit[1];
  Whalf[3] += vhalf * n_unit[2];
}

/**
 * @brief Solve the vacuum Riemann problem and return the fluxes
 */
__attribute__((always_inline)) INLINE static void riemann_solve_vacuum_flux(
    const float* WL, const float* WR, float vL, float vR, float aL, float aR,
    const float* n_unit, const float* vij, float* totflux) {

  float Whalf[5];
  float flux[5][3];
  float vtot[3];
  float rhoe;

  riemann_solve_vacuum(WL, WR, vL, vR, aL, aR, Whalf, n_unit);

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

#endif /* SWIFT_RIEMANN_VACUUM_H */
