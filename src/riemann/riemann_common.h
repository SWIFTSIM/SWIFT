//
// Created by yuyttenh on 28/02/23.
//

#ifndef SWIFT_RIEMANN_COMMON_H
#define SWIFT_RIEMANN_COMMON_H

__attribute__((always_inline)) INLINE static void riemann_flux_from_half_state(
    const float* Whalf, const float* vij, const float* n_unit, float* totflux) {

  float vtot[3];
  vtot[0] = Whalf[1] + vij[0];
  vtot[1] = Whalf[2] + vij[1];
  vtot[2] = Whalf[3] + vij[2];

  float flux[5][3];
  flux[0][0] = Whalf[0] * Whalf[1];
  flux[0][1] = Whalf[0] * Whalf[2];
  flux[0][2] = Whalf[0] * Whalf[3];

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
  float rhoe = Whalf[4] / hydro_gamma_minus_one +
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

#endif  // SWIFT_RIEMANN_COMMON_H
