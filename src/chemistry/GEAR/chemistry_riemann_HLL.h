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
#ifndef SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H
#define SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H

#include "hydro.h"

/**
 * @brief Check if the given input states are vacuum or will generate vacuum
 */
__attribute__((always_inline)) INLINE static int chemistry_riemann_is_vacuum(
    const float rhoL, const float rhoR, float vL, float vR, float aL, float aR) {

  /* vacuum */
  if (!rhoL || !rhoR) return 1;

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
 * @brief Solve the Riemann problem for the RT equations and return the
 * flux at the interface.
 *
 * @param UL left state (radiation energy density, flux)
 * @param UR right state (radiation energy density, flux)
 * @param FLnorm the norm of the radiation flux of the left state
 * @param FRnorm the norm of the radiation flux of the right state
 * @param hyperFluxL the flux of the hyperbolic conservation law of the left
 * state
 * @param hyperFluxR the flux of the hyperbolic conservation law of the right
 * state
 * @param n_unit the unit vector perpendicular to the "intercell" surface.
 * @param flux_half (return) the resulting flux at the interface
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_solve_for_flux(
    const struct part* restrict pi, const struct part* restrict pj, const float UL,
    const float UR, const float n_unit[3], int g,
    const float Anorm, const float F_diff_L[3], const float F_diff_R[3], float* metal_flux) {
										   
  /* Handle pure vacuum */
  if (!UL && !UR) {
    *metal_flux = 0.0f;
    return;
  }

  const float rhoL = pi->rho;
  const float rhoR = pj->rho;
  const float vL[3] = {pi->v[0], pi->v[1], pi->v[2]};
  const float vR[3] = {pj->v[0], pj->v[1], pj->v[2]};
  const float PL = hydro_get_comoving_pressure(pi);
  const float PR = hydro_get_comoving_pressure(pj);
  
  /* STEP 0: obtain velocity in interface frame */
  const float uL = vL[0] * n_unit[0] + vL[1] * n_unit[1] + vL[2] * n_unit[2];
  const float uR = vR[0] * n_unit[0] + vR[1] * n_unit[1] + vR[2] * n_unit[2];
  const float aL = sqrtf(hydro_gamma * PR / rhoL);
  const float aR = sqrtf(hydro_gamma * PR / rhoR);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (chemistry_riemann_is_vacuum(rhoL, rhoR, uL, uR, aL, aR)) {
    *metal_flux = 0.0f;
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = rhoL + rhoR;
  const float abar = aL + aR;
  const float pPVRS =
    0.5f * ((PL + PR) - 0.25f * (uR - uL) * rhobar * abar);
  const float pstar = max(0.f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.0f;
  if (pstar > PL && PL > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
	       (pstar / PL - 1.0f));
  }
  float qR = 1.0f;
  if (pstar > PR && PR > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
	       (pstar / PR - 1.0f));
  }
  const float SLmuL = -aL * qL;
  const float SRmuR = aR * qR;
  const float Sstar =
    (PR - PL + rhoL * uL * SLmuL - rhoR * uR * SRmuR) /
    (rhoL * SLmuL - rhoR * SRmuR);

  const float SL = uL - aL * qL;
  const float SR = uR - aR * qR;

  /* Get the fastet wavespeed */
  const float c_fast = max3(SL, SR, Sstar);

  /* Approximate lambda_plus and lambda_minus */
  const float lambda_plus = max(uL, uR) + c_fast;
  const float lambda_minus = min(uL, UR) - c_fast;

  /* Compute r */
  /* const float nabla_o_q_L[3] = {0.0, 0.0, 0.0}; */
  /* const float nabla_o_q_R[3] = {0, 0, 0};  */
  /* const float grad_otimes_q_star[3] = {0.5*(nabla_o_q_L[0] + nabla_o_q_R[0]), */
  /* 				       0.5*(nabla_o_q_L[1] + nabla_o_q_R[1]), */
  /* 				       0.5*(nabla_o_q_L[2] + nabla_o_q_R[2])}; */

  /* const float K_star_norm = sqrtf(3.0) * (pi->chemistry_data.kappa + pj->chemistry_data.kappa)/2.f; */
  /* const float dx[3] = {pj->x[0] - pi->x[0], pj->x[1] - pi->x[1], pj->x[1] - pi->x[1]}; */
  /* const float dx_norm = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; */
  /* const float r = dx_norm/K_star_norm * (0.5*labs(vR - vL) + Sstar); */
  /* const float r_term = (0.2 + r )/ (0.2 + r + r*r); */
  /* const float norm_term = 0.0; */
  /* const float alpha = r_term * norm_term ; */
  const float alpha = 1.0;

  const float dU = UR - UL;
  const float one_over_dl = 1.f / (lambda_plus - lambda_minus);

  /* Project to reduce to a 1D Problem with 1 quantity */
  const float Flux_L = F_diff_L[0]*n_unit[0] + F_diff_L[1]*n_unit[1] + F_diff_L[2]*n_unit[2];
  const float Flux_R = F_diff_R[0]*n_unit[0] + F_diff_R[1]*n_unit[1] + F_diff_R[2]*n_unit[2];

  /* Compute F_diff_ij^* */
  const float flux_hll  = (lambda_plus*Flux_L - lambda_minus*Flux_R + alpha*lambda_plus*lambda_minus*dU)*one_over_dl;

  if (SL >= 0) {
    *metal_flux = Flux_L;
  } else if (SL <= 0 && SR >= 0) {
    *metal_flux = flux_hll;
  } else /* SR <= 0 */ {
    *metal_flux = Flux_R;
  }

}


#endif /* SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H */
