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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_RIEMANN_HLL_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_RIEMANN_HLL_H

#include "hydro.h"

#define SIGN(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))

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
    const struct part* restrict pi, const struct part* restrict pj, const double UL,
    const double UR, const float WL[5], const float WR[5], const double F_diff_L[3],
    const double F_diff_R[3], const float Anorm, const float n_unit[3], int g,
    double* metal_flux) {

  /* Handle pure vacuum */
  if (!UL && !UR) {
    *metal_flux = 0.0f;
    return;
  }
  
  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const float uR = WR[0] * n_unit[0] + WR[1] * n_unit[1] + WR[3] * n_unit[2];
  const float aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  const float aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (chemistry_riemann_is_vacuum(WL[0], WR[0], uL, uR, aL, aR)) {
    *metal_flux = 0.0f;
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = WL[0] + WR[0];
  const float abar = aL + aR;
  const float pPVRS =
    0.5f * ((WL[4] + WR[4]) - 0.25f * (uR - uL) * rhobar * abar);
  const float pstar = max(0.f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.0f;
  if (pstar > WL[4] && WL[4] > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
	       (pstar / WL[4] - 1.0f));
  }
  float qR = 1.0f;
  if (pstar > WR[4] && WR[4] > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
	       (pstar / WR[4] - 1.0f));
  }
  const float SLmuL = -aL * qL;
  const float SRmuR = aR * qR;
  const float Sstar =
    (WR[4] - WL[4] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
    (WL[0] * SLmuL - WR[0] * SRmuR);

  const float SL = uL - aL * qL;
  const float SR = uR - aR * qR;

  /* Get the fastet wavespeed */
  const float c_fast = max3(SL, SR, Sstar);

  /* Approximate lambda_plus and lambda_minus */
  const float lambda_plus = max(uL, uR) + c_fast;
  const float lambda_minus = min(uL, UR) - c_fast;

  /* Compute alpha */
  const float K_star_norm = 0.5 * sqrtf(3.0) * (pi->chemistry_data.kappa + pj->chemistry_data.kappa);
  const float dx[3] = {pj->x[0] - pi->x[0], pj->x[1] - pi->x[1], pj->x[1] - pi->x[1]};
  const float dx_norm_2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
  const float dx_norm = sqrtf(dx_norm_2);
  const float r = dx_norm/K_star_norm * (0.5*fabs(uR - uL) + Sstar);
  const float r_term = (0.2 + r )/ (0.2 + r + r*r);
  const float norm_term = 1.0 / sqrtf(3.0);
  const float alpha = norm_term * r_term;

  /* Compute F_HLL */
  const double dU = UR - UL;
  const float one_over_dl = 1.f / (lambda_plus - lambda_minus);

  /* Project to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0]*n_unit[0] + F_diff_L[1]*n_unit[1] + F_diff_L[2]*n_unit[2];
  const double Flux_R = F_diff_R[0]*n_unit[0] + F_diff_R[1]*n_unit[1] + F_diff_R[2]*n_unit[2];
  const double F_U = alpha*lambda_plus*lambda_minus*dU*one_over_dl;
  const double F_2 = (lambda_plus*Flux_L - lambda_minus*Flux_R)*one_over_dl;

  /* TODO: psi must be a user-defined parameter */
  const double psi = 0.1;
  const double flux_hll = chemistry_minmod((1+psi)*F_2, F_2 + F_U);

  /* Compute the direct fluxes */
  const double qj = pj->chemistry_data.metal_mass[g] / pj->chemistry_data.geometry.volume;
  const double qi = pi->chemistry_data.metal_mass[g] / pi->chemistry_data.geometry.volume;
  const double dq = qj - qi;
  const double nabla_o_q_dir[3] = {dx[0]* dq / dx_norm_2,
				   dx[1]* dq / dx_norm_2,
				   dx[2]* dq / dx_norm_2};
  const double kappa_mean = 0.5*(pi->chemistry_data.kappa + pj->chemistry_data.kappa);
  const double F_A_left_side[3] = {- kappa_mean * nabla_o_q_dir[0],
				   - kappa_mean * nabla_o_q_dir[1],
				   - kappa_mean * nabla_o_q_dir[2]};
  const double F_A_right_side[3] = { Anorm*dx[0]/dx_norm, Anorm*dx[1]/dx_norm, Anorm*dx[2]/dx_norm };

  const double F_times_A_dir = F_A_left_side[0]*F_A_right_side[0] + F_A_left_side[1]*F_A_right_side[1] + F_A_left_side[2]*F_A_right_side[2];

  /* Get F_HLL * A_ij */
  const double F_HLL_times_A = flux_hll*Anorm;

  /* Now, choose the righ flux to get F_diff_ij^* */
  /* TODO: This must be user-defined *\/ */
  const double epsilon = 0.5;
  if ((SIGN(F_times_A_dir) != SIGN(F_HLL_times_A))
      && fabs(F_times_A_dir) > epsilon*fabs(F_HLL_times_A)) {
    *metal_flux = 0;
  } else {
    *metal_flux = flux_hll;
  }

  /* Simple HLL */
  /* Compute F_diff_ij^* */
  /* if (SL >= 0) { */
  /*   *metal_flux = Flux_L; */
  /* } else if (SL <= 0 && SR >= 0) { */
  /*   *metal_flux = flux_hll; */
  /* } else /\* SR <= 0 *\/ { */
  /*   *metal_flux = Flux_R; */
  /* } */
}


#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_RIEMANN_HLL_H */