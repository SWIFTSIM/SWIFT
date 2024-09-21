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

#include "chemistry_getters.h"
#include "chemistry_slope_limiters_face.h"
#include "hydro.h"

/* TODO: Put this in the sign.h file and define for double. */
#define SIGN(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))

/**
 * @brief Check if the given input states are vacuum or will generate vacuum
 */
__attribute__((always_inline)) INLINE static int chemistry_riemann_is_vacuum(
    const float rhoL, const float rhoR, float vL, float vR, float aL,
    float aR) {

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
 * @brief Solve the Riemann problem for the diffusion equations and return the
 * flux at the interface.
 *
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density)
 * @param UR right diffusion state (metal density)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure)
 * @param F_diff_L The diffusion flux of the left
 * @param F_diff_R The diffusion flux of the right
 * @param Anorm Norm of the face between the left and right particles
 * @param hyperFluxR the flux of the hyperbolic conservation law of the right
 * state
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param metal_flux (return) The resulting flux at the interface
 * @param chem_data Chemistry data.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solve_for_flux(
    const struct part* restrict pi, const struct part* restrict pj,
    const double UL, const double UR, const float WL[5], const float WR[5],
    const double F_diff_L[3], const double F_diff_R[3], const float Anorm,
    const float n_unit[3], int g, double* metal_flux,
    const struct chemistry_global_data* chem_data) {

  /* Handle pure vacuum */
  if (!UL && !UR) {
    *metal_flux = 0.0f;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const float uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* Get the fastet speed of sound */
  /* Think about comoving vs physical units */
  const float c_s_L = hydro_get_comoving_soundspeed(pi);
  const float c_s_R = hydro_get_comoving_soundspeed(pj);
  const float c_fast = max(c_s_L, c_s_R);

  /* Approximate lambda_plus and lambda_minus. Use velocity difference. */
  const float lambda_plus = fabsf(uL - uR) + c_fast;
  const float lambda_minus = -lambda_plus;

  if (lambda_plus == 0.f && lambda_minus == 0.f) {
    *metal_flux = 0.f;
    return;
  }

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] +
                        F_diff_L[2] * n_unit[2];
  const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] +
                        F_diff_R[2] * n_unit[2];

  /* Compute F_HLL */
  const double dU = UR - UL;
  const float one_over_dl = 1.f / (lambda_plus - lambda_minus);
  const double F_2 =
      (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl;
  double F_U = lambda_plus * lambda_minus * dU * one_over_dl;

  if (chem_data->use_hokpins2017_hll_riemann_solver) {
    /****************************************************************************
     * Hopkins 2017 implementation of HLL */
    const float K_star_norm =
        0.5 * sqrtf(3.0) *
        fabsf(pi->chemistry_data.kappa + pj->chemistry_data.kappa);

    /* If the diffusion matrix is null, don't exchange flux. This can happen
       in the first timestep when the density is not yet computed. */
    if (K_star_norm == 0.0) {
      *metal_flux = 0;
      return;
    }

    const float dx[3] = {pj->x[0] - pi->x[0], pj->x[1] - pi->x[1],
                         pj->x[2] - pi->x[2]};
    const float dx_norm_2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const float dx_norm = sqrtf(dx_norm_2);
    const float r = dx_norm / K_star_norm * (0.5 * fabs(uR - uL) + c_fast);
    const float r_term = (0.2 + r) / (0.2 + r + r * r);
    const float norm_term = 1.0 / sqrtf(3.0);
    const float alpha = norm_term * r_term;

    /* Notice that the norm(U) and norm(q) are for U = (metal_1, metal_2,
       etc), the same for norm(grad otimes q) */

    /* Multiply by alpha in Hopkins' HLL version */
    /* This is needed, this ensures we are not overdiffusing */
    F_U *= alpha;

    double flux_hll = 0.0;
    /* The minmod is too restrictive in the diffusion */
    if (chem_data->hll_riemann_solver_psi < 0) {
      flux_hll = chemistry_minmod((1 + chem_data->hll_riemann_solver_psi) * F_2, F_2 + F_U);
    } else {
      flux_hll = F_2 + F_U;
    }

    /* Compute the direct fluxes */
    const double qi = chemistry_part_get_metal_density(pi, g);
    const double qj = chemistry_part_get_metal_density(pj, g);
    const double dq = qj - qi;
    const double nabla_o_q_dir[3] = {
        dx[0] * dq / dx_norm_2, dx[1] * dq / dx_norm_2, dx[2] * dq / dx_norm_2};
    const double kappa_mean =
        0.5 * (pi->chemistry_data.kappa + pj->chemistry_data.kappa);
    const double F_A_left_side[3] = {-kappa_mean * nabla_o_q_dir[0],
                                     -kappa_mean * nabla_o_q_dir[1],
                                     -kappa_mean * nabla_o_q_dir[2]};
    const double F_A_right_side[3] = {Anorm * dx[0] / dx_norm,
                                      Anorm * dx[1] / dx_norm,
                                      Anorm * dx[2] / dx_norm};

    const double F_times_A_dir = F_A_left_side[0] * F_A_right_side[0] +
                                 F_A_left_side[1] * F_A_right_side[1] +
                                 F_A_left_side[2] * F_A_right_side[2];

    /* Get F_HLL * A_ij */
    const double F_HLL_times_A = flux_hll * Anorm;

    /* Now, choose the righ flux to get F_diff_ij^* */
    const double epsilon = chem_data->hll_riemann_solver_epsilon;
    if (F_times_A_dir * F_HLL_times_A < 0 &&
        fabs(F_times_A_dir) > epsilon * fabs(F_HLL_times_A)) {
      *metal_flux = 0;
    } else {
      *metal_flux = flux_hll;
    }
  } else {
    /***************************************************************************
     * Simple HLL (compute F_diff_ij^*) */
    *metal_flux = F_2 + F_U;
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_RIEMANN_HLL_H */
