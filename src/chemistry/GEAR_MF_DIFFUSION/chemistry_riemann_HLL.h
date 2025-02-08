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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_HLL_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_HLL_H

#include "chemistry_riemann_checks.h"
#include "chemistry_riemann_utils.h"
#include "chemistry_gradients.h"
#include "chemistry_struct.h"
#include "hydro.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_riemann_HLL.h
 * @brief File containing functions concerning HLL riemann solver for the
 * anisotropic diffusion equation.
 *
 * */

/**
 * @brief Hokpins (2017) HLL Riemann solver. It works well for parabolic
 * diffusion but not for hyperbolic diffusion when tau >> 1.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param F_diff_L The diffusion flux of the left (in physical units)
 * @param F_diff_R The diffusion flux of the right (in physical units)
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param metal_flux (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solver_hopkins2017_HLL(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL, const double UR,
    const float WL[5], const float WR[5], const double F_diff_L[3],
    const double F_diff_R[3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double *metal_flux) {

  /****************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* Obtain velocity in interface frame */
  const double uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* Get the fastet speed of sound. Use physical soundspeed */
  const double c_s_L = hydro_get_physical_soundspeed(pi, cosmo);
  const double c_s_R = hydro_get_physical_soundspeed(pj, cosmo);
  const double c_fast = max(c_s_L, c_s_R);

  /* Approximate lambda_plus and lambda_minus. Use velocity difference. */
  const double lambda_plus = fabs(uL - uR) + c_fast;
  const double lambda_minus = -lambda_plus;

  if (lambda_plus == 0.f && lambda_minus == 0.f) {
    *metal_flux = 0.f;
    return;
  }

  /****************************************************************************/
  /* Compute the flux artificial diffusion coefficient alpha */
  /* Compute diffusion matrix K_star = 0.5*(KR + KL) */
  double K_star[3][3];
  chemistry_riemann_compute_K_star(pi, pj, chem_data, cosmo, K_star);
  double norm_K_star = chemistry_get_matrix_norm(K_star);

  /* If the diffusion matrix is null, don't exchange flux. This can happen
     in the first timestep when the density is not yet computed. */
  if (norm_K_star == 0.0) {
    *metal_flux = 0;
    return;
  }

  /* Prevent exessively large diffusion coefficients (same as Gizmo) */
  /* chemistry_riemann_prevent_large_K_star(pi, pj, Anorm, cosmo, &norm_K_star,
   * K_star); */

  /* Get U_star. Already in physical units. */
  const double U_star = 0.5 * (UR + UL);

  /* Reconstruct ZR and ZL at the interface. Note that we have reconstructed UL
     and UR in chemistry_gradients_predict(), but not q. We have everything to
     do so now. */
  double ZL = chemistry_get_metal_mass_fraction(pi, m);
  double ZR = chemistry_get_metal_mass_fraction(pj, m);
  chemistry_gradients_predict_Z(pi, pj, m, dx, cosmo, &ZL, &ZR);

  /* Now compute q_star and grad_q_star. Convert the gradient to physical
     units by dividing by a. Z is physical. */
  const double q_star = 0.5 * (ZR + ZL);
  double grad_q_star[3] = {0.5 * cosmo->a_inv *
                               (pi->chemistry_data.gradients.Z[m][0] +
                                pj->chemistry_data.gradients.Z[m][0]),
                           0.5 * cosmo->a_inv *
                               (pi->chemistry_data.gradients.Z[m][1] +
                                pj->chemistry_data.gradients.Z[m][1]),
                           0.5 * cosmo->a_inv *
                               (pi->chemistry_data.gradients.Z[m][2] +
                                pj->chemistry_data.gradients.Z[m][2])};

  /* Define some convenient variables. Convert to physical: add a for the norm
   */
  const float dx_p[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a, dx[2] * cosmo->a};
  const float dx_p_norm_2 =
      sqrtf(dx_p[0] * dx_p[0] + dx_p[1] * dx_p[1] + dx_p[2] * dx_p[2]);
  const float dx_p_norm = sqrtf(dx_p_norm_2);

  /* Now compute alpha to reduce numerical diffusion below physical
     diffusion. */
  const double alpha =
      chemistry_riemann_compute_alpha(c_s_R, c_s_L, uR, uL, dx_p_norm, q_star,
                                      U_star, K_star, norm_K_star, grad_q_star);

  /****************************************************************************/
  /* Now project the fluxes */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] +
                        F_diff_L[2] * n_unit[2];
  const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] +
                        F_diff_R[2] * n_unit[2];

  /****************************************************************************/
  /* Compute variables to determine F_HLL */
  const double dU = UR - UL;
  const double one_over_dl = 1.f / (lambda_plus - lambda_minus);
  const double F_2 =
      (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl;
  const double F_U = lambda_plus * lambda_minus * dU * one_over_dl;
  const double flux_hll = chemistry_riemann_minmod(
      (1 + chem_data->hll_riemann_solver_psi) * F_2, F_2 + alpha * F_U);

  /* Compute the direct fluxes */
  const double qi = chemistry_get_metal_mass_fraction(pi, m);
  const double qj = chemistry_get_metal_mass_fraction(pj, m);
  const double dq = qj - qi;
  const double nabla_o_q_dir[3] = {dx_p[0] * dq / dx_p_norm_2,
                                   dx_p[1] * dq / dx_p_norm_2,
                                   dx_p[2] * dq / dx_p_norm_2};
  const double kappa_mean =
      0.5 * (pi->chemistry_data.kappa + pj->chemistry_data.kappa);

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  const float min_dt =
      (pj->chemistry_data.flux_dt > 0.f)
          ? fminf(pi->chemistry_data.flux_dt, pj->chemistry_data.flux_dt)
          : pi->chemistry_data.flux_dt;
  const double tau_L = pi->chemistry_data.tau;
  const double tau_R = pj->chemistry_data.tau;
  const double F_A_left_side[3] = {
      -min_dt * 0.5 * (F_diff_L[0] / tau_L + F_diff_R[0] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[0],
      -min_dt * 0.5 * (F_diff_L[1] / tau_L + F_diff_R[1] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[1],
      -min_dt * 0.5 * (F_diff_L[2] / tau_L + F_diff_R[2] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[2]};
#else
  const double F_A_left_side[3] = {-kappa_mean * nabla_o_q_dir[0],
                                   -kappa_mean * nabla_o_q_dir[1],
                                   -kappa_mean * nabla_o_q_dir[2]};
#endif

  const double F_A_right_side[3] = {Anorm * dx_p[0] / dx_p_norm,
                                    Anorm * dx_p[1] / dx_p_norm,
                                    Anorm * dx_p[2] / dx_p_norm};

  const double F_times_A_dir = F_A_left_side[0] * F_A_right_side[0] +
                               F_A_left_side[1] * F_A_right_side[1] +
                               F_A_left_side[2] * F_A_right_side[2];

  /* Get F_HLL * A_ij */
  const double F_HLL_times_A = flux_hll * Anorm;

  /* Now, choose the righ flux to get F_diff_ij^* */
  const double epsilon = chem_data->hll_riemann_solver_epsilon;
  if (F_times_A_dir * F_HLL_times_A < 0.0 &&
      fabs(F_times_A_dir) > epsilon * fabs(F_HLL_times_A)) {
    *metal_flux = 0;
  } else {
    *metal_flux = flux_hll;
  }
}

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
/**
 * @brief Hokpins (2017) HLL Riemann solver improved for hyperbolic diffusion.
 *
 * Note: To ensure exact metal mass conservation, the fluxes must be perfectly
 * antisymmetric. Hence, the "normal" HLL is not suited.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param F_diff_L The diffusion flux of the left (in physical units)
 * @param F_diff_R The diffusion flux of the right (in physical units)
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param metal_flux (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solver_hopkins2017_hyperbolic_HLL(
    const float dx[3],  struct part *restrict pi,
    struct part *restrict pj, const double UL, const double UR,
    const float WL[5], const float WR[5], const double F_diff_L[3],
    const double F_diff_R[3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double *metal_flux) {

  /****************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* Obtain velocity in interface frame */
  const double u_L = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double u_R = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];
  const double u_rel = u_R - u_L;

  /* Get the fastet speed of sound. Use physical soundspeed */
  /* const double c_s_L = hydro_get_physical_soundspeed(pi, cosmo); */
  /* const double c_s_R = hydro_get_physical_soundspeed(pj, cosmo); */
  const double c_s_L = chemistry_get_physical_hyperbolic_soundspeed(pi, chem_data, cosmo);
  const double c_s_R = chemistry_get_physical_hyperbolic_soundspeed(pj, chem_data, cosmo);
  const double c_fast = max(c_s_L, c_s_R);
  const double lambda_plus = fabs(u_rel) + c_fast;
  const double lambda_minus = -lambda_plus;

  const double aL = c_s_L;
  const double aR = c_s_R;

  /* Toro's wavespeed */
  /* const double lambda_plus = max(u_L + c_s_L, u_R + c_s_R); */
  /* const double lambda_minus = min(u_L - c_s_L, u_R - c_s_R); */

  /* PVRS wavespeed approximation */
  /* const double u_L = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2]; */
  /* const double u_R = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2]; */
  /* const double rhoLinv = (WL[0] > 0.0f) ? 1.0f / WL[0] : 0.0f; */
  /* const double rhoRinv = (WR[0] > 0.0f) ? 1.0f / WR[0] : 0.0f; */
  /* const double aL = sqrtf(hydro_gamma * WL[4] * rhoLinv); */
  /* const double aR = sqrtf(hydro_gamma * WR[4] * rhoRinv); */

  /* /\* Pressure estimate *\/ */
  /* const double rhobar = WL[0] + WR[0]; */
  /* const double abar = aL + aR; */
  /* const double pPVRS = */
  /*     0.5f * ((WL[4] + WR[4]) - 0.25f * (u_R - u_L) * rhobar * abar); */
  /* const double pstar = max(0.0f, pPVRS); */

  /* /\* Wave speed estimates */
  /*    all these speeds are along the interface normal, since u_L and u_R are *\/ */
  /* double qL = 1.0f; */
  /* if (pstar > WL[4] && WL[4] > 0.0f) { */
  /*   qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma * */
  /*                         (pstar / WL[4] - 1.0f)); */
  /* } */
  /* double qR = 1.0f; */
  /* if (pstar > WR[4] && WR[4] > 0.0f) { */
  /*   qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma * */
  /*                         (pstar / WR[4] - 1.0f)); */
  /* } */
  /* const double lambda_minus = -aL * qL; */
  /* const double lambda_plus = aR * qR; */

  pi->chemistry_data.wavespeed = max(pi->chemistry_data.wavespeed, lambda_plus);
  pj->chemistry_data.wavespeed = max(pj->chemistry_data.wavespeed, lambda_plus);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, u_L, u_R, aL, aR)) {
    *metal_flux = 0.0f;
    return;
  }

  if (lambda_plus == 0.f && lambda_minus == 0.f) {
    *metal_flux = 0.f;
    return;
  }

  /****************************************************************************/
  /* Define some convenient variables. Convert to physical: add a for the norm
   */
  double K_star[3][3];
  chemistry_riemann_compute_K_star(pi, pj, chem_data, cosmo, K_star);
  double norm_K_star = chemistry_get_matrix_norm(K_star);

  /* If the diffusion matrix is null, don't exchange flux. This can happen
     in the first timestep when the density is not yet computed. */
  if (norm_K_star == 0.0) {
    *metal_flux = 0;
    return;
  }

  const float dx_p[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a, dx[2] * cosmo->a};
  const float dx_p_norm_2 =
      sqrtf(dx_p[0] * dx_p[0] + dx_p[1] * dx_p[1] + dx_p[2] * dx_p[2]);
  const float dx_p_norm = sqrtf(dx_p_norm_2);

  /****************************************************************************/
  /* Now project the fluxes */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] +
                        F_diff_L[2] * n_unit[2];
  const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] +
                        F_diff_R[2] * n_unit[2];

  /****************************************************************************/
  /* Compute variables to determine F_HLL */
  const double dU = UR - UL;
  const double one_over_dl = 1.f / (lambda_plus - lambda_minus);
  const double F_2 =
      (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl;
  const double F_U = lambda_plus * lambda_minus * dU * one_over_dl;
  const double flux_hll = F_2 + F_U;

  /* Compute the direct fluxes */
  const double qi = chemistry_get_metal_mass_fraction(pi, m);
  const double qj = chemistry_get_metal_mass_fraction(pj, m);
  const double dq = qj - qi;
  const double nabla_o_q_dir[3] = {dx_p[0] * dq / dx_p_norm_2,
                                   dx_p[1] * dq / dx_p_norm_2,
                                   dx_p[2] * dq / dx_p_norm_2};
  const double kappa_mean =
      0.5 * (pi->chemistry_data.kappa + pj->chemistry_data.kappa);

  const float min_dt =
      (pj->chemistry_data.flux_dt > 0.f)
          ? fminf(pi->chemistry_data.flux_dt, pj->chemistry_data.flux_dt)
          : pi->chemistry_data.flux_dt;
  const double tau_L = pi->chemistry_data.tau;
  const double tau_R = pj->chemistry_data.tau;
  const double F_A_left_side[3] = {
      -min_dt * 0.5 * (F_diff_L[0] / tau_L + F_diff_R[0] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[0],
      -min_dt * 0.5 * (F_diff_L[1] / tau_L + F_diff_R[1] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[1],
      -min_dt * 0.5 * (F_diff_L[2] / tau_L + F_diff_R[2] / tau_R) -
          min_dt * kappa_mean * nabla_o_q_dir[2]};

  const double F_A_right_side[3] = {Anorm * dx_p[0] / dx_p_norm,
                                    Anorm * dx_p[1] / dx_p_norm,
                                    Anorm * dx_p[2] / dx_p_norm};

  const double F_times_A_dir = F_A_left_side[0] * F_A_right_side[0] +
                               F_A_left_side[1] * F_A_right_side[1] +
                               F_A_left_side[2] * F_A_right_side[2];

  /* Get F_HLL * A_ij */
  const double F_HLL_times_A = flux_hll * Anorm;

  /* Now, choose the righ flux to get F_diff_ij^* */
  const double epsilon = chem_data->hll_riemann_solver_epsilon;
  if (F_times_A_dir * F_HLL_times_A < 0.0 &&
      fabs(F_times_A_dir) > epsilon * fabs(F_HLL_times_A)) {
    *metal_flux = 0;
  } else {
    *metal_flux = flux_hll;
  }
}

/**
 * @brief HLL riemann solver.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param F_diff_L The diffusion flux of the left (in physical units)
 * @param F_diff_R The diffusion flux of the right (in physical units)
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param metal_flux (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_solver_HLL(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL, const double UR,
    const float WL[5], const float WR[5], const double F_diff_L[3],
    const double F_diff_R[3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double *metal_flux) {

  /***************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* PVRS wavespeed approximation */
  const double uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];
  const double rhoLinv = (WL[0] > 0.0f) ? 1.0f / WL[0] : 0.0f;
  const double rhoRinv = (WR[0] > 0.0f) ? 1.0f / WR[0] : 0.0f;
  const double aL = sqrtf(hydro_gamma * WL[4] * rhoLinv);
  const double aR = sqrtf(hydro_gamma * WR[4] * rhoRinv);

  /* Pressure estimate */
  const double rhobar = WL[0] + WR[0];
  const double abar = aL + aR;
  const double pPVRS =
      0.5f * ((WL[4] + WR[4]) - 0.25f * (uR - uL) * rhobar * abar);
  const double pstar = max(0.0f, pPVRS);

  /* Wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  double qL = 1.0f;
  if (pstar > WL[4] && WL[4] > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WL[4] - 1.0f));
  }
  double qR = 1.0f;
  if (pstar > WR[4] && WR[4] > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WR[4] - 1.0f));
  }
  const double lambda_minus = -aL * qL;
  const double lambda_plus = aR * qR;

  if (lambda_plus == 0.0 && lambda_minus == 0.0) {
    *metal_flux = 0.0;
  }

  /***************************************************************************/
  /* Now project the fluxes */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] +
                        F_diff_L[2] * n_unit[2];
  const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] +
                        F_diff_R[2] * n_unit[2];

  /***************************************************************************/
  /* Now solve the Riemann problem */

  const double dU = UR - UL;

  const double delta_lambda = lambda_plus - lambda_minus;
  if (fabs(delta_lambda) < 1e-8) {
    *metal_flux = 0.0;
    return;
  }

  const double one_over_dl = 1.f / (lambda_plus - lambda_minus);
  const double F_2 =
      (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl;
  double F_U = lambda_plus * lambda_minus * dU * one_over_dl;

  if (lambda_minus > 0.0) {
    *metal_flux = Flux_L;
  } else if (lambda_minus <= 0.0 && lambda_plus >= 0.0) {
    *metal_flux = F_2 + F_U;
  } else if (lambda_plus < 0.0) {
    *metal_flux = Flux_R;
  }
}
#endif

/**
 * @brief Solve the Riemann problem for the diffusion equations and return the
 * flux at the interface.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param F_diff_L The diffusion flux of the left (in physical units)
 * @param F_diff_R The diffusion flux of the right (in physical units)
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param metal_flux (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solve_for_flux(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL, const double UR,
    const float WL[5], const float WR[5], const double F_diff_L[3],
    const double F_diff_R[3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double *metal_flux) {

  riemann_check_input(WL, WR, UL, UR, n_unit);

  /* Handle pure vacuum */
  if ((!UL && !UR) || (!WL[0] && !WR[0])) {
    *metal_flux = 0.0;
    return;
  }

  /* No conversion to physical needed, everything is physical here */
  /* We probably need a mechanism to switch the Riemann solver when tau << 1.*/
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  if (chem_data->riemann_solver == 1) {
    /* Regular HLL */
    chemistry_riemann_solver_HLL(dx, pi, pj, UL, UR, WL, WR, F_diff_L, F_diff_R,
                                 Anorm, n_unit, m, chem_data, cosmo,
                                 metal_flux);
    return;

  } else if (chem_data->riemann_solver == 2) {
    /* Hopkins Hopkins 2017 implementation of HLL */
    chemistry_riemann_solver_hopkins2017_HLL(dx, pi, pj, UL, UR, WL, WR,
                                             F_diff_L, F_diff_R, Anorm, n_unit,
                                             m, chem_data, cosmo, metal_flux);
    return;
  } else if (chem_data->riemann_solver == 3) {
    chemistry_riemann_solver_hopkins2017_hyperbolic_HLL(dx, pi, pj, UL, UR, WL, WR,
							F_diff_L, F_diff_R, Anorm, n_unit,
							m, chem_data, cosmo, metal_flux);
    return;
  } else {
    return;
  }
#else
  /* Hopkins Hopkins 2017 implementation of HLL */
  chemistry_riemann_solver_hopkins2017_HLL(dx, pi, pj, UL, UR, WL, WR,
					   F_diff_L, F_diff_R, Anorm, n_unit,
					   m, chem_data, cosmo, metal_flux);
#endif

  riemann_check_output(WL, WR, UL, UR, n_unit, metal_flux);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_HLL_H */
