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

#include "chemistry_getters.h"
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
 * The minmod limiter.
 *
 * @param a Left slope
 * @param b Right slope
 */
__attribute__((always_inline)) INLINE static double chemistry_riemann_minmod(
    double a, double b) {
  if (a > 0 && b > 0) {
    return min(a, b);
  } else if (a < 0 && b < 0) {
    return max(a, b);
  } else {
    return 0.0;
  }
}

/**
 * @brief Computes the matrix diffusion matrix at the particles' interface. The
 * matrix is computed as K_\star = 0.5*(K_R + K_L).
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 * @param (return) K_star The diffusion matrix at the interface (in physical
 * units).
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_compute_K_star(const struct part *restrict pi,
                                 const struct part *restrict pj,
                                 const struct chemistry_global_data *chem_data,
                                 const struct cosmology *cosmo,
                                 double K_star[3][3]) {
  double KR[3][3], KL[3][3];
  chemistry_get_physical_matrix_K(pi, chem_data, cosmo, KR);
  chemistry_get_physical_matrix_K(pj, chem_data, cosmo, KL);

  /* Init K_star to 0.0. */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      K_star[i][j] = 0.0;
    }
  }

  /* Compute K_star = 0.5 * (KR + KL) */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      K_star[i][j] += 0.5 * (KR[i][j] + KL[i][j]);
    }
  }
}

/**
 * @brief Limit the diffusion matrix at the interface to aboir excessively
 * large diffusion coefficients.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param Anorm Surface area of the interface (in physical units).
 * @param cosmo The current cosmological model.
 * @param (return) K_star The norm of the diffusion matrix at the interface
 * (in physical units).
 * @param (return) K_star The diffusion matrix at the interface (in physical
 * units).
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_prevent_large_K_star(const struct part *restrict pi,
                                       const struct part *restrict pj,
                                       const float Anorm,
                                       const struct cosmology *cosmo,
                                       double *norm_K_star,
                                       double K_star[3][3]) {

  /* Prevent exessively large diffusion coefficients. This is taken from Gizmo.
   */
  const float min_dt =
      (pj->chemistry_data.flux_dt > 0.f)
          ? fminf(pi->chemistry_data.flux_dt, pj->chemistry_data.flux_dt)
          : pi->chemistry_data.flux_dt;
  const float min_mass = min(hydro_get_mass(pi), hydro_get_mass(pj));
  const float min_h = min(pi->h, pj->h) * cosmo->a * kernel_gamma;

  double mass_flux = Anorm * *norm_K_star / min_h * min_dt / min_mass;
  if (mass_flux > 0.25) {
    /* warning("Mass_flux > 0.25, reducing the diffusion coefficient"); */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K_star[i][j] *= 0.25 / mass_flux;
      }
    }
  }
  *norm_K_star *= 0.25 / mass_flux;
}

/**
 * @brief Compute the artificial diffusivity \alpha from eqn (7-11) in Hopkins
 * 2017. While solving the Rieman problem with HLL/GLF, we do not use \alpha=1
 * because the effective numerical diffusivity can exceed the physical
 * diffusivity \propto \| K \|.
 *
 * @param c_s_R Sound speed of the right state (in physical units).
 * @param c_s_L Sound speed of the left state (in physical units).
 * @param uR Velocity of the right state in the interface frame (in physical
 * units).
 * @param uL Velocity of the left state in the interface frame (in physical
 * units).
 * @param dx_norm Normalized distance between the left and right states (in
 * physical units).
 * @param q_star Metal mass fraction of metal specie "metal" at the interface
 * (in physical units).
 * @param U_star Metal mass density at the interface (in physical units).
 * @param K_star The 3x3 matrix representing anisotropic diffusion at the
 * interface (in physical units).
 * @param norm_K_star The norm of the diffusion matrix at the interface
 * (in physical units).
 * @param grad_q_star Gradient of metal mass fraction at thexs interface (in
 * physical units).
 *
 * @return The artificial diffusivity \alpha (in physical units).
 */
__attribute__((always_inline)) INLINE static double
chemistry_riemann_compute_alpha(const double c_s_R, const double c_s_L,
                                const double uR, const double uL,
                                const double dx_norm, const double q_star,
                                const double U_star, const double K_star[3][3],
                                const double norm_K_star,
                                const double grad_q_star[3]) {
  /* Everything is in physical units here. No need to convert. */

  /* Compute norm(K_star * grad_q_star) */
  double norm_K_star_times_grad_q_star = 0.0;
  double matrix_product = 0.0;

  for (int i = 0; i < 3; ++i) {
    /* Reset the temporary var */
    matrix_product = 0.0;
    for (int j = 0; j < 3; ++j) {
      /* Compute (K_star * grad_q_star)_i = Sum_j (K_star)_ij (grad_q_star)_j */
      matrix_product = K_star[i][j] * grad_q_star[j];
    }
    /* Add the product to the norm squared */
    norm_K_star_times_grad_q_star += matrix_product * matrix_product;
  }
  norm_K_star_times_grad_q_star = sqrtf(norm_K_star_times_grad_q_star);

  const double c_fast_star = 0.5 * (c_s_L + c_s_R);
  const double v_HLL = 0.5 * fabs(uR - uL) + c_fast_star;
  const double r =
      v_HLL * dx_norm * fabs(q_star) / (norm_K_star * fabs(U_star));
  const double r_term = (0.2 + r) / (0.2 + r + r * r);

  const double norm_grad_q_star =
      sqrtf(grad_q_star[0] * grad_q_star[0] + grad_q_star[1] * grad_q_star[1] +
            grad_q_star[2] * grad_q_star[2]);
  double norm_term =
      norm_K_star_times_grad_q_star / (norm_K_star * norm_grad_q_star);

  /* This behaviour is physically not correct. The correct physical behaviour
     is undetermined: we have 0/0 */
  if (norm_grad_q_star == 0.0) {
    norm_term = 1.0;
  }

  double alpha = norm_term * r_term;

  /* Treat pathological cases (physical solutions to the equation) */
  if (U_star == 0.0 || norm_K_star == 0.0) {
    alpha = 0.0;
  }

  return alpha;
}

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

  /* Handle pure vacuum */
  if (!UL && !UR) {
    *metal_flux = 0.0f;
    return;
  }

  /****************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* Obtain velocity in interface frame */
  const float uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const float uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];

  /* Get the fastet speed of sound. Use physical soundspeed */
  const float c_s_L = hydro_get_physical_soundspeed(pi, cosmo);
  const float c_s_R = hydro_get_physical_soundspeed(pj, cosmo);
  const float c_fast = max(c_s_L, c_s_R);

  /* Approximate lambda_plus and lambda_minus. Use velocity difference. */
  const float lambda_plus = fabsf(uL - uR) + c_fast;
  const float lambda_minus = -lambda_plus;

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
  double ZR = chemistry_get_metal_mass_fraction(pi, m);
  double ZL = chemistry_get_metal_mass_fraction(pj, m);
  chemistry_gradients_predict_Z(pi, pj, m, dx, cosmo, &ZR, &ZL);

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
  /* Now compute the HLL flux */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] +
                        F_diff_L[2] * n_unit[2];
  const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] +
                        F_diff_R[2] * n_unit[2];

  /* Compute variables to determine F_HLL */
  const double dU = UR - UL;
  const float one_over_dl = 1.f / (lambda_plus - lambda_minus);
  const double F_2 =
      (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl;
  double F_U = lambda_plus * lambda_minus * dU * one_over_dl;

  /* Multiply by alpha to limit numerical diffusion. */
  F_U *= alpha;

  /****************************************************************************
   * Hopkins 2017 implementation of HLL */
  /* No conversion to physical needed, everything is physical here */
  const double flux_hll = chemistry_riemann_minmod(
      (1 + chem_data->hll_riemann_solver_psi) * F_2, F_2 + F_U);

  /* Compute the direct fluxes */
  const double qi = chemistry_get_metal_mass_fraction(pi, m);
  const double qj = chemistry_get_metal_mass_fraction(pj, m);
  const double dq = qj - qi;
  const double nabla_o_q_dir[3] = {dx_p[0] * dq / dx_p_norm_2,
                                   dx_p[1] * dq / dx_p_norm_2,
                                   dx_p[2] * dq / dx_p_norm_2};
  const double kappa_mean =
      0.5 * (pi->chemistry_data.kappa + pj->chemistry_data.kappa);

#if !defined(GEAR_MF_HYPERBOLIC_DIFFUSION)
  const double F_A_left_side[3] = {-kappa_mean * nabla_o_q_dir[0],
                                   -kappa_mean * nabla_o_q_dir[1],
                                   -kappa_mean * nabla_o_q_dir[2]};
#else
  const double tau_L = pi->chemistry_data.tau;
  const double tau_R = pj->chemistry_data.tau;
  const double F_A_left_side[3] = {- 0.5*(F_diff_L[0]/tau_L + F_diff_R[0]/tau_R) - kappa_mean * nabla_o_q_dir[0],
                                   - 0.5*(F_diff_L[1]/tau_L + F_diff_R[1]/tau_R) - kappa_mean * nabla_o_q_dir[1],
                                   - 0.5*(F_diff_L[2]/tau_L + F_diff_R[2]/tau_R) - kappa_mean * nabla_o_q_dir[2]};
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

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_HLL_H */
