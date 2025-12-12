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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_RIEMANN_HLL_H
#define SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_RIEMANN_HLL_H

#include "../chemistry_gradients.h"
#include "../chemistry_riemann_checks.h"
#include "../chemistry_riemann_utils.h"
#include "../chemistry_struct.h"
#include "../chemistry_properties.h"
#include "hydro.h"
#include "sign.h"
#include "minmax.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_riemann_HLL.h
 * @brief File containing functions concerning HLL riemann solver for the
 * anisotropic parabolic diffusion equation.
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
  const double u_L = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double u_R = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];
  const double u_rel = u_R - u_L;

  /* Get the fastet speed of sound. Use physical soundspeed */
  const double c_s_L = hydro_get_physical_soundspeed(pi, cosmo);
  const double c_s_R = hydro_get_physical_soundspeed(pj, cosmo);
  const double c_fast = max(c_s_L, c_s_R);

  /* Approximate lambda_plus and lambda_minus. Use velocity difference. */
  const double lambda_plus = fabs(u_rel) + c_fast;
  const double lambda_minus = -lambda_plus;

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (chemistry_riemann_is_vacuum(WL, WR, u_L, u_R, c_s_L, c_s_R)) {
    *metal_flux = 0.0f;
    return;
  }

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

  /* Get U_star. Already in physical units. */
  const double U_star = 0.5 * (UR + UL);

  /* Reconstruct qR and qL at the interface. Note that we have reconstructed UL
     and UR in chemistry_gradients_predict(), but not q. We have everything to
     do so now. */
  double qL, qR;
  if (chem_data->diffusion_mode == isotropic_constant) {
    /* For constant isotropic case, U = q = rho*Z.
       This is already predicted at the cell interface, nothing else to do. */
       qL = UL;
       qR = UR;
  } else {
    /* In these cases, U = rho*Z, q = Z */
    qL = chemistry_get_metal_mass_fraction(pi, m);
    qR = chemistry_get_metal_mass_fraction(pj, m);

    chemistry_gradients_predict_Z(pi, pj, m, dx, cosmo, &qL, &qR);
  }
  /* Now compute q_star and grad_q_star. Convert the gradient to physical
     units by dividing by a. Z is physical. */
  const double q_star = 0.5 * (qR + qL);

  /* Get the physical gradients */
  double grad_q_L[3], grad_q_R[3];
  if (chem_data->diffusion_mode == isotropic_constant) {
    /* In this case, q = U = rho*Z, grad q = grad (rho*Z) */
    chemistry_get_metal_density_gradients(pi, m, grad_q_L);
    chemistry_get_metal_density_gradients(pj, m, grad_q_R);

    /* Convert to physical units: multiply by a^-4 (a^-3 for density and a^-1
     * for gradient) */
    const double a_inv_4 =
        cosmo->a_inv * cosmo->a_inv * cosmo->a_inv * cosmo->a_inv;
    grad_q_L[0] *= a_inv_4;
    grad_q_L[1] *= a_inv_4;
    grad_q_L[2] *= a_inv_4;

    grad_q_R[0] *= a_inv_4;
    grad_q_R[1] *= a_inv_4;
    grad_q_R[2] *= a_inv_4;
  } else {
    /* In these cases, U = rho*Z, q = Z, grad q = grad Z */
    chemistry_get_metal_mass_fraction_gradients(pi, m, grad_q_L);
    chemistry_get_metal_mass_fraction_gradients(pj, m, grad_q_R);

    /* Convert to physical units: multiply by a^-1 (Z is already physical) */
    grad_q_L[0] *= cosmo->a_inv;
    grad_q_L[1] *= cosmo->a_inv;
    grad_q_L[2] *= cosmo->a_inv;

    grad_q_R[0] *= cosmo->a_inv;
    grad_q_R[1] *= cosmo->a_inv;
    grad_q_R[2] *= cosmo->a_inv;
  }
  double grad_q_star[3] = {0.5 * (grad_q_L[0] + grad_q_R[0]),
                           0.5 * (grad_q_L[1] + grad_q_R[1]),
                           0.5 * (grad_q_L[2] + grad_q_R[2])};

  /* Define some convenient variables. Convert to physical: add a for the norm
   */

  /* Re check carefully the signs of the flux */
  const float dx_p[3] = {dx[0] * cosmo->a, dx[1] * cosmo->a, dx[2] * cosmo->a};
  const float dx_p_norm_2 = dx_p[0] * dx_p[0] + dx_p[1] * dx_p[1] + dx_p[2] * dx_p[2];
  const float dx_p_norm = sqrtf(dx_p_norm_2);

  /* Now compute alpha to reduce numerical diffusion below physical
     diffusion. */
  const double alpha =
      chemistry_riemann_compute_alpha(c_s_L, c_s_R, u_L, u_R, dx_p_norm, q_star,
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
  double qi, qj;
  if (chem_data->diffusion_mode == isotropic_constant) {
    qi = chemistry_get_physical_metal_density(pi, m, cosmo);
    qj = chemistry_get_physical_metal_density(pj, m, cosmo);
  } else {
    qi = chemistry_get_metal_mass_fraction(pi, m);
    qj = chemistry_get_metal_mass_fraction(pj, m);
  }
  const double dq = qj - qi;
  const double grad_q_dir[3] = {dx_p[0] * dq / dx_p_norm_2,
				dx_p[1] * dq / dx_p_norm_2,
				dx_p[2] * dq / dx_p_norm_2};

  double F_dir_left_side[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      F_dir_left_side[i] += K_star[i][j]*grad_q_dir[j];
    }
  }

  const double F_dir_right_side[3] = {dx_p[0] / dx_p_norm,
				      dx_p[1] / dx_p_norm,
				      dx_p[2] / dx_p_norm};

  /* We don't need to multiply by Anorm since flux_hll is not yet multiplied by
     Anorm */
  const double F_dir = F_dir_left_side[0] * F_dir_right_side[0] +
		       F_dir_left_side[1] * F_dir_right_side[1] +
		       F_dir_left_side[2] * F_dir_right_side[2];

  /* Now, choose the righ flux to get F_diff_ij^* */
  const double epsilon = chem_data->hll_riemann_solver_epsilon;
  if (F_dir * flux_hll < 0.0 &&
      fabs(F_dir) > epsilon * fabs(flux_hll)) {
    *metal_flux = 0.0;
  } else {
#if !defined (GEAR_FVPM_DIFFUSION_FLUX_LIMITER_EXTREMA_PRESERVING)
    *metal_flux = flux_hll;
#else
    /* Diffusion should not create new extrema. We (try to) ensure this with
       the following flux limiter */
    const struct chemistry_part_data *chi = &pi->chemistry_data;
    const struct chemistry_part_data *chj = &pj->chemistry_data;

    const float mindt =
      (chj->flux_dt > 0.f) ? fminf(chi->flux_dt, chj->flux_dt) : chi->flux_dt;

    /* Current metal mass reservoir */
    const double mZi = chi->metal_mass[m];
    const double mZj = chj->metal_mass[m];

    /* Convert the flux to mass */
    const double mZ_exchanged = flux_hll * Anorm * mindt;

    /* Limits */
    const double rhoZ_max = min(chi->limiter.rhoZ[m][1], chj->limiter.rhoZ[m][1]);
    const double rhoZ_min = min(chi->limiter.rhoZ[m][0], chj->limiter.rhoZ[m][0]);

    const float Vi = pi->geometry.volume;
    const float Vj = pj->geometry.volume;
    double Phi_1, Phi_2, Phi;

    if (mZ_exchanged > 0.0) /* i loses mass to j */ {
      /* In this case, we want rhoZi_final >= rhoZ_min.
	 So, (mZi - Phi_1*mZ_exchanged)/Vi >= rhoZ_min
	 <=>  mZi - rhoZ_min*Vi >= Phi_1*mZ_exchanged
	 <=> Phi_1 <= (mZi - rhoZ_min*Vi)/mZ_exchanged.

	 But we also need rhoZj_final <= rhoZ_max.
	 So,  (mZj + Phi_2*mZ_exchanged)/Vj <= rhoZ_max
	 <=>  Phi_2*mZ_exchanged <= rhoZ_max*Vj - mZj
	 <=>  Phi_2 <= (rhoZ_max*Vj - mZj)/mZ_exchanged

	 Finally, we take Phi = min(Phi_1, Phi_2)

	 Note that we multiply the final inequations by a safety factor \in (0,
	 1] to take into account that there are other neighbours to process.
      */
      Phi_1 = (mZi - rhoZ_min * Vi) / mZ_exchanged;
      Phi_2 = (rhoZ_max * Vj - mZj) / mZ_exchanged;
    } else /* j loses mass i */ {
      /* In this case, we want rhoZj_final >= rhoZ_min.
	 So, (mZj - Phi_1*mZ_exchanged)/Vj >= rhoZ_min
	 <=>  mZj - rhoZ_min*Vj >= Phi_1*mZ_exchanged
	 <=> Phi_1 <= (mZj - rhoZ_min*Vj)/mZ_exchanged.

	 But we also need rhoZi_final <= rhoZ_max.
	 So,  (mZi + Phi_2*mZ_exchanged)/Vi <= rhoZ_max
	 <=>  Phi_2*mZ_exchanged <= rhoZ_max*Vi - mZi
	 <=>  Phi_2 <= (rhoZ_max*Vi - mZi)/mZ_exchanged

	 Finally, we take Phi = min(Phi_1, Phi_2)

	 Note that we multiply the final inequations by a safety factor \in (0,
	 1] to take into account that there are other neighbours to process.
      */
      Phi_1 = (mZj - rhoZ_min * Vj) / fabs(mZ_exchanged);
      Phi_2 = (rhoZ_max * Vi - mZi) / fabs(mZ_exchanged);
    }
    Phi = min(Phi_1, Phi_2);
    Phi = min(1.0, GEAR_FVPM_DIFFUSION_EXTREMA_PRESERVING_FLUX_LIMITER_SAFETY_FACTOR*Phi);

    /* Now rescale the flux */
    *metal_flux = Phi * flux_hll;

    if (GEAR_FVPM_DIFFUSION_FLUX_LIMITER_OUTPUT_VERBOSITY == 1 && Phi < 1.0) {
      const float rho_i = hydro_get_comoving_density(pi);
      const float rho_j = hydro_get_comoving_density(pj);
      message(
	  "[%lld, %lld] UL = %e, UR = %e, mZL = %e, mZR = %e, mZi = %e, mZj = "
	  "%e"
	  " | rho_i = %e, rho_j = %e, Vi = %e, Vj = %e | alpha = %e, F_L = %e, "
	  " F_R = %e, F_2 = %e, F_U = %e | flux_hll = %e,"
	  " F_dir = %e | Phi = %e, Phi_1 = %e, Phi_2 = %e | metal_mass_final = "
	  "%e",
	  pi->id, pj->id, UL, UR, UL * Vi, UR * Vj, mZi, mZj, rho_i, rho_j, Vi,
	  Vj, alpha, Flux_L * Anorm * mindt, Flux_R * Anorm * mindt,
	  F_2 * Anorm * mindt, F_U * Anorm * mindt, flux_hll * Anorm * mindt,
	  F_dir * Anorm * mindt, Phi, Phi_1, Phi_2,
	  Phi * flux_hll * Anorm * mindt);
    }
#endif /* GEAR_FVPM_DIFFUSION_FLUX_LIMITER_EXTREMA_PRESERVING */
  }
}

/* This prevents compilation issues with hyperbolic diffusion because of the
   different definitions of the following function and the riemann_check
   functions it calls */
#if !defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
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

  chemistry_riemann_check_input(WL, WR, UL, UR, n_unit);

  /* Handle pure vacuum */
  if ((!UL && !UR) || (!WL[0] && !WR[0])) {
    *metal_flux = 0.0;
    return;
  }

  /* No conversion to physical needed, everything is physical here */
  /* Hopkins Hopkins 2017 implementation of HLL */
  chemistry_riemann_solver_hopkins2017_HLL(dx, pi, pj, UL, UR, WL, WR, F_diff_L,
                                           F_diff_R, Anorm, n_unit, m,
                                           chem_data, cosmo, metal_flux);

  chemistry_riemann_check_output(WL, WR, UL, UR, n_unit, metal_flux);
}
#endif /* CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION */

#endif /* SWIFT_CHEMISTRY_GEAR_MF_PARABIOLIC_DIFFUSION_RIEMANN_HLL_H */
