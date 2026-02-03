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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_UTILS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_UTILS_H

/* Config parameters. */
#include <config.h>

/* Local headers */
#include "chemistry_getters.h"
#include "chemistry_gradients.h"
#include "part.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_riemann_utils.h
 * @brief File containing functions to modularise the Riemann solver
 * computations.
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
  double KL[3][3], KR[3][3];
  chemistry_get_physical_matrix_K(pi, chem_data, cosmo, KL);
  chemistry_get_physical_matrix_K(pj, chem_data, cosmo, KR);

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
chemistry_riemann_compute_alpha(const double c_s_L, const double c_s_R,
                                const double uL, const double uR,
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
      matrix_product += K_star[i][j] * grad_q_star[j];
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
      sqrt(grad_q_star[0] * grad_q_star[0] + grad_q_star[1] * grad_q_star[1] +
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

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
/**
 * @brief Compute the blending factor between the hyperbolic flux and the
 * parabolic flux to reduce numerical diffusion.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 *
 * @return The blending factor.
 */
__attribute__((always_inline)) INLINE static double
chemistry_riemann_compute_hyperbolic_blending_factor(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL[4], const double UR[4],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

  /* Get the correct diffusion driver depending on the diffusion mode */
  double qL, qR;
  if (chem_data->diffusion_mode == isotropic_constant) {
    /* For constant isotropic case, U = q = rho*Z.
       This is already predicted at the cell interface, nothing else to do. */
    qL = UL[0];
    qR = UR[0];
  } else {
    /* In these cases, U = rho*Z, q = Z */
    qL = chemistry_get_metal_mass_fraction(pi, m);
    qR = chemistry_get_metal_mass_fraction(pj, m);

    chemistry_gradients_predict_Z(pi, pj, m, dx, cosmo, &qL, &qR);
  }
  const double q_star = 0.5 * (qL + qR);
  const double U_star = 0.5 * (UR[0] + UL[0]);

  /* Compute diffusion matrix K_star = 0.5*(KR + KL) */
  double K_star[3][3];
  chemistry_riemann_compute_K_star(pi, pj, chem_data, cosmo, K_star);
  const double norm_K_star = chemistry_get_matrix_norm(K_star);
  const double norm_D_star = norm_K_star * q_star / U_star;
  const double delta_x = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const double tau_numerical = delta_x * delta_x / norm_D_star;

  /* Harmonic mean weigths better if we have widely different tau */
  const double tau_L = pi->chemistry_data.tau;
  const double tau_R = pj->chemistry_data.tau;
  const double tau_star = tau_L * tau_R / (tau_L + tau_R);
  const double ratio = tau_star / tau_numerical;

  /* This is a smooth function. When ratio << 1, we are in parabolic diffusion
     regime, then alpha ~ ratio. When ratio >> 1, we are in hyperbolic
     diffusion regime so alpha ~= 1.
     The 0.1 constant is chosen to reach alpha = 1 faster when ratio >= 1.*/
  double alpha = (ratio) / (0.1 + ratio);

  /* Safeguards */
  if (isnan(alpha) || isinf(alpha)) {
    alpha = 1.0;
  }
  /* \alpha \in [0, 1] so enforce this. (It should never happen) */
  if (alpha > 1.0) {
    alpha = 1.0;
  }
  if (alpha < 0.0) {
    alpha = 0.0;
  }

  return alpha;
}
#endif /* CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION */

/**
 * @brief Check if the given input states are vacuum or will generate vacuum.
 *
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param vL Left velocity (in physical units)
 * @param vR Right velocity (in physical units)
 * @param aL Left soundspeed (in physical units)
 * @param aR Right soundspeed (in physical units)
 */
__attribute__((always_inline)) INLINE static int chemistry_riemann_is_vacuum(
    const float *WL, const float *WR, float vL, float vR, float aL, float aR) {

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

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_RIEMANN_UTILS_H */
