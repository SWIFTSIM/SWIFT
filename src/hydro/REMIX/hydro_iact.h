/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_REMIX_HYDRO_IACT_H
#define SWIFT_REMIX_HYDRO_IACT_H

/**
 * @file REMIX/hydro_iact.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2025)
 */

#include "adiabatic_index.h"
#include "const.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "hydro_visc_difn.h"
#include "math.h"
#include "minmax.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r */
  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  hydro_runner_iact_density_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  float wi, wi_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  hydro_runner_iact_nonsym_density_extra_kernel(pi, pj, dx, wi, wi_dx);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute kernel of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute kernel of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  hydro_runner_iact_gradient_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_gradient_extra_visc_difn(pi, pj, dx, wi, wj, wi_dx, wj_dx);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute kernel of pi. */
  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute kernel of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  hydro_runner_iact_nonsym_gradient_extra_kernel(pi, pj, dx, wi, wj, wi_dx,
                                                 wj_dx);
  hydro_runner_iact_nonsym_gradient_extra_visc_difn(pi, pj, dx, wi, wi_dx);
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float rhoi_inv = 1.f / rhoi;
  const float rhoj_inv = 1.f / rhoj;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* Linear-order reproducing kernel gradient term (Sandnes+2025 Eqn. 28) */
  float Gj[3], Gi[3], G_mean[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);

  /* Antisymmetric kernel grad term for conservation of momentum and energy */
  G_mean[0] = 0.5f * (Gi[0] - Gj[0]);
  G_mean[1] = 0.5f * (Gi[1] - Gj[1]);
  G_mean[2] = 0.5f * (Gi[2] - Gj[2]);

  /* Viscous pressures (Sandnes+2025 Eqn. 41) */
  float Qi, Qj;
  float visc_signal_velocity, difn_signal_velocity;
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &difn_signal_velocity, pi,
                  pj, dx);

  /* Pressure terms to be used in evolution equations */
  const float P_i_term = pressurei * rhoi_inv * rhoj_inv;
  const float P_j_term = pressurej * rhoi_inv * rhoj_inv;
  const float Q_i_term = Qi * rhoi_inv * rhoj_inv;
  const float Q_j_term = Qj * rhoi_inv * rhoj_inv;

  /* Use the force Luke! */
  pi->a_hydro[0] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[0];
  pi->a_hydro[1] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[1];
  pi->a_hydro[2] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[2];

  pj->a_hydro[0] +=
      mi * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[0];
  pj->a_hydro[1] +=
      mi * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[1];
  pj->a_hydro[2] +=
      mi * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[2];

  /* v_ij dot kernel gradient term */
  const float dvdotG = (pi->v[0] - pj->v[0]) * G_mean[0] +
                       (pi->v[1] - pj->v[1]) * G_mean[1] +
                       (pi->v[2] - pj->v[2]) * G_mean[2];

  /* Internal energy time derivative */
  pi->u_dt += mj * (P_i_term + Q_i_term) * dvdotG;
  pj->u_dt += mi * (P_j_term + Q_j_term) * dvdotG;

  /* Density time derivative */
  pi->drho_dt += mj * (rhoi * rhoj_inv) * dvdotG;
  pj->drho_dt += mi * (rhoj * rhoi_inv) * dvdotG;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdotG * rhoj_inv;
  pj->force.h_dt -= mi * dvdotG * rhoi_inv;

  const float v_sig = visc_signal_velocity;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

  if ((pi->is_h_max) || (pj->is_h_max)) {
    /* Do not add diffusion or normalising terms if either particle has h=h_max
     */
    return;
  }

  /* Calculate some quantities */
  const float mean_rho = 0.5f * (rhoi + rhoj);
  const float mean_rho_inv = 1.f / mean_rho;
  const float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
  const float mod_G = sqrtf(G_mean[0] * G_mean[0] + G_mean[1] * G_mean[1] +
                            G_mean[2] * G_mean[2]);
  const float v_sig_norm = sqrtf((pi->v[0] - pj->v[0]) * (pi->v[0] - pj->v[0]) +
                                 (pi->v[1] - pj->v[1]) * (pi->v[1] - pj->v[1]) +
                                 (pi->v[2] - pj->v[2]) * (pi->v[2] - pj->v[2]));

  /* Add normalising term to density evolution (Sandnes+2025 Eqn. 51) */
  const float alpha_norm = const_remix_norm_alpha;
  float drho_dt_norm_and_difn_i = alpha_norm * mj * v_sig_norm *
                                  pi->force.vac_switch *
                                  (pi->m0 * rhoi - rhoi) * mod_G * mean_rho_inv;
  float drho_dt_norm_and_difn_j = alpha_norm * mi * v_sig_norm *
                                  pj->force.vac_switch *
                                  (pj->m0 * rhoj - rhoj) * mod_G * mean_rho_inv;

  /* Only include diffusion for same-material particle pair */
  if (pi->mat_id == pj->mat_id) {
    /* Diffusion parameters */
    const float a_difn_rho = const_remix_difn_a_rho;
    const float b_difn_rho = const_remix_difn_b_rho;
    const float a_difn_u = const_remix_difn_a_u;
    const float b_difn_u = const_remix_difn_b_u;

    float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
    hydro_set_u_rho_difn(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi, pj,
                         dx);
    const float v_sig_difn = difn_signal_velocity;

    /* Calculate artificial diffusion of internal energy (Sandnes+2025 Eqn. 42)
     */
    const float du_dt_difn_i = -(a_difn_u + b_difn_u * mean_balsara) * mj *
                               v_sig_difn * (utilde_i - utilde_j) * mod_G *
                               mean_rho_inv;
    const float du_dt_difn_j = -(a_difn_u + b_difn_u * mean_balsara) * mi *
                               v_sig_difn * (utilde_j - utilde_i) * mod_G *
                               mean_rho_inv;

    /* Add artificial diffusion to evolution of internal energy */
    pi->u_dt += du_dt_difn_i;
    pj->u_dt += du_dt_difn_j;

    /* Calculate artificial diffusion of density (Sandnes+2025 Eqn. 43) */
    drho_dt_norm_and_difn_i += -(a_difn_rho + b_difn_rho * mean_balsara) * mj *
                               (rhoi * rhoj_inv) * v_sig_difn *
                               (rhotilde_i - rhotilde_j) * mod_G * mean_rho_inv;
    drho_dt_norm_and_difn_j += -(a_difn_rho + b_difn_rho * mean_balsara) * mi *
                               (rhoj * rhoi_inv) * v_sig_difn *
                               (rhotilde_j - rhotilde_i) * mod_G * mean_rho_inv;
  }

  /* Add normalising term and artificial diffusion to evolution of density */
  pi->drho_dt += drho_dt_norm_and_difn_i;
  pj->drho_dt += drho_dt_norm_and_difn_j;
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float rhoi_inv = 1.f / rhoi;
  const float rhoj_inv = 1.f / rhoj;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* Linear-order reproducing kernel gradient term (Sandnes+2025 Eqn. 28) */
  float Gj[3], Gi[3], G_mean[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);

  /* Antisymmetric kernel grad term for conservation of momentum and energy */
  G_mean[0] = 0.5f * (Gi[0] - Gj[0]);
  G_mean[1] = 0.5f * (Gi[1] - Gj[1]);
  G_mean[2] = 0.5f * (Gi[2] - Gj[2]);

  /* Viscous pressures (Sandnes+2025 Eqn. 41) */
  float Qi, Qj;
  float visc_signal_velocity, difn_signal_velocity;
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &difn_signal_velocity, pi,
                  pj, dx);

  /* Pressure terms to be used in evolution equations */
  const float P_i_term = pressurei * rhoi_inv * rhoj_inv;
  const float P_j_term = pressurej * rhoi_inv * rhoj_inv;
  const float Q_i_term = Qi * rhoi_inv * rhoj_inv;
  const float Q_j_term = Qj * rhoi_inv * rhoj_inv;

  /* Use the force Luke! */
  pi->a_hydro[0] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[0];
  pi->a_hydro[1] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[1];
  pi->a_hydro[2] -=
      mj * (P_i_term + P_j_term + Q_i_term + Q_j_term) * G_mean[2];

  /* v_ij dot kernel gradient term */
  const float dvdotG = (pi->v[0] - pj->v[0]) * G_mean[0] +
                       (pi->v[1] - pj->v[1]) * G_mean[1] +
                       (pi->v[2] - pj->v[2]) * G_mean[2];

  /* Internal energy time derivative */
  pi->u_dt += mj * (P_i_term + Q_i_term) * dvdotG;

  /* Density time derivative */
  pi->drho_dt += mj * (rhoi * rhoj_inv) * dvdotG;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdotG * rhoj_inv;

  const float v_sig = visc_signal_velocity;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

  if ((pi->is_h_max) || (pj->is_h_max)) {
    /* Do not add diffusion or normalising terms if either particle has h=h_max
     */
    return;
  }

  /* Calculate some quantities */
  const float mean_rho = 0.5f * (rhoi + rhoj);
  const float mean_rho_inv = 1.f / mean_rho;
  const float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
  const float mod_G = sqrtf(G_mean[0] * G_mean[0] + G_mean[1] * G_mean[1] +
                            G_mean[2] * G_mean[2]);

  const float v_sig_norm = sqrtf((pi->v[0] - pj->v[0]) * (pi->v[0] - pj->v[0]) +
                                 (pi->v[1] - pj->v[1]) * (pi->v[1] - pj->v[1]) +
                                 (pi->v[2] - pj->v[2]) * (pi->v[2] - pj->v[2]));

  /* Add normalising term to density evolution (Sandnes+2025 Eqn. 51) */
  const float alpha_norm = const_remix_norm_alpha;
  float drho_dt_norm_and_difn_i = alpha_norm * mj * v_sig_norm *
                                  pi->force.vac_switch *
                                  (pi->m0 * rhoi - rhoi) * mod_G * mean_rho_inv;

  /* Only include diffusion for same-material particle pair */
  if (pi->mat_id == pj->mat_id) {
    /* Diffusion parameters */
    const float a_difn_rho = const_remix_difn_a_rho;
    const float b_difn_rho = const_remix_difn_b_rho;
    const float a_difn_u = const_remix_difn_a_u;
    const float b_difn_u = const_remix_difn_b_u;

    float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
    hydro_set_u_rho_difn(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi, pj,
                         dx);
    const float v_sig_difn = difn_signal_velocity;

    /* Calculate artificial diffusion of internal energy (Sandnes+2025 Eqn. 42)
     */
    const float du_dt_difn_i = -(a_difn_u + b_difn_u * mean_balsara) * mj *
                               v_sig_difn * (utilde_i - utilde_j) * mod_G *
                               mean_rho_inv;

    /* Add artificial diffusion to evolution of internal energy */
    pi->u_dt += du_dt_difn_i;

    /* Calculate artificial diffusion of density (Sandnes+2025 Eqn. 43) */
    drho_dt_norm_and_difn_i += -(a_difn_rho + b_difn_rho * mean_balsara) * mj *
                               (rhoi * rhoj_inv) * v_sig_difn *
                               (rhotilde_i - rhotilde_j) * mod_G * mean_rho_inv;
  }

  /* Add normalising term and artificial diffusion to evolution of density */
  pi->drho_dt += drho_dt_norm_and_difn_i;
}

#endif /* SWIFT_REMIX_HYDRO_IACT_H */
