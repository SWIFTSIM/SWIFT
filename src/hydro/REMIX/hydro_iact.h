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
#ifndef SWIFT_PLANETARY_HYDRO_IACT_H
#define SWIFT_PLANETARY_HYDRO_IACT_H

/**
 * @file Planetary/hydro_iact.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2024)
 */

#include "adiabatic_index.h"
#include "const.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "hydro_visc_difn.h"
#include "hydro_strength.h"
#include "material_properties.h"
#include "math.h"
#include "minmax.h"

/**
 * @brief Calculates the stress tensor for force interaction. No strength if
 * either particle is fluid (or not using strength at all)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_pairwise_stress_tensors(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj, const float r,
    const float pressurei, const float pressurej) {

  // Set the default stress tensor with just the pressures, S = -P * I(3)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (i == j) {
            // Only include the pressure if it is positive (i.e. not in tension)
            pairwise_stress_tensor_i[i][i] = -max(pressurei, 0.f);
            pairwise_stress_tensor_j[i][i] = -max(pressurej, 0.f);
        } else {
            pairwise_stress_tensor_i[i][j] = 0.f;
            pairwise_stress_tensor_j[i][j] = 0.f;
        }
    }
  }

#ifdef MATERIAL_STRENGTH
  hydro_set_pairwise_stress_tensors_strength(pairwise_stress_tensor_i, pairwise_stress_tensor_j, pi, pj, r);
#endif /* MATERIAL_STRENGTH */
}

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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pj->n_density += wj;
  pi->N_density++;
  pj->N_density++;
#endif

  hydro_runner_iact_density_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_density_extra_viscosity(pi, pj, dx, wi, wj, wi_dx, wj_dx);
#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_density_extra_strength(pi, pj, dx, wi, wj, wi_dx, wj_dx);
#endif /* MATERIAL_STRENGTH */
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pi->N_density++;
#endif

  hydro_runner_iact_nonsym_density_extra_kernel(pi, pj, dx, wi, wi_dx);
  hydro_runner_iact_nonsym_density_extra_viscosity(pi, pj, dx, wi, wi_dx);
#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_nonsym_density_extra_strength(pi, pj, dx, wi, wi_dx);
#endif /* MATERIAL_STRENGTH */
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
  hydro_runner_iact_gradient_extra_viscosity(pi, pj, dx, wi, wj, wi_dx, wj_dx);
#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_gradient_extra_strength(pi, pj, dx, wi, wj, wi_dx, wj_dx);
#endif /* MATERIAL_STRENGTH */
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

  hydro_runner_iact_nonsym_gradient_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_nonsym_gradient_extra_viscosity(pi, pj, dx, wi, wi_dx);
#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_nonsym_gradient_extra_strength(pi, pj, dx, wi, wi_dx);
#endif /* MATERIAL_STRENGTH */
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
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

  // Linear-order reproducing kernel gradient term
  float Gj[3], Gi[3], G_mean[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);
  for (int i = 0; i < 3; i++) {
    G_mean[i] = 0.5f * (Gi[i] - Gj[i]);
  }

  // Convert the pressures into "stress" tensors for fluids, or set the stress
  // for solid particle pairs with strength
  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  hydro_set_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j,
                                    pi, pj, r, pressurei, pressurej);

  // Viscous pressures
  float Qi, Qj;
  float visc_signal_velocity, difn_signal_velocity;
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &difn_signal_velocity, pi,
                  pj, dx, a, H);

  // Incorporate kernel gradient term into viscous pressures and stress tensors
  float Q_term_i[3], Q_term_j[3];
  float stress_tensor_term_i[3], stress_tensor_term_j[3];
  for (int i = 0; i < 3; i++) {
      Q_term_i[i] = -Qi * G_mean[i];
      Q_term_j[i] = -Qj * G_mean[i];

      stress_tensor_term_i[i] = 0.f;
      stress_tensor_term_j[i] = 0.f;
      for (int j = 0; j < 3; j++) {
          stress_tensor_term_i[i] += pairwise_stress_tensor_i[i][j] * G_mean[j];
          stress_tensor_term_j[i] += pairwise_stress_tensor_j[i][j] * G_mean[j];
      }
  }

  /* Use the force Luke! */
  for (int i = 0; i < 3; i++) {
    pi->a_hydro[i] +=
        mj * (stress_tensor_term_i[i] + stress_tensor_term_j[i] + Q_term_i[i] + Q_term_j[i])
        / (rhoi * rhoj);
    pj->a_hydro[i] -=
        mi * (stress_tensor_term_i[i] + stress_tensor_term_j[i] + Q_term_i[i] + Q_term_j[i])
        / (rhoi * rhoj);
  }

  // Compute dv terms for u and h time derivatives
  float dv_dot_Q_term_i=0.f, dv_dot_Q_term_j=0.f;
  float dv_dot_stress_term_i=0.f, dv_dot_stress_term_j=0.f;
  float dv_dot_G=0.f;
  for (int i=0; i<3; i++) {
      dv_dot_Q_term_i += (pi->v[i] - pj->v[i]) * Q_term_i[i];
      dv_dot_Q_term_j += (pi->v[i] - pj->v[i]) * Q_term_j[i];

      dv_dot_stress_term_i += (pi->v[i] - pj->v[i]) * stress_tensor_term_i[i];
      dv_dot_stress_term_j += (pi->v[i] - pj->v[i]) * stress_tensor_term_j[i];

      dv_dot_G += (pi->v[i] - pj->v[i]) * G_mean[i];
  }

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = -(dv_dot_stress_term_i + dv_dot_Q_term_i) / (rhoi * rhoj);
  float du_dt_j = -(dv_dot_stress_term_j + dv_dot_Q_term_j) / (rhoi * rhoj);

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dv_dot_G / rhoj;
  pj->force.h_dt -= mi * dv_dot_G / rhoi;

  const float v_sig = visc_signal_velocity;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

  pi->drho_dt += mj * (rhoi / rhoj) * dv_dot_G;
  pj->drho_dt += mi * (rhoj / rhoi) * dv_dot_G;

  // Diffusion and kernel normalising term, if particles are not at h=h_max
  if ((!pi->is_h_max) && (!pj->is_h_max)) {
    float mean_rho = 0.5f * (rhoi + rhoj);
    float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
    float mod_G = sqrtf(G_mean[0] * G_mean[0] + G_mean[1] * G_mean[1] +
                        G_mean[2] * G_mean[2]);
    float v_sig_norm = sqrtf((pi->v[0] - pj->v[0]) * (pi->v[0] - pj->v[0]) +
                             (pi->v[1] - pj->v[1]) * (pi->v[1] - pj->v[1]) +
                             (pi->v[2] - pj->v[2]) * (pi->v[2] - pj->v[2]));

    // Kernel normalising term (zero for solid particles with strength)
    const float alpha_norm = diffusion_global.alpha_norm;
    float drho_dt_norm_and_difn_i;
    float drho_dt_norm_and_difn_j;
    if (pi->phase_state == mat_phase_state_fluid) {
      drho_dt_norm_and_difn_i =
        alpha_norm * mj * v_sig_norm * pi->force.vac_switch *
        (pi->m0 * pi->rho_evol - pi->rho_evol) * mod_G / mean_rho;
    } else {
        drho_dt_norm_and_difn_i = 0.f;
    }
    if (pj->phase_state == mat_phase_state_fluid) {
      drho_dt_norm_and_difn_j =
        alpha_norm * mi * v_sig_norm * pj->force.vac_switch *
        (pj->m0 * pj->rho_evol - pj->rho_evol) * mod_G / mean_rho;
    } else {
        drho_dt_norm_and_difn_j = 0.f;
    }

    // Diffusion for same materials
    if (pi->mat_id == pj->mat_id) {
      // Diffusion parameters
      const float a_difn_rho = diffusion_global.a_difn_rho;
      const float b_difn_rho = diffusion_global.b_difn_rho;
      const float a_difn_u = diffusion_global.a_difn_u;
      const float b_difn_u = diffusion_global.b_difn_u;

      // ...
      float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
      hydro_set_u_rho_difn(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi,
                           pj, dx, a, H);
      float v_sig_difn = difn_signal_velocity;
      float du_dt_difn_i = -(a_difn_u + b_difn_u * mean_balsara) * mj *
                           v_sig_difn * (utilde_i - utilde_j) * mod_G /
                           mean_rho;
      float du_dt_difn_j = -(a_difn_u + b_difn_u * mean_balsara) * mi *
                           v_sig_difn * (utilde_j - utilde_i) * mod_G /
                           mean_rho;

      // ...
      pi->u_dt += du_dt_difn_i;
      pj->u_dt += du_dt_difn_j;

      // ...
      drho_dt_norm_and_difn_i += -(a_difn_rho + b_difn_rho * mean_balsara) *
                                 mj * (rhoi / rhoj) * v_sig_difn *
                                 (rhotilde_i - rhotilde_j) * mod_G / mean_rho;
      drho_dt_norm_and_difn_j += -(a_difn_rho + b_difn_rho * mean_balsara) *
                                 mi * (rhoj / rhoi) * v_sig_difn *
                                 (rhotilde_j - rhotilde_i) * mod_G / mean_rho;
    }

    pi->drho_dt += drho_dt_norm_and_difn_i;
    pj->drho_dt += drho_dt_norm_and_difn_j;
  }

#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_force_extra_strength(pi, pj, dx, Gi, Gj);
#endif /* MATERIAL_STRENGTH */

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pj->n_force += wi + wj;
  pi->N_force++;
  pj->N_force++;
#endif
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
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

  // Linear-order reproducing kernel gradient term
  float Gj[3], Gi[3], G_mean[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);
  for (int i = 0; i < 3; i++) {
    G_mean[i] = 0.5f * (Gi[i] - Gj[i]);
  }

  // Convert the pressures into "stress" tensors for fluids, S = -P * I(3),
  // or set the stress for solid particle pairs with strength
  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  hydro_set_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j,
                                    pi, pj, r, pressurei, pressurej);

  // Viscous pressures
  float Qi, Qj;
  float visc_signal_velocity, difn_signal_velocity;
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &difn_signal_velocity, pi,
                  pj, dx, a, H);

  // Incorporate kernel gradient term into viscous pressures and stress tensors
  float Q_term_i[3], Q_term_j[3];
  float stress_tensor_term_i[3], stress_tensor_term_j[3];
  for (int i = 0; i < 3; i++) {
      Q_term_i[i] = -Qi * G_mean[i];
      Q_term_j[i] = -Qj * G_mean[i];

      stress_tensor_term_i[i] = 0.f;
      stress_tensor_term_j[i] = 0.f;
      for (int j = 0; j < 3; j++) {
          stress_tensor_term_i[i] += pairwise_stress_tensor_i[i][j] * G_mean[j];
          stress_tensor_term_j[i] += pairwise_stress_tensor_j[i][j] * G_mean[j];
      }
  }

  /* Use the force Luke! */
  for (int i = 0; i < 3; i++) {
    pi->a_hydro[i] +=
        mj * (stress_tensor_term_i[i] + stress_tensor_term_j[i] + Q_term_i[i] + Q_term_j[i])
        / (rhoi * rhoj);
  }

  // Compute dv terms for u and h time derivatives
  float dv_dot_Q_term_i=0.f;
  float dv_dot_stress_term_i=0.f;
  float dv_dot_G=0.f;
  for (int i=0; i<3; i++) {
      dv_dot_Q_term_i += (pi->v[i] - pj->v[i]) * Q_term_i[i];

      dv_dot_stress_term_i += (pi->v[i] - pj->v[i]) * stress_tensor_term_i[i];

      dv_dot_G += (pi->v[i] - pj->v[i]) * G_mean[i];
  }

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = -(dv_dot_stress_term_i + dv_dot_Q_term_i) / (rhoi * rhoj);

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dv_dot_G / rhoj;

  const float v_sig = visc_signal_velocity;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

  pi->drho_dt += mj * (rhoi / rhoj) * dv_dot_G;

  if ((!pi->is_h_max) && (!pj->is_h_max)) {
    float mean_rho = 0.5f * (rhoi + rhoj);
    float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
    float mod_G = sqrtf(G_mean[0] * G_mean[0] + G_mean[1] * G_mean[1] +
                        G_mean[2] * G_mean[2]);

    float v_sig_norm = sqrtf((pi->v[0] - pj->v[0]) * (pi->v[0] - pj->v[0]) +
                             (pi->v[1] - pj->v[1]) * (pi->v[1] - pj->v[1]) +
                             (pi->v[2] - pj->v[2]) * (pi->v[2] - pj->v[2]));

    const float alpha_norm = diffusion_global.alpha_norm;
    float drho_dt_norm_and_difn_i = 0.f;
    if (pi->phase_state == mat_phase_state_fluid) {
      drho_dt_norm_and_difn_i +=
        alpha_norm * mj * v_sig_norm * pi->force.vac_switch *
        (pi->m0 * pi->rho_evol - pi->rho_evol) * mod_G / mean_rho;
    }

    // Diffusion for same materials
    if (pi->mat_id == pj->mat_id) {
      // Diffusion parameters
      const float a_difn_rho = diffusion_global.a_difn_rho;
      const float b_difn_rho = diffusion_global.b_difn_rho;
      const float a_difn_u = diffusion_global.a_difn_u;
      const float b_difn_u = diffusion_global.b_difn_u;

      // ...
      float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
      hydro_set_u_rho_difn(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi,
                           pj, dx, a, H);
      float v_sig_difn = difn_signal_velocity;
      float du_dt_difn_i = -(a_difn_u + b_difn_u * mean_balsara) * mj *
                           v_sig_difn * (utilde_i - utilde_j) * mod_G /
                           mean_rho;

      // ...
      pi->u_dt += du_dt_difn_i;

      // ...
      drho_dt_norm_and_difn_i += -(a_difn_rho + b_difn_rho * mean_balsara) *
                                 mj * (rhoi / rhoj) * v_sig_difn *
                                 (rhotilde_i - rhotilde_j) * mod_G / mean_rho;
    }

    pi->drho_dt += drho_dt_norm_and_difn_i;
  }

#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_nonsym_force_extra_strength(pi, pj, dx, Gi);
#endif /* MATERIAL_STRENGTH */

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
}

#endif /* SWIFT_PLANETARY_HYDRO_IACT_H */
