/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_MAGMA2_HYDRO_IACT_H
#define SWIFT_MAGMA2_HYDRO_IACT_H

/**
 * @file MAGMA2/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (interaction routines)
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {
  float wi, wj, wi_dx, wj_dx;

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

  adaptive_softening_add_correction_term(pi, ui, hi_inv, mj);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  adaptive_softening_add_correction_term(pj, uj, hj_inv, mi);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Smooth pressure gradient */
  pi->gradients.pressure[0] += faci * pj->u * dx[0];
  pi->gradients.pressure[1] += faci * pj->u * dx[1];
  pi->gradients.pressure[2] += faci * pj->u * dx[2];

  pj->gradients.pressure[0] -= facj * pi->u * dx[0];
  pj->gradients.pressure[1] -= facj * pi->u * dx[1];
  pj->gradients.pressure[2] -= facj * pi->u * dx[2];

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] +=
          mj * dv[i] * dx[k] * wi_dx * r_inv;
      pj->gradients.velocity_tensor_aux[i][k] +=
          mi * dv[i] * dx[k] * wj_dx * r_inv;

      pi->gradients.velocity_tensor_aux_norm[i][k] +=
          mj * dx[i] * dx[k] * wi_dx * r_inv;
      pj->gradients.velocity_tensor_aux_norm[i][k] +=
          mi * dx[i] * dx[k] * wj_dx * r_inv;
    }
  }

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */

  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
  pj->weighted_wcount += mi * r2 * wj_dx * r_inv;
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {
  float wi, wi_dx;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  /* Compute pressure gradient */
  pi->gradients.pressure[0] += faci * pj->u * dx[0];
  pi->gradients.pressure[1] += faci * pj->u * dx[1];
  pi->gradients.pressure[2] += faci * pj->u * dx[2];

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] +=
          mj * dv[i] * dx[k] * wi_dx * r_inv;

      pi->gradients.velocity_tensor_aux_norm[i][k] +=
          mj * dx[i] * dx[k] * wi_dx * r_inv;
    }
  }

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */
  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Get particle properties */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const float rhoi = hydro_get_comoving_density(pi);
  const float rhoj = hydro_get_comoving_density(pj);

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig);

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx, wj, wj_dx;

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  const float shock_ratio_j =
      pi->viscosity.tensor_norm > 0.f
          ? pi->viscosity.shock_indicator / pi->viscosity.tensor_norm
          : 0.f;

  pi->viscosity.shock_limiter += mj * shock_ratio_i * wi;
  pj->viscosity.shock_limiter += mi * shock_ratio_j * wj;

  /* Correction factors for kernel gradients */
  const float rho_inv_i = 1.f / rhoi;
  const float rho_inv_j = 1.f / rhoj;

  pi->weighted_neighbour_wcount += mj * r2 * wi_dx * rho_inv_j * r_inv;
  pj->weighted_neighbour_wcount += mi * r2 * wj_dx * rho_inv_i * r_inv;

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.C[k][i] += mj * rho_inv_j * dx[k] * dx[i] * wi;
      pj->gradients.C[k][i] += mi * rho_inv_i * dx[k] * dx[i] * wj;
    }
  }
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Get particle properties */
  const float mj = hydro_get_mass(pj);
  const float rhoj = hydro_get_comoving_density(pj);

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);

  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  pi->viscosity.shock_limiter += mj * shock_ratio_i * wi;

  /* Correction factors for kernel gradients */

  const float rho_inv_j = 1.f / rhoj;

  pi->weighted_neighbour_wcount += mj * r2 * wi_dx * rho_inv_j * r_inv;

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.C[k][i] += mj * rho_inv_j * dx[k] * dx[i] * wi;
    }
  }
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {
  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  double G_ij[3] = {0., 0., 0.};

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  char C_well_conditioned = 0;
  if (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned) {
    C_well_conditioned = 1;
    double G_i[3] = {0., 0., 0.};
    double G_j[3] = {0., 0., 0.};
    for (int k = 0; k < 3; k++) {
      for (int i = 0; i < 3; i++) {
        /* Note: Negative because dx is (pj-pi) in Rosswog 2020.
         * It is (pj-pi) for both particles. */
        G_i[k] -= pi->gradients.C[k][i] * dx[i] * wi * hi_inv_dim;
        G_j[k] -= pj->gradients.C[k][i] * dx[i] * wj * hj_inv_dim;
      }

      G_ij[k] = 0.5 * (G_i[k] + G_j[k]);
    }
  }

  /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
  const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                             (pi->v[1] - pj->v[1]) * G_ij[1] +
                             (pi->v[2] - pj->v[2]) * G_ij[2];

  /* Compute dv dot dr. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      omega_ij < 0.f
          ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                 mu_ij +
             const_viscosity_beta * mu_ij * mu_ij) /
                (0.5f * rho_ij)
          : 0.f;

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = visc;
  float sph_acc_term = (pressurei + pressurej) / (pi->rho * pj->rho);
  
  float sph_du_term_i = pressurei / (pi->rho * pj->rho);
  float sph_du_term_j = pressurej / (pi->rho * pj->rho);
  float visc_du_term_i = 0.5f * visc_acc_term;
  float visc_du_term_j = visc_du_term_i;

  if (C_well_conditioned) {
    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    visc_du_term_i *= dv_dot_G_ij + a2_Hubble + r2;
    visc_du_term_j *= dv_dot_G_ij + a2_Hubble + r2;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
    pj->force.h_dt -= mi * dv_dot_G_ij;
  }
  else {
    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient =
        0.5f * r_inv * (wi_dr * pi->force.f + wj_dr * pj->force.f);

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;

    sph_du_term_i *= dvdr * kernel_gradient;
    sph_du_term_j *= dvdr * kernel_gradient;
    visc_du_term_i *= dvdr_Hubble * kernel_gradient;
    visc_du_term_j = visc_du_term_i;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dvdr * r_inv * wi_dr;
    pj->force.h_dt -= mi * dvdr * r_inv * wj_dr;
  }

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  pj->a_hydro[0] += mi * acc * G_ij[0];
  pj->a_hydro[1] += mi * acc * G_ij[1];
  pj->a_hydro[2] += mi * acc * G_ij[2];

  /* Get the time derivative for u. */

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term_i;
  const float du_dt_j = sph_du_term_j + visc_du_term_j;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {
  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  double G_ij[3] = {0., 0., 0.};

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  char C_well_conditioned = 0;
  if (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned) {
    C_well_conditioned = 1;
    double G_i[3] = {0., 0., 0.};
    double G_j[3] = {0., 0., 0.};
    for (int k = 0; k < 3; k++) {
      for (int i = 0; i < 3; i++) {
        /* Note: Negative because dx is (pj-pi) in Rosswog 2020 */
        G_i[k] -= pi->gradients.C[k][i] * dx[i] * wi * hi_inv_dim;
        G_j[k] -= pj->gradients.C[k][i] * dx[i] * wj * hj_inv_dim;
      }

      G_ij[k] = 0.5 * (G_i[k] + G_j[k]);
    }
  }

  /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
  const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                             (pi->v[1] - pj->v[1]) * G_ij[1] +
                             (pi->v[2] - pj->v[2]) * G_ij[2];

  /* Compute dv dot dr. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Construct the full viscosity term */
  const float rho_ij = rhoi + rhoj;
  const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const float visc =
      omega_ij < 0.f
          ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                 mu_ij +
             const_viscosity_beta * mu_ij * mu_ij) /
                (0.5f * rho_ij)
          : 0.f;

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = visc;
  float sph_acc_term = (pressurei + pressurej) / (pi->rho * pj->rho);
  float sph_du_term_i = pressurei / (pi->rho * pj->rho);
  float visc_du_term_i = visc;

  if (C_well_conditioned) {
    sph_du_term_i *= dv_dot_G_ij;
    visc_du_term_i *= dv_dot_G_ij + a2_Hubble * r2;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
  }
  else {
    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient =
        0.5f * r_inv * (wi_dr * pi->force.f + wj_dr * pj->force.f);

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term_i *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dvdr * r_inv * wi_dr;
  }

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term_i;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

}

#endif /* SWIFT_MAGMA2_HYDRO_IACT_H */
