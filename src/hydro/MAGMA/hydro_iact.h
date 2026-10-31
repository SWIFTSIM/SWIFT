/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2022 Zhen Xiang (z.xiang@umail.leidenuniv.nl) &
 *                    Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MAGMA_HYDRO_IACT_H
#define SWIFT_MAGMA_HYDRO_IACT_H
//#define MAGMA_USE_SPHENIX_DIFFUSION_PARAMETER
//#define MAGMA_USE_FIRST_ORDER

/**
 * @file MAGMA/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (interaction routines)
 */

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
  float dv[3], curlvr[3];

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;
  pj->viscosity.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Negative because of the change in sign of dx & dv. */
  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];

  /* Compute the D correction matrix */

  for (int i = 0; i < 3; i++) {
    pi->magma.aux_u[i] += pj->mass * (pj->u - pi->u) * wi_dr * dx[i] * r_inv;
    pj->magma.aux_u[i] += pi->mass * (pj->u - pi->u) * wj_dr * dx[i] * r_inv;
    for (int j = 0; j < 3; j++) {
      /* Eq20 without inversing matrix. */
      pi->magma.d_matrix[i][j] -= pj->mass * dx[i] * wi_dr * dx[j] * r_inv;
      pj->magma.d_matrix[i][j] -= pi->mass * dx[i] * wj_dr * dx[j] * r_inv;
      /* Eq19 without multiplying matrix D. */
      pi->magma.aux_v[i][j] +=
          pj->mass * (pj->v[i] - pi->v[i]) * wi_dr * dx[j] * r_inv;
      pj->magma.aux_v[i][j] +=
          pi->mass * (pj->v[i] - pi->v[i]) * wj_dr * dx[j] * r_inv;
    }
  }

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pj->n_density += wj;
  pi->N_density++;
  pj->N_density++;
#endif
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
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float hid_inv = pow_dimension_plus_one(h_inv); /* 1/h^(d+1) */
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->viscosity.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  /* Compute the D correction matrix. */
  for (int i = 0; i < 3; i++) {
    pi->magma.aux_u[i] += pj->mass * (pj->u - pi->u) * wi_dr * dx[i] * r_inv;
    for (int j = 0; j < 3; j++) {
      /* Eq20 without inversing matrix. */
      pi->magma.d_matrix[i][j] -= pj->mass * dx[i] * wi_dr * dx[j] * r_inv;
      /* Eq19 without multiplying matrix D. */
      pi->magma.aux_v[i][j] +=
          pj->mass * (pj->v[i] - pi->v[i]) * wi_dr * dx[j] * r_inv;
    }
  }

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pi->N_density++;
#endif
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

  const float delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;
  pj->diffusion.laplace_u -= pi->mass * delta_u_factor * wj_dx / pi->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const float alpha_i = pi->viscosity.alpha;
  const float alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = max(pi->force.alpha_visc_max_ngb, alpha_j);
  pj->force.alpha_visc_max_ngb = max(pj->force.alpha_visc_max_ngb, alpha_i);

  /* Calculate the correction matrix. */
  for (int i = 0; i < 3; i++) {
    pi->magma.fder_u[i] -= pj->mass * (pj->u - pi->u) * dx[i] * wi / pj->rho;
    pj->magma.fder_u[i] -= pi->mass * (pj->u - pi->u) * dx[i] * wj / pi->rho;
    for (int j = 0; j < 3; j++) {
      /* Eq6 without inversing matrix. */
      pi->magma.c_matrix[i][j] += pj->mass * dx[i] * dx[j] * wi / pj->rho;
      pj->magma.c_matrix[i][j] += pi->mass * dx[i] * dx[j] * wj / pi->rho;

      pi->magma.sder_u[i][j] -= pj->mass *
                                (pj->magma.aux_u[i] - pi->magma.aux_u[i]) *
                                dx[j] * wi / pj->rho;
      pj->magma.sder_u[i][j] -= pi->mass *
                                (pj->magma.aux_u[i] - pi->magma.aux_u[i]) *
                                dx[j] * wj / pi->rho;

      /* Eq18 without multiplying matrix C. */
      pi->magma.fder_v[i][j] -=
          pj->mass * (pj->v[i] - pi->v[i]) * dx[j] * wi / pj->rho;
      pj->magma.fder_v[i][j] -=
          pi->mass * (pj->v[i] - pi->v[i]) * dx[j] * wj / pi->rho;
      for (int k = 0; k < 3; k++) {
        pi->magma.sder_v[i][j][k] -=
            pj->mass * (pj->magma.aux_v[i][j] - pi->magma.aux_v[i][j]) * dx[k] *
            wi / pj->rho;
        pj->magma.sder_v[i][j][k] -=
            pi->mass * (pj->magma.aux_v[i][j] - pi->magma.aux_v[i][j]) * dx[k] *
            wj / pi->rho;
      }
    }
  }

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_gradient += wi;
  pj->n_gradient += wj;
  pi->N_gradient++;
  pj->N_gradient++;
#endif
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

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  const float delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const float alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = max(pi->force.alpha_visc_max_ngb, alpha_j);

  /* Calculate the correction matrix. */
  for (int i = 0; i < 3; i++) {
    pi->magma.fder_u[i] -= pj->mass * (pj->u - pi->u) * dx[i] * wi / pj->rho;
    for (int j = 0; j < 3; j++) {
      /* Eq6 without inversing matrix. */
      pi->magma.c_matrix[i][j] += pj->mass * dx[i] * dx[j] * wi / pj->rho;
      pi->magma.sder_u[i][j] -= pj->mass *
                                (pj->magma.aux_u[i] - pi->magma.aux_u[i]) *
                                dx[j] * wi / pj->rho;
      /* Eq18 without multiplying matrix C. */
      pi->magma.fder_v[i][j] -=
          pj->mass * (pj->v[i] - pi->v[i]) * dx[j] * wi / pj->rho;
      for (int k = 0; k < 3; k++) {
        pi->magma.sder_v[i][j][k] -=
            pj->mass * (pj->magma.aux_v[i][j] - pi->magma.aux_v[i][j]) * dx[k] *
            wi / pj->rho;
      }
    }
  }

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_gradient += wi;
  pi->N_gradient++;
#endif
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
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Get the slope constant h(Eq23). */
  const float Hi = r / hi;
  const float Hj = r / hj;
  const float hij = min(Hi, Hj);
  const float h_crit = 0.5f;

  /* Get the slope constant A(Eq22). */
  float A_i = 0.f, A_j = 0.f, Av_i = 0.f, Av_j = 0.f;
  float A_ij, A_ji, Av_ij, Av_ji;

  for (int i = 0; i < 3; i++) {
    A_i += pi->magma.fder_u[i] * dx[i];
    A_j += pj->magma.fder_u[i] * dx[i];
    for (int j = 0; j < 3; j++) {
      Av_i += pi->magma.fder_v[i][j] * dx[i] * dx[j];
      Av_j += pj->magma.fder_v[i][j] * dx[i] * dx[j];
    }
  }

  if (A_i == 0.f && A_j == 0.f) { /* For smooth internal energy field, we turn
                                     off diffusion term*/
    A_ij = 1.f;
    A_ji = 1.f;
  } else if ((A_i == 0.f && A_j != 0.f) || (A_j == 0.f && A_i != 0.f) ||
             (A_i == -A_j)) { /* For extreme values, we add diffusion term*/
    A_ij = 0.f;
    A_ji = 0.f;
  } else {
    A_ij = A_i / A_j;
    A_ji = A_j / A_i;
  }

  if (Av_i == 0.f &&
      Av_j == 0.f) { /* For smooth velocity field, we turn off viscosity term*/
    Av_ij = 1.f;
    Av_ji = 1.f;
  } else if ((Av_i == 0.f && Av_j != 0.f) || (Av_j == 0.f && Av_i != 0.f) ||
             (Av_i == -Av_j)) { /* For extreme values, we add viscosity term*/
    Av_ij = 0.f;
    Av_ji = 0.f;
  } else {
    Av_ij = Av_i / Av_j;
    Av_ji = Av_j / Av_i;
  }

  /* Compute slope limiter(Eq21). */
  const float Ai_min = 4.f * A_ij / ((1.f + A_ij) * (1.f + A_ij));
  const float Aj_min = 4.f * A_ji / ((1.f + A_ji) * (1.f + A_ji));
  const float Fi_min = min(1.f, Ai_min);
  const float Fj_min = min(1.f, Aj_min);
  const float Fi_max = max(0.f, Fi_min);
  const float Fj_max = max(0.f, Fj_min);

  const float Avi_min = 4.f * Av_ij / ((1.f + Av_ij) * (1.f + Av_ij));
  const float Avj_min = 4.f * Av_ji / ((1.f + Av_ji) * (1.f + Av_ji));
  const float Fvi_min = min(1.f, Avi_min);
  const float Fvj_min = min(1.f, Avj_min);
  const float Fvi_max = max(0.f, Fvi_min);
  const float Fvj_max = max(0.f, Fvj_min);

  float F_ij, F_ji, Fv_ij, Fv_ji;

  if (hij > h_crit) {
    F_ij = Fi_max;
    F_ji = Fj_max;
    Fv_ij = Fvi_max;
    Fv_ji = Fvj_max;
  } else {
    const float e_con = -(hij - h_crit) * (hij - h_crit) * 25.f;
    const float F_exp = expf(e_con);
    F_ij = Fi_max * F_exp;
    F_ji = Fj_max * F_exp;
    Fv_ij = Fvi_max * F_exp;
    Fv_ji = Fvj_max * F_exp;
  }

  /* Compute reconstructed u. */
#if defined(MAGMA_USE_FIRST_ORDER)
  const float para_vis = 0.f;
#else
  const float para_vis = 0.5f;
#endif

  float ui_mid, ui_fder = 0.f, ui_sder = 0.f;
  float uj_mid, uj_fder = 0.f, uj_sder = 0.f;

  for (int i = 0; i < 3; i++) {
    ui_fder -= pi->magma.fder_u[i] * 0.5f * dx[i];
    uj_fder += pj->magma.fder_u[i] * 0.5f * dx[i];
    for (int j = 0; j < 3; j++) {
      ui_sder += pi->magma.sder_u[i][j] * 0.25f * dx[i] * dx[j];
      uj_sder += pj->magma.sder_u[i][j] * 0.25f * dx[i] * dx[j];
    }
  }

  ui_mid = pi->u + F_ij * (ui_fder + para_vis * ui_sder);
  uj_mid = pj->u + F_ji * (uj_fder + para_vis * uj_sder);

  /* Compute reconstructed v. */
  float vi_mid[3], vi_fder[3], vi_sder[3];
  float vj_mid[3], vj_fder[3], vj_sder[3];

  for (int i = 0; i < 3; i++) {
    vi_fder[i] = 0.f;
    vj_fder[i] = 0.f;
    vi_sder[i] = 0.f;
    vj_sder[i] = 0.f;
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      vi_fder[i] -= pi->magma.fder_v[i][j] * 0.5f * dx[j];
      vj_fder[i] += pj->magma.fder_v[i][j] * 0.5f * dx[j];
      for (int k = 0; k < 3; k++) {
        vi_sder[i] += pi->magma.sder_v[i][j][k] * 0.25f * dx[j] * dx[k];
        vj_sder[i] += pj->magma.sder_v[i][j][k] * 0.25f * dx[j] * dx[k];
      }
    }
  }

  for (int i = 0; i < 3; i++) { /* Eq17 */
    vi_mid[i] = pi->v[i] + Fv_ij * (vi_fder[i] + para_vis * vi_sder[i]);
    vj_mid[i] = pj->v[i] + Fv_ji * (vj_fder[i] + para_vis * vj_sder[i]);
  }

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Get de-dimensionalised separations(Eq16). */
  float h_i[3], h_j[3];

  for (int i = 0; i < 3; i++) {
    h_i[i] = dx[i] / hi;
    h_j[i] = -dx[i] / hj;
  }
  const float h_i2 = h_i[0] * h_i[0] + h_i[1] * h_i[1] + h_i[2] * h_i[2];
  const float h_j2 = h_j[0] * h_j[0] + h_j[1] * h_j[1] + h_j[2] * h_j[2];

  /* Compute new dv dot r. */
  const float dvdh_i = (vi_mid[0] - vj_mid[0]) * h_i[0] +
                       (vi_mid[1] - vj_mid[1]) * h_i[1] +
                       (vi_mid[2] - vj_mid[2]) * h_i[2] + a2_Hubble * r2 / hi;
  const float dvdh_j = (vj_mid[0] - vi_mid[0]) * h_j[0] +
                       (vj_mid[1] - vi_mid[1]) * h_j[1] +
                       (vj_mid[2] - vi_mid[2]) * h_j[2] + a2_Hubble * r2 / hj;
  const float dvdh_i2 = dvdh_i / (h_i2 + 0.01f);
  const float dvdh_j2 = dvdh_j / (h_j2 + 0.01f);

  /* compute velocity jump(Eq15). */
  const float mu_i = min(0.f, dvdh_i2);
  const float mu_j = min(0.f, dvdh_j2);

  /* compute viscosity pressure term(Eq14). */
  const float Q_i = fac_mu * rhoi *
                    (-pi->force.soundspeed * mu_i + 2.f * fac_mu * mu_i * mu_i);
  const float Q_j = fac_mu * rhoj *
                    (-pj->force.soundspeed * mu_j + 2.f * fac_mu * mu_j * mu_j);

  /* New gradient functions*/
  float c_matrix_i[3][3], c_matrix_j[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      c_matrix_i[i][j] = pi->magma.c_matrix[i][j];
      c_matrix_j[i][j] = pj->magma.c_matrix[i][j];
    }
  }

  float g_i[3], g_j[3], g_ij[3]; /*Eq4 and Eq5 for gradient functions G. */

  g_i[0] = -(c_matrix_i[0][0] * dx[0] + c_matrix_i[0][1] * dx[1] +
             c_matrix_i[0][2] * dx[2]) *
           wi * hi_inv_dim;
  g_i[1] = -(c_matrix_i[1][0] * dx[0] + c_matrix_i[1][1] * dx[1] +
             c_matrix_i[1][2] * dx[2]) *
           wi * hi_inv_dim;
  g_i[2] = -(c_matrix_i[2][0] * dx[0] + c_matrix_i[2][1] * dx[1] +
             c_matrix_i[2][2] * dx[2]) *
           wi * hi_inv_dim;

  g_j[0] = -(c_matrix_j[0][0] * dx[0] + c_matrix_j[0][1] * dx[1] +
             c_matrix_j[0][2] * dx[2]) *
           wj * hj_inv_dim;
  g_j[1] = -(c_matrix_j[1][0] * dx[0] + c_matrix_j[1][1] * dx[1] +
             c_matrix_j[1][2] * dx[2]) *
           wj * hj_inv_dim;
  g_j[2] = -(c_matrix_j[2][0] * dx[0] + c_matrix_j[2][1] * dx[1] +
             c_matrix_j[2][2] * dx[2]) *
           wj * hj_inv_dim;

  g_ij[0] = 0.5f * (g_i[0] + g_j[0]);
  g_ij[1] = 0.5f * (g_i[1] + g_j[1]);
  g_ij[2] = 0.5f * (g_i[2] + g_j[2]);

  /* Construct the density term */
  const float rho_ij = rhoi + rhoj;

  /* Compute gradient terms */
  const float P_over_rho_ij =
      (pressurei + Q_i + pressurej + Q_j) / (rhoi * rhoj);
  const float P_over_rho_i = (pressurei + Q_i) / (rhoi * rhoj);
  const float P_over_rho_j = (pressurej + Q_j) / (rhoi * rhoj);

  /* Use the force Luke ! */ /* Eq10 */
  pi->a_hydro[0] -= mj * P_over_rho_ij * g_ij[0];
  pi->a_hydro[1] -= mj * P_over_rho_ij * g_ij[1];
  pi->a_hydro[2] -= mj * P_over_rho_ij * g_ij[2];

  pj->a_hydro[0] += mi * P_over_rho_ij * g_ij[0];
  pj->a_hydro[1] += mi * P_over_rho_ij * g_ij[1];
  pj->a_hydro[2] += mi * P_over_rho_ij * g_ij[2];

  /* Get the v*g term*/
  const float dvg_i = (pi->v[0] - pj->v[0]) * g_ij[0] +
                      (pi->v[1] - pj->v[1]) * g_ij[1] +
                      (pi->v[2] - pj->v[2]) * g_ij[2];
  const float dvg_j = (pi->v[0] - pj->v[0]) * g_ij[0] +
                      (pi->v[1] - pj->v[1]) * g_ij[1] +
                      (pi->v[2] - pj->v[2]) * g_ij[2];

  /* Get the time derivative for u(Eq11). */
  const float sph_du_term_i = P_over_rho_i * dvg_i;
  const float sph_du_term_j = P_over_rho_j * dvg_j;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */

#if defined(MAGMA_USE_SPHENIX_DIFFUSION_PARAMETER)
  const float alpha_diff_sphenix =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const float alpha_diff = max(0.f, alpha_diff_sphenix);
#else
  const float alpha_diff = 0.05f;
#endif

  /* Two conductivity signal velocities(Eq26). */
  const float v_diffn = sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij);
  const float v_diffg =
      (sqrtf((vi_mid[0] - vj_mid[0]) * (vi_mid[0] - vj_mid[0]) +
             (vi_mid[1] - vj_mid[1]) * (vi_mid[1] - vj_mid[1]) +
             (vi_mid[2] - vj_mid[2]) * (vi_mid[2] - vj_mid[2])) +
       a2_Hubble * r) *
      fac_mu;
  const float iG = 0.5f;
  const float v_diff = (1.f - iG) * v_diffn + iG * v_diffg; /* Eq25. */
  const float aver_g = 0.5f * sqrtf((g_i[0] + g_j[0]) * (g_i[0] + g_j[0]) +
                                    (g_i[1] + g_j[1]) * (g_i[1] + g_j[1]) +
                                    (g_i[2] + g_j[2]) * (g_i[2] + g_j[2]));
  const float diff_du_term = -2.f * alpha_diff * v_diff * (ui_mid - uj_mid) *
                             aver_g / rho_ij; /* Eq24. */

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + diff_du_term;
  const float du_dt_j = sph_du_term_j - diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

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
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* Get the slope constant h(Eq23). */
  const float Hi = r / hi;
  const float Hj = r / hj;
  const float hij = min(Hi, Hj);
  const float h_crit = 0.5f;

  /* Get the slope constant A(Eq22). */
  float A_i = 0.f, A_j = 0.f, Av_i = 0.f, Av_j = 0.f;
  float A_ij, A_ji, Av_ij, Av_ji;

  for (int i = 0; i < 3; i++) {
    A_i += pi->magma.fder_u[i] * dx[i];
    A_j += pj->magma.fder_u[i] * dx[i];
    for (int j = 0; j < 3; j++) {
      Av_i += pi->magma.fder_v[i][j] * dx[i] * dx[j];
      Av_j += pj->magma.fder_v[i][j] * dx[i] * dx[j];
    }
  }

  if (A_i == 0.f && A_j == 0.f) { /* For smooth internal energy field, we turn
                                     off diffusion term*/
    A_ij = 1.f;
    A_ji = 1.f;
  } else if ((A_i == 0.f && A_j != 0.f) || (A_j == 0.f && A_i != 0.f) ||
             (A_i == -A_j)) { /* For extreme values, we add diffusion term*/
    A_ij = 0.f;
    A_ji = 0.f;
  } else {
    A_ij = A_i / A_j;
    A_ji = A_j / A_i;
  }

  if (Av_i == 0.f &&
      Av_j == 0.f) { /* For smooth velocity field, we turn off viscosity term*/
    Av_ij = 1.f;
    Av_ji = 1.f;
  } else if ((Av_i == 0.f && Av_j != 0.f) || (Av_j == 0.f && Av_i != 0.f) ||
             (Av_i == -Av_j)) { /* For extreme values, we add viscosity term*/
    Av_ij = 0.f;
    Av_ji = 0.f;
  } else {
    Av_ij = Av_i / Av_j;
    Av_ji = Av_j / Av_i;
  }

  /* Compute slope limiter(Eq21). */
  const float Ai_min = 4.f * A_ij / ((1.f + A_ij) * (1.f + A_ij));
  const float Aj_min = 4.f * A_ji / ((1.f + A_ji) * (1.f + A_ji));
  const float Fi_min = min(1.f, Ai_min);
  const float Fj_min = min(1.f, Aj_min);
  const float Fi_max = max(0.f, Fi_min);
  const float Fj_max = max(0.f, Fj_min);

  const float Avi_min = 4.f * Av_ij / ((1.f + Av_ij) * (1.f + Av_ij));
  const float Avj_min = 4.f * Av_ji / ((1.f + Av_ji) * (1.f + Av_ji));
  const float Fvi_min = min(1.f, Avi_min);
  const float Fvj_min = min(1.f, Avj_min);
  const float Fvi_max = max(0.f, Fvi_min);
  const float Fvj_max = max(0.f, Fvj_min);

  float F_ij, F_ji, Fv_ij, Fv_ji;

  if (hij > h_crit) {
    F_ij = Fi_max;
    F_ji = Fj_max;
    Fv_ij = Fvi_max;
    Fv_ji = Fvj_max;
  } else {
    const float e_con = -(hij - h_crit) * (hij - h_crit) * 25.f;
    const float F_exp = expf(e_con);
    F_ij = Fi_max * F_exp;
    F_ji = Fj_max * F_exp;
    Fv_ij = Fvi_max * F_exp;
    Fv_ji = Fvj_max * F_exp;
  }

  /* Compute reconstructed u. */
#if defined(MAGMA_USE_FIRST_ORDER)
  const float para_vis = 0.f;
#else
  const float para_vis = 0.5f;
#endif

  float ui_mid, ui_fder = 0.f, ui_sder = 0.f;
  float uj_mid, uj_fder = 0.f, uj_sder = 0.f;

  for (int i = 0; i < 3; i++) {
    ui_fder -= pi->magma.fder_u[i] * 0.5f * dx[i];
    uj_fder += pj->magma.fder_u[i] * 0.5f * dx[i];
    for (int j = 0; j < 3; j++) {
      ui_sder += pi->magma.sder_u[i][j] * 0.25f * dx[i] * dx[j];
      uj_sder += pj->magma.sder_u[i][j] * 0.25f * dx[i] * dx[j];
    }
  }

  ui_mid = pi->u + F_ij * (ui_fder + para_vis * ui_sder);
  uj_mid = pj->u + F_ji * (uj_fder + para_vis * uj_sder);

  /* Compute reconstructed v. */
  float vi_mid[3], vi_fder[3], vi_sder[3];
  float vj_mid[3], vj_fder[3], vj_sder[3];

  for (int i = 0; i < 3; i++) {
    vi_fder[i] = 0.f;
    vj_fder[i] = 0.f;
    vi_sder[i] = 0.f;
    vj_sder[i] = 0.f;
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      vi_fder[i] -= pi->magma.fder_v[i][j] * 0.5f * dx[j];
      vj_fder[i] += pj->magma.fder_v[i][j] * 0.5f * dx[j];
      for (int k = 0; k < 3; k++) {
        vi_sder[i] += pi->magma.sder_v[i][j][k] * 0.25f * dx[j] * dx[k];
        vj_sder[i] += pj->magma.sder_v[i][j][k] * 0.25f * dx[j] * dx[k];
      }
    }
  }

  for (int i = 0; i < 3; i++) { /* Eq17 */
    vi_mid[i] = pi->v[i] + Fv_ij * (vi_fder[i] + para_vis * vi_sder[i]);
    vj_mid[i] = pj->v[i] + Fv_ji * (vj_fder[i] + para_vis * vj_sder[i]);
  }

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Get de-dimensionalised separations(Eq16). */
  float h_i[3], h_j[3];

  for (int i = 0; i < 3; i++) {
    h_i[i] = dx[i] / hi;
    h_j[i] = -dx[i] / hj;
  }
  const float h_i2 = h_i[0] * h_i[0] + h_i[1] * h_i[1] + h_i[2] * h_i[2];
  const float h_j2 = h_j[0] * h_j[0] + h_j[1] * h_j[1] + h_j[2] * h_j[2];

  /* Compute new dv dot r. */
  const float dvdh_i = (vi_mid[0] - vj_mid[0]) * h_i[0] +
                       (vi_mid[1] - vj_mid[1]) * h_i[1] +
                       (vi_mid[2] - vj_mid[2]) * h_i[2] + a2_Hubble * r2 / hi;
  const float dvdh_j = (vj_mid[0] - vi_mid[0]) * h_j[0] +
                       (vj_mid[1] - vi_mid[1]) * h_j[1] +
                       (vj_mid[2] - vi_mid[2]) * h_j[2] + a2_Hubble * r2 / hj;
  const float dvdh_i2 = dvdh_i / (h_i2 + 0.01f);
  const float dvdh_j2 = dvdh_j / (h_j2 + 0.01f);

  /* compute velocity jump(Eq15). */
  const float mu_i = min(0.f, dvdh_i2);
  const float mu_j = min(0.f, dvdh_j2);

  /* compute viscosity pressure term(Eq14). */
  const float Q_i = fac_mu * rhoi *
                    (-pi->force.soundspeed * mu_i + 2.f * fac_mu * mu_i * mu_i);
  const float Q_j = fac_mu * rhoj *
                    (-pj->force.soundspeed * mu_j + 2.f * fac_mu * mu_j * mu_j);

  /* New gradient functions*/
  float c_matrix_i[3][3], c_matrix_j[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      c_matrix_i[i][j] = pi->magma.c_matrix[i][j];
      c_matrix_j[i][j] = pj->magma.c_matrix[i][j];
    }
  }

  float g_i[3], g_j[3], g_ij[3]; /*Eq4 and Eq5 for gradient functions G. */

  g_i[0] = -(c_matrix_i[0][0] * dx[0] + c_matrix_i[0][1] * dx[1] +
             c_matrix_i[0][2] * dx[2]) *
           wi * hi_inv_dim;
  g_i[1] = -(c_matrix_i[1][0] * dx[0] + c_matrix_i[1][1] * dx[1] +
             c_matrix_i[1][2] * dx[2]) *
           wi * hi_inv_dim;
  g_i[2] = -(c_matrix_i[2][0] * dx[0] + c_matrix_i[2][1] * dx[1] +
             c_matrix_i[2][2] * dx[2]) *
           wi * hi_inv_dim;

  g_j[0] = -(c_matrix_j[0][0] * dx[0] + c_matrix_j[0][1] * dx[1] +
             c_matrix_j[0][2] * dx[2]) *
           wj * hj_inv_dim;
  g_j[1] = -(c_matrix_j[1][0] * dx[0] + c_matrix_j[1][1] * dx[1] +
             c_matrix_j[1][2] * dx[2]) *
           wj * hj_inv_dim;
  g_j[2] = -(c_matrix_j[2][0] * dx[0] + c_matrix_j[2][1] * dx[1] +
             c_matrix_j[2][2] * dx[2]) *
           wj * hj_inv_dim;

  g_ij[0] = 0.5f * (g_i[0] + g_j[0]);
  g_ij[1] = 0.5f * (g_i[1] + g_j[1]);
  g_ij[2] = 0.5f * (g_i[2] + g_j[2]);

  /* Construct the density term */
  const float rho_ij = rhoi + rhoj;

  /* Compute gradient terms */
  const float P_over_rho_ij =
      (pressurei + Q_i + pressurej + Q_j) / (rhoi * rhoj);
  const float P_over_rho_i = (pressurei + Q_i) / (rhoi * rhoj);

  /* Use the force Luke ! */ /* Eq10 */
  pi->a_hydro[0] -= mj * P_over_rho_ij * g_ij[0];
  pi->a_hydro[1] -= mj * P_over_rho_ij * g_ij[1];
  pi->a_hydro[2] -= mj * P_over_rho_ij * g_ij[2];

  /* Get the v*g term*/
  const float dvg_i = (pi->v[0] - pj->v[0]) * g_ij[0] +
                      (pi->v[1] - pj->v[1]) * g_ij[1] +
                      (pi->v[2] - pj->v[2]) * g_ij[2];

  /* Get the time derivative for u(Eq11). */
  const float sph_du_term_i = P_over_rho_i * dvg_i;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */

#if defined(MAGMA_USE_SPHENIX_DIFFUSION_PARAMETER)
  const float alpha_diff_sphenix =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const float alpha_diff = max(0.f, alpha_diff_sphenix);
#else
  const float alpha_diff = 0.05f;
#endif

  /* Two conductivity signal velocities(Eq26). */
  const float v_diffn = sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij);
  const float v_diffg =
      (sqrtf((vi_mid[0] - vj_mid[0]) * (vi_mid[0] - vj_mid[0]) +
             (vi_mid[1] - vj_mid[1]) * (vi_mid[1] - vj_mid[1]) +
             (vi_mid[2] - vj_mid[2]) * (vi_mid[2] - vj_mid[2])) +
       a2_Hubble * r) *
      fac_mu;
  const float iG = 0.5f;
  const float v_diff = (1.f - iG) * v_diffn + iG * v_diffg; /* Eq25. */
  const float aver_g = 0.5f * sqrtf((g_i[0] + g_j[0]) * (g_i[0] + g_j[0]) +
                                    (g_i[1] + g_j[1]) * (g_i[1] + g_j[1]) +
                                    (g_i[2] + g_j[2]) * (g_i[2] + g_j[2]));
  const float diff_du_term = -2.f * alpha_diff * v_diff * (ui_mid - uj_mid) *
                             aver_g / rho_ij; /* Eq24. */

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
}

#endif /* SWIFT_MAGMA_HYDRO_IACT_H */
