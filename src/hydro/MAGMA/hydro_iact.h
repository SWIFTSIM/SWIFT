/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file Minimal/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
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

  runner_iact_nonsym_density(r2, dx, hi, hj, pi, pj, a, H);
  const float dx_rev[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_density(r2, dx_rev, hj, hi, pj, pi, a, H);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * Nothing to do here in this scheme.
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
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  float w;
  kernel_eval(ui, &w);

  /* Get the mass and density. */
  const float mj = pj->mass;
  const float rhoj = pj->rho;

  const float common_term = w * mj / rhoj;

  /* The inverse of the C-matrix.
   * It's symmetric so recall we only store the 6 useful terms. */
  pi->gradient.c_matrix_inv.xx += common_term * dx[0] * dx[0];
  pi->gradient.c_matrix_inv.yy += common_term * dx[1] * dx[1];
  pi->gradient.c_matrix_inv.zz += common_term * dx[2] * dx[2];
  pi->gradient.c_matrix_inv.xy += common_term * dx[0] * dx[1];
  pi->gradient.c_matrix_inv.xz += common_term * dx[0] * dx[2];
  pi->gradient.c_matrix_inv.yz += common_term * dx[1] * dx[2];
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * Nothing to do here in this scheme.
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
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  runner_iact_nonsym_gradient(r2, dx, hi, hj, pi, pj, a, H);
  const float dx_rev[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_gradient(r2, dx_rev, hj, hi, pj, pi, a, H);
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

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension(hi_inv); /* 1/h^d */
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension(hj_inv); /* 1/h^d */
  const float uj = r * hj_inv;
  float wj;
  kernel_eval(uj, &wj);

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi);
  const float P_over_rho2_j = pressurej / (rhoj * rhoj);

  /* Velocity difference */
  const float v_ij[3] = {pi->v[0] - pj->v[0],  /* x */
                         pi->v[1] - pj->v[1],  /* y */
                         pi->v[2] - pj->v[2]}; /* z */

  /* Compute dv dot r. */
  const float dvdr = v_ij[0] * dx[0] + v_ij[1] * dx[1] + v_ij[2] * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute signal velocity */
  const float v_sig = signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Construct the full viscosity term */
  // const float rho_ij = 0.5f * (rhoi + rhoj);
  // const float visc = -0.25f * v_sig * 2.f * mu_ij / rho_ij;

  /* Convolve with the kernel */
  // const float visc_acc_term = 0.;
  // 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* Construct the gradient functions (eq. 4 and 5) */
  float G_i[3], G_j[3];
  sym_matrix_multiply_by_vector(G_i, &pi->force.c_matrix, dx);
  sym_matrix_multiply_by_vector(G_j, &pj->force.c_matrix, dx);

  /* Note we multiply by -1 as dx is (pi - pj) and not (pj - pi) */
  G_i[0] *= -wi * hid_inv;
  G_i[1] *= -wi * hid_inv;
  G_i[1] *= -wi * hid_inv;
  G_j[0] *= -wj * hjd_inv;
  G_j[1] *= -wj * hjd_inv;
  G_j[1] *= -wj * hjd_inv;

  /* Raw fluid acceleration (eq. 2) */
  pi->a_hydro[0] -= mj * (P_over_rho2_i * G_i[0] + P_over_rho2_j * G_j[0]);
  pi->a_hydro[0] -= mj * (P_over_rho2_i * G_i[1] + P_over_rho2_j * G_j[1]);
  pi->a_hydro[0] -= mj * (P_over_rho2_i * G_i[2] + P_over_rho2_j * G_j[2]);

  const float v_ij_dot_G_i =
      v_ij[0] * G_i[0] + v_ij[1] * G_i[1] + v_ij[2] * G_i[2];

  /* Raw change in internal energy (eq. 3) */
  pi->u_dt += P_over_rho2_i * mj * v_ij_dot_G_i;

  /* Get the time derivative for u. */
  // const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  // const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  // const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  //

  /* Get the time derivative for h. */
  // pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  /* TODO: Think about this one */

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
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

  runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj, a, H);
  const float dx_rev[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_force(r2, dx_rev, hj, hi, pj, pi, a, H);
}

#endif /* SWIFT_MAGMA_HYDRO_IACT_H */
