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

  /* Velocity difference */
  const float vij[3] = {pj->v[0] - pi->v[0], pj->v[1] - pi->v[1],
                        pj->v[2] - pi->v[2]};

  /* Internal energy difference */
  const float uij = pj->u - pi->u;

  const float common_term = w * mj / rhoj;

  /* The inverse of the C-matrix. eq. 6
   * It's symmetric so recall we only store the 6 useful terms. */
  pi->gradient.c_matrix_inv.xx += common_term * dx[0] * dx[0];
  pi->gradient.c_matrix_inv.yy += common_term * dx[1] * dx[1];
  pi->gradient.c_matrix_inv.zz += common_term * dx[2] * dx[2];
  pi->gradient.c_matrix_inv.xy += common_term * dx[0] * dx[1];
  pi->gradient.c_matrix_inv.xz += common_term * dx[0] * dx[2];
  pi->gradient.c_matrix_inv.yz += common_term * dx[1] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pi->gradient.gradient_vx[0] -= common_term * vij[0] * dx[0];
  pi->gradient.gradient_vx[1] -= common_term * vij[0] * dx[1];
  pi->gradient.gradient_vx[2] -= common_term * vij[0] * dx[2];

  pi->gradient.gradient_vy[0] -= common_term * vij[1] * dx[0];
  pi->gradient.gradient_vy[1] -= common_term * vij[1] * dx[1];
  pi->gradient.gradient_vy[2] -= common_term * vij[1] * dx[2];

  pi->gradient.gradient_vz[0] -= common_term * vij[2] * dx[0];
  pi->gradient.gradient_vz[1] -= common_term * vij[2] * dx[1];
  pi->gradient.gradient_vz[2] -= common_term * vij[2] * dx[2];

  pi->gradient.gradient_u[0] -= common_term * uij * dx[0];
  pi->gradient.gradient_u[1] -= common_term * uij * dx[1];
  pi->gradient.gradient_u[2] -= common_term * uij * dx[2];
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
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension(hi_inv); /* 1/h^d */
  const float ui = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(ui, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension(hj_inv); /* 1/h^d */
  const float uj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(uj, &wj, &wj_dx);

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

  /* De-dimentionalised distances (eq. 16, recall dx = xi - xj)*/
  const float eta_i[3] = {dx[0] / hi, dx[1] / hi, dx[2] / hi};
  const float eta_j[3] = {-dx[0] / hj, -dx[1] / hj, -dx[2] / hj};

  /* Norms of the eta vectors (eq. 16) */
  const float eta_square_i =
      eta_i[0] * eta_i[0] + eta_i[1] * eta_i[1] + eta_i[2] * eta_i[2];
  const float eta_square_j =
      eta_j[0] * eta_j[0] + eta_j[1] * eta_j[1] + eta_j[2] * eta_j[2];

  /* Reconstructed velocities at the mid-point (before reconstruction) */
  float v_rec_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
  float v_rec_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

#ifndef USE_ZEROTH_ORDER_VELOCITIES

  /* Vectors from the particles to the mid-point */
  const float delta_i[3] = {0.5f * (pj->x[0] - pi->x[0]),
                            0.5f * (pj->x[1] - pi->x[1]),
                            0.5f * (pj->x[2] - pi->x[2])};
  const float delta_j[3] = {-delta_i[0], -delta_i[1], -delta_i[2]};

  /* Terms entering the limiter (eq. 23) */
  const float eta_ij = sqrtf(fminf(eta_square_i, eta_square_j));
  const float eta_crit = 0.5f;

  /* Van Leer limiter fraction (eq. 22) */
  const float A_ij_num = pi->gradient.gradient_vx[0] * dx[0] * dx[0] +
                         pi->gradient.gradient_vx[1] * dx[0] * dx[1] +
                         pi->gradient.gradient_vx[2] * dx[0] * dx[2] +
                         pi->gradient.gradient_vy[0] * dx[1] * dx[0] +
                         pi->gradient.gradient_vy[1] * dx[1] * dx[1] +
                         pi->gradient.gradient_vy[2] * dx[1] * dx[2] +
                         pi->gradient.gradient_vz[0] * dx[2] * dx[0] +
                         pi->gradient.gradient_vz[1] * dx[2] * dx[1] +
                         pi->gradient.gradient_vz[2] * dx[2] * dx[2];

  const float A_ij_den = pj->gradient.gradient_vx[0] * dx[0] * dx[0] +
                         pj->gradient.gradient_vx[1] * dx[0] * dx[1] +
                         pj->gradient.gradient_vx[2] * dx[0] * dx[2] +
                         pj->gradient.gradient_vy[0] * dx[1] * dx[0] +
                         pj->gradient.gradient_vy[1] * dx[1] * dx[1] +
                         pj->gradient.gradient_vy[2] * dx[1] * dx[2] +
                         pj->gradient.gradient_vz[0] * dx[2] * dx[0] +
                         pj->gradient.gradient_vz[1] * dx[2] * dx[1] +
                         pj->gradient.gradient_vz[2] * dx[2] * dx[2];

  const float A_ij = A_ij_num / A_ij_den;

  /* Slope limiter exponential term (eq. 21, right term) */
  const float exp_term =
      eta_ij < eta_crit
          ? expf(-25.f * (eta_ij - eta_crit) * (eta_ij - eta_crit))
          : 1.f;

  /* Van Leer limiter (eq. 21) */
  const float fraction = 4.f * A_ij / ((1.f + A_ij) * (1.f * A_ij));
  const float Phi_ij = fmaxf(0.f, fminf(1.f, fraction)) * exp_term;

  /* Mid-point reconstruction, first order (eq. 17) */
  v_rec_i[0] += Phi_ij * pi->gradient.gradient_vx[0] * delta_i[0];
  v_rec_i[0] += Phi_ij * pi->gradient.gradient_vx[1] * delta_i[1];
  v_rec_i[0] += Phi_ij * pi->gradient.gradient_vx[2] * delta_i[2];
  v_rec_i[1] += Phi_ij * pi->gradient.gradient_vy[0] * delta_i[0];
  v_rec_i[1] += Phi_ij * pi->gradient.gradient_vy[1] * delta_i[1];
  v_rec_i[1] += Phi_ij * pi->gradient.gradient_vy[2] * delta_i[2];
  v_rec_i[2] += Phi_ij * pi->gradient.gradient_vz[0] * delta_i[0];
  v_rec_i[2] += Phi_ij * pi->gradient.gradient_vz[1] * delta_i[1];
  v_rec_i[2] += Phi_ij * pi->gradient.gradient_vz[2] * delta_i[2];

  v_rec_j[0] += Phi_ij * pj->gradient.gradient_vx[0] * delta_j[0];
  v_rec_j[0] += Phi_ij * pj->gradient.gradient_vx[1] * delta_j[1];
  v_rec_j[0] += Phi_ij * pj->gradient.gradient_vx[2] * delta_j[2];
  v_rec_j[1] += Phi_ij * pj->gradient.gradient_vy[0] * delta_j[0];
  v_rec_j[1] += Phi_ij * pj->gradient.gradient_vy[1] * delta_j[1];
  v_rec_j[1] += Phi_ij * pj->gradient.gradient_vy[2] * delta_j[2];
  v_rec_j[2] += Phi_ij * pj->gradient.gradient_vz[0] * delta_j[0];
  v_rec_j[2] += Phi_ij * pj->gradient.gradient_vz[1] * delta_j[1];
  v_rec_j[2] += Phi_ij * pj->gradient.gradient_vz[2] * delta_j[2];

#endif

  /* Difference in velocity at the mid-point */
  const float v_rec_ij[3] = {v_rec_i[0] - v_rec_j[0], v_rec_i[1] - v_rec_j[1],
                             v_rec_i[2] - v_rec_j[2]};

  /* Normalised relative velocity (eq. 15) */
  const float vel_rel_i =
      eta_i[0] * v_rec_ij[0] + eta_i[1] * v_rec_ij[1] + eta_i[2] * v_rec_ij[2];
  const float vel_rel_j = eta_j[0] * -v_rec_ij[0] + eta_j[1] * -v_rec_ij[1] +
                          eta_j[2] * -v_rec_ij[2];

  /* Terms entering the viscosity (eq. 15) */
  const float eps_squared = const_viscosity_epsilon * const_viscosity_epsilon;
  const float mu_i = fminf(0.f, vel_rel_i / (eta_square_i + eps_squared));
  const float mu_j = fminf(0.f, vel_rel_j / (eta_square_j + eps_squared));

  /* Eq. 14 */
  const float Qi =
      1.f * rhoi *
      (-const_viscosity_alpha * ci * mu_i + const_viscosity_beta * mu_i * mu_i);
  const float Qj =
      1.f * rhoj *
      (-const_viscosity_alpha * cj * mu_j + const_viscosity_beta * mu_j * mu_j);

  /* Construct the gradient functions (eq. 4 and 5) */
  float G_i[3], G_j[3];
  sym_matrix_multiply_by_vector(G_i, &pi->force.c_matrix, dx);
  sym_matrix_multiply_by_vector(G_j, &pj->force.c_matrix, dx);

  /* Note we multiply by -1 as dx is (pi - pj) and not (pj - pi) */
  G_i[0] *= -wi * hid_inv;
  G_i[1] *= -wi * hid_inv;
  G_i[2] *= -wi * hid_inv;
  G_j[0] *= -wj * hjd_inv;
  G_j[1] *= -wj * hjd_inv;
  G_j[2] *= -wj * hjd_inv;

#ifdef USE_STANDARD_KERNEL_GRADIENTS
  const float wi_dr = hid_inv * hi_inv * wi_dx;
  const float wj_dr = hjd_inv * hj_inv * wj_dx;
  G_i[0] = wi_dr * r_inv * dx[0];
  G_i[1] = wi_dr * r_inv * dx[1];
  G_i[2] = wi_dr * r_inv * dx[2];
  G_j[0] = wj_dr * r_inv * dx[0];
  G_j[1] = wj_dr * r_inv * dx[1];
  G_j[2] = wj_dr * r_inv * dx[2];
#endif

#ifdef I_LOVE_GASOLINE

  const float G_ij[3] = {0.5f * (G_i[0] + G_j[0]), 0.5f * (G_i[1] + G_j[1]),
                         0.5f * (G_i[2] + G_j[2])};

  const float acc_term = (pressurei + pressurej + Qi + Qj) / (rhoi * rhoj);
  const float du_term = (pressurei + Qi) / (rhoi * rhoj);

  /* Raw fluid acceleration (eq. 10) */
  pi->a_hydro[0] -= mj * acc_term * G_ij[0];
  pi->a_hydro[1] -= mj * acc_term * G_ij[1];
  pi->a_hydro[2] -= mj * acc_term * G_ij[2];

  /* Equivalent of div v */
  const float v_ij_dot_G_ij =
      v_ij[0] * G_ij[0] + v_ij[1] * G_ij[1] + v_ij[2] * G_ij[2];

  /* Raw change in internal energy (eq. 11) */
  pi->u_dt += du_term * mj * v_ij_dot_G_ij;

#else

  /* Compute pressure terms */
  const float P_over_rho2_i = (pressurei + Qi) / (rhoi * rhoi);
  const float P_over_rho2_j = (pressurej + Qj) / (rhoj * rhoj);

  /* Raw fluid acceleration (eq. 2) */
  pi->a_hydro[0] -= mj * (P_over_rho2_i * G_i[0] + P_over_rho2_j * G_j[0]);
  pi->a_hydro[1] -= mj * (P_over_rho2_i * G_i[1] + P_over_rho2_j * G_j[1]);
  pi->a_hydro[2] -= mj * (P_over_rho2_i * G_i[2] + P_over_rho2_j * G_j[2]);

  /* Equivalent of div v */
  const float v_ij_dot_G_i =
      v_ij[0] * G_i[0] + v_ij[1] * G_i[1] + v_ij[2] * G_i[2];

  /* Raw change in internal energy (eq. 3) */
  pi->u_dt += P_over_rho2_i * mj * v_ij_dot_G_i;

#endif

  /* Reconstructed internal energies at the mid-point (before reconstruction) */
  float u_rec_i = pi->u;
  float u_rec_j = pj->u;

#ifndef USE_ZEROTH_ORDER_VELOCITIES

  /* Mid-point reconstruction, first order (eq. 17) */
  u_rec_i += pi->gradient.gradient_u[0] * delta_i[0];
  u_rec_i += pi->gradient.gradient_u[1] * delta_i[1];
  u_rec_i += pi->gradient.gradient_u[2] * delta_i[2];

  u_rec_j += pj->gradient.gradient_u[0] * delta_j[0];
  u_rec_j += pj->gradient.gradient_u[1] * delta_j[1];
  u_rec_j += pj->gradient.gradient_u[2] * delta_j[2];

#endif

  /* Difference in internal energy */
  const float delta_u = u_rec_i - u_rec_j;

  /* Norm of the G vectors */
  const float sum_G[3] = {G_i[0] + G_j[0], G_i[1] + G_j[1], G_i[2] + G_j[2]};
  const float norm_sum_G =
      sqrtf(sum_G[0] * sum_G[0] + sum_G[1] * sum_G[1] + sum_G[2] * sum_G[2]);

  /* Diffusion signal velocity (eq. 26) */
#ifdef GRAVITY_DIFF_VELOCITY
  const float v_sig_u =
      sqrtf(v_ij[0] * v_ij[0] + v_ij[1] * v_ij[1] + v_ij[2] * v_ij[2]);
#else
  const float v_sig_u =
      sqrtf(2. * fabsf(pressurei - pressurej) / (rhoi + rhoj));
#endif

  /* Diffusion term (eq. 24) */
  pi->u_dt += -const_diffusion_alpha * mj * delta_u * v_sig_u * norm_sum_G /
              (rhoi + rhoj);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dx * hi_inv * hid_inv;

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
