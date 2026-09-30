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
 * @file MAGMA/hydro_part.h
 * @brief MAGMA-2 implementation of SPH following Rosswog+2020 (Particle
 * interactions)
 *
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
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
  pi->gradient.c_matrix_inv.yy += common_term * dx[1] * dx[1];
  pi->gradient.c_matrix_inv.xy += common_term * dx[0] * dx[1];
#endif
#if defined(HYDRO_DIMENSION_3D)
  pi->gradient.c_matrix_inv.zz += common_term * dx[2] * dx[2];
  pi->gradient.c_matrix_inv.xz += common_term * dx[0] * dx[2];
  pi->gradient.c_matrix_inv.yz += common_term * dx[1] * dx[2];
#endif

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
  const int use_base_SPH_i = pi->use_base_SPH;
  const int use_base_SPH_j = pj->use_base_SPH;

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

  /* Compute signal velocity (eq. 36) modified to add dimension on the
   * denominator. This is the magnitude of the approach speed (>= 0), zero for
   * receding particles. */
  const float mu_tilde_i =
      -fac_mu * hi * omega_ij /
      (r * r + magma_viscosity.mu_softening * hi * hi);

  /* De-dimentionalised distances (eq. 16, recall dx = xi - xj)*/
  const float eta_i[3] = {dx[0] * hi_inv, dx[1] * hi_inv, dx[2] * hi_inv};
  const float eta_j[3] = {-dx[0] * hj_inv, -dx[1] * hj_inv, -dx[2] * hj_inv};

  /* Norms of the eta vectors (eq. 16) */
  const float eta_square_i =
      eta_i[0] * eta_i[0] + eta_i[1] * eta_i[1] + eta_i[2] * eta_i[2];
  const float eta_square_j =
      eta_j[0] * eta_j[0] + eta_j[1] * eta_j[1] + eta_j[2] * eta_j[2];

  /* Reconstructed velocities at the mid-point (before reconstruction) */
  float v_rec_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
  float v_rec_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

  /* Reconstructed internal energies at the mid-point (before reconstruction) */
  float u_rec_i = pi->u;
  float u_rec_j = pj->u;

#ifdef USE_ZEROTH_ORDER_VELOCITIES
  const int force_zeroth_order = 1;
#else
  const int force_zeroth_order = 0;
#endif

  /* Slope limiter applied to the velocities (0: no reconstruction). Also
   * used for the Hubble flow, see below. */
  float Phi_vel = 0.f;

  /* Reconstruct v and u at the interface unless one of the particles is weird
   */
  if (!use_base_SPH_i && !use_base_SPH_j && !force_zeroth_order) {

    /* Vectors from the particles to the mid-point */
    const float delta_i[3] = {-0.5f * dx[0], -0.5f * dx[1], -0.5f * dx[2]};
    const float delta_j[3] = {-delta_i[0], -delta_i[1], -delta_i[2]};

    /* Terms entering the limiter (eq. 23) */
    const float eta_ij = sqrtf(fminf(eta_square_i, eta_square_j));
    const float eta_crit = magma_viscosity.eta_crit;
    const float width_inv = 1.f / magma_viscosity.limiter_width;

    /* Van Leer limiter fraction (eq. 22) */
    const float A_ij_vel_num = pi->force.gradient_vx[0] * dx[0] * dx[0] +
                               pi->force.gradient_vx[1] * dx[0] * dx[1] +
                               pi->force.gradient_vx[2] * dx[0] * dx[2] +
                               pi->force.gradient_vy[0] * dx[1] * dx[0] +
                               pi->force.gradient_vy[1] * dx[1] * dx[1] +
                               pi->force.gradient_vy[2] * dx[1] * dx[2] +
                               pi->force.gradient_vz[0] * dx[2] * dx[0] +
                               pi->force.gradient_vz[1] * dx[2] * dx[1] +
                               pi->force.gradient_vz[2] * dx[2] * dx[2] +
                               a2_Hubble * r2; /* Hubble flow: a^2 H I */

    const float A_ij_vel_den = pj->force.gradient_vx[0] * dx[0] * dx[0] +
                               pj->force.gradient_vx[1] * dx[0] * dx[1] +
                               pj->force.gradient_vx[2] * dx[0] * dx[2] +
                               pj->force.gradient_vy[0] * dx[1] * dx[0] +
                               pj->force.gradient_vy[1] * dx[1] * dx[1] +
                               pj->force.gradient_vy[2] * dx[1] * dx[2] +
                               pj->force.gradient_vz[0] * dx[2] * dx[0] +
                               pj->force.gradient_vz[1] * dx[2] * dx[1] +
                               pj->force.gradient_vz[2] * dx[2] * dx[2] +
                               a2_Hubble * r2; /* Hubble flow: a^2 H I */

    const float A_ij_vel =
        A_ij_vel_den != 0.f ? A_ij_vel_num / A_ij_vel_den : 0.f;

    /* Slope limiter exponential term (eq. 21, right term) */
    const float delta_eta = (eta_ij - eta_crit) * width_inv;
    const float exp_term =
        eta_ij < eta_crit ? expf(-delta_eta * delta_eta) : 1.f;

    /* Van Leer limiter (eq. 21).
     * Slopes of opposite signs (A <= 0) mean an extremum between the
     * particles: no reconstruction. Also avoids the division by zero at A = -1.
     */
    const float fraction_vel =
        (A_ij_vel > 0.f)
            ? 4.f * A_ij_vel / ((1.f + A_ij_vel) * (1.f + A_ij_vel))
            : 0.f;

    const float Phi_ij_vel = fminf(1.f, fraction_vel) * exp_term;
    Phi_vel = Phi_ij_vel;

    /* Mid-point reconstruction, first order (eq. 17) */
    v_rec_i[0] += Phi_ij_vel * pi->force.gradient_vx[0] * delta_i[0];
    v_rec_i[0] += Phi_ij_vel * pi->force.gradient_vx[1] * delta_i[1];
    v_rec_i[0] += Phi_ij_vel * pi->force.gradient_vx[2] * delta_i[2];
    v_rec_i[1] += Phi_ij_vel * pi->force.gradient_vy[0] * delta_i[0];
    v_rec_i[1] += Phi_ij_vel * pi->force.gradient_vy[1] * delta_i[1];
    v_rec_i[1] += Phi_ij_vel * pi->force.gradient_vy[2] * delta_i[2];
    v_rec_i[2] += Phi_ij_vel * pi->force.gradient_vz[0] * delta_i[0];
    v_rec_i[2] += Phi_ij_vel * pi->force.gradient_vz[1] * delta_i[1];
    v_rec_i[2] += Phi_ij_vel * pi->force.gradient_vz[2] * delta_i[2];

    v_rec_j[0] += Phi_ij_vel * pj->force.gradient_vx[0] * delta_j[0];
    v_rec_j[0] += Phi_ij_vel * pj->force.gradient_vx[1] * delta_j[1];
    v_rec_j[0] += Phi_ij_vel * pj->force.gradient_vx[2] * delta_j[2];
    v_rec_j[1] += Phi_ij_vel * pj->force.gradient_vy[0] * delta_j[0];
    v_rec_j[1] += Phi_ij_vel * pj->force.gradient_vy[1] * delta_j[1];
    v_rec_j[1] += Phi_ij_vel * pj->force.gradient_vy[2] * delta_j[2];
    v_rec_j[2] += Phi_ij_vel * pj->force.gradient_vz[0] * delta_j[0];
    v_rec_j[2] += Phi_ij_vel * pj->force.gradient_vz[1] * delta_j[1];
    v_rec_j[2] += Phi_ij_vel * pj->force.gradient_vz[2] * delta_j[2];

    /* Now, same for the internal energy */

    const float A_ij_u_num = pi->force.gradient_u[0] * dx[0] +
                             pi->force.gradient_u[1] * dx[1] +
                             pi->force.gradient_u[2] * dx[2];

    const float A_ij_u_den = pj->force.gradient_u[0] * dx[0] +
                             pj->force.gradient_u[1] * dx[1] +
                             pj->force.gradient_u[2] * dx[2];

    const float A_ij_u = A_ij_u_den != 0.f ? A_ij_u_num / A_ij_u_den : 0.f;

    /* Van Leer limiter (eq. 21), as above */
    const float fraction_u =
        (A_ij_u > 0.f) ? 4.f * A_ij_u / ((1.f + A_ij_u) * (1.f + A_ij_u)) : 0.f;

    const float Phi_ij_u = fminf(1.f, fraction_u) * exp_term;

    /* Mid-point reconstruction, first order (eq. 17) */
    u_rec_i += Phi_ij_u * pi->force.gradient_u[0] * delta_i[0];
    u_rec_i += Phi_ij_u * pi->force.gradient_u[1] * delta_i[1];
    u_rec_i += Phi_ij_u * pi->force.gradient_u[2] * delta_i[2];

    u_rec_j += Phi_ij_u * pj->force.gradient_u[0] * delta_j[0];
    u_rec_j += Phi_ij_u * pj->force.gradient_u[1] * delta_j[1];
    u_rec_j += Phi_ij_u * pj->force.gradient_u[2] * delta_j[2];

    /* Limiter preventing inversion: the slope limiter above only compares
     * the two slopes, not the actual difference. For steep profiles, the
     * reconstructed difference can change sign (or appear between equal
     * values), which would make the conduction transport energy from the
     * colder to the hotter particle. Use no difference in that case. */
    if ((pi->u - pj->u) * (u_rec_i - u_rec_j) <= 0.f) {
      u_rec_i = u_rec_j = 0.5f * (pi->u + pj->u);
    }
  }

  /* Difference in velocity at the mid-point */
  const float v_rec_ij[3] = {v_rec_i[0] - v_rec_j[0], v_rec_i[1] - v_rec_j[1],
                             v_rec_i[2] - v_rec_j[2]};

  /* Normalised relative velocity (eq. 15) */
  const float vel_rel_i =
      eta_i[0] * v_rec_ij[0] + eta_i[1] * v_rec_ij[1] + eta_i[2] * v_rec_ij[2];
  const float vel_rel_j = eta_j[0] * -v_rec_ij[0] + eta_j[1] * -v_rec_ij[1] +
                          eta_j[2] * -v_rec_ij[2];

  /* Includes the hubble flow term (a^2 H dx, projected on eta = dx / h).
   * The Hubble flow is a linear field: reconstructed to the mid-point like the
   * peculiar velocities (with the same limiter), its difference reduces to
   * (1 - Phi) a^2 H dx. fac_mu converts the internal velocities to the units
   * of the (comoving) sound speed they are combined with in Q. */
  const float Hubble_rec = (1.f - Phi_vel) * a2_Hubble * r2;
  const float vel_rel_Hubble_i = fac_mu * (vel_rel_i + Hubble_rec * hi_inv);
  const float vel_rel_Hubble_j = fac_mu * (vel_rel_j + Hubble_rec * hj_inv);

  /* Construct the gradient functions (eq. 4 and 5) */
  float G_i[3] = {0.f}, G_j[3] = {0.f};
  sym_matrix_multiply_by_vector(G_i, &pi->force.c_matrix, dx);
  sym_matrix_multiply_by_vector(G_j, &pj->force.c_matrix, dx);

  /* Note we multiply by -1 as dx is (pi - pj) and not (pj - pi) */
  G_i[0] *= -wi * hid_inv;
  G_i[1] *= -wi * hid_inv;
  G_i[2] *= -wi * hid_inv;
  G_j[0] *= -wj * hjd_inv;
  G_j[1] *= -wj * hjd_inv;
  G_j[2] *= -wj * hjd_inv;

  const float cos_limit = magma_viscosity.cos_angle_limit;

#ifdef TRADITIONAL_SPH_ACCELERATION_TERM
  /* MI1 uses G_i and G_j individually: test each of them. The averaged G can
   * look fine while one of them is tilted too much or points the wrong way.
   * |G.dx| < cos(limit) |G| r  <=>  angle > limit. */
  const float G_i_norm =
      sqrtf(G_i[0] * G_i[0] + G_i[1] * G_i[1] + G_i[2] * G_i[2]);
  const float G_j_norm =
      sqrtf(G_j[0] * G_j[0] + G_j[1] * G_j[1] + G_j[2] * G_j[2]);
  const float G_i_dot_dx = G_i[0] * dx[0] + G_i[1] * dx[1] + G_i[2] * dx[2];
  const float G_j_dot_dx = G_j[0] * dx[0] + G_j[1] * dx[1] + G_j[2] * dx[2];
  const int G_ij_misaligned = (fabsf(G_i_dot_dx) < cos_limit * G_i_norm * r) ||
                              (fabsf(G_j_dot_dx) < cos_limit * G_j_norm * r);
  const int G_ij_wrong_sign = (G_i_dot_dx > 0.f) || (G_j_dot_dx > 0.f);
#else
  /* Verify that the G vector has the right direction */
  const float G_ij[3] = {0.5f * (G_i[0] + G_j[0]),  /* x */
                         0.5f * (G_i[1] + G_j[1]),  /* y */
                         0.5f * (G_i[2] + G_j[2])}; /* z */

  /* Angle between G and the axis linking the particles */
  const float G_ij_norm =
      sqrtf(G_ij[0] * G_ij[0] + G_ij[1] * G_ij[1] + G_ij[2] * G_ij[2]);
  const float G_ij_dot_dx = G_ij[0] * dx[0] + G_ij[1] * dx[1] + G_ij[2] * dx[2];

  /* MI2 only uses the average of G_i and G_j.
   * |G.dx| < cos(limit) |G| r  <=>  angle > limit. */
  const int G_ij_misaligned = fabsf(G_ij_dot_dx) < cos_limit * G_ij_norm * r;

  /* Check whether the sign of the reconstructed interface normals is wrong */
  const int G_ij_wrong_sign = (G_ij_dot_dx > 0.f);
#endif

  /* if (G_ij_misaligned) */
  /*   warning( */
  /*       "Misaligned dx=[%e %e %e] G_ij=[%e %e %e] G_ij.dx=%e use_SPH_i=%d
   * " */
  /*       "use_SPH_j=%d", */
  /*       dx[0], dx[1], dx[2], G_ij[0], G_ij[1], G_ij[2], G_ij_dot_dx, */
  /*       use_base_SPH_i, use_base_SPH_j); */

  /* if (G_ij_wrong_sign) */
  /*   warning( */
  /* 	    "Wrong sign! dx=[%e %e %e] G_ij=[%e %e %e]  use_SPH_i=%d " */
  /* 	    "use_SPH_j=%d", */
  /* 	    dx[0], dx[1], dx[2], G_ij[0], G_ij[1], G_ij[2], */
  /* 	    use_base_SPH_i, use_base_SPH_j); */

#ifdef USE_STANDARD_KERNEL_GRADIENTS
  const int force_standard_kernel = 1;
#else
  const int force_standard_kernel = 0;
#endif

  /* Default to the traditional SPH gradW term if one of the particles is weird
   */
  if (use_base_SPH_i || use_base_SPH_j || force_standard_kernel ||
      G_ij_misaligned || G_ij_wrong_sign) {

    const float wi_dr = hid_inv * hi_inv * wi_dx;
    const float wj_dr = hjd_inv * hj_inv * wj_dx;
    G_i[0] = wi_dr * r_inv * dx[0];
    G_i[1] = wi_dr * r_inv * dx[1];
    G_i[2] = wi_dr * r_inv * dx[2];
    G_j[0] = wj_dr * r_inv * dx[0];
    G_j[1] = wj_dr * r_inv * dx[1];
    G_j[2] = wj_dr * r_inv * dx[2];
  }

  /* Terms entering the viscosity (eq. 15).
   * Only for pairs that are actually approaching (raw velocities including the
   * Hubble flow, as in the other SPH schemes): the reconstructed velocities
   * can indicate compression while the particles recede, in which case Q
   * would do negative work (Q v_ij . G < 0) and cool the gas. */
  const int pair_approaching = (dvdr_Hubble < 0.f);
  const float eps_squared = magma_viscosity.epsilon * magma_viscosity.epsilon;
  const float mu_i =
      pair_approaching
          ? fminf(0.f, vel_rel_Hubble_i / (eta_square_i + eps_squared))
          : 0.f;
  const float mu_j =
      pair_approaching
          ? fminf(0.f, vel_rel_Hubble_j / (eta_square_j + eps_squared))
          : 0.f;

  /* Relative velocity including the Hubble flow (a^2 H dx) */
  const float v_ij_Hubble[3] = {v_ij[0] + a2_Hubble * dx[0],
                                v_ij[1] + a2_Hubble * dx[1],
                                v_ij[2] + a2_Hubble * dx[2]};

  /* The viscous pressure heats at a rate Q v_ij . G (see the energy
   * equation below). G can be tilted from dx (up to the angle limit parameter),
   * so for shear-dominated pairs v_ij . G can be negative although the pair
   * approaches along dx: Q would then turn internal energy into kinetic
   * energy. Only use Q where its heating term is positive. */
#ifdef TRADITIONAL_SPH_ACCELERATION_TERM
  /* Each particle's Q works through its own gradient function */
  const int Q_dissipative_i =
      (v_ij_Hubble[0] * G_i[0] + v_ij_Hubble[1] * G_i[1] +
       v_ij_Hubble[2] * G_i[2]) > 0.f;
  const int Q_dissipative_j =
      (v_ij_Hubble[0] * G_j[0] + v_ij_Hubble[1] * G_j[1] +
       v_ij_Hubble[2] * G_j[2]) > 0.f;
#else
  /* Both Q work through the averaged gradient function */
  const int Q_dissipative_i =
      (v_ij_Hubble[0] * (G_i[0] + G_j[0]) + v_ij_Hubble[1] * (G_i[1] + G_j[1]) +
       v_ij_Hubble[2] * (G_i[2] + G_j[2])) > 0.f;
  const int Q_dissipative_j = Q_dissipative_i;
#endif

  /* Eq. 14 */
  const float visc_alpha = magma_viscosity.alpha;
  const float visc_beta = magma_viscosity.beta;
  const float Qi =
      Q_dissipative_i
          ? rhoi * (-visc_alpha * ci * mu_i + visc_beta * mu_i * mu_i)
          : 0.f;
  const float Qj =
      Q_dissipative_j
          ? rhoj * (-visc_alpha * cj * mu_j + visc_beta * mu_j * mu_j)
          : 0.f;

#ifdef TRADITIONAL_SPH_ACCELERATION_TERM

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

  /* Same, including the Hubble flow (a^2 H dx) */
  const float v_ij_Hubble_dot_G_i =
      v_ij_dot_G_i +
      a2_Hubble * (dx[0] * G_i[0] + dx[1] * G_i[1] + dx[2] * G_i[2]);

  /* Raw change in internal energy (eq. 3). The comoving internal energy
   * absorbs the P dV work of the Hubble expansion, hence the pressure term
   * uses the peculiar velocity only. The viscous pressure Q is not part of
   * that change of variables: its heating uses the full relative velocity */
  pi->u_dt += mj * (pressurei * v_ij_dot_G_i + Qi * v_ij_Hubble_dot_G_i) /
              (rhoi * rhoi);

#else /* Gasoline-like mixing */

  /* Average of the (possibly replaced) gradient functions */
  const float G_avg[3] = {0.5f * (G_i[0] + G_j[0]), 0.5f * (G_i[1] + G_j[1]),
                          0.5f * (G_i[2] + G_j[2])};

  const float acc_term = (pressurei + pressurej + Qi + Qj) / (rhoi * rhoj);

  /* Raw fluid acceleration (eq. 10) */
  pi->a_hydro[0] -= mj * acc_term * G_avg[0];
  pi->a_hydro[1] -= mj * acc_term * G_avg[1];
  pi->a_hydro[2] -= mj * acc_term * G_avg[2];

  /* Equivalent of div v */
  const float v_ij_dot_G_avg =
      v_ij[0] * G_avg[0] + v_ij[1] * G_avg[1] + v_ij[2] * G_avg[2];

  /* Same, including the Hubble flow (a^2 H dx) */
  const float v_ij_Hubble_dot_G_avg =
      v_ij_dot_G_avg +
      a2_Hubble * (dx[0] * G_avg[0] + dx[1] * G_avg[1] + dx[2] * G_avg[2]);

  /* Raw change in internal energy (eq. 11). The comoving internal energy
   * absorbs the P dV work of the Hubble expansion, hence the pressure term
   * uses the peculiar velocity only. The viscous pressure Q is not part of
   * that change of variables: its heating uses the full relative velocity */
  pi->u_dt += mj * (pressurei * v_ij_dot_G_avg + Qi * v_ij_Hubble_dot_G_avg) /
              (rhoi * rhoj);

#endif

  /* Difference in internal energy */
  const float delta_u = u_rec_i - u_rec_j;

  /* Norm of the G vectors */
  const float sum_G[3] = {G_i[0] + G_j[0], G_i[1] + G_j[1], G_i[2] + G_j[2]};
  const float norm_sum_G =
      sqrtf(sum_G[0] * sum_G[0] + sum_G[1] * sum_G[1] + sum_G[2] * sum_G[2]);

  /* Diffusion signal velocity (eq. 26) */
#ifdef GRAVITY_DIFF_VELOCITY
  /* Norm of the full relative velocity vector (eq. 26, |v_a - v_b|, here
   * without reconstruction), including the Hubble flow, converted to
   * sound-speed units by fac_mu.
   * Note: unlike SPHENIX, which uses only the component along dx,
   * the shear components contribute too. */
  const float v_sig_u = fac_mu * sqrtf(v_ij_Hubble[0] * v_ij_Hubble[0] +
                                       v_ij_Hubble[1] * v_ij_Hubble[1] +
                                       v_ij_Hubble[2] * v_ij_Hubble[2]);
#else
  const float v_sig_u =
      sqrtf(2.f * fabsf(pressurei - pressurej) / (rhoi + rhoj));
#endif

  /* Diffusion term (eq. 24) */
  pi->u_dt += -magma_diffusion.alpha * mj * delta_u * v_sig_u * norm_sum_G /
              (rhoi + rhoj);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dx * hi_inv * hid_inv;

  /* Update the signal velocity. */
  pi->force.mu_tilde = max(pi->force.mu_tilde, mu_tilde_i);

  /* Update the signal speed of the time-step with the neighbour's sound
   * speed: as conservative as the global time-step of the paper, which
   * limits every particle by its hottest neighbour's Courant condition. */
  pi->force.c_sig = max(pi->force.c_sig, cj);
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
