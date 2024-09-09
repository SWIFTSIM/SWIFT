/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch (optional).
 * No thermal conduction term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "const.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

//### temp
#include "../REMIX/hydro_strength.h"

/**
 * @brief Calculates the stress tensor for force interaction. No strength if
 * either particle is fluid (or not using strength at all)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_set_pairwise_stress_tensors(float pairwise_stress_tensor_i[3][3],
                                  float pairwise_stress_tensor_j[3][3],
                                  const struct part *restrict pi,
                                  const struct part *restrict pj, const float r,
                                  const float pressurei,
                                  const float pressurej) {

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
  hydro_set_pairwise_stress_tensors_strength(
      pairwise_stress_tensor_i, pairwise_stress_tensor_j, pi, pj, r);
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

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

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;
  pj->density.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];

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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

#ifdef MATERIAL_STRENGTH
  hydro_runner_iact_nonsym_density_extra_strength(pi, pj, dx, wi, wi_dx);
#endif /* MATERIAL_STRENGTH */
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
    const float H) {}

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
    const float H) {}

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

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each other? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float v_sig = signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

  float Gi[3], Gj[3], G_mean[3];
  for (int i = 0; i < 3; i++) {
    Gi[i] = wi_dr * f_ij * dx[i] * r_inv;
    Gj[i] = -wj_dr * f_ji * dx[i] * r_inv;
    G_mean[i] = 0.5f * (Gi[i] - Gj[i]);
  }

  // Convert the pressures into "stress" tensors for fluids, S = -P * I(3),
  // or set the stress for solid particle pairs with strength
  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  hydro_set_pairwise_stress_tensors(pairwise_stress_tensor_i,
                                    pairwise_stress_tensor_j, pi, pj, r,
                                    pressurei, pressurej);
    
  float stress_tensor_term_i[3], stress_tensor_term_j[3];
  for (int i = 0; i < 3; i++) {
    stress_tensor_term_i[i] = 0.f;
    stress_tensor_term_j[i] = 0.f;
    for (int j = 0; j < 3; j++) {
      stress_tensor_term_i[i] += pairwise_stress_tensor_i[i][j] * Gi[j];
      stress_tensor_term_j[i] += pairwise_stress_tensor_j[i][j] * Gj[j];
    }
  }

  /* Construct acceleration terms */
  float visc_acc_term[3], sph_acc_term[3], adapt_soft_acc_term[3];
  for (int i = 0; i < 3; i++) {  
     sph_acc_term[i] = stress_tensor_term_i[i] / (rhoi * rhoi) -  stress_tensor_term_j[i] / (rhoj * rhoj);
     visc_acc_term[i] = -visc * G_mean[i];
     adapt_soft_acc_term[i] =
         -adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv) * dx[i];
  }

  /* Use the force Luke ! */
  for (int i = 0; i < 3; i++) {  
    pi->a_hydro[i] += mj * (sph_acc_term[i] + visc_acc_term[i] + adapt_soft_acc_term[i]);
    pj->a_hydro[i] -= mi * (sph_acc_term[i] + visc_acc_term[i] + adapt_soft_acc_term[i]);
  }

  float dv_dot_stress_term_i = 0.f, dv_dot_stress_term_j = 0.f;
  float dv_dot_visc = 0.f;
  float dv_dot_G_i = 0.f;
  float dv_dot_G_j = 0.f;  
  for (int i = 0; i < 3; i++) {
    dv_dot_stress_term_i += (pi->v[i] - pj->v[i]) * stress_tensor_term_i[i];
    dv_dot_stress_term_j += -(pi->v[i] - pj->v[i]) * stress_tensor_term_j[i];
    dv_dot_visc += (pi->v[i] - pj->v[i]) * visc_acc_term[i];
    dv_dot_G_i += (pi->v[i] - pj->v[i]) * Gi[i];
    dv_dot_G_j += (pj->v[i] - pi->v[i]) * Gj[i];
  }

  /* Get the time derivative for u. */
  const float sph_du_term_i = -(dv_dot_stress_term_i) / (rhoi * rhoi);
  const float sph_du_term_j = -(dv_dot_stress_term_j) / (rhoj * rhoj);

  /* Viscosity term */
  const float visc_du_term = -0.5f * dv_dot_visc;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

#if defined(MATERIAL_STRENGTH)    
  pi->drho_dt += mj * dv_dot_G_i;
  pj->drho_dt += mi * dv_dot_G_j;
    
  hydro_runner_iact_force_extra_strength(pi, pj, dx, Gi, Gj);
#endif /* MATERIAL_STRENGTH */    
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
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each other? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

  float Gi[3], Gj[3], G_mean[3];
  for (int i = 0; i < 3; i++) {
    Gi[i] = wi_dr * f_ij * dx[i] * r_inv;
    Gj[i] = -wj_dr * f_ji * dx[i] * r_inv;
    G_mean[i] = 0.5f * (Gi[i] - Gj[i]);
  }

  // Convert the pressures into "stress" tensors for fluids, S = -P * I(3),
  // or set the stress for solid particle pairs with strength
  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  hydro_set_pairwise_stress_tensors(pairwise_stress_tensor_i,
                                    pairwise_stress_tensor_j, pi, pj, r,
                                    pressurei, pressurej);
    
  float stress_tensor_term_i[3], stress_tensor_term_j[3];
  for (int i = 0; i < 3; i++) {
    stress_tensor_term_i[i] = 0.f;
    stress_tensor_term_j[i] = 0.f;
    for (int j = 0; j < 3; j++) {
      stress_tensor_term_i[i] += pairwise_stress_tensor_i[i][j] * Gi[j];
      stress_tensor_term_j[i] += pairwise_stress_tensor_j[i][j] * Gj[j];
    }
  }

  /* Construct acceleration terms */
  float visc_acc_term[3], sph_acc_term[3], adapt_soft_acc_term[3];
  for (int i = 0; i < 3; i++) {  
     sph_acc_term[i] = stress_tensor_term_i[i] / (rhoi * rhoi) -  stress_tensor_term_j[i] / (rhoj * rhoj);
     visc_acc_term[i] = -visc * G_mean[i];
     adapt_soft_acc_term[i] =
         -adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv) * dx[i];
  }

  /* Use the force Luke ! */
  for (int i = 0; i < 3; i++) {  
    pi->a_hydro[i] += mj * (sph_acc_term[i] + visc_acc_term[i] + adapt_soft_acc_term[i]);
  }

  float dv_dot_stress_term_i = 0.f;
  float dv_dot_visc = 0.f;
  float dv_dot_G_i = 0.f;
  for (int i = 0; i < 3; i++) {
    dv_dot_stress_term_i += (pi->v[i] - pj->v[i]) * stress_tensor_term_i[i];
    dv_dot_visc += (pi->v[i] - pj->v[i]) * visc_acc_term[i];
    dv_dot_G_i += (pi->v[i] - pj->v[i]) * Gi[i];
  }

  /* Get the time derivative for u. */
  const float sph_du_term_i = -(dv_dot_stress_term_i) / (rhoi * rhoi);

  /* Viscosity term */
  const float visc_du_term = -0.5f * dv_dot_visc;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

#if defined(MATERIAL_STRENGTH)    
  pi->drho_dt += mj * dv_dot_G_i;

  hydro_runner_iact_nonsym_force_extra_strength(pi, pj, dx, Gi);
#endif /* MATERIAL_STRENGTH */       
}

#endif /* SWIFT_PLANETARY_HYDRO_IACT_H */
