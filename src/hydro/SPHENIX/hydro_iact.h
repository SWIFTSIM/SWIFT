/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SPHENIX_HYDRO_IACT_H
#define SWIFT_SPHENIX_HYDRO_IACT_H

/**
 * @file SPHENIX/hydro_iact.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added SPHENIX physics (Borrow 2020) (interaction routines)
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "fvpm_geometry.h"
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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  float wi, wj, wi_dx, wj_dx;

  const double r = sqrt(r2);

  /* Get the masses. */
  const double mi = pi->mass;
  const double mj = pj->mass;

  /* Compute density of pi. */
  const double hi_inv = 1.0 / hi;
  const double ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
  adaptive_softening_add_correction_term(pi, ui, hi_inv, mj);

  /* Collect data for FVPM matrix construction */
  fvpm_accumulate_geometry_and_matrix(pi, wi, dx);
  fvpm_update_centroid_left(pi, dx, wi);

  /* Compute density of pj. */
  const double hj_inv = 1.0 / hj;
  const double uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);
  adaptive_softening_add_correction_term(pj, uj, hj_inv, mi);

  /* Collect data for FVPM matrix construction */
  fvpm_accumulate_geometry_and_matrix(pj, wj, dx);
  fvpm_update_centroid_right(pj, dx, wj);

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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, const struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  float wi, wi_dx;

  /* Get the masses. */
  const double mj = pj->mass;

  /* Get r and r inverse. */
  const double r = sqrt(r2);

  const double h_inv = 1.0 / hi;
  const double ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  /* Collect data for FVPM matrix construction */
  fvpm_accumulate_geometry_and_matrix(pi, wi, dx);
  fvpm_update_centroid_left(pi, dx, wi);

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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  double dv[3], curlvr[3];

  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  const double mi = pi->mass;
  const double mj = pj->mass;

  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx, wj, wj_dx;

  const double ui = r / hi;
  const double uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* Cosmology terms for the signal velocity */
  const double fac_mu = pow_three_gamma_minus_five_over_two(a);
  const double a2_Hubble = a * a * H;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const double dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  /* Add Hubble flow */

  const double dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const double omega_ij = fmin(dvdr_Hubble, 0.0);
  const double mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const double new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Update if we need to */
  pi->viscosity.v_sig = fmax(pi->viscosity.v_sig, new_v_sig);
  pj->viscosity.v_sig = fmax(pj->viscosity.v_sig, new_v_sig);

  /* Now we need to compute the div terms */
  const double faci = mj * f_ij * wi_dx * r_inv;
  const double facj = mi * f_ji * wj_dx * r_inv;

  pi->viscosity.div_v -= faci * dvdr;
  pj->viscosity.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->viscosity.rot_v[0] += faci * curlvr[0];
  pi->viscosity.rot_v[1] += faci * curlvr[1];
  pi->viscosity.rot_v[2] += faci * curlvr[2];

  /* Negative because of the change in sign of dx & dv. */
  pj->viscosity.rot_v[0] += facj * curlvr[0];
  pj->viscosity.rot_v[1] += facj * curlvr[1];
  pj->viscosity.rot_v[2] += facj * curlvr[2];

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  const double delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;
  pj->diffusion.laplace_u -= pi->mass * delta_u_factor * wj_dx / pi->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const double alpha_i = pi->viscosity.alpha;
  const double alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = fmax(pi->force.alpha_visc_max_ngb, alpha_j);
  pj->force.alpha_visc_max_ngb = fmax(pj->force.alpha_visc_max_ngb, alpha_i);

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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  double dv[3], curlvr[3];

  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  const double mj = pj->mass;

  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const double ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;

  /* Cosmology terms for the signal velocity */
  const double fac_mu = pow_three_gamma_minus_five_over_two(a);
  const double a2_Hubble = a * a * H;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const double dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  /* Add Hubble flow */

  const double dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const double omega_ij = fmin(dvdr_Hubble, 0.0);
  const double mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const double new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Update if we need to */
  pi->viscosity.v_sig = fmax(pi->viscosity.v_sig, new_v_sig);

  /* Now we need to compute the div terms */
  const double faci = mj * f_ij * wi_dx * r_inv;

  pi->viscosity.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->viscosity.rot_v[0] += faci * curlvr[0];
  pi->viscosity.rot_v[1] += faci * curlvr[1];
  pi->viscosity.rot_v[2] += faci * curlvr[2];

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  const double delta_u_factor = (pi->u - pj->u) * r_inv;
  pi->diffusion.laplace_u += pj->mass * delta_u_factor * wi_dx / pj->rho;

  /* Set the maximal alpha from the previous step over the neighbours
   * (this is used to limit the diffusion in hydro_prepare_force) */
  const double alpha_j = pj->viscosity.alpha;
  pi->force.alpha_visc_max_ngb = fmax(pi->force.alpha_visc_max_ngb, alpha_j);

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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  /* Cosmological factors entering the EoMs */
  const double fac_mu = pow_three_gamma_minus_five_over_two(a);
  const double a2_Hubble = a * a * H;

  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mj = pj->mass;
  const double mi = pi->mass;

  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  const double pressurei = pi->force.pressure;
  const double pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const double hj_inv = 1.0 / hj;
  const double hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const double xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const double wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const double dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const double dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const double omega_ij = fmin(dvdr_Hubble, 0.0);
  const double mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const double v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* Construct the full viscosity term */
  const double rho_ij = rhoi + rhoj;
  const double alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const double visc = -0.5 * alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const double visc_acc_term =
      0.5 * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* Compute gradient terms */
  const double P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const double P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term */
  const double sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Adaptive softening acceleration term */
  const double adapt_soft_acc_term =
      adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv);

  /* Assemble the acceleration */
  const double acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const double sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const double sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const double visc_du_term = 0.5 * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */
  const double alpha_diff =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const double v_diff = alpha_diff * 0.5 *
                       (sqrt(2.0 * fabs(pressurei - pressurej) / rho_ij) +
                        fabs(fac_mu * r_inv * dvdr_Hubble));
  /* wi_dx + wj_dx / 2 is F_ij */
  const double diff_du_term =
      v_diff * (pi->u - pj->u) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

  /* Assemble the energy equation term */
  const double du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;
  const double du_dt_j = sph_du_term_j + visc_du_term - diff_du_term;

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
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const double r2, const float dx[3], const double hi, const double hj,
    struct part* restrict pi, const struct part* restrict pj, const double mu_0,
    const double a, const double H) {

  /* Cosmological factors entering the EoMs */
  const double fac_mu = pow_three_gamma_minus_five_over_two(a);
  const double a2_Hubble = a * a * H;

  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mi = pi->mass;
  const double mj = pj->mass;

  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  const double pressurei = pi->force.pressure;
  const double pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const double hj_inv = 1.0 / hj;
  const double hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const double xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const double wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  const double dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Includes the hubble flow term; not used for du/dt */
  const double dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const double omega_ij = fmin(dvdr_Hubble, 0.0);
  const double mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const double v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* Construct the full viscosity term */
  const double rho_ij = rhoi + rhoj;
  const double alpha = pi->viscosity.alpha + pj->viscosity.alpha;
  const double visc = -0.5 * alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const double visc_acc_term =
      0.5 * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* Compute gradient terms */
  const double P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const double P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term */
  const double sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Adaptive softening acceleration term */
  const double adapt_soft_acc_term =
      adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv);

  /* Assemble the acceleration */
  const double acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const double sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const double visc_du_term = 0.5 * visc_acc_term * dvdr_Hubble;

  /* Diffusion term */
  /* Combine the alpha_diff into a pressure-based switch -- this allows the
   * alpha from the highest pressure particle to dominate, so that the
   * diffusion limited particles always take precedence - another trick to
   * allow the scheme to work with thermal feedback. */
  const double alpha_diff =
      (pressurei * pi->diffusion.alpha + pressurej * pj->diffusion.alpha) /
      (pressurei + pressurej);
  const double v_diff = alpha_diff * 0.5 *
                       (sqrt(2.0 * fabs(pressurei - pressurej) / rho_ij) +
                        fabs(fac_mu * r_inv * dvdr_Hubble));
  /* wi_dx + wj_dx / 2 is F_ij */
  const double diff_du_term =
      v_diff * (pi->u - pj->u) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

  /* Assemble the energy equation term */
  const double du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
}

#endif /* SWIFT_SPHENIX_HYDRO_IACT_H */
