/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_VECTOR_POTENTIAL_MHD_H
#define SWIFT_VECTOR_POTENTIAL_MHD_H

#include "hydro.h"
#include "mhd_parameters.h"
#include "space.h"

#include <float.h>

/**
 * @brief Returns the magnetic energy contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_energy(
    const struct part *p, const struct xpart *xp, const float mu_0,
    const float a) {

  const float afact = pow(a, 2.f * mhd_comoving_factor);
  // convert to physical
  const float B2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  return 0.5f * afact * B2 / mu_0 * p->mass / p->rho;
}
/**
 * @brief Returns the magnetic field squared contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static float mhd_get_Bms(
    const struct part *p, const struct xpart *xp) {

  const float B2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  return B2;
}
/**
 * @brief Returns the magnetic field divergence of a particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */

__attribute__((always_inline)) INLINE static float mhd_get_magnetic_divergence(
    const struct part *p, const struct xpart *xp) {

  return p->mhd_data.divB;
}

/**
 * @brief Returns the magnetic helicity contained in the particle.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_magnetic_helicity(
    const struct part *p, const struct xpart *xp) {

  return p->mhd_data.APred[0] * p->mhd_data.BPred[0] +
         p->mhd_data.APred[1] * p->mhd_data.BPred[1] +
         p->mhd_data.APred[2] * p->mhd_data.BPred[2];
}

__attribute__((always_inline)) INLINE static float mhd_get_cross_helicity(
    const struct part *p, const struct xpart *xp) {

  return (p->v[0] * p->mhd_data.BPred[0] + p->v[1] * p->mhd_data.BPred[1] +
          p->v[2] * p->mhd_data.BPred[2]);
}

/**
 * @brief Returns the magnetic field divergence error of the particle.
 *
 * This is (div B) / (B / h) and is hence dimensionless.
 *
 * @param p the #part.
 * @param xp the #xpart.
 */
__attribute__((always_inline)) INLINE static float mhd_get_divB_error(
    const struct part *p, const struct xpart *xp) {

  const float B2 = p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                   p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                   p->mhd_data.BPred[2] * p->mhd_data.BPred[2];
  const float error =
      B2 != 0.0f ? fabs(p->mhd_data.divB) * p->h / sqrtf(B2) : 0.0f;
  return error;
}

/**
 * @brief Computes the MHD time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float mhd_compute_timestep(
    const struct part *p, const struct xpart *xp,
    const struct hydro_props *hydro_properties, const struct cosmology *cosmo,
    const float mu_0) {

  const float afac_divB = cosmo->a * cosmo->a;
  // pow(cosmo->a, -mhd_comoving_factor - 0.5f);
  const float afac_resistive = cosmo->a * cosmo->a;

  float dt_divB =
      p->mhd_data.divB != 0.0f
          ? afac_divB * hydro_properties->CFL_condition *
                sqrtf(p->rho / (p->mhd_data.divB * p->mhd_data.divB) * mu_0)
          : FLT_MAX;
  const float resistive_eta = p->mhd_data.resistive_eta;
  const float dt_eta = resistive_eta != 0.0f
                           ? afac_resistive * hydro_properties->CFL_condition *
                                 p->h * p->h / resistive_eta
                           : FLT_MAX;

  return fminf(dt_eta, dt_divB);
}

/**
 * @brief Compute Alfven speed
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_Alfven_speed(const struct part *restrict p, const float mu_0) {

  /* Recover some data */
  const float rho = p->rho;

  /* B squared */
  const float B2 = (p->mhd_data.BPred[0] * p->mhd_data.BPred[0] +
                    p->mhd_data.BPred[1] * p->mhd_data.BPred[1] +
                    p->mhd_data.BPred[2] * p->mhd_data.BPred[2]);

  /* Square of Alfven speed */
  const float vA2 = B2 / (mu_0 * rho);

  return sqrtf(vA2);
}

/**
 * @brief Compute magnetosonic speed
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_magnetosonic_speed(const struct part *restrict p) {

  /* Compute square of fast magnetosonic speed */
  const float cs = hydro_get_comoving_soundspeed(p);
  const float cs2 = cs * cs;

  const float vA = p->mhd_data.Alfven_speed;
  const float vA2 = vA * vA;

  const float cms2 = cs2 + vA2;

  return sqrtf(cms2);
}

/**
 * @brief Compute fast magnetosonic wave phase veolcity
 */
__attribute__((always_inline)) INLINE static float
mhd_get_comoving_fast_magnetosonic_wave_phase_velocity(
    const float dx[3], const struct part *restrict p, const float a,
    const float mu_0) {

  /* Get r and 1/r. */
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float rho = p->rho;
  float B[3];
  B[0] = p->mhd_data.BPred[0];
  B[1] = p->mhd_data.BPred[1];
  B[2] = p->mhd_data.BPred[2];

  /* B dot r. */
  const float Br = B[0] * dx[0] + B[1] * dx[1] + B[2] * dx[2];
  const float permeability_inv = 1.0f / mu_0;

  /* Compute effective sound speeds */
  const float cs = p->force.soundspeed;
  const float cs2 = cs * cs;
  const float c_ms = mhd_get_comoving_magnetosonic_speed(p);
  const float c_ms2 = c_ms * c_ms;
  const float projection_correction = c_ms2 * c_ms2 - 4.0f * permeability_inv *
                                                          cs2 * Br * r_inv *
                                                          Br * r_inv / rho;

  const float v_fmsw2 = 0.5f * (c_ms2 + sqrtf(projection_correction));

  return sqrtf(v_fmsw2);
}

/**
 * @brief Compute the MHD signal velocity between two gas particles
 *
 * This is eq. (131) of Price D., JCoPh, 2012, Vol. 231, Issue 3
 *
 * Warning ONLY to be called just after preparation of the force loop.
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float mhd_signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta,
    const float a, const float mu_0) {

  const float v_sigi = mhd_get_comoving_magnetosonic_speed(pi);
  const float v_sigj = mhd_get_comoving_magnetosonic_speed(pj);

  const float v_sig = v_sigi + v_sigj - beta * mu_ij;

  return v_sig;
}

/**
 * @brief Adapts signal velocity to change in drifted physical internal energy of a particle at feedback events
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void
mhd_set_drifted_physical_internal_energy(struct part *p) {

  /* Re-set MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p);
  hydro_set_signal_velocity(p, fmaxf(hydro_get_signal_velocity(p), 2.0f * cms));
}

/**
 * @brief Correct the signal velocity of the particle partaking in
 * supernova (kinetic) feedback based on the velocity kick the particle receives
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param dv_phys The velocity kick received by the particle expressed in
 * physical units (note that dv_phys must be positive or equal to zero)
 */
__attribute__((always_inline)) INLINE static void
mhd_set_v_sig_based_on_velocity_kick(struct part *p,
                                       const struct cosmology *cosmo,
                                       const float dv_phys) {

  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Fast magnetosonic speed */
  const float cms = mhd_get_comoving_magnetosonic_speed(p);
  
  /* Update the signal velocity */
  hydro_set_signal_velocity(p, fmaxf(2.f * cms, hydro_get_signal_velocity(p) +
                                                    const_viscosity_beta * dv));
}

/**
 * @brief Returns the Gauge Scalar Phi evolution
 * time the particle. Gauge all variables in full step
 *
 * @param p The particle of interest
 * @param Gauge Gauge
 */
__attribute__((always_inline)) INLINE static float mhd_get_dGau_dt(
    const struct part *restrict p, const struct cosmology *c, const float mu0) {

  const float Gauge = p->mhd_data.Gau;
  const float v_sig = 0.5f * hydro_get_signal_velocity(p);
  const float afac1 = pow(c->a, 2.f * c->a_factor_sound_speed);
  const float afac2 = pow(c->a, (c->a_factor_sound_speed + 1.f));

  /* Hyperbolic term */
  const float Source_Term = 2.0f * afac1 * p->mhd_data.divA * (v_sig * v_sig);
  /* Parabolic evolution term */
  const float Damping_Term = 1.0f * afac2 * v_sig * Gauge / p->h;
  /* Density change term */
  const float DivV_Term = 0.0 * hydro_get_div_v(p) * Gauge;
  /* Cosmological term */
  const float Hubble_Term = (2.f + mhd_comoving_factor) * c->H * Gauge;

  return (-Source_Term - Damping_Term - DivV_Term - Hubble_Term) * 0.f * c->a *
         c->a;
}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density loop over neighbours. Typically, all fields of the
 * density sub-structure of a particle get zeroed in here.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_init_part(
    struct part *p) {
  
  zero_sym_matrix(&p->mhd_data.dens.d_mat_inv);
  for (int i = 0; i < 3; ++i) 
    for (int j = 0; j < 3; ++j){ 
      p->mhd_data.dens.Mat_b[i][j] = 0.f;
      p->mhd_data.grad.Mat_da[i][j] = 0.f;
      }
}

/**
 * @brief Finishes the density calculation.
 *
 * Multiplies the density and number of neighbours by the appropiate constants
 * and add the self-contribution term.
 * Additional quantities such as velocity gradients will also get the final
 * terms added to them here.
 *
 * Also adds/multiplies the cosmological terms if need be.
 *
 * @param p The particle to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_density(
    struct part *p, const struct cosmology *cosmo) {
  
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension_plus_one(h_inv); /* 1/h^(d+1) */

  /* Finish the construction of the inverse of the c-matrix by
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 6; ++i) {
    p->mhd_data.dens.d_mat_inv.elements[i] *= h_inv_dim;
  }
  /* Finish the construction of the inverse of the A gradient
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 3; ++i) 
    for (int j = 0; j < 3; ++j) {
      p->mhd_data.dens.Mat_b[i][j] *= h_inv_dim;
      p->mhd_data.grad.Mat_da[i][j] *= h_inv_dim;
      }
  /* Invert the c-matrix */
  float d_mat_temp[3][3];
  get_matrix_from_sym_matrix(d_mat_temp, &p->mhd_data.dens.d_mat_inv);
  int res = invert_dimension_by_dimension_matrix(d_mat_temp);
  if (res) {
    sym_matrix_print(&p->mhd_data.dens.d_mat_inv);
    error("Error inverting D matrix");
  }
  const float g_b[3][3] = {
  {p->mhd_data.dens.Mat_b[0][0], p->mhd_data.dens.Mat_b[0][1], p->mhd_data.dens.Mat_b[0][2]},
  {p->mhd_data.dens.Mat_b[1][0], p->mhd_data.dens.Mat_b[1][1], p->mhd_data.dens.Mat_b[1][2]},
  {p->mhd_data.dens.Mat_b[2][0], p->mhd_data.dens.Mat_b[2][1], p->mhd_data.dens.Mat_b[2][2]}};
  
  const float g_da[3][3] = {
  {p->mhd_data.grad.Mat_da[0][0], p->mhd_data.grad.Mat_da[0][1], p->mhd_data.grad.Mat_da[0][2]},
  {p->mhd_data.grad.Mat_da[1][0], p->mhd_data.grad.Mat_da[1][1], p->mhd_data.grad.Mat_da[1][2]},
  {p->mhd_data.grad.Mat_da[2][0], p->mhd_data.grad.Mat_da[2][1], p->mhd_data.grad.Mat_da[2][2]}};
  
  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++){ 
         p->mhd_data.dens.Mat_b[i][j] = 0.f;
         p->mhd_data.grad.Mat_da[i][j] = 0.f;
      for (int k = 0; k < 3; k++){ 
         p->mhd_data.dens.Mat_b[i][j] += d_mat_temp[j][k] * g_b[i][k];
         p->mhd_data.grad.Mat_da[i][j] += d_mat_temp[j][k] * g_da[i][k];
	 }
  }
  
  for (int i = 0; i < 3; i++) 
     p->mhd_data.BPred[i] =
         p->mhd_data.dens.Mat_b[(i+2)%3][(i+1)%3] - p->mhd_data.dens.Mat_b[(i+1)%3][(i+2)%3];
  p->mhd_data.divA = p->mhd_data.dens.Mat_b[0][0] + p->mhd_data.dens.Mat_b[1][1] +
                     p->mhd_data.dens.Mat_b[2][2];

  get_sym_matrix_from_matrix(&p->mhd_data.grad.d_mat, d_mat_temp);

}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const float mu_0) {

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after mhd_prepare_gradient.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_reset_gradient(
    struct part *p) {

  zero_sym_matrix(&p->mhd_data.grad.c_mat_inv);
  for (int i = 0; i < 3; ++i) 
    for (int j = 0; j < 3; ++j){ 
      p->mhd_data.grad.Mat_b[i][j] = 0.f;
      p->mhd_data.grad.Mat_F[i][j] = 0.f;
      }
  /* Curl B*/
  for (int k = 0; k < 3; k++) p->mhd_data.curl_B[k] = 0.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->mhd_data.grad_B_tensor[i][j] = 0.0f;
    }
  }

  /* Initialise MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p);
  hydro_set_signal_velocity(p, 2.0f * cms);

  /* SPH error*/
  p->mhd_data.mean_SPH_err = 0.f;
  for (int k = 0; k < 3; k++) {
    p->mhd_data.mean_grad_SPH_err[k] = 0.f;
  }
}

/**
 * @brief Finishes the gradient calculation.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void mhd_end_gradient(
    struct part *p, const float mu_0) {

  /* Add self contribution */
  p->mhd_data.mean_SPH_err += p->mass * kernel_root;
  /* Finish SPH_1 calculation*/
  p->mhd_data.mean_SPH_err *= pow_dimension(1.f / (p->h)) / p->rho;
  
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Finish the construction of the inverse of the c-matrix by
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 6; ++i) {
    p->mhd_data.grad.c_mat_inv.elements[i] *= h_inv_dim;
  }
  /* Finish the construction of the inverse of the A gradient
   * multiplying in the factors of h coming from W */
  for (int i = 0; i < 3; ++i) 
    for (int j = 0; j < 3; ++j) {
      p->mhd_data.grad.Mat_b[i][j] *= h_inv_dim;
      p->mhd_data.grad.Mat_F[i][j] *= h_inv_dim;
      }
  /* Invert the c-matrix */
  float c_mat_temp[3][3];
  get_matrix_from_sym_matrix(c_mat_temp, &p->mhd_data.grad.c_mat_inv);
  int res = invert_dimension_by_dimension_matrix(c_mat_temp);
  if (res) {
    sym_matrix_print(&p->mhd_data.grad.c_mat_inv);
    error("Error inverting C matrix");
  }
  
  const float g_b[3][3] = {
  {p->mhd_data.grad.Mat_b[0][0], p->mhd_data.grad.Mat_b[0][1], p->mhd_data.grad.Mat_b[0][2]},
  {p->mhd_data.grad.Mat_b[1][0], p->mhd_data.grad.Mat_b[1][1], p->mhd_data.grad.Mat_b[1][2]},
  {p->mhd_data.grad.Mat_b[2][0], p->mhd_data.grad.Mat_b[2][1], p->mhd_data.grad.Mat_b[2][2]}};
  
  const float g_F[3][3] = {
  {p->mhd_data.grad.Mat_F[0][0], p->mhd_data.grad.Mat_F[0][1], p->mhd_data.grad.Mat_F[0][2]},
  {p->mhd_data.grad.Mat_F[1][0], p->mhd_data.grad.Mat_F[1][1], p->mhd_data.grad.Mat_F[1][2]},
  {p->mhd_data.grad.Mat_F[2][0], p->mhd_data.grad.Mat_F[2][1], p->mhd_data.grad.Mat_F[2][2]}};

  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++){ 
         p->mhd_data.grad.Mat_b[i][j] = 0.f;
         p->mhd_data.grad.Mat_F[i][j] = 0.f;
      for (int k = 0; k < 3; k++){ 
         p->mhd_data.grad.Mat_b[i][j] += c_mat_temp[j][k] * g_b[i][k];
         p->mhd_data.grad.Mat_F[i][j] += c_mat_temp[j][k] * g_F[i][k];
	 }
  }
  
  for (int i = 0; i < 3; i++) 
     p->mhd_data.JPred[i] =
         p->mhd_data.grad.Mat_b[(i+2)%3][(i+1)%3] - p->mhd_data.grad.Mat_b[(i+1)%3][(i+2)%3];
  p->mhd_data.divB = p->mhd_data.grad.Mat_b[0][0] + p->mhd_data.grad.Mat_b[1][1] +
                     p->mhd_data.grad.Mat_b[2][2];

  get_sym_matrix_from_matrix(&p->mhd_data.force.c_mat, c_mat_temp);
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * In the desperate case where a particle has no neighbours (likely because
 * of the h_max ceiling), set the particle fields to something sensible to avoid
 * NaNs in the next calculations.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_part_has_no_neighbours(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo) {

  p->mhd_data.divB = 0.0f;
  p->mhd_data.curl_B[0] = 0.0f;
  p->mhd_data.curl_B[1] = 0.0f;
  p->mhd_data.curl_B[2] = 0.0f;
}

/**
 * @brief Prepare a particle for the force calculation.
 *
 * This function is called in the ghost task to convert some quantities coming
 * from the density loop over neighbours into quantities ready to be used in the
 * force loop over neighbours. Quantities are typically read from the density
 * sub-structure and written to the force sub-structure.
 * Examples of calculations done here include the calculation of viscosity term
 * constants, thermal conduction terms, hydro conversions, etc.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 * @param cosmo The current cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 */
__attribute__((always_inline)) INLINE static void mhd_prepare_force(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float dt_alpha,
    const float mu_0) {

  const float h = p->h;
  float B[3];
  B[0] = p->mhd_data.BPred[0];
  B[1] = p->mhd_data.BPred[1];
  B[2] = p->mhd_data.BPred[2];

  const float B2 = B[0] * B[0] + B[1] * B[1] + B[2] * B[2];
  const float normB = sqrtf(B2);

  float grad_B_mean_square = 0.0f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      grad_B_mean_square +=
          p->mhd_data.grad_B_tensor[i][j] * p->mhd_data.grad_B_tensor[i][j];
    }
  }

  const float alpha_AR_max = 2.0;

  p->mhd_data.alpha_AR =
      normB ? fminf(alpha_AR_max, h * sqrtf(grad_B_mean_square) / normB) : 0.0f;

  /* Sets Induction equation */
  for (int i = 0; i < 3; i++) {
      p->mhd_data.dAdt[i] = - p->mhd_data.resistive_eta * p->mhd_data.JPred[i];
      for (int k = 0; k < 3; k++) 
         p->mhd_data.dAdt[i] += p->mhd_data.grad.Mat_da[k][i];
  }
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_reset_acceleration(
    struct part *restrict p) {

  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] = 0.0f;
  }
  /* Save induction sources*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.Adv_A_source[k] = 0.0f;
    p->mhd_data.Diff_A_source[k] = 0.0f;
    p->mhd_data.Delta_A[k] = 0.0f;
  }
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model
 */
__attribute__((always_inline)) INLINE static void mhd_reset_predicted_values(
    struct part *p, const struct xpart *xp, const struct cosmology *cosmo,
    const float mu_0) {

  /* Re-set the predicted magnetic flux densities */
  p->mhd_data.APred[0] = xp->mhd_data.Afull[0];
  p->mhd_data.APred[1] = xp->mhd_data.Afull[1];
  p->mhd_data.APred[2] = xp->mhd_data.Afull[2];

  p->mhd_data.Gau = xp->mhd_data.Gaufull;

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);

  /* Re-set MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p);
  hydro_set_signal_velocity(p, fmaxf(hydro_get_signal_velocity(p), 2.0f * cms));
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Note the different time-step sizes used for the different quantities as they
 * include cosmological factors.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 * @param mu_0 The vacuum magnetic permeability.
 */
__attribute__((always_inline)) INLINE static void mhd_predict_extra(
    struct part *p, const struct xpart *xp, const float dt_drift,
    const float dt_therm, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, const float mu_0) {

  /* Predict the VP magnetic field */
  p->mhd_data.APred[0] += p->mhd_data.dAdt[0] * dt_therm;
  p->mhd_data.APred[1] += p->mhd_data.dAdt[1] * dt_therm;
  p->mhd_data.APred[2] += p->mhd_data.dAdt[2] * dt_therm;

  p->mhd_data.Gau += p->mhd_data.Gau_dt * dt_therm;

  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);

  /* Initialise MHD signal velocity */
  const float cms = mhd_get_comoving_magnetosonic_speed(p);
  hydro_set_signal_velocity(p, fmaxf(hydro_get_signal_velocity(p), 2.0f * cms));
}

/**
 * @brief Finishes the force calculation.
 *
 * Multiplies the force and accelerations by the appropiate constants
 * and add the self-contribution term. In most cases, there is little
 * to do here.
 *
 * Cosmological terms are also added/multiplied here.
 *
 * @param p The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void mhd_end_force(
    struct part *p, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float mu_0) {

  /* Get time derivative of Dedner scalar */
  p->mhd_data.Gau_dt = mhd_get_dGau_dt(p, cosmo, mu_0);

  /* Hubble expansion contribution to induction equation */
  float a_fac = (2.f + mhd_comoving_factor) * cosmo->a * cosmo->a * cosmo->H;
  p->mhd_data.dAdt[0] -= a_fac * p->mhd_data.APred[0];
  p->mhd_data.dAdt[1] -= a_fac * p->mhd_data.APred[1];
  p->mhd_data.dAdt[2] -= a_fac * p->mhd_data.APred[2];
  /*
  for (int i = 0; i < 3; i++) 
      for (int k = 0; k < 3; k++) 
         p->a_hydro[i] += p->mhd_data.grad.Mat_F[i][k];
 */ 
  /* Save forces*/
  for (int k = 0; k < 3; k++) {
    p->mhd_data.tot_mag_F[k] *= p->mass;
  }
  for (int i = 0; i < 3; i++) {
    p->mhd_data.tot_mag_F[i] = 0.f;
      for (int k = 0; k < 3; k++) 
         p->mhd_data.tot_mag_F[i] -= p->mhd_data.grad.Mat_F[i][k];
  }
}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantities are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void mhd_kick_extra(
    struct part *p, struct xpart *xp, const float dt_therm, const float dt_grav,
    const float dt_hydro, const float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {

  /* Integrate the magnetic flux density forward in time */
  xp->mhd_data.Afull[0] += p->mhd_data.dAdt[0] * dt_therm;
  xp->mhd_data.Afull[1] += p->mhd_data.dAdt[1] * dt_therm;
  xp->mhd_data.Afull[2] += p->mhd_data.dAdt[2] * dt_therm;

  xp->mhd_data.Gaufull += p->mhd_data.Gau_dt * dt_therm;
}

/**
 * @brief Converts MHD quantities of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy in the case
 * of hydro for instance.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 */
__attribute__((always_inline)) INLINE static void mhd_convert_quantities(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props, const float mu_0) {

  const float a_fact = pow(cosmo->a, -mhd_comoving_factor - 1.f);
  /* Set Restitivity Eta */
  p->mhd_data.resistive_eta = hydro_props->mhd.mhd_eta;

  p->mhd_data.APred[0] *= a_fact;
  p->mhd_data.APred[1] *= a_fact;
  p->mhd_data.APred[2] *= a_fact;
  /* Full Step */
  xp->mhd_data.Afull[0] = p->mhd_data.APred[0];
  xp->mhd_data.Afull[1] = p->mhd_data.APred[1];
  xp->mhd_data.Afull[2] = p->mhd_data.APred[2];

  /* Instantiate Alfven speed */
  p->mhd_data.Alfven_speed = mhd_get_comoving_Alfven_speed(p, mu_0);
}

/**
 * @brief Initialises the particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void mhd_first_init_part(
    struct part *restrict p, struct xpart *restrict xp,
    const struct mhd_global_data *mhd_data, const double Lsize) {
  xp->mhd_data.Gaufull = 0.0f;
  p->mhd_data.Gau = 0.0f;
  p->mhd_data.divB = 0.0f;

  mhd_reset_acceleration(p);
  mhd_init_part(p);
}

#endif /* SWIFT_VECTOR_POTENTIAL_MHD_H */
