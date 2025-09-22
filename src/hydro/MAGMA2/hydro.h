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
#ifndef SWIFT_MAGMA2_HYDRO_H
#define SWIFT_MAGMA2_HYDRO_H

/**
 * @file MAGMA2/hydro.h
 * @brief Density-Energy non-conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (Non-neighbour loop
 *        equations)
 */

#include "adiabatic_index.h"
#include "approx_math.h"
#include "cosmology.h"
#include "dimension.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "hydro_properties.h"
#include "hydro_space.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "pressure_floor.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <float.h>

/**
 * @brief Returns the comoving internal energy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not internal energy, this function computes the internal
 * energy from the thermodynamic variable.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp) {
  return xp->u_full;
}

/**
 * @brief Returns the physical internal energy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not internal energy, this function computes the internal
 * energy from the thermodynamic variable and converts it to
 * physical coordinates.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy(const struct part *restrict p,
                                   const struct xpart *restrict xp,
                                   const struct cosmology *cosmo) {
  return xp->u_full * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_internal_energy(const struct part *restrict p) {
  return p->u;
}

/**
 * @brief Returns the physical internal energy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_internal_energy(const struct part *restrict p,
                                           const struct cosmology *cosmo) {
  return p->u * cosmo->a_factor_internal_energy;
}

/**
 * @brief Returns the comoving pressure of a particle
 *
 * Computes the pressure based on the particle's properties.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_pressure(
    const struct part *restrict p) {
  return gas_pressure_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the physical pressure of a particle
 *
 * Computes the pressure based on the particle's properties and
 * convert it to physical coordinates.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_pressure(
    const struct part *restrict p, const struct cosmology *cosmo) {
  return cosmo->a_factor_pressure * hydro_get_comoving_pressure(p);
}

/**
 * @brief Returns the comoving entropy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not entropy, this function computes the entropy from
 * the thermodynamic variable.
 *
 * @param p The particle of interest
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part *restrict p, const struct xpart *restrict xp) {
  return gas_entropy_from_internal_energy(p->rho, xp->u_full);
}

/**
 * @brief Returns the physical entropy of a particle at the last
 * time the particle was kicked.
 *
 * For implementations where the main thermodynamic variable
 * is not entropy, this function computes the entropy from
 * the thermodynamic variable and converts it to
 * physical coordinates.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_entropy(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo) {
  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(p->rho, xp->u_full);
}

/**
 * @brief Returns the comoving entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_comoving_entropy(const struct part *restrict p) {
  return gas_entropy_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the physical entropy of a particle drifted to the
 * current time.
 *
 * @param p The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_drifted_physical_entropy(const struct part *restrict p,
                                   const struct cosmology *cosmo) {
  /* Note: no cosmological conversion required here with our choice of
   * coordinates. */
  return gas_entropy_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the comoving sound speed of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_soundspeed(const struct part *restrict p) {
  return gas_soundspeed_from_internal_energy(p->rho, p->u);
}

/**
 * @brief Returns the physical sound speed of a particle
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_soundspeed(const struct part *restrict p,
                              const struct cosmology *cosmo) {
  return cosmo->a_factor_sound_speed * hydro_get_comoving_soundspeed(p);
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_density(
    const struct part *restrict p) {
  return p->rho;
}

/**
 * @brief Returns the comoving density of a particle.
 *
 * @param p The particle of interest
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_density(
    const struct part *restrict p, const struct cosmology *cosmo) {
  return cosmo->a3_inv * p->rho;
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part *restrict p) {
  return p->mass;
}

/**
 * @brief Sets the mass of a particle
 *
 * @param p The particle of interest
 * @param m The mass to set.
 */
__attribute__((always_inline)) INLINE static void hydro_set_mass(
    struct part *restrict p, float m) {
  p->mass = m;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy_dt(const struct part *restrict p) {
  return p->u_dt;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest
 * @param cosmo Cosmology data structure
 */
__attribute__((always_inline)) INLINE static float
hydro_get_physical_internal_energy_dt(const struct part *restrict p,
                                      const struct cosmology *cosmo) {
  return p->u_dt * cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_comoving_internal_energy_dt(struct part *restrict p, float du_dt) {
  p->u_dt = du_dt;
}

/**
 * @brief Returns the time derivative of internal energy of a particle
 *
 * We assume a constant density.
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param du_dt The new time derivative of the internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy_dt(struct part *restrict p,
                                      const struct cosmology *cosmo,
                                      float du_dt) {
  p->u_dt = du_dt / cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the physical entropy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param entropy The physical entropy
 */
__attribute__((always_inline)) INLINE static void hydro_set_physical_entropy(
    struct part *p, struct xpart *xp, const struct cosmology *cosmo,
    const float entropy) {
  /* Note there is no conversion from physical to comoving entropy */
  const float comoving_entropy = entropy;
  xp->u_full = gas_internal_energy_from_entropy(p->rho, comoving_entropy);
}

/**
 * @brief Sets the physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param xp The extended particle data.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_physical_internal_energy(struct part *p, struct xpart *xp,
                                   const struct cosmology *cosmo,
                                   const float u) {
  xp->u_full = u / cosmo->a_factor_internal_energy;
}

/**
 * @brief Sets the drifted physical internal energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param pressure_floor The properties of the pressure floor.
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_set_drifted_physical_internal_energy(
    struct part *p, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, const float u) {
  /* There is no need to use the floor here as this function is called in the
   * feedback, so the new value of the internal energy should be strictly
   * higher than the old value. */

  p->u = u / cosmo->a_factor_internal_energy;

  /* Now recompute the extra quantities */

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.soundspeed = soundspeed;
  p->force.pressure = pressure_including_floor;

  /* Signal speed */
  const float v_sig = const_viscosity_alpha_prefactor * soundspeed;

  p->dt_min = min(p->dt_min, p->h_min / v_sig);
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
hydro_set_v_sig_based_on_velocity_kick(struct part *p,
                                       const struct cosmology *cosmo,
                                       const float dv_phys) {
  /* Compute the velocity kick in comoving coordinates */
  const float dv = dv_phys / cosmo->a_factor_sound_speed;

  /* Sound speed */
  const float soundspeed = hydro_get_comoving_soundspeed(p);

  /* Signal speed */
  const float v_sig_sound = const_viscosity_alpha_prefactor * soundspeed;
  const float v_sig_kick = const_viscosity_beta_prefactor * dv;
  const float v_sig = v_sig_sound + v_sig_kick;

  p->dt_min = min(p->dt_min, p->h_min / v_sig);
}

/**
 * @brief Update the value of the viscosity alpha for the scheme.
 *
 * @param p the particle of interest
 * @param alpha the new value for the viscosity coefficient.
 */
__attribute__((always_inline)) INLINE static void hydro_set_viscosity_alpha(
    struct part *restrict p, float alpha) { }

/**
 * @brief Update the value of the diffusive coefficients to the
 *        feedback reset value for the scheme.
 *
 * @param p the particle of interest
 */
__attribute__((always_inline)) INLINE static void
hydro_diffusive_feedback_reset(struct part *restrict p) {
  /* Set the viscosity to the max, and the diffusion to the min */
  hydro_set_viscosity_alpha(p, const_viscosity_alpha);
}

/**
 * @brief Computes the hydro time-step of a given particle
 *
 * This function returns the time-step of a particle given its hydro-dynamical
 * state. A typical time-step calculation would be the use of the CFL condition.
 *
 * @param p Pointer to the particle data
 * @param xp Pointer to the extended particle data
 * @param hydro_properties The SPH parameters
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float hydro_compute_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct hydro_props *restrict hydro_properties,
    const struct cosmology *restrict cosmo) {

  if (p->dt_min == 0.f) return FLT_MAX;

  /* Use full kernel support and physical time */
  const float conv = kernel_gamma * cosmo->a / cosmo->a_factor_sound_speed;

  /* CFL condition */
  const float dt_cfl = 2. * hydro_properties->CFL_condition * conv * p->dt_min;

  /* Do not allow more than 0.25 * |u|/|du/dt| per step */
  const float dt_u = 
      (p->u_dt_cond != 0.) ? 0.25 * p->u / fabs(p->u_dt_cond) : FLT_MAX;

#ifdef MAGMA2_DEBUG_CHECKS
  if (dt_u < dt_cfl) {
    message("dt_u < dt_cfl for pid=%lld u=%g u_dt_cond=%g dt_min=%g conv=%g "
            "dt_cfl=%g dt_u=%g",
            p->id, p->u, p->u_dt_cond, p->dt_min, conv, dt_cfl, dt_u);
  }
#endif

  return fmin(dt_cfl, dt_u);
}

/**
 * @brief Compute the signal velocity between two gas particles
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float hydro_signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta) {

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float c_ij = 0.5f * (ci + cj);

  const float v_sig_alpha = const_viscosity_alpha_prefactor * c_ij;
  const float v_sig_beta = const_viscosity_beta_prefactor * mu_ij;
  const float v_sig = v_sig_alpha - fmin(v_sig_beta, 0.f);

  return v_sig;
}

/**
 * @brief Returns the physical velocity divergence.
 *
 * @brief p The particle
 */
__attribute__((always_inline)) INLINE static float hydro_get_div_v(
    const struct part *restrict p) {

  return p->gradients.velocity_tensor[0][0] +
         p->gradients.velocity_tensor[1][1] +
         p->gradients.velocity_tensor[2][2];
}

/**
 * @brief Returns the physical velocity divergence.
 *
 * @brief p The particle
 */
__attribute__((always_inline)) INLINE static float hydro_get_physical_div_v(
    const struct part *restrict p, const struct cosmology* cosmo) {

  const float div_v = hydro_get_div_v(p);
  return div_v * cosmo->a2_inv + hydro_dimension * cosmo->H;
}

/**
 * @brief Does some extra hydro operations once the actual physical time step
 * for the particle is known.
 *
 * @param p The particle to act upon.
 * @param dt Physical time step of the particle during the next step.
 */
__attribute__((always_inline)) INLINE static void hydro_timestep_extra(
    struct part *p, float dt) {}

/**
 * @brief Operations performed when a particle gets removed from the
 * simulation volume.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void hydro_remove_part(
    const struct part *p, const struct xpart *xp, const double time) {}

/**
 * @brief Prepares a particle for the density calculation.
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the various density loop over neighbours. Typically, all fields of the
 * density sub-structure of a particle get zeroed in here.
 *
 * @param p The particle to act upon
 * @param hs #hydro_space containing hydro specific space information.
 */
__attribute__((always_inline)) INLINE static void hydro_init_part(
    struct part *restrict p, const struct hydro_space *hs) {
  p->density.wcount = 0.f;
  p->density.wcount_dh = 0.f;
  p->rho = 0.f;
  p->density.rho_dh = 0.f;
#ifdef MAGMA2_DEBUG_CHECKS
  p->debug.num_ngb = 0;
  p->debug.v_sig_visc_max = 0.;
  p->debug.v_sig_cond_max = 0.;
#endif

#ifdef hydro_props_use_adiabatic_correction
  p->gradients.adiabatic_f_numerator = 0.;
#endif

  p->gradients.wcount = 0.;
  p->gradients.u_well_conditioned = 1;
  p->gradients.D_well_conditioned = 1;

  p->gradients.du_min = FLT_MAX;
  p->gradients.du_max = -FLT_MAX;
  p->gradients.kernel_size = FLT_MIN;

  /* These must be zeroed before the density loop */
  for (int i = 0; i < 3; i++) {
    p->gradients.u_aux[i] = 0.;
    p->gradients.u_aux_norm[i] = 0.;

    p->gradients.dv_min[i] = FLT_MAX;
    p->gradients.dv_max[i] = -FLT_MAX;

    for (int j = 0; j < 3; j++) {
      p->gradients.velocity_tensor_aux[i][j] = 0.;
      p->gradients.velocity_tensor_aux_norm[i][j] = 0.;
    }
  }
}

/**
 * @brief Computes the condition number of a matrix A
 *
 *
 * @param A The matrix to compute the condition number of.
 */
__attribute__((always_inline)) INLINE static 
double condition_number(gsl_matrix *A) {

  gsl_matrix *A_copy = gsl_matrix_alloc(3, 3);
  gsl_matrix_memcpy(A_copy, A);

  gsl_vector *S = gsl_vector_alloc(3);
  gsl_vector *work = gsl_vector_alloc(3);
  gsl_matrix *V = gsl_matrix_alloc(3, 3);

  gsl_linalg_SV_decomp(A_copy, V, S, work);

  double s_max = gsl_vector_get(S, 0);
  double s_min = gsl_vector_get(S, 2);

  gsl_matrix_free(A_copy);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);

  return (s_min != 0.) ? s_max / s_min : const_condition_number_upper_limit;
}

/**
 * @brief Vector dot product of two 3D vectors.
 *
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @param result The result of the dot product.
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_vec3_vec3_dot(const hydro_real_t *restrict vec_a, 
                                 const hydro_real_t *restrict vec_b) {

  return vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2];
}

/**
 * @brief Norm of two 3D vectors.
 *
 *
 * @param vec The vector.
 * @param result The result of the norm.
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_vec3_norm(const hydro_real_t *restrict vec) {

  return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

/**
 * @brief Unit vector of the given vector.
 *
 *
 * @param vec The vector.
 * @param result The unit vector of vec
 */
__attribute__((always_inline)) INLINE static 
void hydro_vec3_unit(const hydro_real_t *restrict vec,
                     hydro_real_t *restrict result) {

  result[0] = 0.;
  result[1] = 0.;
  result[2] = 0.;

  const hydro_real_t vec_norm = hydro_vec3_norm(vec);
  if (vec_norm > 0.) {
    result[0] = vec[0] / vec_norm;
    result[1] = vec[1] / vec_norm;
    result[2] = vec[2] / vec_norm;
  }
}

/**
 * @brief The Frobenius inner product of two matrices.
 *
 *
 * @param mat The matrix to contract with the vector.
 * @param vec The vector to contract with the matrix.
 * @param result The result of the contraction.
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_mat3x3_mat3x3_dot(const hydro_real_t (*restrict mat_a)[3], 
                                     const hydro_real_t (*restrict mat_b)[3]) {

  return mat_a[0][0] * mat_b[0][0] +
         mat_a[0][1] * mat_b[0][1] +
         mat_a[0][2] * mat_b[0][2] +
         mat_a[1][0] * mat_b[1][0] +
         mat_a[1][1] * mat_b[1][1] +
         mat_a[1][2] * mat_b[1][2] +
         mat_a[2][0] * mat_b[2][0] +
         mat_a[2][1] * mat_b[2][1] +
         mat_a[2][2] * mat_b[2][2];
}

/**
 * @brief Contracts the last index of matrix mat with a vector vec and stores in
 *        result.
 *
 *
 * @param mat The matrix to contract with the vector.
 * @param vec The vector to contract with the matrix.
 * @param result The result of the contraction.
 */
__attribute__((always_inline)) INLINE static void hydro_mat3x3_vec3_dot(
    const hydro_real_t (*restrict mat)[3], 
    const hydro_real_t *restrict vec, 
    hydro_real_t *restrict result) {

  result[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2];
  result[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2];
  result[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2];
}

/**
 * @brief Contracts the last two indices of the tensor tensor with a matrix
 *        mat and stores in result. Form: mat^T * tensor * mat.
 *
 *
 * @param tensor The tensor to contract with the matrix.
 * @param mat The matrix to contract with the tensor.
 * @param result The result of the contraction.
 */
__attribute__((always_inline)) INLINE static 
void hydro_tensor3x3x3_matrix3x3_dot(
    const hydro_real_t (*restrict tensor)[3][3],
    const hydro_real_t (*restrict mat)[3],
    hydro_real_t *restrict result) {

  result[0] = tensor[0][0][0] * mat[0][0] +
              tensor[0][0][1] * mat[0][1] +
              tensor[0][0][2] * mat[0][2] +
              tensor[0][1][0] * mat[1][0] +
              tensor[0][1][1] * mat[1][1] +
              tensor[0][1][2] * mat[1][2] +
              tensor[0][2][0] * mat[2][0] +
              tensor[0][2][1] * mat[2][1] +
              tensor[0][2][2] * mat[2][2];
  
  result[1] = tensor[1][0][0] * mat[0][0] +
              tensor[1][0][1] * mat[0][1] +
              tensor[1][0][2] * mat[0][2] +
              tensor[1][1][0] * mat[1][0] +
              tensor[1][1][1] * mat[1][1] +
              tensor[1][1][2] * mat[1][2] +
              tensor[1][2][0] * mat[2][0] +
              tensor[1][2][1] * mat[2][1] +
              tensor[1][2][2] * mat[2][2];

  result[2] = tensor[2][0][0] * mat[0][0] +
              tensor[2][0][1] * mat[0][1] +
              tensor[2][0][2] * mat[0][2] +
              tensor[2][1][0] * mat[1][0] +
              tensor[2][1][1] * mat[1][1] +
              tensor[2][1][2] * mat[1][2] +
              tensor[2][2][0] * mat[2][0] +
              tensor[2][2][1] * mat[2][1] +
              tensor[2][2][2] * mat[2][2];
}

/**
 * @brief Constructs the outer product of two 3D vectors.
 *
 *
 * @param vec_a The first vector.
 * @param vec_b The second vector.
 * @param result The result of the outer product.
 */
__attribute__((always_inline)) INLINE static void hydro_vec3_vec3_outer(
    const hydro_real_t *restrict vec_a, 
    const hydro_real_t *restrict vec_b, 
    hydro_real_t (*restrict result)[3]) {

  result[0][0] = vec_a[0] * vec_b[0];
  result[0][1] = vec_a[0] * vec_b[1];
  result[0][2] = vec_a[0] * vec_b[2];

  result[1][0] = vec_a[1] * vec_b[0];
  result[1][1] = vec_a[1] * vec_b[1];
  result[1][2] = vec_a[1] * vec_b[2];

  result[2][0] = vec_a[2] * vec_b[0];
  result[2][1] = vec_a[2] * vec_b[1];
  result[2][2] = vec_a[2] * vec_b[2];
}

/**
 * @brief Limit gradients to variation across the kernel.
 *
 *
 * @param df_raw Raw value
 * @param df_reconstructed Reconstructed value
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_scalar_minmod(const hydro_real_t df_raw, 
                                 const hydro_real_t df_reconstructed) {

  if (df_raw * df_reconstructed <= 0.) return 0.;

  return (fabs(df_raw) < fabs(df_reconstructed)) ? df_raw : df_reconstructed;
}

/**
 * @brief Limit gradients to variation across the kernel.
 *
 *
 * @param df_reconstructed Reconstructed estimate of df
 * @param df_raw Particle estimate of df
 * @param df_min_i Minimum value of df_raw across the kernel.
 * @param df_max_i Maximum value of df_raw across the kernel.
 * @param df_min_j Minimum value of df_raw across the kernel.
 * @param df_max_j Maximum value of df_raw across the kernel
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_scalar_minmod_limiter(const hydro_real_t df_reconstructed,
                                         const hydro_real_t df_raw,
                                         const hydro_real_t df_min_i,
                                         const hydro_real_t df_max_i,
                                         const hydro_real_t df_min_j,
                                         const hydro_real_t df_max_j) {

#ifdef hydro_props_use_strict_minmod_limiter
  hydro_real_t df = hydro_scalar_minmod(df_raw, df_reconstructed);

  const hydro_real_t lo = fmax(df_min_i, -df_max_j);
  const hydro_real_t hi = fmin(df_max_i, -df_min_j);

  if (lo > hi) return df_raw;

  if (df < lo) df = lo;
  if (df > hi) df = hi;

  return df;
#else
  return df_reconstructed;
#endif
}

/**
 * @brief Limit variations across the kernel (for 3-vectors)
 *
 *
 * @param dv_reconstructed Reconstructed estimate of df
 * @param dv_raw Particle estimate of df
 * @param dv_min_i Minimum value of df_raw across the kernel.
 * @param dv_max_i Maximum value of df_raw across the kernel.
 * @param dv_min_j Minimum value of df_raw across the kernel.
 * @param dv_max_j Maximum value of df_raw across the kernel
 * @param dv_ij The vector to return.
 */
__attribute__((always_inline)) INLINE static 
void hydro_vec_minmod_limiter(const hydro_real_t *restrict dv_reconstructed,
                              const hydro_real_t *restrict dv_raw,
                              const hydro_real_t *restrict dv_min_i,
                              const hydro_real_t *restrict dv_max_i,
                              const hydro_real_t *restrict dv_min_j,
                              const hydro_real_t *restrict dv_max_j,
                              hydro_real_t *restrict dv_ij) {

  dv_ij[0] = hydro_scalar_minmod_limiter(dv_reconstructed[0], dv_raw[0],
                                         dv_min_i[0], dv_max_i[0], 
                                         dv_min_j[0], dv_max_j[0]);
  dv_ij[1] = hydro_scalar_minmod_limiter(dv_reconstructed[1], dv_raw[1],
                                         dv_min_i[1], dv_max_i[1], 
                                         dv_min_j[1], dv_max_j[1]);
  dv_ij[2] = hydro_scalar_minmod_limiter(dv_reconstructed[2], dv_raw[2],
                                         dv_min_i[2], dv_max_i[2], 
                                         dv_min_j[2], dv_max_j[2]);

  /* If any of the components are exactly zero, we return a zero vector */
  if (dv_ij[0] == 0. || dv_ij[1] == 0. || dv_ij[2] == 0.) {
    dv_ij[0] = 0.;
    dv_ij[1] = 0.;
    dv_ij[2] = 0.;
  }
}

/**
 * @brief Limit gradients to variation across the kernel.
 *
 *
 * @param df_min Minimum value of delta f across the kernel.
 * @param df_max Maximum value of delta f across the kernel.
 * @param kernel_size Interaction distance across the kernel.
 * @param grad The vector gradient to slope limit.
 */
__attribute__((always_inline)) INLINE static
void hydro_vec_slope_limiter(const hydro_real_t df_min,
                             const hydro_real_t df_max,
                             const hydro_real_t kernel_size,
                             hydro_real_t *restrict grad) {

#ifdef hydro_props_use_extra_slope_limiter
  const hydro_real_t grad_norm = hydro_vec3_norm(grad);
  const hydro_real_t length = const_grad_overshoot_length * kernel_size;

  /* Nothing to do if there is no gradient or no look-ahead distance */
  if (grad_norm > 0. && length > 0.) {
    const hydro_real_t df_min_abs = (df_min < 0.) ? fabs(df_min) : 0.;
    const hydro_real_t df_max_abs = (df_max > 0.) ? fabs(df_max) : 0.;
    const hydro_real_t df_abs_min = fmin(df_min_abs, df_max_abs);

    hydro_real_t bound = df_abs_min;
    
#ifdef const_grad_overshoot_tolerance
    const hydro_real_t tolerance = const_grad_overshoot_tolerance;
    if (tolerance > 0.) {
      const hydro_real_t df_abs_max = fmax(df_min_abs, df_max_abs);
      const hydro_real_t extra = tolerance * df_abs_max;
      const hydro_real_t cap = 
          (df_abs_min + extra < df_abs_max) ? df_abs_min + extra : df_abs_max;
      bound = cap;
    }
#endif

    const hydro_real_t limiter = 
        (bound > 0.) ? (bound / (length * grad_norm)) : 0.;

    if (limiter < 1.) {
      grad[0] *= limiter;
      grad[1] *= limiter;
      grad[2] *= limiter;
    }
  }
#endif

}

/**
 * @brief Computes Phi_ab from Rosswog 2020 21 for a field given A. This is the 
 * van Leer 1974 slope limiting procedure.
 *
 *
 * @param A_ij The ratio of the gradients of the two particles.
 * @param eta_i The normed smoothing length of the first particle.
 * @param eta_j The normed smoothing length of the second particle.
 */
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_van_leer_phi(const hydro_real_t A_ij,
                                const hydro_real_t eta_i, 
                                const hydro_real_t eta_j) {

  hydro_real_t phi_raw = (4. * A_ij) / ((1. + A_ij) * (1. + A_ij));
  phi_raw = fmin(1., phi_raw);
  phi_raw = fmax(0., phi_raw);

  /* η_ab and η_crit damping */
  const hydro_real_t eta_ij = min(eta_i, eta_j);

  hydro_real_t damping = 1.;
  if (eta_ij <= const_slope_limiter_eta_crit) {
    const hydro_real_t diff = 
        (eta_ij - const_slope_limiter_eta_crit) / const_slope_limiter_eta_fold;
    damping = exp(-diff * diff);
  }

  phi_raw *= damping;

  /* Handle any edge cases */
  phi_raw = fmin(1., phi_raw);
  phi_raw = fmax(0., phi_raw);

  return phi_raw;
}

/**
 * @brief Computes A_ab from Rosswog 2020 Eq. 22 for a scalar field.
 *
 *
 * @param grad_i The gradient of the quantity for the first particle.
 * @param grad_j The gradient of the quantity for the second particle.
 * @param dx The distance vector between the two particles ( ri - rj ).
 */
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_scalar_van_leer_A(const hydro_real_t *restrict grad_i,
                                     const hydro_real_t *restrict grad_j,
                                     const hydro_real_t *restrict dx) {

  const hydro_real_t grad_dot_x_i = hydro_vec3_vec3_dot(grad_i, dx);
  const hydro_real_t grad_dot_x_j = hydro_vec3_vec3_dot(grad_j, dx);

  /* Fall-back to raw estimates */
  if (grad_dot_x_i * grad_dot_x_j <= 0.) return 0.;

  /* Regularize denominator */
  if (fabs(grad_dot_x_j) < 1.e-10) return 0.;

  const hydro_real_t A_ij = grad_dot_x_i / grad_dot_x_j;

  return A_ij;
}

/**
 * @brief Computes A_ab from Rosswog 2020 Eq. 22 for a vector field.
 *
 *
 * @param grad_i The gradient tensor for the first particle.
 * @param grad_j The gradient tensor for the second particle.
 * @param dx The distance vector between the two particles ( ri - rj ).
 */
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_vector_van_leer_A(const hydro_real_t (*restrict grad_i)[3],
                                     const hydro_real_t (*restrict grad_j)[3],
                                     const hydro_real_t *restrict dx) {

  hydro_real_t delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(dx, dx, delta_ij);

  const hydro_real_t grad_dot_x_x_i = hydro_mat3x3_mat3x3_dot(grad_i, delta_ij);
  const hydro_real_t grad_dot_x_x_j = hydro_mat3x3_mat3x3_dot(grad_j, delta_ij);

  /* Fall-back to raw estimates */
  if (grad_dot_x_x_i * grad_dot_x_x_j <= 0.) return 0.;

  /* Regularize denominator */
  if (fabs(grad_dot_x_x_j) < 1.e-10) return 0.;

  const hydro_real_t A_ij = grad_dot_x_x_i / grad_dot_x_x_j;

  return A_ij;
}

/**
 * @brief Computes Phi_ab from Rosswog 2020 21 for a scalar field. This is the 
 * van Leer 1974 slope limiting procedure.
 *
 *
 * @param grad_i The gradient of the quantity for the first particle.
 * @param grad_j The gradient of the quantity for the second particle.
 * @param dx The distance vector between the two particles ( ri - rj ).
 * @param eta_i The normed smoothing length of the first particle.
 * @param eta_j The normed smoothing length of the second particle.
 */
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_scalar_van_leer_phi(const hydro_real_t *restrict grad_i,
                                       const hydro_real_t *restrict grad_j,
                                       const hydro_real_t *restrict dx,
                                       const hydro_real_t eta_i, 
                                       const hydro_real_t eta_j) {

  const hydro_real_t A_ij = hydro_scalar_van_leer_A(grad_i, grad_j, dx);

  return hydro_van_leer_phi(A_ij, eta_i, eta_j);
}

/**
 * @brief Computes Phi_ab from Rosswog 2020 21 for a vector field. This is the 
 * van Leer 1974 slope limiting procedure.
 *
 *
 * @param grad_i The gradient of the quantity for the first particle.
 * @param grad_j The gradient of the quantity for the second particle.
 * @param dx The distance vector between the two particles ( ri - rj ).
 * @param eta_i The normed smoothing length of the first particle.
 * @param eta_j The normed smoothing length of the second particle.
 */
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_vector_van_leer_phi(const hydro_real_t (*restrict grad_i)[3],
                                       const hydro_real_t (*restrict grad_j)[3],
                                       const hydro_real_t *restrict dx,
                                       const hydro_real_t eta_i, 
                                       const hydro_real_t eta_j) {

  const hydro_real_t A_ij = hydro_vector_van_leer_A(grad_i, grad_j, dx);

  return hydro_van_leer_phi(A_ij, eta_i, eta_j);
}

/**
 * @brief Reconstructs the scalar field at the mid-point between to the 
 *        two particles to second order.
 *
 *
 * @param phi The slope limiter value from the van Leer 1974 scheme.
 * @param dx The distance vector between the two particles ( ri - rj ).
 * @param f The field at the particle.
 * @param grad The gradient of the field at the particle.
 * @param hess The Hessian of the field at the particle.
 * @param f_reconstructed The reconstructed field at the mid-point (2nd order).
 */
__attribute__((always_inline)) INLINE static void 
hydro_scalar_second_order_reconstruction(const hydro_real_t phi,
                                         const hydro_real_t *restrict dx,
                                         const hydro_real_t f, 
                                         const hydro_real_t *restrict grad,
                                         const hydro_real_t (*restrict hess)[3],
                                         hydro_real_t *f_reconstructed) {

  /* Midpoint from Equation 17 Rosswog 2020 */
  const hydro_real_t midpoint[3] = {0.5 * dx[0],
                                    0.5 * dx[1],
                                    0.5 * dx[2]};
  hydro_real_t delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(midpoint, midpoint, delta_ij);

  const hydro_real_t df_dx_dot_r = hydro_vec3_vec3_dot(grad, midpoint);
  const hydro_real_t df_dx_dx_dot_r2 = hydro_mat3x3_mat3x3_dot(hess, delta_ij);

  /* Apply limited slope reconstruction */
  *f_reconstructed = f + phi * (df_dx_dot_r + 0.5 * df_dx_dx_dot_r2);
}

/**
 * @brief Reconstructs the vector field at the mid-point between to the 
 *        two particles to second order.
 *
 *
 * @param phi The slope limiter value from the van Leer 1974 scheme.
 * @param dx The distance vector between the two particles ( ri - rj ).
 * @param f The vector field at the particle.
 * @param grad The gradient of the vector field at the particle.
 * @param hess The Hessian of the vector field at the particle.
 * @param f_reconstructed The reconstructed vector field at the mid-point (2nd order).
 */
__attribute__((always_inline)) INLINE static void 
hydro_vector_second_order_reconstruction(const hydro_real_t phi,
                                         const hydro_real_t *restrict dx,
                                         const hydro_real_t *restrict f, 
                                         const hydro_real_t (*restrict grad)[3],
                                         const hydro_real_t (*restrict hess)[3][3],
                                         hydro_real_t *restrict f_reconstructed) {

  /* Midpoint from Equation 17 Rosswog 2020 */
  const hydro_real_t midpoint[3] = {0.5 * dx[0],
                                    0.5 * dx[1],
                                    0.5 * dx[2]};
  hydro_real_t delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(midpoint, midpoint, delta_ij);

  hydro_real_t df_dx_dot_r[3] = {0};
  hydro_real_t df_dx_dx_dot_r2[3] = {0};

  hydro_mat3x3_vec3_dot(grad, midpoint, df_dx_dot_r);
  hydro_tensor3x3x3_matrix3x3_dot(hess, delta_ij, df_dx_dx_dot_r2);

  /* Apply limited slope reconstruction */
  f_reconstructed[0] = f[0] + phi * (df_dx_dot_r[0] + 0.5 * df_dx_dx_dot_r2[0]);
  f_reconstructed[1] = f[1] + phi * (df_dx_dot_r[1] + 0.5 * df_dx_dx_dot_r2[1]);
  f_reconstructed[2] = f[2] + phi * (df_dx_dot_r[2] + 0.5 * df_dx_dx_dot_r2[2]);
}

/**
 * @brief Get the balsara limiter for viscosity.
 *
 *
 * @param p Particle p
 * @param cosmo The cosmology structure
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_get_balsara_limiter(const struct part *restrict p,
                                       const struct cosmology* cosmo) {
  
#ifdef hydro_props_use_balsara_limiter
  const hydro_real_t fac_B = cosmo->a_factor_Balsara_eps;
  hydro_real_t balsara = 1.;

  /* Can't trust velocity_tensor when having ill-conditioned matrices */
  if (p->gradients.C_well_conditioned && p->gradients.D_well_conditioned) {
    const hydro_real_t div_v_phys = hydro_get_physical_div_v(p, cosmo);

    const hydro_real_t curl_v[3] = {
        p->gradients.velocity_tensor[2][1] - p->gradients.velocity_tensor[1][2],
        p->gradients.velocity_tensor[0][2] - p->gradients.velocity_tensor[2][0],
        p->gradients.velocity_tensor[1][0] - p->gradients.velocity_tensor[0][1]
    };

    const hydro_real_t curl_v_norm_phys = 
        hydro_vec3_norm(curl_v) * cosmo->a2_inv;
    const hydro_real_t Balsara_eps = 1.e-4 * fac_B * p->force.soundspeed / p->h;
    balsara = 
        fabs(div_v_phys) / (fabs(div_v_phys) + curl_v_norm_phys + Balsara_eps);
    balsara = min(1., balsara);
    balsara = max(0., balsara);
  }

  return balsara;
#else
  return 1.;
#endif
}

/**
 * @brief Computes the average kernel gradient between pi and pj
 *
 *
 * @param pi Particle pi
 * @param pj Particle pj
 * @param G_i Kernel gradient at pi
 * @param G_j Kernel gradient at pj
 * @param G_ij Kernel gradient average to fill
 */
__attribute__((always_inline)) INLINE static 
void hydro_get_average_kernel_gradient(const struct part *restrict pi, 
                                       const struct part *restrict pj,
                                       const hydro_real_t *restrict G_i,
                                       const hydro_real_t *restrict G_j,
                                       hydro_real_t *restrict G_ij) {

#if (hydro_props_kernel_gradient_weighting == 0)
  G_ij[0] = 0.5 * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 1)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  G_ij[0] = 0.5 * (f_i * G_i[0] + f_j * G_j[0]);
  G_ij[1] = 0.5 * (f_i * G_i[1] + f_j * G_j[1]);
  G_ij[2] = 0.5 * (f_i * G_i[2] + f_j * G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 2)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t f_ij = 0.5 * (f_i + f_j);
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 3)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  hydro_real_t f_ij = 0.;
  const hydro_real_t f_sum = f_i + f_j;
  if (f_sum > 0.) f_ij = 2. * f_i * f_j / f_sum;

  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 4)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t f_ij = sqrt(f_i * f_j);
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 5)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t rho_sum = pi->rho + pj->rho;
  const hydro_real_t f_ij = (pi->rho * f_i + pj->rho * f_j) / rho_sum;
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 6)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t rho_f_sum = pi->rho * f_i + pj->rho * f_j;
  const hydro_real_t f_ij = 2. * pi->rho * f_i * pj->rho * f_j / rho_f_sum;
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 7)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t P_sum = pi->force.pressure + pj->force.pressure;
  const hydro_real_t f_ij = 
      (pi->force.pressure * f_i + pj->force.pressure * f_j) / P_sum;
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 8)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t P_f_sum = 
      pi->force.pressure * f_i + pj->force.pressure * f_j;
  const hydro_real_t f_ij =
      2. * pi->force.pressure * f_i * pj->force.pressure * f_j / P_f_sum;
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#elif (hydro_props_kernel_gradient_weighting == 9)
  const hydro_real_t f_i = pi->force.f;
  const hydro_real_t f_j = pj->force.f;
  const hydro_real_t mi = hydro_get_mass(pi);
  const hydro_real_t mj = hydro_get_mass(pj);
  const hydro_real_t rhoi_inv = 1. / hydro_get_comoving_density(pi);
  const hydro_real_t rhoj_inv = 1. / hydro_get_comoving_density(pj);
  const hydro_real_t Vi = mi * rhoi_inv;
  const hydro_real_t Vj = mj * rhoj_inv;

  const hydro_real_t f_ij = 0.5 * (Vi * f_i + Vj * f_j) / (Vi + Vj);
  G_ij[0] = 0.5 * f_ij * (G_i[0] + G_j[0]);
  G_ij[1] = 0.5 * f_ij * (G_i[1] + G_j[1]);
  G_ij[2] = 0.5 * f_ij * (G_i[2] + G_j[2]);
#else
  error("Invalid hydro_props_kernel_gradient_weighting value: %d",
        hydro_props_kernel_gradient_weighting);
#endif
}

/**
 * @brief Computes the viscosity acceleration and mu_ij terms.
 *
 *
 * @param pi Particle pi
 * @param pj Particle pj
 * @param dv Second order reconstructed velocity difference
 * @param dx Distance vector
 * @param r Distance vector norm
 * @param fac_mu Cosmological factor for mu_ij
 * @param a2_Hubble Hubble flow
 * @param mu_i The velocity jump indicator at pi to fill
 * @param mu_j The velocity jump indicator at pj to fill
 */ 
__attribute__((always_inline)) INLINE static
hydro_real_t hydro_get_visc_acc_term(const struct part *restrict pi,
                                     const struct part *restrict pj,
                                     const hydro_real_t *restrict dv,
                                     const hydro_real_t *restrict dx,
                                     const hydro_real_t fac_mu,
                                     const hydro_real_t a2_Hubble) {

  const hydro_real_t rhoi = pi->rho;
  const hydro_real_t rhoj = pj->rho;

  const hydro_real_t bi = pi->gradients.balsara;
  const hydro_real_t bj = pj->gradients.balsara;

  /* Viscosity direction */
  hydro_real_t dx_hat[3] = {0., 0., 0.};
  hydro_vec3_unit(dx, dx_hat);

  /* Separation and velocity along the gradient vector */
  const hydro_real_t dv_Hubble[3] = {dv[0] + a2_Hubble * dx[0],
                                     dv[1] + a2_Hubble * dx[1],
                                     dv[2] + a2_Hubble * dx[2]};
  const hydro_real_t dv_dot_dx_hat = hydro_vec3_vec3_dot(dv_Hubble, dx_hat);
  const hydro_real_t conv = (dv_dot_dx_hat < 0.) ? fac_mu : 0.;

  /* Must be a converging flow */
  if (conv == 0.) return 0.;

  const hydro_real_t mu_ij = conv * dv_dot_dx_hat;

#ifdef hydro_props_viscosity_weighting_type
#if (hydro_props_viscosity_weighting_type == 0)

  /* Each particle gets its own Q and then weighted by density */
  const hydro_real_t rhoij_inv = 1.0 / (rhoi * rhoj);

  const hydro_real_t q_i_alpha = 
      -const_viscosity_alpha * pi->force.soundspeed * mu_ij;
  const hydro_real_t q_i_beta = const_viscosity_beta  * mu_ij * mu_ij;
  const hydro_real_t Q_i = rhoi * (q_i_alpha + q_i_beta);

  const hydro_real_t q_j_alpha = 
      -const_viscosity_alpha * pj->force.soundspeed * mu_ij;
  const hydro_real_t q_j_beta = const_viscosity_beta  * mu_ij * mu_ij;
  const hydro_real_t Q_j = rhoj * (q_j_alpha + q_j_beta);

  return (bi * Q_i + bj * Q_j) * rhoij_inv;

#elif (hydro_props_viscosity_weighting_type == 1)

  /* Each particle has the same Q but is density weighted */
  const hydro_real_t b_ij = 0.5 * (bi + bj);
  const hydro_real_t rhoij_inv = 1. / (rhoi * rhoj);
  const hydro_real_t c_ij = 0.5 * (pi->force.soundspeed + pj->force.soundspeed);
  const hydro_real_t q_ij_alpha = -const_viscosity_alpha * c_ij * mu_ij;
  const hydro_real_t q_ij_beta = const_viscosity_beta  * mu_ij * mu_ij;
  const hydro_real_t Q_ij = (rhoi + rhoj) * (q_ij_alpha + q_ij_beta);

  return b_ij * Q_ij * rhoij_inv;

#elif (hydro_props_viscosity_weighting_type == 2)

  /* Particles average symmetrically and arithmetically */
  const hydro_real_t b_ij = 0.5 * (bi + bj);
  const hydro_real_t c_ij = 0.5 * (pi->force.soundspeed + pj->force.soundspeed);
  const hydro_real_t q_ij_alpha = -const_viscosity_alpha * c_ij * mu_ij;
  const hydro_real_t q_ij_beta = const_viscosity_beta  * mu_ij * mu_ij;
  const hydro_real_t q_ij = 2. * (q_ij_alpha + q_ij_beta) / (rhoi + rhoj);

  return b_ij * q_ij;

#endif
#else
  error("Unknown compiled hydro_props_viscosity_weighting_type value: %d\n"
        "Valid values are 0, 1, or 2.\n",
        (int)hydro_props_viscosity_weighting_type);
#endif
}

/**
 * @brief Gets the signal velocity for the conduction interaction based on mu_ij.
 *
 *
 * @param dx The separation vector
 * @param pi Particle pi
 * @param pj Particle pj
 * @param mu_i The velocity jump indicator at pi
 * @param mu_j The velocity jump indicator at pj
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_get_cond_signal_velocity(const float *restrict dx,
                                            const struct part *restrict pi,
                                            const struct part *restrict pj,
                                            const float mu_i,
                                            const float mu_j,
                                            const float beta) {


  const float mu_ij = 0.5f * (mu_i + mu_j);
  return hydro_signal_velocity(dx, pi, pj, mu_ij, beta);
}

/**
 * @brief Gets the signal velocity for the viscous interaction based on mu_ij.
 *
 *
 * @param dx The separation vector
 * @param pi Particle pi
 * @param pj Particle pj
 * @param mu_i The velocity jump indicator at pi
 * @param mu_j The velocity jump indicator at pj
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_get_visc_signal_velocity(const float *restrict dx,
                                            const struct part *restrict pi,
                                            const struct part *restrict pj,
                                            const float mu_i,
                                            const float mu_j,
                                            const float beta) {

#ifdef hydro_props_viscosity_weighting_type
#if (hydro_props_viscosity_weighting_type == 0)
  const float dx_ji[3] = {-dx[0], -dx[1], -dx[2]};
  const float v_sig_visc_i = 
      hydro_signal_velocity(dx, pi, pj, mu_i, beta);
  const float v_sig_visc_j = 
      hydro_signal_velocity(dx_ji, pj, pi, mu_j, beta);
  return 0.5 * (v_sig_visc_i + v_sig_visc_j);
#else
  const float mu_ij = 0.5f * (mu_i + mu_j);
  return hydro_signal_velocity(dx, pi, pj, mu_ij, beta);
#endif
#else
  error("No viscosity weighting type defined at compile time.");
#endif
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
__attribute__((always_inline)) INLINE static void hydro_end_density(
    struct part *restrict p, const struct cosmology *cosmo) {
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  p->rho += p->mass * kernel_root;
  p->density.rho_dh -= hydro_dimension * p->mass * kernel_root;
  p->density.wcount += kernel_root;
  p->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  p->rho *= h_inv_dim;
  p->density.rho_dh *= h_inv_dim_plus_one;
  p->density.wcount *= h_inv_dim;
  p->density.wcount_dh *= h_inv_dim_plus_one;

  /* Need this for correct dh/dt */
  p->gradients.wcount = p->density.wcount;

  p->gradients.u_aux[0] *= h_inv_dim_plus_one;
  p->gradients.u_aux[1] *= h_inv_dim_plus_one;
  p->gradients.u_aux[2] *= h_inv_dim_plus_one;

  p->gradients.u_aux_norm[0] *= h_inv_dim_plus_one;
  p->gradients.u_aux_norm[1] *= h_inv_dim_plus_one;
  p->gradients.u_aux_norm[2] *= h_inv_dim_plus_one;

  if (p->gradients.u_aux_norm[0] != 0. &&
      p->gradients.u_aux_norm[1] != 0. &&
      p->gradients.u_aux_norm[2] != 0.) {
    p->gradients.u_aux[0] /= p->gradients.u_aux_norm[0];
    p->gradients.u_aux[1] /= p->gradients.u_aux_norm[1];
    p->gradients.u_aux[2] /= p->gradients.u_aux_norm[2];

    hydro_vec_slope_limiter(p->gradients.du_min, p->gradients.du_max,
                            p->gradients.kernel_size,
                            p->gradients.u_aux);
  }
  else {
    p->gradients.u_well_conditioned = 0;
#ifdef MAGMA2_DEBUG_CHECKS
    p->debug.u_aux[0] = p->gradients.u_aux[0];
    p->debug.u_aux[1] = p->gradients.u_aux[1];
    p->debug.u_aux[2] = p->gradients.u_aux[2];

    p->debug.u_aux_norm[0] = p->gradients.u_aux_norm[0];
    p->debug.u_aux_norm[1] = p->gradients.u_aux_norm[1];
    p->debug.u_aux_norm[2] = p->gradients.u_aux_norm[2];

    p->debug.u_ill_conditioned_count++;
#endif

    p->gradients.u_aux[0] = 0.;
    p->gradients.u_aux[1] = 0.;
    p->gradients.u_aux[2] = 0.;

    p->gradients.u_aux_norm[0] = 0.;
    p->gradients.u_aux_norm[1] = 0.;
    p->gradients.u_aux_norm[2] = 0.;
  }

  double aux_norm[3][3] = {0};
  for (int i = 0; i < 3; i++) {
    p->gradients.velocity_tensor_aux[i][0] *= h_inv_dim_plus_one;
    p->gradients.velocity_tensor_aux[i][1] *= h_inv_dim_plus_one;
    p->gradients.velocity_tensor_aux[i][2] *= h_inv_dim_plus_one;
    p->gradients.velocity_tensor_aux_norm[i][0] *= h_inv_dim_plus_one;
    p->gradients.velocity_tensor_aux_norm[i][1] *= h_inv_dim_plus_one;
    p->gradients.velocity_tensor_aux_norm[i][2] *= h_inv_dim_plus_one;
    aux_norm[i][0] = p->gradients.velocity_tensor_aux_norm[i][0];
    aux_norm[i][1] = p->gradients.velocity_tensor_aux_norm[i][1];
    aux_norm[i][2] = p->gradients.velocity_tensor_aux_norm[i][2];
  }

  /* Invert the aux_norm matrix */
  gsl_matrix_view A_view = gsl_matrix_view_array((double *)aux_norm, 3, 3);
  gsl_matrix *A = &A_view.matrix;

  double cond = condition_number(A);
  if (cond < const_condition_number_upper_limit) {
    gsl_matrix *A_inv = gsl_matrix_alloc(3, 3);
    gsl_permutation *p_perm = gsl_permutation_alloc(3);
    int signum;

    gsl_linalg_LU_decomp(A, p_perm, &signum);
    gsl_linalg_LU_invert(A, p_perm, A_inv);

    for (int i = 0; i < 3; i++) {
      p->gradients.velocity_tensor_aux_norm[i][0] = gsl_matrix_get(A_inv, i, 0);
      p->gradients.velocity_tensor_aux_norm[i][1] = gsl_matrix_get(A_inv, i, 1);
      p->gradients.velocity_tensor_aux_norm[i][2] = gsl_matrix_get(A_inv, i, 2);
    }

    gsl_matrix_free(A_inv);
    gsl_permutation_free(p_perm);

    /* aux_norm is equation 20 in Rosswog 2020, finalize the gradient in 19 */
    hydro_real_t aux_matrix[3][3] = {0};
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        /**
         * The indices of aux_norm and velocity_gradient_aux are important.
         * aux_norm j: dx vector direction.
         *          k: kernel gradient direction
         *
         * velocity_gradient_aux i: dv vector direction
         *                       k: kernel gradient direction
         *
         * aux_matrix j: spatial derivative direction
         *            i: velocity direction
         */
        for (int k = 0; k < 3; k++) {
          aux_matrix[j][i] += p->gradients.velocity_tensor_aux_norm[j][k] * 
                              p->gradients.velocity_tensor_aux[i][k];
        }
      }
    }

    /* Copy over the matrices for use later */

    /* D. Rennehan: For some reason, memcpy does not work here? Could it 
     * be because of the union in the particle struct? */
    for (int j = 0; j < 3; j++) {
      hydro_vec_slope_limiter(p->gradients.dv_min[j], p->gradients.dv_max[j],
                              p->gradients.kernel_size,
                              aux_matrix[j]);

      p->gradients.velocity_tensor_aux[j][0] = aux_matrix[j][0];
      p->gradients.velocity_tensor_aux[j][1] = aux_matrix[j][1];
      p->gradients.velocity_tensor_aux[j][2] = aux_matrix[j][2];
    }
  }
  else {

    p->gradients.D_well_conditioned = 0;
#ifdef MAGMA2_DEBUG_CHECKS
    p->debug.D_ill_conditioned_count++;
#endif

    /* Ensure no crazy gradients later */
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
#ifdef MAGMA2_DEBUG_CHECKS
        p->debug.velocity_tensor_aux[i][j] =
            p->gradients.velocity_tensor_aux[i][j];
        p->debug.velocity_tensor_aux_norm[i][j] =
            p->gradients.velocity_tensor_aux_norm[i][j];
#endif

        p->gradients.velocity_tensor_aux[i][j] = 0.;
        p->gradients.velocity_tensor_aux_norm[i][j] = 0.;
      }
    }
  }
}

/**
 * @brief Prepare a particle for the gradient calculation.
 *
 * This function is called after the density loop and before the gradient loop.
 *
 * We use it to set the physical timestep for the particle and to copy the
 * actual velocities, which we need to boost our interfaces during the flux
 * calculation. We also initialize the variables used for the time step
 * calculation.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 * @param hydro_props Hydrodynamic properties.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_gradient(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {

#ifdef hydro_props_use_adiabatic_correction
  /* Prepare the denominator for the adiabatic correction term */
  p->gradients.adiabatic_f_denominator = 0.;
#endif

  /* Compute the sound speed  */
  const float pressure = hydro_get_comoving_pressure(p);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  
  /* ------ Compute the kernel correction for SPH gradients ------ */


  /* Compute the "grad h" term */
  const float common_factor = p->h * hydro_dimension_inv / p->density.wcount;
  float grad_W_term = -1.f;

  /* Ignore changing-kernel effects when h ~= h_max */
  if (p->h > 0.9999f * hydro_props->h_max) {
    warning("h ~ h_max for particle with ID %lld (h: %g)", p->id, p->h);
  } 
  else {
    grad_W_term = common_factor * p->density.wcount_dh;

    if (grad_W_term < -0.9999f) {
      /* if we get here, we either had very small neighbour contributions
         (which should be treated as a no neighbour case in the ghost) or
         a very weird particle distribution (e.g. particles sitting on
         top of each other). Either way, we cannot use the normal
         expression, since that would lead to overflow or excessive round
         off and cause excessively high accelerations in the force loop */
      grad_W_term = -1.f;
      warning(
          "grad_W_term very small for particle with ID %lld (h: %g, wcount: "
          "%g, wcount_dh: %g)",
          p->id, p->h, p->density.wcount, p->density.wcount_dh);
    } 
  }

  /* Update variables. */
  p->force.f = (grad_W_term > -1.f) ? 1.f / (1.f + grad_W_term) : 1.f;

}

/**
 * @brief Resets the variables that are required for a gradient calculation.
 *
 * This function is called after hydro_prepare_gradient.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_gradient(
    struct part *restrict p) {

  p->gradients.C_well_conditioned = 1;

  for (int i = 0; i < 3; i++) {
    p->gradients.u[i] = 0.;
    p->gradients.u_hessian[i][0] = 0.;
    p->gradients.u_hessian[i][1] = 0.;
    p->gradients.u_hessian[i][2] = 0.;

    p->gradients.correction_matrix[i][0] = 0.;
    p->gradients.correction_matrix[i][1] = 0.;
    p->gradients.correction_matrix[i][2] = 0.;

    p->gradients.velocity_tensor[i][0] = 0.;
    p->gradients.velocity_tensor[i][1] = 0.;
    p->gradients.velocity_tensor[i][2] = 0.;

    for (int j = 0; j < 3; j++) {
      p->gradients.velocity_hessian[i][j][0] = 0.;
      p->gradients.velocity_hessian[i][j][1] = 0.;
      p->gradients.velocity_hessian[i][j][2] = 0.; 
    }
  }
}

/**
 * @brief Finishes the gradient calculation.
 *
 * Just a wrapper around hydro_gradients_finalize, which can be an empty method,
 * in which case no gradients are used.
 *
 * This method also initializes the force loop variables.
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_end_gradient(
    struct part *p) {

  /* Calculate smoothing length powers */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  p->gradients.u[0] *= h_inv_dim;
  p->gradients.u[1] *= h_inv_dim;
  p->gradients.u[2] *= h_inv_dim;

  /* Temporary double for GSL */
  double correction_matrix[3][3] = {0};

  /* Apply correct normalisation */
  for (int k = 0; k < 3; k++) {
    p->gradients.correction_matrix[k][0] *= h_inv_dim;
    p->gradients.correction_matrix[k][1] *= h_inv_dim;
    p->gradients.correction_matrix[k][2] *= h_inv_dim;

    correction_matrix[k][0] = p->gradients.correction_matrix[k][0];
    correction_matrix[k][1] = p->gradients.correction_matrix[k][1];
    correction_matrix[k][2] = p->gradients.correction_matrix[k][2];

    p->gradients.u_hessian[k][0] *= h_inv_dim;
    p->gradients.u_hessian[k][1] *= h_inv_dim;
    p->gradients.u_hessian[k][2] *= h_inv_dim;

    p->gradients.velocity_tensor[k][0] *= h_inv_dim;
    p->gradients.velocity_tensor[k][1] *= h_inv_dim;
    p->gradients.velocity_tensor[k][2] *= h_inv_dim;

    for (int j = 0; j < 3; j++) {
      p->gradients.velocity_hessian[k][j][0] *= h_inv_dim;
      p->gradients.velocity_hessian[k][j][1] *= h_inv_dim;
      p->gradients.velocity_hessian[k][j][2] *= h_inv_dim;
    }
  }

  /* Invert the p->gradients.correction_matrix[3][3] matrix */
  gsl_matrix_view C_view = 
      gsl_matrix_view_array((double *)correction_matrix, 3, 3);
  gsl_matrix *C = &C_view.matrix;

  const double cond = condition_number(C);
  if (cond < const_condition_number_upper_limit) {
    gsl_matrix *C_inv = gsl_matrix_alloc(3, 3);
    gsl_permutation *p_perm = gsl_permutation_alloc(3);
    int signum;

    gsl_linalg_LU_decomp(C, p_perm, &signum);
    gsl_linalg_LU_invert(C, p_perm, C_inv);

    for (int k = 0; k < 3; k++) {
      p->gradients.correction_matrix[k][0] = gsl_matrix_get(C_inv, k, 0);
      p->gradients.correction_matrix[k][1] = gsl_matrix_get(C_inv, k, 1);
      p->gradients.correction_matrix[k][2] = gsl_matrix_get(C_inv, k, 2);
    }

    gsl_matrix_free(C_inv);
    gsl_permutation_free(p_perm);
  }
  else {
#ifdef MAGMA2_DEBUG_CHECKS
    for (int k = 0; k < 3; k++) {
      p->debug.correction_matrix[k][0] = p->gradients.correction_matrix[k][0];
      p->debug.correction_matrix[k][1] = p->gradients.correction_matrix[k][1];
      p->debug.correction_matrix[k][2] = p->gradients.correction_matrix[k][2];
    }
    
    p->debug.C_ill_conditioned_count++;
#endif

    /* Ill-condition matrix, revert back to normal SPH gradients */
    p->gradients.C_well_conditioned = 0;
    for (int k = 0; k < 3; k++) {
      p->gradients.correction_matrix[k][0] = 0.;
      p->gradients.correction_matrix[k][1] = 0.;
      p->gradients.correction_matrix[k][2] = 0.;
    }
  }

  /* Contract the correction matrix with the internal energy gradient and with 
   * the velocity tensor */
  hydro_real_t u_gradient[3] = {0};
  hydro_real_t u_hessian[3][3] = {0};
  hydro_real_t velocity_tensor[3][3] = {0};
  hydro_real_t velocity_hessian[3][3][3] = {0};
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      u_gradient[j] += p->gradients.correction_matrix[j][i] *
                       p->gradients.u[i];

      for (int k = 0; k < 3; k++) {
        u_hessian[j][i] += p->gradients.correction_matrix[j][k] *
                           p->gradients.u_hessian[i][k];
        velocity_tensor[j][i] += p->gradients.correction_matrix[j][k] *
                                 p->gradients.velocity_tensor[i][k];

        for (int m = 0; m < 3; m++) {
          velocity_hessian[j][i][k] += p->gradients.correction_matrix[j][m] * 
                                       p->gradients.velocity_hessian[j][i][m];
        }
      }
    }
  }

  /* Slope limiter for internal energy */
  hydro_vec_slope_limiter(p->gradients.du_min, p->gradients.du_max,
                          p->gradients.kernel_size,
                          u_gradient);

  /* Copy back over to the particle for later */
  p->gradients.u[0] = u_gradient[0];
  p->gradients.u[1] = u_gradient[1];
  p->gradients.u[2] = u_gradient[2];

  for (int j = 0; j < 3; j++) {
    p->gradients.u_hessian[j][0] = u_hessian[j][0];
    p->gradients.u_hessian[j][1] = u_hessian[j][1];
    p->gradients.u_hessian[j][2] = u_hessian[j][2];

    hydro_vec_slope_limiter(p->gradients.dv_min[j], p->gradients.dv_max[j],
                            p->gradients.kernel_size,
                            velocity_tensor[j]);

    p->gradients.velocity_tensor[j][0] = velocity_tensor[j][0];
    p->gradients.velocity_tensor[j][1] = velocity_tensor[j][1];
    p->gradients.velocity_tensor[j][2] = velocity_tensor[j][2];

    p->gradients.velocity_hessian[j][0][0] = velocity_hessian[j][0][0];
    p->gradients.velocity_hessian[j][0][1] = velocity_hessian[j][0][1];
    p->gradients.velocity_hessian[j][0][2] = velocity_hessian[j][0][2];

    p->gradients.velocity_hessian[j][1][0] = velocity_hessian[j][1][0];
    p->gradients.velocity_hessian[j][1][1] = velocity_hessian[j][1][1];
    p->gradients.velocity_hessian[j][1][2] = velocity_hessian[j][1][2];

    p->gradients.velocity_hessian[j][2][0] = velocity_hessian[j][2][0];
    p->gradients.velocity_hessian[j][2][1] = velocity_hessian[j][2][1];
    p->gradients.velocity_hessian[j][2][2] = velocity_hessian[j][2][2];
  }
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
__attribute__((always_inline)) INLINE static void hydro_part_has_no_neighbours(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo) {
  /* Some smoothing length multiples. */
  const float h = p->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  warning(
      "Gas particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      p->id, h, p->density.wcount);

  /* Re-set problematic values */
  p->rho = p->mass * kernel_root * h_inv_dim;
  p->h_min = FLT_MAX;
  p->dt_min = FLT_MAX;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.f;
  p->density.wcount_dh = 0.f;

#ifdef MAGMA2_DEBUG_CHECKS
  p->debug.num_ngb = 0;
  p->debug.v_sig_visc_max = 0;
  p->debug.v_sig_cond_max = 0;
#endif
  p->gradients.C_well_conditioned = 0;
  p->gradients.D_well_conditioned = 0;
  p->gradients.u_well_conditioned = 0;
  p->gradients.du_min = 0.;
  p->gradients.du_max = 0.;
  p->gradients.kernel_size = (hydro_real_t)h;

  for (int i = 0; i < 3; i++) {
    p->gradients.u[i] = 0.;
    p->gradients.u_aux[i] = 0.;
    p->gradients.u_aux_norm[i] = 0.;

    p->gradients.dv_min[i] = 0.;
    p->gradients.dv_max[i] = 0.;

    for (int j = 0; j < 3; j++) {
      p->gradients.correction_matrix[i][j] = 0.;
      p->gradients.velocity_tensor[i][j] = 0.;
      p->gradients.velocity_tensor_aux[i][j] = 0.;
      p->gradients.velocity_tensor_aux_norm[i][j] = 0.;
      p->gradients.u_hessian[i][j] = 0.;

      for (int k = 0; k < 3; k++) {
        p->gradients.velocity_hessian[i][j][k] = 0.;
      }
    }
  }
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
 * @param pressure_floor The properties of the pressure floor.
 * @param dt_alpha The time-step used to evolve non-cosmological quantities such
 *                 as the artificial viscosity.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_prepare_force(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor, const float dt_alpha,
    const float dt_therm) {

  /* First estimates for the timestepping. Missing the kernel_gamma factors
   * for now, but will be added at the end of the force loop. */
  p->h_min = FLT_MAX;
  p->dt_min = FLT_MAX;
  p->gradients.balsara = hydro_get_balsara_limiter(p, cosmo);

#ifdef hydro_props_use_adiabatic_correction
  const hydro_real_t prev_f = p->force.f;
  const hydro_real_t rho_inv = 1. / hydro_get_comoving_density(p);
  /* Finish final kernel gradient correction factor */
  if (p->gradients.adiabatic_f_denominator != 0.) {
    const hydro_real_t kernel_r2 = p->gradients.adiabatic_f_numerator;
    const hydro_real_t weighted_kernel_r2_inv = 
        1. / p->gradients.adiabatic_f_denominator;
    p->force.f = rho_inv * kernel_r2 * weighted_kernel_r2_inv;
  }
  else {
    p->force.f = 1.;
  }

#ifdef MAGMA2_DEBUG_CHECKS
  if (p->force.f > 100.) {
    warning("Final force factor is very high for particle with ID %lld"
            " (prev f: %g, f: %g, numerator: %g, denominator: %g, rho_inv: %g)",
            p->id, prev_f, p->force.f, p->gradients.adiabatic_f_numerator,
            p->gradients.adiabatic_f_denominator, rho_inv);
  }
#endif
#endif
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * Resets all hydro acceleration and time derivative fields in preparation
 * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_reset_acceleration(
    struct part *restrict p) {
  /* Reset the acceleration. */
  p->a_hydro[0] = 0.0f;
  p->a_hydro[1] = 0.0f;
  p->a_hydro[2] = 0.0f;

  /* Reset the time derivatives. */
  p->u_dt = 0.0f;
  p->u_dt_cond = 0.0f;
  p->force.h_dt = 0.0f;
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 * @param cosmo The cosmological model.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp,
    const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor) {
  /* Re-set the predicted velocities */
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];

  /* Re-set the entropy */
  p->u = xp->u_full;

  /* Compute the sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Signal speed */
  const float v_sig = const_viscosity_alpha_prefactor * soundspeed;

  /* Update the signal velocity, if we need to. */
  p->dt_min = min(p->dt_min, p->h_min / v_sig);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * Additional hydrodynamic quantites are drifted forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * Note the different time-step sizes used for the different quantities as they
 * include cosmological factors.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @param dt_drift The drift time-step for positions.
 * @param dt_therm The drift time-step for thermal quantities.
 * @param dt_kick_grav The time-step for gravity quantities.
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor_props The properties of the entropy floor.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra(
    struct part *restrict p, const struct xpart *restrict xp, float dt_drift,
    float dt_therm, float dt_kick_grav, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct pressure_floor_props *pressure_floor) {

  /* Predict the internal energy */
  p->u += p->u_dt * dt_therm;

  const float h_inv = 1.f / p->h;

  /* Predict smoothing length */
  const float w1 = p->force.h_dt * h_inv * dt_drift;
  if (fabsf(w1) < 0.2f)
    p->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    p->h *= expf(w1);

  /* Predict density and weighted pressure */
  const float w2 = -hydro_dimension * w1;
  if (fabsf(w2) < 0.2f) {
    const float expf_approx =
        approx_expf(w2); /* 4th order expansion of exp(w) */
    p->rho *= expf_approx;
  } else {
    const float expf_exact = expf(w2);
    p->rho *= expf_exact;
  }

  /* Check against entropy floor - explicitly do this after drifting the
   * density as this has a density dependence. */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  p->u = max(p->u, floor_u);
  p->u = max(p->u, min_u);

  /* Compute the new sound speed */
  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;

  /* Signal speed */
  const float v_sig = const_viscosity_alpha_prefactor * soundspeed;

  /* Update signal velocity if we need to */
  p->dt_min = min(p->dt_min, p->h_min / v_sig);
}

/**
 * @brief Returns the sum term for the dh/dt calculation
 *
 *
 * @param m The particle mass of the neighbour
 * @param rho_inv The inverse density of the neighbour
 * @param r_inv The inverse distance between particles
 * @param w_dr The kernel gradient for this particle
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_get_h_dt_sum(const hydro_real_t dv_dot_dx,
                                const hydro_real_t dv_dot_G,
                                const hydro_real_t m, 
                                const hydro_real_t rho_inv,
                                const hydro_real_t r_inv,
                                const hydro_real_t w_dr) {

  hydro_real_t dvdx = 0.;
#ifdef hydro_props_dh_dt_estimator_type
#if (hydro_props_dh_dt_estimator_type == 0)
  const hydro_real_t grad = r_inv * w_dr;
  const hydro_real_t wt = m * rho_inv;
  dvdx = dv_dot_dx;
#elif (hydro_props_dh_dt_estimator_type == 1)
  const hydro_real_t grad = r_inv * w_dr;
  const hydro_real_t wt = 1.;
  dvdx = dv_dot_dx;
#elif (hydro_props_dh_dt_estimator_type == 2)
  const hydro_real_t grad = 1.;
  const hydro_real_t wt = m * rho_inv;
  dvdx = dv_dot_G;
#else
  error("Compiled with an unknown dh/dt estimator type");
#endif
#else
  error("Must compile with hydro_props_dh_dt_estimator_type.");
#endif

  return wt * dvdx * grad;
}

/**
 * @brief Returns the normalization for the dh/dt sum
 *
 *
 * @param p The particle
 */
__attribute__((always_inline)) INLINE static 
hydro_real_t hydro_get_h_dt_norm(struct part *restrict p) {

  hydro_real_t renormalization = 1.;

#ifdef hydro_props_dh_dt_estimator_type
#if (hydro_props_dh_dt_estimator_type == 1)
  renormalization = p->force.f / p->gradients.wcount;
#endif
#else
  error("Must compile with hydro_props_dh_dt_estimator_type.");
#endif

  return renormalization * p->h * hydro_dimension_inv;
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
__attribute__((always_inline)) INLINE static void hydro_end_force(
    struct part *restrict p, const struct cosmology *cosmo) {

  p->force.h_dt *= hydro_get_h_dt_norm(p);
}

/**
 * @brief Kick the additional variables
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 * @param dt_grav The time-step for this kick (for gravity quantities).
 * @param dt_grav_mesh The time-step for this kick (mesh gravity).
 * @param dt_hydro The time-step for this kick (for hydro quantities).
 * @param dt_kick_corr The time-step for this kick (for gravity corrections).
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme
 * @param floor_props The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void hydro_kick_extra(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm,
    float dt_grav, float dt_grav_mesh, float dt_hydro, float dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props) {
  /* Integrate the internal energy forward in time */
  const float delta_u = p->u_dt * dt_therm;

  /* Do not decrease the energy by more than a factor of 2*/
  xp->u_full = max(xp->u_full + delta_u, 0.5f * xp->u_full);

  /* Check against entropy floor */
  const float floor_A = entropy_floor(p, cosmo, floor_props);
  const float floor_u = gas_internal_energy_from_entropy(p->rho, floor_A);

  /* Check against absolute minimum */
  const float min_u =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;

  /* Take highest of both limits */
  const float energy_min = max(min_u, floor_u);

  if (xp->u_full < energy_min) {
    xp->u_full = energy_min;
    p->u_dt = 0.f;
  }
}

/**
 * @brief Converts hydro quantity of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 * This can be used to convert internal energy into entropy for instance.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param pressure_floor The properties of the pressure floor.
 */
__attribute__((always_inline)) INLINE static void hydro_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct pressure_floor_props *pressure_floor) {
  /* Convert the physcial internal energy to the comoving one. */
  /* u' = a^(3(g-1)) u */
  const float factor = 1.f / cosmo->a_factor_internal_energy;
  p->u *= factor;
  xp->u_full = p->u;

  /* Apply the minimal energy limit */
  const float min_comoving_energy =
      hydro_props->minimal_internal_energy / cosmo->a_factor_internal_energy;
  if (xp->u_full < min_comoving_energy) {
    xp->u_full = min_comoving_energy;
    p->u = min_comoving_energy;
    p->u_dt = 0.f;
  }

  /* Set the initial value of the artificial viscosity based on the non-variable
     schemes for safety */

  const float pressure = gas_pressure_from_internal_energy(p->rho, p->u);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;
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
__attribute__((always_inline)) INLINE static void hydro_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {
  p->time_bin = 0;

  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
  xp->u_full = p->u;

  p->gradients.C_well_conditioned = 1;
  p->gradients.D_well_conditioned = 1;
  p->gradients.u_well_conditioned = 1;
#ifdef MAGMA2_DEBUG_CHECKS
  p->debug.N_force_low_order_grad = 0;
  p->debug.N_force_high_order_grad = 0;
  p->debug.C_ill_conditioned_count = 0;
  p->debug.D_ill_conditioned_count = 0;
  p->debug.u_ill_conditioned_count = 0;
  for (int i = 0; i < 3; i++) {
    p->debug.correction_matrix[i][0] = 0.;
    p->debug.correction_matrix[i][1] = 0.;
    p->debug.correction_matrix[i][2] = 0.;

    p->debug.velocity_tensor_aux[i][0] = 0.;
    p->debug.velocity_tensor_aux[i][1] = 0.;
    p->debug.velocity_tensor_aux[i][2] = 0.;

    p->debug.velocity_tensor_aux_norm[i][0] = 0.;
    p->debug.velocity_tensor_aux_norm[i][1] = 0.;
    p->debug.velocity_tensor_aux_norm[i][2] = 0.;
  }
#endif

  hydro_reset_acceleration(p);
  hydro_init_part(p, NULL);
}

/**
 * @brief Overwrite the initial internal energy of a particle.
 *
 * Note that in the cases where the thermodynamic variable is not
 * internal energy but gets converted later, we must overwrite that
 * field. The conversion to the actual variable happens later after
 * the initial fake time-step.
 *
 * @param p The #part to write to.
 * @param u_init The new initial internal energy.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_init_internal_energy(struct part *p, float u_init) {
  p->u = u_init;
}

#endif /* SWIFT_MAGMA2_HYDRO_H */
