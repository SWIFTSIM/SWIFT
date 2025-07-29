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
 * @brief Density-Energy conservative implementation of SPH,
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

  p->v_sig_max = max(p->v_sig_max, soundspeed);
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

  /* Update the signal velocity */
  p->v_sig_max =
      max(soundspeed, p->v_sig_max + const_viscosity_beta * dv);
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

  /* CFL condition */
  return hydro_properties->CFL_condition * p->dt_min;
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

  return 0.5f * (ci + cj - beta * mu_ij);
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
  p->num_ngb = 0;

  /* These must be zeroed before the density loop */
  for (int i = 0; i < 3; i++) {
    p->gradients.u_aux[i] = 0.f;
    p->gradients.u_aux_norm[i] = 0.f;

    for (int j = 0; j < 3; j++) {
      p->gradients.velocity_tensor_aux[i][j] = 0.f;
      p->gradients.velocity_tensor_aux_norm[i][j] = 0.f;
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
void hydro_vec3_vec3_dot(const float *restrict vec_a, 
                         const float *restrict vec_b, 
                         float *result) {

  *result = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2];
}

/**
 * @brief The scalar triple product of vector vec and matrix mat.
 *        That is vec^T * mat * vec.
 *
 *
 * @param mat The matrix to contract with the vector.
 * @param vec The vector to contract with the matrix.
 * @param result The result of the contraction.
 */
__attribute__((always_inline)) INLINE static 
void hydro_mat3x3_mat3x3_dot(const float (*restrict mat_a)[3], 
                             const float (*restrict mat_b)[3], 
                             float *result) {

  *result = mat_a[0][0] * mat_b[0][0] +
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
    const float (*restrict mat)[3], 
    const float *restrict vec, 
    float *restrict result) {

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
    const float (*restrict tensor)[3][3],
    const float (*restrict mat)[3],
    float *restrict result) {

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

  result[3] = tensor[2][0][0] * mat[0][0] +
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
    const float *restrict vec_a, 
    const float *restrict vec_b, 
    float (*restrict result)[3]) {

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
 * @brief Computes Phi_ab from Rosswog 2020 21 for a field given A. This is the 
 * van Leer 1974 slope limiting procedure.
 *
 *
 * @param A_ij The ratio of the gradients of the two particles.
 * @param eta_i The normed smoothing length of the first particle.
 * @param eta_j The normed smoothing length of the second particle.
 * @param num_ngb The number of neighbours in the scheme.
 */
__attribute__((always_inline)) INLINE static
float hydro_van_leer_phi(const float A_ij,
                         const float eta_i, 
                         const float eta_j) {

  float phi_raw = (4.f * A_ij) / ((1.f + A_ij) * (1.f + A_ij));
  phi_raw = fminf(1.f, fmaxf(0.f, phi_raw));

  /* η_ab and η_crit damping */
  const float eta_ij = fminf(eta_i, eta_j);

  float damping = 1.f;
  if (eta_ij <= const_slope_limiter_eta_crit) {
    const float diff = 
        (eta_ij - const_slope_limiter_eta_crit) / const_slope_limiter_eta_fold;
    damping = expf(-diff * diff);
  }

  return phi_raw * damping;
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
float hydro_scalar_van_leer_A(const float *restrict grad_i,
                              const float *restrict grad_j,
                              const float *restrict dx) {

  float grad_dot_x_i = 0.f;
  float grad_dot_x_j = 0.f;
  hydro_vec3_vec3_dot(grad_i, dx, &grad_dot_x_i);
  hydro_vec3_vec3_dot(grad_j, dx, &grad_dot_x_j);

  /* Regularize denominator */
  if (fabsf(grad_dot_x_j) < FLT_EPSILON) return 0.f;

  const float A_ij = grad_dot_x_i / grad_dot_x_j;

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
float hydro_vector_van_leer_A(const float (*restrict grad_i)[3],
                              const float (*restrict grad_j)[3],
                              const float *restrict dx) {

  float delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(dx, dx, delta_ij);

  float grad_dot_x_x_i = 0.f;
  float grad_dot_x_x_j = 0.f;
  hydro_mat3x3_mat3x3_dot(grad_i, delta_ij, &grad_dot_x_x_i);
  hydro_mat3x3_mat3x3_dot(grad_j, delta_ij, &grad_dot_x_x_j);

  /* Regularize denominator */
  if (fabsf(grad_dot_x_x_j) < FLT_EPSILON) return 0.f;

  const float A_ij = grad_dot_x_x_i / grad_dot_x_x_j;

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
 * @param num_ngb The number of neighbours in the scheme.
 */
__attribute__((always_inline)) INLINE static
float hydro_scalar_van_leer_phi(const float *restrict grad_i,
                                const float *restrict grad_j,
                                const float *restrict dx,
                                const float eta_i, 
                                const float eta_j) {

  const float A_ij = hydro_scalar_van_leer_A(grad_i, grad_j, dx);

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
 * @param num_ngb The number of neighbours in the scheme.
 */
__attribute__((always_inline)) INLINE static
float hydro_vector_van_leer_phi(const float (*restrict grad_i)[3],
                                const float (*restrict grad_j)[3],
                                const float *restrict dx,
                                const float eta_i, 
                                const float eta_j) {

  const float A_ij = hydro_vector_van_leer_A(grad_i, grad_j, dx);

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
__attribute__((always_inline)) INLINE static
void hydro_scalar_second_order_reconstruction(const float phi,
                                              const float *restrict dx,
                                              const float f, 
                                              const float *restrict grad,
                                              const float (*restrict hess)[3],
                                              float *f_reconstructed) {

  /* Midpoint from Equation 17 Rosswog 2020 */
  const float midpoint[3] = {0.5f * dx[0],
                             0.5f * dx[1],
                             0.5f * dx[2]};
  float delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(midpoint, midpoint, delta_ij);

  float df_dx_dot_r = 0.f;
  float df_dx_dx_dot_r2 = 0.f;

  hydro_vec3_vec3_dot(grad, midpoint, &df_dx_dot_r);
  hydro_mat3x3_mat3x3_dot(hess, delta_ij, &df_dx_dx_dot_r2);

  /* Apply limited slope reconstruction */
  *f_reconstructed = f + phi * (df_dx_dot_r + 0.5f * df_dx_dx_dot_r2);
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
__attribute__((always_inline)) INLINE static
void hydro_vector_second_order_reconstruction(const float phi,
                                              const float *restrict dx,
                                              const float *restrict f, 
                                              const float (*restrict grad)[3],
                                              const float (*restrict hess)[3][3],
                                              float *restrict f_reconstructed) {

  /* Midpoint from Equation 17 Rosswog 2020 */
  const float midpoint[3] = {0.5f * dx[0],
                             0.5f * dx[1],
                             0.5f * dx[2]};
  float delta_ij[3][3] = {0};
  hydro_vec3_vec3_outer(midpoint, midpoint, delta_ij);

  float df_dx_dot_r[3] = {0};
  float df_dx_dx_dot_r2[3] = {0};

  hydro_mat3x3_vec3_dot(grad, midpoint, df_dx_dot_r);
  hydro_tensor3x3x3_matrix3x3_dot(hess, delta_ij, df_dx_dx_dot_r2);

  /* Apply limited slope reconstruction */
  f_reconstructed[0] = f[0] + phi * (df_dx_dot_r[0] + 0.5f * df_dx_dx_dot_r2[0]);
  f_reconstructed[1] = f[1] + phi * (df_dx_dot_r[1] + 0.5f * df_dx_dx_dot_r2[1]);
  f_reconstructed[2] = f[2] + phi * (df_dx_dot_r[2] + 0.5f * df_dx_dx_dot_r2[2]);
}

/**
 e @brief Finishes the density calculation.
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

  p->gradients.u_aux[0] *= h_inv_dim_plus_one;
  p->gradients.u_aux[1] *= h_inv_dim_plus_one;
  p->gradients.u_aux[2] *= h_inv_dim_plus_one;

  p->gradients.u_aux_norm[0] *= h_inv_dim_plus_one;
  p->gradients.u_aux_norm[1] *= h_inv_dim_plus_one;
  p->gradients.u_aux_norm[2] *= h_inv_dim_plus_one;

  if (fabsf(p->gradients.u_aux_norm[0]) > FLT_EPSILON &&
      fabsf(p->gradients.u_aux_norm[1]) > FLT_EPSILON &&
      fabsf(p->gradients.u_aux_norm[2]) > FLT_EPSILON) {
    p->gradients.u_aux[0] /= p->gradients.u_aux_norm[0];
    p->gradients.u_aux[1] /= p->gradients.u_aux_norm[1];
    p->gradients.u_aux[2] /= p->gradients.u_aux_norm[2];
  }
  else {
    p->gradients.u_aux[0] = 0.f;
    p->gradients.u_aux[1] = 0.f;
    p->gradients.u_aux[2] = 0.f;

    p->gradients.u_aux_norm[0] = 0.f;
    p->gradients.u_aux_norm[1] = 0.f;
    p->gradients.u_aux_norm[2] = 0.f;
  }

  double aux_norm[3][3];
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
    float aux_matrix[3][3] = {0};
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
      p->gradients.velocity_tensor_aux[j][0] = aux_matrix[j][0];
      p->gradients.velocity_tensor_aux[j][1] = aux_matrix[j][1];
      p->gradients.velocity_tensor_aux[j][2] = aux_matrix[j][2];
    }
  }
  else {

#ifdef MAGMA2_DEBUG_CHECKS
    p->debug.D_well_conditioned = 0;
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

        p->gradients.velocity_tensor_aux[i][j] = 0.f;
        p->gradients.velocity_tensor_aux_norm[i][j] = 0.f;
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

  /* Compute the sound speed  */
  const float pressure = hydro_get_comoving_pressure(p);
  const float pressure_including_floor =
      pressure_floor_get_comoving_pressure(p, pressure_floor, pressure, cosmo);
  const float soundspeed =
      gas_soundspeed_from_pressure(p->rho, pressure_including_floor);

  /* Update variables. */
  p->force.pressure = pressure_including_floor;
  p->force.soundspeed = soundspeed;
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
#ifdef MAGMA2_DEBUG_CHECKS
  p->debug.D_well_conditioned = 1;
#endif

  for (int i = 0; i < 3; i++) {
    p->gradients.u[i] = 0.f;
    p->gradients.u_hessian[i][0] = 0.f;
    p->gradients.u_hessian[i][1] = 0.f;
    p->gradients.u_hessian[i][2] = 0.f;

    p->gradients.correction_matrix[i][0] = 0.f;
    p->gradients.correction_matrix[i][1] = 0.f;
    p->gradients.correction_matrix[i][2] = 0.f;

    p->gradients.velocity_tensor[i][0] = 0.f;
    p->gradients.velocity_tensor[i][1] = 0.f;
    p->gradients.velocity_tensor[i][2] = 0.f;

    for (int j = 0; j < 3; j++) {
      p->gradients.velocity_hessian[i][j][0] = 0.f;
      p->gradients.velocity_hessian[i][j][1] = 0.f;
      p->gradients.velocity_hessian[i][j][2] = 0.f; 
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

  /* Normalize correctly with the smoothing length */
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

  double cond = condition_number(C);
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
      p->gradients.correction_matrix[k][0] = 0.f;
      p->gradients.correction_matrix[k][1] = 0.f;
      p->gradients.correction_matrix[k][2] = 0.f;
    }
  }

  /* Contract the correction matrix with the internal energy gradient and with 
   * the velocity tensor */
  float u_gradient[3] = {0};
  float u_hessian[3][3] = {0};
  float velocity_tensor[3][3] = {0};
  float velocity_hessian[3][3][3] = {0};
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

  /* Copy back over to the particle for later */
  p->gradients.u[0] = u_gradient[0];
  p->gradients.u[1] = u_gradient[1];
  p->gradients.u[2] = u_gradient[2];

  for (int j = 0; j < 3; j++) {
    p->gradients.u_hessian[j][0] = u_hessian[j][0];
    p->gradients.u_hessian[j][1] = u_hessian[j][1];
    p->gradients.u_hessian[j][2] = u_hessian[j][2];

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
  p->h_min = 0.f;
  p->dt_min = 0.f;
  p->v_sig_max = 0.f;
  p->density.wcount = kernel_root * h_inv_dim;
  p->density.rho_dh = 0.f;
  p->density.wcount_dh = 0.f;

  p->num_ngb = 0;
  p->gradients.C_well_conditioned = 0;
  for (int i = 0; i < 3; i++) {
    p->gradients.u[i] = 0.f;
    p->gradients.u_aux[i] = 0.f;
    p->gradients.u_aux_norm[i] = 0.f;

    for (int j = 0; j < 3; j++) {
      p->gradients.correction_matrix[i][j] = 0.f;
      p->gradients.velocity_tensor[i][j] = 0.f;
      p->gradients.velocity_tensor_aux[i][j] = 0.f;
      p->gradients.velocity_tensor_aux_norm[i][j] = 0.f;
      p->gradients.u_hessian[i][j] = 0.f;

      for (int k = 0; k < 3; k++) {
        p->gradients.velocity_hessian[i][j][k] = 0.f;
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
  p->h_min = p->h;
  p->v_sig_max = p->force.soundspeed;
  p->dt_min = p->h_min / p->v_sig_max;
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

  /* Update the signal velocity, if we need to. */
  p->v_sig_max = max(p->v_sig_max, soundspeed);
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

  /* Update signal velocity if we need to */
  p->v_sig_max = max(p->v_sig_max, soundspeed);
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
  
  const float rho_inv = hydro_get_comoving_density(p);  
  p->force.h_dt *= p->h * rho_inv * hydro_dimension_inv;

  /* dt_min is in physical units, and requires the kernel_gamma factor for h */
  p->dt_min *= kernel_gamma * cosmo->a / cosmo->a_factor_sound_speed;

#ifdef MAGMA2_DEBUG_CHECKS
  if (p->force.h_dt > 1.e10f) {
    warning(
        "Particle %lld has a very large h_dt value (%g). This may be due to "
        "a very low density or a very high smoothing length."
        "rho_inv = %g, h = %g, hydro_dimension_inv = %g",
        p->id, p->force.h_dt,
        rho_inv, p->h, hydro_dimension_inv);
  }
#endif
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
#ifdef MAGMA2_DEBUG_CHECKS
  p->debug.C_ill_conditioned_count = 0;
  p->debug.D_well_conditioned = 1;
  p->debug.D_ill_conditioned_count = 0;
  for (int i = 0; i < 3; i++) {
    p->debug.correction_matrix[i][0] = 0.f;
    p->debug.correction_matrix[i][1] = 0.f;
    p->debug.correction_matrix[i][2] = 0.f;

    p->debug.velocity_tensor_aux[i][0] = 0.f;
    p->debug.velocity_tensor_aux[i][1] = 0.f;
    p->debug.velocity_tensor_aux[i][2] = 0.f;

    p->debug.velocity_tensor_aux_norm[i][0] = 0.f;
    p->debug.velocity_tensor_aux_norm[i][1] = 0.f;
    p->debug.velocity_tensor_aux_norm[i][2] = 0.f;
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
