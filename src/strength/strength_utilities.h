/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_STRENGTH_UTILITIES_H
#define SWIFT_STRENGTH_UTILITIES_H

/**
 * @file strength/strength_utilities.h
 * @brief Utilities used throughout the material strength scheme.
 */

#include "math.h"
#include "symmetric_matrix.h"

/**
 * @brief Computes the J_2 invariant of the deviatoric stress tensor.
 *
 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 */
__attribute__((always_inline)) INLINE static float strength_compute_stress_tensor_J_2(
    const struct sym_matrix deviatoric_stress_tensor) {

  // ### Does j_2 need to be decreased by a factor of (1 - damage)^2 for B&A?

  /* Compute J_2 invariant. */
  return 0.5f * deviatoric_stress_tensor.xx * deviatoric_stress_tensor.xx +
         0.5f * deviatoric_stress_tensor.yy * deviatoric_stress_tensor.yy +
         0.5f * deviatoric_stress_tensor.zz * deviatoric_stress_tensor.zz +
         deviatoric_stress_tensor.xy * deviatoric_stress_tensor.xy +
         deviatoric_stress_tensor.xz * deviatoric_stress_tensor.xz +
         deviatoric_stress_tensor.yz * deviatoric_stress_tensor.yz;
}

/**
 * @brief Computes the strain rate tensor.
 *
 * @param strain_rate_tensor The strain rate tensor to be computed.
 * @param dv The velocity gradient dv/dr.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_strain_rate_tensor(float strain_rate_tensor[3][3], const float dv[3][3]) {

  /* Compute strain rate tensor elements. */
  strain_rate_tensor[0][0] = dv[0][0];
  strain_rate_tensor[1][1] = dv[1][1];
  strain_rate_tensor[2][2] = dv[2][2];
  strain_rate_tensor[0][1] = 0.5f * (dv[0][1] + dv[1][0]);
  strain_rate_tensor[0][2] = 0.5f * (dv[0][2] + dv[2][0]);
  strain_rate_tensor[1][0] = 0.5f * (dv[1][0] + dv[0][1]);
  strain_rate_tensor[1][2] = 0.5f * (dv[1][2] + dv[2][1]);
  strain_rate_tensor[2][0] = 0.5f * (dv[2][0] + dv[0][2]);
  strain_rate_tensor[2][1] = 0.5f * (dv[2][1] + dv[1][2]);
}

/**
 * @brief Computes the rotation rate tensor.
 *
 * @param rotation_rate_tensor The rotation rate tensor to be computed.
 * @param dv The velocity gradient dv/dr.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_rotation_rate_tensor(float rotation_rate_tensor[3][3], const float dv[3][3]) {

  /* Compute rotation rate tensor elements. */
  rotation_rate_tensor[0][0] = 0.f;
  rotation_rate_tensor[1][1] = 0.f;
  rotation_rate_tensor[2][2] = 0.f;
  rotation_rate_tensor[0][1] = 0.5f * (dv[1][0] - dv[0][1]);
  rotation_rate_tensor[0][2] = 0.5f * (dv[2][0] - dv[0][2]);
  rotation_rate_tensor[1][0] = 0.5f * (dv[0][1] - dv[1][0]);
  rotation_rate_tensor[1][2] = 0.5f * (dv[2][1] - dv[1][2]);
  rotation_rate_tensor[2][0] = 0.5f * (dv[0][2] - dv[2][0]);
  rotation_rate_tensor[2][1] = 0.5f * (dv[1][2] - dv[2][1]);
}

/**
 * @brief Computes the rotation contribution for transforming a tensor into the co-rotating frame.
 *
 * This function calculates the term M*R - R*M, where R is the rotation rate tensor,
 * and M is the tensor being rotated. This term accounts for the apparent change
 * of the tensor due to rotation of the reference frame.
 *
 * Note: Papers often make errors in the signs in this equation. For the correct
 *       equation, see Dienes 1979 for a detailed derivation, which leads to
 *       the final expression in Eqn. 4.8.
 *
 * @param rotation_term The rotation term to be computed.
 * @param rotation_rate_tensor The rotation rate tensor.
 * @param M The tensor to rotate.
 */
__attribute__((always_inline)) INLINE static void strength_compute_rotation_term(float rotation_term[3][3],
const float rotation_rate_tensor[3][3], const float M[3][3]) {

  /* Compute rotation term elements. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;
      for (int k = 0; k < 3; k++) {
        rotation_term[i][j] += M[i][k] * rotation_rate_tensor[k][j] -
                               rotation_rate_tensor[i][k] * M[k][j];
      }
    }
  }
}

/**
 * @brief Evolves a strain tensor over a time-step.
 *
 * This function updates the strain tensor according to the combined effects
 * of the strain rate and rotation contributions. After updating, the function
 * enforces symmetry of the strain tensor.
 *
 * @param strain_tensor The strain tensor to be evolved.
 * @param strain_rate_tensor The strain rate tensor.
 * @param rotation_term The rotation term.
 * @param dt The time-step over which to evolve the strain tensor.
 */
__attribute__((always_inline)) INLINE static void strength_evolve_strain_tensor(float strain_tensor[3][3],
const float strain_rate_tensor[3][3], const float rotation_term[3][3], const float dt) {

  /* Evolve strain tensor. */
  float evolved_strain_tensor[3][3];
  memcpy(evolved_strain_tensor, strain_tensor, 3 * 3 * sizeof(float));
  evolved_strain_tensor[0][0] += (strain_rate_tensor[0][0] + rotation_term[0][0]) * dt;
  evolved_strain_tensor[0][1] += (strain_rate_tensor[0][1] + rotation_term[0][1]) * dt;
  evolved_strain_tensor[0][2] += (strain_rate_tensor[0][2] + rotation_term[0][2]) * dt;
  evolved_strain_tensor[1][0] += (strain_rate_tensor[1][0] + rotation_term[1][0]) * dt;
  evolved_strain_tensor[1][1] += (strain_rate_tensor[1][1] + rotation_term[1][1]) * dt;
  evolved_strain_tensor[1][2] += (strain_rate_tensor[1][2] + rotation_term[1][2]) * dt;
  evolved_strain_tensor[2][0] += (strain_rate_tensor[2][0] + rotation_term[2][0]) * dt;
  evolved_strain_tensor[2][1] += (strain_rate_tensor[2][1] + rotation_term[2][1]) * dt;
  evolved_strain_tensor[2][2] += (strain_rate_tensor[2][2] + rotation_term[2][2]) * dt;

  /* Enforce symmetry of strain tensor. */
  strain_tensor[0][0] = evolved_strain_tensor[0][0];
  strain_tensor[1][1] = evolved_strain_tensor[1][1];
  strain_tensor[2][2] = evolved_strain_tensor[2][2];
  strain_tensor[0][1] = 0.5f * (evolved_strain_tensor[0][1] + evolved_strain_tensor[1][0]);
  strain_tensor[0][2] = 0.5f * (evolved_strain_tensor[0][2] + evolved_strain_tensor[2][0]);
  strain_tensor[1][0] = 0.5f * (evolved_strain_tensor[1][0] + evolved_strain_tensor[0][1]);
  strain_tensor[1][2] = 0.5f * (evolved_strain_tensor[1][2] + evolved_strain_tensor[2][1]);
  strain_tensor[2][0] = 0.5f * (evolved_strain_tensor[2][0] + evolved_strain_tensor[0][2]);
  strain_tensor[2][1] = 0.5f * (evolved_strain_tensor[2][1] + evolved_strain_tensor[1][2]);
}


/**
 * @brief Evolves a rotation tensor over a time-step.
 *
 * This function updates the rotation tensor according to
 *
 *     d(rotation_tensor)/dt = rotation_rate_tensor * rotation_tensor.
 *
 * For sufficiently small time steps, the rotation rate tensor is approximately
 * constant, allowing the solution to be expressed via
 *
 *     rotation_tensor(t + dt) = exp(rotation_rate_tensor * dt) * rotation_tensor(t).
 *
 * The matrix exponential can be expressed using Rodrigues' rotation formula:
 *
 *     exp(theta * N) = I + sin(theta) * N + (1 - cos(theta)) * N^2,
 *
 * where `theta` is the magnitude of the angular velocity vector derived from the
 * rotation_rate_tensor, and `N` is the corresponding unit skew-symmetric matrix.
 *
 * NOTE: The rotation tensor should remain approximately orthogonal, though
 * small numerical errors can accumulate over time, meaning the tensor could in
 * principle deviate from representing a strict rotation.
 *
 * @param rotation_tensor The rotation tensor to be evolved.
 * @param rotation_rate_tensor The rotation rate tensor.
 * @param dt The time-step over which to evolve the rotation tensor.
 */
__attribute__((always_inline)) INLINE static void
strength_evolve_rotation_tensor(float rotation_tensor[3][3], const float rotation_rate_tensor[3][3], const float dt) {

  /* Angular velocity vector from rotation rate tensor. */
  float omega[3];
  omega[0] = 0.5f * (rotation_rate_tensor[2][1] - rotation_rate_tensor[1][2]);
  omega[1] = 0.5f * (rotation_rate_tensor[0][2] - rotation_rate_tensor[2][0]);
  omega[2] = 0.5f * (rotation_rate_tensor[1][0] - rotation_rate_tensor[0][1]);
  const float omega_mag = sqrtf(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]);

  /* Rotation angle. */
  const float theta = omega_mag * dt;

  /* Return to avoid division by 0. */
  // ### Test what number to use in this if statement
  if (theta < 1e-8) {
    return;
  }

  /* Unit vector. */
  const float n[3] = {omega[0]/omega_mag, omega[1]/omega_mag, omega[2]/omega_mag};

  /* Build skew-symmetric matrix of n. */
  float N[3][3] = {
    {  0.f, -n[2],  n[1]},
    { n[2],   0.f, -n[0]},
    {-n[1],  n[0],   0.f}
  };

  /* Compute N^2. */
  float N2[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      N2[i][j] = 0.f;
      for (int k = 0; k < 3; k++) {
        N2[i][j] += N[i][k] * N[k][j];
      }
    }
  }

  /* Compute Rodrigues' rotation formula expression for exp(theta): I + sinθ N + (1 - cosθ) N^2. */
  float matrix_exp_theta[3][3] = {{0.f}};
  for (int i = 0; i < 3; i++) {
    matrix_exp_theta[i][i] += 1.f;
    for (int j = 0; j < 3; j++) {
          matrix_exp_theta[i][j] += sinf(theta) * N[i][j] + (1.f - cosf(theta)) * N2[i][j];
    }
  }

  /* Update rotation_tensor. */
  float prev_rotation_tensor[3][3];
  memcpy(prev_rotation_tensor, rotation_tensor, 3 * 3 * sizeof(float));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
    rotation_tensor[i][j] = 0.f;
      for (int k = 0; k < 3; k++) {
        rotation_tensor[i][j] += matrix_exp_theta[i][k] * prev_rotation_tensor[k][j];
      }
    }
  }
}
#endif /* SWIFT_STRENGTH_UTILITIES_H */
