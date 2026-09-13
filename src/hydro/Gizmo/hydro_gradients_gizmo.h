/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
#ifndef SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H
#define SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H

#include "fvpm_geometry.h"
#include "hydro_getters.h"
#include "hydro_setters.h"

__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  hydro_part_reset_gradients(p);

  hydro_slope_limit_cell_init(p);

#ifdef SWIFT_DEBUG_CHECKS
  /* Item 9b: reset once per gradient loop, not once per h-iteration. */
  p->geometry.psi_c_sum[0] = 0.0f;
  p->geometry.psi_c_sum[1] = 0.0f;
  p->geometry.psi_c_sum[2] = 0.0f;
  p->geometry.psi_c_abs_sum[0] = 0.0f;
  p->geometry.psi_c_abs_sum[1] = 0.0f;
  p->geometry.psi_c_abs_sum[2] = 0.0f;
#endif
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_collect(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  float wi, wj, wi_dx, wj_dx;
  float Bi[3][3];
  float Bj[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* Compute kernel of pi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  float wiBidx[3];
  if (fvpm_part_geometry_well_behaved(pi)) {
    wiBidx[0] = wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    wiBidx[1] = wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    wiBidx[2] = wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    wiBidx[0] = norm * dx[0];
    wiBidx[1] = norm * dx[1];
    wiBidx[2] = norm * dx[2];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Item 9b: role-aware Sum_j psitilde^c_j(x_i) = 0 check, independent of the
     gradient computation above (which still reads raw matrix_E). Uses only
     pi's own matrix_E_centred/centroid_offset and its own dx role, mirroring
     wiBidx, so the check does not validate itself through the same path
     first_moment/c_i were built from. */
  {
    const float dxi[3] = {dx[0] + pi->geometry.centroid_offset[0],
                          dx[1] + pi->geometry.centroid_offset[1],
                          dx[2] + pi->geometry.centroid_offset[2]};
    const float(*Bic)[3] = (const float(*)[3])pi->geometry.matrix_E_centred;
    for (int k = 0; k < 3; k++) {
      const float psi =
          -wi * (Bic[k][0] * dxi[0] + Bic[k][1] * dxi[1] + Bic[k][2] * dxi[2]);
      pi->geometry.psi_c_sum[k] += psi;
      pi->geometry.psi_c_abs_sum[k] += fabsf(psi);
    }
  }
#endif

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  drho_i[0] = dW[0] * wiBidx[0];
  drho_i[1] = dW[0] * wiBidx[1];
  drho_i[2] = dW[0] * wiBidx[2];

  dvx_i[0] = dW[1] * wiBidx[0];
  dvx_i[1] = dW[1] * wiBidx[1];
  dvx_i[2] = dW[1] * wiBidx[2];
  dvy_i[0] = dW[2] * wiBidx[0];
  dvy_i[1] = dW[2] * wiBidx[1];
  dvy_i[2] = dW[2] * wiBidx[2];
  dvz_i[0] = dW[3] * wiBidx[0];
  dvz_i[1] = dW[3] * wiBidx[1];
  dvz_i[2] = dW[3] * wiBidx[2];

  dP_i[0] = dW[4] * wiBidx[0];
  dP_i[1] = dW[4] * wiBidx[1];
  dP_i[2] = dW[4] * wiBidx[2];

  hydro_part_update_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);

  hydro_slope_limit_cell_collect(pi, pj, r);

  /* Compute kernel of pj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  float wjBjdx[3];
  if (fvpm_part_geometry_well_behaved(pj)) {

    wjBjdx[0] = wj * (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    wjBjdx[1] = wj * (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    wjBjdx[2] = wj * (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

  } else {
    const float norm = -wj_dx * r_inv;
    wjBjdx[0] = norm * dx[0];
    wjBjdx[1] = norm * dx[1];
    wjBjdx[2] = norm * dx[2];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Item 9b, pj's own role: uses only pj's matrix_E_centred/centroid_offset.
     No leading minus: d_ji - c_j = dx - c_j = dxj directly (plan section
     0.1), unlike the i-role where d_ij - c_i = -(dx + c_i) supplies the
     minus sign naturally. Copying the i-role's "-w*(B.dx)" pattern here
     would be the sign of A2's outer face-assembly term, not of this
     particle's own row-sum contribution. */
  {
    const float dxj[3] = {dx[0] - pj->geometry.centroid_offset[0],
                          dx[1] - pj->geometry.centroid_offset[1],
                          dx[2] - pj->geometry.centroid_offset[2]};
    const float(*Bjc)[3] = (const float(*)[3])pj->geometry.matrix_E_centred;
    for (int k = 0; k < 3; k++) {
      const float psi =
          wj * (Bjc[k][0] * dxj[0] + Bjc[k][1] * dxj[1] + Bjc[k][2] * dxj[2]);
      pj->geometry.psi_c_sum[k] += psi;
      pj->geometry.psi_c_abs_sum[k] += fabsf(psi);
    }
  }
#endif

  float drho_j[3], dvx_j[3], dvy_j[3], dvz_j[3], dP_j[3];

  /* Compute gradients for pj */
  /* there is no sign difference w.r.t. eqn. (6) because dx is now what we
   * want it to be */
  drho_j[0] = dW[0] * wjBjdx[0];
  drho_j[1] = dW[0] * wjBjdx[1];
  drho_j[2] = dW[0] * wjBjdx[2];

  dvx_j[0] = dW[1] * wjBjdx[0];
  dvx_j[1] = dW[1] * wjBjdx[1];
  dvx_j[2] = dW[1] * wjBjdx[2];
  dvy_j[0] = dW[2] * wjBjdx[0];
  dvy_j[1] = dW[2] * wjBjdx[1];
  dvy_j[2] = dW[2] * wjBjdx[2];
  dvz_j[0] = dW[3] * wjBjdx[0];
  dvz_j[1] = dW[3] * wjBjdx[1];
  dvz_j[2] = dW[3] * wjBjdx[2];

  dP_j[0] = dW[4] * wjBjdx[0];
  dP_j[1] = dW[4] * wjBjdx[1];
  dP_j[2] = dW[4] * wjBjdx[2];

  hydro_part_update_gradients(pj, drho_j, dvx_j, dvy_j, dvz_j, dP_j);

  hydro_slope_limit_cell_collect(pj, pi, r);

  /* Collect data for FVPM Face are checks */
  fvpm_accumulate_total_face_area_vector_and_norm(pi, pj, r2, dx, hi, hj, 1);
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_nonsym_collect(const float r2, const float dx[3],
                               const float hi, const float hj,
                               struct part *restrict pi,
                               struct part *restrict pj) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  float Bi[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
    }
  }
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  float wiBidx[3];
  if (fvpm_part_geometry_well_behaved(pi)) {
    wiBidx[0] = wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    wiBidx[1] = wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    wiBidx[2] = wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    wiBidx[0] = norm * dx[0];
    wiBidx[1] = norm * dx[1];
    wiBidx[2] = norm * dx[2];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Item 9b: role-aware Sum_j psitilde^c_j(x_i) = 0 check, non-symmetric
     variant. Only pi is updated here; pj gets its own contribution from the
     interaction in which it plays the pi role. */
  {
    const float dxi[3] = {dx[0] + pi->geometry.centroid_offset[0],
                          dx[1] + pi->geometry.centroid_offset[1],
                          dx[2] + pi->geometry.centroid_offset[2]};
    const float(*Bic)[3] = (const float(*)[3])pi->geometry.matrix_E_centred;
    for (int k = 0; k < 3; k++) {
      const float psi =
          -wi * (Bic[k][0] * dxi[0] + Bic[k][1] * dxi[1] + Bic[k][2] * dxi[2]);
      pi->geometry.psi_c_sum[k] += psi;
      pi->geometry.psi_c_abs_sum[k] += fabsf(psi);
    }
  }
#endif

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  drho_i[0] = dW[0] * wiBidx[0];
  drho_i[1] = dW[0] * wiBidx[1];
  drho_i[2] = dW[0] * wiBidx[2];

  dvx_i[0] = dW[1] * wiBidx[0];
  dvx_i[1] = dW[1] * wiBidx[1];
  dvx_i[2] = dW[1] * wiBidx[2];
  dvy_i[0] = dW[2] * wiBidx[0];
  dvy_i[1] = dW[2] * wiBidx[1];
  dvy_i[2] = dW[2] * wiBidx[2];
  dvz_i[0] = dW[3] * wiBidx[0];
  dvz_i[1] = dW[3] * wiBidx[1];
  dvz_i[2] = dW[3] * wiBidx[2];

  dP_i[0] = dW[4] * wiBidx[0];
  dP_i[1] = dW[4] * wiBidx[1];
  dP_i[2] = dW[4] * wiBidx[2];

  hydro_part_update_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);

  hydro_slope_limit_cell_collect(pi, pj, r);

  /* Collect data for FVPM Face are checks */
  fvpm_accumulate_total_face_area_vector_and_norm(pi, pj, r2, dx, hi, hj, 0);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  /* add kernel normalization to gradients */
  const float volume = p->geometry.volume;
  const float h = p->h;
  const float h_inv = 1.0f / h;
  const float ihdim = pow_dimension(h_inv);

  float norm;
  if (fvpm_part_geometry_well_behaved(p)) {
    norm = ihdim;
  } else {
    const float ihdimp1 = pow_dimension_plus_one(h_inv);
    norm = ihdimp1 * volume;
  }

  hydro_part_normalise_gradients(p, norm);

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H */
