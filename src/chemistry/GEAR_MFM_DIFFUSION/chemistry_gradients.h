/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GRADIENTS_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GRADIENTS_H

#include "chemistry_getters.h"
#include "chemistry_setters.h"
#include "chemistry_slope_limiters_cell.h"
#include "chemistry_slope_limiters_face.h"
#include "chemistry_unphysical.h"
#include "kernel_hydro.h"

/**
 * @file src/chemistry/GEAR_MFM_diffusion/chemistry_gradients.h
 * @brief Main header file for the GEAR MFM diffusion scheme gradients
 */

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_init(
    struct part *p) {

  chemistry_part_reset_gradients(p);
  chemistry_slope_limit_cell_init(p);
}

/**
 * @brief Gradient calculations done during the gradient loop
 *
 * We compute \nabla \otimes q.
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float Bi[3][3];
  float Bj[3][3];

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = chi->geometry.matrix_E[k][l];
      Bj[k][l] = chj->geometry.matrix_E[k][l];
    }
  }

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float qi = r * hi_inv;
  /* Note: factor 1/omega for psi is swallowed in matrix */
  kernel_deval(qi, &wi, &wi_dx);

  /* Compute kernel of pj */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float qj = r * hj_inv;
  kernel_deval(qj, &wj, &wj_dx);

  /*****************************************/
  /* Compute psi tilde */
  float psii_tilde[3];
  if (chemistry_geometry_well_behaved(pi)) {
    psii_tilde[0] =
        wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    psii_tilde[1] =
        wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    psii_tilde[2] =
        wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    psii_tilde[0] = norm * dx[0];
    psii_tilde[1] = norm * dx[1];
    psii_tilde[2] = norm * dx[2];
  }

  float psij_tilde[3];
  if (chemistry_geometry_well_behaved(pj)) {
    psij_tilde[0] =
        wi * (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    psij_tilde[1] =
        wi * (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    psij_tilde[2] =
        wi * (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
  } else {
    const float norm = -wj_dx * r_inv;
    psij_tilde[0] = norm * dx[0];
    psij_tilde[1] = norm * dx[1];
    psij_tilde[2] = norm * dx[2];
  }

  /*****************************************/
  /* Update diffusion gradients */
  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {
    const double Zi = chemistry_get_metal_mass_fraction(pi, g);
    const double Zj = chemistry_get_metal_mass_fraction(pj, g);
    const double dZ = Zi - Zj;

    /* First do pi (i.e. \grad n = \nabla \otimes q = \grad U) */
    double dF_i[3];

    /* There is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dF_i[0] = dZ * psii_tilde[0];
    dF_i[1] = dZ * psii_tilde[1];
    dF_i[2] = dZ * psii_tilde[2];

    chemistry_part_update_diffusion_gradients(pi, g, dF_i);

    /* Now do the gradients of pj */
    double dF_j[3];

    /* We don't need a sign change here: both the dx and the dU
     * should switch their sign, resulting in no net change */
    dF_j[0] = dZ * psij_tilde[0];
    dF_j[1] = dZ * psij_tilde[1];
    dF_j[2] = dZ * psij_tilde[2];

    chemistry_part_update_diffusion_gradients(pj, g, dF_j);
  }

  /*****************************************/
  /* Update velocity gradients */
  const float dv[5] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Compute velocity gradients for pi */
  float dvx_i[3], dvy_i[3], dvz_i[3];

  dvx_i[0] = dv[0] * psii_tilde[0];
  dvx_i[1] = dv[0] * psii_tilde[1];
  dvx_i[2] = dv[0] * psii_tilde[2];
  dvy_i[0] = dv[1] * psii_tilde[0];
  dvy_i[1] = dv[1] * psii_tilde[1];
  dvy_i[2] = dv[1] * psii_tilde[2];
  dvz_i[0] = dv[2] * psii_tilde[0];
  dvz_i[1] = dv[2] * psii_tilde[1];
  dvz_i[2] = dv[2] * psii_tilde[2];

  chemistry_part_update_hydro_gradients(pi, dvx_i, dvy_i, dvz_i);

  /* Compute velocity gradients for pj */
  float dvx_j[3], dvy_j[3], dvz_j[3];

  dvx_j[0] = dv[0] * psij_tilde[0];
  dvx_j[1] = dv[0] * psij_tilde[1];
  dvx_j[2] = dv[0] * psij_tilde[2];
  dvy_j[0] = dv[1] * psij_tilde[0];
  dvy_j[1] = dv[1] * psij_tilde[1];
  dvy_j[2] = dv[1] * psij_tilde[2];
  dvz_j[0] = dv[2] * psij_tilde[0];
  dvz_j[1] = dv[2] * psij_tilde[1];
  dvz_j[2] = dv[2] * psij_tilde[2];

  chemistry_part_update_hydro_gradients(pj, dvx_j, dvy_j, dvz_j);

  /*****************************************/
  /* Collect the cell's min and max for the slope limiter. */
  chemistry_slope_limit_cell_collect(pi, pj, r);
  chemistry_slope_limit_cell_collect(pj, pi, r);
}

/**
 * @brief Gradient calculations done during the gradient loop: non-symmetric
 * version
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
chemistry_gradients_nonsym_collect(float r2, const float *dx, float hi,
                                   float hj, struct part *restrict pi,
                                   struct part *restrict pj) {

  struct chemistry_part_data *chi = &pi->chemistry_data;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float Bi[3][3];

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = chi->geometry.matrix_E[k][l];
    }
  }

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float qi = r * hi_inv;
  /* factor 1/omega for psi is swallowed in matrix */
  kernel_deval(qi, &wi, &wi_dx);

  /*****************************************/
  /* Compute psi tilde */
  float psii_tilde[3];
  if (chemistry_geometry_well_behaved(pi)) {
    psii_tilde[0] =
        wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    psii_tilde[1] =
        wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    psii_tilde[2] =
        wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    psii_tilde[0] = norm * dx[0];
    psii_tilde[1] = norm * dx[1];
    psii_tilde[2] = norm * dx[2];
  }

  /*****************************************/
  /* Update diffusion gradients */
  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {
    const double Zi = chemistry_get_metal_mass_fraction(pi, g);
    const double Zj = chemistry_get_metal_mass_fraction(pj, g);
    const double dZ = Zi - Zj;

    double dF_i[3];

    /* Compute gradients for pi */
    /* There is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dF_i[0] = dZ * psii_tilde[0];
    dF_i[1] = dZ * psii_tilde[1];
    dF_i[2] = dZ * psii_tilde[2];

    chemistry_part_update_diffusion_gradients(pi, g, dF_i);
  }

  /*****************************************/
  /* Update velocity gradients */
  const float dv[5] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Compute velocity gradients for pi */
  float dvx_i[3], dvy_i[3], dvz_i[3];

  dvx_i[0] = dv[0] * psii_tilde[0];
  dvx_i[1] = dv[0] * psii_tilde[1];
  dvx_i[2] = dv[0] * psii_tilde[2];
  dvy_i[0] = dv[1] * psii_tilde[0];
  dvy_i[1] = dv[1] * psii_tilde[1];
  dvy_i[2] = dv[1] * psii_tilde[2];
  dvz_i[0] = dv[2] * psii_tilde[0];
  dvz_i[1] = dv[2] * psii_tilde[1];
  dvz_i[2] = dv[2] * psii_tilde[2];

  chemistry_part_update_hydro_gradients(pi, dvx_i, dvy_i, dvz_i);

  /*****************************************/
  /* Collect the cell's min and max for the slope limiter. */
  chemistry_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_finalise(
    struct part *p) {

  /* add kernel normalization to gradients */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  float norm;
  if (chemistry_geometry_well_behaved(p)) {
    const float hinvdim = pow_dimension(h_inv);
    norm = hinvdim;
  } else {
    const float hinvdimp1 = pow_dimension_plus_one(h_inv);
    norm = hinvdimp1 * chemistry_get_volume(p);
  }

  /* Normalise the gradients */
  chemistry_part_normalise_gradients(p, norm);

  /* Limit the cell gradients */
  chemistry_slope_limit_cell(p);
}

/**
 * @brief Extrapolate the given gradient over the given distance.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static double
chemistry_gradients_extrapolate(const double gradient[3], const float dx[3]) {
  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Extrapolate the given gradient over the given distance. Float version.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float
chemistry_gradients_extrapolate_float(const float gradient[3],
                                      const float dx[3]) {
  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 *
 * Only reconstruct U_R and U_L. We do not use a linear reconstuction for
 * nabla_otimes_q_L/R, we simply use nabla_otimes_q_L/R = nabla_otimes_q_i/j,
 * i.e. first order reconstruction.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param Ui (return) Resulting predicted and limited diffusion state of
 * particle i
 * @param Uj (return) Resulting predicted and limited diffusion state of
 * particle j
 * @param group which metal to use
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, double *Ui,
    double *Uj, int group, const float *dx, const float r,
    const float xij_i[3]) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  chemistry_get_diffusion_state_vector(pi, group, Ui);
  chemistry_get_diffusion_state_vector(pj, group, Uj);
  /* No need to check unphysical state here:
   * they haven't been touched since the call
   * to chemistry_end_density() */

  /* Get grad U = grad (rho*Z) = Z*grad_rho + rho*grad_Z */
  const float Delta_rho = max(pi->rho, pj->rho) - min(pi->rho, pj->rho);
  const float grad_rho[3] = {Delta_rho * dx[0] / (r * r),
                             Delta_rho * dx[1] / (r * r),
                             Delta_rho * dx[2] / (r * r)};
  double dF_i[3];
  double dF_j[3];
  chemistry_get_diffusion_gradients(pi, group, grad_rho, dF_i);
  chemistry_get_diffusion_gradients(pj, group, grad_rho, dF_j);

  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8)
   * Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z) */
  double dUi = chemistry_gradients_extrapolate(dF_i, xij_i);
  double dUj = chemistry_gradients_extrapolate(dF_j, xij_j);

  /* Apply the slope limiter at this interface */
  chemistry_slope_limit_face(Ui, Uj, &dUi, &dUj, xij_i, xij_j, r);

  *Ui += dUi;
  *Uj += dUj;

  /* Check we have physical masses and that we are not overshooting the
     particle's mass */
  const double m_Zi_old = *Ui * chi->geometry.volume;
  const double m_Zj_old = *Uj * chj->geometry.volume;
  double m_Zi = m_Zi_old;
  double m_Zj = m_Zj_old;

  /* Check and correct unphysical extrapolated states */
  chemistry_check_unphysical_state(&m_Zi, /*m_old=*/0.f, hydro_get_mass(pi),
                                   /*callloc=*/1);
  chemistry_check_unphysical_state(&m_Zj, /*m_old=*/0.f, hydro_get_mass(pj),
                                   /*callloc=*/1);

  /* If the new masses have been changed, update the state vectors */
  if (m_Zi != m_Zi_old) {
    *Ui = m_Zi / chemistry_get_volume(pi);
  }
  if (m_Zj != m_Zj_old) {
    *Uj = m_Zj / chemistry_get_volume(pj);
  }
}

/**
 * @brief Velocity gradients reconstruction. Predict the value at point x_ij
 * given current values at particle positions and gradients at particle
 * positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param Wi (return) Resulting predicted and limited state of particle i.
 * @param Wj (return) Resulting predicted and limited state of particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 */
__attribute__((always_inline)) INLINE static void
chemistry_gradients_predict_hydro(struct part *restrict pi,
                                  struct part *restrict pj, float hi, float hj,
                                  const float dx[3], float r,
                                  const float xij_i[3], float Wi[5],
                                  float Wj[5]) {

  /* Perform gradient reconstruction in space and time */
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8) */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  float dvx_i[3], dvy_i[3], dvz_i[3];
  float dvx_j[3], dvy_j[3], dvz_j[3];
  chemistry_get_hydro_gradients(pi, dvx_i, dvy_i, dvz_i);
  chemistry_get_hydro_gradients(pj, dvx_j, dvy_j, dvz_j);

  float dvi[3];
  dvi[0] = chemistry_gradients_extrapolate_float(dvx_i, xij_i);
  dvi[1] = chemistry_gradients_extrapolate_float(dvy_i, xij_i);
  dvi[2] = chemistry_gradients_extrapolate_float(dvz_i, xij_i);

  float dvj[3];
  dvj[0] = chemistry_gradients_extrapolate_float(dvx_j, xij_j);
  dvj[1] = chemistry_gradients_extrapolate_float(dvy_j, xij_j);
  dvj[2] = chemistry_gradients_extrapolate_float(dvz_j, xij_j);

  /* Apply the slope limiter at this interface */
  chemistry_slope_limit_face_hydro(Wi, Wj, dvi, dvj, xij_i, xij_j, r);

  Wi[1] += dvi[0];
  Wi[2] += dvi[1];
  Wi[3] += dvi[2];

  Wj[1] += dvj[0];
  Wj[2] += dvj[1];
  Wj[3] += dvj[2];
}

#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_GRADIENTS_H */
