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

/* #include "hydro_slope_limiters.h" */
#include "chemistry_getters.h"
#include "chemistry_setters.h"
#include "chemistry_slope_limiter.h"
#include "chemistry_unphysical.h"
#include "kernel_hydro.h"

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_init(
    struct part *p) {

  chemistry_part_reset_gradients(p);
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * We compute grad \otimes q.
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

  /* Compute psi tilde */
  float psii_tilde[3];
  if (chemistry_part_geometry_well_behaved(pi)) {
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
  if (chemistry_part_geometry_well_behaved(pj)) {
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

  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {

    double Ui, Uj;
    chemistry_part_get_diffusion_state_vector(pi, g, &Ui);
    chemistry_part_get_diffusion_state_vector(pj, g, &Uj);
    const double dU = Ui - Uj;

    /* First to the gradients of pi (i.e. grad n = nabla otimes q = grad U) */
    double dF_i[3];

    /* Compute gradients for pi */
    /* there is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dF_i[0] = dU * psii_tilde[0];
    dF_i[1] = dU * psii_tilde[1];
    dF_i[2] = dU * psii_tilde[2];

    chemistry_part_update_gradients(pi, g, dF_i);

    /* Now do the gradients of pj */
    double dF_j[3];

    /* We don't need a sign change here: both the dx and the dU
     * should switch their sign, resulting in no net change */
    dF_j[0] = dU * psij_tilde[0];
    dF_j[1] = dU * psij_tilde[1];
    dF_j[2] = dU * psij_tilde[2];

    chemistry_part_update_gradients(pj, g, dF_j);
  }
}

/**
 * @brief Gradient calculations done during the neighbour loop: non-symmetric
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

  /* Compute psi tilde */
  float psii_tilde[3];
  if (chemistry_part_geometry_well_behaved(pi)) {
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

  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {

    double Ui, Uj;
    chemistry_part_get_diffusion_state_vector(pi, g, &Ui);
    chemistry_part_get_diffusion_state_vector(pj, g, &Uj);
    const double dU = Ui - Uj;

    double dF_i[3];

    /* Compute gradients for pi */
    /* there is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dF_i[0] = dU * psii_tilde[0];
    dF_i[1] = dU * psii_tilde[1];
    dF_i[2] = dU * psii_tilde[2];

    chemistry_part_update_gradients(pi, g, dF_i);
  }
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
  if (chemistry_part_geometry_well_behaved(p)) {
    const float hinvdim = pow_dimension(h_inv);
    norm = hinvdim;
  } else {
    const float hinvdimp1 = pow_dimension_plus_one(h_inv);
    const float volume = p->chemistry_data.geometry.volume;
    norm = hinvdimp1 * volume;
  }

  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {
    chemistry_part_normalise_gradients(p, g, norm);
  }
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
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 *
 * Only reconstruct U_R and U_L. We do not use a liner reconstuction for
 * nabla_otimes_q_L/R, we simply use nabla_otimes_q_L/R = nabla_otimes_q_i/j,
 * i.e. first order reconstruction.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param Ui (return) Resulting predicted and limited diffusion state of
 * particle i
 * @param Uj (return) Resulting predicted and limited diffusion state of
 * particle j
 * @param group which photon group to use
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, double *Ui,
    double *Uj, int group, const float *dx, const float r,
    const float xij_i[3]) {

  chemistry_part_get_diffusion_state_vector(pi, group, Ui);
  chemistry_part_get_diffusion_state_vector(pj, group, Uj);
  /* No need to check unphysical state here:
   * they haven't been touched since the call
   * to rt_injection_update_photon_density */

  double dF_i[3];
  double dF_j[3];
  chemistry_part_get_gradients(pi, group, dF_i);
  chemistry_part_get_gradients(pj, group, dF_j);

  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8)
   * Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L */
  double dUi = chemistry_gradients_extrapolate(dF_i, xij_i);
  double dUj = chemistry_gradients_extrapolate(dF_j, xij_j);

  /* Apply the slope limiter at this interface */
  chemistry_slope_limit_face(Ui, Uj, &dUi, &dUj, xij_i, xij_j, r);

  *Ui += dUi;
  *Uj += dUj;

  const double m_Zi_old = *Ui*pi->chemistry_data.geometry.volume;
  const double m_Zj_old = *Uj*pj->chemistry_data.geometry.volume;
  double m_Zi = *Ui*pi->chemistry_data.geometry.volume;
  double m_Zj = *Uj*pj->chemistry_data.geometry.volume;

  /* Check and correct unphysical extrapolated states */
  chemistry_check_unphysical_state(&m_Zi, /*m_old=*/0.f, pi->mass, /*callloc=*/1);
  chemistry_check_unphysical_state(&m_Zj, /*m_old=*/0.f, pj->mass, /*callloc=*/1);

  /* If the new masses have been changed, update the state vectors */
  if (m_Zi != m_Zi_old) {
    *Ui = m_Zi / pi->chemistry_data.geometry.volume;
  }
  if (m_Zj != m_Zj_old) {
    *Uj = m_Zj / pj->chemistry_data.geometry.volume;
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_GRADIENTS_H */
