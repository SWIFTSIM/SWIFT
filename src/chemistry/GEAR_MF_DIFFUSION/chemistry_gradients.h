/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GRADIENTS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GRADIENTS_H

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
 * We compute \nabla Z.
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
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
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
  if (fvpm_part_geometry_well_behaved(pi)) {
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
  if (fvpm_part_geometry_well_behaved(pj)) {
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

    chemistry_part_update_metal_mass_fraction_gradients(pi, g, dF_i);

    /* Now do the gradients of pj */
    double dF_j[3];

    /* We don't need a sign change here: both the dx and the dU
     * should switch their sign, resulting in no net change */
    dF_j[0] = dZ * psij_tilde[0];
    dF_j[1] = dZ * psij_tilde[1];
    dF_j[2] = dZ * psij_tilde[2];

    chemistry_part_update_metal_mass_fraction_gradients(pj, g, dF_j);
  }

  /*****************************************/
  /* Update velocity gradients */
  const float dv[3] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};
  const float vi_tilde[3] = {chi->filtered.rho_v[0] / chi->filtered.rho,
                             chi->filtered.rho_v[1] / chi->filtered.rho,
                             chi->filtered.rho_v[2] / chi->filtered.rho};
  const float vj_tilde[3] = {chj->filtered.rho_v[0] / chj->filtered.rho,
                             chj->filtered.rho_v[1] / chj->filtered.rho,
                             chj->filtered.rho_v[2] / chj->filtered.rho};
  const float dv_tilde[3] = {vi_tilde[0] - vj_tilde[0],
                             vi_tilde[1] - vj_tilde[1],
                             vi_tilde[2] - vj_tilde[2]};

  /* Compute velocity gradients for pi */
  float dvx_i[3], dvy_i[3], dvz_i[3], dvx_tilde_i[3], dvy_tilde_i[3],
      dvz_tilde_i[3];

  dvx_i[0] = dv[0] * psii_tilde[0];
  dvx_i[1] = dv[0] * psii_tilde[1];
  dvx_i[2] = dv[0] * psii_tilde[2];
  dvy_i[0] = dv[1] * psii_tilde[0];
  dvy_i[1] = dv[1] * psii_tilde[1];
  dvy_i[2] = dv[1] * psii_tilde[2];
  dvz_i[0] = dv[2] * psii_tilde[0];
  dvz_i[1] = dv[2] * psii_tilde[1];
  dvz_i[2] = dv[2] * psii_tilde[2];

  dvx_tilde_i[0] = dv_tilde[0] * psii_tilde[0];
  dvx_tilde_i[1] = dv_tilde[0] * psii_tilde[1];
  dvx_tilde_i[2] = dv_tilde[0] * psii_tilde[2];
  dvy_tilde_i[0] = dv_tilde[1] * psii_tilde[0];
  dvy_tilde_i[1] = dv_tilde[1] * psii_tilde[1];
  dvy_tilde_i[2] = dv_tilde[1] * psii_tilde[2];
  dvz_tilde_i[0] = dv_tilde[2] * psii_tilde[0];
  dvz_tilde_i[1] = dv_tilde[2] * psii_tilde[1];
  dvz_tilde_i[2] = dv_tilde[2] * psii_tilde[2];

  chemistry_part_update_hydro_gradients(pi, dvx_i, dvy_i, dvz_i, dvx_tilde_i,
                                        dvy_tilde_i, dvz_tilde_i);

  /* Compute velocity gradients for pj */
  float dvx_j[3], dvy_j[3], dvz_j[3], dvx_tilde_j[3], dvy_tilde_j[3],
      dvz_tilde_j[3];

  dvx_j[0] = dv[0] * psij_tilde[0];
  dvx_j[1] = dv[0] * psij_tilde[1];
  dvx_j[2] = dv[0] * psij_tilde[2];
  dvy_j[0] = dv[1] * psij_tilde[0];
  dvy_j[1] = dv[1] * psij_tilde[1];
  dvy_j[2] = dv[1] * psij_tilde[2];
  dvz_j[0] = dv[2] * psij_tilde[0];
  dvz_j[1] = dv[2] * psij_tilde[1];
  dvz_j[2] = dv[2] * psij_tilde[2];

  dvx_tilde_j[0] = dv_tilde[0] * psij_tilde[0];
  dvx_tilde_j[1] = dv_tilde[0] * psij_tilde[1];
  dvx_tilde_j[2] = dv_tilde[0] * psij_tilde[2];
  dvy_tilde_j[0] = dv_tilde[1] * psij_tilde[0];
  dvy_tilde_j[1] = dv_tilde[1] * psij_tilde[1];
  dvy_tilde_j[2] = dv_tilde[1] * psij_tilde[2];
  dvz_tilde_j[0] = dv_tilde[2] * psij_tilde[0];
  dvz_tilde_j[1] = dv_tilde[2] * psij_tilde[1];
  dvz_tilde_j[2] = dv_tilde[2] * psij_tilde[2];

  chemistry_part_update_hydro_gradients(pj, dvx_j, dvy_j, dvz_j, dvx_tilde_j,
                                        dvy_tilde_j, dvz_tilde_j);

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
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float Bi[3][3];

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
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
  if (fvpm_part_geometry_well_behaved(pi)) {
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

    chemistry_part_update_metal_mass_fraction_gradients(pi, g, dF_i);
  }

  /*****************************************/
  /* Update velocity gradients */
  const float dv[3] = {pi->v[0] - pj->v[0], pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};
  const float vi_tilde[3] = {chi->filtered.rho_v[0] / chi->filtered.rho,
                             chi->filtered.rho_v[1] / chi->filtered.rho,
                             chi->filtered.rho_v[2] / chi->filtered.rho};
  const float vj_tilde[3] = {chj->filtered.rho_v[0] / chj->filtered.rho,
                             chj->filtered.rho_v[1] / chj->filtered.rho,
                             chj->filtered.rho_v[2] / chj->filtered.rho};
  const float dv_tilde[3] = {vi_tilde[0] - vj_tilde[0],
                             vi_tilde[1] - vj_tilde[1],
                             vi_tilde[2] - vj_tilde[2]};

  /* Compute velocity gradients for pi */
  float dvx_i[3], dvy_i[3], dvz_i[3], dvx_tilde_i[3], dvy_tilde_i[3],
      dvz_tilde_i[3];

  dvx_i[0] = dv[0] * psii_tilde[0];
  dvx_i[1] = dv[0] * psii_tilde[1];
  dvx_i[2] = dv[0] * psii_tilde[2];
  dvy_i[0] = dv[1] * psii_tilde[0];
  dvy_i[1] = dv[1] * psii_tilde[1];
  dvy_i[2] = dv[1] * psii_tilde[2];
  dvz_i[0] = dv[2] * psii_tilde[0];
  dvz_i[1] = dv[2] * psii_tilde[1];
  dvz_i[2] = dv[2] * psii_tilde[2];

  dvx_tilde_i[0] = dv_tilde[0] * psii_tilde[0];
  dvx_tilde_i[1] = dv_tilde[0] * psii_tilde[1];
  dvx_tilde_i[2] = dv_tilde[0] * psii_tilde[2];
  dvy_tilde_i[0] = dv_tilde[1] * psii_tilde[0];
  dvy_tilde_i[1] = dv_tilde[1] * psii_tilde[1];
  dvy_tilde_i[2] = dv_tilde[1] * psii_tilde[2];
  dvz_tilde_i[0] = dv_tilde[2] * psii_tilde[0];
  dvz_tilde_i[1] = dv_tilde[2] * psii_tilde[1];
  dvz_tilde_i[2] = dv_tilde[2] * psii_tilde[2];

  chemistry_part_update_hydro_gradients(pi, dvx_i, dvy_i, dvz_i, dvx_tilde_i,
                                        dvy_tilde_i, dvz_tilde_i);

  /*****************************************/
  /* Collect the cell's min and max for the slope limiter. */
  chemistry_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_finalise(
  struct part *p, const struct chemistry_global_data* cd) {

  /* add kernel normalization to gradients */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  float norm;
  if (fvpm_part_geometry_well_behaved(p)) {
    const float hinvdim = pow_dimension(h_inv);
    norm = hinvdim;
  } else {
    const float hinvdimp1 = pow_dimension_plus_one(h_inv);
    norm = hinvdimp1 * p->geometry.volume;
  }

  /* Normalise the gradients */
  chemistry_part_normalise_gradients(p, norm);

  /* Limit the cell gradients */
  chemistry_slope_limit_cell(p, cd);
}

/**
 * @brief Extrapolate the given gradient over the given distance. Double
 * version.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static double
chemistry_gradients_extrapolate_double(const double gradient[3],
                                       const float dx[3]) {
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
 * @brief Metal density gradients reconstruction. Predict the value at point
 * x_ij given current values at particle positions and gradients at particle
 * positions.
 *
 * Only reconstruct U_R and U_L. We do not use a linear reconstuction for
 * grad Z_L/R, we simply use grad Z_L/R = grad Z_i/j, i.e. first order
 * reconstruction.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param metal Metal specie to update
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 * @param Ui (return) Resulting predicted and limited diffusion state of
 * particle i
 * @param Uj (return) Resulting predicted and limited diffusion state of
 * particle j
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const float dx[3], const float r, const float xij_i[3], double *Ui,
    double *Uj) {

  const double mi = hydro_get_mass(pi);
  const double mj = hydro_get_mass(pj);

  *Ui = chemistry_get_metal_mass_fraction(pi, metal);
  *Uj = chemistry_get_metal_mass_fraction(pj, metal);
  /* No need to check unphysical state here: they haven't been touched since
     the call to chemistry_end_density() */

  double m_Zi_not_extrapolated = *Ui * mi;
  double m_Zj_not_extrapolated = *Uj * mj;

  double dF_i[3];
  double dF_j[3];
  chemistry_get_metal_mass_fraction_gradients(pi, metal, dF_i);
  chemistry_get_metal_mass_fraction_gradients(pj, metal, dF_j);

  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8)
   * Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z) */
  double dUi = chemistry_gradients_extrapolate_double(dF_i, xij_i);
  double dUj = chemistry_gradients_extrapolate_double(dF_j, xij_j);

  /* Apply the slope limiter at this interface */
  chemistry_slope_limit_face(Ui, Uj, &dUi, &dUj, xij_i, xij_j, r);

  *Ui += dUi;
  *Uj += dUj;

  /* Convert to density */
  *Ui *= mi/pi->geometry.volume;
  *Uj *= mj/pj->geometry.volume;

  /* Check we have physical masses and that we are not overshooting the
     particle's mass */
  double m_Zi = *Ui * pi->geometry.volume; /* extrapolated masses */
  double m_Zj = *Uj * pj->geometry.volume;

  chemistry_check_unphysical_state(&m_Zi, m_Zi_not_extrapolated, mi,
                                   /*callloc=*/1, /*element*/ metal);
  chemistry_check_unphysical_state(&m_Zj, m_Zj_not_extrapolated, mj,
                                   /*callloc=*/1, /*element*/ metal);

  /* If the new masses have been changed, do not extrapolate, use 0th order
     reconstruction and update the state vectors */
  if (m_Zi == m_Zi_not_extrapolated) {
    *Ui = m_Zi / pi->geometry.volume;
  }
  if (m_Zj == m_Zj_not_extrapolated) {
    *Uj = m_Zj / pj->geometry.volume;
  }
}

/**
 * @brief Gradients reconstruction of the metal mass fraction of specie
 * "metal". Predict the value at point x_ij given current values at particle
 * positions and gradients at particle positions.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param metal Metal specie to update.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param cosmo The current cosmological model.
 * @param Zi (return) Particle i's metal mass fraction of metal specie "metal"
 * (in physical units).
 * @param Zj (return) Particle j's metal mass fraction of metal specie "metal"
 * (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict_Z(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const float dx[3], const struct cosmology *cosmo, double *Zi, double *Zj) {

  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

  double grad_Z_i[3], grad_Z_j[3];
  chemistry_get_metal_mass_fraction_gradients(pi, metal, grad_Z_i);
  chemistry_get_metal_mass_fraction_gradients(pj, metal, grad_Z_j);

  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8)
   * Do it this way in case dx contains periodicity corrections already */
  const float xfac = -pi->h / (pi->h + pj->h);
  const float xij_i[3] = {xfac * dx[0], xfac * dx[1], xfac * dx[2]};
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of Z_i and Z_j */
  double dZi = chemistry_gradients_extrapolate_double(grad_Z_i, xij_i);
  double dZj = chemistry_gradients_extrapolate_double(grad_Z_j, xij_j);

  /* Apply the slope limiter at this interface */
  chemistry_slope_limit_face(Zi, Zj, &dZi, &dZj, xij_i, xij_j, r);

  /* Check that we do not have unphysical values */
  if (*Zi > 1) {
    *Zi = 1;
  } else if (*Zi < 0) {
    *Zi = 0;
  }

  if (*Zi > 1) {
    *Zj = 1;
  } else if (*Zj < 0) {
    *Zj = 0;
  }

  /* Pay attention here to convert this gradient to physical units... Z is
     always physical. */
  *Zi += dZi * cosmo->a_inv;
  *Zj += dZj * cosmo->a_inv;
}

/**
 * @brief Velocity gradients reconstruction. Predict the value at point x_ij
 * given current values at particle positions and gradients at particle
 * positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param float r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 * @param Wi (return) Resulting predicted and limited state of particle i.
 * @param Wj (return) Resulting predicted and limited state of particle j.
 */
__attribute__((always_inline)) INLINE static void
chemistry_gradients_predict_hydro(struct part *restrict pi,
                                  struct part *restrict pj, const float dx[3],
                                  const float r, const float xij_i[3],
                                  float Wi[5], float Wj[5]) {

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

  /* Note: We do not reconstruct v_tilde at the interface since it is not used
     during the Riemann problem. */
}

#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_GRADIENTS_H */
