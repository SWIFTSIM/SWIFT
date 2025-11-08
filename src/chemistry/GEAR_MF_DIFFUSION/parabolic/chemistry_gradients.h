/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_PARABOLIC_DIFFUSION_H
#define SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_PARABOLIC_DIFFUSION_H

#include "../chemistry_gradients_extrapolate.h"

/**
 * @file src/chemistry/GEAR_MFM_diffusion/chemistry_gradients.h
 * @brief Header file for parabolic diffusion gradient functions.
 */

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
    const float dx[3], const float r, const float xij_i[3], double Ui[4],
    double Uj[4]) {

  /* Metal density */
  Ui[0] = chemistry_get_comoving_metal_density(pi, metal);
  Uj[0] = chemistry_get_comoving_metal_density(pj, metal);

  /* Hyperbolic flux. Note for parabolic diffusion this is 0. */
  Ui[1] = 0.0;
  Ui[2] = 0.0;
  Ui[3] = 0.0;

  Uj[1] = 0.0;
  Uj[2] = 0.0;
  Uj[3] = 0.0;
  /* No need to check unphysical state here: they haven't been touched since
     the call to chemistry_end_density() */

  double grad_rhoZ_i[3];
  double grad_rhoZ_j[3];
  chemistry_get_metal_density_gradients(pi, metal, grad_rhoZ_i);
  chemistry_get_metal_density_gradients(pj, metal, grad_rhoZ_j);

  /* Compute interface position (relative to pj, since we don't need the
     actual position) eqn. (8)
     Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z) */
  double dUi = chemistry_gradients_extrapolate_double(grad_rhoZ_i, xij_i);
  double dUj = chemistry_gradients_extrapolate_double(grad_rhoZ_j, xij_j);

  chemistry_slope_limit_face(&Ui[0], &Uj[0], &dUi, &dUj, xij_i, xij_j, r);

  Ui[0] += dUi;
  Uj[0] += dUj;

  /* Check we have physical masses and that we are not overshooting the
     particle's mass */
  const double mi = hydro_get_mass(pi);
  const double mj = hydro_get_mass(pj);
  const double m_Zi_not_extrapolated =
      chemistry_get_metal_mass_fraction(pi, metal) * mi;
  const double m_Zj_not_extrapolated =
      chemistry_get_metal_mass_fraction(pj, metal) * mj;
  double m_Zi = Ui[0] * mi / hydro_get_comoving_density(pi);
  double m_Zj = Uj[0] * mj / hydro_get_comoving_density(pj);

  chemistry_check_unphysical_state(&m_Zi, m_Zi_not_extrapolated, mi,
                                   /*callloc=*/1, /*element*/ metal, pi->id);
  chemistry_check_unphysical_state(&m_Zj, m_Zj_not_extrapolated, mj,
                                   /*callloc=*/1, /*element*/ metal, pj->id);

  /* If the new masses have been changed, do not extrapolate, use 0th order
     reconstruction and update the state vectors */
  if (m_Zi == m_Zi_not_extrapolated) {
    Ui[0] = m_Zi_not_extrapolated * hydro_get_comoving_density(pi) / mi;
  }
  if (m_Zj == m_Zj_not_extrapolated) {
    Uj[0] = m_Zj_not_extrapolated * hydro_get_comoving_density(pj) / mj;
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

#endif /* SWIFT_CHEMISTRY_GRADIENTS_GEAR_MF_PARABOLIC_DIFFUSION_H */
