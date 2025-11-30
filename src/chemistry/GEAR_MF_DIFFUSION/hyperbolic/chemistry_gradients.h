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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GRADIENTS_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GRADIENTS_H

#include "../chemistry_getters.h"
#include "../chemistry_setters.h"
#include "../chemistry_slope_limiters_cell.h"
#include "../chemistry_slope_limiters_face.h"
#include "../chemistry_unphysical.h"

/**
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param metal Metal specie to update
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i.
 * @param cosmo The #cosmology.
 * @param Ui (return) Resulting predicted and limited diffusion state of
 * particle i (in physical units).
 * @param Uj (return) Resulting predicted and limited diffusion state of
 * particle j (in physical units).
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const float dx[3], const float r, const float xij_i[3],
    const struct cosmology *cosmo, double Ui[4], double Uj[4]) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get the metal density (comoving) and the hyperbolic flux (physical) */
  Ui[0] = chemistry_get_comoving_metal_density(pi, metal);
  Ui[1] = pi->chemistry_data.flux[metal][0];
  Ui[2] = pi->chemistry_data.flux[metal][1];
  Ui[3] = pi->chemistry_data.flux[metal][2];

  Uj[0] = chemistry_get_comoving_metal_density(pj, metal);
  Uj[1] = pj->chemistry_data.flux[metal][0];
  Uj[2] = pj->chemistry_data.flux[metal][1];
  Uj[3] = pj->chemistry_data.flux[metal][2];
  /* No need to check unphysical state here: they haven't been touched since
     the call to chemistry_end_density() */

  double grad_rhoZ_i[3];
  double grad_rhoZ_j[3];
  chemistry_get_metal_density_gradients(pi, metal, grad_rhoZ_i);
  chemistry_get_metal_density_gradients(pj, metal, grad_rhoZ_j);

  /* Cell limit the metal density gradients now. The fluxes were already cell
     limited. */
  const double shoot_tol = 0.0;
  const double alpha_i = chemistry_slope_limit_quantity(
      /*gradient=*/ grad_rhoZ_i,
      /*maxr=    */ chi->limiter.maxr,
      /*value=   */ chemistry_get_comoving_metal_density(pi, metal),
      /*valmin=  */ chi->limiter.rhoZ[metal][0],
      /*valmax=  */ chi->limiter.rhoZ[metal][1],
      /*condition_number*/ pi->geometry.condition_number,
      /*pos_preserve*/ 1,
      /*shoot_tol*/ shoot_tol);
  const double alpha_j = chemistry_slope_limit_quantity(
      /*gradient=*/ grad_rhoZ_j,
      /*maxr=    */ chj->limiter.maxr,
      /*value=   */ chemistry_get_comoving_metal_density(pj, metal),
      /*valmin=  */ chj->limiter.rhoZ[metal][0],
      /*valmax=  */ chj->limiter.rhoZ[metal][1],
      /*condition_number*/ pj->geometry.condition_number,
      /*pos_preserve*/ 1,
      /*shoot_tol*/ shoot_tol);
  chemistry_slope_limit_quantity_apply(grad_rhoZ_i, alpha_i);
  chemistry_slope_limit_quantity_apply(grad_rhoZ_j, alpha_j);

  /* Since the fluxes are stored in physical units, we need to convert the
     gradients to physical units */
  double dFx_i[3], dFy_i[3], dFz_i[3];
  double dFx_j[3], dFy_j[3], dFz_j[3];
  chemistry_get_hyperbolic_flux_gradients(pi, metal, dFx_i, dFy_i, dFz_i);
  chemistry_get_hyperbolic_flux_gradients(pj, metal, dFx_j, dFy_j, dFz_j);

  /* Compute interface position (relative to pj, since we don't need the
     actual position) eqn. (8)
     Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z)
     Note for the flux extrapolation: grad_p = a^{-1} grad_c and xij_i/j_p = a
     * xij_i/j_p. Hence, the scale factors compensate and dUi/j[1-3] end up in
     phyiscal units.*/
  double dUi[4], dUj[4];
  dUi[0] = chemistry_gradients_extrapolate_double(grad_rhoZ_i, xij_i);
  dUi[1] = chemistry_gradients_extrapolate_double(dFx_i, xij_i);
  dUi[2] = chemistry_gradients_extrapolate_double(dFy_i, xij_i);
  dUi[3] = chemistry_gradients_extrapolate_double(dFz_i, xij_i);

  dUj[0] = chemistry_gradients_extrapolate_double(grad_rhoZ_j, xij_j);
  dUj[1] = chemistry_gradients_extrapolate_double(dFx_j, xij_j);
  dUj[2] = chemistry_gradients_extrapolate_double(dFy_j, xij_j);
  dUj[3] = chemistry_gradients_extrapolate_double(dFz_j, xij_j);

  chemistry_slope_limit_face(Ui, Uj, dUi, dUj, xij_i, xij_j, r);

  Ui[0] += dUi[0];
  Ui[1] += dUi[1];
  Ui[2] += dUi[2];
  Ui[3] += dUi[3];

  Uj[0] += dUj[0];
  Uj[1] += dUj[1];
  Uj[2] += dUj[2];
  Uj[3] += dUj[3];

  /* Check that we have physical masses and that we are not overshooting the
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

  /* Check that we have physical fluxes */
  double flux_i[3] = {Ui[1], Ui[2], Ui[3]};
  double flux_j[3] = {Uj[1], Uj[2], Uj[3]};
  chemistry_check_unphysical_diffusion_flux(flux_i);
  chemistry_check_unphysical_diffusion_flux(flux_j);

  /* If the new masses have been changed, do not extrapolate, use 0th order
     reconstruction and update the state vectors */
  if (m_Zi == m_Zi_not_extrapolated) {
    Ui[0] = m_Zi_not_extrapolated * hydro_get_comoving_density(pi) / mi;
  }
  if (m_Zj == m_Zj_not_extrapolated) {
    Uj[0] = m_Zj_not_extrapolated * hydro_get_comoving_density(pj) / mj;
  }
  Ui[1] = flux_i[0];
  Ui[2] = flux_i[1];
  Ui[3] = flux_i[2];
  Uj[1] = flux_j[0];
  Uj[2] = flux_j[1];
  Uj[3] = flux_j[2];

  /* Convert Ui[0] and Uj[0] (metal density) to physical units */
  Ui[0] *= cosmo->a3_inv;
  Uj[0] *= cosmo->a3_inv;
  /* The fluxes are already in physical units. No conversion needed */
}

#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_HYPERBOLIC_GRADIENTS_H */
