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

#include "chemistry_getters.h"
#include "chemistry_setters.h"
#include "chemistry_slope_limiters_cell.h"
#include "chemistry_slope_limiters_face.h"
#include "chemistry_unphysical.h"



/**
 * @brief Flux gradients reconstruction. Predict the value at point
 * x_ij given current values at particle positions and gradients at particle
 * positions.
 *
 * @param pi Particle i
 * @param pj Particle j
 * @param metal Metal specie to update
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 * @param Ui (return) Resulting predicted and limited diffusion flux state of
 * particle i
 * @param Uj (return) Resulting predicted and limited diffusion flux state of
 * particle j
 */
__attribute__((always_inline)) INLINE static void chemistry_gradients_predict(
    const struct part *restrict pi, const struct part *restrict pj, int metal,
    const float dx[3], const float r, const float xij_i[3], double Ui[3],
    double Uj[3]) {

  // TODO: I think that here we want to group the gradient update of the metal
  // density, so that it is consistent...
  // Find a modular way to call this (without cyclic dependencies)
  // Similarly with the extrapolate functions...
  // Probably move the chemistry_predict to a diffusion submodule...  
/* chemistry_gradients_predict( */
/*     const struct part *restrict pi, const struct part *restrict pj, int metal, */
/*     const float dx[3], const float r, const float xij_i[3], double *Ui, */
/*     double *Uj);   */

  chemistry_get_physical_flux(pi, metal, Ui);
  chemistry_get_physical_flux(pj, metal, Uj);
  /* No need to check unphysical state here: they haven't been touched since
     the call to chemistry_end_density() */

  double dFx_i[3], dFy_i[3], dFz_i[3];
  double dFx_j[3], dFy_j[3], dFz_j[3];
  chemistry_get_flux_gradients(pi, metal, dFx_i, dFy_i, dFz_i);
  chemistry_get_flux_gradients(pj, metal, dFx_j, dFy_j, dFz_j);

  /* Compute interface position (relative to pj, since we don't need the
     actual position) eqn. (8)
     Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  /* Linear reconstruction of U_R and U_L (rho*Z) */
  double dUi[3];
  dUi[0] = chemistry_gradients_extrapolate_double(dFx_i, xij_i);
  dUi[1] = chemistry_gradients_extrapolate_double(dFy_i, xij_i);
  dUi[2] = chemistry_gradients_extrapolate_double(dFz_i, xij_i);

  double dUj[3];
  dUj[0] = chemistry_gradients_extrapolate_double(dFx_j, xij_j);
  dUj[1] = chemistry_gradients_extrapolate_double(dFy_j, xij_j);
  dUj[2] = chemistry_gradients_extrapolate_double(dFz_j, xij_j);  

  chemistry_slope_limit_face_flux(Ui, Uj, dUi, dUj, xij_i, xij_j, r);

  Ui[0] += dUi[0];
  Ui[1] += dUi[1];
  Ui[2] += dUi[2];

  Uj[0] += dUj[0];
  Uj[1] += dUj[1];
  Uj[2] += dUj[2];

  /* Check we have physical fluxes */
  /* TODO */
}

  
#endif /* SWIFT_CHEMISTRY_GEAR_CHEMISTRY_HYPERBOLIC_GRADIENTS_H */
