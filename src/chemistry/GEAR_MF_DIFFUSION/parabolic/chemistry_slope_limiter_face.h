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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_SLOPE_LIMITERS_FACE_H
#define SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_SLOPE_LIMITERS_FACE_H

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/parabolic/chemistry_slope_limiters_face.h
 * @brief Face slope limiter for the parabolic diffusion scheme.
 **/

/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Ui Chemistry variables of particle i.
 * @param Uj Chemistry variables of particle j.
 * @param dUi Difference between the chemistry variables of particle i at the
 * position of particle i and at the interface position.
 * @param dUj Difference between the chemistry variables of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_face(
    double *Ui, double *Uj, double *dUi, double *dUj, const float xij_i[3],
    const float *xij_j, float r) {
  /* chemistry_limiter_minmod(dUi, dUj); */

  /* The Gizmo slope limiter works better. */
  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  *dUi = chemistry_slope_limit_face_quantity_double(
      Ui[0], Uj[0], Ui[0] + dUi[0], xij_i_norm, r_inv, 1);

  *dUj = chemistry_slope_limit_face_quantity_double(
      Uj[0], Ui[0], Uj[0] + dUj[0], xij_j_norm, r_inv, 1);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_PARABOLIC_DIFFUSION_SLOPE_LIMITERS_FACE_H */
