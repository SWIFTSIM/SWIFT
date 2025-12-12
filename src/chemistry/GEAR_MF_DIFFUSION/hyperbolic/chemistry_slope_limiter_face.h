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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SLOPE_LIMITER_FACE_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SLOPE_LIMITER_FACE_H

/**
 * @file
 *src/chemistry/GEAR_MF_DIFFUSION/hyperbolic/chemistry_slope_limiters_face.h
 * @brief Face slope limiter for the hyperbolic diffusion scheme.
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
    double Ui[4], double Uj[4], double dUi[4], double dUj[4],
    const float xij_i[3], const float xij_j[3], float r) {

  /* The Gizmo slope limiter works even better. */
  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  for (int i = 0; i < 4; i++) {
    /* chemistry_limiter_minmod(&dUi[i], &dUj[i]); */

    /* For hyperbolic diffusion, all these slope limiter perform better than
       minmod and reduce numerical diffusion */
    /* const double alphai = chemistry_limiter_mc(dUi[i], dUj[i]); */
    /* const double alphaj = chemistry_limiter_mc(dUj[i], dUi[i]); */
    /* dUi[i] *= alphai; */
    /* dUj[i] *= alphaj; */

    /* const double alphai = chemistry_limiter_superbee(dUi[i], dUj[i]); */
    /* const double alphaj = chemistry_limiter_superbee(dUj[i], dUi[i]); */
    /* dUi[i] *= alphai; */
    /* dUj[i] *= alphaj; */

    /* const double alphai = chemistry_limiter_vanLeer(dUi[i], dUj[i]); */
    /* const double alphaj = chemistry_limiter_vanLeer(dUj[i], dUi[i]); */
    /* dUi[i] *= alphai; */
    /* dUj[i] *= alphaj; */

    /* const double alphai = chemistry_limiter_koren(dUi[i], dUj[i]); */
    /* const double alphaj = chemistry_limiter_koren(dUj[i], dUi[i]); */
    /* dUi[i] *= alphai; */
    /* dUj[i] *= alphaj; */

    /* Do NOT use the positivity preserving mode: it creates anisotropy */
    dUi[i] = chemistry_slope_limit_face_quantity_double(
        Ui[i], Uj[i], Ui[i] + dUi[i], xij_i_norm, r_inv, 0);

    dUj[i] = chemistry_slope_limit_face_quantity_double(
        Uj[i], Ui[i], Uj[i] + dUj[i], xij_j_norm, r_inv, 0);
  }
}

/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * Scalar version.
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
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_face_scalar(double *Ui, double *Uj, double *dUi,
                                  double *dUj, const float xij_i[3],
                                  const float *xij_j, float r) {

  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  *dUi = chemistry_slope_limit_face_quantity_double(
      Ui[0], Uj[0], Ui[0] + dUi[0], xij_i_norm, r_inv, 0);

  *dUj = chemistry_slope_limit_face_quantity_double(
      Uj[0], Ui[0], Uj[0] + dUj[0], xij_j_norm, r_inv, 0);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SLOPE_LIMITER_FACE_H */
