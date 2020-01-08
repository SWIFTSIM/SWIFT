/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_GRADIENTS_H
#define SWIFT_GIZMO_HYDRO_GRADIENTS_H

#include "hydro_getters.h"
#include "hydro_slope_limiters.h"
#include "hydro_unphysical.h"

#if defined(GRADIENTS_SPH)

#define HYDRO_GRADIENT_IMPLEMENTATION "SPH gradients (Price 2012)"
#include "hydro_gradients_sph.h"

#elif defined(GRADIENTS_GIZMO)

#define HYDRO_GRADIENT_IMPLEMENTATION "GIZMO gradients (Hopkins 2015)"
#include "hydro_gradients_gizmo.h"

#else

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/**
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part* p) {}

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
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj) {}

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
hydro_gradients_nonsym_collect(float r2, const float* dx, float hi, float hj,
                               struct part* restrict pi,
                               struct part* restrict pj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {}

#endif

/**
 * @brief Extrapolate the given gradient over the given distance.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float hydro_gradients_extrapolate(
    const float* gradient, const float* dx) {

  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, float hi, float hj,
    const float* dx, float r, const float* xij_i, float* Wi, float* Wj) {

  /* perform gradient reconstruction in space and time */
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8) */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  float gradWi[6][3], gradWj[6][3];
  hydro_part_get_gradients(pi, gradWi[0], gradWi[1], gradWi[2], gradWi[3],
                           gradWi[4], gradWi[5]);
  hydro_part_get_gradients(pj, gradWj[0], gradWj[1], gradWj[2], gradWj[3],
                           gradWj[4], gradWj[5]);

  float dWi[6], dWj[6];
  for (int i = 0; i < 6; ++i) {
    dWi[i] = hydro_gradients_extrapolate(gradWi[i], xij_i);
    dWj[i] = hydro_gradients_extrapolate(gradWj[i], xij_j);
  }

  /* Apply the slope limiter at this interface */
  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  for (int i = 0; i < 6; ++i) {
    Wi[i] += dWi[i];
    Wj[i] += dWj[i];
  }

  gizmo_check_physical_quantities("density", "pressure", Wi[0], Wi[1], Wi[2],
                                  Wi[3], Wi[4]);
  gizmo_check_physical_quantities("density", "pressure", Wj[0], Wj[1], Wj[2],
                                  Wj[3], Wj[4]);
}

#endif /* SWIFT_GIZMO_HYDRO_GRADIENTS_H */
