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
#ifndef SWIFT_PLANETARY_GIZMO_HYDRO_GRADIENTS_H
#define SWIFT_PLANETARY_GIZMO_HYDRO_GRADIENTS_H

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

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];
  float drho_j[3], dvx_j[3], dvy_j[3], dvz_j[3], dP_j[3];
  hydro_part_get_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);
  hydro_part_get_gradients(pj, drho_j, dvx_j, dvy_j, dvz_j, dP_j);

  float dWi[5];
  dWi[0] = hydro_gradients_extrapolate(drho_i, xij_i);
  dWi[1] = hydro_gradients_extrapolate(dvx_i, xij_i);
  dWi[2] = hydro_gradients_extrapolate(dvy_i, xij_i);
  dWi[3] = hydro_gradients_extrapolate(dvz_i, xij_i);
  dWi[4] = hydro_gradients_extrapolate(dP_i, xij_i);

  float dWj[5];
  dWj[0] = hydro_gradients_extrapolate(drho_j, xij_j);
  dWj[1] = hydro_gradients_extrapolate(dvx_j, xij_j);
  dWj[2] = hydro_gradients_extrapolate(dvy_j, xij_j);
  dWj[3] = hydro_gradients_extrapolate(dvz_j, xij_j);
  dWj[4] = hydro_gradients_extrapolate(dP_j, xij_j);

  /* Apply the slope limiter at this interface */
  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  Wi[0] += dWi[0];
  Wi[1] += dWi[1];
  Wi[2] += dWi[2];
  Wi[3] += dWi[3];
  Wi[4] += dWi[4];

  Wj[0] += dWj[0];
  Wj[1] += dWj[1];
  Wj[2] += dWj[2];
  Wj[3] += dWj[3];
  Wj[4] += dWj[4];

  gizmo_check_physical_quantities("density", "pressure", Wi[0], Wi[1], Wi[2],
                                  Wi[3], Wi[4]);
  gizmo_check_physical_quantities("density", "pressure", Wj[0], Wj[1], Wj[2],
                                  Wj[3], Wj[4]);
}

#endif /* SWIFT_PLANETARY_GIZMO_HYDRO_GRADIENTS_H */
