/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_SLOPE_LIMITERS_CELL_GEAR_H
#define SWIFT_RT_SLOPE_LIMITERS_CELL_GEAR_H

#include "rt_parameters.h"

/**
 * @file src/rt/GEAR/rt_slope_limiters_cell.h
 * @brief File containing routines concerning the cell slope
 * limiter for the GEAR RT scheme. (= first slope limiting step
 * that limits gradients such that they don't predict new extrema
 * at neighbour particle's positions)
 *
 * Used by the "gizmo" and "cell" choices of GEARRT:slope_limiter.
 */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_init(
    struct part *p) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.limiter[g].energy_density[0] = FLT_MAX;
    p->rt_data.limiter[g].energy_density[1] = -FLT_MAX;
    for (int k = 0; k < 3; k++) {
      p->rt_data.limiter[g].flux[k][0] = FLT_MAX;
      p->rt_data.limiter[g].flux[k][1] = -FLT_MAX;
    }
  }
  p->rt_data.limiter_maxr = -FLT_MAX;
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop. We need the min and max value of densities of conserved
 * quantities amongst all neighbours of particle pi.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param g index of photon group
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_collect(
    struct part *restrict pi, const struct part *restrict pj, const int g,
    const float r) {

  if (rt_params.slope_limiter == rt_slope_limiter_minmod) return;

  struct rt_part_data *rtdi = &pi->rt_data;
  const struct rt_part_data *rtdj = &pj->rt_data;

  rtdi->limiter[g].energy_density[0] = min(rtdj->radiation[g].energy_density,
                                           rtdi->limiter[g].energy_density[0]);
  rtdi->limiter[g].energy_density[1] = max(rtdj->radiation[g].energy_density,
                                           rtdi->limiter[g].energy_density[1]);
  for (int k = 0; k < 3; k++) {
    rtdi->limiter[g].flux[k][0] =
        min(rtdj->radiation[g].flux[k], rtdi->limiter[g].flux[k][0]);
    rtdi->limiter[g].flux[k][1] =
        max(rtdj->radiation[g].flux[k], rtdi->limiter[g].flux[k][1]);
  }
  rtdi->limiter_maxr = max(r, rtdi->limiter_maxr);
}

/**
 * @brief Slope-limit the given quantity. Result will be written directly
 * to float gradient[3].
 *
 * The bounds include the particle's own value, so the limiting factor stays
 * in [0, 1] and is 0 at a local extremum.
 *
 * @param gradient the gradient of the quantity
 * @param maxr maximal distance to any neighbour of the particle
 * @param value the current value of the quantity
 * @param valmin the minimal value amongst all neighbours of the quantity
 * @param valmax the maximal value amongst all neighbours of the quantity
 * @param condition_number condition number of the particle's geometry matrix
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_quantity(
    float gradient[3], const float maxr, const float value, const float valmin,
    const float valmax, const float condition_number) {

  /* Square in double: the squares of RT gradients in internal units fall
   * outside the float range. */
  const double gradtrue = sqrt((double)gradient[0] * gradient[0] +
                               (double)gradient[1] * gradient[1] +
                               (double)gradient[2] * gradient[2]) *
                          maxr;
  if (gradtrue <= 0.0) return;

  const double gradtrue_inv = 1.0 / gradtrue;
  const double gradmax = (double)max(valmax, value) - value;
  const double gradmin = value - (double)min(valmin, value);

  /* Allow twice the bound for a well-conditioned geometry matrix */
  const double beta_2 =
      (condition_number > const_gizmo_max_condition_number)
          ? const_gizmo_max_condition_number / condition_number
          : 1.0;
  const double beta = max(1.0, 2.0 * beta_2);

  const double min_temp = min(gradmax, gradmin) * gradtrue_inv * beta;
  const double alpha = min(1.0, min_temp);
  if (alpha < 1.0) {
    gradient[0] *= alpha;
    gradient[1] *= alpha;
    gradient[2] *= alpha;
  }
}

/**
 * @brief Slope limit cell gradients.
 * This is done in rt_gradients_finalise().
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell(
    struct part *p) {

  if (rt_params.slope_limiter == rt_slope_limiter_minmod) return;

  /* No neighbours: maxr is still -FLT_MAX, and scaling it can overflow. */
  if (p->rt_data.limiter_maxr <= 0.f) return;

  /* The stored gradients carry the factor 1 / gradient_scale_inv. */
  const float maxr = p->rt_data.limiter_maxr * p->rt_data.gradient_scale_inv;
  const float N_cond = p->geometry.condition_number;
  struct rt_part_data *rtd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rt_slope_limit_quantity(rtd->gradient[g].energy_density, maxr,
                            rtd->radiation[g].energy_density,
                            rtd->limiter[g].energy_density[0],
                            rtd->limiter[g].energy_density[1], N_cond);
    for (int k = 0; k < 3; k++) {
      rt_slope_limit_quantity(
          rtd->gradient[g].flux[k], maxr, rtd->radiation[g].flux[k],
          rtd->limiter[g].flux[k][0], rtd->limiter[g].flux[k][1], N_cond);
    }
  }
}
#endif /* SWIFT_RT_SLOPE_LIMITERS_CELL_GEAR_H */
