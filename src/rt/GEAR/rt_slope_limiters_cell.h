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

/**
 * @file src/rt/GEAR/rt_slope_limiters_cell.h
 * @brief File containing routines concerning the cell slope
 * limiter for the GEAR RT scheme. (= fist slope limiting step
 * that limits gradients such that they don't predict new extrema
 * at neighbour praticle's positions )
 *
 * The Gizmo-style slope limiter doesn't help for RT problems for
 * now, so nothing in this file should actually be called.
 * */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_init(
    struct part* p) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.limiter[g].energy_density[0] = FLT_MAX;
    p->rt_data.limiter[g].energy_density[1] = -FLT_MAX;
    p->rt_data.limiter[g].flux[0][0] = FLT_MAX;
    p->rt_data.limiter[g].flux[0][1] = -FLT_MAX;
    p->rt_data.limiter[g].flux[1][0] = FLT_MAX;
    p->rt_data.limiter[g].flux[1][1] = -FLT_MAX;
    p->rt_data.limiter[g].flux[2][0] = FLT_MAX;
    p->rt_data.limiter[g].flux[2][1] = -FLT_MAX;
  }

  /* just use the hydro one */
  /* p->limiter.maxr = -FLT_MAX; */
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop. We need the min and max value of densities of conserved
 * quantities amongst all neighbours of particle pi.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param g index of photon group
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_collect(
    struct part* restrict pi, struct part* restrict pj, int g) {

  struct rt_part_data* rtdi = &pi->rt_data;
  struct rt_part_data* rtdj = &pj->rt_data;

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  rtdi->limiter[g].energy_density[0] = min(rtdj->radiation[g].energy_density,
                                           rtdi->limiter[g].energy_density[0]);
  rtdi->limiter[g].energy_density[1] = max(rtdj->radiation[g].energy_density,
                                           rtdi->limiter[g].energy_density[1]);

  rtdi->limiter[g].flux[0][0] =
      min(rtdj->radiation[g].flux[0], rtdi->limiter[g].flux[0][0]);
  rtdi->limiter[g].flux[0][1] =
      max(rtdj->radiation[g].flux[0], rtdi->limiter[g].flux[0][1]);
  rtdi->limiter[g].flux[1][0] =
      min(rtdj->radiation[g].flux[1], rtdi->limiter[g].flux[1][0]);
  rtdi->limiter[g].flux[1][1] =
      max(rtdj->radiation[g].flux[1], rtdi->limiter[g].flux[1][1]);
  rtdi->limiter[g].flux[2][0] =
      min(rtdj->radiation[g].flux[2], rtdi->limiter[g].flux[2][0]);
  rtdi->limiter[g].flux[2][1] =
      max(rtdj->radiation[g].flux[2], rtdi->limiter[g].flux[2][1]);
  /* just use the hydro one */
  /* pi->limiter.maxr = max(r, pi->limiter.maxr); */
}

/**
 * @brief Slope-limit the given quantity. Result will be written directly
 * to float gradient[3].
 *
 * @param gradient the gradient of the quantity
 * @param maxr maximal distance to any neighbour of the particle
 * @param value the current value of the quantity
 * @param valmin the minimal value amongst all neighbours of the quantity
 * @param valmax the maximal value amongst all neighbours of the quantity
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_quantity(
    float gradient[3], const float maxr, const float value, const float valmin,
    const float valmax) {

  float gradtrue = sqrtf(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                         gradient[2] * gradient[2]);
  if (gradtrue != 0.0f) {
    gradtrue *= maxr;
    const float gradtrue_inv = 1.0f / gradtrue;
    const float gradmax = valmax - value;
    const float gradmin = valmin - value;
    const float beta = 1.f; /* TODO: test for best value here. For now, take
                                stability over diffusivity. */
    const float min_temp =
        min(gradmax * gradtrue_inv, gradmin * gradtrue_inv) * beta;
    const float alpha = min(1.f, min_temp);
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
    struct part* p) {

  const float maxr = p->limiter.maxr;
  struct rt_part_data* rtd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].energy_density,
                            /*maxr=    */ maxr,
                            /*value=   */ rtd->radiation[g].energy_density,
                            /*valmin=  */ rtd->limiter[g].energy_density[0],
                            /*valmax=  */ rtd->limiter[g].energy_density[1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[0],
                            /*maxr=    */ maxr,
                            /*value=   */ rtd->radiation[g].flux[0],
                            /*valmin=  */ rtd->limiter[g].flux[0][0],
                            /*valmax=  */ rtd->limiter[g].flux[0][1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[1],
                            /*maxr=    */ maxr,
                            /*value=   */ rtd->radiation[g].flux[1],
                            /*valmin=  */ rtd->limiter[g].flux[1][0],
                            /*valmax=  */ rtd->limiter[g].flux[1][1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[2],
                            /*maxr=    */ maxr,
                            /*value=   */ rtd->radiation[g].flux[2],
                            /*valmin=  */ rtd->limiter[g].flux[2][0],
                            /*valmax=  */ rtd->limiter[g].flux[2][1]);
  }
}
#endif /* SWIFT_RT_SLOPE_LIMITERS_CELL_GEAR_H */
