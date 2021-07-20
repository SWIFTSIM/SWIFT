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
 * limiter for the GEAR RT scheme */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_init(
    struct part* p) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.limiter[g].energy[0] = FLT_MAX;
    p->rt_data.limiter[g].energy[1] = -FLT_MAX;
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
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_cell_collect(
    struct part* pi, struct part* pj) {

  struct rt_part_data* rtdi = &pi->rt_data;
  struct rt_part_data* rtdj = &pj->rt_data;

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  for (int g = 0; g < RT_NGROUPS; g++) {
    rtdi->limiter[g].energy[0] =
        min(rtdj->density[g].energy, rtdi->limiter[g].energy[0]);
    rtdi->limiter[g].energy[1] =
        max(rtdj->density[g].energy, rtdi->limiter[g].energy[1]);

    rtdi->limiter[g].flux[0][0] =
        min(rtdj->density[g].flux[0], rtdi->limiter[g].flux[0][0]);
    rtdi->limiter[g].flux[0][1] =
        max(rtdj->density[g].flux[0], rtdi->limiter[g].flux[0][1]);
    rtdi->limiter[g].flux[1][0] =
        min(rtdj->density[g].flux[1], rtdi->limiter[g].flux[1][0]);
    rtdi->limiter[g].flux[1][1] =
        max(rtdj->density[g].flux[1], rtdi->limiter[g].flux[1][1]);
    rtdi->limiter[g].flux[2][0] =
        min(rtdj->density[g].flux[2], rtdi->limiter[g].flux[2][0]);
    rtdi->limiter[g].flux[2][1] =
        max(rtdj->density[g].flux[2], rtdi->limiter[g].flux[2][1]);
  }

  /* just use the hydro one */
  /* pi->limiter.maxr = max(r, pi->limiter.maxr); */
}

/**
 * @brief Slope-limit the given quantity. Result will be written directly
 * to float* gradient.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_quantity(
    float* gradient, const float maxr, const float value, const float valmin,
    const float valmax) {

  float gradtrue = sqrtf(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                         gradient[2] * gradient[2]);
  if (gradtrue != 0.0f) {
    gradtrue *= maxr;
    const float gradmax = valmax - value;
    const float gradmin = value - valmin;
    const float gradtrue_inv = 1.0f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
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
    float test = rtd->gradient[g].energy[0];
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].energy, maxr,
                            /*value=*/rtd->density[g].energy,
                            /*valmin=*/rtd->limiter[g].energy[0],
                            /*valmax=*/rtd->limiter[g].energy[1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[0], maxr,
                            /*value=*/rtd->density[g].flux[0],
                            /*valmin=*/rtd->limiter[g].flux[0][0],
                            /*valmax=*/rtd->limiter[g].flux[0][1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[1], maxr,
                            /*value=*/rtd->density[g].flux[1],
                            /*valmin=*/rtd->limiter[g].flux[1][0],
                            /*valmax=*/rtd->limiter[g].flux[1][1]);
    rt_slope_limit_quantity(/*gradient=*/rtd->gradient[g].flux[2], maxr,
                            /*value=*/rtd->density[g].flux[2],
                            /*valmin=*/rtd->limiter[g].flux[2][0],
                            /*valmax=*/rtd->limiter[g].flux[2][1]);
  }
}
#endif /* SWIFT_RT_SLOPE_LIMITERS_CELL_GEAR_H */
