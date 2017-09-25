/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H

/**
 * @file PressureEnergy/hydro_iact.h
 * @brief Pressure-Energy implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (16), (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;

  /* Get r and r inverse. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float h_inv = 1.0f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= hydro_dimension * wi + ui * wi_dx;
}

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* MISSING ! */
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {
}

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {
}

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H */
