/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GRADIENTS_GEAR_H
#define SWIFT_RT_GRADIENTS_GEAR_H

/**
 * @file src/rt/GEAR/rt_gradients.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * gradients
 */

/**
 * @brief symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pi->rt_data.debug_injection_done != 1)
    error(
        "Trying to do symmetric iact gradient when finalise injection count is "
        "%d ID %lld",
        pi->rt_data.debug_injection_done, pi->id);

  if (pj->rt_data.debug_injection_done != 1)
    error(
        "Trying to do symmetric iact gradient when finalise injection count is "
        "%d ID %lld",
        pj->rt_data.debug_injection_done, pj->id);

  pi->rt_data.debug_calls_iact_gradient_interaction += 1;

  pj->rt_data.debug_calls_iact_gradient_interaction += 1;
#endif
}

/**
 * @brief Non-symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_nonsym_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pi->rt_data.debug_injection_done != 1)
    error(
        "Trying to do nonsym iact gradients when finalise injection count is "
        "%d ID %lld",
        pi->rt_data.debug_injection_done, pi->id);

  pi->rt_data.debug_calls_iact_gradient_interaction += 1;
#endif
}

#endif /* SWIFT_RT_GRADIENT_GEAR_H */
