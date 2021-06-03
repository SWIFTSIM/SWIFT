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
#ifndef SWIFT_RT_THERMOCHEMISTRY_GEAR_H
#define SWIFT_RT_THERMOCHEMISTRY_GEAR_H

/**
 * @file src/rt/GEAR/rt_thermochemistry.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * thermochemistry related functions.
 */

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 */

__attribute__((always_inline)) INLINE static void rt_do_thermochemistry(
    struct part *restrict p) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (!p->rt_data.debug_injection_done)
    error("Trying to do thermochemistry when injection step hasn't been done");
  if (!p->rt_data.debug_gradients_done)
    error("Trying to do thermochemistry when gradient step hasn't been done");
  if (!p->rt_data.debug_transport_done)
    error("Trying to do thermochemistry when transport step hasn't been done");

  p->rt_data.debug_thermochem_done += 1;
#endif
}

#endif /* SWIFT_RT_THERMOCHEMISTRY_DEBUG_H */
