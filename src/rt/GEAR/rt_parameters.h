/*******************************************************************************
 * This file is part of SWIFT.
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

#ifndef SWIFT_GEAR_RT_PARAMETERS_H
#define SWIFT_GEAR_RT_PARAMETERS_H

/**
 * @file src/rt/GEAR/rt_parameters.h
 * @brief Global RT parameters.
 */

extern struct rt_parameters rt_params;

/**
 * @brief Slope limiters for the reconstruction of the radiation state at the
 * particle interfaces.
 */
enum rt_slope_limiter_type {
  /*! minmod of the two extrapolations to the interface (default). */
  rt_slope_limiter_minmod = 0,
  /*! Cell limiter, then the Hopkins (2015) face limiter, as for Gizmo hydro. */
  rt_slope_limiter_gizmo,
  /*! Cell limiter only, no face limiter. */
  rt_slope_limiter_cell,
};

/**
 * Some global RT related parameters.
 */
struct rt_parameters {
  float reduced_speed_of_light;
  float reduced_speed_of_light_inverse;

  /*! Slope limiter, one of #rt_slope_limiter_type. */
  int slope_limiter;
};

#endif /* SWIFT_GEAR_RT_PARAMETERS_H */
