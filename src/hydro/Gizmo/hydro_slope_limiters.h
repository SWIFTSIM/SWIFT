/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_HYDRO_SLOPE_LIMITERS_H
#define SWIFT_HYDRO_SLOPE_LIMITERS_H

#define PER_FACE_LIMITER
#define CELL_WIDE_LIMITER

#ifdef PER_FACE_LIMITER
#include "hydro_slope_limiters_face.h"
#else

__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, float *xij_i, float *xij_j,
    float r) {}

#endif

#ifdef CELL_WIDE_LIMITER
#include "hydro_slope_limiters_cell.h"
#else

__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part *p) {}

__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj, float r) {}

__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part *p) {}

#endif

#endif  // SWIFT_HYDRO_SLOPE_LIMITERS_H
