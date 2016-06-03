/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_IACT_GRAV_H
#define SWIFT_RUNNER_IACT_GRAV_H

/* Includes. */
#include "const.h"
#include "kernel_gravity.h"
#include "vector.h"

/**
 * @brief Gravity potential
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav(
    float r2, float *dx, struct gpart *pi, struct gpart *pj) {}

/**
 * @brief Gravity potential (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_grav(
    float *R2, float *Dx, struct gpart **pi, struct gpart **pj) {}

#endif /* SWIFT_RUNNER_IACT_GRAV_H */
