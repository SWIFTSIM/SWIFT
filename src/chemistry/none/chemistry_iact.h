/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_NONE_CHEMISTRY_IACT_H
#define SWIFT_NONE_CHEMISTRY_IACT_H

/**
 * @file none/chemistry_iact.h
 * @brief Density computation
 */

#include "cache.h"
#include "minmax.h"
#include "chemistry_struct.h"

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric version)
 *
 * @param r2 Distance squared between particles
 * @param dx Distance between particles
 * @param hi Smoothing length of i
 * @param hj Smoothing length of j
 * @param pi #part i
 * @param pj #part j
 * @param chem_data Chemistry informations
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
const struct chemistry_data *chem_data) {}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric version)
 *
 * @param r2 Distance squared between particles
 * @param dx Distance between particles
 * @param hi Smoothing length of i
 * @param hj Smoothing length of j
 * @param pi #part i
 * @param pj #part j
 * @param chem_data Chemistry informations
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
float r2, float *dx, float hi, float hj, struct part *pi, const struct part *pj,
const struct chemistry_data *chem_data) {}



#endif /* SWIFT_NONE_CHEMISTRY_IACT_H */
