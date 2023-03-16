/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef  SWIFT_PLANETARY_HYDRO_MISC_UTILS_H
#define  SWIFT_PLANETARY_HYDRO_MISC_UTILS_H

/**
 * @file Planetary/hydro_misc_utils.h
 * @brief Miscellaneous hydro utilities for various scheme options.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Sets standard SPH or GDF rho^2-like denominator factors.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_rho_factors(
    float *rho_factor_i, float *rho_factor_j, const struct part *restrict pi, const struct part *restrict pj) {
    
#ifdef PLANETARY_GDF
  *rho_factor_i = pi->rho * pj->rho;
  *rho_factor_j = *rho_factor_i;
#else
  *rho_factor_i = pi->rho * pi->rho;
  *rho_factor_j = pj->rho * pj->rho;
#endif

}



#endif /*  SWIFT_PLANETARY_HYDRO_MISC_UTILS_H */
