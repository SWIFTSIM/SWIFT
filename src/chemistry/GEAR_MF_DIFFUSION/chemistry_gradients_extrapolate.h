/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_DIFFUSION_H
#define SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_DIFFUSION_H

/**
 * @file src/chemistry/GEAR_MFM_diffusion/chemistry_gradients.h
 * @brief Header file for common gradient extrapolation functions.
 */

/**
 * @brief Extrapolate the given gradient over the given distance. Double
 * version.SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_DIFFUSION_H
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static double
chemistry_gradients_extrapolate_double(const double gradient[3],
                                       const float dx[3]) {
  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Extrapolate the given gradient over the given distance. Float version.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float
chemistry_gradients_extrapolate_float(const float gradient[3],
                                      const float dx[3]) {
  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

#endif /* SWIFT_CHEMISTRY_GRADIENTS_EXTRAPOLATE_GEAR_MF_DIFFUSION_H */
