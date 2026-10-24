/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_STAR_GAS_WEIGHTS_GEAR_H
#define SWIFT_STAR_GAS_WEIGHTS_GEAR_H

/**
 * @file src/feedback/GEAR/star_gas_weights.h
 * @brief Weights that split a star's output between its gas neighbours.
 */

#include "inline.h"
#include "kernel_hydro.h"

/**
 * @brief Volume weight `(m / rho) W(r, h)` of a gas neighbour of a star, in
 * double.
 *
 * Summed over the neighbours, the weights sample the kernel-weighted volume
 * around the star, so each direction gets its geometric share whatever the
 * gas density on that side. The factor `h^-d`, common to all neighbours of
 * one star, is left out.
 *
 * @param u Distance over the star's smoothing length.
 * @param mass Mass of the gas particle.
 * @param density Density of the gas particle, in the same frame as @p mass.
 * @return The weight, 0 beyond the kernel support or for a non-positive
 * density.
 */
__attribute__((always_inline)) INLINE static double star_gas_volume_weight(
    const double u, const double mass, const double density) {
  if (!(density > 0.0)) return 0.0;
  double W;
  kernel_eval_double(u, &W);
  return mass / density * W;
}

#endif /* SWIFT_STAR_GAS_WEIGHTS_GEAR_H */
