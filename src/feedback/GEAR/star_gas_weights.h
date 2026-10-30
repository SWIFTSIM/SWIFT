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

#include <math.h>

/*! Distance over the star's kernel support below which a gas particle
    receives the star's light. The weights are zero there, so the runner's
    float admission of a pair at the support edge cannot change a share. */
#define STAR_GAS_WEIGHT_EDGE 0.95

/*! The light of a cone is fully delivered once its weight sum reaches this
    fraction of the sum a uniform gas gives; below it the rest escapes. */
#define STAR_GAS_CONE_FLOOR_FRACTION 0.5

/**
 * @brief Weight `(m / rho) (1 - q / q_edge)` of a gas neighbour of a star, in
 * double.
 *
 * Summed over the neighbours of one emission cone, the weights sample the
 * volume of the cone whatever the gas density. It is zero at `q = q_edge`, so
 * it is continuous across the receiver boundary.
 *
 * @param q Distance over the star's kernel support radius.
 * @param mass Mass of the gas particle.
 * @param density Density of the gas particle, in the same frame as @p mass.
 * @return The weight, 0 at and beyond #STAR_GAS_WEIGHT_EDGE or for a
 * non-positive density.
 */
__attribute__((always_inline)) INLINE static double star_gas_tent_weight(
    const double q, const double mass, const double density) {
  if (!(density > 0.0) || !(q < STAR_GAS_WEIGHT_EDGE)) return 0.0;
  return mass / density * (1.0 - q / STAR_GAS_WEIGHT_EDGE);
}

/**
 * @brief Weight sum of one emission cone at which its light is fully delivered.
 *
 * It is #STAR_GAS_CONE_FLOOR_FRACTION times the share of the cone of the
 * integral of the weight over a uniform gas, `(pi / 3) (q_edge H)^3 / n_cones`.
 *
 * @param support Comoving kernel support radius `H` of the star.
 * @param n_cones Number of emission cones.
 * @return The weight sum, 0 for a support that is not positive.
 */
__attribute__((always_inline)) INLINE static double star_gas_cone_floor(
    const double support, const int n_cones) {
  if (!(support > 0.0)) return 0.0;
  const double R = STAR_GAS_WEIGHT_EDGE * support;
  return STAR_GAS_CONE_FLOOR_FRACTION * (M_PI / 3.0) * R * R * R / n_cones;
}

/**
 * @brief Fraction of the light of the star that escapes through one cone: its
 * `1 / n_cones` minus the shares of its receivers.
 *
 * @param cone_sum Weight sum of the cone.
 * @param floor The cone floor, see star_gas_cone_floor().
 * @param n_cones Number of emission cones.
 * @return The fraction, `1 / n_cones` for an empty cone, 0 at or above the
 * floor.
 */
__attribute__((always_inline)) INLINE static double star_gas_cone_escape(
    const double cone_sum, const double floor, const int n_cones) {
  if (!(floor > cone_sum)) return cone_sum > 0.0 ? 0.0 : 1.0 / n_cones;
  return (1.0 - cone_sum / floor) / n_cones;
}

/**
 * @brief Fraction of the light of the star that escapes: the sum over the
 * cones.
 *
 * @param cone_sums Weight sum of each cone.
 * @param floor The cone floor, see star_gas_cone_floor().
 * @param n_cones Number of emission cones.
 * @return The fraction, 1 for a star with no receiver.
 */
__attribute__((always_inline)) INLINE static double star_gas_escape_fraction(
    const double cone_sums[], const double floor, const int n_cones) {
  double escape = 0.0;
  for (int c = 0; c < n_cones; c++)
    escape += star_gas_cone_escape(cone_sums[c], floor, n_cones);
  return escape;
}

/**
 * @brief Fraction of the light of the star that one receiver gets.
 *
 * The light that does not escape, `1 - sum_c escape_c`, is shared by all
 * receivers of the star in proportion to their weights, so a deficit of
 * receivers in one cone does not raise the share of the receivers there.
 *
 * @param weight Weight of the receiver.
 * @param cone_sums Weight sum of each cone.
 * @param floor The cone floor, see star_gas_cone_floor().
 * @param n_cones Number of emission cones.
 * @return The fraction, 0 if no cone holds a weight.
 */
__attribute__((always_inline)) INLINE static double star_gas_global_share(
    const double weight, const double cone_sums[], const double floor,
    const int n_cones) {
  double total = 0.0;
  for (int c = 0; c < n_cones; c++) total += cone_sums[c];
  if (!(total > 0.0)) return 0.0;
  const double delivered =
      1.0 - star_gas_escape_fraction(cone_sums, floor, n_cones);
  return delivered > 0.0 ? delivered * weight / total : 0.0;
}

#endif /* SWIFT_STAR_GAS_WEIGHTS_GEAR_H */
