/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_STRENGTH_INTERACTION_H
#define SWIFT_STRENGTH_INTERACTION_H

/**
 * @file strength/strength_interaction.h
 * @brief Conditions for pairs of particles to interact without strength.
 */

/**
 * @brief Returns whether two particles interact without strength.
 *
 * Interactions are strengthless e.g if either particle is a fluid, i.e. for
 * fluid--fluid and fluid--solid pairs.
 *
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static int
strength_is_strengthless_interaction(const struct part *restrict pi,
                                     const struct part *restrict pj) {

  /* Interactions with fluid particles */
  if ((pi->phase != mat_phase_solid) || (pj->phase != mat_phase_solid)) {
    return 1;
  }

#ifdef STRENGTH_OBJECT_IDS
  /* Interactions between particles in different objects */
  if (pi->strength_data.object_id != pj->strength_data.object_id) {
    return 1;
  }
#endif

  return 0;
}

/**
 * @brief Returns whether two particles are on either side of the interface of
 * a solid i.e. a strengthless interaction that involves a solid particle.
 *
 * @param pi First particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static int strength_is_solid_interface(
    const struct part *restrict pi, const struct part *restrict pj) {

  return strength_is_strengthless_interaction(pi, pj) &&
         ((pi->phase == mat_phase_solid) || (pj->phase == mat_phase_solid));
}

#endif /* SWIFT_STRENGTH_INTERACTION_H */