/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwn Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_GEAR_CHEMISTRY_SETTERS_H
#define SWIFT_GEAR_CHEMISTRY_SETTERS_H

/**
 * @brief Set the gradients for the given particle to zero.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_part_reset_gradients(
    struct part* restrict p) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {      
    chd->gradients[i].nabla_otimes_q[0] = 0.0f;
    chd->gradients[i].nabla_otimes_q[1] = 0.0f;
    chd->gradients[i].nabla_otimes_q[2] = 0.0f;
  }
}

/**
 * @brief Set the gradients for the given particle to the given values.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_part_set_gradients(
   struct part* restrict p, int g, const double gradF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients[g].nabla_otimes_q[0] = gradF[0];
  chd->gradients[g].nabla_otimes_q[1] = gradF[1];
  chd->gradients[g].nabla_otimes_q[2] = gradF[2];
}

/**
 * @brief Update the gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle
 * @param g photon group index to update (0 <= g < RT_NGROUPS)
 * @param dE energy density gradient
 * @param dF gradient of the flux component
 */
__attribute__((always_inline)) INLINE static void chemistry_part_update_gradients(
    struct part *restrict p, int g, double dF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients[g].nabla_otimes_q[0] += dF[0];
  chd->gradients[g].nabla_otimes_q[1] += dF[1];
  chd->gradients[g].nabla_otimes_q[2] += dF[2];
}

/**
 * @brief Normalise the gradients for the given particle with the given
 * normalisation factor.
 *
 * @param p Particle.
 * @param g photon group index to update (0 <= g < RT_NGROUPS)
 * @param norm Normalisation factor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_normalise_gradients(struct part* restrict p, int g, const float norm) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients[g].nabla_otimes_q[0] *= norm;
  chd->gradients[g].nabla_otimes_q[1] *= norm;
  chd->gradients[g].nabla_otimes_q[2] *= norm;
}


#endif /* SWIFT_GEAR_CHEMISTRY_SETTERS_H */
