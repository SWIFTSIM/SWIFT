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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SETTERS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"
#include "hydro.h"
#include "kernel_hydro.h"

/**
 * @brief Set the gradients for the given particle to zero.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_gradients(struct part *restrict p) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chd->gradients.nabla_otimes_q[i][0] = 0.0f;
    chd->gradients.nabla_otimes_q[i][1] = 0.0f;
    chd->gradients.nabla_otimes_q[i][2] = 0.0f;
  }

  chd->gradients.v[0][0] = 0.0f;
  chd->gradients.v[0][1] = 0.0f;
  chd->gradients.v[0][2] = 0.0f;
  chd->gradients.v[1][0] = 0.0f;
  chd->gradients.v[1][1] = 0.0f;
  chd->gradients.v[1][2] = 0.0f;
  chd->gradients.v[2][0] = 0.0f;
  chd->gradients.v[2][1] = 0.0f;
  chd->gradients.v[2][2] = 0.0f;
}

/**
 * @brief Set the gradients for the given particle to the given values.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_set_diffusion_gradients(struct part *restrict p, int metal,
                                       const double gradF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients.nabla_otimes_q[metal][0] = gradF[0];
  chd->gradients.nabla_otimes_q[metal][1] = gradF[1];
  chd->gradients.nabla_otimes_q[metal][2] = gradF[2];
}

/**
 * @brief Update the diffusion gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle
 * @param metal metal specie index to update (0 <= metal <
 * GEAR_CHEMISTRY_ELEMENT_COUNT)
 * @param dF gradient of the diffusion flux
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_diffusion_gradients(struct part *restrict p, int metal,
                                          double dF[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  /* Now this is grad Z (and not grad rho_Z) */
  chd->gradients.nabla_otimes_q[metal][0] += dF[0];
  chd->gradients.nabla_otimes_q[metal][1] += dF[1];
  chd->gradients.nabla_otimes_q[metal][2] += dF[2];
}

/**
 * @brief Update the velocity gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle
 * @param dvx x velocity gradient contribution.
 * @param dvy y velocity gradient contribution.
 * @param dvz z velocity gradient contribution.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_hydro_gradients(struct part *restrict p, float dvx[3],
                                      float dvy[3], float dvz[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradients.v[0][0] += dvx[0];
  chd->gradients.v[0][1] += dvx[1];
  chd->gradients.v[0][2] += dvx[2];
  chd->gradients.v[1][0] += dvy[0];
  chd->gradients.v[1][1] += dvy[1];
  chd->gradients.v[1][2] += dvy[2];
  chd->gradients.v[2][0] += dvz[0];
  chd->gradients.v[2][1] += dvz[1];
  chd->gradients.v[2][2] += dvz[2];
}

/**
 * @brief Normalise the gradients for the given particle with the given
 * normalisation factor.
 *
 * @param p Particle.
 * @param norm Normalisation factor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_normalise_gradients(struct part *restrict p, const float norm) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chd->gradients.nabla_otimes_q[i][0] *= norm;
    chd->gradients.nabla_otimes_q[i][1] *= norm;
    chd->gradients.nabla_otimes_q[i][2] *= norm;
  }

  chd->gradients.v[0][0] *= norm;
  chd->gradients.v[0][1] *= norm;
  chd->gradients.v[0][2] *= norm;
  chd->gradients.v[1][0] *= norm;
  chd->gradients.v[1][1] *= norm;
  chd->gradients.v[1][2] *= norm;
  chd->gradients.v[2][0] *= norm;
  chd->gradients.v[2][1] *= norm;
  chd->gradients.v[2][2] *= norm;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SETTERS_H */
