/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SETTERS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"
#include "hydro.h"
#include "kernel_hydro.h"

/**
 * @brief Update the flux gradients for the given particle to the given values.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param gradF Metal mass fraction gradient (of size 3) to set.
 * @param dFx gradient of the x direction flux component
 * @param dFy gradient of the y direction flux component
 * @param dFz gradient of the z direction flux component
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_flux_gradients(struct part *restrict p, int metal,
                                     float dFx[3], float dFy[3], float dFz[3]) {

  struct chemistry_part_data *chd = &p->chemistry_data;

  chd->gradient[metal].F_diff[0][0] += dFx[0];
  chd->gradient[metal].F_diff[0][1] += dFx[1];
  chd->gradient[metal].F_diff[0][2] += dFx[2];

  chd->gradient[metal].F_diff[1][0] += dFy[0];
  chd->gradient[metal].F_diff[1][1] += dFy[1];
  chd->gradient[metal].F_diff[1][2] += dFy[2];

  chd->gradient[metal].F_diff[2][0] += dFz[0];
  chd->gradient[metal].F_diff[2][1] += dFz[1];
  chd->gradient[metal].F_diff[2][2] += dFz[2];
}
#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_SETTERS_H */
