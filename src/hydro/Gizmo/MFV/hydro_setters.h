/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_MFV_HYDRO_SETTERS_H
#define SWIFT_GIZMO_MFV_HYDRO_SETTERS_H

#include "const.h"

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_fluxes_left(
    struct part *restrict p, const float *fluxes, const float *dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->conserved.flux.mass -= fluxes[0];
  p->conserved.flux.momentum[0] -= fluxes[1];
  p->conserved.flux.momentum[1] -= fluxes[2];
  p->conserved.flux.momentum[2] -= fluxes[3];
  p->conserved.flux.energy -= fluxes[4];

#ifndef GIZMO_TOTAL_ENERGY
  const float ekin = 0.5f * (p->primitives.v[0] * p->primitives.v[0] +
                             p->primitives.v[1] * p->primitives.v[1] +
                             p->primitives.v[2] * p->primitives.v[2]);
  p->conserved.flux.energy += fluxes[1] * p->primitives.v[0];
  p->conserved.flux.energy += fluxes[2] * p->primitives.v[1];
  p->conserved.flux.energy += fluxes[3] * p->primitives.v[2];
  p->conserved.flux.energy -= fluxes[0] * ekin;
#endif
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_update_fluxes_right(struct part *restrict p, const float *fluxes,
                               const float *dx) {

  p->gravity.mflux[0] += fluxes[0] * dx[0];
  p->gravity.mflux[1] += fluxes[0] * dx[1];
  p->gravity.mflux[2] += fluxes[0] * dx[2];

  p->conserved.flux.mass += fluxes[0];
  p->conserved.flux.momentum[0] += fluxes[1];
  p->conserved.flux.momentum[1] += fluxes[2];
  p->conserved.flux.momentum[2] += fluxes[3];
  p->conserved.flux.energy += fluxes[4];

#ifndef GIZMO_TOTAL_ENERGY
  const float ekin = 0.5f * (p->primitives.v[0] * p->primitives.v[0] +
                             p->primitives.v[1] * p->primitives.v[1] +
                             p->primitives.v[2] * p->primitives.v[2]);
  p->conserved.flux.energy -= fluxes[1] * p->primitives.v[0];
  p->conserved.flux.energy -= fluxes[2] * p->primitives.v[1];
  p->conserved.flux.energy -= fluxes[3] * p->primitives.v[2];
  p->conserved.flux.energy += fluxes[0] * ekin;
#endif
}

/**
 * @brief Set the gradients for the given particle to zero.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_gradients(
    struct part *restrict p) {

  p->primitives.gradients.rho[0] = 0.0f;
  p->primitives.gradients.rho[1] = 0.0f;
  p->primitives.gradients.rho[2] = 0.0f;

  p->primitives.gradients.v[0][0] = 0.0f;
  p->primitives.gradients.v[0][1] = 0.0f;
  p->primitives.gradients.v[0][2] = 0.0f;
  p->primitives.gradients.v[1][0] = 0.0f;
  p->primitives.gradients.v[1][1] = 0.0f;
  p->primitives.gradients.v[1][2] = 0.0f;
  p->primitives.gradients.v[2][0] = 0.0f;
  p->primitives.gradients.v[2][1] = 0.0f;
  p->primitives.gradients.v[2][2] = 0.0f;

  p->primitives.gradients.P[0] = 0.0f;
  p->primitives.gradients.P[1] = 0.0f;
  p->primitives.gradients.P[2] = 0.0f;
}

/**
 * @brief Update the gradients for the given particle with the given
 * contributions.
 *
 * @param p Particle.
 * @param drho Density gradient contribution.
 * @param dvx x velocity gradient contribution.
 * @param dvy y velocity gradient contribution.
 * @param dvz z velocity gradient contribution.
 * @param dP Pressure gradient contribution.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_gradients(
    struct part *restrict p, const float *drho, const float *dvx,
    const float *dvy, const float *dvz, const float *dP) {

  p->primitives.gradients.rho[0] += drho[0];
  p->primitives.gradients.rho[1] += drho[1];
  p->primitives.gradients.rho[2] += drho[2];

  p->primitives.gradients.v[0][0] += dvx[0];
  p->primitives.gradients.v[0][1] += dvx[1];
  p->primitives.gradients.v[0][2] += dvx[2];
  p->primitives.gradients.v[1][0] += dvy[0];
  p->primitives.gradients.v[1][1] += dvy[1];
  p->primitives.gradients.v[1][2] += dvy[2];
  p->primitives.gradients.v[2][0] += dvz[0];
  p->primitives.gradients.v[2][1] += dvz[1];
  p->primitives.gradients.v[2][2] += dvz[2];

  p->primitives.gradients.P[0] += dP[0];
  p->primitives.gradients.P[1] += dP[1];
  p->primitives.gradients.P[2] += dP[2];
}

/**
 * @brief Normalise the gradients for the given particle with the given
 * normalisation factor.
 *
 * @param p Particle.
 * @param norm Normalisation factor.
 */
__attribute__((always_inline)) INLINE static void
hydro_part_normalise_gradients(struct part *restrict p, const float norm) {

  p->primitives.gradients.rho[0] *= norm;
  p->primitives.gradients.rho[1] *= norm;
  p->primitives.gradients.rho[2] *= norm;

  p->primitives.gradients.v[0][0] *= norm;
  p->primitives.gradients.v[0][1] *= norm;
  p->primitives.gradients.v[0][2] *= norm;
  p->primitives.gradients.v[1][0] *= norm;
  p->primitives.gradients.v[1][1] *= norm;
  p->primitives.gradients.v[1][2] *= norm;
  p->primitives.gradients.v[2][0] *= norm;
  p->primitives.gradients.v[2][1] *= norm;
  p->primitives.gradients.v[2][2] *= norm;

  p->primitives.gradients.P[0] *= norm;
  p->primitives.gradients.P[1] *= norm;
  p->primitives.gradients.P[2] *= norm;
}

#endif /* SWIFT_GIZMO_MFV_HYDRO_SETTERS_H */
