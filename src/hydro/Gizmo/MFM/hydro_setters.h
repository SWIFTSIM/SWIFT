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
#ifndef SWIFT_GIZMO_MFM_HYDRO_SETTERS_H
#define SWIFT_GIZMO_MFM_HYDRO_SETTERS_H

/**
 * @brief Reset the fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_fluxes(
    struct part* restrict p) {

  p->flux.momentum[0] = 0.0f;
  p->flux.momentum[1] = 0.0f;
  p->flux.momentum[2] = 0.0f;
  p->flux.energy = 0.0f;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * @param p Particle.
 * @param fluxes Fluxes accross the interface.
 * @param dx Distance between the particles that share the interface.
 */
__attribute__((always_inline)) INLINE static void hydro_part_update_fluxes_left(
    struct part* restrict p, const float* fluxes, const float* dx) {

  p->flux.momentum[0] -= fluxes[1];
  p->flux.momentum[1] -= fluxes[2];
  p->flux.momentum[2] -= fluxes[3];
  p->flux.energy -= fluxes[4];

#ifndef GIZMO_TOTAL_ENERGY
  p->flux.energy += fluxes[1] * p->v[0];
  p->flux.energy += fluxes[2] * p->v[1];
  p->flux.energy += fluxes[3] * p->v[2];
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
hydro_part_update_fluxes_right(struct part* restrict p, const float* fluxes,
                               const float* dx) {

  p->flux.momentum[0] += fluxes[1];
  p->flux.momentum[1] += fluxes[2];
  p->flux.momentum[2] += fluxes[3];
  p->flux.energy += fluxes[4];

#ifndef GIZMO_TOTAL_ENERGY
  p->flux.energy -= fluxes[1] * p->v[0];
  p->flux.energy -= fluxes[2] * p->v[1];
  p->flux.energy -= fluxes[3] * p->v[2];
#endif
}

#endif /* SWIFT_GIZMO_MFM_HYDRO_SETTERS_H */
