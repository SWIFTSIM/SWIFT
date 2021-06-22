/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Coypright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_GETTERS_H
#define SWIFT_GEAR_RT_GETTERS_H

/**
 * @file src/rt/GEAR/rt_getters.h
 * @brief Getter functions for GEAR RT scheme to keep code clean and lean
 */

/**
 * @brief Get a 4-element state vector Q containing the density photon
 * quantities for a specific photon group
 *
 * @param p Particle.
 * @param group Index of photon group
 * @param Q Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void rt_part_get_density_vector(
    const struct part* restrict p, int group, float Q[4]) {

  Q[0] = p->rt_data.density[group].energy;
  Q[1] = p->rt_data.density[group].flux[0];
  Q[2] = p->rt_data.density[group].flux[1];
  Q[3] = p->rt_data.density[group].flux[2];
}

/**
 * @brief Get the gradients of energy and fluxes for a given photon group
 *
 * @param p Particle.
 * @param group Index of photon group
 * @param Q Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void rt_part_get_gradients(
    const struct part* restrict p, int group, float dE[3], float dFx[3],
    float dFy[3], float dFz[3]) {

  dE[0] = p->rt_data.gradient[group].energy[0];
  dE[1] = p->rt_data.gradient[group].energy[1];
  dE[2] = p->rt_data.gradient[group].energy[2];

  dFx[0] = p->rt_data.gradient[group].flux[0][0];
  dFx[1] = p->rt_data.gradient[group].flux[0][1];
  dFx[2] = p->rt_data.gradient[group].flux[0][2];

  dFy[0] = p->rt_data.gradient[group].flux[1][0];
  dFy[1] = p->rt_data.gradient[group].flux[1][1];
  dFy[2] = p->rt_data.gradient[group].flux[1][2];

  dFz[0] = p->rt_data.gradient[group].flux[2][0];
  dFz[1] = p->rt_data.gradient[group].flux[2][1];
  dFz[2] = p->rt_data.gradient[group].flux[2][2];
}

#endif
