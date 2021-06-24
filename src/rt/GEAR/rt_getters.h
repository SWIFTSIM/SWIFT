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

#include "rt_parameters.h"

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

/**
 * @brief compute the pressure tensor for a given state U
 *
 * @param U the state (photon energy, photon energy flux) to use
 * @param flux the resulting flux F(U) of the hyperbolic conservation law
 */
__attribute__((always_inline)) INLINE static void rt_get_pressure_tensor(const float U[4], float pressure_tensor[3][3]){
 
  const double c_red = rt_params.reduced_speed_of_light;
  const float normF = sqrtf(U[1]*U[1] + U[2]*U[2] + U[3] * U[3]);
  float f = normF / (c_red * U[0]);
  const float f2 = f * f;
  /* Because we use reduced speed of light, we may find f > 1,
   * which will lead to negative values in sqrt. Handle this. */
  const float rootterm = max((4.f - 3.f * f2), 0.f);
  const float chi = (3.f + 4.f * f2) / (5.f + 2. * sqrtf(rootterm));
  const float normF_inv = 1.f/normF;
  const float n[3] = {U[1]*normF_inv, U[2]*normF_inv, U[3]*normF_inv};

  const float temp = 0.5f * (3.f * chi - 1.f);
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      pressure_tensor[i][j] = temp * n[i] * n[j];
    }
  }

  const float temp2 = 0.5f * (1.f - chi);
  pressure_tensor[0][0] += temp2;
  pressure_tensor[1][1] += temp2;
  pressure_tensor[2][2] += temp2;

  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      pressure_tensor[i][j] *= U[0];
    }
  }
}

/**
 * @brief compute the flux of the hyperbolic conservation law for a given
 * state U
 *
 * @param U the state (photon energy, photon energy flux) to use
 * @param flux the resulting flux F(U) of the hyperbolic conservation law
 */
__attribute__((always_inline)) INLINE static void rt_get_hyperbolic_flux(const float U[4], float flux[4][3]){

  flux[0][0] = U[1];
  flux[0][1] = U[2];
  flux[0][2] = U[3];

  float pressure_tensor[3][3];
  rt_get_pressure_tensor(U, pressure_tensor);

  const float c_red = rt_params.reduced_speed_of_light;
  const float c2 = c_red * c_red;
  flux[1][0] = pressure_tensor[0][0] * c2;
  flux[1][1] = pressure_tensor[0][1] * c2;
  flux[1][2] = pressure_tensor[0][2] * c2;
  flux[2][0] = pressure_tensor[1][0] * c2;
  flux[2][1] = pressure_tensor[1][1] * c2;
  flux[2][2] = pressure_tensor[1][2] * c2;
  flux[3][0] = pressure_tensor[2][0] * c2;
  flux[3][1] = pressure_tensor[2][1] * c2;
  flux[3][2] = pressure_tensor[2][2] * c2;
}

/**
 * @brief Compute the time-step length for an RT step.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the gravity kick (internal units).
 */
__attribute__((always_inline)) INLINE static double rt_get_part_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology *cosmo) {

  if (with_cosmology) {
    error("GEAR RT with cosmology not implemented yet! :(");
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}
#endif
