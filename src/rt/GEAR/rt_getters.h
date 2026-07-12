/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
 * @brief Get the comoving radiation energy densities of a particle.
 *
 * @param p Particle.
 * @param E (return) Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void
rt_part_get_comoving_radiation_energy_density(const struct part *restrict p,
                                              float E[RT_NGROUPS]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    E[g] = p->rt_data.radiation[g].energy_density;
  }
}

/**
 * @brief Get the physical radiation energy densities of a particle
 *
 * @param p Particle.
 * @param E (return) Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void
rt_part_get_physical_radiation_energy_density(const struct part *restrict p,
                                              float E[RT_NGROUPS],
                                              const struct cosmology *cosmo) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    E[g] = cosmo->a3_inv * p->rt_data.radiation[g].energy_density;
  }
}

/**
 * @brief Get a 4-element state vector U containing the radiation energy
 * density and fluxes for a specific photon group.
 *
 * @param p Particle.
 * @param group Index of photon group
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void
rt_part_get_radiation_state_vector(const struct part *restrict p, int group,
                                   float U[4]) {

  U[0] = p->rt_data.radiation[group].energy_density;
  U[1] = p->rt_data.radiation[group].flux[0];
  U[2] = p->rt_data.radiation[group].flux[1];
  U[3] = p->rt_data.radiation[group].flux[2];
}

/**
 * @brief Get the gradients of energy density and fluxes for a given photon
 * group
 *
 * @param p Particle.
 * @param group Index of photon group
 * @param dE (return) Array to write energy density gradient into
 * @param dFx (return) Array to write flux x component gradient into
 * @param dFy (return) Array to write flux y component gradient into
 * @param dFz (return) Array to write flux z component gradient into
 */
__attribute__((always_inline)) INLINE static void rt_part_get_gradients(
    const struct part *restrict p, int group, float dE[3], float dFx[3],
    float dFy[3], float dFz[3]) {

  dE[0] = p->rt_data.gradient[group].energy_density[0];
  dE[1] = p->rt_data.gradient[group].energy_density[1];
  dE[2] = p->rt_data.gradient[group].energy_density[2];

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
 * @brief compute the pressure tensor for a given radiation state U
 *
 * @param U the state (radiation energy density, radiation flux) to use
 * @param Fnorm the norm of the radiation flux
 * @param pressure_tensor (return) 3x3 array to write resulting Eddington
 * pressure tensor into
 */
__attribute__((always_inline)) INLINE static void rt_get_pressure_tensor(
    const float U[4], const float Fnorm, float pressure_tensor[3][3]) {

  /* We may encounter zero flux even with nonzero energy.
   * Also even with nonzero flux, the norm may round down
   * to exactly zero, so exit early if that is the case. */
  if (Fnorm == 0.f) {
    const float diagonal_element = U[0] / 3.f;
    pressure_tensor[0][0] = diagonal_element;
    pressure_tensor[0][1] = 0.f;
    pressure_tensor[0][2] = 0.f;
    pressure_tensor[1][0] = 0.f;
    pressure_tensor[1][1] = diagonal_element;
    pressure_tensor[1][2] = 0.f;
    pressure_tensor[2][0] = 0.f;
    pressure_tensor[2][1] = 0.f;
    pressure_tensor[2][2] = diagonal_element;
    return;
  }

  /* f mustn't be > 1. This may happen because of the reduced speed of light.
   * Energy density U[0] is nonzero at this point, or this function wouldn't
   * have been called. */
  const float f =
      min(1.f, rt_params.reduced_speed_of_light_inverse * Fnorm / U[0]);
  const float f2 = f * f;
  const float rootterm = 4.f - 3.f * f2;
  const float chi = (3.f + 4.f * f2) / (5.f + 2.f * sqrtf(rootterm));

  /* get unit vector n */
  const float Fnorm_inv = 1.f / Fnorm;
  const float n[3] = {U[1] * Fnorm_inv, U[2] * Fnorm_inv, U[3] * Fnorm_inv};

  const float temp = 0.5f * (3.f * chi - 1.f);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pressure_tensor[i][j] = temp * n[i] * n[j];
    }
  }

  const float temp2 = 0.5f * (1.f - chi);
  pressure_tensor[0][0] += temp2;
  pressure_tensor[1][1] += temp2;
  pressure_tensor[2][2] += temp2;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pressure_tensor[i][j] *= U[0];
    }
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pressure_tensor[0][0] != pressure_tensor[0][0]) {
    message(
        "Found NaNs in pressure tensor: 1/c %.3e |"
        " |F| %.3e f %.3e f^2 %.3e root %.3e chi %.3e |"
        " n %.3e %.3e %.3e | U %.3e %.3e %.3e |"
        " temp %.3e %.3e",
        rt_params.reduced_speed_of_light_inverse, Fnorm, f, f2, rootterm, chi,
        n[0], n[1], n[2], U[1], U[2], U[3], temp, temp2);
  }
#endif
}

/**
 * @brief compute the flux of the hyperbolic conservation law for a given
 * state U
 *
 * @param U the state (radiation energy density, radiation flux) to use
 * @param Fnorm the norm of the radiation flux
 * @param hypflux (return) resulting flux F(U) of the hyperbolic conservation
 * law
 */
__attribute__((always_inline)) INLINE static void rt_get_hyperbolic_flux(
    const float U[4], const float Fnorm, float hypflux[4][3]) {

  if (U[0] == 0.f) {
    /* At this point, the state U has been corrected to not contain
     * unphysical values. If we encounter this situation, it means
     * that the fluxes are zero as well, meaning that when we compute
     * 1/|F| we get infinities. So skip this. The pressure tensor is
     * P_ij = D_ij * E_i anyway. */

    for (int i = 0; i < 4; i++) {
      hypflux[i][0] = 0.f;
      hypflux[i][1] = 0.f;
      hypflux[i][2] = 0.f;
    }
    return;
  }

  hypflux[0][0] = U[1];
  hypflux[0][1] = U[2];
  hypflux[0][2] = U[3];

  float pressure_tensor[3][3];
  rt_get_pressure_tensor(U, Fnorm, pressure_tensor);

  const float c_red = rt_params.reduced_speed_of_light;
  const float c2 = c_red * c_red;
  hypflux[1][0] = pressure_tensor[0][0] * c2;
  hypflux[1][1] = pressure_tensor[0][1] * c2;
  hypflux[1][2] = pressure_tensor[0][2] * c2;
  hypflux[2][0] = pressure_tensor[1][0] * c2;
  hypflux[2][1] = pressure_tensor[1][1] * c2;
  hypflux[2][2] = pressure_tensor[1][2] * c2;
  hypflux[3][0] = pressure_tensor[2][0] * c2;
  hypflux[3][1] = pressure_tensor[2][1] * c2;
  hypflux[3][2] = pressure_tensor[2][2] * c2;
}

#endif
