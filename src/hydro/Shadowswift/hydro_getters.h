/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *                             Yolan Uyttenhove (Yolan.Uyttenhove@UGent.be)
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

#ifndef SWIFTSIM_HYDRO_SHADOWSWIFT_GETTERS_H
#define SWIFTSIM_HYDRO_SHADOWSWIFT_GETTERS_H

#include "hydro_part.h"

/**
 * @brief Get a 6-element state vector Q containing the conserved hydrodynamic
 * variables.
 *
 * @param p Particle.
 * @param Q Pointer to the array in which the result needs to be stored (of size
 * 6 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_part_get_conserved_variables(const struct part* restrict p, float* Q) {

  Q[0] = p->conserved.mass;
  Q[1] = p->conserved.momentum[0];
  Q[2] = p->conserved.momentum[1];
  Q[3] = p->conserved.momentum[2];
  Q[4] = p->conserved.energy;
  Q[5] = p->conserved.entropy;
}

/**
 * @brief Get the fluxes for the given particle.
 *
 * @param p Particle.
 * @param flux Fluxes for the particle (array of size 6 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_fluxes(
    const struct part* restrict p, float* flux) {

  flux[0] = p->flux.mass;
  flux[1] = p->flux.momentum[0];
  flux[2] = p->flux.momentum[1];
  flux[3] = p->flux.momentum[2];
  flux[4] = p->flux.energy;
  flux[5] = p->flux.entropy;
}

/**
 * @brief Returns the comoving internal energy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float
hydro_get_comoving_internal_energy(const struct part* restrict p,
                                   const struct xpart* restrict xp) {
  if (p->rho > 0.0f)
    return gas_internal_energy_from_pressure(p->rho, p->P);
  else
    return 0.f;
}

/**
 * @brief Returns the comoving entropy of a particle
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float hydro_get_comoving_entropy(
    const struct part* restrict p, const struct xpart* restrict xp) {
  if (p->rho > 0.0f) {
    return gas_entropy_from_pressure(p->rho, p->P);
  } else {
    return 0.f;
  }
}

/**
 * @brief Returns the mass of a particle
 *
 * @param p The particle of interest
 */
__attribute__((always_inline)) INLINE static float hydro_get_mass(
    const struct part* restrict p) {
  return p->conserved.mass;
}

#endif  // SWIFTSIM_HYDRO_SHADOWSWIFT_GETTERS_H
