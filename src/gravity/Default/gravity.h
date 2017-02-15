/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016 Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_GRAVITY_H
#define SWIFT_DEFAULT_GRAVITY_H

#include <float.h>
#include "minmax.h"

/**
 * @brief Computes the gravity time-step of a given particle due to self-gravity
 *
 * @param gp Pointer to the g-particle data.
 */
__attribute__((always_inline)) INLINE static float
gravity_compute_timestep_self(const struct gpart* const gp) {

  const float ac2 = gp->a_grav[0] * gp->a_grav[0] +
                    gp->a_grav[1] * gp->a_grav[1] +
                    gp->a_grav[2] * gp->a_grav[2];

  const float ac = (ac2 > 0.f) ? sqrtf(ac2) : FLT_MIN;

  const float dt = sqrtf(2.f * const_gravity_eta * gp->epsilon / ac);

  return dt;
}

/**
 * @brief Prepares a g-particle for the gravity calculation
 *
 * Zeroes all the relevant arrays in preparation for the sums taking place in
 * the variaous tasks
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void gravity_init_gpart(
    struct gpart* gp) {

  /* Zero the acceleration */
  gp->a_grav[0] = 0.f;
  gp->a_grav[1] = 0.f;
  gp->a_grav[2] = 0.f;
}

/**
 * @brief Finishes the gravity calculation.
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param gp The particle to act upon
 * @param const_G Newton's constant in internal units
 */
__attribute__((always_inline)) INLINE static void gravity_end_force(
    struct gpart* gp, float const_G) {

  /* Let's get physical... */
  gp->a_grav[0] *= const_G;
  gp->a_grav[1] *= const_G;
  gp->a_grav[2] *= const_G;
}

/**
 * @brief Kick the additional variables
 *
 * @param gp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void gravity_kick_extra(
    struct gpart* gp, float dt) {}

/**
 * @brief Initialises the g-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void gravity_first_init_gpart(
    struct gpart* gp) {

  gp->time_bin = 0;
  gp->epsilon = 0.;  // MATTHIEU

  gravity_init_gpart(gp);
}

#endif /* SWIFT_DEFAULT_GRAVITY_H */
