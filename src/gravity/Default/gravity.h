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
#include "gravity_properties.h"
#include "minmax.h"

/**
 * @brief Returns the mass of a particle
 *
 * @param gp The particle of interest
 */
__attribute__((always_inline)) INLINE static float gravity_get_mass(
    const struct gpart* restrict gp) {

  return gp->mass;
}

/**
 * @brief Returns the softening of a particle
 *
 * @param gp The particle of interest
 */
__attribute__((always_inline)) INLINE static float gravity_get_softening(
    const struct gpart* restrict gp) {

  return gp->epsilon;
}

/**
 * @brief Returns the potential of a particle
 *
 * @param gp The particle of interest
 */
__attribute__((always_inline)) INLINE static float gravity_get_potential(
    const struct gpart* restrict gp) {

  return gp->potential;
}

/**
 * @brief Computes the gravity time-step of a given particle due to self-gravity
 *
 * We use Gadget-2's type 0 time-step criterion.
 *
 * @param gp Pointer to the g-particle data.
 * @param grav_props Constants used in the gravity scheme.
 */
__attribute__((always_inline)) INLINE static float
gravity_compute_timestep_self(const struct gpart* const gp,
                              const struct gravity_props* restrict grav_props) {

  const float ac2 = gp->a_grav[0] * gp->a_grav[0] +
                    gp->a_grav[1] * gp->a_grav[1] +
                    gp->a_grav[2] * gp->a_grav[2];

  const float ac_inv = (ac2 > 0.f) ? 1.f / sqrtf(ac2) : FLT_MAX;

  const float epsilon = gravity_get_softening(gp);

  /* Note that 0.66666667 = 2. (from Gadget) / 3. (Plummer softening) */
  const float dt = sqrtf(0.66666667f * grav_props->eta * epsilon * ac_inv);

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
  gp->potential = 0.f;

#ifdef SWIFT_DEBUG_CHECKS
  gp->num_interacted = 0;
#endif
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
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param gp The particle.
 */
__attribute__((always_inline)) INLINE static void
gravity_reset_predicted_values(struct gpart* gp) {}

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
  gp->epsilon = 0.f;

  gravity_init_gpart(gp);
}

/**
 * @brief Initialises the softening of the g-particles
 *
 * @param gp The particle to act upon.
 * @param grav_props The properties of the gravity scheme.
 */
__attribute__((always_inline)) INLINE static void gravity_init_softening(
    struct gpart* gp, const struct gravity_props* grav_props) {

  /* Note 3 is the Plummer-equivalent correction */
  gp->epsilon = grav_props->epsilon;
}

#endif /* SWIFT_DEFAULT_GRAVITY_H */
