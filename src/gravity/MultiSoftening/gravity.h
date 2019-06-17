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

/* Local includes. */
#include "cosmology.h"
#include "gravity_properties.h"
#include "kernel_gravity.h"
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
 * @brief Returns the current co-moving softening of a particle
 *
 * @param gp The particle of interest
 * @param grav_props The global gravity properties.
 */
__attribute__((always_inline)) INLINE static float gravity_get_softening(
    const struct gpart* gp, const struct gravity_props* restrict grav_props) {

  return grav_props->epsilon_cur;
}

/**
 * @brief Add a contribution to this particle's potential.
 *
 * Here we do nothing as this version does not accumulate potential.
 *
 * @param gp The particle.
 * @param pot The contribution to add.
 */
__attribute__((always_inline)) INLINE static void
gravity_add_comoving_potential(struct gpart* restrict gp, float pot) {}

/**
 * @brief Returns the comoving potential of a particle.
 *
 * This returns 0 as this flavour of gravity does not compute the
 * particles' potential.
 *
 * @param gp The particle of interest
 */
__attribute__((always_inline)) INLINE static float
gravity_get_comoving_potential(const struct gpart* restrict gp) {

  return 0.f;
}

/**
 * @brief Returns the physical potential of a particle
 *
 * This returns 0 as this flavour of gravity does not compute the
 * particles' potential.
 *
 * @param gp The particle of interest.
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static float
gravity_get_physical_potential(const struct gpart* restrict gp,
                               const struct cosmology* cosmo) {

  return 0.f;
}

/**
 * @brief Computes the gravity time-step of a given particle due to self-gravity
 *
 * We use Gadget-2's type 0 time-step criterion.
 *
 * @param gp Pointer to the g-particle data.
 * @param a_hydro The accelerations coming from the hydro scheme (can be 0).
 * @param grav_props Constants used in the gravity scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float
gravity_compute_timestep_self(const struct gpart* const gp,
                              const float a_hydro[3],
                              const struct gravity_props* restrict grav_props,
                              const struct cosmology* cosmo) {

  /* Get physical acceleration (gravity contribution) */
  float a_phys_x = gp->a_grav[0] * cosmo->a_factor_grav_accel;
  float a_phys_y = gp->a_grav[1] * cosmo->a_factor_grav_accel;
  float a_phys_z = gp->a_grav[2] * cosmo->a_factor_grav_accel;

  /* Get physical acceleration (hydro contribution) */
  a_phys_x += a_hydro[0] * cosmo->a_factor_hydro_accel;
  a_phys_y += a_hydro[1] * cosmo->a_factor_hydro_accel;
  a_phys_z += a_hydro[2] * cosmo->a_factor_hydro_accel;

  const float ac2 =
      a_phys_x * a_phys_x + a_phys_y * a_phys_y + a_phys_z * a_phys_z;

  const float ac_inv = (ac2 > 0.f) ? 1.f / sqrtf(ac2) : FLT_MAX;

  const float epsilon = gravity_get_softening(gp, grav_props);

  const float dt = sqrtf(2. * kernel_gravity_softening_plummer_equivalent_inv *
                         cosmo->a * grav_props->eta * epsilon * ac_inv);

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

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gp->potential_PM = 0.f;
  gp->a_grav_PM[0] = 0.f;
  gp->a_grav_PM[1] = 0.f;
  gp->a_grav_PM[2] = 0.f;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  gp->num_interacted = 0;
  gp->initialised = 1;
#endif
}

/**
 * @brief Finishes the gravity calculation.
 *
 * Multiplies the forces and accelerations by the appropiate constants.
 * Applies cosmological correction for periodic BCs.
 *
 * No need to apply the potential normalisation correction for periodic
 * BCs here since we do not compute the potential.
 *
 * @param gp The particle to act upon
 * @param const_G Newton's constant in internal units.
 * @param potential_normalisation Term to be added to all the particles.
 * @param periodic Are we using periodic BCs?
 */
__attribute__((always_inline)) INLINE static void gravity_end_force(
    struct gpart* gp, float const_G, const float potential_normalisation,
    const int periodic) {

  /* Let's get physical... */
  gp->a_grav[0] *= const_G;
  gp->a_grav[1] *= const_G;
  gp->a_grav[2] *= const_G;

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  gp->potential_PM *= const_G;
  gp->a_grav_PM[0] *= const_G;
  gp->a_grav_PM[1] *= const_G;
  gp->a_grav_PM[2] *= const_G;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  gp->initialised = 0; /* Ready for next step */
#endif
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
 * @param grav_props The global properties of the gravity calculation.
 */
__attribute__((always_inline)) INLINE static void gravity_first_init_gpart(
    struct gpart* gp, const struct gravity_props* grav_props) {

  gp->time_bin = 0;

  gravity_init_gpart(gp);
}

#endif /* SWIFT_DEFAULT_GRAVITY_H */
