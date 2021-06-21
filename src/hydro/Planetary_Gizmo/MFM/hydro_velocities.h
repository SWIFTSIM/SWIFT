/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_PLANETARY_GIZMO_MFM_HYDRO_VELOCITIES_H
#define SWIFT_PLANETARY_GIZMO_MFM_HYDRO_VELOCITIES_H

#ifdef GIZMO_FIX_PARTICLES
#error "Fixed particles are not allowed for GIZMO MFM!"
#endif

#ifdef GIZMO_STEER_MOTION
#error "Steering particle movement is not allowed for GIZMO MFM!"
#endif

/**
 * @brief Initialize the GIZMO particle velocities before the start of the
 * actual run based on the initial value of the primitive velocity.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_init(
    struct part* restrict p, struct xpart* restrict xp) {

  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
}

/**
 * @brief Set the particle velocity field that will be used to deboost fluid
 * velocities during the force loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_prepare_force(struct part* restrict p,
                               const struct xpart* restrict xp) {}

/**
 * @brief Set the variables that will be used to update the smoothing length
 * during the drift (these will depend on the movement of the particles).
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_end_force(
    struct part* restrict p) {

  /* Add normalization to h_dt. */
  p->force.h_dt *= p->h * hydro_dimension_inv;
}

/**
 * @brief Set the velocity of a GIZMO particle, based on the values of its
 * primitive variables and the geometry of its mesh-free "cell".
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_set(
    struct part* restrict p, struct xpart* restrict xp) {

  /* Set the velocities: */
  /* We first set the particle velocity */
  if (p->conserved.mass > 0.0f && p->rho > 0.0f) {

    const float inverse_mass = 1.0f / p->conserved.mass;

    /* Normal case: set particle velocity to fluid velocity. */
    xp->v_full[0] = p->conserved.momentum[0] * inverse_mass;
    xp->v_full[1] = p->conserved.momentum[1] * inverse_mass;
    xp->v_full[2] = p->conserved.momentum[2] * inverse_mass;

  } else {
    /* Vacuum particles have no fluid velocity. */
    xp->v_full[0] = 0.0f;
    xp->v_full[1] = 0.0f;
    xp->v_full[2] = 0.0f;
  }

  if (p->gpart) {
    p->gpart->v_full[0] = xp->v_full[0];
    p->gpart->v_full[1] = xp->v_full[1];
    p->gpart->v_full[2] = xp->v_full[2];
  }
}

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_reset_centroids(struct part* restrict p) {}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param p Particle.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_normalise_centroid(struct part* restrict p,
                                    const float wcount) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param p Particle (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_update_centroid_left(struct part* restrict p, const float* dx,
                                      const float w) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param p Particle (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_update_centroid_right(struct part* restrict p, const float* dx,
                                       const float w) {}

#endif /* SWIFT_PLANETARY_GIZMO_MFM_HYDRO_VELOCITIES_H */
