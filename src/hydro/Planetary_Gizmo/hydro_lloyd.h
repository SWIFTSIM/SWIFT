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
#ifndef SWIFT_PLANETARY_GIZMO_HYDRO_LLOYD_H
#define SWIFT_PLANETARY_GIZMO_HYDRO_LLOYD_H

#ifdef GIZMO_LLOYD_ITERATION

#if !defined(GIZMO_MFV_SPH)
#error "LLoyd's algorithm only works if you use GIZMO MFV!"
#endif

/*! @brief Set the time step for each particle to the CFL condition. */
#define hydro_gizmo_lloyd_skip_timestep(CFL_condition) return CFL_condition;

/*! @brief Skip the drift for Lloyd's algorithm. */
#define hydro_gizmo_lloyd_skip_drift() return;

/**
 * @brief Set the initial values for the primitive and conserved variables and
 * the particle velocity to safe values for use with Lloyd's algorithm.
 *
 * @param W Primitive variable array (size 5 or more).
 * @param Q Conserved variable array (size 5 or more).
 * @param v Particle velocity array (size 3 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_initialize_particle(float *W, float *Q, float *v) {

  W[0] = 1.0f;
  W[1] = 0.0f;
  W[2] = 0.0f;
  W[3] = 0.0f;
  W[4] = 1.0f;

  Q[0] = 1.0f;
  Q[1] = 0.0f;
  Q[2] = 0.0f;
  Q[3] = 0.0f;
  Q[4] = 1.0f;

  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = 0.0f;
}

/**
 * @brief Reset the primitive variables after they have been updated to make
 * sure they still have safe values for use with Lloyd's algorithm.
 *
 * @param W Primitive variable array (size 5 or more).
 */
__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_reset_primitive_variables(float *W) {

  W[0] = 1.0f;
  W[1] = 0.0f;
  W[2] = 0.0f;
  W[3] = 0.0f;
  W[4] = 1.0f;
}

/**
 * @brief Reset the gradients to 0 after they have been computed, since we don't
 * use them for Lloyd's algorithm.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_reset_gradients(struct part *restrict p) {
  hydro_gradients_init(p);
}

/**
 * @brief Reset the force variables to safe values for use with Lloyd's
 * algorithm.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gizmo_lloyd_end_force(
    struct part *restrict p) {

  p->force.h_dt = 0.0f;
}

/**
 * @brief Reset the conserved variables after the flux exchange and "kick" the
 * particle towards its centroid using the given kick time step.
 *
 * @param p Particle.
 * @param xp Extended particle data.
 * @param dt Kick time step.
 */
__attribute__((always_inline)) INLINE static void hydro_gizmo_lloyd_kick(
    struct part *restrict p, struct xpart *restrict xp, const float dt) {

  /* reset conserved variables to safe values */
  p->conserved.mass = 1.;
  p->conserved.momentum[0] = 0.;
  p->conserved.momentum[1] = 0.;
  p->conserved.momentum[2] = 0.;
  p->conserved.energy = 1.;

  /* set the particle velocities to the Lloyd velocities */
  /* note that centroid is the relative position of the centroid w.r.t. the
     particle position (position - centroid) */
  xp->v_full[0] = -p->geometry.centroid[0] / dt;
  xp->v_full[1] = -p->geometry.centroid[1] / dt;
  xp->v_full[2] = -p->geometry.centroid[2] / dt;
  p->v[0] = xp->v_full[0];
  p->v[1] = xp->v_full[1];
  p->v[2] = xp->v_full[2];
}

#else /* no GIZMO_LLOYD_ITERATION */

#define hydro_gizmo_lloyd_skip_timestep(CFL_condition)
#define hydro_gizmo_lloyd_skip_drift()

__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_initialize_particle(float *W, float *Q, float *v) {}

__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_reset_primitive_variables(float *W) {}

__attribute__((always_inline)) INLINE static void
hydro_gizmo_lloyd_reset_gradients(struct part *restrict p) {}

__attribute__((always_inline)) INLINE static void hydro_gizmo_lloyd_end_force(
    struct part *restrict p) {}

__attribute__((always_inline)) INLINE static void hydro_gizmo_lloyd_kick(
    struct part *restrict p, struct xpart *restrict xp, const float dt) {}

#endif /* GIZMO_LLOYD_ITERATION */

#endif /* SWIFT_PLANETARY_GIZMO_HYDRO_LLOYD_H */
