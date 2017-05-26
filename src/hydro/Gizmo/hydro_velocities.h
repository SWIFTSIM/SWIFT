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
#ifndef SWIFT_HYDRO_VELOCITIES_H
#define SWIFT_HYDRO_VELOCITIES_H

/**
 * @brief Initialize the GIZMO particle velocities before the start of the
 * actual run based on the initial value of the primitive velocity.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_init(
    struct part* restrict p, struct xpart* restrict xp) {

#ifdef GIZMO_FIX_PARTICLES
  p->v[0] = 0.;
  p->v[1] = 0.;
  p->v[2] = 0.;
#else
  p->v[0] = p->primitives.v[0];
  p->v[1] = p->primitives.v[1];
  p->v[2] = p->primitives.v[2];
#endif

  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
}

/**
 * @brief Set the particle velocity field that will be used to deboost fluid
 * velocities during the force loop.
 *
 * @param p The particle to act upon.
 * @param xp The extended particel data to act upon.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_prepare_force(struct part* restrict p,
                               const struct xpart* restrict xp) {

#ifndef GIZMO_FIX_PARTICLES
  p->force.v_full[0] = xp->v_full[0];
  p->force.v_full[1] = xp->v_full[1];
  p->force.v_full[2] = xp->v_full[2];
#endif
}

/**
 * @brief Set the variables that will be used to update the smoothing length
 * during the drift (these will depend on the movement of the particles).
 *
 * @param p The particle to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_end_force(
    struct part* restrict p) {

#ifdef GIZMO_FIX_PARTICLES
  /* disable the smoothing length update, since the smoothing lengths should
     stay the same for all steps (particles don't move) */
  p->force.h_dt = 0.0f;
#else
  /* Add normalization to h_dt. */
  p->force.h_dt *= p->h * hydro_dimension_inv;
#endif
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

/* We first set the particle velocity. */
#ifdef GIZMO_FIX_PARTICLES

  p->v[0] = 0.;
  p->v[1] = 0.;
  p->v[2] = 0.;

#else  // GIZMO_FIX_PARTICLES

  if (p->conserved.mass > 0. && p->primitives.rho > 0.) {
    /* Normal case: set particle velocity to fluid velocity. */
    p->v[0] = p->conserved.momentum[0] / p->conserved.mass;
    p->v[1] = p->conserved.momentum[1] / p->conserved.mass;
    p->v[2] = p->conserved.momentum[2] / p->conserved.mass;

#ifdef GIZMO_STEER_MOTION

    /* Add a correction to the velocity to keep particle positions close enough
       to
       the centroid of their mesh-free "cell". */
    /* The correction term below is the same one described in Springel (2010).
     */
    float ds[3];
    ds[0] = p->geometry.centroid[0];
    ds[1] = p->geometry.centroid[1];
    ds[2] = p->geometry.centroid[2];
    const float d = sqrtf(ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]);
    const float R = get_radius_dimension_sphere(p->geometry.volume);
    const float eta = 0.25;
    const float etaR = eta * R;
    const float xi = 1.;
    const float soundspeed =
        sqrtf(hydro_gamma * p->primitives.P / p->primitives.rho);
    /* We only apply the correction if the offset between centroid and position
       is
       too large. */
    if (d > 0.9 * etaR) {
      float fac = xi * soundspeed / d;
      if (d < 1.1 * etaR) {
        fac *= 5. * (d - 0.9 * etaR) / etaR;
      }
      p->v[0] -= ds[0] * fac;
      p->v[1] -= ds[1] * fac;
      p->v[2] -= ds[2] * fac;
    }

#endif  // GIZMO_STEER_MOTION
  } else {
    /* Vacuum particles have no fluid velocity. */
    p->v[0] = 0.;
    p->v[1] = 0.;
    p->v[2] = 0.;
  }

#endif  // GIZMO_FIX_PARTICLES

  /* Now make sure all velocity variables are up to date. */
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];

  if (p->gpart) {
    p->gpart->v_full[0] = p->v[0];
    p->gpart->v_full[1] = p->v[1];
    p->gpart->v_full[2] = p->v[2];
  }
}

#endif  // SWIFT_HYDRO_VELOCITIES_H
