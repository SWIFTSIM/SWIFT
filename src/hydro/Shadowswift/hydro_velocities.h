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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_VELOCITIES_H
#define SWIFT_SHADOWSWIFT_HYDRO_VELOCITIES_H

/**
 * @brief Initialize the ShadowSWIFT particle velocities before the start of the
 * actual run based on the initial value of the primitive velocity.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_init(
    struct part* restrict p, struct xpart* restrict xp) {

#ifdef SHADOWSWIFT_FIX_PARTICLES
  xp->v_full[0] = 0.0f;
  xp->v_full[1] = 0.0f;
  xp->v_full[2] = 0.0f;
#else
  xp->v_full[0] = p->v[0];
  xp->v_full[1] = p->v[1];
  xp->v_full[2] = p->v[2];
#endif
}

/**
 * @brief Set the velocities based on the particles momentum.
 *
 * Velocities near vacuum are linearly suppressed.
 */
__attribute__((always_inline)) INLINE static void
hydro_set_velocity_from_momentum(const float* restrict momentum,
                                 float inverse_mass, float rho,
                                 float* restrict /*return*/ velocity) {
  if (rho < 0) {
    /* Suppress velocity linearly near vacuum */
    const float fac = rho * 1e20f;
    velocity[0] = fac * momentum[0] * inverse_mass;
    velocity[1] = fac * momentum[1] * inverse_mass;
    velocity[2] = fac * momentum[2] * inverse_mass;
  } else {
    /* Normal case: update fluid velocity and set particle velocity
     * accordingly. */
    velocity[0] = momentum[0] * inverse_mass;
    velocity[1] = momentum[1] * inverse_mass;
    velocity[2] = momentum[2] * inverse_mass;
  }
}

/**
 * @brief Applies a hydrodynamical half kick to the genererator velocity.
 *
 * Note: The fluid velocity is updated in the drift, we apply a hydro kick to
 * the generator velocity.
 * The gravity half kick is applied elsewhere.
 *
 * @param p The #part to kick
 * @param xp THe #xpart
 * @param dt The kick timestep for thermodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void
hydro_generator_velocity_half_kick(struct part* p, struct xpart* xp, float dt) {
  if (p->rho > 0.) {
    float rho_inv = 1.f / p->rho;
    xp->v_full[0] -= dt * p->gradients.P[0] * rho_inv;
    xp->v_full[1] -= dt * p->gradients.P[1] * rho_inv;
    xp->v_full[2] -= dt * p->gradients.P[2] * rho_inv;
  }
}

/**
 * @brief Set the velocity of a ShadowSWIFT particle, based on the values of its
 * primitive variables and the geometry of its voronoi cell.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param dt The hydrodynamical time-step of the particle.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_set(
    struct part* restrict p, struct xpart* restrict xp, float dt) {

  /* We first get the particle velocity. */
  float v[3] = {0.f, 0.f, 0.f};
  float fluid_v[3] = {0.f, 0.f, 0.f};
  int fix_particle = 0;

#ifdef SHADOWSWIFT_FIX_PARTICLES
  fix_particle = 1;
#elif defined(SWIFT_BOUNDARY_PARTICLES)
  if (p->id < SWIFT_BOUNDARY_PARTICLES) {
    fix_particle = 1;
  }
#endif

  if (!fix_particle && p->conserved.mass > 0.0f && p->rho > 0.0f) {

    /* Normal case: use (kicked) fluid velocity.
     * NOTE (yuytenh, 2025): At this point, xp->v_full will have recieved all
     * the necessary hydro and gravity half kicks, while p->v and p->v_part_full
     * will still be at the value from the last full timestep. So we can
     * calculate the kicked fluid velocity as follows: */
    fluid_v[0] = p->v[0] + xp->v_full[0] - p->v_part_full[0];
    fluid_v[1] = p->v[1] + xp->v_full[1] - p->v_part_full[1];
    fluid_v[2] = p->v[2] + xp->v_full[2] - p->v_part_full[2];

#ifdef SHADOWSWIFT_STEER_MOTION
    /* Add a correction to the velocity to keep particle positions close enough
       to the centroid of their voronoi cell. */
    /* The correction term below is based on the one from Springel (2010). */
    float ds[3];
    ds[0] = p->geometry.centroid[0];
    ds[1] = p->geometry.centroid[1];
    ds[2] = p->geometry.centroid[2];
    const float d = sqrtf(ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]);
    const float R = get_radius_dimension_sphere(p->geometry.volume);
    const float eta = 0.25f;
    const float etaR = eta * R;
    const float xi = 1.0f;
    const float soundspeed = hydro_get_comoving_soundspeed(p);
    /* We only apply the correction if the offset between centroid and position
       is too large, or if the distance to the nearest face is smaller than the
       distance to the centroid. */
    if (d > 0.9f * etaR || d > p->geometry.min_face_dist) {
      float fac = xi * soundspeed / d;
      float fac_dt = 0.5f * xi / dt;
      if (fac_dt < fac) {
        fac = fac_dt;
      } else {
#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS
        /* In very cold flows, the sound speed may be significantly slower than
         * the actual speed of the particles, rendering this scheme ineffective.
         * In this case, use a criterion based on the timestep instead */
        if (25.f * soundspeed * soundspeed <
            p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]) {
          fac = fac_dt;
        }
#endif
      }
      if (d < 1.1f * etaR) {
        fac *= 5.0f * (d - 0.9f * etaR) / etaR;
      }
      v[0] = fluid_v[0] + ds[0] * fac;
      v[1] = fluid_v[1] + ds[1] * fac;
      v[2] = fluid_v[2] + ds[2] * fac;
    }
#endif  // SHADOWSWIFT_STEER_MOTION
  }

  if (p->rho < 0.) {
    float fac = p->rho * 1e-10;
    v[0] *= fac;
    v[1] *= fac;
    v[2] *= fac;
  }

  /* Now make sure all velocity variables are up-to-date with the new generator
   * velocity. */
  xp->v_full[0] = v[0];
  xp->v_full[1] = v[1];
  xp->v_full[2] = v[2];
  p->v_part_full[0] = v[0];
  p->v_part_full[1] = v[1];
  p->v_part_full[2] = v[2];
  p->v_rel_full[0] = fluid_v[0] - p->v_part_full[0];
  p->v_rel_full[1] = fluid_v[1] - p->v_part_full[1];
  p->v_rel_full[2] = fluid_v[2] - p->v_part_full[2];

  if (p->gpart) {
    p->gpart->v_full[0] = v[0];
    p->gpart->v_full[1] = v[1];
    p->gpart->v_full[2] = v[2];
  }
}

#endif /* SWIFT_SHADOWSWIFT_HYDRO_VELOCITIES_H */
