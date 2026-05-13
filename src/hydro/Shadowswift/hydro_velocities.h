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
#include <grackle.h>

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
 */
__attribute__((always_inline)) INLINE static void
hydro_set_velocity_from_momentum(const float* restrict momentum,
                                 float inverse_mass, float rho,
                                 const struct hydro_props* hydro_properties,
                                 float* restrict /*return*/ velocity) {
  if (rho < hydro_properties->epsilon_rho) {
    /* Suppress velocity linearly near vacuum */
    const float fac = rho * hydro_properties->epsilon_rho_inv;
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
 * @brief Calculates gradient of cold steering.
 *
 * Note: Only called with SHADOWSWIFT_STEERING_COLD_FLOWS_GRADIENT
 *
 * @param soundspeed
 * @param vchar_dt the maximum cold steering corrective velocity
 * @param Mach1 minimum threshold for cold steering
 * @param Mach2 threshold for using maximum steering
 * @param Mach_part mach of current particle
 *
 * Returns point along linear gradient that scales strength of velocity steering
 * from soundspeed to vchar_dt linearly between Mach values
 */
__attribute__((always_inline)) INLINE static float
hydro_velocities_steering_gradient(
  const float soundspeed, const float vchar_dt,
  const float Mach1, const float Mach2, const float Mach_part) {

  if (Mach1 == Mach2) {
    error("Mach1 and Mach2 Equal");
  }

  /* Set steering linear gradient (y = mx + c, find m) */
  float steering_gradient = (vchar_dt - soundspeed) / (Mach2 - Mach1);

  /* Fix c */
  float steering_constant = soundspeed - (vchar_dt - soundspeed) /
                                          (Mach2 / Mach1 - 1);

  /* Set Vchar */
  return steering_gradient * Mach_part + steering_constant;

}

/**
 * @brief Set the velocity of a ShadowSWIFT particle, based on the values of its
 * primitive variables and the geometry of its voronoi cell.
 *
 * @param p The particle to act upon.
 * @param xp The extended particle data to act upon.
 * @param hydro_properties The hydro_props struct
 * @param dt The hydrodynamical time-step of the particle.
 */
__attribute__((always_inline)) INLINE static void hydro_velocities_set(
    struct part* restrict p, struct xpart* restrict xp,
    const struct hydro_props* hydro_properties, float dt) {
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
    v[0] = fluid_v[0];
    v[1] = fluid_v[1];
    v[2] = fluid_v[2];

#ifdef SHADOWSWIFT_STEER_MOTION

#ifdef CENTER_OF_MASS_DENSITY
    /* Using the centre of mass not the centroid for steering */
    float ds[3];
    ds[0] = p->geometry.center_of_mass[0];
    ds[1] = p->geometry.center_of_mass[1];
    ds[2] = p->geometry.center_of_mass[2];
#else

    float ds[3];
    ds[0] = p->geometry.centroid[0];
    ds[1] = p->geometry.centroid[1];
    ds[2] = p->geometry.centroid[2];
#endif

    const float d = sqrtf(ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]);
    const float soundspeed = hydro_get_comoving_soundspeed(p);

#ifdef SHADOWSWIFT_STEERING_FACEANGLE_FLOWS
    /* Add a correction to the velocity to move generator closer to centre of
     * cell based on geometric arguments made in Vogelsberger 2012 2.2.2 (i)
     * where we consider the maximum angle of the face to the generator and
     * apply corrections if the cell is irregular.
     *
     * If this is not defined then it defaults onto previous Springel 2010
     * description, with the possibility to further enable cold steering.
     *
     * This is however not recommended, see Vogelsberger 2012 2.2.2 and
     * Weinberger 2020 (Arepo public release paper) 5.1 for more justifications
     * why not. The argument of face angles is allegedly more robust, it is
     * recommended that this feature remain enabled by default. */

    /* Angle Steering Parameters (AREPO Defaults) */
    const float max_angle = p->geometry.max_face_angle;
    const float beta = 2.25f; // Controls how steep angle is before steering. Lower -> more cells steered. AREPO = 2.25 or 2.0
    const float f_shaping_speed = 0.5; // Steering aggressiveness factor. AREPO 0.5
    const float vchar_dt = dt > 0. ? d / dt : 0.;// Timestep based correction
    float vchar = soundspeed; // Determines cold steering corrective velocity. Default is soundspeed.

    /* Prefactor for steering, catch 1/d errors now before inf velocity */
    float fac = d > 0. ? f_shaping_speed / d : 0.;

    /* Additionally impose that we use the timestep based steering if it is less
     * aggressive than soundspeed steering, and if the cell is very deformed */
    if ((vchar_dt < vchar) || (max_angle > 2 * beta)) {
      vchar = vchar_dt;
    }

    if (max_angle > 0.75 * beta) {
      /* Enter bottom two options of Weinberger 2020 Eq. 39
       * and apply default factor f_shaping */

      /* Fancy settings for vchar */
#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS
      /* Enables much more aggressive cold steering if the cell
       * has high Mach (assumed cold) */
      if (soundspeed <= 0.f) {
        vchar = 0;
      } else {
        const float Mach_low = 5.f;
        const float Mach_part =
            sqrtf(fluid_v[0] * fluid_v[0] + fluid_v[1] * fluid_v[1] +
                  fluid_v[2] * fluid_v[2]) /
            soundspeed;
        if (Mach_part > Mach_low) {

#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS_GRADIENT  // Apply cold steering
          /* Controls for maximum values for steering gradient */
          const float Mach_high = 10.f;

          /* Apply a smoother steering gradient if between Mach low and high */
          if (p->rho < 12.f) {
            // Activate the donkey steering
            if (Mach_part < Mach_high) {
              vchar = hydro_velocities_steering_gradient(
                  soundspeed, vchar_dt, Mach_low, Mach_high, Mach_part);
            } else {
              /* Mach is higher than max, just set to max steering */
              vchar = vchar_dt;
            }
          }

#else
          /* Default cold steering */
          vchar = vchar_dt;
#endif
        }
      }
#endif

      if (max_angle < beta) {
        /* Enter second option of Weinberger 2020 Eq. 39
         * and apply correction
         * This can only be < 1, therefore if we do not enter this criteria
         * it automatically applies maximum correction */
        fac *= (max_angle - 0.75 * beta) / (0.25 * beta);

      }  // If not, then hits third criteria.
    } // End Angle Criteria

    /* Set characteristic speed into factor */
    fac *= vchar;

    /* Apply velocity corrections */
    v[0] += ds[0] * fac;
    v[1] += ds[1] * fac;
    v[2] += ds[2] * fac;

/* Enter else correspondning to SHADOWSWIFT_STEERING_FACEANGLE_FLOWS
 * being off, but SHADOWSWIFT_STEER_MOTION being on*/
#else

    /* Add a correction to the velocity to keep particle positions close enough
       to the centroid of their voronoi cell. */
    /* The correction term below is based on the one from Springel (2010). */
    const float R = get_radius_dimension_sphere(p->geometry.volume);
    const float eta =
        0.25f;  // Roundness criterion - lower = rounder. Default 0.25
    const float etaR = eta * R;
    const float xi = 1.0f;  // Corrective velocity strength - should prevent
                            // shells. Default 1.0
    /* We only apply the correction if the offset between centroid and position
       is too large, or if the distance to the nearest face is smaller than the
       distance to the centroid. */
    if (d > 0.9f * etaR || d > p->geometry.min_face_dist) {
      float fac = xi * soundspeed / d;
      float fac_dt = dt > 0. ? 0.5f * xi / dt : 0.;
      if (fac_dt < fac) {
        fac = fac_dt;
      } else {
#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS
        /* In very cold flows, the sound speed may be significantly slower than
         * the actual speed of the particles, rendering this scheme ineffective.
         * In this case, use a criterion based on the timestep instead
         *
         * NOTE: A huge variety of criteria have been attempted here, including
         * temperature, density, soundspeed, a mix, and some brute force ones.
         * It is highly recommended to use the default face angle steering,
         * and use this with mach 5 - 50 as settings. Less is too inclusive,
         *
         */

        if (2500.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
                                                   fluid_v[1] * fluid_v[1] +
                                                   fluid_v[2] * fluid_v[2]) {
          fac = fac_dt;
        }
#endif
      }
      if (d < 1.1f * etaR) {
        fac *= 5.0f * (d - 0.9f * etaR) / etaR;
      }
      v[0] += ds[0] * fac;
      v[1] += ds[1] * fac;
      v[2] += ds[2] * fac;
    }
#endif  // Face Angle Flows
#endif  // Steering Motion
  }

  if (p->rho < hydro_properties->epsilon_rho) {
    /* Linearly suppress particle velocity near vacuum */
    float fac = p->rho * hydro_properties->epsilon_rho_inv;
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
