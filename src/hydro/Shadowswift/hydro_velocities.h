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
 * @param vchar the characteristic velocity governming steering strength
 * @param soundspeed
 * @param vchar_dt the maximum cold steering corrective velocity
 * @param Mach1 minimum threshold for cold steering
 * @param Mach2 threshold for using maximum steering
 * @param Mach_part mach of current particle
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_steering_gradient(
  float* restrict vchar, float soundspeed, float vchar_dt, float Mach1,
  float Mach2, float Mach_part) {
  /* Set steering linear gradient (y = mx + c, find m) */
  float steering_gradient = (vchar_dt - soundspeed) / (Mach2 - Mach1);

  /* Fix c */
  float steering_constant = soundspeed - (vchar_dt - soundspeed) /
                                          (Mach2 / Mach1 - 1);

  /* Set Vchar */
  *vchar = steering_gradient * Mach_part + steering_constant;
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

    float ds[3];
    ds[0] = p->geometry.centroid[0];
    ds[1] = p->geometry.centroid[1];
    ds[2] = p->geometry.centroid[2];
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
    float beta = 2.25f;
    float f_shaping_speed = 0.5f;
    float vchar = soundspeed; // Determines cold steering aggresssiveness
    float vchar_dt = d / dt; // Timestep based correction

#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS
    /* Enables much more aggressive cold steering if the cell
     * has high Mach (assumed cold) */
    if (soundspeed <= 0.f) {
      vchar = 0;
    } else {
      const float Mach_low = 5.f;
      const float Mach_part = sqrtf(fluid_v[0] * fluid_v[0] +
                                fluid_v[1] * fluid_v[1] +
                                fluid_v[2] * fluid_v[2]) / soundspeed;
      if (Mach_part > Mach_low) {
        const float Mach_high = 15.f;
        // Apply cold steering
#ifdef SHADOWSWIFT_STEERING_COLD_FLOWS_GRADIENT
        /* Apply a smoother gradient to steering if between Mach low and high */
        if (Mach_part < Mach_high) {
          hydro_velocities_steering_gradient(&vchar, soundspeed, vchar_dt,
            Mach_low, Mach_high, Mach_part);
        } else {
          /* Mach is higher than max, just set to max steering */
          vchar = vchar_dt;
        }
#else
        /* Default cold steering */
          vchar = vchar_dt;
#endif
        }
    }
#endif

    /* Additionally impose that we use the timestep based steering if
     * it is less aggressive than soundspeed steering */
    if (vchar_dt < vchar) {
      vchar = vchar_dt;
    }

    if (max_angle > 0.75 * beta) {
      /* Enter bottom two options of Weinberger 2020 Eq. 39
       * and apply default factor f_shaping */
      float fac = f_shaping_speed * vchar / d;

      if (max_angle < beta) {
        /* Enter second option of Weinberger 2020 Eq. 39
         * and apply correction*/
        fac *= (max_angle - 0.75 * beta) / (0.25 * beta);

      } // Else, no further factors needed to apply third option

      /* Apply velocity corrections */
      v[0] += ds[0] * fac;
      v[1] += ds[1] * fac;
      v[2] += ds[2] * fac;
    }
/* Enter else correspondning strictly to SHADOWSWIFT_STEERING_FACEANGLE_FLOWS * being off, but
 * SHADOWSWIFT_STEER_MOTION being on*/
#else

    /* Add a correction to the velocity to keep particle positions close enough
       to the centroid of their voronoi cell. */
    /* The correction term below is based on the one from Springel (2010). */
    const float R = get_radius_dimension_sphere(p->geometry.volume);
    const float eta = 0.25f; // Roundness criterion - lower = rounder. Default 0.25
    const float etaR = eta * R;
    const float xi = 1.0f; // Corrective velocity strength - should prevent shells. Default 1.0
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
         * In this case, use a criterion based on the timestep instead */

        // /* Gather all variables for cold steering criteria */
        // const float *g = xp->a_grav;
        //
        // float Egrav = p->conserved.mass * sqrtf(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) *
        //               hydro_get_comoving_psize(p);
        //
        // const float thermal_energy = xp->u_full * p->conserved.mass;
        //
        // /* NOT SAFE, CAN have 0 soundspeed */
        // const float mach_cold = sqrtf(fluid_v[0] * fluid_v[0] +
        //                                  fluid_v[1] * fluid_v[1] +
        //                                  fluid_v[2] * fluid_v[2]) / soundspeed;
        //
        // const float Ekin_cold = 0.5f * p->conserved.mass * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        // /* Apply cold steering if needed */
        // /* Demand cold and high mach, if very high mach always */
        // if (((thermal_energy < 1e-1 * Egrav) &&
        //           (mach_cold > 10)) || mach_cold > 50) {
        //   // message("here! Thermal = %g, Ekin + thermal = %g, Egrav = %g, mach_cold = %g, dens = %g",
        //   //   thermal_energy, p->timestepvars.Ekin + thermal_energy, Egrav, mach_cold, p->rho);
        //
        //   float div_v = p->gradients.v[0][0] +
        //                p->gradients.v[1][1] +
        //                p->gradients.v[2][2];
        //
        //   fac = fac_dt;
        // }

        /* If cold (and high M) -> No steering
         * If gravitationally dominant (and high M) -> No steering
         * Else -> Steering (for hotter less gravitationally excited gas) */

        // /* Steer if gas satisfies:
        //  * Cold
        //  * Gravitationally indifferent (not largely perturbed by g)
        //  * Mach is above 5.
        //  * Else: Do Not Cold Steer */
        // if ((thermal_energy < 1e-2 * (Ekin_cold + thermal_energy) &&
        //      thermal_energy > 1e-2 * Egrav) && mach_cold > 5) {
        //
        //   message("here! Thermal = %g, Ekinmax + thermal = %g, Ekin_cold = %g, Egrav = %g, mach_cold = %g, dens = %g, mass = %g, massflux = %g",
        //     thermal_energy, p->timestepvars.Ekin + thermal_energy, Ekin_cold, Egrav, mach_cold, p->rho, p->conserved.mass, p->flux.mass);
        //
        //
        //
        //   fac = fac_dt;
        // } // Not so sure this works... Egrav term still isnt enough. Plenty of cold unperturbed gas at centres


        /* Physical constants
         * Can't access properly but if we assume neutral minimum internal
         * energy initialisation, this works to get m_h * (gamma-1) / kB
         */
        const double mp_gammaminusone_over_kb =
          hydro_properties->minimal_temperature /
            (hydro_properties->minimal_internal_energy *
              hydro_properties->mu_neutral);



        /* Gas properties */
        const double T_transition = hydro_properties->hydrogen_ionization_temperature;
        const double mu_neutral = hydro_properties->mu_neutral;
        const double mu_ionised = hydro_properties->mu_ionised;

        /* Particle temperature */
        const double u = xp->u_full;

        /* Temperature over mean molecular weight */
        const double T_over_mu = u * mp_gammaminusone_over_kb;
        double temperature;

        /* Are we above or below the HII -> HI transition? */
        if (T_over_mu > (T_transition + 1.) / mu_ionised) {
          temperature = T_over_mu * mu_ionised;
        }
        else if (T_over_mu < (T_transition - 1.) / mu_neutral) {
          temperature = T_over_mu * mu_neutral;
        }
        else {
          temperature = T_transition;
        }

        // /* Donkey steering */
        // if (((p->rho < 5) &&
        //   (2500.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
        //                                       fluid_v[1] * fluid_v[1] +
        //                                       fluid_v[2] * fluid_v[2])) ||
        //   ((temperature > 10 * hydro_properties->hydrogen_ionization_temperature)
        //     && (25.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
        //                                           fluid_v[1] * fluid_v[1] +
        //                                           fluid_v[2] * fluid_v[2]))) {
        //
        //   // message("here! Thermal = %g, Ekinmax + thermal = %g, Ekin_cold = %g, Egrav = %g, mach_cold = %g, dens = %g, mass = %g, massflux = %g",
        //   //   thermal_energy, p->timestepvars.Ekin + thermal_energy, Ekin_cold, Egrav, mach_cold, p->rho, p->conserved.mass, p->flux.mass);
        //
        //
        //
        //   fac = fac_dt;
        //      }

        // /* Donkey steering, no density criteria needed */
        // if (((2500.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
        //                                       fluid_v[1] * fluid_v[1] +
        //                                       fluid_v[2] * fluid_v[2])) ||
        //   ((temperature > 10 * hydro_properties->hydrogen_ionization_temperature)
        //     && (25.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
        //                                           fluid_v[1] * fluid_v[1] +
        //                                           fluid_v[2] * fluid_v[2]))) {
        //
        //   fac = fac_dt;
        //                                           }



        // /* If mflux/m > 1e-2 (and low M) -> No Steering
        //  * If mflux/m < 1e-2 (and high M) -> Steering
        //  * Basically, little growth + high mach -> steer. The rest is fine */
        // if ((( (p->flux.mass / p->conserved.mass) < 1e-3) && mach_cold > 50) ) {
        //
        //   // message("here! Thermal = %g, Ekin_cold = %g, Egrav = %g, mach_cold = %g, dens = %g, mass = %g, massflux = %g, potential E = %g",
        //   //   thermal_energy, Ekin_cold, Egrav, mach_cold, p->rho, p->conserved.mass, p->flux.mass, p->gpart->potential * p->conserved.mass);
        //
        //   // if (thermal_energy > 1e-5 * (p->timestepvars.Ekin + thermal_energy)) {
        //   //   message("Ekin_max = %g, p.mach = %g,", p->timestepvars.Ekin, p->timestepvars.mach_number);
        //   // }
        //
        //   fac = fac_dt;
        //
        // } /* SEEMINGLY LOVES TO APPLY TO EXACT WRONG PARTICLES!, actually in tests it seems fine? IDK why it does that. In restarts from Mach50 it was good? Redo. Mach50 AND
        // bad idea sadly SADLY THIS LEADS TO NaNs and CRASHES WITH EXPLODING SHELLS!*/






        // /* Apply cold steering if needed */
        // /* Demand cold and high mach, if very high mach always */
        // if (((thermal_energy < 1e-3 * (p->timestepvars.Ekin + thermal_energy) ||
        //           thermal_energy < 1e-3 * Egrav) &&
        //           mach_cold < 8) || mach_cold > 50) {
        //   // message("here! Thermal = %g, Ekin + thermal = %g, Egrav = %g, mach_cold = %g",
        //   //   thermal_energy, p->timestepvars.Ekin + thermal_energy, Egrav, mach_cold);
        //   fac = fac_dt;
        //
        //           }

        // /* Gravitational Potential Based Steering:
        //  * If negative Potential -> Assume collapsing, No steering
        //  * If positive and cold (High M) -> Assume void and steer
        //  * Unsure how this will work for non-cosmo runs whose
        //  * potentials are always negative........... */
        // if ((100.f * soundspeed * soundspeed < fluid_v[0] * fluid_v[0] +
        //                                          fluid_v[1] * fluid_v[1] +
        //                                          fluid_v[2] * fluid_v[2])
        //                                          && p->gpart->potential > 0.f) {
        //
        //   fac = fac_dt;
        // }





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
#endif // Face Angle Flows
#endif // Steering Motion

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
