/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_BH_IACT_H
#define SWIFT_EAGLE_BH_IACT_H

/* Local includes */
#include "black_holes_parameters.h"
#include "gravity.h"
#include "hydro.h"
#include "random.h"
#include "space.h"
#include "timestep_sync_part.h"
#include "tracers.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas, not updated).
 * @param xpj The extended data of the second particle (not updated).
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_density(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct bpart *bi, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct gravity_props *grav_props,
                                  const integertime_t ti_current) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the number of neighbours */
  bi->num_ngbs += 1;

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Contribution to the BH gas density */
  bi->rho_gas += mj * wi;

  /* Contribution to the total neighbour mass */
  bi->ngb_mass += mj;

  /* Neighbour sounds speed */
  const float cj = hydro_get_comoving_soundspeed(pj);

  /* Contribution to the smoothed sound speed */
  bi->sound_speed_gas += mj * cj * wi;

  /* Neighbour peculiar drifted velocity */
  const float vj[3] = {pj->v[0], pj->v[1], pj->v[2]};

  /* Contribution to the smoothed velocity */
  bi->velocity_gas[0] += mj * vj[0] * wi;
  bi->velocity_gas[1] += mj * vj[1] * wi;
  bi->velocity_gas[2] += mj * vj[2] * wi;

  /* Contribution to the circular valocity */
  const float dv[3] = {bi->v[0] - vj[0], bi->v[1] - vj[1], bi->v[2] - vj[2]};
  bi->circular_velocity_gas[0] += mj * wi * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->circular_velocity_gas[1] += mj * wi * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->circular_velocity_gas[2] += mj * wi * (dx[0] * dv[1] - dx[1] * dv[0]);

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Swallowing interaction between two particles (non-symmetric).
 *
 * Function used to flag the gas particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct bpart *bi, struct part *pj,
                                  struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const struct gravity_props *grav_props,
                                  const integertime_t ti_current) {

  float wi;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Start by checking the repositioning criteria */

  /* (Square of) Max repositioning distance allowed based on the softening */
  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      const_max_repositioning_distance_ratio *
      const_max_repositioning_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* This gas neighbour is close enough that we can consider it's potential
     for repositioning */
  if (r2 < max_dist_repos2) {

    /* Compute relative peculiar velocity between the two BHs
     * Recall that in SWIFT v is (v_pec * a) */
    const float delta_v[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                              bi->v[2] - pj->v[2]};
    const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                     delta_v[2] * delta_v[2];

    const float v2_pec = v2 * cosmo->a2_inv;

    /* Check the velocity criterion */
    if (v2_pec < 0.25f * bi->sound_speed_gas * bi->sound_speed_gas) {

      const float potential = pj->black_holes_data.potential;

      /* Is the potential lower? */
      if (potential < bi->reposition.min_potential) {

        /* Store this as our new best */
        bi->reposition.min_potential = potential;
        bi->reposition.delta_x[0] = -dx[0];
        bi->reposition.delta_x[1] = -dx[1];
        bi->reposition.delta_x[2] = -dx[2];
      }
    }
  }

  /* Is the BH hungry? */
  if (bi->subgrid_mass > bi->mass) {

    /* Probability to swallow this particle
     * Recall that in SWIFT the SPH kernel is recovered by computing
     * kernel_eval() and muliplying by (1/h^d) */
    const float prob =
        (bi->subgrid_mass - bi->mass) * hi_inv_dim * wi / bi->rho_gas;

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                            random_number_BH_swallow);

    /* Are we lucky? */
    if (rand < prob) {

      /* This particle is swallowed by the BH with the largest ID of all the
       * candidates wanting to swallow it */
      if (pj->black_holes_data.swallow_id < bi->id) {

        message("BH %lld wants to swallow gas particle %lld", bi->id, pj->id);

        pj->black_holes_data.swallow_id = bi->id;

      } else {

        message(
            "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
            "swallow id=%lld)",
            bi->id, pj->id, pj->black_holes_data.swallow_id);
      }
    }
  }
}

/**
 * @brief Swallowing interaction between two BH particles (non-symmetric).
 *
 * Function used to flag the BH particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param bj Second particle (black hole)
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_swallow(const float r2, const float *dx,
                                 const float hi, const float hj,
                                 struct bpart *bi, struct bpart *bj,
                                 const struct cosmology *cosmo,
                                 const struct gravity_props *grav_props,
                                 const integertime_t ti_current) {

  /* Compute relative peculiar velocity between the two BHs
   * Recall that in SWIFT v is (v_pec * a) */
  const float delta_v[3] = {bi->v[0] - bj->v[0], bi->v[1] - bj->v[1],
                            bi->v[2] - bj->v[2]};
  const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                   delta_v[2] * delta_v[2];

  const float v2_pec = v2 * cosmo->a2_inv;

  /* (Square of) Max repositioning distance allowed based on the softening */
  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      const_max_repositioning_distance_ratio *
      const_max_repositioning_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* This gas neighbour is close enough that we can consider it's potential
     for repositioning */
  if (r2 < max_dist_repos2) {

    /* Check the velocity criterion */
    if (v2_pec < 0.25f * bi->sound_speed_gas * bi->sound_speed_gas) {

      const float potential = bj->reposition.potential;

      /* Is the potential lower? */
      if (potential < bi->reposition.min_potential) {

        /* Store this as our new best */
        bi->reposition.min_potential = potential;
        bi->reposition.delta_x[0] = -dx[0];
        bi->reposition.delta_x[1] = -dx[1];
        bi->reposition.delta_x[2] = -dx[2];
      }
    }
  }

  /* Find the most massive of the two BHs */
  float M = bi->subgrid_mass;
  float h = hi;
  if (bj->subgrid_mass > M) {
    M = bj->subgrid_mass;
    h = hj;
  }

  /* (Square of) Max swallowing distance allowed based on the softening */
  const float max_dist_merge2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      const_max_merging_distance_ratio * const_max_merging_distance_ratio *
      grav_props->epsilon_baryon_cur * grav_props->epsilon_baryon_cur;

  const float G_Newton = grav_props->G_Newton;

  /* The BH with the smaller mass will be merged onto the one with the
   * larger mass.
   * To avoid rounding issues, we additionally check for IDs if the BHs
   * have the exact same mass. */
  if ((bj->subgrid_mass < bi->subgrid_mass) ||
      (bj->subgrid_mass == bi->subgrid_mass && bj->id < bi->id)) {

    /* Merge if gravitationally bound AND if within max distance
     * Note that we use the kernel support here as the size and not just the
     * smoothing length */
    if (v2_pec < G_Newton * M / (kernel_gamma * h) && (r2 < max_dist_merge2)) {

      /* This particle is swallowed by the BH with the largest ID of all the
       * candidates wanting to swallow it */
      if ((bj->merger_data.swallow_mass < bi->subgrid_mass) ||
          (bj->merger_data.swallow_mass == bi->subgrid_mass &&
           bj->merger_data.swallow_id < bi->id)) {

        message("BH %lld wants to swallow BH particle %lld", bi->id, bj->id);

        bj->merger_data.swallow_id = bi->id;
        bj->merger_data.swallow_mass = bi->subgrid_mass;

      } else {

        message(
            "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
            "swallow id=%lld)",
            bi->id, bj->id, bj->merger_data.swallow_id);
      }
    }
  }
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_feedback(const float r2, const float *dx,
                                   const float hi, const float hj,
                                   const struct bpart *bi, struct part *pj,
                                   struct xpart *xpj,
                                   const struct cosmology *cosmo,
                                   const struct gravity_props *grav_props,
                                   const integertime_t ti_current) {

  /* Get the heating probability */
  const float prob = bi->to_distribute.AGN_heating_probability;

  /* Are we doing some feedback? */
  if (prob > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                            random_number_BH_feedback);

    /* Are we lucky? */
    if (rand < prob) {

      /* Compute new energy of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = bi->to_distribute.AGN_delta_u;
      const double u_new = u_init + delta_u;

      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Mark this particle has having been heated by AGN feedback */
      tracers_after_black_holes_feedback(xpj);

      /* message( */
      /*     "We did some AGN heating! id %llu BH id %llu probability " */
      /*     " %.5e  random_num %.5e du %.5e du/ini %.5e", */
      /*     pj->id, bi->id, prob, rand, delta_u, delta_u / u_init); */

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_force < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_force[si->num_ngb_force] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_force;
#endif
}

#endif /* SWIFT_EAGLE_BH_IACT_H */
