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
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_density(
    const float r2, const float *dx, const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props, const integertime_t ti_current,
    const double time) {

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

  /* Neighbour sound speed */
  const float cj = hydro_get_comoving_soundspeed(pj);

  /* Contribution to the smoothed sound speed */
  bi->sound_speed_gas += mj * cj * wi;

  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);

  /* Contribution to the smoothed internal energy */
  bi->internal_energy_gas += mj * uj * wi;

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {pj->v[0] - bi->v[0], pj->v[1] - bi->v[1],
                       pj->v[2] - bi->v[2]};

  /* Contribution to the smoothed velocity (gas w.r.t. black hole) */
  bi->velocity_gas[0] += mj * dv[0] * wi;
  bi->velocity_gas[1] += mj * dv[1] * wi;
  bi->velocity_gas[2] += mj * dv[2] * wi;

  /* Contribution to the specific angular momentum of gas, which is later
   * converted to the circular velocity at the smoothing length */
  bi->spec_angular_momentum_gas[0] += mj * wi * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->spec_angular_momentum_gas[1] += mj * wi * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->spec_angular_momentum_gas[2] += mj * wi * (dx[0] * dv[1] - dx[1] * dv[0]);

  if (bh_props->multi_phase_bondi) {
    /* Contribution to BH accretion rate
     *
     * i) Calculate denominator in Bondi formula */
    const double gas_v_phys[3] = {dv[0] * cosmo->a_inv, dv[1] * cosmo->a_inv,
                                  dv[2] * cosmo->a_inv};
    const double gas_v_norm2 = gas_v_phys[0] * gas_v_phys[0] +
                               gas_v_phys[1] * gas_v_phys[1] +
                               gas_v_phys[2] * gas_v_phys[2];

    const double gas_c_phys = cj * cosmo->a_factor_sound_speed;
    const double gas_c_phys2 = gas_c_phys * gas_c_phys;
    const double denominator2 = gas_v_norm2 + gas_c_phys2;

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure that the denominator is strictly positive */
    if (denominator2 <= 0)
      error(
          "Invalid denominator for BH particle %lld and gas particle "
          "%lld in Bondi rate calculation.",
          bi->id, pj->id);
#endif
    const double denominator_inv = 1. / sqrt(denominator2);

    /* ii) Contribution of gas particle to the BH accretion rate
     *     (without constant pre-factor)
     *     N.B.: rhoj is the weighted contribution to BH gas density. */
    const float rhoj = mj * wi * cosmo->a3_inv;
    bi->accretion_rate +=
        rhoj * denominator_inv * denominator_inv * denominator_inv;

  } /* End of accretion contribution calculation */

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
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(const float r2, const float *dx,
                                  const float hi, const float hj,
                                  struct bpart *bi, struct part *pj,
                                  struct xpart *xpj, const int with_cosmology,
                                  const struct cosmology *cosmo,
                                  const struct gravity_props *grav_props,
                                  const struct black_holes_props *bh_props,
                                  const integertime_t ti_current,
                                  const double time) {

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
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Is this gas neighbour close enough that we can consider its potential
     for repositioning? */
  if (r2 < max_dist_repos2) {

    /* Flag to check whether neighbour is slow enough to be considered
     * as repositioning target. Always true if velocity cut is switched off. */
    int neighbour_is_slow_enough = 1;
    if (bh_props->with_reposition_velocity_threshold) {

      /* Compute relative peculiar velocity between the two BHs
       * Recall that in SWIFT v is (v_pec * a) */
      const float delta_v[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                                bi->v[2] - pj->v[2]};
      const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                       delta_v[2] * delta_v[2];
      const float v2_pec = v2 * cosmo->a2_inv;

      /* Compute the maximum allowed velocity */
      float v2_max = bh_props->max_reposition_velocity_ratio *
                     bh_props->max_reposition_velocity_ratio *
                     bi->sound_speed_gas * bi->sound_speed_gas;

      /* If desired, limit the value of the threshold (v2_max) to be no
       * smaller than a user-defined value */
      if (bh_props->min_reposition_velocity_threshold > 0) {
        const float v2_min_thresh =
            bh_props->min_reposition_velocity_threshold *
            bh_props->min_reposition_velocity_threshold;
        v2_max = max(v2_max, v2_min_thresh);
      }

      /* Is the neighbour too fast to jump to? */
      if (v2_pec >= v2_max) neighbour_is_slow_enough = 0;
    }

    if (neighbour_is_slow_enough) {
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

  /* Check if the BH needs to be fed. If not, we're done here */
  const float bh_mass_deficit = bi->subgrid_mass - bi->mass_at_start_of_step;
  if (bh_mass_deficit <= 0) return;

  if (bh_props->use_nibbling) {

    /* If we do nibbling, things are quite straightforward. We transfer
     * the mass and all associated quantities right here. */

    if (bh_props->epsilon_r == 1) return;

    const float bi_mass_orig = bi->mass;
    const float pj_mass_orig = hydro_get_mass(pj);

    /* Don't nibble from particles that are too small already */
    if (pj_mass_orig < bh_props->min_gas_mass_for_nibbling) return;

    /* Next line is equivalent to w_ij * m_j / Sum_j (w_ij * m_j) */
    const float particle_weight = hi_inv_dim * wi * pj_mass_orig / bi->rho_gas;
    float nibble_mass = bh_mass_deficit * particle_weight;

    /* We radiated away some of the accreted mass, so need to take slightly
     * more from the gas than the BH gained */
    const float excess_fraction = 1.0 / (1.0 - bh_props->epsilon_r);

    /* Need to check whether nibbling would push gas mass below minimum
     * allowed mass */
    float new_gas_mass = pj_mass_orig - nibble_mass * excess_fraction;
    if (new_gas_mass < bh_props->min_gas_mass_for_nibbling) {
      new_gas_mass = bh_props->min_gas_mass_for_nibbling;
      nibble_mass = (pj_mass_orig - bh_props->min_gas_mass_for_nibbling) /
                    excess_fraction;
    }

    /* Transfer (dynamical) mass from the gas particle to the BH */
    bi->mass += nibble_mass;
    hydro_set_mass(pj, new_gas_mass);

    /* Add the angular momentum of the accreted gas to the BH total.
     * Note no change to gas here. The cosmological conversion factors for
     * velocity (a^-1) and distance (a) cancel out, so the angular momentum
     * is already in physical units. */
    const float dv[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                         bi->v[2] - pj->v[2]};
    bi->swallowed_angular_momentum[0] +=
        nibble_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
    bi->swallowed_angular_momentum[1] +=
        nibble_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
    bi->swallowed_angular_momentum[2] +=
        nibble_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

    /* Update the BH momentum and velocity. Again, no change to gas here. */
    const float bi_mom[3] = {bi_mass_orig * bi->v[0] + nibble_mass * pj->v[0],
                             bi_mass_orig * bi->v[1] + nibble_mass * pj->v[1],
                             bi_mass_orig * bi->v[2] + nibble_mass * pj->v[2]};

    bi->v[0] = bi_mom[0] / bi->mass;
    bi->v[1] = bi_mom[1] / bi->mass;
    bi->v[2] = bi_mom[2] / bi->mass;
    bi->gpart->v_full[0] = bi->v[0];
    bi->gpart->v_full[1] = bi->v[1];
    bi->gpart->v_full[2] = bi->v[2];

    /* Update the BH and also gas metal masses */
    struct chemistry_bpart_data *bi_chem = &bi->chemistry_data;
    struct chemistry_part_data *pj_chem = &pj->chemistry_data;
    chemistry_transfer_part_to_bpart(
        bi_chem, pj_chem, nibble_mass * excess_fraction,
        nibble_mass * excess_fraction / pj_mass_orig);

  } else { /* ends nibbling section, below comes swallowing */

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
  } /* ends section for swallowing */
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
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_swallow(const float r2, const float *dx,
                                 const float hi, const float hj,
                                 struct bpart *bi, struct bpart *bj,
                                 const struct cosmology *cosmo,
                                 const struct gravity_props *grav_props,
                                 const struct black_holes_props *bh_props,
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
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Is this BH neighbour close enough that we can consider its potential
     for repositioning? */
  if (r2 < max_dist_repos2) {

    /* Flag to check whether neighbour is slow enough to be considered
     * as repositioning target. Always true if velocity cut switched off */
    int neighbour_is_slow_enough = 1;
    if (bh_props->with_reposition_velocity_threshold) {

      /* Compute the maximum allowed velocity */
      float v2_max = bh_props->max_reposition_velocity_ratio *
                     bh_props->max_reposition_velocity_ratio *
                     bi->sound_speed_gas * bi->sound_speed_gas;

      /* If desired, limit the value of the threshold (v2_max) to be no
       * smaller than a user-defined value */
      if (bh_props->min_reposition_velocity_threshold > 0) {
        const float v2_min_thresh =
            bh_props->min_reposition_velocity_threshold *
            bh_props->min_reposition_velocity_threshold;
        v2_max = max(v2_max, v2_min_thresh);
      }

      /* Is the neighbour too fast to jump to? */
      if (v2_pec >= v2_max) neighbour_is_slow_enough = 0;
    }

    if (neighbour_is_slow_enough) {
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

  /* (Square of) max swallowing distance allowed based on the softening */
  const float max_dist_merge2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_merging_distance_ratio *
      bh_props->max_merging_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

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

    /* Maximum velocity difference between BHs allowed to merge */
    float v2_threshold;

    if (bh_props->merger_threshold_type == 0) {

      /* 'Old-style' merger threshold using circular velocity at the
       * edge of the more massive BH's kernel */
      v2_threshold = G_Newton * M / (kernel_gamma * h);
    } else {

      /* Arguably better merger threshold using the escape velocity at
       * the distance of the lower-mass BH */
      const float r_12 = sqrt(r2);

      if ((bh_props->merger_threshold_type == 1) &&
          (r_12 < grav_props->epsilon_baryon_cur)) {

        /* If BHs are within softening range, take this into account */
        const float w_grav =
            kernel_grav_pot_eval(r_12 / grav_props->epsilon_baryon_cur);
        const float r_mod = w_grav / grav_props->epsilon_baryon_cur;
        v2_threshold = 2.f * G_Newton * M / (r_mod);

      } else {

        /* Standard formula if BH interactions are not softened */
        v2_threshold = 2.f * G_Newton * M / (r_12);
      }
    } /* Ends sections for different merger thresholds */

    if ((v2_pec < v2_threshold) && (r2 < max_dist_merge2)) {

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
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_feedback(const float r2, const float *dx,
                                   const float hi, const float hj,
                                   const struct bpart *bi, struct part *pj,
                                   struct xpart *xpj, const int with_cosmology,
                                   const struct cosmology *cosmo,
                                   const struct gravity_props *grav_props,
                                   const struct black_holes_props *bh_props,
                                   const integertime_t ti_current,
                                   const double time) {

  /* Get the heating probability */
  const float prob = bi->to_distribute.AGN_heating_probability;

  /* Are we doing some feedback? */
  if (prob > 0.f) {

    /* Draw a random number (Note mixing both IDs) */
    const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                            random_number_BH_feedback);

    /* Are we lucky? */
    if (rand < prob) {

      /* Compute new energy per unit mass of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = bi->to_distribute.AGN_delta_u;
      const double u_new = u_init + delta_u;

      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Store the feedback energy */
      const double delta_energy = delta_u * hydro_get_mass(pj);
      tracers_after_black_holes_feedback(xpj, with_cosmology, cosmo->a, time,
                                         delta_energy);

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
