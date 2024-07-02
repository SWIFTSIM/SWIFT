/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 * Copyright (c) 2022 Filip Husko (filip.husko@durham.ac.uk)
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
#ifndef SWIFT_SPIN_JET_BH_IACT_H
#define SWIFT_SPIN_JET_BH_IACT_H

/* Local includes */
#include "black_holes_parameters.h"
#include "engine.h"
#include "equation_of_state.h"
#include "gravity.h"
#include "gravity_iact.h"
#include "hydro.h"
#include "random.h"
#include "rays.h"
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
 * @param step The current time-step.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_density(
    const float r2, const float *dx, const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  float wi, wi_dx;

  /* Get r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the number of neighbours */
  bi->num_ngbs++;

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Contribution to the BH gas density */
  bi->rho_gas += mj * wi;

  /* Contribution to the total neighbour mass */
  bi->ngb_mass += mj;

  /* Neighbour's sound speed */
  float cj;
  if (bh_props->use_subgrid_gas_properties && engine_current_step >= 0) {
    const float pressure_j = hydro_get_comoving_pressure(pj);
    const float subgrid_dens = cooling_get_subgrid_density(pj, xpj);
    cj = gas_soundspeed_from_pressure(
        subgrid_dens * cosmo->a * cosmo->a * cosmo->a, pressure_j);
  } else {
    cj = hydro_get_comoving_soundspeed(pj);
  }

  /* Contribution to the smoothed sound speed */
  bi->sound_speed_gas += mj * cj * wi;

  if (cj * cosmo->a_factor_sound_speed > bh_props->sound_speed_hot_gas_min) {
    bi->sound_speed_gas_hot += mj * cj * wi;
    bi->rho_gas_hot += mj * wi;
  }

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we do include a Hubble term) */
  const float dv[3] = {pj->v[0] - bi->v[0], pj->v[1] - bi->v[1],
                       pj->v[2] - bi->v[2]};

  const float a = cosmo->a;
  const float H = cosmo->H;
  const float a2H = a * a * H;

  /* Calculate the velocity with the Hubble flow */
  const float v_plus_H_flow[3] = {a2H * dx[0] + dv[0], a2H * dx[1] + dv[1],
                                  a2H * dx[2] + dv[2]};

  /* Contribution to the smoothed velocity (gas w.r.t. black hole) */
  bi->velocity_gas[0] += mj * dv[0] * wi;
  bi->velocity_gas[1] += mj * dv[1] * wi;
  bi->velocity_gas[2] += mj * dv[2] * wi;

  /* Contribution to the specific angular momentum of gas, which is later
   * converted to the circular velocity at the smoothing length */
  bi->spec_angular_momentum_gas[0] -= mj * wi * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->spec_angular_momentum_gas[1] -= mj * wi * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->spec_angular_momentum_gas[2] -= mj * wi * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Contribution to the smoothed squared relative velocity (for dispersion)
   * We will convert this to actual dispersion later. */
  const float norm_v2 = v_plus_H_flow[0] * v_plus_H_flow[0] +
                        v_plus_H_flow[1] * v_plus_H_flow[1] +
                        v_plus_H_flow[2] * v_plus_H_flow[2];

  bi->velocity_dispersion_gas += norm_v2 * wi * mj;

  if (bh_props->use_multi_phase_bondi) {
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

  /* Need to compute gas vorticity around the black hole, in analogy to
   * calculation for Balsara switch in hydro part */

  /* Factor to make sure we get curl, not angular momentum */
  const float faci = mj * wi_dx * r_inv;

  /* Compute dv cross r */
  const float v_cross_r[3] = {dv[1] * dx[2] - dv[2] * dx[1],
                              dv[2] * dx[0] - dv[0] * dx[2],
                              dv[0] * dx[1] - dv[1] * dx[0]};

  bi->curl_v_gas[0] += faci * v_cross_r[0];
  bi->curl_v_gas[1] += faci * v_cross_r[1];
  bi->curl_v_gas[2] += faci * v_cross_r[2];

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif

  /* Gas particle id */
  const long long gas_id = pj->id;

  /* Choose AGN feedback model */
  switch (bh_props->feedback_model) {
    case AGN_isotropic_model: {
      /* Compute arc lengths in AGN isotropic feedback and collect
       * relevant data for later use in the feedback_apply loop */

      /* Loop over rays */
      for (int i = 0; i < spinjet_blackhole_number_of_rays; i++) {

        /* We generate two random numbers that we use
        to randomly select the direction of the ith ray */

        /* Random number in [0, 1[ */
        const double rand_theta = random_unit_interval_part_ID_and_index(
            bi->id, i, ti_current,
            random_number_isotropic_AGN_feedback_ray_theta);

        /* Random number in [0, 1[ */
        const double rand_phi = random_unit_interval_part_ID_and_index(
            bi->id, i, ti_current,
            random_number_isotropic_AGN_feedback_ray_phi);

        /* Compute arc length */
        ray_minimise_arclength(dx, r, bi->rays + i,
                               /*ray_type=*/ray_feedback_thermal, gas_id,
                               rand_theta, rand_phi, mj, /*ray_ext=*/NULL,
                               /*v=*/NULL);
      }
      break;
    }
    case AGN_minimum_distance_model: {
      /* Compute the size of the array that we want to sort. If the current
       * function is called for the first time (at this time-step for this BH),
       * then bi->num_ngbs = 1 and there is nothing to sort. Note that the
       * maximum size of the sorted array cannot be larger then the maximum
       * number of rays. */
      const int arr_size = min(bi->num_ngbs, spinjet_blackhole_number_of_rays);

      /* Minimise separation between the gas particles and the BH. The rays
       * structs with smaller ids in the ray array will refer to the particles
       * with smaller distances to the BH. */
      ray_minimise_distance(r, bi->rays, arr_size, gas_id, mj);
      break;
    }
  }

  const int arr_size_jet = min(bi->num_ngbs, spinjet_blackhole_number_of_rays);

  /* Scalar product of the spin vector and position vector of the gas particle
     relative to the BH */
  float cosine_theta = -dx[0] * bi->angular_momentum_direction[0] -
                       dx[1] * bi->angular_momentum_direction[1] -
                       dx[2] * bi->angular_momentum_direction[2];

  /* Norm of the scalar product (ang. mom. dir. is already normalized) */
  const float norm = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

  /* Apply the norm to find the cosine of the angle between the two vectors,
     if norm is 0 then manually set to small value */
  if (norm > 0.) {
    cosine_theta = cosine_theta / norm;
  } else {
    cosine_theta = 0.001;
  }

  /* Define a variable which will be used in ray sorting to make sure that
     jet-kicked particles always end up at the end of the rays*/
  float ray_jet_correction = 0.;

  /* Define variables that we will minimize. Two because we have two rays. */
  float quantity_to_minimize = 0.;
  float quantity_to_minimize_pos = 0.;

  /* Calculate relative velocity of particle and BH, to be used to see if
     this particle was recently kicked */
  const float relative_velocity =
      sqrtf((bi->v[0] - pj->v[0]) * (bi->v[0] - pj->v[0]) +
            (bi->v[1] - pj->v[1]) * (bi->v[1] - pj->v[1]) +
            (bi->v[2] - pj->v[2]) * (bi->v[2] - pj->v[2])) *
      cosmo->a_inv;

  /* Choose AGN jet feedback model. Here we calculate the quantities to
     minimize, depending on the model. We calculate two numbers, one for each
     side of the BH. Particles are prioritized to be kicked from the 'same'
     hemisphere as defined by the BH spin vector. However, there may be cases
     in which one hemisphere is empty, so particles from the other side are
     used. If the particle is on the 'wrong' side of the BH, relative to the
     spin vector, we still put it in each ray, but they are pushed to the end
     (modulo other particles that were already kicked). We push them to the end
     so that they can still be used if the other hemisphere is empty. The way
     we push them to the end of the ray is by multiplying whatever quantity is
     being minimized by arbitrary numbers (large or small, depending on what is
     being minimized). This way, these particles never compete with the ones
     that are in the correct hemisphere of the BH. */
  switch (bh_props->jet_feedback_model) {
    case AGN_jet_minimum_distance_model: {

      /* Check if the relative velocity is a significant fraction of the jet
         launching velocity. If it is, set the ray correction variable to some
         arbitrarily large value. */
      if (relative_velocity > 0.8 * bi->v_jet) {
        ray_jet_correction = 1e11 * bi->h;
      }

      /* In this case we minimize using particle separations from the BH, with
         the order closest --> farthest.  */
      if (cosine_theta < 0) {
        quantity_to_minimize = r + ray_jet_correction;
      } else {
        quantity_to_minimize_pos = r + ray_jet_correction;
      }
      break;
    }
    case AGN_jet_maximum_distance_model: {

      /* Check if the relative velocity is a significant fraction of the jet
         launching velocity. If it is, set the ray correction variable to some
         arbitrarily large value. */
      if (relative_velocity > 0.8 * bi->v_jet) {
        ray_jet_correction = 1e13 * 1. / bi->h;
      }

      /* In this case we minimize using particle separations from the BH, with
         the order farthest --> closest  */
      if (cosine_theta < 0) {
        quantity_to_minimize = r_inv + ray_jet_correction;
      } else {
        quantity_to_minimize_pos = r_inv + ray_jet_correction;
      }
      break;
    }
    case AGN_jet_spin_axis_model: {

      /* Check if the relative velocity is a significant fraction of the jet
         launching velocity. If it is, set the ray correction variable to some
         arbitrarily large value. */
      if (relative_velocity > 0.8 * bi->v_jet) {
        ray_jet_correction = 1e3;
      }

      /* In this case we minimize using the angle (cosine) between the position
         vector of the particle (relative to the BH) and the spin vector of the
         BH. I.e. we launch particles along the spin axis, regardless of the
         distances from the BH. */
      if (cosine_theta < 0) {
        quantity_to_minimize = cosine_theta + ray_jet_correction;
      } else {
        quantity_to_minimize_pos = -cosine_theta + ray_jet_correction;
      }
      break;
    }
    case AGN_jet_minimum_density_model: {

      /* Check if the relative velocity is a significant fraction of the jet
         launching velocity. If it is, set the ray correction variable to some
         arbitrarily large value. */
      if (relative_velocity > 0.8 * bi->v_jet) {
        ray_jet_correction = 1e15 * pj->rho;
      }

      /* In this case we minimize using particle densities, i.e. we target the
         low-density gas.  */
      if (cosine_theta < 0) {
        quantity_to_minimize = pj->rho + ray_jet_correction;
      } else {
        quantity_to_minimize_pos = pj->rho + ray_jet_correction;
      }
      break;
    }
  }

  /* Do the actual minimization. */
  if (cosine_theta < 0) {
    ray_minimise_distance(quantity_to_minimize, bi->rays_jet, arr_size_jet,
                          gas_id, pj->mass);
  } else {
    ray_minimise_distance(quantity_to_minimize_pos, bi->rays_jet_pos,
                          arr_size_jet, gas_id, pj->mass);
  }
}

/**
 * @brief Swallowing interaction between two particles (non-symmetric).
 *
 * Function used to identify the gas particle that this BH may move towards.
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
 * @param step The current time-step.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_repos(
    const float r2, const float *dx, const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  float wi;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
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

      float potential = pj->black_holes_data.potential;

      if (bh_props->correct_bh_potential_for_repositioning) {
        /* Let's not include the contribution of the BH
         * itself to the potential of the gas particle */

        /* Note: This assumes the BH and gas have the same
         * softening, which is currently true */
        const float eps = gravity_get_softening(bi->gpart, grav_props);
        const float eps2 = eps * eps;
        const float eps_inv = 1.f / eps;
        const float eps_inv3 = eps_inv * eps_inv * eps_inv;
        const float BH_mass = bi->mass;

        /* Compute the Newtonian or truncated potential the BH
         * exherts onto the gas particle */
        float dummy, pot_ij;
        runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, BH_mass, &dummy,
                                 &pot_ij);

        /* Deduct the BH contribution */
        potential -= pot_ij * grav_props->G_Newton;
      }

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
 * @param step The current time-step.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(
    const float r2, const float *dx, const float hi, const float hj,
    struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  float wi;

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Check if the BH needs to be fed. If not, we're done here */
  const float bh_mass_deficit = bi->subgrid_mass - bi->mass_at_start_of_step;
  if (bh_mass_deficit <= 0) return;

  if (bh_props->use_nibbling) {

    /* If we do nibbling, things are quite straightforward. We transfer
     * the mass and all associated quantities right here. */

    const float bi_mass_orig = bi->mass;
    const float pj_mass_orig = hydro_get_mass(pj);

    /* Don't nibble from particles that are too small already */
    if (pj_mass_orig < bh_props->min_gas_mass_for_nibbling) return;

    /* Next line is equivalent to w_ij * m_j / Sum_j (w_ij * m_j) */
    const float particle_weight = hi_inv_dim * wi * pj_mass_orig / bi->rho_gas;
    float nibble_mass = bh_mass_deficit * particle_weight;

    /* We radiated away some of the accreted mass, so need to take slightly
     * more from the gas than the BH gained */
    const float excess_fraction = 1.0 / (1.0 - bi->radiative_efficiency);

    /* Need to check whether nibbling would push gas mass below minimum
     * allowed mass */
    float new_gas_mass = pj_mass_orig - nibble_mass * excess_fraction;
    if (new_gas_mass < bh_props->min_gas_mass_for_nibbling) {
      new_gas_mass = bh_props->min_gas_mass_for_nibbling;
      nibble_mass = (pj_mass_orig - bh_props->min_gas_mass_for_nibbling) /
                    excess_fraction;
    }

    /* Correct for nibbling the particle mass that is stored in rays */
    for (int i = 0; i < spinjet_blackhole_number_of_rays; i++) {
      if (bi->rays[i].id_min_length == pj->id) bi->rays[i].mass = new_gas_mass;
      if (bi->rays_jet[i].id_min_length == pj->id) {
        bi->rays_jet[i].mass = new_gas_mass;
      }
      if (bi->rays_jet_pos[i].id_min_length == pj->id) {
        bi->rays_jet_pos[i].mass = new_gas_mass;
      }
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

    const float nibbled_mass = nibble_mass * excess_fraction;
    const float nibbled_fraction = nibbled_mass / pj_mass_orig;

    /* Update the BH and also gas metal masses */
    struct chemistry_bpart_data *bi_chem = &bi->chemistry_data;
    struct chemistry_part_data *pj_chem = &pj->chemistry_data;
    chemistry_transfer_part_to_bpart(bi_chem, pj_chem, nibbled_mass,
                                     nibbled_fraction);

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
 * Function used to identify the BH particle that this BH may move towards.
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
runner_iact_nonsym_bh_bh_repos(const float r2, const float *dx, const float hi,
                               const float hj, struct bpart *bi,
                               struct bpart *bj, const struct cosmology *cosmo,
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

      float potential = bj->reposition.potential;

      if (bh_props->correct_bh_potential_for_repositioning) {

        /* Let's not include the contribution of the BH i
         * to the potential of the BH j */

        const float eps = gravity_get_softening(bi->gpart, grav_props);
        const float eps2 = eps * eps;
        const float eps_inv = 1.f / eps;
        const float eps_inv3 = eps_inv * eps_inv * eps_inv;
        const float BH_mass = bi->mass;

        /* Compute the Newtonian or truncated potential the BH
         * exherts onto the gas particle */
        float dummy, pot_ij;
        runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, BH_mass, &dummy,
                                 &pot_ij);

        /* Deduct the BH contribution */
        potential -= pot_ij * grav_props->G_Newton;
      }

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

    /* Maximum velocity difference between BHs allowed to merge */
    float v2_threshold;

    if (bh_props->merger_threshold_type == BH_mergers_circular_velocity) {

      /* 'Old-style' merger threshold using circular velocity at the
       * edge of the more massive BH's kernel (note: we are using the kernel
       * support radius here and not just the smoothing length). */
      v2_threshold = G_Newton * M / (kernel_gamma * h);
    } else {

      /* Arguably better merger threshold using the escape velocity at
       * the distance between the BHs */

      if (bh_props->merger_threshold_type == BH_mergers_escape_velocity) {
        /* Standard formula (not softening BH interactions) */
        v2_threshold = 2.f * G_Newton * M / sqrt(r2);
      } else if (bh_props->merger_threshold_type ==
                 BH_mergers_dynamical_escape_velocity) {
        /* General two-body escape velocity based on dynamical masses */
        v2_threshold = 2.f * G_Newton * (bi->mass + bj->mass) / sqrt(r2);
      } else {
        error("Unexpected BH merger threshold type!");
        v2_threshold = 0.f;
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
 * @param step The current time-step.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_feedback(
    const float r2, const float *dx, const float hi, const float hj,
    const struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  /* Number of energy injections per BH per time-step */
  const int num_energy_injections_per_BH =
      bi->to_distribute.AGN_number_of_energy_injections;

  /* Are we doing some feedback? */
  if (num_energy_injections_per_BH > 0) {

    /* Number of energy injections that have reached this gas particle */
    int num_of_energy_inj_received_by_gas = 0;

    /* Find out how many rays (= energy injections) this gas particle
     * has received */
    for (int i = 0; i < num_energy_injections_per_BH; i++) {
      if (pj->id == bi->rays[i].id_min_length)
        num_of_energy_inj_received_by_gas++;
    }

    /* If the number of received rays is non-zero, inject
     * AGN energy in thermal form */
    if (num_of_energy_inj_received_by_gas > 0) {

      /* Compute new energy per unit mass of this particle
       * The energy the particle receives is proportional to the number of rays
       * (num_of_energy_inj_received_by_gas) to which the particle was found to
       * be closest. */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float delta_u = bi->to_distribute.AGN_delta_u *
                            (float)num_of_energy_inj_received_by_gas;
      const double u_new = u_init + delta_u;

      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, /*pfloor=*/NULL,
                                                 u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Store the feedback energy */
      const double delta_energy = delta_u * hydro_get_mass(pj);
      tracers_after_black_holes_feedback(pj, xpj, with_cosmology, cosmo->a,
                                         time, delta_energy);

      /* message( */
      /*     "We did some AGN heating! id %llu BH id %llu probability " */
      /*     " %.5e  random_num %.5e du %.5e du/ini %.5e", */
      /*     pj->id, bi->id, prob, rand, delta_u, delta_u / u_init); */

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }
  }

  /* Number of jet injections per BH per time-step */
  const int num_jet_injections_per_BH =
      bi->to_distribute.AGN_number_of_jet_injections;

  /* Are we doing some jet feedback? */
  if (num_jet_injections_per_BH > 0) {

    /* Number of jet injections that have reached this gas particle */
    int num_of_jet_inj_received_by_gas = 0;

    /* Define a variable to assign a velocity kick direction depending
       on which side of the BH smoothing kernel the particle is */
    float direction = 0.;

    /* Find out if this gas particle has received any jet injections (rays).
    Loop through num_jet_injections divided by 2 because of two sets of rays */
    for (int i = 0; i < num_jet_injections_per_BH / 2; i++) {
      if (pj->id == bi->rays_jet[i].id_min_length) {

        num_of_jet_inj_received_by_gas++;

        /*This particle is in the 'negative' hemisphere (pointing away from the
          spin vector of the BH), so it receives a negative kick direction */
        direction = -1.;
      }
    }

    for (int i = 0; i < num_jet_injections_per_BH / 2; i++) {
      if (pj->id == bi->rays_jet_pos[i].id_min_length) {

        num_of_jet_inj_received_by_gas++;

        /* This particle is in the 'positive' hemisphere (pointing in the
           direction of the spin vector of the BH), so it receives a positive
           kick direction */
        direction = 1.;
      }
    }

    /* If the number of received rays is non-zero, inject
     * AGN jet energy as a kinetic kick */
    if (num_of_jet_inj_received_by_gas > 0) {

      /* Get the kinetic energy per unit mass */
      const float delta_u_jet = bi->to_distribute.AGN_delta_u_jet *
                                (float)num_of_jet_inj_received_by_gas;

      /* Get the (physical) kick velocity, and convert to code units */
      const float vel_kick = sqrtf(2. * delta_u_jet) * cosmo->a;

      /* Compute velocity kick direction using the previously generated
       * jet direction.*/
      float vel_kick_direction[3];

      /* Include the -1./1. factor (direction) which accounts for kicks in the
       * opposite direction of the spin vector */
      vel_kick_direction[0] = direction * bi->jet_direction[0];
      vel_kick_direction[1] = direction * bi->jet_direction[1];
      vel_kick_direction[2] = direction * bi->jet_direction[2];

      /* Get the initial velocity in the frame of the black hole */
      const float v_init[3] = {xpj->v_full[0] - bi->v[0],
                               xpj->v_full[1] - bi->v[1],
                               xpj->v_full[2] - bi->v[2]};

      /* We compute this final velocity by requiring that the final energy and
       * the inital one differ by the energy received by the particle, i.e.
       *
       *        (pi + delta_pi)^2 / (2m) - pi^2 / (2m) = u,
       *
       * u here being the energy per unit mass received by the particle. pi is
       * the initial momentum, and the momenta terms are expressed in vector
       * form. The equation, if expressed in terms of velocities, amounts to
       *
       *   norm(delta_v) + 2 * norm(delta_v) * norm(v_i) * cos_theta_v = v_k^2.
       *
       * Here, delta_v is the change in velocity which we wish to apply, v_i is
       * the initial velocity, cos_theta_v the cosine of the angle between the
       * two and v_k^2 is the vel_kick term computed from the energy received
       * by the particle. The delta_v applied to the particle will differ in
       * norm from the parameter v_j used for jet feedback for two reasons:
       * 1) v_k is slightly different from v_j if the particle mass is not
       * equal to the mean neighbour mass, and 2) the presence of the initial
       * velocity means we need to increase the magnitude of the velocity by
       * less than v_j in order to increase its energy by (1/2)mv_j^2. We solve
       * the above quadratic equation for the norm(delta_v). We begin by
       * calculating norm(v_i) * cos_theta_v, which is the initial velocity
       * projected onto the velocity kick direction. */
      const float v_init_proj = v_init[0] * vel_kick_direction[0] +
                                v_init[1] * vel_kick_direction[1] +
                                v_init[2] * vel_kick_direction[2];
      const float delta_v =
          sqrtf(max(0., v_init_proj * v_init_proj + vel_kick * vel_kick)) -
          v_init_proj;

      /* Calculate final velocity by adding delta_v in the direction of the kick
       */
      xpj->v_full[0] += delta_v * vel_kick_direction[0];
      xpj->v_full[1] += delta_v * vel_kick_direction[1];
      xpj->v_full[2] += delta_v * vel_kick_direction[2];

#ifdef SWIFT_DEBUG_CHECKS
      message(
          "Black hole with id %lld kicked particle with id %lld , with a final "
          "velocity of (%f, %f, %f).",
          bi->id, pj->id, xpj->v_full[0], xpj->v_full[1], xpj->v_full[2]);
#endif

      /* Store the jet energy and other variables of interest */
      const double delta_energy_jet = delta_u_jet * hydro_get_mass(pj);
      tracers_after_jet_feedback(pj, xpj, with_cosmology, cosmo->a, time,
                                 delta_energy_jet, vel_kick, bi->accretion_mode,
                                 bi->id);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Update the signal velocity */
      hydro_set_v_sig_based_on_velocity_kick(pj, cosmo,
                                             sqrtf(2. * delta_u_jet));

      /* Synchronize particle on the time-line */
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

#endif /* SWIFT_SPIN_JET_BH_IACT_H */
