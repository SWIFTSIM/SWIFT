/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2021 Edo Altamura (edoardo.altamura@manchester.ac.uk)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_SIMBA_BH_IACT_H
#define SWIFT_SIMBA_BH_IACT_H

/* Local includes */
#include "black_holes_parameters.h"
#include "entropy_floor.h"
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
 * @param bi First particle (black hole).
 * @param sj Second particle (stars, not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_stars_density(const float r2, const float dx[3],
                                    struct bpart *bi, const struct spart *sj) { }

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param bi First particle (black hole).
 * @param sj Second particle (stars, not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_stars_bulge(const float r2, const float dx[3],
                                  struct bpart *bi, const struct spart *sj) { }

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
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time, const double time_base) {

  /* Ignore decoupled winds in density computation */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  float wi, wi_dx;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
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

  /* Contribution to the smoothed sound speed */
  const float cj = hydro_get_comoving_soundspeed(pj);
  bi->sound_speed_gas += mj * wi * cj;

  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);

  /* Contribution to the smoothed internal energy */
  bi->internal_energy_gas += mj * uj * wi;

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {pj->v[0] - bi->v[0], pj->v[1] - bi->v[1],
                       pj->v[2] - bi->v[2]};

  /* Account for hot and cold gas surrounding the SMBH */
  const float Tj =
      uj * cosmo->a_factor_internal_energy / bh_props->temp_to_u_factor;
  int is_hot_gas = 0;
  /* Check whether we are close to the entropy floor. If we are, we
   * classify the gas as cold regardless of temperature */
  if (Tj > bh_props->environment_temperature_cut) {
    const float T_EoS = entropy_floor_temperature(pj, cosmo, floor_props);
    if (Tj > T_EoS * bh_props->fixed_T_above_EoS_factor) {
      is_hot_gas = 1;
    }
  }

  /* Star forming gas is never considered "hot" */
  if (pj->sf_data.SFR > 0.f) is_hot_gas = 0;

  if (is_hot_gas) {
    bi->hot_gas_mass += mj;
    bi->hot_gas_internal_energy += mj * uj; /* Not kernel weighted */
  } else {
    bi->cold_gas_mass += mj;
  }

  /* Sum up cold disk mass corotating relative to ang mom computed so far.  This is not fully
   * accurate but it is convenient and probably not too bad */
  const double Lx = mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  const double Ly = mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  const double Lz = mj * (dx[0] * dv[1] - dx[1] * dv[0]);
  const double proj = Lx * bi->angular_momentum_gas[0] + Ly * bi->angular_momentum_gas[1] + Lz * bi->angular_momentum_gas[2];
  if ((proj > 0.f) && (is_hot_gas == 0)) bi->cold_disk_mass += mj;

  /* Gas angular momentum in kernel */
  bi->angular_momentum_gas[0] += mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->angular_momentum_gas[1] += mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->angular_momentum_gas[2] += mj * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Contribution to the specific angular momentum of gas, which is later
   * converted to the circular velocity at the smoothing length */
  bi->circular_velocity_gas[0] -= mj * wi * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->circular_velocity_gas[1] -= mj * wi * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->circular_velocity_gas[2] -= mj * wi * (dx[0] * dv[1] - dx[1] * dv[0]);

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

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif

}

/**
 * @brief Repositioning interaction between two particles (non-symmetric).
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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_repos(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time, const double time_base) {

  /* Ignore decoupled wind particles */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  float wi;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

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
                     bi->sound_speed_gas * bi->sound_speed_gas *
                     cosmo->a_factor_sound_speed * cosmo->a_factor_sound_speed;

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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time, const double time_base) {

  /* IMPORTANT: Do not even consider wind particles for accretion/feedback */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bi->group_data.mass <= 0.f) return;
  
  /* If there is no gas, skip */
  if (bi->rho_gas <= 0.f) return;

  float wi;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Get particle time-step */
  double dt;
  if (with_cosmology) {
    const integertime_t ti_step = get_integer_timestep(bi->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(ti_current - 1, bi->time_bin);

    dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
  } else {
    dt = get_timestep(bi->time_bin, time_base);
  }

  /* Probability to swallow this particle */
  float prob = -1.f;
  float f_accretion = bi->f_accretion;
  if (f_accretion <= 0.f) return;

  /* Radiation was already accounted for in bi->subgrid_mass
   * so if is is bigger than bi->mass we can simply
   * flag particles to eat and satisfy the mass constraint.
   *
   * If bi->subgrid_mass < bi->mass then there is a problem,
   * but we use the full accretion rate to select particles
   * and then don't actually take anything away from them.
   * The bi->mass variable is decreased previously to account
   * for the radiative losses.
   */
  const float mass_deficit = bi->subgrid_mass - (bi->mass_at_start_of_step +
                                                 bi->mass_accreted_this_step);
  if (mass_deficit >= 0.f) {
    /* Don't nibble from particles that are too small already */
    if (hydro_get_mass(pj) < bh_props->min_gas_mass_for_nibbling) return;

    prob = (mass_deficit / bi->f_accretion) * hi_inv_dim * wi / bi->rho_gas;
  } else {
    prob = 0.f;
  }

  if (bi->subgrid_mass - bi->mass < 0.f) {
    prob = ((1.f - bi->f_accretion) / bi->f_accretion) * bi->accretion_rate *
           dt * (hi_inv_dim * wi / bi->rho_gas);
    /* We do NOT accrete when subgrid_mass < physical_mass
     * but we still kick.
     */
    f_accretion = 0.f;
  }

  /* Draw a random number (Note mixing both IDs) */
  const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                          random_number_BH_swallow);

  /* Are we lucky? */
  if (rand < prob) {

    if (f_accretion > 0.f) {
      const float bi_mass_orig = bi->mass;
      const float pj_mass_orig = hydro_get_mass(pj);
      const float nibbled_mass = f_accretion * pj_mass_orig;
      const float new_gas_mass = pj_mass_orig - nibbled_mass;
      /* Don't go below the minimum for stability */
      if (new_gas_mass < bh_props->min_gas_mass_for_nibbling) return;

      bi->mass += nibbled_mass;
      hydro_set_mass(pj, new_gas_mass);
      bi->mass_accreted_this_step += nibbled_mass;

      /* Add the angular momentum of the accreted gas to the BH total.
       * Note no change to gas here. The cosmological conversion factors for
       * velocity (a^-1) and distance (a) cancel out, so the angular momentum
       * is already in physical units. */
      const float dv[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                           bi->v[2] - pj->v[2]};
      bi->swallowed_angular_momentum[0] +=
          nibbled_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      bi->swallowed_angular_momentum[1] +=
          nibbled_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      bi->swallowed_angular_momentum[2] +=
          nibbled_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      /* Update the BH momentum and velocity. Again, no change to gas here. */
      const float bi_mom[3] = {
          bi_mass_orig * bi->v[0] + nibbled_mass * pj->v[0],
          bi_mass_orig * bi->v[1] + nibbled_mass * pj->v[1],
          bi_mass_orig * bi->v[2] + nibbled_mass * pj->v[2]};

      bi->v[0] = bi_mom[0] / bi->mass;
      bi->v[1] = bi_mom[1] / bi->mass;
      bi->v[2] = bi_mom[2] / bi->mass;

      /* Update the BH and also gas metal masses */
      struct chemistry_bpart_data *bi_chem = &bi->chemistry_data;
      struct chemistry_part_data *pj_chem = &pj->chemistry_data;
      chemistry_transfer_part_to_bpart(bi_chem, pj_chem, nibbled_mass,
                                       nibbled_mass / pj_mass_orig);

    } else { /* When f_accretion <= 0.f, but bi->f_accretion > 0.f */
      bi->mass_accreted_this_step += bi->f_accretion * hydro_get_mass(pj);
    }

    /* This particle is swallowed by the BH with the largest ID of all the
     * candidates wanting to swallow it */
    if (pj->black_holes_data.swallow_id < bi->id) {
      pj->black_holes_data.swallow_id = bi->id;
    } else {
      message(
          "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
          "swallow id=%lld)",
          bi->id, pj->id, pj->black_holes_data.swallow_id);
    }
  }
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
runner_iact_nonsym_bh_bh_repos(const float r2, const float dx[3],
                               const float hi, const float hj, struct bpart *bi,
                               const struct bpart *bj,
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
                     bi->sound_speed_gas * bi->sound_speed_gas *
                     cosmo->a_factor_sound_speed * cosmo->a_factor_sound_speed;

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
runner_iact_nonsym_bh_bh_swallow(const float r2, const float dx[3],
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

    /* Merge if gravitationally bound AND if within max distance
     * Note that we use the kernel support here as the size and not just the
     * smoothing length */

    /* Maximum velocity difference between BHs allowed to merge */
    float v2_threshold;

    if (bh_props->merger_threshold_type == BH_mergers_circular_velocity) {

      /* 'Old-style' merger threshold using circular velocity at the
       * edge of the more massive BH's kernel */
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
        /* Cannot happen! */
#ifdef SWIFT_DEBUG_CHECKS
        error("Invalid choice of BH merger threshold type");
#endif
        v2_threshold = 0.f;
      }
    } /* Ends sections for different merger thresholds */

    if ((v2_pec < v2_threshold) && (r2 < max_dist_merge2)) {

      /* This particle is swallowed by the BH with the largest mass of all the
       * candidates wanting to swallow it (we use IDs to break ties)*/
      if ((bj->merger_data.swallow_mass < bi->subgrid_mass) ||
          (bj->merger_data.swallow_mass == bi->subgrid_mass &&
           bj->merger_data.swallow_id < bi->id)) {

#ifdef SWIFT_DEBUG_CHECKS
        message("BH %lld wants to swallow BH particle %lld", bi->id, bj->id);
#endif

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
runner_iact_nonsym_bh_gas_feedback(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time, const double time_base) {

  /* This shouldn't happen, but just be sure anyway */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bi->group_data.mass <= 0.f) return;

  /* No energy left; no feedback */
  if (bi->energy_reservoir <= 0.f) return;

  /* No distance, no feedback */
  if (r2 <= 0.f) return;

  /* Do X-ray feedback first */
  if (pj->black_holes_data.swallow_id != bi->id) {
    /* We were not lucky for kick, but we might be lucky for X-ray feedback */
    if (bi->v_kick > bh_props->xray_heating_velocity_threshold) {
      const float group_gas_mass = bi->group_data.mass -
                                   bi->group_data.stellar_mass;

      const float f_gas = group_gas_mass / bi->group_data.mass;

      float f_rad_loss = bh_props->xray_radiation_loss * 
                         (bh_props->xray_f_gas_limit - f_gas) / 
                         bh_props->xray_f_gas_limit;
      if (f_rad_loss > bh_props->xray_radiation_loss) {
        f_rad_loss = bh_props->xray_radiation_loss;
      }

      if (f_rad_loss < 0.f) return;

      /* Get particle time-step */
      double dt;
      if (with_cosmology) {
        const integertime_t ti_step = get_integer_timestep(bi->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current - 1, bi->time_bin);

        dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
      } else {
        dt = get_timestep(bi->time_bin, time_base);
      }
  
      /* Hydrogen number density (X_H * rho / m_p) [cm^-3] */
      const float n_H_cgs =
          hydro_get_physical_density(pj, cosmo) * bh_props->rho_to_n_cgs;
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const float T_gas_cgs =
          u_init / (bh_props->temp_to_u_factor * bh_props->T_K_to_int);

      double du_xray_phys = black_holes_compute_xray_feedback(
          bi, pj, bh_props, cosmo, dx, dt, n_H_cgs, T_gas_cgs);

      /* Limit the amount of heating BEFORE dividing to avoid numerical
       * instability */
      if (du_xray_phys > bh_props->xray_maximum_heating_factor * u_init) {
        du_xray_phys = bh_props->xray_maximum_heating_factor * u_init;
      }

      /* Account for X-rays lost due to radiation */
      du_xray_phys *= f_rad_loss;

      /* Conserve energy: reduce available BH energy */
      const double dE_this_step = du_xray_phys * pj->mass;
      if (bi->energy_reservoir > dE_this_step) {
	      bi->energy_reservoir -= dE_this_step;
      }
      else {
	      /* Not enough energy left; reduce feedback energy */
        du_xray_phys = bi->energy_reservoir / pj->mass;
	      bi->energy_reservoir = 0.f;
      }

      /* Look for cold dense gas. Then push it. */

      /* Check whether we are close to the entropy floor. If we are, we
       * classify the gas as cold regardless of temperature.
       * All star forming gas is considered cold.
       */
      const float T_EoS_cgs = entropy_floor_temperature(pj, cosmo, floor_props)
                                  / bh_props->T_K_to_int;
      if ((n_H_cgs > bh_props->xray_heating_n_H_threshold_cgs &&
            (T_gas_cgs < bh_props->xray_heating_T_threshold_cgs ||
                T_gas_cgs < T_EoS_cgs * bh_props->fixed_T_above_EoS_factor)) ||
            pj->sf_data.SFR > 0.f) {
	      /* compute kick velocity */
        const float dv_phys = 2.f * sqrtf(
                                  bh_props->xray_kinetic_fraction * 
                                  du_xray_phys);
        /* Push gas radially */
        const float r = sqrtf(r2);
        const float dv_comoving = dv_phys * cosmo->a;
        const float prefactor = dv_comoving / r;
        pj->v_full[0] += prefactor * dx[0];
        pj->v_full[1] += prefactor * dx[1];
        pj->v_full[2] += prefactor * dx[2];

        /* Update the signal velocity of the particle based on the velocity
         * kick. */
        hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, dv_phys);

#ifdef SIMBA_DEBUG_CHECKS
        message("BH_XRAY_KICK: bid=%lld, pid=%lld, %g km/s", 
                bi->id, pj->id,
                dv_phys / bh_props->kms_to_internal);
#endif
      }
      /* If we are not in ISM, then dump heat */
      du_xray_phys *= (1. - bh_props->xray_kinetic_fraction);

      const double u_new = u_init + du_xray_phys;

      /* Do the energy injection. */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      if (bh_props->xray_shutoff_cooling) {
        /* u_init is physical so cs_physical is physical */
        const double cs_physical = gas_soundspeed_from_internal_energy(pj->rho, u_new);
  
        /* a_factor_sound_speed converts cs_physical to internal (comoving) units) */
        pj->feedback_data.cooling_shutoff_delay_time = maxf(
              cosmo->a_factor_sound_speed * (pj->h / cs_physical),
              dt); /* BH timestep as a lower limit */
      }

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);

#ifdef SIMBA_DEBUG_CHECKS
      const double T_gas_final_cgs = 
          u_new / (bh_props->temp_to_u_factor * bh_props->T_K_to_int);
      message("BH_XRAY_HEAT: bid=%lld, pid=%lld, T_init %g K, T_new %g K, T_new/T_init=%g, dt_shutoff=%g Myr",
              bi->id, pj->id,
              T_gas_cgs,
              T_gas_final_cgs,
              T_gas_final_cgs / T_gas_cgs,
              pj->feedback_data.cooling_shutoff_delay_time * bh_props->time_to_Myr);
#endif

    }
  } 
  else { /* Below is swallow_id = id for particle/bh */

    /* Save gas density and entropy before feedback */
    tracers_before_black_holes_feedback(pj, xpj, cosmo->a);

    /* Kick along the angular momentum axis of gas in the kernel */
    float norm =
        sqrtf(bi->angular_momentum_gas[0] * bi->angular_momentum_gas[0] +
              bi->angular_momentum_gas[1] * bi->angular_momentum_gas[1] +
              bi->angular_momentum_gas[2] * bi->angular_momentum_gas[2]);

    /* No norm, no wind */
    if (norm <= 0.f) return;

    /* Kick particle */
    norm = 1.f / norm;
    const double random_number =
        random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
    const float dirsign = (random_number > 0.5) ? 1.f : -1.f;
    double dv = bi->v_kick * cosmo->a * dirsign;
    const double dE_this_step = 0.5f * pj->mass * dv * dv;
    if (bi->energy_reservoir > dE_this_step) {
      bi->energy_reservoir -= 0.5f * pj->mass * dv * dv;
    }
    else if (bi->energy_reservoir > 0.f) {
      dv = sqrtf(2.f * bi->energy_reservoir / pj->mass);
      bi->energy_reservoir = 0.f;
    }

    pj->v_full[0] += dv * bi->angular_momentum_gas[0] * norm;
    pj->v_full[1] += dv * bi->angular_momentum_gas[1] * norm;
    pj->v_full[2] += dv * bi->angular_momentum_gas[2] * norm;

#ifdef SIMBA_DEBUG_CHECKS
    message("BH_KICK: bid=%lld kicking pid=%lld, v_kick=%g km/s",
       bi->id, pj->id, bi->v_kick / bh_props->kms_to_internal);
#endif

    /* Set delay time */
    pj->feedback_data.decoupling_delay_time =
        bh_props->wind_decouple_time_factor *
        cosmology_get_time_since_big_bang(cosmo, cosmo->a);

    /* Update the signal velocity of the particle based on the velocity kick. */
    hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, bi->v_kick);

    /* If we have a jet, we heat! */
    if (bi->v_kick >= bh_props->jet_heating_velocity_threshold) {
#ifdef SIMBA_DEBUG_CHECKS
      message("BH_KICK_JET: bid=%lld kicking pid=%lld at v_kick=%g km/s",
        bi->id, pj->id, bi->v_kick / bh_props->kms_to_internal);
#endif

      float new_Tj = 0.f;
      /* Use the halo Tvir? */
      if (bh_props->scale_jet_temperature_with_mass) {
        const float mass_scaling = 
            bh_props->mass_to_solar_mass / bh_props->jet_temperature_mass_norm; 
        new_Tj = bh_props->jet_temperature *
                 powf(bi->subgrid_mass * mass_scaling, 2.0 / 3.0);
      } else {
        new_Tj = bh_props->jet_temperature; /* K */
      }

      /* Simba scales with velocity */
      new_Tj *= (bi->v_kick * bi->v_kick) /
                (bh_props->jet_velocity * bh_props->jet_velocity);

      /* Treat the jet temperature as an upper limit, in case v_kick > v_jet */
      if (new_Tj > bh_props->jet_temperature)
        new_Tj = bh_props->jet_temperature;

#ifdef SIMBA_DEBUG_CHECKS
      message("BH_KICK_JET_HEAT: bid=%lld heating pid=%lld to T=%g K",
        bi->id, pj->id, new_Tj);
#endif

      /* Compute new energy per unit mass of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      const double u_new = bh_props->temp_to_u_factor * new_Tj;

      /* Don't decrease the gas temperature if it's already hotter */
      if (u_new > u_init) {
        /* We are overwriting the internal energy of the particle */
        hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
        hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

        const double delta_energy = (u_new - u_init) * hydro_get_mass(pj);

        tracers_after_black_holes_feedback(pj, xpj, with_cosmology, cosmo->a,
                                           time, delta_energy);
      }
      pj->feedback_data.number_of_times_decoupled += 100000;
    }
    else {
      pj->feedback_data.number_of_times_decoupled += 1000;
    }

    /* Wind cannot be star forming */
    if (pj->sf_data.SFR > 0.f) {

      /* Record the current time as an indicator of when this particle was last
        star-forming. */
      if (with_cosmology) {
        pj->sf_data.SFR = -cosmo->a;
      } else {
        pj->sf_data.SFR = -time;
      }

    }

    /* Impose maximal viscosity */
    hydro_diffusive_feedback_reset(pj);

    /* Synchronize the particle on the timeline */
    timestep_sync_part(pj);

    /* IMPORTANT: The particle MUST NOT be swallowed. 
     * We are taking a f_accretion from each particle, and then
     * kicking the rest. We used the swallow marker as a temporary
     * passer in order to remember which particles have been "nibbled"
     * so that we can kick them out.
     */
    black_holes_mark_part_as_not_swallowed(&pj->black_holes_data);
  }
}

#endif /* SWIFT_SIMBA_BH_IACT_H */
