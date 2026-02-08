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
#ifndef SWIFT_OBSIDIAN_BH_IACT_H
#define SWIFT_OBSIDIAN_BH_IACT_H

/* Local includes */
#include "black_holes_parameters.h"
#include "black_holes_properties.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "gravity.h"
#include "gravity_iact.h"
#include "hydro.h"
#include "random.h"
#include "space.h"
#include "timestep_sync_part.h"
#include "tracers.h"

// #define OBSIDIAN_DEBUG_CHECKS

/**
 * @brief Set direction of black hole feedback kick
 *
 * @param bi Black hole particle.
 * @param pj Gas particle being kicked.
 * @param ti_current Current integer time value (for random numbers).
 * @param dir_flag Flag to choose direction: 0=rendom, 1=L_gas, 2=L_BH.
 * @param dir Direction of kick (returned).
 */
__attribute__((always_inline)) INLINE static float
black_hole_set_kick_direction(const struct bpart *bi, const struct part *pj,
                              const integertime_t ti_current,
                              const int dir_flag, float *dir) {

  float kick_dir = 1.f;
  double random_number = 1.;
  switch (dir_flag) {
    /* Isotropic */
    case 0: {
      const double random_for_theta =
          random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
      const double random_for_phi =
          random_unit_interval(bi->id, ti_current, random_number_BH_feedback);

      const float theta = acosf(2.f * random_for_theta - 1.f);
      const float phi = 2.f * M_PI * random_for_phi;

      dir[0] = sinf(theta) * cosf(phi);
      dir[1] = sinf(theta) * sinf(phi);
      dir[2] = cosf(theta);
      break;
    }

    /* Along the angular momentum vector of the gas */
    case 1:
      dir[0] = bi->angular_momentum_gas[0];
      dir[1] = bi->angular_momentum_gas[1];
      dir[2] = bi->angular_momentum_gas[2];
      random_number =
          random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
      kick_dir = (random_number > 0.5) ? 1.f : -1.f;
      break;

    /* Along the accreted angular momentum of the gas */
    case 2:
      dir[0] =
          bi->accreted_angular_momentum[0] + bi->swallowed_angular_momentum[0];
      dir[1] =
          bi->accreted_angular_momentum[1] + bi->swallowed_angular_momentum[1];
      dir[2] =
          bi->accreted_angular_momentum[2] + bi->swallowed_angular_momentum[2];
      random_number =
          random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
      kick_dir = (random_number > 0.5) ? 1.f : -1.f;
      break;

    /* Outwards from BH*/
    case 3:
      dir[0] = pj->x[0] - bi->x[0];
      dir[1] = pj->x[1] - bi->x[1];
      dir[2] = pj->x[2] - bi->x[2];
      break;

    default:
      error("dir_flag=%d but must be 0, 1, or 2", dir_flag);
      break;
  }

  return kick_dir;
}

/**
 * @brief Categorise gas temperature as hot or cold or neither for BH accretion
 *
 * @param bi First particle (black hole).
 * @param pj Second particle (gas, not updated).
 * @param cosmo The cosmological model.
 * @param bh_props The properties of the BH scheme
 */
__attribute__((always_inline)) INLINE static int
black_hole_gas_hot_or_cold(const struct bpart *bi, const struct part *pj,
    const struct cosmology *cosmo, const struct black_holes_props *bh_props)
{
  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);

  /* Determine if it is hot or cold gas surrounding the SMBH */
  int gas_state = 0;
  const float Tj =
      uj * cosmo->a_factor_internal_energy / bh_props->temp_to_u_factor;
  if (Tj > bh_props->environment_temperature_cut && pj->sf_data.SFR <= 0.f) {
    gas_state = 1;
  }
  /* Cold gas must be cool and SF'ing */
  if (Tj < bh_props->cold_gas_temperature_cut && pj->sf_data.SFR > 0.) {
    gas_state = -1;
  }

  return gas_state;
}
		
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
    const integertime_t ti_current, const double time) {

  /* Get the total gravitationally interacting mass within the kernel */
  const float mj = hydro_get_mass(pj);

  /* Compute total mass that contributes to the dynamical time */
  bi->gravitational_ngb_mass += mj;
  bi->num_gravitational_ngbs += 1;

  float wi, wi_dx;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Contribution to the BH gas density */
  bi->rho_gas += mj * wi;

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the smoothed sound speed */
  const float cj = hydro_get_comoving_soundspeed(pj);
  bi->sound_speed_gas += mj * wi * cj;

  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);

  /* Contribution to the smoothed internal energy */
  bi->internal_energy_gas += mj * uj * wi;

  /* weighting for feedback */
  bi->kernel_wt_sum += mj * wi;

  /* Contribution to the number of neighbours */
  bi->num_ngbs += 1;

  /* Contribution to the total neighbour mass */
  bi->ngb_mass += mj;

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {pj->v[0] - bi->v[0], pj->v[1] - bi->v[1],
                       pj->v[2] - bi->v[2]};

  /* Classify gas as hot or cold for accretion */
  const int gas_temperature_state = black_hole_gas_hot_or_cold(bi, pj, cosmo, bh_props);

  if (gas_temperature_state == 1) {
    bi->hot_gas_mass += mj;
    bi->hot_gas_internal_energy += mj * uj; /* Not kernel weighted */
  } 
  else if (gas_temperature_state == -1) {
#if COOLING_GRACKLE_MODE >= 2
    /* With subgrid ISM model, only allow H2 component to be accreted */
    bi->cold_gas_mass += mj * pj->cooling_data.subgrid_fcold * pj->sf_data.H2_fraction;
#else
    bi->cold_gas_mass += mj;
#endif
    bi->gas_SFR += max(pj->sf_data.SFR, 0.);
    //if (bi->subgrid_mass * bh_props->mass_to_solar_mass > 1.e10) message("BH_SFR bid=%lld pid=%lld mj=%g psfr=%g nH=%g T=%g fH2=%g totSFR=%g", bi->id, pj->id, mj * bh_props->mass_to_solar_mass, max(pj->sf_data.SFR, 0.) * bh_props->mass_to_solar_mass / bh_props->time_to_yr, pj->rho * bh_props->conv_factor_density_to_cgs / 1.673e-24, pj->u * cosmo->a_factor_internal_energy / (bh_props->T_K_to_int * bh_props->temp_to_u_factor), pj->sf_data.H2_fraction, bi->gas_SFR * bh_props->mass_to_solar_mass / bh_props->time_to_yr);
  }

  const float L_x = mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  const float L_y = mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  const float L_z = mj * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Gas angular momentum in kernel */
  bi->angular_momentum_gas[0] -= L_x;
  bi->angular_momentum_gas[1] -= L_y;
  bi->angular_momentum_gas[2] -= L_z;

  /* Contribution to the velocity (gas w.r.t. black hole) */
  bi->velocity_gas[0] += mj * dv[0];
  bi->velocity_gas[1] += mj * dv[1];
  bi->velocity_gas[2] += mj * dv[2];

  /* Contribution to the specific angular momentum of gas, which is later
   * converted to the circular velocity */
  bi->circular_velocity_gas[0] -= L_x;
  bi->circular_velocity_gas[1] -= L_y;
  bi->circular_velocity_gas[2] -= L_z;

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
    const integertime_t ti_current, const double time) {

  /* Ignore decoupled wind particles for repositioning BH */
  if (pj->decoupled) return;

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
      bh_props->max_reposition_distance_ratio * bi->h * bi->h;

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
    const integertime_t ti_current, const double time) {

  /* Collect information about galaxy that the particle belongs to */
  float galaxy_mstar = bi->galaxy_data.stellar_mass;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (galaxy_mstar <= 0.f) return;

  /* If there is no gas, skip */
  if (bi->ngb_mass <= 0.f) return;

  /* Black hole should only accrete gas within the same galaxy */
  if (galaxy_mstar != pj->galaxy_data.stellar_mass) return;

  /* Gas particle must be bound to BH kernel mass to be swallowed */
  const float r = sqrtf(r2);
  const float dv[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                       bi->v[2] - pj->v[2]};
  const float u_kinetic = 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
  const float u_potential =
      bh_props->const_newton_G * (bi->ngb_mass + bi->mass) / r;

  if (u_kinetic > u_potential) return;

  /* Compute the kernel function */
  float wi;
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Sum up cold disk mass corotating relative to total angular momentum.
   * This can't be used to compute accretion since it is not summable
   * prior to this swallow routine; we need an intermediate loop
   * if we want to use this to set accretion ala Simba.*/
   
  const float mj = hydro_get_mass(pj);

  /* Classify gas as hot or cold for accretion */
  const int gas_temperature_state = black_hole_gas_hot_or_cold(bi, pj, cosmo, bh_props);

  const float Lx = mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  const float Ly = mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  const float Lz = mj * (dx[0] * dv[1] - dx[1] * dv[0]);
  const float proj = Lx * bi->angular_momentum_gas[0] +
                     Ly * bi->angular_momentum_gas[1] +
                     Lz * bi->angular_momentum_gas[2];
  if ((proj > 0.f) && gas_temperature_state == -1) {
    bi->cold_disk_mass += mj;
  }

  /* Probability to swallow this particle */
  float prob = -1.f;
  float f_accretion = bi->f_accretion;

  /* No accretion? Nothing to do */
  if (f_accretion <= 0.f) return;

  const float pj_mass_orig = mj;
  const float nibbled_mass = f_accretion * pj_mass_orig;

  /* Normalize the weights */
  const float kernel_wt =
      (bi->kernel_wt_sum > 0.f) ? wi / bi->kernel_wt_sum : 0.f;

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
  const float mass_deficit = bi->subgrid_mass - bi->mass_at_start_of_step;
  if (mass_deficit > 0.f) {
    /* Don't nibble from particles that are too small already */
    if (mj < bh_props->min_gas_mass_for_nibbling) return;

    /* Just enough to satisfy M_dot,inflow */
    prob = (mass_deficit / f_accretion) * kernel_wt;
  } else {
    /* Do not grow the physical mass, only kick */
    f_accretion = 0.f;

    /* Check the accretion reservoir and if it has passed the limit */
    if (bi->unresolved_mass_reservoir > 0.f) {
      prob = bi->unresolved_mass_reservoir * kernel_wt;
    }
  }

  /* Draw a random number (Note mixing both IDs) */
  const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                          random_number_BH_swallow);
  float new_gas_mass = pj_mass_orig;
  /* Are we lucky? */
  if (rand < prob) {

    if (f_accretion > 0.f) {
      const float bi_mass_orig = bi->mass;
      new_gas_mass = pj_mass_orig - nibbled_mass;
      /* Don't go below the minimum for stability */
      if (new_gas_mass < bh_props->min_gas_mass_for_nibbling) return;

      bi->mass += nibbled_mass;
      hydro_set_mass(pj, new_gas_mass);

      /* Add the angular momentum of the accreted gas to the BH total.
       * Note no change to gas here. The cosmological conversion factors for
       * velocity (a^-1) and distance (a) cancel out, so the angular momentum
       * is already in physical units. */
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

      /* TODO: Spoke to Matthieu about this, it is a bug cannot assign here. */
      bi->v[0] = bi_mom[0] / bi->mass;
      bi->v[1] = bi_mom[1] / bi->mass;
      bi->v[2] = bi_mom[2] / bi->mass;

      /* Update the BH and also gas metal masses */
      struct chemistry_bpart_data *bi_chem = &bi->chemistry_data;
      struct chemistry_part_data *pj_chem = &pj->chemistry_data;
      chemistry_transfer_part_to_bpart(bi_chem, pj_chem, nibbled_mass,
                                       nibbled_mass / pj_mass_orig);
    }

    /* This particle is swallowed by the BH with the largest ID of all the
     * candidates wanting to swallow it */
    if (pj->black_holes_data.swallow_id < bi->id) {
      /* Handle the ADAF heating separately */
      if (bi->state != BH_states_adaf) {
        pj->black_holes_data.swallow_id = bi->id;
      }

      /* Keep track of unresolved mass kicks */
      if (mass_deficit <= 0.f) {
        bi->unresolved_mass_kicked_this_step += new_gas_mass;
      }
    } else {
      message(
          "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
          "swallow id=%lld)",
          bi->id, pj->id, pj->black_holes_data.swallow_id);
    }
  }

  /* When there is zero mass loading the ADAF mode heats the entire kernel
   * so all of the weights are required in the sum. */
  if (bi->state == BH_states_adaf) {
    /* Zero mass loading implies entire
    kernel heating */
    if (bh_props->adaf_wind_mass_loading == 0.f) {
      const float adaf_wt = new_gas_mass * wi;
      bi->adaf_wt_sum += adaf_wt;

      /* Normalized later */
      bi->adaf_energy_used_this_step += bi->adaf_energy_to_dump * adaf_wt;
    } else {
      /* --- The entire kernel is NOT heated here ---- */

      /* Heating is based on specific energy and the mass loading */
      if (bi->adaf_energy_to_dump > 0.f && bh_props->adaf_wind_speed > 0.f) {
        const float adaf_v2 =
            bh_props->adaf_wind_speed * bh_props->adaf_wind_speed;
        const float adaf_mass_to_heat = 2.f * bi->adaf_energy_to_dump / adaf_v2;

        /* Bernoulli trial P = M_adaf * wi / sum(mj * wj) */
        const float adaf_heat_prob = adaf_mass_to_heat * kernel_wt;

        /* Draw a random number (Note mixing both IDs) */
        const float adaf_rand = random_unit_interval(
            bi->id + pj->id, ti_current, random_number_BH_swallow);

        /* Identify ADAF heating particles by their ID, and only heat those! */
        if (adaf_rand < adaf_heat_prob) {
          if (pj->black_holes_data.adaf_id < bi->id) {
            pj->black_holes_data.adaf_id = bi->id;

            /* New normalization for ADAF heating. Use new gas mass since that
             * is what is in the feedback loop. */
            const float adaf_wt = new_gas_mass * wi;
            bi->adaf_wt_sum += adaf_wt;
            /* Will be normalized at the end when reducing the reservoir */
            bi->adaf_energy_used_this_step += bi->adaf_energy_to_dump * adaf_wt;
          } else {
            message(
                "BH %lld wants to heat particle %lld BUT CANNOT (old "
                "adaf_id=%lld)",
                bi->id, pj->id, pj->black_holes_data.adaf_id);
          }
        }
      }
    }
  }

  /* Check jet reservoir regardless of the state, allows simultaneous
   * ADAF heating and a kinetic jet. */
  if (bi->jet_mass_reservoir >= bh_props->jet_minimum_reservoir_mass) {

    /* Make sure there is enough gas to kick */
    if (bi->ngb_mass < bh_props->jet_minimum_reservoir_mass) return;

    float jet_prob = bi->jet_mass_reservoir * kernel_wt;
    const float rand_jet = random_unit_interval(bi->id + pj->id, ti_current,
                                                random_number_BH_kick);

    /* Here the particle is also identified to be kicked out as a jet */
    if (rand_jet < jet_prob) {

      /* If we also are accreting above, the mass loss is already taken
       * into account */

      if (pj->black_holes_data.jet_id < bi->id) {
        bi->jet_mass_kicked_this_step += new_gas_mass;
        pj->black_holes_data.jet_id = bi->id;
      } else {
        message(
            "BH %lld wants to kick jet particle %lld BUT CANNOT (old "
            "jet_id=%lld)",
            bi->id, pj->id, pj->black_holes_data.jet_id);
      }
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
      bh_props->max_reposition_distance_ratio * bi->h * bi->h;

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
      bh_props->max_merging_distance_ratio * bi->h * bi->h;

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

#ifdef OBSIDIAN_DEBUG_CHECKS
        //if (bi->mass * bh_props->mass_to_solar_mass > 1.e9) message("BH_MERGER: z=%g bid=%lld to swallow bid=%lld: MBHi=%g MBHj=%g Mgali=%g Mgalj=%g", cosmo->z, bi->id, bj->id, bi->subgrid_mass * bh_props->mass_to_solar_mass, bj->subgrid_mass * bh_props->mass_to_solar_mass, bi->galaxy_data.stellar_mass * bh_props->mass_to_solar_mass, bj->galaxy_data.stellar_mass * bh_props->mass_to_solar_mass);
#endif

        bj->merger_data.swallow_id = bi->id;
        bj->merger_data.swallow_mass = bi->subgrid_mass;

      } else {

        message(
            "BH %lld wants to swallow bh particle %lld BUT CANNOT (old "
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
    const struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time) {

  /* Gas particle must be bound to BH kernel mass to have feedback
   * (avoids fast-moving decoupled particles) */
  const float r = sqrtf(r2);
  const float dv[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                       bi->v[2] - pj->v[2]};
  const float u_kinetic = 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
  const float u_potential =
      bh_props->const_newton_G * (bi->ngb_mass + bi->mass) / r;

  if (u_kinetic > u_potential) return;

  /* Collect information about galaxy that the particle belongs to */
  const float galaxy_mstar = bi->galaxy_data.stellar_mass;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (galaxy_mstar <= 0.f) return;

  /* A black hole should have gas surrounding it. */
  if (bi->ngb_mass <= 0.f) return;

  /* Need time-step for decoupling and ADAF heating */
  double dt;
  if (with_cosmology) {
    const integertime_t ti_step = get_integer_timestep(bi->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(ti_current - 1, bi->time_bin);

    dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
  } else {
    error("Kiara BH model can only be run with cosmology.");
    dt = 0.;
  }

  /* Save gas density and entropy before feedback */
  tracers_before_black_holes_feedback(pj, xpj, cosmo->a);

  float v_kick = bi->v_kick; /* PHYSICAL */

  const float bh_mass_msun = bi->subgrid_mass * bh_props->mass_to_solar_mass;

  /* by-eye fit from Fig 6 of Zhang+2023 (2305.06803) */
  double halo_mass = 1.e12 * pow(bh_mass_msun * 1.e-7, 0.75);
  if (halo_mass < 6.3e11) halo_mass = 6.3e11;

  /* In internal temperature units */
  const double T_vir = 9.52e7 * pow(halo_mass * 1.e-15, 0.6666) *
                       bh_props->T_K_to_int * (1. + cosmo->z);

  /* In the swallow loop the particle was marked as a kick particle */
  const int swallow_flag = (pj->black_holes_data.swallow_id == bi->id);

  /* Can get reset if there is a adaf_kick_factor > 0 and heating */
  int adaf_kick_flag = 0;

  /* Initialize heating */
  const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
  double E_heat = 0.;
  double E_inject = 0.;
  double u_new = u_init;
  double T_new = u_new / bh_props->temp_to_u_factor;

  int adaf_energy_flag =
      (bi->state == BH_states_adaf && bi->adaf_energy_to_dump > 0.f);
  int adaf_heat_flag =
      (bi->state == BH_states_adaf && pj->black_holes_data.adaf_id == bi->id);

  /* In the case of non-zero mass loading, require only certain particles
   * to be heated to satisfy M_dot,ADAF  = psi_ADAF * M_dot,acc */
  if (bh_props->adaf_wind_mass_loading > 0.f) {
    adaf_heat_flag = (adaf_energy_flag && adaf_heat_flag);
  } else {
    /* In this case, all of the kernel is heated with adaf_energy_to_dump */
    adaf_heat_flag = adaf_energy_flag;
  }

  /* In the swallow loop the particle was marked as a jet particle */
  int jet_flag = (pj->black_holes_data.jet_id == bi->id);

  /* Compute ramp-up in energy above ADAF mass limit */
  float jet_ramp = black_hole_compute_jet_energy_ramp(bi, cosmo, bh_props);

  /* ADAF heating: Only heat this particle if it is NOT a jet particle */
  if (adaf_heat_flag && !jet_flag) {

    /* compute kernel weights */
    float wj;
    kernel_eval(sqrtf(r2) / hi, &wj);
    const float mj = hydro_get_mass(pj);

    /* Below is equivalent to
     * E_inject_i = E_ADAF * (w_j * m_j) / Sum(w_i * mi) */
    E_inject = bi->adaf_energy_to_dump * mj * wj / bi->adaf_wt_sum;

    /* Initialise heat energy to injection energy */
    E_heat = E_inject;

    /* Heat and/or kick the particle */
    if (E_inject > 0.f) {

      const double n_H_cgs =
          hydro_get_physical_density(pj, cosmo) * bh_props->rho_to_n_cgs;
      const double T_gas_cgs =
          u_init / (bh_props->temp_to_u_factor * bh_props->T_K_to_int);
      const double T_EoS_cgs =
          entropy_floor_temperature(pj, cosmo, floor_props) /
          bh_props->T_K_to_int;

      /* Check whether we are close to the entropy floor or SF/ing. If we are,
       * we classify the gas as cold regardless of temperature. */
      if ((n_H_cgs > bh_props->adaf_heating_n_H_threshold_cgs &&
           (T_gas_cgs < bh_props->adaf_heating_T_threshold_cgs ||
            T_gas_cgs < T_EoS_cgs * bh_props->fixed_T_above_EoS_factor)) ||
          pj->sf_data.SFR > 0.f) {

        /* Kick with some fraction of the energy, if desired */
        if (bh_props->adaf_kick_factor > 0.f) {

          /* Compute kick velocity */
          double E_kick = bh_props->adaf_kick_factor * E_inject;
          v_kick = sqrt(2. * E_kick / mj);

          /* Apply ramp-up in kick velocity above ADAF mass limit */
          const float adaf_max_speed =
              bh_props->adaf_wind_speed * sqrtf(jet_ramp);

          /* Limit kick energy if velocity exceeds max */
          if (v_kick > adaf_max_speed) {
            v_kick = adaf_max_speed;
            E_kick = 0.5 * mj * v_kick * v_kick;
          }

          /* Reduce energy available to heat */
          E_heat = E_inject - E_kick;

          /* Later will apply velocity as if it was flagged to swallow */
          adaf_kick_flag = 1;

        } /* adaf_kick_factor > 0 */

      } /* If in ISM */

      /* Heat gas with remaining energy, if any */
      if (E_heat > 0.) {

        /* Compute new energy per unit mass of this particle */
        u_new = u_init + E_heat / mj;

        /* New temperature */
        T_new = u_new / bh_props->temp_to_u_factor;

        /* Limit heating.  There can sometimes be VERY large amounts of
         * energy to deposit */
        if (bh_props->adaf_maximum_temperature > 0.f) {
          if (T_new > bh_props->adaf_maximum_temperature) {
            u_new =
                bh_props->adaf_maximum_temperature * bh_props->temp_to_u_factor;
          }
        } else {
          const float T_max = fabs(bh_props->adaf_maximum_temperature) * T_vir;
          if (T_new > T_max) {
            u_new = T_max * bh_props->temp_to_u_factor;
            T_new = T_max;
          }
        }

        /* Reset in case clipped at the upper limit */
        E_heat = (u_new - u_init) * mj;

        /* Heat particle: We are overwriting the internal energy of the
         * particle */
        hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
        hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

#if COOLING_GRACKLE_MODE >= 2
	/* Take AGN-heated gas out of subgrid ISM mode */
	pj->cooling_data.subgrid_temp = 0.f;
	pj->cooling_data.subgrid_dens = hydro_get_physical_density(pj, cosmo);
	pj->cooling_data.subgrid_fcold = 0.f;
#endif

        /* Shut off cooling for some time, if desired */
        if (bh_props->adaf_cooling_shutoff_factor > 0.f) {

          /* u_init is physical so cs_physical is physical */
          const double u_com = u_new / cosmo->a_factor_internal_energy;
          const double cs = gas_soundspeed_from_internal_energy(pj->rho, u_com);

          const float h_phys = kernel_gamma * pj->h * cosmo->a;
          const float cs_physical = cs * cosmo->a_factor_sound_speed;
          const float dt_sound_phys = h_phys / cs_physical;

          /* a_factor_sound_speed converts cs_physical to comoving units,
           * twice the BH timestep as a lower limit */
          pj->feedback_data.cooling_shutoff_delay_time =
              bh_props->adaf_cooling_shutoff_factor * min(dt_sound_phys, dt);
        }

      } /* E_heat > 0 */

    } /* E_inject > 0 */

  } /* non-jet ADAF mode */

  /* ----- If particle is marked as a jet, do jet feedback ----- */

  /* Heat the particle and set kinetic kick information if jet particle */
  if (jet_flag) {

    /* Set jet velocity, accounting for energy ramp-up */
    v_kick = black_hole_compute_jet_velocity(bi, cosmo, bh_props);
    v_kick *= sqrtf(jet_ramp);

    /* Heat jet particle */
    float new_Tj = bh_props->jet_temperature;

    /* Use the halo T_vir? */
    if (bh_props->jet_temperature < 0.f) {
      new_Tj = fabs(bh_props->jet_temperature) * T_vir;
    }

    /* Compute new energy per unit mass of this particle */
    u_new = new_Tj * bh_props->temp_to_u_factor;

    /* Only increase the gas temperature if it's below the target T */
    if (u_new > u_init) {
      /* account for energy ramp-up */
      u_new = jet_ramp * (u_new - u_init) + u_init;
      /* We are overwriting the internal energy of the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

      const double delta_energy = (u_new - u_init) * hydro_get_mass(pj);
      E_heat += delta_energy;

      tracers_after_black_holes_feedback(pj, xpj, with_cosmology, cosmo->a,
                                         time, delta_energy);
    }

  } /* jet_flag */

#ifdef OBSIDIAN_DEBUG_CHECKS
  float pj_vel_norm = FLT_MAX;
#endif

  /* Flagged if it is a jet particle, marked to swallow (i.e. kick) or
   * if there was an ADAF kick because of energy splitting. */
  int flagged_to_kick =
      jet_flag || (swallow_flag && !adaf_heat_flag) || adaf_kick_flag;

  /* Kick the particle if is was tagged only */
  if (v_kick > 0.f && flagged_to_kick) {

    /* Set direction of launch: 0=random, 1=L_gas, 2=L_BH, 3=outwards */
    float dir[3] = {0.f, 0.f, 0.f};
    int dir_flag = 0;
    if (jet_flag) {
      dir_flag = bh_props->jet_launch_dir;
    } else if (adaf_heat_flag) {
      dir_flag = bh_props->adaf_wind_dir;
    } else if (bi->state == BH_states_quasar) {
      dir_flag = bh_props->quasar_wind_dir;
      if (bi->radiative_luminosity > bh_props->quasar_luminosity_thresh &&
          bh_props->quasar_luminosity_thresh > 0.f) {
        dir_flag = 3;  // outwards blowout above threshold luminosity
      }
    } else if (bi->state == BH_states_slim_disk) {
      dir_flag = bh_props->slim_disk_wind_dir;
    } else {
      warning(
          "Cannot determine wind direction (BH state=%d) for v_kick=%g, "
          "setting to random",
          bi->state, v_kick / bh_props->kms_to_internal);
    }

    float dirsign =
        black_hole_set_kick_direction(bi, pj, ti_current, dir_flag, dir);

    /* Do the kick */
    const float norm =
        sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);

    if (norm > 0.f) {
      const float prefactor = v_kick * cosmo->a * dirsign / norm;

#ifdef OBSIDIAN_DEBUG_CHECKS
      pj_vel_norm = sqrtf(xpj->v_full[0] * xpj->v_full[0] +
                          xpj->v_full[1] * xpj->v_full[1] +
                          xpj->v_full[2] * xpj->v_full[2]);
#endif

      xpj->v_full[0] += prefactor * dir[0];
      xpj->v_full[1] += prefactor * dir[1];
      xpj->v_full[2] += prefactor * dir[2];

#ifdef OBSIDIAN_DEBUG_CHECKS
      const float v_mag = sqrtf(xpj->v_full[0] * xpj->v_full[0] +
                                xpj->v_full[1] * xpj->v_full[1] +
                                xpj->v_full[2] * xpj->v_full[2]);

      if (prefactor * norm > 1.e3 * v_mag) {
        warning(
            "LARGE KICK! z=%g id=%lld dv=%g vkick=%g vadaf=%g vjet=%g v=%g "
            "(%g,%g,%g) dir=%g,%g,%g",
            cosmo->z, pj->id, prefactor, v_kick, bh_props->adaf_wind_speed,
            bh_props->jet_velocity, v_mag, xpj->v_full[0], xpj->v_full[1],
            xpj->v_full[2], dir[0], dir[1], dir[2]);
      }
#endif

      /* Update the signal velocity of the particle based
       * on the PHYSICAL velocity kick. */
      hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, v_kick);
      pj->chemistry_data.diffusion_coefficient = 0.f;

      float f_decouple = 0.f;
      switch (bi->state) {
        case BH_states_adaf:
          f_decouple = bh_props->adaf_decouple_time_factor;
          break;
        case BH_states_quasar:
          f_decouple = bh_props->quasar_decouple_time_factor;
          break;
        case BH_states_slim_disk:
          f_decouple = bh_props->slim_disk_decouple_time_factor;
          break;
      }

      /* Hubble time in internal units */
      const double t_H = cosmology_get_time_since_big_bang(cosmo, cosmo->a);

      /* Set delay time to at least the time-step*/
      pj->feedback_data.decoupling_delay_time = dt + f_decouple * t_H;
      pj->decoupled = 1;

      /* Count number of decouplings */
      if (jet_flag) {
        if (bh_props->jet_decouple_time_factor > 0.f) {
          pj->feedback_data.decoupling_delay_time =
              dt + bh_props->jet_decouple_time_factor * t_H;
        } else {
          pj->feedback_data.decoupling_delay_time = 0.f;
        }
        pj->feedback_data.number_of_times_decoupled += 100000;
      } else {
        if ((bh_props->slim_disk_decouple_time_factor > 0.f &&
             bi->state == BH_states_slim_disk) ||
            (bh_props->quasar_decouple_time_factor > 0.f &&
             bi->state == BH_states_quasar)) {
          pj->feedback_data.decoupling_delay_time =
              dt + bh_props->slim_disk_decouple_time_factor * t_H;
        } else {
          pj->feedback_data.decoupling_delay_time = 0.f;
        }
        pj->feedback_data.number_of_times_decoupled += 1000;
      }
    } else {
      v_kick = 0.f;
    }
  }

  /* This particle was touched by BH feedback, so reset some variables */
  if ((v_kick > 0.f && flagged_to_kick) || E_heat > 0.f) {
    /* set SFR=0 for BH feedback particle */
    if (pj->sf_data.SFR > 0.f) {
      /* Record the current time as an indicator of when this particle was last
        star-forming. */
      if (with_cosmology) {
        pj->sf_data.SFR = -cosmo->a;
      } else {
        pj->sf_data.SFR = -time;
      }
    }

    /* Destroy all H2 and put into HI */
    xpj->cooling_data.HI_frac += xpj->cooling_data.HM_frac +
                                 xpj->cooling_data.H2I_frac +
                                 xpj->cooling_data.H2II_frac;
    xpj->cooling_data.HM_frac = 0.f;
    xpj->cooling_data.H2I_frac = 0.f;
    xpj->cooling_data.H2II_frac = 0.f;

    /* Only take it out of ISM mode if it was kicked */
    if (v_kick > 0.f && flagged_to_kick) {
      /* Take particle out of subgrid ISM mode */
      pj->cooling_data.subgrid_temp = 0.f;
      pj->cooling_data.subgrid_dens = hydro_get_physical_density(pj, cosmo);
      pj->cooling_data.subgrid_fcold = 0.f;
    }

    /* Destroy all dust in ADAF-"touched" gas and the jet */
    if (jet_flag || E_inject > 0.) {
#if COOLING_GRACKLE_MODE >= 2
      const float old_dust_mass = pj->cooling_data.dust_mass;
      pj->cooling_data.dust_mass = 0.f;
      float new_Z_total = 0.f;
      pj->chemistry_data.metal_mass_fraction_total = 0.f;
      for (int elem = chemistry_element_He; elem < chemistry_element_count;
           ++elem) {
        const float old_metal_mass_elem =
            pj->chemistry_data.metal_mass_fraction[elem] * hydro_get_mass(pj);
        const float old_dust_mass_elem =
            pj->cooling_data.dust_mass_fraction[elem] * old_dust_mass;

        pj->chemistry_data.metal_mass_fraction[elem] =
            (old_metal_mass_elem + old_dust_mass_elem) / hydro_get_mass(pj);

        if (elem != chemistry_element_H && elem != chemistry_element_He) {
          new_Z_total += pj->chemistry_data.metal_mass_fraction[elem];
        }

        pj->cooling_data.dust_mass_fraction[elem] = 0.f;
      }

      pj->chemistry_data.metal_mass_fraction_total = new_Z_total;
#endif
    }

    /* Impose maximal viscosity */
    hydro_diffusive_feedback_reset(pj);

    /* Synchronize the particle on the timeline */
    timestep_sync_part(pj);

#ifdef OBSIDIAN_DEBUG_CHECKS
    if (E_heat > 0.f) {
      message(
          "BH_HEAT_ADAF: z=%g bid=%lld pid=%lld mbh=%g Msun u=%g T=%g K "
          "Tvir=%g K",
          cosmo->z, bi->id, pj->id, bh_mass_msun, pj->u,
          T_new / bh_props->T_K_to_int, T_vir / bh_props->T_K_to_int);
    }

    switch (bi->state) {
      case BH_states_quasar:
        message(
            "BH_KICK_QSO: z=%g bid=%lld mbh=%g Msun v_kick=%g km/s "
            "v_kick/v_part=%g T=%g K",
            cosmo->z, bi->id, bh_mass_msun, v_kick / bh_props->kms_to_internal,
            v_kick * cosmo->a / pj_vel_norm,
            hydro_get_physical_internal_energy(pj, xpj, cosmo) /
                (bh_props->T_K_to_int * bh_props->temp_to_u_factor));
        break;
      case BH_states_slim_disk:
        message("BH_KICK_SLIM: z=%g bid=%lld mbh=%g Msun v_kick=%g km/s T=%g K",
                cosmo->z, bi->id, bh_mass_msun,
                v_kick / bh_props->kms_to_internal,
                hydro_get_physical_internal_energy(pj, xpj, cosmo) /
                    (bh_props->T_K_to_int * bh_props->temp_to_u_factor));
        break;
      case BH_states_adaf:
        if (jet_flag) {
          message(
              "BH_KICK_JET: z=%g bid=%lld mbh=%g Msun v_kick=%g km/s "
              "v_kick/v_part=%g T=%g",
              cosmo->z, bi->id, bh_mass_msun,
              v_kick / bh_props->kms_to_internal,
              v_kick * cosmo->a / pj_vel_norm,
              hydro_get_physical_internal_energy(pj, xpj, cosmo) /
                  (bh_props->T_K_to_int * bh_props->temp_to_u_factor));
        } else {
          message(
              "BH_KICK_ADAF: z=%g bid=%lld pid=%lld mbh=%g Msun "
              "v_kick=%g km/s "
              "v_kick/v_part=%g u=%g T=%g",
              cosmo->z, bi->id, pj->id, bh_mass_msun,
              v_kick / bh_props->kms_to_internal,
              v_kick * cosmo->a / pj_vel_norm, pj->u,
              hydro_get_physical_internal_energy(pj, xpj, cosmo) /
                  (bh_props->T_K_to_int * bh_props->temp_to_u_factor));
        }
        break;
    }
#endif
  }

  if (swallow_flag) {
    /* IMPORTANT: The particle MUST NOT be swallowed.
     * We are taking a f_accretion from each particle, and then
     * kicking the rest. We used the swallow marker as a temporary
     * passer in order to remember which particles have been "nibbled"
     * so that we can kick them out.
     */
    black_holes_mark_part_as_not_swallowed(&pj->black_holes_data);
  }
}

#endif /* SWIFT_OBSIDIAN_BH_IACT_H */
