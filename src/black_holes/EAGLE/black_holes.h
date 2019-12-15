/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_BLACK_HOLES_H
#define SWIFT_EAGLE_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "cosmology.h"
#include "dimension.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"

/* Standard includes */
#include <float.h>

/**
 * @brief Computes the gravity time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart* const bp) {

  return FLT_MAX;
}

/**
 * @brief Initialises the b-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param bp The particle to act upon
 * @param props The properties of the black holes model.
 */
__attribute__((always_inline)) INLINE static void black_holes_first_init_bpart(
    struct bpart* bp, const struct black_holes_props* props) {

  bp->time_bin = 0;
  bp->subgrid_mass = bp->mass;
  bp->total_accreted_mass = 0.f;
  bp->accretion_rate = 0.f;
  bp->formation_time = -1.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->last_high_Eddington_fraction_scale_factor = -1.f;
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart* bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_density[i] = -1;
  bp->num_ngb_density = 0;
#endif

  bp->density.wcount = 0.f;
  bp->density.wcount_dh = 0.f;
  bp->rho_gas = 0.f;
  bp->sound_speed_gas = 0.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->circular_velocity_gas[0] = 0.f;
  bp->circular_velocity_gas[1] = 0.f;
  bp->circular_velocity_gas[2] = 0.f;
  bp->ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * The fields do not get predicted but we move the BH to its new position
 * if a new one was calculated in the repositioning loop.
 *
 * @param bp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void black_holes_predict_extra(
    struct bpart* restrict bp, float dt_drift) {

  /* Are we doing some repositioning? */
  if (bp->reposition.min_potential != FLT_MAX) {

#ifdef SWIFT_DEBUG_CHECKS
    if (bp->reposition.delta_x[0] == -FLT_MAX ||
        bp->reposition.delta_x[1] == -FLT_MAX ||
        bp->reposition.delta_x[2] == -FLT_MAX) {
      error("Something went wrong with the new repositioning position");
    }

    const double dx = bp->reposition.delta_x[0];
    const double dy = bp->reposition.delta_x[1];
    const double dz = bp->reposition.delta_x[2];
    const double d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d > 1.01 * kernel_gamma * bp->h)
      error("Repositioning BH beyond the kernel support!");
#endif

    /* Move the black hole */
    bp->x[0] += bp->reposition.delta_x[0];
    bp->x[1] += bp->reposition.delta_x[1];
    bp->x[2] += bp->reposition.delta_x[2];

    /* Move its gravity properties as well */
    bp->gpart->x[0] += bp->reposition.delta_x[0];
    bp->gpart->x[1] += bp->reposition.delta_x[1];
    bp->gpart->x[2] += bp->reposition.delta_x[2];

    /* Store the delta position */
    bp->x_diff[0] -= bp->reposition.delta_x[0];
    bp->x_diff[1] -= bp->reposition.delta_x[1];
    bp->x_diff[2] -= bp->reposition.delta_x[2];

    /* Reset the reposition variables */
    bp->reposition.delta_x[0] = -FLT_MAX;
    bp->reposition.delta_x[1] = -FLT_MAX;
    bp->reposition.delta_x[2] = -FLT_MAX;
    bp->reposition.min_potential = FLT_MAX;
  }
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param bp The particle.
 */
__attribute__((always_inline)) INLINE static void
black_holes_reset_predicted_values(struct bpart* bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart* bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart* bp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  bp->sound_speed_gas *= h_inv_dim;
  bp->velocity_gas[0] *= h_inv_dim;
  bp->velocity_gas[1] *= h_inv_dim;
  bp->velocity_gas[2] *= h_inv_dim;
  bp->circular_velocity_gas[0] *= h_inv_dim;
  bp->circular_velocity_gas[1] *= h_inv_dim;
  bp->circular_velocity_gas[2] *= h_inv_dim;

  const float rho_inv = 1.f / bp->rho_gas;

  /* Finish the calculation by undoing the mass smoothing */
  bp->sound_speed_gas *= rho_inv;
  bp->velocity_gas[0] *= rho_inv;
  bp->velocity_gas[1] *= rho_inv;
  bp->velocity_gas[2] *= rho_inv;
  bp->circular_velocity_gas[0] *= rho_inv * h_inv;
  bp->circular_velocity_gas[1] *= rho_inv * h_inv;
  bp->circular_velocity_gas[2] *= rho_inv * h_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_bpart_has_no_neighbours(struct bpart* bp,
                                    const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a gas particle.
 *
 * @param bp The #bpart to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_part(
    struct bpart* bp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float BH_mass = bp->mass;

  /* Increase the dynamical mass of the BH. */
  bp->mass += gas_mass;
  bp->gpart->mass += gas_mass;

  /* Physical velocity difference between the particles */
  const float dv[3] = {(bp->v[0] - p->v[0]) * cosmo->a_inv,
                       (bp->v[1] - p->v[1]) * cosmo->a_inv,
                       (bp->v[2] - p->v[2]) * cosmo->a_inv};

  /* Physical distance between the particles */
  const float dx[3] = {(bp->x[0] - p->x[0]) * cosmo->a,
                       (bp->x[0] - p->x[0]) * cosmo->a,
                       (bp->x[0] - p->x[0]) * cosmo->a};

  /* Collect the swallowed angular momentum */
  bp->swallowed_angular_momentum[0] +=
      gas_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
  bp->swallowed_angular_momentum[1] +=
      gas_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
  bp->swallowed_angular_momentum[2] +=
      gas_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Update the BH momentum */
  const float BH_mom[3] = {BH_mass * bp->v[0] + gas_mass * p->v[0],
                           BH_mass * bp->v[1] + gas_mass * p->v[1],
                           BH_mass * bp->v[2] + gas_mass * p->v[2]};

  bp->v[0] = BH_mom[0] / bp->mass;
  bp->v[1] = BH_mom[1] / bp->mass;
  bp->v[2] = BH_mom[2] / bp->mass;
  bp->gpart->v_full[0] = bp->v[0];
  bp->gpart->v_full[1] = bp->v[1];
  bp->gpart->v_full[2] = bp->v[2];

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_add_part_to_bpart(bp_chem, p_chem, gas_mass);

  /* This BH lost a neighbour */
  bp->num_ngbs--;
  bp->ngb_mass -= gas_mass;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a BH particle.
 *
 * @param bpi The #bpart to update.
 * @param bpj The #bpart that is swallowed.
 * @param cosmo The current cosmological model.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param props The properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_bpart(
    struct bpart* bpi, const struct bpart* bpj, const struct cosmology* cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props* props) {

  /* Get the current dynamical masses */
  const float bpi_dyn_mass = bpi->mass;
  const float bpj_dyn_mass = bpj->mass;

  /* Is this merger ratio above the threshold for recording? */
  const double merger_ratio = bpj->subgrid_mass / bpi->subgrid_mass;
  if (merger_ratio > props->major_merger_threshold) {
    if (with_cosmology) {
      bpi->last_major_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_major_merger_time = time;
    }
  } else if (merger_ratio > props->minor_merger_threshold) {
    if (with_cosmology) {
      bpi->last_minor_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_minor_merger_time = time;
    }
  }

  /* Increase the masses of the BH. */
  bpi->mass += bpj->mass;
  bpi->gpart->mass += bpj->mass;
  bpi->subgrid_mass += bpj->subgrid_mass;

  /* Update the BH momentum */
  const float BH_mom[3] = {bpi_dyn_mass * bpi->v[0] + bpj_dyn_mass * bpj->v[0],
                           bpi_dyn_mass * bpi->v[1] + bpj_dyn_mass * bpj->v[1],
                           bpi_dyn_mass * bpi->v[2] + bpj_dyn_mass * bpj->v[2]};

  bpi->v[0] = BH_mom[0] / bpi->mass;
  bpi->v[1] = BH_mom[1] / bpi->mass;
  bpi->v[2] = BH_mom[2] / bpi->mass;
  bpi->gpart->v_full[0] = bpi->v[0];
  bpi->gpart->v_full[1] = bpi->v[1];
  bpi->gpart->v_full[2] = bpi->v[2];

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bpi_chem = &bpi->chemistry_data;
  const struct chemistry_bpart_data* bpj_chem = &bpj->chemistry_data;
  chemistry_add_bpart_to_bpart(bpi_chem, bpj_chem);

  /* Update the energy reservoir */
  bpi->energy_reservoir += bpj->energy_reservoir;

  /* Add up all the BH seeds */
  bpi->cumulative_number_seeds += bpj->cumulative_number_seeds;

  /* We had another merger */
  bpi->number_of_mergers++;
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantites
 * required for the feedback loop.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const double time, const int with_cosmology, const double dt) {

  if (dt == 0.) return;

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd = props->f_Edd;
  const double f_Edd_recording = props->f_Edd_recording;
  const double epsilon_r = props->epsilon_r;
  const double epsilon_f = props->epsilon_f;
  const double num_ngbs_to_heat = props->num_ngbs_to_heat;
  const double delta_T = props->AGN_delta_T_desired;
  const double delta_u = delta_T * props->temp_to_u_factor;
  const double alpha_visc = props->alpha_visc;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (all internal units)
   * Note: for the velocities this means peculiar velocities */
  const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
  const double gas_c_phys = bp->sound_speed_gas * cosmo->a_factor_sound_speed;
  const double gas_v_peculiar[3] = {bp->velocity_gas[0] * cosmo->a_inv,
                                    bp->velocity_gas[1] * cosmo->a_inv,
                                    bp->velocity_gas[2] * cosmo->a_inv};

  const double bh_v_peculiar[3] = {bp->v[0] * cosmo->a_inv,
                                   bp->v[1] * cosmo->a_inv,
                                   bp->v[2] * cosmo->a_inv};

  const double gas_v_circular[3] = {
      bp->circular_velocity_gas[0] * cosmo->a_inv,
      bp->circular_velocity_gas[1] * cosmo->a_inv,
      bp->circular_velocity_gas[2] * cosmo->a_inv};

  /* Difference in peculiar velocity between the gas and the BH
   * Note that there is no need for a Hubble flow term here. We are
   * computing the gas velocity at the position of the black hole. */
  const double v_diff_peculiar[3] = {gas_v_peculiar[0] - bh_v_peculiar[0],
                                     gas_v_peculiar[1] - bh_v_peculiar[1],
                                     gas_v_peculiar[2] - bh_v_peculiar[2]};
  const double v_diff_norm2 = v_diff_peculiar[0] * v_diff_peculiar[0] +
                              v_diff_peculiar[1] * v_diff_peculiar[1] +
                              v_diff_peculiar[2] * v_diff_peculiar[2];

  /* Norm of the circular velocity of the gas around the BH */
  const double tangential_velocity2 = gas_v_circular[0] * gas_v_circular[0] +
                                      gas_v_circular[1] * gas_v_circular[1] +
                                      gas_v_circular[2] * gas_v_circular[2];
  const double tangential_velocity = sqrt(tangential_velocity2);

  /* We can now compute the Bondi accretion rate (internal units) */
  const double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double denominator2 = v_diff_norm2 + gas_c_phys2;
  const double denominator_inv = 1. / sqrt(denominator2);
  double Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                      denominator_inv * denominator_inv * denominator_inv;

  /* Compute the reduction factor from Rosas-Guevara et al. (2015) */
  const double Bondi_radius = G * BH_mass / gas_c_phys2;
  const double Bondi_time = Bondi_radius / gas_c_phys;
  const double r_times_v_tang = Bondi_radius * tangential_velocity;
  const double r_times_v_tang_3 =
      r_times_v_tang * r_times_v_tang * r_times_v_tang;
  const double viscous_time = 2. * M_PI * r_times_v_tang_3 /
                              (1e-6 * alpha_visc * G * G * BH_mass * BH_mass);
  const double f_visc = max(Bondi_time / viscous_time, 1.);

  /* Limit the Bondi rate by the Bondi viscuous time ratio */
  Bondi_rate *= f_visc;

  /* Compute the Eddington rate (internal units) */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (epsilon_r * c * sigma_Thomson);

  /* Shall we record this high rate? */
  if (Bondi_rate > f_Edd_recording * Eddington_rate) {
    if (with_cosmology) {
      bp->last_high_Eddington_fraction_scale_factor = cosmo->a;
    } else {
      bp->last_high_Eddington_fraction_time = time;
    }
  }

  /* Limit the accretion rate to the Eddington fraction */
  const double accr_rate = min(Bondi_rate, f_Edd * Eddington_rate);
  bp->accretion_rate = accr_rate;

  /* Factor in the radiative efficiency */
  const double mass_rate = (1. - epsilon_r) * accr_rate;
  const double luminosity = epsilon_r * accr_rate * c * c;

  /* Integrate forward in time */
  bp->subgrid_mass += mass_rate * dt;
  bp->total_accreted_mass += mass_rate * dt;
  bp->energy_reservoir += luminosity * epsilon_f * dt;

  /* Energy required to have a feedback event
   * Note that we have subtracted the particles we swallowed from the ngb_mass
   * and num_ngbs accumulators. */
  const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);
  const double E_feedback_event = num_ngbs_to_heat * delta_u * mean_ngb_mass;

  /* Are we doing some feedback? */
  if (bp->energy_reservoir > E_feedback_event) {

    /* Default probability of heating */
    double target_prob = bp->energy_reservoir / (delta_u * bp->ngb_mass);

    /* Calculate the change in internal energy of the gas particles that get
     * heated. Adjust the prbability if needed. */
    double gas_delta_u;
    double prob;
    if (target_prob <= 1.) {

      /* Normal case */
      prob = target_prob;
      gas_delta_u = delta_u;

    } else {

      /* Special case: we need to adjust the energy irrespective of the
       * desired deltaT to ensure we inject all the available energy. */

      prob = 1.;
      gas_delta_u = bp->energy_reservoir / bp->ngb_mass;
    }

    /* Store all of this in the black hole for delivery onto the gas. */
    bp->to_distribute.AGN_heating_probability = prob;
    bp->to_distribute.AGN_delta_u = gas_delta_u;

    /* Decrement the energy in the reservoir by the mean expected energy */
    const double energy_used = bp->energy_reservoir / max(prob, 1.);
    bp->energy_reservoir -= energy_used;

  } else {

    /* Flag that we don't want to heat anyone */
    bp->to_distribute.AGN_heating_probability = 0.f;
    bp->to_distribute.AGN_delta_u = 0.f;
  }
}

/**
 * @brief Finish the calculation of the new BH position.
 *
 * Here, we check that the BH should indeed be moved in the next drift.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo) {

  const float potential = gravity_get_comoving_potential(bp->gpart);

  /* Is the potential lower (i.e. the BH is at the bottom already)
   * OR is the BH massive enough that we don't reposition? */
  if (potential < bp->reposition.min_potential ||
      bp->subgrid_mass > props->max_reposition_mass) {

    /* No need to reposition */
    bp->reposition.min_potential = FLT_MAX;
    bp->reposition.delta_x[0] = -FLT_MAX;
    bp->reposition.delta_x[1] = -FLT_MAX;
    bp->reposition.delta_x[2] = -FLT_MAX;
  }
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on black hole, therefore no need to use
 * it.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_reset_feedback(
    struct bpart* restrict bp) {

  bp->to_distribute.AGN_heating_probability = 0.f;
  bp->to_distribute.AGN_delta_u = 0.f;

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart* bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct part* p) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time */
  bp->formation_scale_factor = cosmo->a;

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

  /* We haven't accreted anything yet */
  bp->total_accreted_mass = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;

  /* Initial metal masses */
  const float gas_mass = hydro_get_mass(p);
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_bpart_from_part(bp_chem, p_chem, gas_mass);

  /* No swallowed angular momentum */
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;

  /* Last time this BH had a high Eddington fraction */
  bp->last_high_Eddington_fraction_scale_factor = -1.f;

  /* Last time of mergers */
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;

  /* First initialisation */
  black_holes_init_bpart(bp);

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

#endif /* SWIFT_EAGLE_BLACK_HOLES_H */
