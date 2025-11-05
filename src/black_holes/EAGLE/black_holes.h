/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "cooling.h"
#include "cooling/PS2020/cooling_tables.h"
#include "cosmology.h"
#include "dimension.h"
#include "exp10.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"
#include "rays.h"

/* Standard includes */
#include <float.h>
#include <math.h>

/**
 * @brief Computes the time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart* const bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo) {

  /* Gather some physical constants (in internal units) */
  const double c = constants->const_speed_light_c;

  /* Compute instantaneous energy supply rate to the BH energy reservoir
   * which is proportional to the BH mass accretion rate */
  const double Energy_rate =
      props->epsilon_r * props->epsilon_f * bp->accretion_rate * c * c;

  /* BHs not accreting can't estimate their time-step length based on the
   * the time-scale for heating */
  if (Energy_rate == 0.) {
    return FLT_MAX;
  }

  /* Average particle mass in BH's kernel */
  const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);

  /* Get the AGN heating temperature that is tored in this BH */
  const double AGN_delta_T = bp->AGN_delta_T;

  /* Without multiplying by mean_ngb_mass we'd get energy per unit mass */
  const double E_heat = AGN_delta_T * props->temp_to_u_factor * mean_ngb_mass;

  /* Compute average time between heating events for the given accretion
   * rate. The time is multiplied by the number of Ngbs to heat because
   * if more particles are heated at once then the time between different
   * AGN feedback events increases proportionally. */
  const double dt_heat = E_heat * bp->num_ngbs_to_heat / Energy_rate;

  /* The new timestep of the BH cannot be smaller than the miminum allowed
   * time-step */
  const double bh_timestep = max(dt_heat, props->time_step_min);

  return bh_timestep;
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
  if (props->use_subgrid_mass_from_ics == 0) {
    bp->subgrid_mass = bp->mass;
  } else if (props->with_subgrid_mass_check && bp->subgrid_mass <= 0) {
    error(
        "Black hole %lld has a subgrid mass of %f (internal units).\n"
        "If this is because the ICs do not contain a 'SubgridMass' data "
        "set, you should set the parameter "
        "'EAGLEAGN:use_subgrid_mass_from_ics' to 0 to initialize the "
        "black hole subgrid masses to the corresponding dynamical masses.\n"
        "If the subgrid mass is intentionally set to this value, you can "
        "disable this error by setting 'EAGLEAGN:with_subgrid_mass_check' "
        "to 0.",
        bp->id, bp->subgrid_mass);
  }
  bp->total_accreted_mass = 0.f;
  bp->accretion_rate = 0.f;
  bp->formation_time = -1.f;
  bp->energy_reservoir = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->number_of_time_steps = 0;
  bp->last_high_Eddington_fraction_scale_factor = -1.f;
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;
  bp->accreted_angular_momentum[0] = 0.f;
  bp->accreted_angular_momentum[1] = 0.f;
  bp->accreted_angular_momentum[2] = 0.f;
  bp->last_repos_vel = 0.f;
  bp->num_ngbs_to_heat = props->num_ngbs_to_heat; /* Filler value */
  bp->dt_heat = 0.f;
  bp->AGN_number_of_AGN_events = 0;
  bp->AGN_number_of_energy_injections = 0;
  bp->AGN_cumulative_energy = 0.f;

  /* Set the initial targetted heating temperature, used for the
   * BH time step determination */
  bp->AGN_delta_T = props->AGN_delta_T_desired;
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
  bp->internal_energy_gas = 0.f;
  bp->rho_subgrid_gas = -1.f;
  bp->sound_speed_subgrid_gas = -1.f;
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
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->f_visc = FLT_MAX;
  bp->accretion_boost_factor = -FLT_MAX;
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */

  /* Reset the rays carried by this BH */
  ray_init(bp->rays, eagle_blackhole_number_of_rays);
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

    /* Count the jump */
    bp->number_of_repositions++;
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

  /* --- Finish the calculation by inserting the missing h factors --- */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  const float rho_inv = 1.f / bp->rho_gas;

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  bp->internal_energy_gas *= h_inv_dim * rho_inv;
  bp->velocity_gas[0] *= h_inv_dim * rho_inv;
  bp->velocity_gas[1] *= h_inv_dim * rho_inv;
  bp->velocity_gas[2] *= h_inv_dim * rho_inv;
  bp->circular_velocity_gas[0] *= h_inv_dim * rho_inv;
  bp->circular_velocity_gas[1] *= h_inv_dim * rho_inv;
  bp->circular_velocity_gas[2] *= h_inv_dim * rho_inv;

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv) */
  bp->circular_velocity_gas[0] *= h_inv;
  bp->circular_velocity_gas[1] *= h_inv;
  bp->circular_velocity_gas[2] *= h_inv;
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

  warning(
      "BH particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      bp->id, bp->h, bp->density.wcount);

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;

  bp->velocity_gas[0] = FLT_MAX;
  bp->velocity_gas[1] = FLT_MAX;
  bp->velocity_gas[2] = FLT_MAX;

  bp->internal_energy_gas = -FLT_MAX;
}

/**
 * @brief Return the current instantaneous accretion rate of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accretion_rate(const struct bpart* bp) {
  return bp->accretion_rate;
}

/**
 * @brief Return the total accreted gas mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accreted_mass(const struct bpart* bp) {
  return bp->total_accreted_mass;
}

/**
 * @brief Return the subgrid mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_subgrid_mass(const struct bpart* bp) {
  return bp->subgrid_mass;
}

/**
 * @brief Return the current bolometric luminosity of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart* bp,
                                      const struct phys_const* constants) {
  return 0.;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart* bp, const struct phys_const* constants) {
  return 0.;
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
                       (bp->x[1] - p->x[1]) * cosmo->a,
                       (bp->x[2] - p->x[2]) * cosmo->a};

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

  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "BH %lld swallowing gas particle %lld "
      "(Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      bp->id, p->id, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_add_part_to_bpart(bp_chem, p_chem, gas_mass);

  /* This BH swallowed a gas particle */
  bp->number_of_gas_swallows++;
  bp->number_of_direct_gas_swallows++;

  /* This BH lost a neighbour */
  bp->num_ngbs--;
  bp->ngb_mass -= gas_mass;

  /* The ray(s) should not point to the no-longer existing particle */
  ray_reset_part_id(bp->rays, eagle_blackhole_number_of_rays, p->id);
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
 * @param constants The physical constants in internal units.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_bpart(
    struct bpart* bpi, const struct bpart* bpj, const struct cosmology* cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props* props, const struct phys_const* constants) {

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

  /* Collect the swallowed angular momentum */
  bpi->swallowed_angular_momentum[0] += bpj->swallowed_angular_momentum[0];
  bpi->swallowed_angular_momentum[1] += bpj->swallowed_angular_momentum[1];
  bpi->swallowed_angular_momentum[2] += bpj->swallowed_angular_momentum[2];

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

  /* Add up all the gas particles we swallowed */
  bpi->number_of_gas_swallows += bpj->number_of_gas_swallows;

  /* Add the subgrid angular momentum that we swallowed */
  bpi->accreted_angular_momentum[0] += bpj->accreted_angular_momentum[0];
  bpi->accreted_angular_momentum[1] += bpj->accreted_angular_momentum[1];
  bpi->accreted_angular_momentum[2] += bpj->accreted_angular_momentum[2];

  /* We had another merger */
  bpi->number_of_mergers++;
}

/**
 * @brief Computes the temperature increase delta_T for black hole feedback.
 *
 * This is calculated as delta_T = min(max(dT_crit, T_gas), dT_num, dT_max):
 * dT_crit = f_crit * critical temperature for suppressing numerical losses
 * T_gas = f_gas * temperature of ambient gas
 * dT_num = temperature increase affordable if N particles should be heated
 * dT_max = maximum allowed energy increase.
 *
 * @param bp The #bpart.
 * @param props The properties of the black hole model.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double black_hole_feedback_delta_T(
    const struct bpart* bp, const struct black_holes_props* props,
    const struct cosmology* cosmo) {

  /* If we do not want a variable delta T, we can stop right here. */
  if (!props->use_variable_delta_T) return props->AGN_delta_T_desired;

  if (bp->internal_energy_gas < 0)
    error("Attempting to compute feedback energy for BH without neighbours.");

  /* Black hole properties */
  const double n_gas_phys = bp->rho_gas * cosmo->a3_inv * props->rho_to_n_cgs;
  const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);
  const double T_gas = bp->internal_energy_gas *
                       cosmo->a_factor_internal_energy /
                       props->temp_to_u_factor;

  /* Calculate base line delta T from BH subgrid mass. The assumption is that
   * the BH mass scales (via halo mass) with the virial temperature, so that
   * this aims for delta T > T_vir. */
  double delta_T = props->AGN_delta_T_mass_norm *
                   pow((bp->subgrid_mass / props->AGN_delta_T_mass_reference),
                       props->AGN_delta_T_mass_exponent);

  /* If desired, also make sure that delta T is not below the numerically
   * critical temperature or that of the ambient gas */
  if (props->AGN_with_locally_adaptive_delta_T) {

    /* Critical temperature for numerical efficiency, based on equation 18
     * of Dalla Vecchia & Schaye (2012) */
    const double T_crit =
        3.162e7 * pow(n_gas_phys * 0.1, 0.6666667) *
        pow(mean_ngb_mass * props->mass_to_solar_mass * 1e-6, 0.33333333);
    delta_T = max(delta_T, T_crit * props->AGN_delta_T_crit_factor);

    /* Delta_T should be (at least) some multiple of the local gas T */
    delta_T = max(delta_T, T_gas * props->AGN_delta_T_background_factor);
  }

  /* Respect the limits */
  delta_T = max(delta_T, props->AGN_delta_T_min);
  return min(delta_T, props->AGN_delta_T_max);
}

/**
 * @brief Computes the energy reservoir threshold for AGN feedback.
 *
 * If adaptive, this is proportional to the accretion rate, with an
 * asymptotic upper limit.
 *
 * @param bp The #bpart.
 * @param props The properties of the black hole model.
 */
__attribute__((always_inline)) INLINE static double
black_hole_energy_reservoir_threshold(struct bpart* bp,
                                      const struct black_holes_props* props) {

  /* If we want a constant threshold, this is short and sweet. */
  if (!props->use_adaptive_energy_reservoir_threshold)
    return props->num_ngbs_to_heat;

  double num_to_heat = props->nheat_alpha *
                       (bp->accretion_rate / props->nheat_maccr_normalisation);

  /* Impose smooth truncation of num_to_heat towards props->nheat_limit */
  if (num_to_heat > props->nheat_alpha) {
    const double coeff_b = 1. / (props->nheat_limit - props->nheat_alpha);
    const double coeff_a = exp(coeff_b * props->nheat_alpha) / coeff_b;
    num_to_heat = props->nheat_limit - coeff_a * exp(-coeff_b * num_to_heat);
  }

  bp->num_ngbs_to_heat = num_to_heat;
  return num_to_heat;
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantities
 * required for the feedback loop.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy fllor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 * @param ti_begin The time at which the step begun (ti_current).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {

  /* Record that the black hole has another active time step */
  bp->number_of_time_steps++;

  if (dt == 0. || bp->rho_gas == 0.) return;

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
  const int with_angmom_limiter = props->with_angmom_limiter;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (all internal units).
   * Note: for the velocities this means peculiar velocities */
  double gas_c_phys = bp->sound_speed_gas * cosmo->a_factor_sound_speed;
  double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double gas_v_circular[3] = {
      bp->circular_velocity_gas[0] * cosmo->a_inv,
      bp->circular_velocity_gas[1] * cosmo->a_inv,
      bp->circular_velocity_gas[2] * cosmo->a_inv};

  /* Norm of the circular velocity of the gas around the BH */
  const double tangential_velocity2 = gas_v_circular[0] * gas_v_circular[0] +
                                      gas_v_circular[1] * gas_v_circular[1] +
                                      gas_v_circular[2] * gas_v_circular[2];
  const double tangential_velocity = sqrt(tangential_velocity2);

  /* We can now compute the Bondi accretion rate (internal units) */
  double Bondi_rate;

  if (props->use_multi_phase_bondi) {

    /* In this case, we are in 'multi-phase-Bondi' mode -- otherwise,
     * the accretion_rate is still zero (was initialised to this) */
    const float hi_inv = 1.f / bp->h;
    const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */

    Bondi_rate = bp->accretion_rate *
                 (4. * M_PI * G * G * BH_mass * BH_mass * hi_inv_dim);
  } else {

    /* Standard approach: compute accretion rate for all gas simultaneously */

    /* Convert velocities to physical frame
     * Note: velocities are already in black hole frame. */
    const double gas_v_phys[3] = {bp->velocity_gas[0] * cosmo->a_inv,
                                  bp->velocity_gas[1] * cosmo->a_inv,
                                  bp->velocity_gas[2] * cosmo->a_inv};

    const double gas_v_norm2 = gas_v_phys[0] * gas_v_phys[0] +
                               gas_v_phys[1] * gas_v_phys[1] +
                               gas_v_phys[2] * gas_v_phys[2];

    if (props->use_subgrid_bondi) {

      /* Use subgrid rho and c for Bondi model */

      /* Construct basic properties of the gas around the BH in
         physical coordinates */
      const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
      const double gas_u_phys =
          bp->internal_energy_gas * cosmo->a_factor_internal_energy;
      const double gas_P_phys =
          gas_pressure_from_internal_energy(gas_rho_phys, gas_u_phys);

      /* Assume primordial abundance and solar metallicity pattern
       * (Yes, that is inconsitent but makes no difference) */
      const double logZZsol = 0.;
      const double XH = 0.75;
      float abundance_ratio[colibre_cooling_N_elementtypes];
      for (int i = 0; i < colibre_cooling_N_elementtypes; ++i)
        abundance_ratio[i] = 1.f;

      /* Get the gas temperature */
      const float gas_T = cooling_get_temperature_from_gas(
          constants, cosmo, cooling, gas_rho_phys, logZZsol, XH, gas_u_phys,
          /*HII_region=*/0);
      const float log10_gas_T = log10f(gas_T);

      /* Internal energy on the entropy floor */
      const double P_EOS = entropy_floor_gas_pressure(gas_rho_phys, bp->rho_gas,
                                                      cosmo, floor_props);
      const double u_EOS =
          gas_internal_energy_from_pressure(gas_rho_phys, P_EOS);
      const double u_EOS_max = u_EOS * exp10(cooling->dlogT_EOS);

      const float log10_u_EOS_max_cgs =
          log10f(u_EOS_max * cooling->internal_energy_to_cgs + FLT_MIN);

      /* Compute the subgrid density assuming pressure
       * equilibirum if on the entropy floor */
      const double rho_sub = compute_subgrid_property(
          cooling, constants, floor_props, cosmo, gas_rho_phys, logZZsol, XH,
          gas_P_phys, log10_gas_T, log10_u_EOS_max_cgs, /*HII_region=*/0,
          abundance_ratio, 0.f, cooling_compute_subgrid_density);

      /* Record what we used */
      bp->rho_subgrid_gas = rho_sub;

      /* And the subgrid sound-speed */
      const float c_sub = gas_soundspeed_from_pressure(rho_sub, gas_P_phys);

      /* Also update the sound-speed to use in the angular momentum limiter */
      gas_c_phys = c_sub;
      gas_c_phys2 = c_sub * c_sub;

      /* Record what we used */
      bp->sound_speed_subgrid_gas = c_sub;

      /* Now, compute the Bondi rate based on the normal velocities and
       * the subgrid density and sound-speed */
      const double denominator2 = gas_v_norm2 + gas_c_phys2;
#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure that the denominator is strictly positive */
      if (denominator2 <= 0)
        error(
            "Invalid denominator for black hole particle %lld in Bondi rate "
            "calculation.",
            bp->id);
#endif
      const double denominator_inv = 1. / sqrt(denominator2);
      Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * rho_sub *
                   denominator_inv * denominator_inv * denominator_inv;

    } else {

      /* Use dynamical rho and c for Bondi model */

      const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;

      const double denominator2 = gas_v_norm2 + gas_c_phys2;
#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure that the denominator is strictly positive */
      if (denominator2 <= 0)
        error(
            "Invalid denominator for black hole particle %lld in Bondi rate "
            "calculation.",
            bp->id);
#endif

      const double denominator_inv = 1. / sqrt(denominator2);
      Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                   denominator_inv * denominator_inv * denominator_inv;
    }
  }

  /* Compute the boost factor from Booth & Schaye (2009) */
  if (props->with_boost_factor) {
    const double XH = 0.75;
    const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
    const double n_H = gas_rho_phys * XH / proton_mass;
    const double boost_ratio = n_H / props->boost_n_h_star;

    double boost_factor = 1.;
    if (props->boost_alpha_only) {
      boost_factor = props->boost_alpha;
    } else {

      /* Booth & Schaye (2009), eq. 4 */
      if (n_H > props->boost_n_h_star) {
        boost_factor = pow(boost_ratio, props->boost_beta);
      }
    }
    Bondi_rate *= boost_factor;
    bp->accretion_boost_factor = boost_factor;
  } else {
    bp->accretion_boost_factor = 1.;
  }

  /* Compute the reduction factor from Rosas-Guevara et al. (2015) */
  if (with_angmom_limiter) {
    const double Bondi_radius = G * BH_mass / gas_c_phys2;
    const double Bondi_time = Bondi_radius / gas_c_phys;
    const double r_times_v_tang = Bondi_radius * tangential_velocity;
    const double r_times_v_tang_3 =
        r_times_v_tang * r_times_v_tang * r_times_v_tang;
    const double viscous_time =
        2. * M_PI * r_times_v_tang_3 /
        (1e-6 * props->alpha_visc * G * G * BH_mass * BH_mass);

    const double f_visc = min(Bondi_time / viscous_time, 1.);
    bp->f_visc = f_visc;

    /* Limit the accretion rate by the Bondi-to-viscous time ratio */
    Bondi_rate *= f_visc;
  } else {
    bp->f_visc = 1.0;
  }

  /* Compute the Eddington rate (internal units) */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (epsilon_r * c * sigma_Thomson);

  /* Should we record this time as the most recent high accretion rate? */
  if (Bondi_rate > f_Edd_recording * Eddington_rate) {
    if (with_cosmology) {
      bp->last_high_Eddington_fraction_scale_factor = cosmo->a;
    } else {
      bp->last_high_Eddington_fraction_time = time;
    }
  }

  /* Limit the accretion rate to a fraction of the Eddington rate */
  const double accr_rate = min(Bondi_rate, f_Edd * Eddington_rate);
  bp->accretion_rate = accr_rate;
  bp->eddington_fraction = Bondi_rate / Eddington_rate;

  /* Factor in the radiative efficiency */
  const double mass_rate = (1. - epsilon_r) * accr_rate;
  const double luminosity = epsilon_r * accr_rate * c * c;

  /* Integrate forward in time */
  bp->subgrid_mass += mass_rate * dt;
  bp->total_accreted_mass += mass_rate * dt;
  bp->energy_reservoir += luminosity * epsilon_f * dt;

  if (props->use_nibbling && bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative losses, so we must decrease the particle mass
     * in proprtion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -= epsilon_r * accr_rate * dt;
    if (bp->mass < 0)
      error("Black hole %lld reached negative mass (%g). Trouble ahead...",
            bp->id, bp->mass);
  }

  /* Increase the subgrid angular momentum according to what we accreted
   * Note that this is already in physical units, a factors from velocity and
   * radius cancel each other. Also, the circular velocity contains an extra
   * smoothing length factor that we undo here. */
  bp->accreted_angular_momentum[0] +=
      bp->circular_velocity_gas[0] * mass_rate * dt / bp->h;
  bp->accreted_angular_momentum[1] +=
      bp->circular_velocity_gas[1] * mass_rate * dt / bp->h;
  bp->accreted_angular_momentum[2] +=
      bp->circular_velocity_gas[2] * mass_rate * dt / bp->h;

  /* Below we compute energy required to have a feedback event(s)
   * Note that we have subtracted the particles we swallowed from the ngb_mass
   * and num_ngbs accumulators. */

  /* Now find the temperature increase for a possible feedback event */
  const double delta_T = black_hole_feedback_delta_T(bp, props, cosmo);
  bp->AGN_delta_T = delta_T;
  double delta_u = delta_T * props->temp_to_u_factor;
  const double delta_u_ref =
      props->AGN_use_nheat_with_fixed_dT
          ? props->AGN_delta_T_desired * props->temp_to_u_factor
          : delta_u;

  /* Energy required to have a feedback event
   * Note that we have subtracted the particles we swallowed from the ngb_mass
   * and num_ngbs accumulators. */
  const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);
  const double E_feedback_event =
      num_ngbs_to_heat * delta_u_ref * mean_ngb_mass;

  /* Compute and store BH accretion-limited time-step */
  if (luminosity > 0.) {
    const float dt_acc = delta_u * mean_ngb_mass * props->num_ngbs_to_heat /
                         (luminosity * props->epsilon_f);
    bp->dt_heat = max(dt_acc, props->time_step_min);
  } else {
    bp->dt_heat = FLT_MAX;
  }

  /* Are we doing some feedback? */
  if (bp->energy_reservoir > E_feedback_event) {

    int number_of_energy_injections;

    /* How are we doing feedback? */
    if (props->AGN_deterministic) {

      number_of_energy_injections =
          (int)(bp->energy_reservoir / (delta_u * mean_ngb_mass));

    } else {

      /* Probability of heating. */
      const double prob = bp->energy_reservoir / (delta_u * bp->ngb_mass);

      /* Compute the number of energy injections based on probability */
      if (prob < 1.) {

        /* Initialise counter of energy injections */
        number_of_energy_injections = 0;

        /* How many AGN energy injections will we get?
         *
         * Note that we use the particles here to draw random numbers. This does
         * not mean that the 'lucky' particles here are the ones that will be
         * heated. If we get N lucky particles, we will use the first N random
         * ray directions in the isotropic case or the first N closest particles
         * in the other modes. */
        for (int i = 0; i < bp->num_ngbs; i++) {
          const double rand = random_unit_interval_part_ID_and_index(
              bp->id, i, ti_begin, random_number_BH_feedback);

          /* Increase the counter if we are lucky */
          if (rand < prob) number_of_energy_injections++;
        }

      } else {

        /* We want to use up all energy avaliable in the reservoir. Therefore,
         * number_of_energy_injections is > or = props->num_ngbs_to_heat */
        number_of_energy_injections =
            (int)(bp->energy_reservoir / (delta_u * mean_ngb_mass));
      }
    }

    /* Maximum number of energy injections allowed */
    const int N_energy_injections_allowed =
        min(eagle_blackhole_number_of_rays, bp->num_ngbs);

    /* If there are more energy-injection events than min(the number of Ngbs in
     * the kernel, maximum number of rays) then lower the number of events &
     * proportionally increase the energy per event */
    if (number_of_energy_injections > N_energy_injections_allowed) {

      /* Increase the thermal energy per event */
      const double alpha_thermal = (double)number_of_energy_injections /
                                   (double)N_energy_injections_allowed;

      delta_u *= alpha_thermal;

      /* Lower the maximum number of events to the max allowed value */
      number_of_energy_injections = N_energy_injections_allowed;
    }

    /* Compute how much energy will be deposited onto the gas */
    /* Note that it will in general be different from E_feedback_event if
     * gas particles are of different mass. */
    double Energy_deposited = 0.0;

    /* Count the number of unsuccessful energy injections (e.g., if the particle
     * that the BH wants to heat has been swallowed and thus no longer exists)
     */
    int N_unsuccessful_energy_injections = 0;

    for (int i = 0; i < number_of_energy_injections; i++) {

      /* If the gas particle that the BH wants to heat has just been swallowed
       * by the same BH, increment the counter of unsuccessful injections. If
       * the particle has not been swallowed by the BH, increase the energy that
       * will later be subtracted from the BH's energy reservoir. */
      if (bp->rays[i].id_min_length != -1)
        Energy_deposited += delta_u * bp->rays[i].mass;
      else
        N_unsuccessful_energy_injections++;
    }

    /* Store all of this in the black hole for delivery onto the gas. */
    bp->to_distribute.AGN_delta_u = delta_u;
    bp->to_distribute.AGN_number_of_energy_injections =
        number_of_energy_injections;

    /* Subtract the deposited energy from the BH energy reservoir. Note
     * that in the stochastic case, the resulting value might be negative.
     * This happens when (due to the probabilistic nature of the model) the
     * BH injects more energy than it actually has in the reservoir. */
    bp->energy_reservoir -= Energy_deposited;

    /* Total number successful energy injections at this time-step. In each
     * energy injection, a certain gas particle from the BH's kernel gets
     * heated. (successful = the particle(s) that is going to get heated by
     * this BH has not been swallowed by the same BH). */
    const int N_successful_energy_injections =
        number_of_energy_injections - N_unsuccessful_energy_injections;

    /* Increase the number of energy injections the black hole has heated so
     * far. Note that in the isotropic model, a gas particle may receive AGN
     * energy several times at the same time-step. In this case, the number of
     * particles heated at this time-step for this BH will be smaller than the
     * total number of energy injections for this BH. */
    bp->AGN_number_of_energy_injections += N_successful_energy_injections;

    /* Increase the number of AGN events the black hole has had so far.
     * If the BH does feedback, the number of AGN events is incremented by one.
     */
    bp->AGN_number_of_AGN_events += N_successful_energy_injections > 0;

    /* Update the total (cumulative) energy used for gas heating in AGN feedback
     * by this BH */
    bp->AGN_cumulative_energy += Energy_deposited;

    /* Store the time/scale factor when the BH last did AGN feedback */
    if (N_successful_energy_injections) {
      if (with_cosmology) {
        bp->last_AGN_event_scale_factor = cosmo->a;
      } else {
        bp->last_AGN_event_time = time;
      }
    }

  } else {

    /* Flag that we don't want to heat anyone */
    bp->to_distribute.AGN_number_of_energy_injections = 0;
    bp->to_distribute.AGN_delta_u = 0.f;
  }
}

/**
 * @brief Computes the (maximal) repositioning speed for a black hole.
 *
 * Calculated as upsilon * (m_BH / m_ref) ^ beta_m * (n_H_BH / n_ref) ^ beta_n
 * where m_BH = BH subgrid mass, n_H_BH = physical gas density around BH
 * and upsilon, m_ref, beta_m, n_ref, and beta_n are parameters.
 *
 * @param bp The #bpart.
 * @param props The properties of the black hole model.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_repositioning_speed(const struct bpart* restrict bp,
                                    const struct black_holes_props* props,
                                    const struct cosmology* cosmo) {

  const double n_gas_phys = bp->rho_gas * cosmo->a3_inv * props->rho_to_n_cgs;
  const double v_repos =
      props->reposition_coefficient_upsilon *
      pow(bp->subgrid_mass / props->reposition_reference_mass,
          props->reposition_exponent_mass) *
      pow(n_gas_phys / props->reposition_reference_n_H,
          props->reposition_exponent_n_H);

  /* Make sure the repositioning is not back-firing... */
  if (v_repos < 0)
    error(
        "BH %lld wants to reposition at negative speed (%g U_V). Do you "
        "think you are being funny? No-one is laughing.",
        bp->id, v_repos);

  return v_repos;
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
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the temp
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const double dt, const integertime_t ti_begin) {

  /* First check: did we find any eligible neighbour particle to jump to? */
  if (bp->reposition.min_potential != FLT_MAX) {

    /* Record that we have a (possible) repositioning situation */
    bp->number_of_reposition_attempts++;

    /* Is the potential lower (i.e. the BH is at the bottom already)
     * OR is the BH massive enough that we don't reposition? */
    const float potential = gravity_get_comoving_potential(bp->gpart);
    if (potential < bp->reposition.min_potential ||
        bp->subgrid_mass > props->max_reposition_mass) {

      /* No need to reposition */
      bp->reposition.min_potential = FLT_MAX;
      bp->reposition.delta_x[0] = -FLT_MAX;
      bp->reposition.delta_x[1] = -FLT_MAX;
      bp->reposition.delta_x[2] = -FLT_MAX;

    } else if (props->set_reposition_speed) {

      /* If we are re-positioning, move the BH a fraction of delta_x, so
       * that we have a well-defined re-positioning velocity (repos_vel
       * cannot be negative). */
      double repos_vel = black_holes_get_repositioning_speed(bp, props, cosmo);

      /* Convert target reposition velocity to a fractional reposition
       * along reposition.delta_x */
      const double dx = bp->reposition.delta_x[0];
      const double dy = bp->reposition.delta_x[1];
      const double dz = bp->reposition.delta_x[2];
      const double d = sqrt(dx * dx + dy * dy + dz * dz);

      /* Exclude the pathological case of repositioning by zero distance */
      if (d > 0) {
        double repos_frac = repos_vel * dt / d;

        /* We should never get negative repositioning fractions... */
        if (repos_frac < 0)
          error("Wanting to reposition by negative fraction (%g)?", repos_frac);

        /* ... but fractions > 1 can occur if the target velocity is high.
         * We do not want this, because it could lead to overshooting the
         * actual potential minimum. */
        if (repos_frac > 1) {
          repos_frac = 1.;
          repos_vel = repos_frac * d / dt;
        }

        bp->last_repos_vel = (float)repos_vel;
        bp->reposition.delta_x[0] *= repos_frac;
        bp->reposition.delta_x[1] *= repos_frac;
        bp->reposition.delta_x[2] *= repos_frac;
      }

      /* ends section for fractional repositioning */
    } else {

      /* We _should_ reposition, but not fractionally. Here, we will
       * reposition exactly on top of another gas particle - which
       * could cause issues, so we add on a small fractional offset
       * of magnitude 0.001 h in the reposition delta. */

      /* Generate three random numbers in the interval [-0.5, 0.5[; id,
       * id**2, and id**3 are required to give unique random numbers (as
       * random_unit_interval is completely reproducible). */
      const float offset_dx =
          random_unit_interval(bp->id, ti_begin, random_number_BH_reposition) -
          0.5f;
      const float offset_dy =
          random_unit_interval(bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;
      const float offset_dz =
          random_unit_interval(bp->id * bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;

      const float length_inv =
          1.0f / sqrtf(offset_dx * offset_dx + offset_dy * offset_dy +
                       offset_dz * offset_dz);

      const float norm = 0.001f * bp->h * length_inv;

      bp->reposition.delta_x[0] += offset_dx * norm;
      bp->reposition.delta_x[1] += offset_dy * norm;
      bp->reposition.delta_x[2] += offset_dz * norm;
    }
  } /* ends section if we found eligible repositioning target(s) */
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

  bp->to_distribute.AGN_delta_u = 0.f;
  bp->to_distribute.AGN_number_of_energy_injections = 0;

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

/**
 * @brief Store the gravitational potential of a black hole by copying it from
 * its #gpart friend.
 *
 * @param bp The black hole particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_bpart(struct bpart* bp, const struct gpart* gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->gpart != gp) error("Copying potential to the wrong black hole!");
#endif

  bp->reposition.potential = gp->potential;
}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The black hole data of a gas particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_part(struct black_holes_part_data* p_data,
                                    const struct gpart* gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 * @param xp The #xpart that became a black hole.
 * @param ti_current the current time on the time-line.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart* bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct part* p, const struct xpart* xp,
    const integertime_t ti_current) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time and density */
  bp->formation_scale_factor = cosmo->a;
  bp->formation_gas_density = hydro_get_physical_density(p, cosmo);

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

  /* We haven't accreted anything yet */
  bp->total_accreted_mass = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_time_steps = 0;

  /* Initialise the energy reservoir threshold to the constant default */
  bp->num_ngbs_to_heat = props->num_ngbs_to_heat; /* Filler value */

  /* We haven't repositioned yet, nor attempted it */
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->last_repos_vel = 0.f;

  /* Copy over the splitting struct */
  bp->split_data = xp->split_data;

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

  /* Set the initial targetted heating temperature, used for the
   * BH time step determination */
  bp->AGN_delta_T = props->AGN_delta_T_desired;

  /* First initialisation */
  black_holes_init_bpart(bp);

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

#endif /* SWIFT_EAGLE_BLACK_HOLES_H */
