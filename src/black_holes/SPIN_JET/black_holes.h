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

#ifndef SWIFT_SPIN_JET_BLACK_HOLES_H
#define SWIFT_SPIN_JET_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_spin.h"
#include "black_holes_struct.h"
#include "cooling_properties.h"
#include "cosmology.h"
#include "dimension.h"
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

  /* Do something if and only if the accretion rate is non-zero */
  if (bp->accretion_rate > 0.f) {

    /* Gather some physical constants (all in internal units) */
    const double c = constants->const_speed_light_c;

    /* Compute instantaneous energy supply rate to the BH energy reservoir
     * which is proportional to the BH mass accretion rate */
    const double Energy_rate =
        (bp->radiative_efficiency * props->epsilon_f + bp->wind_efficiency) *
        bp->accretion_rate * c * c;

    /* Compute instantaneous jet energy supply rate to the BH jet reservoir
     * which is proportional to the BH mass accretion rate */
    const double Jet_rate =
        bp->jet_efficiency * props->eps_f_jet * bp->accretion_rate * c * c;

    /* Compute average heating energy in AGN feedback, both thermal and jet */

    /* Average particle mass in BH's kernel */
    const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);
    /* Without multiplying by mean_ngb_mass we'd get energy per unit mass */
    const double E_heat = bp->delta_T * props->temp_to_u_factor * mean_ngb_mass;
    const double E_jet = 0.5 * mean_ngb_mass * bp->v_jet * bp->v_jet;

    /* Compute average time between energy injections for the given accretion
     * rate. The time is multiplied by the number of Ngbs to heat because
     * if more particles are heated at once then the time between different
     * AGN feedback events increases proportionally. */
    const double dt_heat = E_heat * props->num_ngbs_to_heat / Energy_rate;

    /* Similar for jets, with N_jet the target number of particles to kick. */
    const double dt_jet = E_jet * props->N_jet / Jet_rate;

    /* Pick the shorter of the two feedback time steps */
    double bh_timestep = min(dt_jet, dt_heat);

    /* See if the spin evolution time step is even smaller */
    bh_timestep = min(bh_timestep, bp->dt_ang_mom);

    /* The final timestep of the BH cannot be smaller than the miminum allowed
     * time-step */
    bh_timestep = max(bh_timestep, props->time_step_min);

    return bh_timestep;
  }
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
  if (props->use_subgrid_mass_from_ics == 0) {
    bp->subgrid_mass = bp->mass;

  } else if (props->with_subgrid_mass_check && bp->subgrid_mass <= 0)
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
  bp->last_AGN_event_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;
  bp->accreted_angular_momentum[0] = 0.f;
  bp->accreted_angular_momentum[1] = 0.f;
  bp->accreted_angular_momentum[2] = 0.f;
  bp->dt_heat = 0.f;
  bp->AGN_number_of_AGN_events = 0;
  bp->AGN_number_of_energy_injections = 0;
  bp->AGN_cumulative_energy = 0.f;
  bp->jet_efficiency = 0.1f;
  bp->radiative_efficiency = 0.1f;
  bp->accretion_disk_angle = 0.01f;
  bp->accretion_mode = BH_thin_disc;
  bp->eddington_fraction = 0.01f;
  bp->jet_reservoir = 0.f;
  bp->total_jet_energy = 0.f;
  bp->wind_efficiency = 0.f;
  bp->wind_energy = 0.f;
  bp->radiated_energy = 0.f;
  bp->accretion_efficiency = 1.f;
  bp->dt_jet = 0.f;
  bp->dt_ang_mom = 0.f;
  bp->AGN_number_of_AGN_jet_events = 0;
  bp->AGN_number_of_jet_injections = 0;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->accreted_mass_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->thermal_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->jet_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->wind_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->radiated_energy_by_mode[i] = 0.f;
  bp->jet_direction[0] = bp->angular_momentum_direction[0];
  bp->jet_direction[1] = bp->angular_momentum_direction[1];
  bp->jet_direction[2] = bp->angular_momentum_direction[2];

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
  bp->rho_gas_hot = 0.f;
  bp->sound_speed_gas = 0.f;
  bp->sound_speed_gas_hot = 0.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->spec_angular_momentum_gas[0] = 0.f;
  bp->spec_angular_momentum_gas[1] = 0.f;
  bp->spec_angular_momentum_gas[2] = 0.f;
  bp->curl_v_gas[0] = 0.f;
  bp->curl_v_gas[1] = 0.f;
  bp->curl_v_gas[2] = 0.f;
  bp->velocity_dispersion_gas = 0.f;
  bp->ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->f_visc = FLT_MAX;
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */
  bp->rho_gas_hot = 0.f;
  bp->sound_speed_gas_hot = 0.f;

  /* Reset the rays carried by this BH */
  ray_init(bp->rays, spinjet_blackhole_number_of_rays);
  ray_init(bp->rays_jet, spinjet_blackhole_number_of_rays);
  ray_init(bp->rays_jet_pos, spinjet_blackhole_number_of_rays);
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
  bp->rho_gas_hot *= h_inv_dim;
  const float rho_inv = 1.f / bp->rho_gas;

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  if (bp->rho_gas_hot > 0.)
    bp->sound_speed_gas_hot *= h_inv_dim / bp->rho_gas_hot;
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  if (bp->rho_gas_hot > 0.) {
    bp->sound_speed_gas_hot *= h_inv_dim / bp->rho_gas_hot;
  }
  bp->velocity_gas[0] *= h_inv_dim * rho_inv;
  bp->velocity_gas[1] *= h_inv_dim * rho_inv;
  bp->velocity_gas[2] *= h_inv_dim * rho_inv;
  bp->velocity_dispersion_gas *= h_inv_dim * rho_inv;
  bp->spec_angular_momentum_gas[0] *= h_inv_dim * rho_inv;
  bp->spec_angular_momentum_gas[1] *= h_inv_dim * rho_inv;
  bp->spec_angular_momentum_gas[2] *= h_inv_dim * rho_inv;

  /* ... and for the curl, we also need to divide by an extra h factor */
  bp->curl_v_gas[0] *= h_inv_dim_plus_one * rho_inv;
  bp->curl_v_gas[1] *= h_inv_dim_plus_one * rho_inv;
  bp->curl_v_gas[2] *= h_inv_dim_plus_one * rho_inv;

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv) */
  bp->circular_velocity_gas[0] = bp->spec_angular_momentum_gas[0] * h_inv;
  bp->circular_velocity_gas[1] = bp->spec_angular_momentum_gas[1] * h_inv;
  bp->circular_velocity_gas[2] = bp->spec_angular_momentum_gas[2] * h_inv;

  /* Calculate (actual) gas velocity dispersion. Currently, the variable
   * 'velocity_dispersion_gas' holds <v^2> instead. */
  const double speed_gas2 = bp->velocity_gas[0] * bp->velocity_gas[0] +
                            bp->velocity_gas[1] * bp->velocity_gas[1] +
                            bp->velocity_gas[2] * bp->velocity_gas[2];

  bp->velocity_dispersion_gas -= speed_gas2;
  bp->velocity_dispersion_gas = sqrt(fabs(bp->velocity_dispersion_gas));
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

  bp->velocity_gas[0] = FLT_MAX;
  bp->velocity_gas[1] = FLT_MAX;
  bp->velocity_gas[2] = FLT_MAX;

  bp->velocity_dispersion_gas = FLT_MAX;

  bp->curl_v_gas[0] = FLT_MAX;
  bp->curl_v_gas[1] = FLT_MAX;
  bp->curl_v_gas[2] = FLT_MAX;
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
  const double c = constants->const_speed_light_c;
  return bp->accretion_rate * bp->radiative_efficiency * c * c;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart* bp, const struct phys_const* constants) {
  const double c = constants->const_speed_light_c;
  return bp->accretion_rate * bp->jet_efficiency * c * c;
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
  ray_reset_part_id(bp->rays, spinjet_blackhole_number_of_rays, p->id);
  ray_reset_part_id(bp->rays_jet, spinjet_blackhole_number_of_rays, p->id);
  ray_reset_part_id(bp->rays_jet_pos, spinjet_blackhole_number_of_rays, p->id);
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

  /* Evolve the black hole spin. */
  black_hole_merger_spin_evolve(bpi, bpj, constants);

  /* Increase the masses of the BH. */
  bpi->mass += bpj->mass;
  bpi->gpart->mass += bpj->mass;
  bpi->subgrid_mass += bpj->subgrid_mass;

  /* We need to see if the new spin is positive or negative, by implementing
     King et al. (2005) condition of (counter-)alignment. */
  const float gas_spec_ang_mom_norm = sqrtf(
      bpi->spec_angular_momentum_gas[0] * bpi->spec_angular_momentum_gas[0] +
      bpi->spec_angular_momentum_gas[1] * bpi->spec_angular_momentum_gas[1] +
      bpi->spec_angular_momentum_gas[2] * bpi->spec_angular_momentum_gas[2]);

  float dot_product = -1.;
  if (gas_spec_ang_mom_norm > 0.) {
    dot_product = 1. / gas_spec_ang_mom_norm *
                  (bpi->spec_angular_momentum_gas[0] *
                       bpi->angular_momentum_direction[0] +
                   bpi->spec_angular_momentum_gas[1] *
                       bpi->angular_momentum_direction[1] +
                   bpi->spec_angular_momentum_gas[2] *
                       bpi->angular_momentum_direction[2]);
  } else {
    dot_product = 0.;
  }

  if (black_hole_angular_momentum_magnitude(bpi, constants) * dot_product <
      -0.5 * black_hole_warp_angular_momentum(bpi, constants, props)) {
    bpi->spin = -1. * fabsf(bpi->spin);
  } else {
    bpi->spin = fabsf(bpi->spin);
  }

  /* Update various quantities with new spin */
  black_hole_select_accretion_mode(bpi, props);
  bpi->accretion_disk_angle = dot_product;
  bpi->radiative_efficiency = black_hole_radiative_efficiency(bpi, props);
  bpi->jet_efficiency = black_hole_jet_efficiency(bpi, props);
  bpi->wind_efficiency = black_hole_wind_efficiency(bpi, props);

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

  /* Update the jet reservoir */
  bpi->jet_reservoir += bpj->jet_reservoir;

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
 * @brief Compute the accretion rate of the black hole and all the quantites
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
 * @param ti_begin Integer time value at the beginning of timestep
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

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->num_ngbs <= 0) {
    error(
        "The number of BH neighbours is %d, despite the fact that the gas "
        " density in the BH kernel is non-zero.",
        bp->num_ngbs);
  }
#endif

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (all internal units).
   * Note: for the velocities this means peculiar velocities */
  const double gas_c_phys = bp->sound_speed_gas * cosmo->a_factor_sound_speed;
  const double gas_c_phys2 = gas_c_phys * gas_c_phys;
  const double gas_v_circular[3] = {
      bp->circular_velocity_gas[0] * cosmo->a_inv,
      bp->circular_velocity_gas[1] * cosmo->a_inv,
      bp->circular_velocity_gas[2] * cosmo->a_inv};

  /* Norm of the circular velocity of the gas around the BH */
  const double tangential_velocity2 = gas_v_circular[0] * gas_v_circular[0] +
                                      gas_v_circular[1] * gas_v_circular[1] +
                                      gas_v_circular[2] * gas_v_circular[2];
  const double tangential_velocity = sqrt(tangential_velocity2);

  /* We can now compute the accretion rate (internal units) */
  double accr_rate;

  if (props->use_multi_phase_bondi) {

    /* In this case, we are in 'multi-phase-Bondi' mode -- otherwise,
     * the accretion_rate is still zero (was initialised to this) */
    const float hi_inv = 1.f / bp->h;
    const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
    accr_rate = bp->accretion_rate *
                (4. * M_PI * G * G * BH_mass * BH_mass * hi_inv_dim);
  } else {

    /* Standard approach: compute accretion rate for all gas simultaneously.
     *
     * Convert the quantities we gathered to physical frame (all internal
     * units). Note: velocities are already in black hole frame. */
    const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
    const double gas_v_phys[3] = {bp->velocity_gas[0] * cosmo->a_inv,
                                  bp->velocity_gas[1] * cosmo->a_inv,
                                  bp->velocity_gas[2] * cosmo->a_inv};

    const double gas_v_norm2 = gas_v_phys[0] * gas_v_phys[0] +
                               gas_v_phys[1] * gas_v_phys[1] +
                               gas_v_phys[2] * gas_v_phys[2];

    /* In the Bondi-Hoyle-Lyttleton formula, the bulk flow of gas is
     * added to the sound speed in quadrature. Treated separately (below)
     * in the Krumholz et al. (2006) prescription */
    const double denominator2 =
        props->use_krumholz ? gas_c_phys2 : gas_v_norm2 + gas_c_phys2;
#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure that the denominator is strictly positive */
    if (denominator2 <= 0)
      error(
          "Invalid denominator for BH particle %lld in Bondi rate "
          "calculation.",
          bp->id);
#endif
    const double denominator_inv = 1. / sqrt(denominator2);
    accr_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                denominator_inv * denominator_inv * denominator_inv;

    if (props->use_krumholz) {

      /* Compute the additional correction factors from Krumholz+06,
       * accounting for bulk flow and turbulence of ambient gas. */
      const double lambda = 1.1;
      const double gas_v_dispersion =
          bp->velocity_dispersion_gas * cosmo->a_inv;
      const double mach_turb = gas_v_dispersion / gas_c_phys;
      const double mach_bulk = sqrt(gas_v_norm2) / gas_c_phys;
      const double mach2 = mach_turb * mach_turb + mach_bulk * mach_bulk;
      const double m1 = 1. + mach2;
      const double mach_factor =
          sqrt((lambda * lambda + mach2) / (m1 * m1 * m1 * m1));
      accr_rate *= mach_factor;
    }

    if (props->with_krumholz_vorticity) {

      /* Change the accretion rate to equation (3) of Krumholz et al. (2006)
       * by adding a vorticity-dependent term in inverse quadrature */

      /* Convert curl to vorticity in physical units */
      const double gas_curlv_phys[3] = {bp->curl_v_gas[0] * cosmo->a2_inv,
                                        bp->curl_v_gas[1] * cosmo->a2_inv,
                                        bp->curl_v_gas[2] * cosmo->a2_inv};
      const double gas_vorticity = sqrt(gas_curlv_phys[0] * gas_curlv_phys[0] +
                                        gas_curlv_phys[1] * gas_curlv_phys[1] +
                                        gas_curlv_phys[2] * gas_curlv_phys[2]);

      const double Bondi_radius = G * BH_mass / gas_c_phys2;
      const double omega_star = gas_vorticity * Bondi_radius / gas_c_phys;
      const double f_omega_star = 1.0 / (1.0 + pow(omega_star, 0.9));
      const double mdot_omega = 4. * M_PI * gas_rho_phys * G * G * BH_mass *
                                BH_mass * denominator_inv * denominator_inv *
                                denominator_inv * 0.34 * f_omega_star;

      const double accr_rate_inv = 1. / accr_rate;
      const double mdot_omega_inv = 1. / mdot_omega;
      accr_rate = 1. / sqrt(accr_rate_inv * accr_rate_inv +
                            mdot_omega_inv * mdot_omega_inv);
    } /* ends calculation of vorticity addition to Krumholz prescription */

  } /* ends section without multi-phase accretion */

  /* Compute the boost factor from Booth, Schaye (2009) */
  if (props->with_boost_factor) {
    const double gas_rho_phys = bp->rho_gas * cosmo->a3_inv;
    const double n_H = gas_rho_phys * 0.75 / proton_mass;
    const double boost_ratio = n_H / props->boost_n_h_star;
    const double boost_factor =
        max(pow(boost_ratio, props->boost_beta), props->boost_alpha);
    accr_rate *= boost_factor;
  }
  /* Compute the reduction factor from Rosas-Guevara et al. (2015) */
  if (props->with_angmom_limiter) {
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
    accr_rate *= f_visc;
  } else {
    bp->f_visc = 1.0;
  }

  /* Compute the Eddington rate (internal units). The radiative efficiency
     is taken to be equal to the Novikov-Thorne (1973) one, and thus varies
     with spin between 4% and 40%. */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass /
      (eps_Novikov_Thorne(bp->spin) * c * sigma_Thomson);

  /* Should we record this time as the most recent high accretion rate? */
  if (accr_rate > props->f_Edd_recording * Eddington_rate) {
    if (with_cosmology) {
      bp->last_high_Eddington_fraction_scale_factor = cosmo->a;
    } else {
      bp->last_high_Eddington_fraction_time = time;
    }
  }

  /* Limit the accretion rate to a fraction of the Eddington rate */
  bp->eddington_fraction = accr_rate / Eddington_rate;
  accr_rate = min(accr_rate, props->f_Edd * Eddington_rate);

  /* Include the effects of the accretion efficiency */
  bp->accretion_rate = accr_rate * bp->accretion_efficiency;
  bp->eddington_fraction *= bp->accretion_efficiency;

  /* Change mode based on new accretion rates if necessary, and update
     all important quantities. */
  black_hole_select_accretion_mode(bp, props);
  bp->accretion_efficiency =
      black_hole_accretion_efficiency(bp, props, constants, cosmo);
  bp->accretion_rate = accr_rate * bp->accretion_efficiency;
  bp->eddington_fraction = bp->accretion_rate / Eddington_rate;

  bp->radiative_efficiency = black_hole_radiative_efficiency(bp, props);
  bp->jet_efficiency = black_hole_jet_efficiency(bp, props);
  bp->wind_efficiency = black_hole_wind_efficiency(bp, props);

  /* Define feedback-related quantities that we will update and need later on */
  double luminosity = 0.;
  double jet_power = 0.;

  /* How much mass will be consumed over this time step? */
  double delta_m_0 = bp->accretion_rate * dt;

  /* Norm of the specific angular momentum, will be needed in a moment. */
  float spec_ang_mom_norm = sqrtf(max(
      0.,
      bp->spec_angular_momentum_gas[0] * bp->spec_angular_momentum_gas[0] +
          bp->spec_angular_momentum_gas[1] * bp->spec_angular_momentum_gas[1] +
          bp->spec_angular_momentum_gas[2] * bp->spec_angular_momentum_gas[2]));

  /* Cosine of the angle between the spin vector and the specific angular
     momentum vector of gas around the BH. */
  float dot_product = -1.;
  if (spec_ang_mom_norm > 0.) {
    dot_product =
        1. / spec_ang_mom_norm *
        (bp->spec_angular_momentum_gas[0] * bp->angular_momentum_direction[0] +
         bp->spec_angular_momentum_gas[1] * bp->angular_momentum_direction[1] +
         bp->spec_angular_momentum_gas[2] * bp->angular_momentum_direction[2]);
  } else {
    dot_product = 0.;
  }

  /* Decide if accretion is prograde (spin positive) or retrograde
     (spin negative) based on condition from King et al. (2005) */
  if ((black_hole_angular_momentum_magnitude(bp, constants) * dot_product <
       -0.5 * black_hole_warp_angular_momentum(bp, constants, props)) &&
      (fabsf(bp->spin) > 0.01)) {
    bp->spin = -1. * fabsf(bp->spin);
  } else {
    bp->spin = fabsf(bp->spin);
  }

  /* Calculate how many warp increments the BH will swallow over this time
     step */
  double n_i = 0.;
  if (bp->accretion_rate > 0.) {
    n_i = delta_m_0 / black_hole_warp_mass(bp, constants, props);
  }

  /* Update the angular momentum vector of the BH based on how many
     increments of warp angular momenta have been consumed. */
  if (spec_ang_mom_norm > 0.) {

    /* If spin is at its floor value of 0.01, we immediately redirect
       the spin in the direction of the accreting gas */
    if (fabsf(bp->spin) <= 0.01) {
      bp->angular_momentum_direction[0] =
          bp->spec_angular_momentum_gas[0] / spec_ang_mom_norm;
      bp->angular_momentum_direction[1] =
          bp->spec_angular_momentum_gas[1] / spec_ang_mom_norm;
      bp->angular_momentum_direction[2] =
          bp->spec_angular_momentum_gas[2] / spec_ang_mom_norm;
    } else {
      const double j_warp =
          black_hole_warp_angular_momentum(bp, constants, props);
      const double j_BH = black_hole_angular_momentum_magnitude(bp, constants);
      const double ang_mom_total[3] = {
          bp->angular_momentum_direction[0] * j_BH +
              n_i * bp->spec_angular_momentum_gas[0] / spec_ang_mom_norm *
                  j_warp,
          bp->angular_momentum_direction[1] * j_BH +
              n_i * bp->spec_angular_momentum_gas[1] / spec_ang_mom_norm *
                  j_warp,
          bp->angular_momentum_direction[2] * j_BH +
              n_i * bp->spec_angular_momentum_gas[2] / spec_ang_mom_norm *
                  j_warp};

      /* Modulus of the new J_BH */
      const double modulus = sqrt(ang_mom_total[0] * ang_mom_total[0] +
                                  ang_mom_total[1] * ang_mom_total[1] +
                                  ang_mom_total[2] * ang_mom_total[2]);

      if (modulus > 0.) {
        bp->angular_momentum_direction[0] = ang_mom_total[0] / modulus;
        bp->angular_momentum_direction[1] = ang_mom_total[1] / modulus;
        bp->angular_momentum_direction[2] = ang_mom_total[2] / modulus;
      }
    }
  }

  /* Check if we are fixing the jet along the z-axis. */
  if (props->fix_jet_direction) {
    bp->angular_momentum_direction[0] = 0.;
    bp->angular_momentum_direction[1] = 0.;
    bp->angular_momentum_direction[2] = 1;
  }

  bp->jet_direction[0] = bp->angular_momentum_direction[0];
  bp->jet_direction[1] = bp->angular_momentum_direction[1];
  bp->jet_direction[2] = bp->angular_momentum_direction[2];

  float spin_final = -1.;
  /* Calculate the change in the BH spin */
  if (bp->subgrid_mass > 0.) {
    spin_final = bp->spin + delta_m_0 / bp->subgrid_mass *
                                black_hole_spinup_rate(bp, constants, props);
  } else {
    error(
        "Black hole with id %lld tried to evolve spin with zero "
        "(or less) subgrid mass. ",
        bp->id);
  }

  /* Make sure that the spin does not shoot above 1 or below -1. Physically
     this wouldn't happen because the spinup function goes to 0 at a=1, but
     numerically it may happen due to finite increments. If this happens,
     many spin-related quantities begin to diverge. We also want to avoid
     spin equal to zero, or very close to it. The black hole time steps
     become shorter and shorter as the BH approaches spin 0, so we simply
     'jump' through 0 instead, choosing a small value of 0.01 (in magnitude)
     as a floor for the spin. The spin will always jump to +0.01 since the
     spinup function will always make spin go from negative to positive at
     these small values. */
  if (spin_final > 0.998) {
    spin_final = 0.998;
  } else if (spin_final < -0.998) {
    spin_final = -0.998;
  } else if (fabsf(spin_final) < 0.01) {
    spin_final = 0.01;
  }

  /* Update the spin */
  bp->spin = spin_final;

  /* Update other quantities */
  bp->accretion_disk_angle = dot_product;

  /* Update efficiencies given the new spin */
  bp->radiative_efficiency = black_hole_radiative_efficiency(bp, props);
  bp->jet_efficiency = black_hole_jet_efficiency(bp, props);
  bp->wind_efficiency = black_hole_wind_efficiency(bp, props);

  /* Final jet power at the end of the step */
  jet_power = bp->jet_efficiency * bp->accretion_rate * c * c;

  /* Final luminosity at the end of the step */
  luminosity = bp->radiative_efficiency * bp->accretion_rate * c * c;

  /* The amount of mass the BH is actually swallowing, including the
     effects of the updated efficiencies. */
  const double delta_m_real =
      delta_m_0 * (1. - black_hole_radiative_efficiency(bp, props) -
                   black_hole_jet_efficiency(bp, props) -
                   black_hole_wind_efficiency(bp, props));

  /* Increase the reservoir */
  bp->jet_reservoir += delta_m_0 * c * c * props->eps_f_jet *
                       black_hole_jet_efficiency(bp, props);
  bp->energy_reservoir +=
      delta_m_0 * c * c *
      (props->epsilon_f * black_hole_radiative_efficiency(bp, props) +
       black_hole_wind_efficiency(bp, props));

  /* Update the masses */
  bp->subgrid_mass += delta_m_real;
  bp->total_accreted_mass += delta_m_real;

  /* Update the total accreted masses split by accretion mode of the BHs */
  bp->accreted_mass_by_mode[bp->accretion_mode] += delta_m_real;

  /* Update total energies launched as radiation/winds, and by mode */
  bp->radiated_energy += delta_m_0 * c * c * bp->radiative_efficiency;
  bp->radiated_energy_by_mode[bp->accretion_mode] +=
      delta_m_0 * c * c * bp->radiative_efficiency;
  bp->wind_energy_by_mode[bp->accretion_mode] +=
      delta_m_0 * c * c * bp->wind_efficiency;
  bp->wind_energy += delta_m_0 * c * c * bp->wind_efficiency;

  if (bp->subgrid_mass < 0.) {
    warning(
        "Black hole %lld has reached a negative mass (%f) due"
        " to jet spindown.",
        bp->id, bp->subgrid_mass);
    bp->subgrid_mass = props->subgrid_seed_mass;
  }

  /* Increase the subgrid angular momentum according to what we accreted
   * (already in physical units, a factors from velocity and radius cancel) */
  bp->accreted_angular_momentum[0] +=
      bp->spec_angular_momentum_gas[0] * delta_m_real;
  bp->accreted_angular_momentum[1] +=
      bp->spec_angular_momentum_gas[1] * delta_m_real;
  bp->accreted_angular_momentum[2] +=
      bp->spec_angular_momentum_gas[2] * delta_m_real;

  if (props->use_nibbling && bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative and jet losses, so we must decrease the particle mass
     * in proprtion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -=
        (bp->radiative_efficiency + bp->jet_efficiency + bp->wind_efficiency) *
        bp->accretion_rate * dt;
    if (bp->mass < 0)
      error(
          "Black hole %lld has reached a negative mass (%f). This is "
          "not a great situation, so I am stopping.",
          bp->id, bp->mass);
  }

  /* Below we compute energy required to have a feedback event(s)
   * Note that we have subtracted the particles we swallowed from the ngb_mass
   * and num_ngbs accumulators. */

  /* Update the heating temperature of the BH */
  bp->delta_T = black_hole_feedback_delta_T(bp, props, cosmo, constants);

  /* Mean gas particle mass in the BH's kernel */
  const double mean_ngb_mass = bp->ngb_mass / ((double)bp->num_ngbs);
  /* Energy per unit mass corresponding to the temperature jump delta_T */
  double delta_u = bp->delta_T * props->temp_to_u_factor;
  /* Number of energy injections at this time-step (will be computed below) */
  int number_of_energy_injections;
  /* Average total energy needed to heat the target number of Ngbs */
  const double E_feedback_event =
      delta_u * mean_ngb_mass * props->num_ngbs_to_heat;

  /* Compute and store BH heating-limited time-step */
  if (luminosity > 0.) {
    const float dt_acc = delta_u * mean_ngb_mass * props->num_ngbs_to_heat /
                         (luminosity * props->epsilon_f);
    bp->dt_heat = max(dt_acc, props->time_step_min);
  } else {
    bp->dt_heat = FLT_MAX;
  }

  /* Compute and store BH jet-limited time-step */
  bp->v_jet = black_hole_feedback_dv_jet(bp, props, cosmo, constants);
  const float V_jet = bp->v_jet;
  if (jet_power > 0.) {
    const float dt_acc2 =
        0.5 * mean_ngb_mass * V_jet * V_jet / (jet_power * props->eps_f_jet);
    bp->dt_jet = max(dt_acc2, props->time_step_min);
  } else {
    bp->dt_jet = FLT_MAX;
  }

  /* Are we doing some feedback? */
  if (bp->energy_reservoir > E_feedback_event) {

    /* Probability of heating. Relevant only for the stochastic model. */
    double prob = bp->energy_reservoir / (delta_u * bp->ngb_mass);

    /* Compute the number of energy injections based on probability if and
     * only if we are using the stochastic (i.e. not deterministic) model
     * and the probability prob < 1. */
    if (prob < 1. && !props->AGN_deterministic) {

      /* Initialise counter of energy injections */
      number_of_energy_injections = 0;

      /* How many AGN energy injections will we get? */
      for (int i = 0; i < bp->num_ngbs; i++) {
        const double rand = random_unit_interval_part_ID_and_index(
            bp->id, i, ti_begin, random_number_BH_feedback);
        /* Increase the counter if we are lucky */
        if (rand < prob) number_of_energy_injections++;
      }
    }
    /* Deterministic model or prob >= 1. */
    else {
      /* We want to use up all energy available in the reservoir. Therefore,
       * number_of_energy_injections is > or = props->num_ngbs_to_heat */
      number_of_energy_injections =
          (int)(bp->energy_reservoir / (delta_u * mean_ngb_mass));
    }

    /* Maximum number of energy injections allowed */
    const int N_energy_injections_allowed =
        min(spinjet_blackhole_number_of_rays, bp->num_ngbs);

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
     * If the BH does feedback, the number of AGN events is incremented by one
     */
    bp->AGN_number_of_AGN_events += N_successful_energy_injections > 0;

    /* Update the total (cumulative) energy used for gas heating in AGN feedback
     * by this BH */
    bp->AGN_cumulative_energy += Energy_deposited;

    /* Update the total (cumulative) energies of radiative (thermal) feedback
       split by accretion mode of the BHs */
    bp->thermal_energy_by_mode[bp->accretion_mode] += Energy_deposited;

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
    bp->to_distribute.AGN_delta_u = 0.f;
    bp->to_distribute.AGN_number_of_energy_injections = 0;
  }

  /* Calculate energy required for a jet event (with N_jet particles kicked) */
  const double E_jet_event = 0.5 * V_jet * V_jet * mean_ngb_mass * props->N_jet;

  /* Are we doing some feedback? */
  if (bp->jet_reservoir > E_jet_event) {

    int number_of_jet_injections;
    double delta_u_jet = 0.5 * V_jet * V_jet;

    number_of_jet_injections =
        (int)(bp->jet_reservoir / (0.5 * V_jet * V_jet * mean_ngb_mass));

    /* Limit the number of injections to 2 x N_jet. The jet feedback time
       steps will try to target so that the BH ejects N_jet particles every
       time step, but excesses may build up. Allowing for more than N_jet
       kicks lets these excesses get released.*/
    if (number_of_jet_injections > 2 * props->N_jet) {
      number_of_jet_injections = 2 * props->N_jet;
    } else {
      number_of_jet_injections = props->N_jet;
    }

    /* Compute how much jet energy will be deposited onto the gas */
    /* Note that it will in general be different from E_feedback_event if
     * gas particles are of different mass. */
    double Jet_energy_deposited = 0.0;

    /* Count the number of unsuccessful energy injections (e.g., if the particle
     * that the BH wants to heat has been swallowed and thus no longer exists)
     */
    int N_unsuccessful_jet_injections = 0;

    /* Store all of this in the black hole for delivery onto the gas. */
    bp->to_distribute.AGN_delta_u_jet = delta_u_jet;
    bp->to_distribute.AGN_number_of_jet_injections = number_of_jet_injections;

    /* Check for situations where the gas has been swallowed (not used in the
     * nibbling scheme, or for situations where one of the BH hemispheres is
     * empty. In this case we do not perform any jet kicks at all. */
    for (int i = 0; i < number_of_jet_injections / 2; i++) {

      /* If the gas particle that the BH wants to heat has just been swallowed
       * by the same BH (or if there are no particles in either of the jet
       * jet rays, i.e. on either side of the BH hemisphere), increment the
       * counter of unsuccessful injections. If the particle has not been
       * swallowed by the BH and both hemispheres are populated, the energy
       * to be deposited will be subracted from the jet reservoir. */

      if ((bp->rays_jet[i].id_min_length != -1) &&
          (bp->rays_jet_pos[i].id_min_length != -1)) {

        /* Neither hemisphere of the BH is empty, nor have any particles been
         * swallowed, so we can safely deposity the jet energy. */
        Jet_energy_deposited += delta_u_jet * bp->rays_jet[i].mass;
        Jet_energy_deposited += delta_u_jet * bp->rays_jet_pos[i].mass;

      } else {

        /* Track how many unsuccessful jet injections happened. */
        N_unsuccessful_jet_injections += 2;

        /* Reduce the number of jet injections to be handed out. */
        bp->to_distribute.AGN_number_of_jet_injections -= 2;
      }
    }

    /* Subtract the deposited energy from the BH jet energy reservoir. Note
     * that in the stochastic case, the resulting value might be negative.
     * This happens when (due to the probabilistic nature of the model) the
     * BH injects more energy than it actually has in the reservoir. */
    bp->jet_reservoir -= Jet_energy_deposited;

    /* Total number successful jet energy injections at this time-step. In each
     * energy injection, a certain gas particle from the BH's kernel gets
     * heated. (successful = the particle(s) that is going to get heated by
     * this BH has not been swallowed by the same BH). */
    const int N_successful_jet_injections =
        number_of_jet_injections - N_unsuccessful_jet_injections;

    /* If we will be kicking something, choose a new direction for the jets.
     * This is done by randomly generating a unit vector within a cone of a
     * given opening angle around the BH spin vector. */
    if (N_successful_jet_injections > 0) {
      random_direction_in_cone(
          bp->id, ti_begin, random_number_BH_kick, props->opening_angle,
          bp->angular_momentum_direction, bp->jet_direction);
    }

    /* Increase the number of jet energy injections the black hole has kicked
     * so far.  */
    bp->AGN_number_of_jet_injections += N_successful_jet_injections;

    /* Increase the number of AGN jet events the black hole has had so far.
     * If the BH does feedback, the number of AGN events is incremented by one.
     */
    bp->AGN_number_of_AGN_jet_events += N_successful_jet_injections > 0;

    /* Update the total (cumulative) energy used for gas heating in AGN jet
      feedback by this BH */
    bp->total_jet_energy += Jet_energy_deposited;

    /* Update the total (cumulative) energies of jet feedback split by
       accretion mode of the BHs */
    bp->jet_energy_by_mode[bp->accretion_mode] += Jet_energy_deposited;

    /* Store the time/scale factor when the BH last did AGN jet  feedback */
    if (N_successful_jet_injections) {
      if (with_cosmology) {
        bp->last_AGN_jet_event_scale_factor = cosmo->a;
      } else {
        bp->last_AGN_jet_event_time = time;
      }
    }
  } else {

    /* Flag that we don't want to kick anyone */
    bp->to_distribute.AGN_number_of_jet_injections = 0;
    bp->to_distribute.AGN_delta_u_jet = 0.f;
  }

  /* Decide the accretion mode of the BH, based on the new spin and Eddington
   * fraction */
  black_hole_select_accretion_mode(bp, props);
  bp->accretion_efficiency =
      black_hole_accretion_efficiency(bp, props, constants, cosmo);
  bp->accretion_rate = accr_rate * bp->accretion_efficiency;
  bp->eddington_fraction = bp->accretion_rate / Eddington_rate;

  /* The accretion/feedback mode is now possibly different, so the feedback
     efficiencies need to be updated. This is important for computing the
     next time-step of the BH. */
  bp->radiative_efficiency = black_hole_radiative_efficiency(bp, props);
  bp->jet_efficiency = black_hole_jet_efficiency(bp, props);
  bp->wind_efficiency = black_hole_wind_efficiency(bp, props);

  /* Calculate a BH angular momentum evolution time step. Two conditions are
     used, one ensures that the BH spin changes by a small amount over the
     next time-step, while the other ensures that the spin is not wildly
     redirected between two time-steps */
  if ((fabsf(bp->spin) > 0.01) && (bp->accretion_rate > 0.)) {

    /* How much do we want to allow the spin to change over the next time-
       step? If the spin is large (above 0.1 in magnitude), we allow it
       to only change by 1% of its current value. If it is small, we instead
       use a fixed increment of 0.01. This is to prevent the code from
       becoming very slow. */
    const float spin_magnitude = fabsf(bp->spin);
    float epsilon_spin = 0.01;
    if (spin_magnitude > 0.1) {
      epsilon_spin = 0.01 * spin_magnitude;
    }
    float dt_ang_mom = epsilon_spin /
                       fabsf(black_hole_spinup_rate(bp, constants, props)) *
                       bp->subgrid_mass / bp->accretion_rate;

    /* We now compute the angular-momentum (direction) time-step. We allow
       the anticipated warp angular momentum along the direction perpendicular
       to the current spin vector to be 10% of the current BH angular momentum,
       which should correspond to a change of direction of around 5 degrees. We
       also apply this only if the spin is not low. */
    if (spin_magnitude > 0.1) {

      /* Angular momentum direction of gas in kernel along the direction of
         the current spin vector */

      if (spec_ang_mom_norm > 0.) {
        const float cosine = (bp->spec_angular_momentum_gas[0] *
                                  bp->angular_momentum_direction[0] +
                              bp->spec_angular_momentum_gas[1] *
                                  bp->angular_momentum_direction[1] +
                              bp->spec_angular_momentum_gas[2] *
                                  bp->angular_momentum_direction[2]) /
                             spec_ang_mom_norm;
        /* Compute sine, i.e. the componenent perpendicular to that. */
        const float sine = fmaxf(0., sqrtf(1. - cosine * cosine));

        const float dt_redirection =
            0.1 * black_hole_warp_mass(bp, constants, props) *
            black_hole_angular_momentum_magnitude(bp, constants) /
            (bp->accretion_rate *
             black_hole_warp_angular_momentum(bp, constants, props) * sine);
        dt_ang_mom = min(dt_ang_mom, dt_redirection);
      }
    }

    bp->dt_ang_mom = dt_ang_mom;
  } else {
    bp->dt_ang_mom = FLT_MAX;
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
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the step
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
       * that we have a well-defined re-positioning velocity. We have
       * checked already that reposition_coefficient_upsilon is positive. */
      const float repos_vel =
          props->reposition_coefficient_upsilon *
          pow(bp->subgrid_mass / constants->const_solar_mass,
              props->reposition_exponent_xi);

      const double dx = bp->reposition.delta_x[0];
      const double dy = bp->reposition.delta_x[1];
      const double dz = bp->reposition.delta_x[2];
      const double d = sqrt(dx * dx + dy * dy + dz * dz);

      /* Convert target reposition velocity to a fractional reposition
       * along reposition.delta_x */

      /* Exclude the pathological case of repositioning by zero distance */
      if (d > 0) {
        double repos_frac = repos_vel * dt / d;

        /* We should never get negative repositioning fractions... */
        if (repos_frac < 0)
          error("Wanting to reposition by negative fraction (%g)?", repos_frac);

        /* ... but fractions > 1 can occur if the target velocity is high.
         * We do not want this, because it could lead to overshooting the
         * actual potential minimum. */
        if (repos_frac > 1) repos_frac = 1.;

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
 */
INLINE static void black_holes_create_from_gas(
    struct bpart* bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct part* p, const struct xpart* xp,
    const integertime_t ti_current) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time */
  bp->formation_scale_factor = cosmo->a;

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

  /* Small initial spin magnitude */
  bp->spin = props->seed_spin;

  /* Generate a random unit vector for the spin direction */
  const float rand_cos_theta =
      2. *
      (0.5 - random_unit_interval(bp->id, ti_current, random_number_BH_spin));
  const float rand_sin_theta =
      sqrtf(max(0., (1. - rand_cos_theta) * (1. + rand_cos_theta)));
  const float rand_phi =
      2. * M_PI *
      random_unit_interval(bp->id * bp->id, ti_current, random_number_BH_spin);

  bp->angular_momentum_direction[0] = rand_sin_theta * cos(rand_phi);
  bp->angular_momentum_direction[1] = rand_sin_theta * sin(rand_phi);
  bp->angular_momentum_direction[2] = rand_cos_theta;

  /* Point the jets in the same direction to begin with */
  bp->jet_direction[0] = bp->angular_momentum_direction[0];
  bp->jet_direction[1] = bp->angular_momentum_direction[1];
  bp->jet_direction[2] = bp->angular_momentum_direction[2];

  bp->jet_efficiency = 0.1f;
  bp->radiative_efficiency = 0.1f;
  bp->accretion_disk_angle = 0.01f;
  bp->accretion_mode = BH_thin_disc;
  bp->eddington_fraction = 0.01f;
  bp->jet_reservoir = 0.f;
  bp->total_jet_energy = 0.f;
  bp->wind_efficiency = 0.f;
  bp->wind_energy = 0.f;
  bp->radiated_energy = 0.f;
  bp->accretion_efficiency = 1.f;
  bp->dt_jet = 0.f;
  bp->dt_ang_mom = 0.f;
  bp->delta_T = black_hole_feedback_delta_T(bp, props, cosmo, constants);
  bp->v_jet = black_hole_feedback_dv_jet(bp, props, cosmo, constants);
  bp->AGN_number_of_AGN_jet_events = 0;
  bp->AGN_number_of_jet_injections = 0;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->accreted_mass_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->thermal_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->jet_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->wind_energy_by_mode[i] = 0.f;
  for (int i = 0; i < BH_accretion_modes_count; ++i)
    bp->radiated_energy_by_mode[i] = 0.f;

  /* We haven't accreted anything yet */
  bp->total_accreted_mass = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_time_steps = 0;

  /* We haven't repositioned yet, nor attempted it */
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;

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

  /* First initialisation */
  black_holes_init_bpart(bp);

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

#endif /* SWIFT_SPIN_JET_BLACK_HOLES_H */
