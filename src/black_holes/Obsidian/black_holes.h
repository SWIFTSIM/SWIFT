/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_OBSIDIAN_BLACK_HOLES_H
#define SWIFT_OBSIDIAN_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "cooling.h"
#include "cosmology.h"
#include "dimension.h"
#include "exp10.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"
#include "star_formation.h"

/* Standard includes */
#include <float.h>
#include <gsl/gsl_poly.h>
#include <math.h>

/**
 * @brief How much of the feedback actually couples to the medium?
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param BH_state The current state of the black hole.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double get_black_hole_coupling(
    const struct bpart *const bp, const struct black_holes_props *props,
    const struct cosmology *cosmo, const struct phys_const *phys_const) {
  const int BH_state = bp->state;
  switch (BH_state) {
    case BH_states_adaf: {
      float scaling = 1.f;
      if (props->adaf_coupling < 0.f) {
        min(pow(1. + cosmo->z, props->adaf_z_scaling), 1.);
      }
      return fabs(props->adaf_coupling) * scaling;
      break;
    }
    case BH_states_quasar: {
      float quasar_coupling = fabs(props->quasar_coupling);
      const double c = phys_const->const_speed_light_c;
      const double luminosity =
          bp->radiative_efficiency * bp->accretion_rate * c * c;
      const float mass_limit = get_black_hole_adaf_mass_limit(bp, props, cosmo);
      if (luminosity > props->quasar_luminosity_thresh &&
          props->quasar_luminosity_thresh > 0.f &&
          bp->subgrid_mass * props->mass_to_solar_mass > mass_limit) {
        quasar_coupling = fmin(
            quasar_coupling * luminosity / props->quasar_luminosity_thresh, 1.);
        // message("BOOST: z=%g id=%lld MBH=%g LBH=%g boost=%g qcoupling=%g",
        // cosmo->z, bp->id, bp->subgrid_mass * props->mass_to_solar_mass,
        // luminosity * props->conv_factor_energy_rate_to_cgs, luminosity /
        // props->quasar_luminosity_thresh, quasar_coupling);
      }
      return quasar_coupling;
      break;
    }
    case BH_states_slim_disk:
      return fabs(props->slim_disk_coupling);
      break;
    default:
      error("Invalid black hole state.");
      return 0.;
      break;
  }
}

/**
 * @brief Computes the radiative efficiency in the slim disk mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static double
get_black_hole_slim_disk_efficiency(const struct black_holes_props *props,
                                    const double f_Edd) {
  if (f_Edd <= 0.) return 0.;
  const double R = 1. / f_Edd;
  /* Efficiency from Lupi et al. (2014),
   * super eddington accretion and feedback */
  return (R / 16.) * props->A_sd *
         (0.985 / (R + (5. / 8.) * props->B_sd) +
          0.015 / (R + (5. / 8.) * props->C_sd));
}

/**
 * @brief Computes the radiative efficiency in the ADAF mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static double
get_black_hole_adaf_efficiency(const struct black_holes_props *props,
                               const double f_Edd) {
  return props->epsilon_r * f_Edd; /* scales with M_dot,BH */
}

/**
 * @brief Chooses and calls the proper radiative efficiency function for the
 *        state.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd The accretion rate over the Eddington rate.
 * @param BH_state The current state of the BH.
 */
__attribute__((always_inline)) INLINE static double
get_black_hole_radiative_efficiency(const struct black_holes_props *props,
                                    const double f_Edd, const int BH_state) {
  switch (BH_state) {
    case BH_states_adaf:
      return get_black_hole_adaf_efficiency(props, f_Edd);
    case BH_states_quasar:
      return props->epsilon_r;
    case BH_states_slim_disk:
      return get_black_hole_slim_disk_efficiency(props, f_Edd);
    default:
      error("Invalid black hole state.");
      break;
  }

  return 0.;
}

/**
 * @brief Compute the wind launch speed for this feedback step.
 *
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const (in internal units).
 * @param bp The black hole particle.
 */
__attribute__((always_inline)) INLINE static double get_black_hole_wind_speed(
    const struct black_holes_props *props, const struct phys_const *phys_const,
    const struct bpart *bp) {

  if (bp->accretion_rate < 0.f || bp->m_dot_inflow < 0.f) return 0.f;

  float v_kick = 0.f;
  if (props->quasar_wind_speed < 0.f || props->slim_disk_wind_speed < 0.f) {
    const float subgrid_mass_Msun =
        bp->subgrid_mass * props->mass_to_solar_mass;

    if (bp->subgrid_mass > props->subgrid_seed_mass) {
      const float min_BH_mass_Msun =
          props->minimum_black_hole_mass_v_kick * props->mass_to_solar_mass;
      const float dlog10_BH_mass =
          log10f(subgrid_mass_Msun) - log10f(min_BH_mass_Msun);
      v_kick = 500.f + (500.f / 3.f) * dlog10_BH_mass;

      /* Sometimes can get very small leading to huge mass loadings */
      if (v_kick < props->minimum_v_kick_km_s) {
        v_kick = props->minimum_v_kick_km_s;
      }

      v_kick *= props->kms_to_internal;
    }
  }

  switch (bp->state) {
    case BH_states_adaf:
      return fabs(props->adaf_wind_speed);
      break;
    case BH_states_quasar:
      if (props->quasar_wind_speed < 0.f && v_kick > 0.f) {
        return v_kick;
      } else {
        return fabs(props->quasar_wind_speed);
      }
      break;
    case BH_states_slim_disk:
      if (props->slim_disk_wind_speed < 0.f && v_kick > 0.f) {
        return v_kick;
      } else {
        return fabs(props->slim_disk_wind_speed);
      }
      break;
    default:
      error("Invalid black hole state.");
      return 0.f;
      break;
  }
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const (in internal units).
 * @param cosmo The current cosmological model.
 * @param bp The black hole particle.
 * @param m_dot_inflow_m_dot_edd M_dot,inflow scaled to M_dot,Edd for the BH.
 */
__attribute__((always_inline)) INLINE static double
get_black_hole_upper_mdot_medd(const struct black_holes_props *props,
                               const struct phys_const *phys_const,
                               const struct cosmology *cosmo,
                               const struct bpart *const bp,
                               const double m_dot_inflow_m_dot_edd) {

  if (m_dot_inflow_m_dot_edd <= 0.) return 0.;

  double x1, x2, x3;
  double a3, a2, a1, a0;

  double phi = props->slim_disk_phi;
  if (props->slim_disk_wind_speed < 0.f) {
    const float v_kick = get_black_hole_wind_speed(props, phys_const, bp);

    if (v_kick > 0.f) {
      /* Set the slim disk mass loading to be continuous at the
       * eta upper boundary. Compute the phi term to solve for the
       * accretion fraction.
       * WARNING: Using the quasar_wind_momentum_flux because that is the
       * default way of using the model. Any changes to that require a
       * change here. */
      const double c_over_v = phys_const->const_speed_light_c / v_kick;

      /* TODO: Compute once at the beginning of the simulation */
      double mom_flux_times_epsilon_sd =
          props->quasar_wind_momentum_flux *
          get_black_hole_coupling(bp, props, cosmo, phys_const);
      phi = mom_flux_times_epsilon_sd * c_over_v;
    }
  }

  int num_roots;

  a3 = ((5. * 5.) / (8. * 8.)) * props->B_sd * props->C_sd;
  a2 =
      (5. / 8.) *
      ((props->B_sd + props->C_sd) +
       (phi / 16.) * props->A_sd * (0.015 * props->B_sd + 0.985 * props->C_sd) -
       (5. / 8.) * props->B_sd * props->C_sd * m_dot_inflow_m_dot_edd);
  a1 = 1. + (phi / 16.) * props->A_sd -
       (5. / 8.) * (props->B_sd + props->C_sd) * m_dot_inflow_m_dot_edd;
  a0 = -m_dot_inflow_m_dot_edd;

  a2 /= a3;
  a1 /= a3;
  a0 /= a3;

  num_roots = gsl_poly_solve_cubic(a2, a1, a0, &x1, &x2, &x3);
  if (num_roots == 1) {
    if (x1 >= 0.) {
      return x1;
    } else {
      if (m_dot_inflow_m_dot_edd > 1.e-6)
      warning(
          "num_roots=1 m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g "
          "a1=%g a0=%g",
          m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
      return 0.;
    }
  }
  if (x3 >= 0.) {
    return x3;
  } else {
    if (m_dot_inflow_m_dot_edd > 1.e-6)
    warning(
        "num_roots=0 m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g a1=%g "
        "a0=%g",
        m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
    return 0.;
  }

  return 0.;
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const (in internal units).
 * @param m_dot_inflow M_dot,inflow in internal units.
 * @param BH_mass The subgrid mass of the BH in internal units.
 * @param BH_state The current state of the BH.
 * @param Eddington_rate M_dot,Edd in internal units.
 */
__attribute__((always_inline)) INLINE static double
get_black_hole_accretion_factor(const struct black_holes_props *props,
                                const struct phys_const *phys_const,
                                const struct cosmology *cosmo,
                                const struct bpart *const bp,
                                const double Eddington_rate) {

  const double m_dot_inflow = bp->m_dot_inflow;
  const double BH_mass = bp->subgrid_mass;
  const int BH_state = bp->state;

  if (m_dot_inflow <= 0. || BH_mass <= 0.) return 0.;

  switch (BH_state) {
    case BH_states_adaf:
      return props->adaf_f_accretion;
      break;
    case BH_states_quasar: {
      float v_kick = 0.f;
      float f_accretion = 0.f;
      if (props->quasar_wind_speed < 0.f) {
        /* Save computation by only computing when specified by the user */
        v_kick = get_black_hole_wind_speed(props, phys_const, bp);
        if (v_kick > 0.) {
          const double c = phys_const->const_speed_light_c;
          const double c_over_v = c / v_kick;
          double quasar_coupling =
              get_black_hole_coupling(bp, props, cosmo, phys_const);
          double quasar_wind_mass_loading = props->quasar_wind_momentum_flux *
                                            quasar_coupling * props->epsilon_r *
                                            c_over_v;
          f_accretion = 1.f / (1.f + quasar_wind_mass_loading);
        }
      }

      if (f_accretion > 0.f) {
        return f_accretion;
      } else {
        return props->quasar_f_accretion;
      }
      break;
    }
    case BH_states_slim_disk: {
      /* This is the FRACTION of the total so divide by M_dot,inflow */
      const double f_edd = m_dot_inflow / Eddington_rate;
      double mdot_medd =
          get_black_hole_upper_mdot_medd(props, phys_const, cosmo, bp, f_edd);
      return mdot_medd * Eddington_rate / m_dot_inflow;
      break;
    }
    default:
      error("Invalid black hole state.");
      return 0.;
      break;
  }
}

/**
 * @brief Computes the time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart *const bp, const struct black_holes_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo) {

  /* Allow for finer timestepping if necessary! */
  float dt_accr = FLT_MAX;
  float dt_overall = FLT_MAX;
  float dt_kick = FLT_MAX;
  const float min_subgrid_mass = props->minimum_black_hole_mass_unresolved;

  /* Only limit when in the resolved feedback regime */
  if (bp->accretion_rate > 0.f && bp->subgrid_mass > min_subgrid_mass) {
    dt_accr = props->dt_accretion_factor * bp->mass / bp->accretion_rate;

    if (bp->state == BH_states_adaf && bp->jet_mass_loading > 0.f) {
      dt_kick = bp->ngb_mass / (bp->jet_mass_loading * bp->accretion_rate);
    } else {
      if (bp->f_accretion > 0.f) {
        /* Make sure that the wind mass does not exceed the kernel gas mass */
        const float psi = (1.f - bp->f_accretion) / bp->f_accretion;
        dt_kick = bp->ngb_mass / (psi * bp->accretion_rate);
      }
    }

    dt_overall = min(dt_kick, dt_accr);
  }

  if (dt_overall < props->time_step_min) {
    message(
        "Warning! BH_TIMESTEP_LOW: id=%lld (%g Myr) is below time_step_min (%g "
        "Myr).",
        bp->id, dt_overall * props->time_to_Myr,
        props->time_step_min * props->time_to_Myr);
  }

  return max(dt_overall, props->time_step_min);
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
    struct bpart *bp, const struct black_holes_props *props) {

  bp->time_bin = 0;
  if (props->use_subgrid_mass_from_ics == 0) {
    bp->subgrid_mass = bp->mass;
  } else if (props->with_subgrid_mass_check && bp->subgrid_mass <= 0) {
    error(
        "Black hole %lld has a subgrid mass of %f (internal units).\n"
        "If this is because the ICs do not contain a 'SubgridMass' data "
        "set, you should set the parameter "
        "'ObsidianAGN:use_subgrid_mass_from_ics' to 0 to initialize the "
        "black hole subgrid masses to the corresponding dynamical masses.\n"
        "If the subgrid mass is intentionally set to this value, you can "
        "disable this error by setting 'ObsidianAGN:with_subgrid_mass_check' "
        "to 0.",
        bp->id, bp->subgrid_mass);
  }
  bp->total_accreted_mass = 0.f;
  bp->accretion_disk_mass = 0.f;
  bp->accretion_rate = 0.f;
  bp->bondi_accretion_rate = 0.f;
  bp->formation_time = -1.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->number_of_time_steps = 0;
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;
  bp->accreted_angular_momentum[0] = 0.f;
  bp->accreted_angular_momentum[1] = 0.f;
  bp->accreted_angular_momentum[2] = 0.f;
  bp->last_repos_vel = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->delta_energy_this_timestep = 0.f;
  bp->state = BH_states_slim_disk;
  bp->radiative_efficiency = 0.f;
  bp->f_accretion = 0.f;
  bp->m_dot_inflow = 0.f;
  bp->cold_disk_mass = 0.f;
  bp->jet_mass_reservoir = 0.f;
  /* Default to the original value at fixed jet_velocity */
  bp->jet_mass_loading = props->jet_mass_loading;
  bp->jet_mass_kicked_this_step = 0.f;
  bp->adaf_energy_to_dump = 0.f;
  bp->adaf_energy_used_this_step = 0.f;
}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart *bp) {

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
  bp->gas_SFR = 0.f;
  bp->hot_gas_mass = 0.f;
  bp->cold_gas_mass = 0.f;
  bp->cold_disk_mass = 0.f;
  bp->hot_gas_internal_energy = 0.f;
  bp->sound_speed_subgrid_gas = -1.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->circular_velocity_gas[0] = 0.f;
  bp->circular_velocity_gas[1] = 0.f;
  bp->circular_velocity_gas[2] = 0.f;
  bp->angular_momentum_gas[0] = 0.f;
  bp->angular_momentum_gas[1] = 0.f;
  bp->angular_momentum_gas[2] = 0.f;
  bp->stellar_mass = 0.f;
  bp->stellar_bulge_mass = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->ngb_mass = 0.f;
  bp->gravitational_ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->num_gravitational_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f;       /* Optionally accumulated ngb-by-ngb */
  bp->bondi_accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */
  bp->m_dot_inflow = 0.f;               /* reset accretion rate */
  bp->kernel_wt_sum = 0.f;

  /* update the reservoir */
  bp->jet_mass_reservoir -= bp->jet_mass_kicked_this_step;
  bp->jet_mass_kicked_this_step = 0.f;
  if (bp->jet_mass_reservoir < 0.f) {
    bp->jet_mass_reservoir = 0.f; /* reset reservoir if used up */
  }
  /* update the unresolved reservoir */
  bp->unresolved_mass_reservoir -= bp->unresolved_mass_kicked_this_step;
  bp->unresolved_mass_kicked_this_step = 0.f;
  if (bp->unresolved_mass_reservoir < 0.f) {
    bp->unresolved_mass_reservoir = 0.f;
  }
  /* update the adaf energy reservoir */
  if (bp->adaf_wt_sum > 0.f) {
    const double adaf_energy_used =
        bp->adaf_energy_used_this_step / bp->adaf_wt_sum;
    bp->adaf_energy_to_dump -= adaf_energy_used;
    bp->adaf_wt_sum = 0.f;
    bp->adaf_energy_used_this_step = 0.f;
    if (bp->adaf_energy_to_dump < 0.f) {
      bp->adaf_energy_to_dump = 0.f;
    }
  } else {
    bp->adaf_energy_used_this_step = 0.f;
  }
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
    struct bpart *restrict bp, float dt_drift) {

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
black_holes_reset_predicted_values(struct bpart *bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart *bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart *bp, const struct cosmology *cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* --- Finish the calculation by inserting the missing h factors --- */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  float rho_inv = 1.f;
  if (bp->rho_gas > 0.f) rho_inv = 1.f / bp->rho_gas;

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  bp->internal_energy_gas *= h_inv_dim * rho_inv;

  /* Non-weighted (no decoupled winds) properties below.
   * All mass-weighted quantities are for the hot & cold gas */
  float m_hot_inv = 1.f;
  if (bp->hot_gas_mass > 0.f) m_hot_inv /= bp->hot_gas_mass;
  /* Or the total mass */
  float m_tot_inv = 1.f;
  if (bp->ngb_mass > 0.f) m_tot_inv /= bp->ngb_mass;

  bp->hot_gas_internal_energy *= m_hot_inv;
  bp->velocity_gas[0] *= m_tot_inv;
  bp->velocity_gas[1] *= m_tot_inv;
  bp->velocity_gas[2] *= m_tot_inv;
  bp->circular_velocity_gas[0] *= m_tot_inv;
  bp->circular_velocity_gas[1] *= m_tot_inv;
  bp->circular_velocity_gas[2] *= m_tot_inv;

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv). It is now a VELOCITY.
   */
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
black_holes_bpart_has_no_neighbours(struct bpart *bp,
                                    const struct cosmology *cosmo) {

  // warning(
  //     "BH particle with ID %lld treated as having no neighbours (h: %g, "
  //     "wcount: %g).",
  //     bp->id, bp->h, bp->density.wcount);

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
  bp->hot_gas_internal_energy = -FLT_MAX;
}

/**
 * @brief Return the current instantaneous accretion rate of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accretion_rate(const struct bpart *bp) {
  return bp->accretion_rate;
}

/**
 * @brief Return the total accreted gas mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accreted_mass(const struct bpart *bp) {
  return bp->total_accreted_mass;
}

/**
 * @brief Return the subgrid mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_subgrid_mass(const struct bpart *bp) {
  return bp->subgrid_mass;
}

/**
 * @brief Return the current bolometric luminosity of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart *bp,
                                      const struct phys_const *phys_const) {
  const double c = phys_const->const_speed_light_c;
  return bp->accretion_rate * bp->radiative_efficiency * c * c;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart *bp, const struct phys_const *phys_const) {
  const double c = phys_const->const_speed_light_c;
  /* accretion_rate is M_dot,acc from the paper */
  return bp->radiative_efficiency * bp->accretion_rate * c * c;
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
    struct bpart *bp, const struct part *p, const struct xpart *xp,
    const struct cosmology *cosmo) {

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
  struct chemistry_bpart_data *bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data *p_chem = &p->chemistry_data;
  chemistry_add_part_to_bpart(bp_chem, p_chem, gas_mass);

  /* This BH swallowed a gas particle */
  bp->number_of_gas_swallows++;
  bp->number_of_direct_gas_swallows++;

  /* This BH lost a neighbour */
  bp->num_ngbs--;
  bp->num_gravitational_ngbs--;
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
    struct bpart *bpi, const struct bpart *bpj, const struct cosmology *cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props *props,
    const struct phys_const *phys_const) {

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
  struct chemistry_bpart_data *bpi_chem = &bpi->chemistry_data;
  const struct chemistry_bpart_data *bpj_chem = &bpj->chemistry_data;
  chemistry_add_bpart_to_bpart(bpi_chem, bpj_chem);

  /* Update the energy reservoir */
  bpi->jet_mass_reservoir += bpj->jet_mass_reservoir;

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
 * @brief Function to generate a random number from a Gaussian distribution.
 * @param mu Mean of Gaussian
 * @param sigma Standard deviation of Gaussian
 * @param u1 Random number in (0,1)
 * @param u2 Random number in (0,1)
 */
__attribute__((always_inline)) INLINE static float gaussian_random_number(
    float mu, float sigma, double u1, double u2) {
  double mag, z0, z1;

  /* Apply the Box-Muller transform */
  mag = sigma * sqrt(-2.0 * log(u1));
  z0 = mag * cos(2.0 * M_PI * u2) + mu;
  z1 = mag * sin(2.0 * M_PI * u2) + mu;
  if (u1 + u2 < 1.f) {
    return z0;
  }
  return z1;
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantities
 * required for the feedback loop.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy fllor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 * @param ti_begin The time at which the step begun (ti_current).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart *restrict bp, const struct black_holes_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct entropy_floor_properties *floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {

  /* Record that the black hole has another active time step */
  bp->number_of_time_steps++;

  if (dt == 0. || bp->rho_gas == 0. || bp->h == 0.) return;

  /* Collect information about galaxy that the particle belongs to */
  const float galaxy_mstar = bp->galaxy_data.stellar_mass;
  const float galaxy_mgas = bp->galaxy_data.gas_mass;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (galaxy_mstar <= 0.f) return;

  /* Gather some physical phys_const (all in internal units) */
  const double G = phys_const->const_newton_G;
  const double c = phys_const->const_speed_light_c;
  const double proton_mass = phys_const->const_proton_mass;
  const double sigma_Thomson = phys_const->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd_maximum = props->f_Edd_maximum;

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->subgrid_mass)) error("subgrid_mass nan");
#endif
  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Compute the Eddington rate (internal units).
   * IMPORTANT: epsilon_r = 0.1 is the SET value for the Eddington rate.
   * It is assumed in all of the derivations for the state model.
   */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (0.1 * c * sigma_Thomson);

  const double bh_h = kernel_gamma * bp->h;
  const double bh_h_inv = 1. / bh_h;
  const double volume_bh_inv =
      (3. / (4. * M_PI)) * bh_h_inv * bh_h_inv * bh_h_inv;
  double gas_rho = 0.;
  if (props->bondi_use_all_gas) {
    gas_rho = bp->rho_gas;
  } else {
    gas_rho = bp->hot_gas_mass * volume_bh_inv;
  }

  const double gas_rho_phys = gas_rho * cosmo->a3_inv;

  /* We can now compute the Bondi accretion rate (internal units)
   * D. Rennehan: In Simba, we only consider the hot gas within
   * the kernel for the Bondi rate, and the cold gas using the
   * torque accretion estimator.
   */
  double Bondi_rate = 0.;

  /* Use all gas within the kernel or only the hot gas*/
  double gas_internal_energy = 0.;
  if (props->bondi_use_all_gas) {
    if (bp->internal_energy_gas > 0.) {
      gas_internal_energy = bp->internal_energy_gas;
    }
  } else {
    if (bp->hot_gas_internal_energy > 0.) {
      gas_internal_energy = bp->hot_gas_internal_energy;
    }
  }

  /* Check if there is hot/any gas in the kernel */
  if (gas_internal_energy > 0.) {
    double gas_c =
        gas_soundspeed_from_internal_energy(gas_rho, gas_internal_energy);

    if (gas_c > 0.) {
      const double gas_c_phys_inv = 1. / (cosmo->a_factor_sound_speed * gas_c);

      float BH_mass_bondi = BH_mass;
      if (BH_mass_bondi > props->bondi_BH_mass_cap) 
	BH_mass_bondi = props->bondi_BH_mass_cap;

      Bondi_rate = 4. * M_PI * G * G * BH_mass_bondi * BH_mass_bondi * gas_rho_phys *
                   gas_c_phys_inv * gas_c_phys_inv * gas_c_phys_inv;

      /* In the case of standard Bondi, we limit it to the Eddington rate */
      Bondi_rate =
          fmin(Bondi_rate, props->f_Edd_Bondi_maximum * Eddington_rate);
    }
  }

  /* The accretion rate estimators give Mdot,inflow
   * (Mdot,BH = f_acc * Mdot,inflow) */
  const double bondi_accr_rate = props->bondi_alpha * Bondi_rate;

  /* Compute the torque-limited accretion rate */
  double torque_accr_rate = 0.;

  double f_corr_stellar = 10.;
  if (galaxy_mgas > 0.) {
    f_corr_stellar = min(galaxy_mstar / galaxy_mgas, f_corr_stellar);
  }

  /* Torque accretion rate based on some fraction of gas near BH
   * falling in on dynamical time. (This is the default.)
   * Here the accretion rate is only based on Mgas / tdyn.
   * We do not use the DM mass to compute tdyn since it probably
   * doesn't contribute much near the core of the system. We also
   * assume that the gas fraction is constant in the galaxy in
   * order to compute Mstar within the kernel of the black hole.
   * Therefore, Mdot = Mgas / tdyn = Mgas / sqrt(3pi/(32 G rho))
   * and rho = (Mgas + Mstar + Mdm) / (4pi h^3 / 3) where
   * Mstar = Mgas / fgas, Mdm = 0. Therefore,
   * rho = 3 * ((1 + fgas) / fgas) * Mgas / (4 * pi * h^3)
   * and
   * Mdot = Mgas * sqrt(32 * G * 3 * ((1 + fgas) / fgas) * Mgas)) /
   *    sqrt(3 * pi * 4 * pi * h^3)
   *      = sqrt(96 * G * ((1 + fgas) / fgas) * Mgas^3) /
   *    sqrt(12 * pi^2 * h^3)
   *      = (1 / pi) * sqrt(8 * G * ((1 + fgas) / fgas) * (Mgas / h)^3))
   */
  float tdyn_inv = FLT_MAX;
  const float potential = fabs(gravity_get_comoving_potential(bp->gpart));
  /* Includes dynamical mass of the BH */
  switch (props->dynamical_time_calculation_method) {
    /* Assume gas fraction is the same in the kernel and outside */
    case 0: {
      /* Compute correction to total dynamical mass around
       * BH contributed by stars */
      const float m_star_gal = galaxy_mstar;
      const float m_gas_cold_gal = galaxy_mgas;
      const float m_gas_bh = bp->gravitational_ngb_mass;
      const float m_bh = bp->mass;

      /* Compute stellar mass assuming a constant cold gas fraction in the
       * entire galaxy. If m_gas_cold_bh is zero it doesn't matter since
       * the BH won't acrrete in the torque mode anyway. */
      const float m_gas_cold_bh = bp->cold_gas_mass;
      float m_star_bh = 0.;
      if (m_gas_cold_gal > 0.) {
        m_star_bh = fmin(m_star_gal / m_gas_cold_gal, 10.f) * m_gas_cold_bh;
      }

      /* Have to assume baryon dominance within the kernel */
      const float rho_est = (m_star_bh + m_gas_bh + m_bh) * volume_bh_inv;

      /* Inverse physical dynamical time */
      tdyn_inv = sqrt(32. * G * rho_est * cosmo->a3_inv / (3. * M_PI));
      break;
    }

    /* Assume BH potential */
    case 1:
      if (potential >= 0.f) {
        tdyn_inv = (sqrt(potential) / bh_h) * cosmo->a2_inv;
      }
      break;

    /* Assume dynamical time from the kernel mass */
    case 2: {
      /* do not have gravity_props here */
      const float hsml = kernel_gamma * bh_h;
      const float volume = (4. * M_PI / 3.) * hsml * hsml * hsml;
      const float rho = bp->mass / volume;
      tdyn_inv = sqrt(32. * G * rho * cosmo->a3_inv / (3. * M_PI));
      break;
    }

    default:
      error("Unknown dynamical time calculation method %d",
            props->dynamical_time_calculation_method);
      break;
  }

  /* Compute inverse of accretion time into BH */
  tdyn_inv /= props->dynamical_time_factor;

  /* Limit by max dynamical time, with z=0 value scaled by H0/H */
  if (props->dynamical_time_max > 0.) {
    const float t_inv_min = cosmo->H / (props->dynamical_time_max * cosmo->H0);
    tdyn_inv = fmax(tdyn_inv, t_inv_min);
  }

  /* Create a spread in accretion times, with minimum at
   * free-fall time=0.5*tdyn */
  const float tdyn_sigma = props->tdyn_sigma;
  if (tdyn_sigma > 0.f) {
    const double ran1 =
        random_unit_interval(bp->id, ti_begin, random_number_BH_swallow);
    const double ran2 =
        random_unit_interval(bp->id, ti_begin, random_number_BH_swallow);
    const float gaussian_random =
        gaussian_random_number(0.f, tdyn_sigma, ran1, ran2);
    tdyn_inv /= 0.5 * (1.f + fabs(gaussian_random));
  }

  float f_suppress = 1.f;

  const float corot_gas_mass =
      bp->cold_gas_mass - 2. * (bp->cold_gas_mass - bp->cold_disk_mass);
  float torque_norm = fabs(props->torque_accretion_norm);
  if (props->torque_accretion_norm < 0.f && bp->subgrid_mass * props->mass_to_solar_mass > 1.e9) {
    torque_norm *= exp(-bp->subgrid_mass * props->mass_to_solar_mass * 1.e-9) * 2.71828; // suppress torque accr
  }
  if (torque_norm > 0.f) {
    switch (props->torque_accretion_method) {
      case 0:
        if (galaxy_mgas > 0.) {
          torque_accr_rate =
              torque_norm * bp->cold_disk_mass * tdyn_inv;
        }
        break;

      case 1:
        if (corot_gas_mass > 0. && bp->cold_gas_mass > 0.) {
          torque_accr_rate =
              torque_norm * corot_gas_mass * tdyn_inv;
        }
        break;

      case 2:
        if (corot_gas_mass > 0. && bp->cold_gas_mass > 0.) {
          const float m_disk = bp->cold_gas_mass * f_corr_stellar;
          const float f_disk = corot_gas_mass / bp->cold_gas_mass;

          const float r0 = bh_h * cosmo->a * (props->length_to_parsec * 0.01);

          const float alpha = 5.;
          const float mass_to_1e9solar = props->mass_to_solar_mass * 1.0e-9;
          const float mass_to_1e8solar = props->mass_to_solar_mass * 1.0e-8;

          const float f0 =
              0.31 * f_disk * f_disk * pow(m_disk * mass_to_1e9solar, -1. / 3.);
          const float f_gas = corot_gas_mass / m_disk;
          const float mass_in_1e8solar = BH_mass * mass_to_1e8solar;

          torque_accr_rate = torque_norm * alpha *
                             corot_gas_mass * mass_to_1e9solar *
                             powf(f_disk, 5. / 2.) *
                             powf(mass_in_1e8solar, 1. / 6.) *
                             powf(r0, -3. / 2.) / (1. + f0 / f_gas);
          torque_accr_rate *= (props->time_to_yr / props->mass_to_solar_mass);
        }
        break;

      case 3:
        if (galaxy_mgas > 0.) {
          torque_accr_rate =
              torque_norm * bp->cold_gas_mass * tdyn_inv;
        }
        break;

      default:
        error("Unknown torque_accretion_method=%d",
              props->torque_accretion_method);
        break;
    }

    /* Do suppression of BH growth */
    switch (props->suppress_growth) {
      case 1: {
        const double r0 = bh_h * cosmo->a * props->length_to_parsec;
        const double sigma_eff = f_corr_stellar * bp->ngb_mass *
                                 props->mass_to_solar_mass / (M_PI * r0 * r0);
	
        f_suppress = sigma_eff / (sigma_eff + 3000.);
        break;
      }

      case 2:
      case 6:
      case 7: {
        double m_suppress = fabs(props->bh_characteristic_suppression_mass);
        if (props->bh_characteristic_suppression_mass < 0) {
          m_suppress *= cosmo->a;
        }

        f_suppress =
            1. - exp(-BH_mass * props->mass_to_solar_mass / m_suppress);
        break;
      }

      case 4:
      case 5: {
        /* compute mass loading factor from SF feedback,
         * should be same as used in feedback_mass_loading_factor()
         */
        const float galaxy_stellar_mass = galaxy_mstar;
        const float eta_norm = props->FIRE_eta_normalization;
        const float eta_break = props->FIRE_eta_break;
        const float eta_lower_slope = props->FIRE_eta_lower_slope;
        const float eta_upper_slope = props->FIRE_eta_upper_slope;
        const float eta_lower_slope_EOR = props->FIRE_eta_lower_slope;
        const float eta_minmass = props->minimum_galaxy_stellar_mass;
        const float eta_suppress = props->wind_eta_suppression_redshift;

        const double eta = feedback_mass_loading_factor(
            cosmo, galaxy_stellar_mass, eta_minmass, eta_norm, eta_break,
            eta_lower_slope, eta_upper_slope, eta_lower_slope_EOR,
            eta_suppress);

        if (bp->cold_gas_mass * tdyn_inv > 0.f) {
          /* star formation efficiency, frac of gas converted
           * to stars per tdyn */
          float sf_eff = props->suppression_sf_eff;
          if (sf_eff < 0.f) {
            /* SF efficiency within BH kernel. Cap at cloud-scale SFE from
             * Leroy+25 */
            sf_eff = fmin(bp->gas_SFR / (tdyn_inv * bp->cold_gas_mass),
                          fabs(sf_eff));
          }

          /* Suppresses accretion by factor accounting for mass
           * lost in outflow over accretion time. ODE:
           * dM/dt = -eta * sf_eff * M / tdyn */
          f_suppress = exp(-eta * sf_eff);
          // message("BH_SUPPRESS: z=%g id=%lld M*=%g eta=%g eff=%g tfac=%g
          // fsupp=%g", cosmo->z, bp->id, galaxy_stellar_mass *
          // props->mass_to_solar_mass, eta, sf_eff, t_accrete * tdyn_inv,
          // exp(-eta * sf_eff * t_accrete * tdyn_inv));
        }
        break;
      }

      default:
        break;
    }
  }

  /* Apply suppression factor to torque accretion */
  torque_accr_rate *= f_suppress;

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bondi_accr_rate)) error("bondi_accr_rate nan");
  if (isnan(torque_accr_rate)) error("torque_accr_rate nan");
#endif

  /* Right now this is M_dot,inflow. We will multiply by
   * f_accretion later to make it M_dot,acc */
  bp->accretion_rate = bondi_accr_rate + torque_accr_rate;

  /* We will use eddington_fraction_lower_boundary and
   * eddington_fraction_upper_boundary to divide up the accretion rate
   * in three regimes.
   *
   * In order to switch out of a regime (i.e. a state), it is necessary
   * for the true accretion rate (compared to Eddington rate) to switch
   * over the boundary. Therefore, before we switch a state we must calculate
   * what the previous state predicts the true accretion rate onto the SMBH is,
   * and then update the state if it crosses a boundary.
   */

  /* We need to store the full M_dot,inflow rate to calculate the
   * fraction at high accretion rate */
  bp->m_dot_inflow = bp->accretion_rate;
  const double f_accretion = get_black_hole_accretion_factor(
      props, phys_const, cosmo, bp, Eddington_rate);
  double predicted_mdot_medd =
      bp->accretion_rate * f_accretion / Eddington_rate;
  const float my_adaf_mass_limit =
      get_black_hole_adaf_mass_limit(bp, props, cosmo);

  /* Switch between states depending on f_edd */
  switch (bp->state) {
    case BH_states_adaf:
      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
        break;
      }
      if (predicted_mdot_medd > props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case ADAF */
    case BH_states_quasar:
      if (BH_mass > my_adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
      }

      break; /* end case quasar */
    case BH_states_slim_disk:
      if (BH_mass > my_adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd < props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case slim disk */
    default:
      error("Invalid black hole state.");
      break;
  }

  /* This depends on the new state */
  bp->f_accretion = get_black_hole_accretion_factor(props, phys_const, cosmo,
                                                    bp, Eddington_rate);
#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->f_accretion)) error("f_accretion nan");
#endif
  if (bp->f_accretion <= 1.e-10f) bp->f_accretion = 0.f;

  /* bp->accretion_rate is M_dot,acc in Rennehan+24 */
  bp->accretion_rate *= bp->f_accretion;

  /* Track Bondi accretion separately for diagnostics (remainder is torque) */
  bp->bondi_accretion_rate = bondi_accr_rate * bp->f_accretion;

  if (!props->bondi_use_all_gas) {
    /* Now we can Eddington limit. */
    bp->accretion_rate =
        min(bp->accretion_rate, f_Edd_maximum * Eddington_rate);
  }

  /* All accretion is done, now we can set the eddington fraction */
  bp->eddington_fraction = bp->accretion_rate / Eddington_rate;

  /* Get the new radiative efficiency based on the new state */
  bp->radiative_efficiency = get_black_hole_radiative_efficiency(
      props, bp->eddington_fraction, bp->state);
#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(bp->radiative_efficiency)) error("radiative_efficiency nan");
#endif
  if (bp->radiative_efficiency < 1.e-10f) bp->radiative_efficiency = 0.f;

  double mass_rate = 0.;
  const double luminosity =
      bp->radiative_efficiency * bp->accretion_rate * c * c;

  /* Factor in the radiative efficiency, don't subtract
   * jet BZ efficiency (spin is fixed) */
  mass_rate = (1. - bp->radiative_efficiency) * bp->accretion_rate;

  /* This is used for jet feedback later */
  bp->radiative_luminosity = luminosity;

#ifdef OBSIDIAN_DEBUG_CHECKS
  if (isnan(dt)) error("dt nan");
#endif
  /* Integrate forward in time */
  double delta_mass = mass_rate * dt;

  /* If desired we put mass into accretion disk which feeds BH on some
   * frac of tdyn
   */
  if (tdyn_inv > 0.f) {
    /* Add accreted mass into a reservoir representing BH accretion disk */
    bp->accretion_disk_mass += delta_mass;

    /* Compute mass that will actually go into BH */
    delta_mass = bp->accretion_disk_mass * (1. - exp(-dt * tdyn_inv));

    /* This mass gets removed from the accretion disk */
    if (bp->accretion_disk_mass > delta_mass) {
      bp->accretion_disk_mass -= delta_mass;
    } else {
      delta_mass = bp->accretion_disk_mass;
      bp->accretion_disk_mass = 0.;
    }

    /* Recompute accretion rate based on the reservoir change */
    bp->accretion_rate = delta_mass / (dt * (1. - bp->radiative_efficiency));
  }

  bp->subgrid_mass += delta_mass;
  bp->total_accreted_mass += delta_mass;

  /* Note: bp->subgrid_mass has been integrated, so avoid BH_mass variable */
  if (bp->state == BH_states_adaf && BH_mass > my_adaf_mass_limit) {

    /* ergs to dump in a kernel-weighted fashion */
    if (props->adaf_wind_mass_loading == 0.f) {
      if (bp->subgrid_mass < my_adaf_mass_limit) {
        bp->adaf_energy_to_dump = 0.f;
      }
      /*else if (bp->subgrid_mass < 1.5f * my_adaf_mass_limit) {
        bp->adaf_energy_to_dump *=
            4.f * powf(bp->subgrid_mass / my_adaf_mass_limit - 1.f, 2.f);
      }*/
      else {
        bp->adaf_energy_to_dump =
            get_black_hole_coupling(bp, props, cosmo, phys_const) *
            props->adaf_disk_efficiency * bp->accretion_rate * c * c * dt;
      }
    } else {
      const float adaf_v2 = props->adaf_wind_speed * props->adaf_wind_speed;
      const float mass_this_step =
          props->adaf_wind_mass_loading * bp->accretion_rate * dt;
      bp->adaf_energy_to_dump += 0.5f * mass_this_step * adaf_v2;
    }
  }

  if (bp->state == BH_states_adaf ||
      (props->slim_disk_jet_active && bp->state == BH_states_slim_disk) ||
      (bp->radiative_luminosity > props->lum_thresh_always_jet &&
      props->lum_thresh_always_jet > 0.f)) {

    float jet_velocity = black_hole_compute_jet_velocity(bp, cosmo, props);

    /* If there is a variable jet velocity we must recalculate the mass loading
     */
    if (jet_velocity != props->jet_velocity) {
      const double c_over_v = phys_const->const_speed_light_c / jet_velocity;

      if (props->jet_loading_type == BH_jet_momentum_loaded) {
        bp->jet_mass_loading = props->jet_efficiency * c_over_v;
      } else if (props->jet_loading_type == BH_jet_mixed_loaded) {
        const double energy_loading =
            2. * props->jet_efficiency * pow(c_over_v, 2.);
        const double momentum_loading = props->jet_efficiency * c_over_v;

        /* Divide the contribution between energy and momentum loading */
        const double energy_term = props->jet_frac_energy * energy_loading;
        const double momentum_term =
            (1. - props->jet_frac_energy) * momentum_loading;

        bp->jet_mass_loading = energy_term + momentum_term;
      } else {
        bp->jet_mass_loading = 2. * props->jet_efficiency * pow(c_over_v, 2.);
      }

      /* Psi_jet*M_dot,acc*dt is the total mass expected in the jet this step */
      bp->jet_mass_reservoir += bp->jet_mass_loading * bp->accretion_rate * dt;
    } else {
      bp->jet_mass_reservoir +=
          props->jet_mass_loading * bp->accretion_rate * dt;
    }
  }

  if (bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative losses, so we must decrease the particle mass
     * in proportion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -= bp->radiative_efficiency * bp->accretion_rate * dt;

    if (bp->mass < 0) {
      error("Black hole %lld reached negative mass (%g). Trouble ahead...",
            bp->id, bp->mass);
    }

    /* Make sure not to destroy low mass galaxies */
    if (bp->subgrid_mass > props->minimum_black_hole_mass_unresolved &&
        bp->state != BH_states_adaf) {
      /* Make sure if many mergers have driven up the dynamical mass at low
       * subgrid mass, that we still kick out particles! */
      const float psi = (1.f - bp->f_accretion) / bp->f_accretion;
      bp->unresolved_mass_reservoir += psi * bp->accretion_rate * dt;
    }
  }

  /* Increase the subgrid angular momentum according to what we accreted
   * Note that this is already in physical units, a factors from velocity and
   * radius cancel each other. Also, the circular velocity contains an extra
   * smoothing length factor that we undo here. */
  const double m_times_r = (mass_rate * dt) * bp->h;
  /* L = m * r * v */
  bp->accreted_angular_momentum[0] += m_times_r * bp->circular_velocity_gas[0];
  bp->accreted_angular_momentum[1] += m_times_r * bp->circular_velocity_gas[1];
  bp->accreted_angular_momentum[2] += m_times_r * bp->circular_velocity_gas[2];

  /* Keep v_kick physical, there are a lot of comparisons */
  bp->v_kick = get_black_hole_wind_speed(props, phys_const, bp);

  /* This is always true in the ADAF mode; only heating happens */
  if (bp->state == BH_states_adaf) bp->v_kick = 0.f;

#ifdef OBSIDIAN_DEBUG_CHECKS
  const float galaxy_sfr = bp->galaxy_data.stellar_mass * bp->galaxy_data.specific_sfr;
  tdyn_inv = (tdyn_inv > 0.f) ? tdyn_inv : FLT_MIN;
  message(
      "BH_ACC: z=%g bid=%lld ms=%g dms=%g sfr=%g mbh=%g dmbh=%g state=%d "
      "torque=%g bondi=%g fEdd=%g facc=%g fsupp=%g mcold=%g mhot=%g mdisk=%g"
      " tin=%g vkick=%g dmass=%g radeff=%g mres=%g tdyn=%g",
      cosmo->z, bp->id, galaxy_mstar * props->mass_to_solar_mass,
      galaxy_sfr * dt * props->mass_to_solar_mass,
      galaxy_sfr * props->mass_to_solar_mass / props->time_to_yr,
      bp->subgrid_mass * props->mass_to_solar_mass,
      delta_mass * props->mass_to_solar_mass, bp->state,
      torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
      bondi_accr_rate * props->mass_to_solar_mass / props->time_to_yr,
      bp->eddington_fraction, bp->f_accretion,
      1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass /
               fabs(props->bh_characteristic_suppression_mass) * cosmo->a),
      bp->cold_gas_mass * props->mass_to_solar_mass,
      bp->hot_gas_mass * props->mass_to_solar_mass,
      corot_gas_mass * props->mass_to_solar_mass, props->time_to_Myr / tdyn_inv,
      bp->v_kick / props->kms_to_internal, delta_mass, bp->radiative_efficiency,
      bp->accretion_disk_mass, (1.f / tdyn_inv) * props->time_to_Myr);

  message(
      "BH_STATES: id=%lld, new_state=%d, predicted_mdot_medd=%g, "
      "eps_r=%g, f_Edd=%g, f_acc=%g, "
      "luminosity=%g, accr_rate=%g Msun/yr, coupling=%g, v_kick=%g km/s, "
      "jet_mass_reservoir=%g Msun unresolved_reservoir=%g Msun "
      "jet_mass_loading=%g",
      bp->id, bp->state, predicted_mdot_medd, bp->radiative_efficiency,
      bp->eddington_fraction, bp->f_accretion,
      bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs,
      bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr,
      get_black_hole_coupling(bp, props, cosmo, phys_const),
      bp->v_kick / props->kms_to_internal,
      bp->jet_mass_reservoir * props->mass_to_solar_mass,
      bp->unresolved_mass_reservoir * props->mass_to_solar_mass,
      bp->jet_mass_loading);
#endif

#define OBSIDIAN_BH_DETAILS
#ifdef OBSIDIAN_BH_DETAILS
  const float galaxy_sfr = bp->galaxy_data.stellar_mass * bp->galaxy_data.specific_sfr;
  printf(
      "BH_DETAILS "
      "z=%2.12f bid=%lld galM*=%g galSFR=%g"
      " Mdyn=%g MBH=%g h=%g Mres=%g BHAR=%g Bondi=%g torque=%g dt=%g dM=%g"
      " nH=%g Thot=%g SFR=%g mngb=%g "
      " mhot=%g mcold=%g m*=%g mdisk=%g mcorot=%g "
      " x=%2.10f y=%2.10f z=%2.10f "
      " vx=%2.7f vy=%2.7f vz=%2.7f "
      " Lgasx=%g Lgasy=%g Lgasz=%g  Lbhx=%g Lbhy=%g Lbhz=%g"
      " Lrad=%g state=%d facc=%g eff=%g"
      " fedd=%g madaf=%g mngb=%g tdyn=%g fsupp=%g\n",
      cosmo->z, bp->id, 
      galaxy_mstar * props->mass_to_solar_mass,
      galaxy_sfr * props->mass_to_solar_mass / props->time_to_yr,
      bp->mass * props->mass_to_solar_mass,
      bp->subgrid_mass * props->mass_to_solar_mass,
      bp->h * cosmo->a * props->length_to_parsec / 1.0e3f,
      bp->jet_mass_reservoir * props->mass_to_solar_mass,
      bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr,
      bp->bondi_accretion_rate * props->mass_to_solar_mass / props->time_to_yr,
      (bp->accretion_rate - bp->bondi_accretion_rate) * props->mass_to_solar_mass / props->time_to_yr,
      dt * props->time_to_Myr,
      bp->accretion_rate * dt * props->mass_to_solar_mass,
      (bp->rho_gas * cosmo->a3_inv) * props->rho_to_n_cgs,
      bp->hot_gas_internal_energy * cosmo->a_factor_internal_energy /
          (props->T_K_to_int * props->temp_to_u_factor),
      bp->gas_SFR * props->mass_to_solar_mass / props->time_to_yr,
      bp->ngb_mass * props->mass_to_solar_mass,
      bp->hot_gas_mass * props->mass_to_solar_mass,
      bp->cold_gas_mass * props->mass_to_solar_mass, 
      bp->stellar_mass * props->mass_to_solar_mass, 
      bp->cold_disk_mass * props->mass_to_solar_mass, 
      (bp->cold_gas_mass - 2. * (bp->cold_gas_mass - bp->cold_disk_mass)) * props->mass_to_solar_mass, 
      bp->x[0] * cosmo->a * props->length_to_parsec / 1.0e3f,
      bp->x[1] * cosmo->a * props->length_to_parsec / 1.0e3f,
      bp->x[2] * cosmo->a * props->length_to_parsec / 1.0e3f,
      bp->v[0] * cosmo->a_inv / props->kms_to_internal,
      bp->v[1] * cosmo->a_inv / props->kms_to_internal,
      bp->v[2] * cosmo->a_inv / props->kms_to_internal,
      bp->angular_momentum_gas[0], bp->angular_momentum_gas[1],
      bp->angular_momentum_gas[2],
      bp->swallowed_angular_momentum[0], bp->swallowed_angular_momentum[1],
      bp->swallowed_angular_momentum[2],
      bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs,
      bp->state, bp->f_accretion, bp->radiative_efficiency,
      bp->eddington_fraction, my_adaf_mass_limit * props->mass_to_solar_mass,
      bp->gravitational_ngb_mass * props->mass_to_solar_mass,
      1.f / tdyn_inv * 1.e-6 * props->time_to_yr, f_suppress);
#endif
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
black_holes_get_repositioning_speed(const struct bpart *restrict bp,
                                    const struct black_holes_props *props,
                                    const struct cosmology *cosmo) {

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
 * @param phys_const The physical phys_const (in internal units).
 * @param cosmo The cosmological model.
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the temp
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart *restrict bp, const struct black_holes_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
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
    struct bpart *restrict bp) {

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
black_holes_store_potential_in_bpart(struct bpart *bp, const struct gpart *gp) {

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
black_holes_store_potential_in_part(struct black_holes_part_data *p_data,
                                    const struct gpart *gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param phys_const The physical phys_const in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 * @param xp The #xpart that became a black hole.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart *bp, const struct black_holes_props *props,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct part *p, const struct xpart *xp,
    const integertime_t ti_current) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time and density */
  bp->formation_scale_factor = cosmo->a;
  bp->formation_gas_density = hydro_get_physical_density(p, cosmo);

  /* Copy FoF galaxy data from spawning particle */
  bp->galaxy_data.stellar_mass = p->galaxy_data.stellar_mass;
  bp->galaxy_data.gas_mass = p->galaxy_data.gas_mass;
  bp->galaxy_data.specific_sfr = p->galaxy_data.specific_sfr;

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

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
  bp->last_repos_vel = 0.f;

  /* Copy over the splitting struct */
  bp->split_data = xp->split_data;

  /* Initial metal masses */
  const float gas_mass = hydro_get_mass(p);
  struct chemistry_bpart_data *bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data *p_chem = &p->chemistry_data;
  chemistry_bpart_from_part(bp_chem, p_chem, gas_mass);

  /* No swallowed angular momentum */
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;

  /* Last time of mergers */
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;

  /* First initialisation */
  black_holes_init_bpart(bp);

  bp->state = BH_states_slim_disk;

  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

#endif /* SWIFT_OBSIDIAN_BLACK_HOLES_H */
