/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIMBA_BLACK_HOLE_PART_H
#define SWIFT_SIMBA_BLACK_HOLE_PART_H

/*! The total number of rays used in AGN feedback */
#define eagle_blackhole_number_of_rays FEEDBACK_NR_RAYS_AGN

#include "black_holes_struct.h"
#include "chemistry_struct.h"
#ifdef WITH_FOF_GALAXIES
#include "fof_struct.h"
#endif
#include "particle_splitting_struct.h"
#include "rays_struct.h"
#include "timeline.h"

/**
 * @brief Particle fields for the black hole particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct bpart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Particle velocity. */
  float v[3];

  /*! Black hole mass */
  float mass;

  /*! Black hole mass at the start of each step, prior to any nibbling */
  float mass_at_start_of_step;

  /* Particle cutoff radius. */
  float h;

  /*! Particle time bin */
  timebin_t time_bin;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Union for the formation time and formation scale factor */
  union {

    /*! Formation time */
    float formation_time;

    /*! Formation scale factor */
    float formation_scale_factor;
  };

  /*! Physical density of the converted part (internal units) */
  float formation_gas_density;

  /*! Subgrid mass of the black hole */
  float subgrid_mass;

  /*! The mass marked for accretion this timestep */
  float mass_accreted_this_step;

  /*! Total accreted mass of the black hole (not including mass merged in
   * from other black holes) */
  float total_accreted_mass;

  /*! Energy reservoir for feedback */
  float energy_reservoir;

  /*! Instantaneous accretion rate */
  float accretion_rate;

  /*! Density of the gas surrounding the black hole. */
  float rho_gas;

  /*! Internal energy of the gas surrounding the black hole. */
  float internal_energy_gas;

  /*! The mass of hot gas surrounding the black hole */
  float hot_gas_mass;

  /*! The mass of cold gas surrounding the black hole */
  float cold_gas_mass;

  /*! The mass-weighted internal energy surrounding the black hole (unsmoothed)
   */
  float hot_gas_internal_energy;

  /*! Smoothed sound speed of the gas surrounding the black hole. */
  float sound_speed_gas;

  /*! Subgrid physical density of the gas (updated when using the subgrid Bondi
   * model) */
  float rho_subgrid_gas;

  /*! Subgrid physical sound speed of the gas (updated when using the subgrid
   * Bondi model) */
  float sound_speed_subgrid_gas;

  /*! Circular velocity of the gas around the black hole at the smoothing
   * radius (calculated as j_gas / h_BH, where j is specific ang. mom.) */
  float circular_velocity_gas[3];

  /*! Total mass of the gas neighbours. */
  float ngb_mass;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Number of seeds in this BH (i.e. itself + the merged ones) */
  int cumulative_number_seeds;

  /*! Total number of BH merger events (i.e. not including all progenies) */
  int number_of_mergers;

  /*! Total number of gas particles swallowed (including particles swallowed
   * by merged-in black holes) */
  int number_of_gas_swallows;

  /*! Total number of gas particles swallowed (excluding particles swallowed
   * by merged-in black holes) */
  int number_of_direct_gas_swallows;

  /*! Total number of times the black hole has been repositioned (excluding
   * repositionings of merged-in black holes) */
  int number_of_repositions;

  /*! Total number of times a black hole attempted repositioning (including
   * cases where it was aborted because the black hole was already at a
   * lower potential than all eligible neighbours) */
  int number_of_reposition_attempts;

  /* Velocity of most recent reposition jump */
  float last_repos_vel;

  /*! Total number of time steps in which the black hole was active. */
  int number_of_time_steps;

  /*! Total (physical) angular momentum accumulated by swallowing particles */
  float swallowed_angular_momentum[3];

  /*! Accretion boost factor */
  float accretion_boost_factor;

  /*! Instantaneous temperature increase for feedback */
  float AGN_delta_T;

  /*! Instantaneous energy reservoir threshold (num-to-heat) */
  float num_ngbs_to_heat;

  /*! Eddington fractions */
  float eddington_fraction;

  /*! Wind velocity kick */
  float v_kick;

  /*! Fraction of Mdot,inflow that should be accreted, the rest is a wind */
  float f_accretion;

  /*! Angular momentum of the gas */
  float angular_momentum_gas[3];

  /*! Specific angular momentum of the stars within the kernel */
  float specific_angular_momentum_stars[3];

  /*! Co-rotating cold gas mass within the kernel (estimate) */
  float cold_disk_mass;

  /*! Mass in accretion disk from which BH accretes */
  float accr_disk_mass;

  /*! Bulge mass of stars within the kernel (twice the counter-rotating mass) */
  float stellar_bulge_mass;

  /*! The mass of stars within the kernel */
  float stellar_mass;

  /*! The radiative luminosity of the black hole */
  float radiative_luminosity;

  /*! Integer (cumulative) number of energy injections in AGN feedback. At a
   * given time-step, an AGN-active BH may produce multiple energy injections.
   * The number of energy injections is equal to or more than the number of
   * particles heated by the BH during this time-step. */
  int AGN_number_of_energy_injections;

  /*! Integer (cumulative) number of AGN events. If a BH does feedback at a
   * given time-step, the number of its AGN events is incremented by 1. Each
   * AGN event may have multiple energy injections. */
  int AGN_number_of_AGN_events;

  /* Total energy injected into the gas in AGN feedback by this BH */
  float AGN_cumulative_energy;

  /*! BH feedback-limited time-step */
  float dt_heat;

  /*! BH accretion-limited time-step */
  float dt_accr;

  /*! How much energy has been given away in this timestep? */
  float delta_energy_this_timestep;

  /*! Union for the last AGN event time and the last AGN event scale factor */
  union {

    /*! Last AGN event time */
    float last_AGN_event_time;

    /*! Last AGN event scale-factor */
    float last_AGN_event_scale_factor;
  };

  /*! Union for the last high Eddington ratio point in time */
  union {

    /*! Last time the BH had a a high Eddington fraction */
    float last_high_Eddington_fraction_time;

    /*! Last scale factor the BH had a a high Eddington fraction */
    float last_high_Eddington_fraction_scale_factor;
  };

  /*! Union for the last minor merger point in time */
  union {

    /*! Last time the BH had a a high Eddington fraction */
    float last_minor_merger_time;

    /*! Last scale factor the BH had a a high Eddington fraction */
    float last_minor_merger_scale_factor;
  };

  /*! Union for the last major merger point in time */
  union {

    /*! Last time the BH had a a high Eddington fraction */
    float last_major_merger_time;

    /*! Last scale factor the BH had a a high Eddington fraction */
    float last_major_merger_scale_factor;
  };

  /*! Properties used in the feedback loop to distribute to gas neighbours. */
  struct {

    /*! Number of energy injections per time-step */
    int AGN_number_of_energy_injections;

    /*! Change in energy from SNII feedback energy injection */
    float AGN_delta_u;

  } to_distribute;

  struct {

    /*! Gravitational potential copied from the #gpart. */
    float potential;

    /*! Value of the minimum potential across all neighbours. */
    float min_potential;

    /*! Delta position to apply after the reposition procedure */
    double delta_x[3];

  } reposition;

  /*! Splitting structure */
  struct particle_splitting_data split_data;

  /*! Chemistry information (e.g. metal content at birth, swallowed metal
   * content, etc.) */
  struct chemistry_bpart_data chemistry_data;

  /*! Black holes merger information (e.g. merging ID) */
  struct black_holes_bpart_data merger_data;

#ifdef WITH_FOF_GALAXIES
  /*! Additional data used by the FoF */
  struct group_data group_data;
#endif

  /*! Isotropic AGN feedback information */
  struct ray_data rays[eagle_blackhole_number_of_rays];

  /*! Tracer structure */
  struct tracers_bpart_data tracers_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES];
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_SIMBA_BLACK_HOLE_PART_H */
