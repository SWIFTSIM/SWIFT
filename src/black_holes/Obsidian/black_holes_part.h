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
#ifndef SWIFT_OBSIDIAN_BLACK_HOLE_PART_H
#define SWIFT_OBSIDIAN_BLACK_HOLE_PART_H

#include "black_holes_struct.h"
#include "chemistry_struct.h"
#ifdef WITH_FOF_GALAXIES
#include "fof_struct.h"
#endif
#include "particle_splitting_struct.h"
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
  struct gpart *gpart;

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

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char depth_h;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

#ifdef WITH_FOF_GALAXIES
  /*! Struct for host galaxy information */
  struct fof_galaxy_data galaxy_data;

  /*! Distance to host galaxy center of mass */
  float galactocentric_radius;
#endif

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

  /*! Total accreted mass of the black hole (not including mass merged in
   * from other black holes) */
  float total_accreted_mass;

  /*! Instantaneous accretion rate */
  float accretion_rate;

  /*! Instantaneous Bondi component of accretion rate */
  float bondi_accretion_rate;

  /*! Density of the gas surrounding the black hole. */
  float rho_gas;

  /*! Internal energy of the gas surrounding the black hole. */
  float internal_energy_gas;

  /*! The mass of hot gas surrounding the black hole */
  float hot_gas_mass;

  /*! The mass of cold gas surrounding the black hole */
  float cold_gas_mass;

  /*! The total SFR of cold gas surrounding the black hole */
  float gas_SFR;

  /*! The current state of the black hole */
  int state;

  /*! The radiative efficiency associated with M_dot,bh */
  float radiative_efficiency;

  /*! The large scale accretion rate onto the black hole */
  float m_dot_inflow;

  /*! The amount of jet energy available */
  float jet_mass_reservoir;

  /*! The amount of jet mass kicked this time step */
  float jet_mass_kicked_this_step;

  /*! The mass loading in the jet for the variable velocity scheme */
  float jet_mass_loading;

  /*! The amount of unresolved mass available to kick */
  float unresolved_mass_reservoir;

  /*! The amount of unresolved mass kicked this step */
  float unresolved_mass_kicked_this_step;

  /*! Energy to dump this step via the ADAF hot-wind, kernel-weighted */
  float adaf_energy_to_dump;

  /*! Energy injected this step in the ADAF mode */
  float adaf_energy_used_this_step;

  /*! sum(mi * wi) weights for accretion/feedback */
  float kernel_wt_sum;

  /*! sum(mi * wi) weights for adaf heating */
  float adaf_wt_sum;

  /*! The mass of cold disk around the black hole */
  float cold_disk_mass;

  /*! Mass in accretion disk from which BH accretes */
  float accretion_disk_mass;

  /*! The mass-weighted internal energy surrounding the black hole (unsmoothed)
   */
  float hot_gas_internal_energy;

  /*! Smoothed sound speed of the gas surrounding the black hole. */
  float sound_speed_gas;

  /*! Total gravitational gas mass within the kernel */
  float gravitational_ngb_mass;

  /*! Subgrid physical sound speed of the gas (updated when using the subgrid
   * Bondi model) */
  float sound_speed_subgrid_gas;

  /*! Smoothed velocity of the gas surrounding the black hole,
   * in the frame of the black hole (internal units) */
  float velocity_gas[3];

  /*! The real angular momentum of the gas in the kernel */
  float angular_momentum_gas[3];

  /*! Circular velocity of the gas around the black hole at the smoothing
   * radius (calculated as j_gas / h_BH, where j is specific ang. mom.) */
  float circular_velocity_gas[3];

  /*! Total mass of the gas neighbours. */
  float ngb_mass;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Integer number of gravitational neighbors */
  int num_gravitational_ngbs;

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

  /*! Total (physical) angular momentum accumulated from subgrid accretion */
  float accreted_angular_momentum[3];

  /*! Eddington fractions */
  float eddington_fraction;

  /*! Wind velocity kick */
  float v_kick;

  /*! Fraction of Mdot,inflow that should be accreted, the rest is a wind */
  float f_accretion;

  /*! Bulge mass of stars within the kernel (twice the counter-rotating mass) */
  float stellar_bulge_mass;

  /*! The mass of stars within the kernel */
  float stellar_mass;

  /*! The radiative luminosity of the black hole */
  float radiative_luminosity;

  /*! How much energy has been given away in this timestep? */
  float delta_energy_this_timestep;

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

    /*! Last scale factor the BH had a major merger */
    float last_major_merger_scale_factor;
  };

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

#endif /* SWIFT_OBSIDIAN_BLACK_HOLE_PART_H */
