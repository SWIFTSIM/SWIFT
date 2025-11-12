/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_GEAR_SINK_PART_H
#define SWIFT_GEAR_SINK_PART_H

#include "timeline.h"

#define sink_need_unique_id 1

/**
 * @brief Particle fields for the sink particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct sink {

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

  /* Particle smoothing length, or r_cut/kernel_gamma if using a fixed cutoff*/
  float h;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Sink particle mass */
  float mass;

  /*! Sink target mass. In Msun. */
  float target_mass_Msun;

  /* Mass of the IMF this sinks is currently affected to. In internal units. */
  double mass_IMF;

  /*! Integer number of neighbours */
  int num_ngbs;

  /*! Mass of the sink before starting the star spawning loop */
  float mass_tot_before_star_spawning;

  /*! Sink target stellar type */
  enum stellar_type target_type;

  /*! Union for the birth time and birth scale factor */
  union {

    /*! Birth time */
    float birth_time;

    /*! Birth scale factor */
    float birth_scale_factor;
  };

  struct {

    /*! Minimal gas smoothing length */
    float minimal_h_gas;

    /*! Density of the gas surrounding the sink. */
    float rho_gas;

    /*! Smoothed sound speed of the gas surrounding the sink. */
    float sound_speed_gas;

    /*! Smoothed velocity of the gas surrounding the sink, in the frame of the
      sink (internal units) */
    float velocity_gas[3];

    /*! Minimal t_c between all sink neighbours */
    float minimal_sink_t_c;

    /*! Minimal dynamical time between all sink neighbours */
    float minimal_sink_t_dyn;

    /*! Total mass that passes all criteria before the accretion limit */
    float mass_eligible_swallow;

    /*! Swallowed mass during this timestep */
    float mass_swallowed;
  } to_collect;

  /*! Particle time bin */
  timebin_t time_bin;

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char depth_h;

  /*! Number of stars spawned by this sink */
  int n_stars;

  /*! Total (physical) angular momentum accumulated by swallowing particles */
  float swallowed_angular_momentum[3];

  /*! Total number of sink merger events (including sink swallowed
   * by merged-in sinks) */
  int number_of_sink_swallows;

  /*! Total number of sink merger events (excluding sink swallowed
   * by merged-in sinks) */
  int number_of_direct_sink_swallows;

  /*! Total number of gas particles swallowed (including particles swallowed
   * by merged-in sinks) */
  int number_of_gas_swallows;

  /*! Total number of gas particles swallowed (excluding particles swallowed
   * by merged-in sinks) */
  int number_of_direct_gas_swallows;

  /*! Flag to determine if a sink has already changed its IMF from pop III to
     pop II. */
  int has_IMF_changed_from_popIII_to_popII;

  /*! Chemistry information (e.g. metal content at birth, swallowed metal
   * content, etc.) */
  struct chemistry_sink_data chemistry_data;

  /*! sink merger information (e.g. merging ID) */
  struct sink_sink_data merger_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_SINKS
  /*! Number of interactions in merger SELF and PAIR */
  int num_ngb_merger;

  /*! List of interacting particles in merger SELF and PAIR */
  long long ids_ngbs_merger[MAX_NUM_OF_NEIGHBOURS_SINKS];

  /*! Number of interactions in compute formation SELF and PAIR */
  int num_ngb_formation;

  /*! List of interacting particles in compute formation SELF and PAIR */
  long long ids_ngbs_formation[MAX_NUM_OF_NEIGHBOURS_SINKS];

  /*! Number of interactions in compute formation SELF and PAIR */
  int num_ngb_accretion;

  /*! List of interacting particles in compute formation SELF and PAIR */
  long long ids_ngbs_accretion[MAX_NUM_OF_NEIGHBOURS_SINKS];
#endif

#ifdef SWIFT_SINK_DENSITY_CHECKS

  /* Integer number of neighbours in the density loop */
  int N_check_density;

  /* Exact integer number of neighbours in the density loop */
  int N_check_density_exact;

  /*! Has this particle interacted with any unhibited neighbour? */
  char inhibited_check_exact;

  float n_check;

  float n_check_exact;

  float rho_check;

  /*! Exact value of the density field obtained via brute-force loop */
  float rho_check_exact;

#endif
} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_GEAR_SINK_PART_H */
