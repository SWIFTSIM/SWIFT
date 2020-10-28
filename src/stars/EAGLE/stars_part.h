/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STAR_PART_H
#define SWIFT_EAGLE_STAR_PART_H

/* Some standard headers. */
#include <stdlib.h>

/* Read additional aubgrid models */
#include "chemistry_struct.h"
#include "feedback_struct.h"
#include "particle_splitting_struct.h"
#include "star_formation_struct.h"
#include "tracers_struct.h"

/**
 * @brief Particle fields for the star particles.
 *
 * All quantities related to gravity are stored in the associate #gpart.
 */
struct spart {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between current position and position at last tree rebuild. */
  float x_diff_sort[3];

  /*! Particle velocity. */
  float v[3];

  /*! Star mass */
  float mass;

  /*! Particle smoothing length. */
  float h;

  float rho;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  /*! Union for the birth time and birth scale factor */
  union {

    /*! Birth time */
    float birth_time;

    /*! Birth scale factor */
    float birth_scale_factor;
  };

  /*! Scale-factor / time at which this particle last did enrichment */
  float last_enrichment_time;

  /*! Initial star mass */
  float mass_init;

  /*! Total number of SNII injection events this star performed */
  int number_of_SNII_events;

  /*! Feedback energy fraction */
  float f_E;

  /*! The physical birth density */
  float birth_density;

  /*! The birth temperature */
  float birth_temperature;

  /*! Total number of (expected) feedback heating events so far */
  float number_of_heating_events;

  /*! Star formation struct */
  struct star_formation_spart_data sf_data;

  /*! Feedback structure */
  struct feedback_spart_data feedback_data;

  /*! Tracer structure */
  struct tracers_xpart_data tracers_data;

  /*! Chemistry structure */
  struct chemistry_spart_data chemistry_data;

  /*! Splitting structure */
  struct particle_splitting_data split_data;

  /*! Particle time bin */
  timebin_t time_bin;

  /*! Number of time-steps since the last enrichment step */
  char count_since_last_enrichment;

#ifdef WITH_LOGGER
  /* Additional data for the particle logger */
  struct logger_part_data logger_data;
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  /* Integer number of neighbours in the density loop */
  int N_density;

  /* Exact integer number of neighbours in the density loop */
  int N_density_exact;

  /*! Has this particle interacted with any unhibited neighbour? */
  char inhibited_exact;

  /*! Exact value of the density field obtained via brute-force loop */
  float rho_exact;
#endif

#ifdef DEBUG_INTERACTIONS_STARS

  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_STARS];

  /*! Number of interactions in the feedback SELF and PAIR */
  int num_ngb_feedback;

  /*! List of interacting particles in the feedback SELF and PAIR */
  long long ids_ngbs_feedback[MAX_NUM_OF_NEIGHBOURS_STARS];
#endif

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Contains all the constants and parameters of the stars scheme
 */
struct stars_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weighted number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /*! Are we overwriting the stars' birth time read from the ICs? */
  int overwrite_birth_time;

  /*! Are we overwriting the stars' birth density read from the ICs? */
  int overwrite_birth_density;

  /*! Are we overwriting the stars' birth temperature read from the ICs? */
  int overwrite_birth_temperature;

  /*! Value to set birth time of stars read from ICs */
  float spart_first_init_birth_time;

  /*! Value to set birth density of stars read from ICs */
  float spart_first_init_birth_density;

  /*! Value to set birth temperature of stars read from ICs */
  float spart_first_init_birth_temperature;

  /*! Maximal time-step length of young stars (internal units) */
  double max_time_step_young;

  /*! Maximal time-step length of old stars (internal units) */
  double max_time_step_old;

  /*! Age threshold for the young/old transition (internal units) */
  double age_threshold;

  /*! Age threshold for the transition to unlimited time-step size (internal
   * units) */
  double age_threshold_unlimited;
};

#endif /* SWIFT_EAGLE_STAR_PART_H */
