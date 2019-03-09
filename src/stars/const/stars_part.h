/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_STAR_PART_H
#define SWIFT_DEFAULT_STAR_PART_H

/* Some standard headers. */
#include <stdlib.h>

/* Read chemistry */
#include "chemistry_struct.h"

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

  /*! Initial star mass */
  float mass_init;

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

  struct {
    /* Change in smoothing length over time. */
    float h_dt;

  } feedback;

  struct {
    /* Mass of ejecta */
    float mass;

    /* Mass fractions of ejecta */
    struct chemistry_part_data chemistry_data;

    float ejecta_specific_thermal_energy;

    /* Number of type 1a SNe per unit mass */
    float num_SNIa;

  } to_distribute;

  /* kernel normalisation factor (equivalent to metalweight_norm in
   * eagle_enrich.c:811, TODO: IMPROVE COMMENT) */
  float omega_normalisation_inv;

  /* total mass of neighbouring gas particles */
  float ngb_mass;

  /*! Tracer structure */
  struct tracers_xpart_data tracers_data;

  /*! Chemistry structure */
  struct chemistry_part_data chemistry_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_STARS
  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_STARS];

  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;
#endif

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Contains all the constants and parameters of the stars scheme
 */
struct stars_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal smoothing length */
  float h_max;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  /* Flag to switch between continuous and stochastic heating */
  int continuous_heating;

  /* Fraction of energy in SNIa (Note: always set to 1 in EAGLE, so may be not
   * necessary) */
  float SNIa_energy_fraction;

  /* Desired temperature increase due to supernovae */
  float SNe_deltaT_desired;

  /* Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /* Energy released by one supernova */
  float total_energy_SNe;

  /* Temperature due to SNe (corresponding to units_factor1 in EAGLE) */
  float SNe_temperature;

  /* Timescale for feedback (used only for testing in const feedback model) */
  float feedback_timescale;

  /* Number of supernovae per solar mass (used only for testing in const
   * feedback model) */
  float sn_per_msun;

  /* Solar mass (used only for testing in const feedback model) */
  float const_solar_mass;

  /* Flag for testing energy injection */
  int const_feedback_energy_testing;

  // CHANGE THIS TO BE CONSISTENT WITH RAND MAX USED IN STAR FORMATION
  double inv_rand_max;
};

#endif /* SWIFT_DEFAULT_STAR_PART_H */
