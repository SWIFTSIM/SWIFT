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

/* Read chemistry */
#include "chemistry_struct.h"
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
  struct gpart *gpart;

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

  /*! Particle smoothing length. */
  float h;

  /*! Density of the gas surrounding the star. */
  float rho_gas;

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

    /* Total metal mass released */
    float total_metal_mass;

    /* Total mass released by element */
    float metal_mass[chemistry_element_count];

    /*! Total mass released due to SNIa */
    float mass_from_SNIa;

    /*! Total metal mass released due to SNIa */
    float metal_mass_from_SNIa;

    /*! Total iron mass released due to SNIa */
    float Fe_mass_from_SNIa;

    /*! Total mass released due to SNII */
    float mass_from_SNII;

    /*! Total metal mass released due to SNII */
    float metal_mass_from_SNII;

    /*! Total mass released due to AGB */
    float mass_from_AGB;

    /*! Total metal mass released due to AGB */
    float metal_mass_from_AGB;

    /* Number of type Ia SNe per unit mass */
    float num_SNIa;

    /* Number of type II SNe per unit mass */
    float num_SNII;

    /* Number of SNe in timestep  */
    float num_SNe;

    /* Energy change due to thermal and kinetic energy of ejecta */
    float d_energy;

    /* Probability for heating neighbouring gas particles */
    float heating_probability;

  } to_distribute;

  /* Normalisation factor for density weight fraction for feedback (equivalent
   * to metalweight_norm in EAGLE, see eagle_enrich.c:811) */
  float density_weighted_frac_normalisation_inv;

  /* total mass (unweighted) of neighbouring gas particles */
  float ngb_mass;

  /*! Union for the birth time and birth scale factor */
  union {

    /*! Birth time */
    float birth_time;

    /*! Birth scale factor */
    float birth_scale_factor;
  };

  /*! Birth density */
  float birth_density;

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
  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS_STARS];

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS_STARS];
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
};

#endif /* SWIFT_EAGLE_STAR_PART_H */
