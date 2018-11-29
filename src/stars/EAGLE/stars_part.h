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
#ifndef SWIFT_EAGLE_STAR_PART_H
#define SWIFT_EAGLE_STAR_PART_H

/* Some standard headers. */
#include <stdlib.h>

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

  /*! Particle velocity. */
  float v[3];

  /*! Star mass */
  float mass;

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

  // should these go somewhere else? perhaps their own struct?
  float mass_from_agb;
  float metals_from_agb;
  float mass_from_snii;
  float metals_from_snii;
  float mass_from_snia;
  float metals_from_snia;
  float iron_from_snia;
  float metal_mass_released;
  float metals_released[chemistry_element_count];

  float num_snia;

  float age_Gyr;

} SWIFT_STRUCT_ALIGN;

struct yield_table {
  // Names coincide with EAGLE, change to be more descriptive?
  float *mass;
  float *metallicity;
  float *SPH;
  float *yield;
  float *ejecta_SPH;
  float *ejecta;
  float *total_metals_SPH;
  float *total_metals;
  int N_Z;
};

struct lifetime_table {
  /* number of elements, mass, and initial metallicity bins */
  int n_mass, n_z; 

  /* mass bins[N_MASS] (check units, give more descriptive name) */
  float *mass;    

  /* metallicity bins[N_Z] (total metallicity/total mass) */
  float *metallicity; 

  /* size [N_MASS*N_Z] with index_2d(imass,iz) */
  float **dyingtime; 
};


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

  // restructure tables to be neater.
  struct yield_table yield_AGB;
  struct yield_table yield_SNII;

  float *yield_SNIa_SPH;
  float yield_SNIa_total_metals_SPH;
  float *yields_SNIa;
  
  char **SNIa_element_names;
  char **SNII_element_names;
  char **AGB_element_names;

  int SNIa_n_elements;
  int SNII_n_mass;
  int SNII_n_elements;
  int SNII_n_z;
  int AGB_n_mass;
  int AGB_n_elements;
  int AGB_n_z;

  char IMF_Model[10];
  float IMF_Exponent;
  float *imf_by_number;
  float *imf_by_number1;
  float *imf_mass_bin;
  float *imf_mass_bin_log10;

  int stellar_lifetime_flag; // 0 for padovani & matteucci 1993, 1 for maeder & meynet 1989, 2 for portinari et al. 1998. 

  // Table of lifetime values, should these be somewhere else?
  struct lifetime_table lifetimes;

  int SNIa_mode;
  float SNIa_efficiency;
  float SNIa_timescale;

  int SNIa_mass_transfer;
  int SNII_mass_transfer;
  int AGB_mass_transfer;

  char *yield_table_path;
};

#endif /* SWIFT_EAGLE_STAR_PART_H */
