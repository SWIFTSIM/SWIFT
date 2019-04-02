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

    /* Total metal mass released per unit initial stellar mass */
    float total_metal_mass;

    /* Total mass released by element per unit initial stellar mass */
    float metal_mass[chemistry_element_count];
    
    /*! Mass coming from SNIa */
    float mass_from_SNIa;

    /*! Fraction of total gas mass in metals coming from SNIa */
    float metal_mass_from_SNIa;

    /*! Mass coming from AGB */
    float mass_from_AGB;

    /*! Fraction of total gas mass in metals coming from AGB */
    float metal_mass_from_AGB;

    /*! Mass coming from SNII */
    float mass_from_SNII;

    /*! Fraction of total gas mass in metals coming from SNII */
    float metal_mass_from_SNII;

    /*! Fraction of total gas mass in Iron coming from SNIa */
    float Fe_mass_from_SNIa;

    /* Number of type Ia SNe per unit mass */
    float num_SNIa;

    /* Number of type II SNe per unit mass */
    float num_SNII;

    /* Number of SNe in timestep  */
    float num_SNe;

  } to_distribute;

  /* Normalisation factor for density weight fraction for feedback (equivalent to metalweight_norm in EAGLE, see eagle_enrich.c:811) */
  float density_weight_frac_normalisation_inv;

  /* total mass of neighbouring gas particles */
  float ngb_mass;

  /*! Union for the birth time and birht scale factor */
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
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {
  // insert comments to differentiate between sph and non-sph fields
  double *mass;
  double *metallicity;
  double *yield_IMF_resampled;
  double *yield;
  double *ejecta_IMF_resampled;
  double *ejecta;
  double *total_metals_IMF_resampled;
  double *total_metals;
};

/**
 * @brief Stores tables to determine stellar lifetimes
 */
struct lifetime_table {
  /* number of elements, mass, and initial metallicity bins */
  int n_mass;
  int n_z;

  /* table of masses */
  double *mass;

  /* table of metallicities */
  double *metallicity;

  /* table of lifetimes depending on mass an metallicity */
  double **dyingtime;
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

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;

  struct {
    /* Flag to switch between continuous and stochastic heating */
    int continuous_heating;

    /* Desired temperature increase due to supernovae */
    float SNe_deltaT_desired;

    /* Conversion factor from temperature to internal energy */
    float temp_to_u_factor;

    /* Energy released by one supernova */
    float total_energy_SNe;

    /* Kinetic energy of SN ejecta per unit mass (check name with Richard)*/
    float ejecta_specific_thermal_energy;

    /* Solar mass */
    float const_solar_mass;

    /* Flag for testing energy injection */
    int const_feedback_energy_testing;

    /* Yield tables for AGB and SNII  */
    struct yield_table yield_AGB;
    struct yield_table yield_SNII;

    /* Array of adjustment factors for SNII  */
    double *typeII_factor;

    /* Yield tables for SNIa  */
    double *yield_SNIa_IMF_resampled;
    double yield_SNIa_total_metals_IMF_resampled;
    double *yields_SNIa;

    /* Parameters to SNIa enrichment model  */
    int SNIa_mode;
    float SNIa_efficiency;
    float SNIa_timescale;

    /* Mass transfer due to enrichment  */
    int SNIa_mass_transfer;
    int SNII_mass_transfer;
    int AGB_mass_transfer;

    /* Arrays for elements being tracked */
    char **SNIa_element_names;
    char **SNII_element_names;
    char **AGB_element_names;

    /* Element name string length */
    int element_name_length;

    /* Dimensions of arrays in yield tables */
    int SNIa_n_elements;
    int SNII_n_mass;
    int SNII_n_elements;
    int SNII_n_z;
    int AGB_n_mass;
    int AGB_n_elements;
    int AGB_n_z;

    /* log10 of max and min allowable masses for SNII and SNIa in msun */
    float log10_SNII_min_mass_msun;
    float log10_SNII_max_mass_msun;
    float log10_SNIa_max_mass_msun;

    /* Array of mass bins for yield calculations */
    double *yield_mass_bins;

    /* Parameters for IMF  */
    char IMF_Model[10];
    float IMF_Exponent;
    float *imf;
    float *imf_mass_bin;
    float *imf_mass_bin_log10;
    int n_imf_mass_bins;
    float imf_max_mass_msun;
    float imf_min_mass_msun;
    float log10_imf_min_mass_msun;
    float log10_imf_max_mass_msun;

    /* Table of lifetime values */
    struct lifetime_table lifetimes;

    /* Location of yield tables */
    char yield_table_path[50];

    /* number of type II supernovae per solar mass */
    float num_SNII_per_msun;

    /* wind delay time for SNII */
    float SNII_wind_delay;
  } feedback;

};

#endif /* SWIFT_EAGLE_STAR_PART_H */
