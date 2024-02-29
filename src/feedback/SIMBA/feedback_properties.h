/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIMBA_FEEDBACK_PROPERTIES_H
#define SWIFT_SIMBA_FEEDBACK_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "chemistry.h"
#include "hydro_properties.h"


/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {

  /*! Yield table mass bins */
  double *mass;

  /*! Yield table metallicity (metal mass fractions) bins */
  double *metallicity;

  /*! Array to store yield table (individual metals produced by the star)
     resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /*! Array to store yield table being read in */
  double *yield;

  /*! Array to store table of ejecta (metals alredy in the stars that are
    ejected) resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /*! Array to store table of ejecta being read in */
  double *ejecta;

  /*! Array to store table of total mass released ( metals produced by the star)
    resampled by IMF mass bins */
  double *total_metals_IMF_resampled;

  /* Array to store table of total mass released being read in */
  double *total_metals;
};

/**
 * @brief Stores tables to determine stellar lifetimes. Used for calculation of
 * IMF
 */
struct lifetime_table {

  /* table of masses */
  double *mass;

  /* table of metallicities */
  double *metallicity;

  /* table of lifetimes depending on mass an metallicity */
  double **dyingtime;
};

/**
 * @brief Functional form of the SNIa delay time distribution.
 */
enum eagle_feedback_SNIa_DTD {

  /*! Power-law with slope -1 */
  eagle_feedback_SNIa_DTD_power_law = 1,

  /*! Exponential model (EAGLE default) */
  eagle_feedback_SNIa_DTD_exponential = 2
};

/**
 * @brief Properties of the SIMBA feedback model.
 */
struct feedback_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we doing AGB enrichment? */
  int with_AGB_enrichment;

  /*! Are we doing SNII enrichment? */
  int with_SNII_enrichment;

  /*! Are we doing SNIa enrichment? */
  int with_SNIa_enrichment;

  /*! Are we doing SNII feedback? */
  int with_SNII_feedback;

  /*! Are we doing SNIa feedback? */
  int with_SNIa_feedback;

  /* ------------ Yield tables    ----------------- */

  /* Yield tables for AGB and SNII  */
  struct yield_table yield_AGB;
  struct yield_table yield_SNII;

  /* Arrays of yield tables for SNIa */
  double *yield_SNIa_IMF_resampled;
  double yield_SNIa_total_metals_IMF_resampled;
  double *yields_SNIa;

  /* Arrays for names of elements being tracked for each enrichment channel */
  char **SNIa_element_names;
  char **SNII_element_names;
  char **AGB_element_names;

  /* Array of mass bins for yield calculations */
  double *yield_mass_bins;

  /* Location of yield tables */
  char yield_table_path[200];

  /* ------------- Lifetime tracks   --------------- */

  /* Table of lifetime values */
  struct lifetime_table lifetimes;

  /* ------------- SNII parameters    --------------- */

  /* Array of adjustment factors for SNII  */
  float SNII_yield_factor[chemistry_element_count];

  /* ------------- SNIa parameters    --------------- */

  /* What delay time distribution are we using? */
  enum eagle_feedback_SNIa_DTD SNIa_DTD;

  /*! Normalisation of the SNIa DTD in the exponential model */
  float SNIa_DTD_exp_norm;

  /*! Time-scale of the SNIa decay function in the exponential model in
   * Giga-years */
  float SNIa_DTD_exp_timescale_Gyr;

  /*! Inverse of time-scale of the SNIa decay function in the exponential model
   * in Giga-years */
  float SNIa_DTD_exp_timescale_Gyr_inv;

  /*! Normalisation of the SNIa DTD in the power-law model */
  float SNIa_DTD_power_law_norm;

  /*! Stellar age below which no SNIa explode in Giga-years */
  float SNIa_DTD_delay_Gyr;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNIa_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNIa;

  /* ------------- AGB parameters    ---------------- */

  /*! Specific kinetic energy injected from AGB ejectas (in internal units). */
  float AGB_ejecta_specific_kinetic_energy;

  /* ------------- Conversion factors --------------- */

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! The mass of the sun in g */
  double solar_mass_in_g;

  /*! Conversion factor from internal mass unit to solar mass */
  double solar_mass_to_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /*! Conversion factor from km/s to cm/s */
  float kms_to_cms;

  /* ------------- Parameters for IMF --------------- */

  /*! Array to store calculated IMF */
  double *imf;

  /*! Arrays to store IMF mass bins */
  double *imf_mass_bin;

  /*! Arrays to store IMF mass bins (log10)*/
  double *imf_mass_bin_log10;

  /*! Minimal stellar mass considered by the IMF (in solar masses) */
  double imf_min_mass_msun;

  /*! Maximal stellar mass considered by the IMF (in solar masses) */
  double imf_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered by the IMF (in solar masses)
   */
  double log10_imf_max_mass_msun;

  /* ------------ SNe feedback properties ------------ */

  /*! Minimal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_min_mass_msun;

  /*! Maximal stellar mass considered for SNII feedback (in solar masses) */
  double SNII_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered for SNII feedback (in solar
   * masses) */
  double log10_SNII_max_mass_msun;

  /*! Number of type II supernovae per solar mass */
  float num_SNII_per_msun;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNII_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNII;

  /* ------------ Enrichment sampling properties ------------ */

  /*! Star age above which the enrichment will be downsampled (in internal
   * units) */
  double stellar_evolution_age_cut;

  /*! Number of time-steps in-between two enrichment events */
  int stellar_evolution_sampling_rate;

  /* ------------ Kinetic feedback properties --------------- */

  /* Velocity normalization */
  float FIRE_velocity_normalization;

  /* FIRE velocity slope */
  float FIRE_velocity_slope;

  /* Normalization for the mass loading curve */
  float FIRE_eta_normalization;

  /* The location (in internal mass units) where the break in the 
   * mass loading curve occurs */
  float FIRE_eta_break;

  /* The power-law slope of eta below FIRE_eta_break */
  float FIRE_eta_lower_slope;

  /* The power-law slope of eta above FIRE_eta_break */
  float FIRE_eta_upper_slope;

  /* Are we suppressing stellar feedback at high-z? */
  int early_wind_suppression_enabled;

  /* The minimum stellar mass normalization at high-z */
  float early_stellar_mass_norm;

  /* The scale factor when the suppression becomes negligible */
  float early_wind_suppression_scale_factor;

  /* The intensity of stellar feedback suppression at high-z */
  float early_wind_suppression_slope;

  /* The minimum galaxy stellar mass in internal units */
  float minimum_galaxy_stellar_mass;

  /* Added scatter to the wind velocities */
  float kick_velocity_scatter;

  /*! max decoupling time is (this factor) * current Hubble time */
  float wind_decouple_time_factor;

  /*! The internal energy corresponding to the cold gas temperature */
  float cold_wind_internal_energy;

  /* ------------ Common conversion factors --------------- */

  /*! Factor to convert km/s to internal units */
  float kms_to_internal;

  /*! Convert internal units to kpc */
  float length_to_kpc;

  /*! Convert internal time to Myr */
  float time_to_Myr;
};

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);

#endif /* SWIFT_SIMBA_FEEDBACK_PROPERTIES_H */
