/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_PROPERTIES_H
#define SWIFT_EAGLE_FEEDBACK_PROPERTIES_H

#include "chemistry.h"
#include "hydro_properties.h"

/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {

  /* Yield table mass bins */
  double *mass;

  /* Yield table metallicity bins */
  double *metallicity;

  /* Array to store yield table resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /* Array to store yield table being read in */
  double *yield;

  /* Array to store table of ejecta resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /* Array to store table of ejecta being read in */
  double *ejecta;

  /* Array to store table of total mass released resampled by IMF mass bins */
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
 * @brief Properties of the EAGLE feedback model.
 */
struct feedback_props {

  /* ------------ Main operation modes ------------- */

  /*! Are we doing SNe feedback? */
  int with_SNe_feedback;

  /*! Are we doing AGB enrichment? */
  int with_AGB_enrichment;

  /*! Are we doing SNII enrichment? */
  int with_SNII_enrichment;

  /*! Are we doing SNIa enrichment? */
  int with_SNIa_enrichment;

  /* ------------ Energy injection ----------------- */

  /* Kinetic energy of SN ejecta per unit mass (check name with Richard)*/
  float ejecta_specific_thermal_energy;

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

  /*! Efficiency of the SNIa model */
  float SNIa_efficiency;

  /*! Time-scale of the SNIa decay function in Giga-years */
  float SNIa_timescale_Gyr;

  /*! Inverse of time-scale of the SNIa decay function in Giga-years */
  float SNIa_timescale_Gyr_inv;

  /*! Log 10 of the maximal mass used for SNIa feedback (in solar masses) */
  float log10_SNIa_max_mass_msun;

  /* ------------- Conversion factors --------------- */

  /*! Conversion factor from internal mass unit to solar mass */
  double mass_to_solar_mass;

  /*! Conversion factor from internal mass unit to solar mass */
  double solar_mass_to_mass;

  /*! Conversion factor from density in internal units to Hydrogen number
   * density in cgs */
  double rho_to_n_cgs;

  /*! Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /* ------------- Parameters for IMF --------------- */

  /*! Array to store calculated IMF */
  float *imf;

  /*! Arrays to store IMF mass bins */
  float *imf_mass_bin;

  /*! Arrays to store IMF mass bins (log10)*/
  float *imf_mass_bin_log10;

  /*! Minimal stellar mass considered by the IMF (in solar masses) */
  float imf_min_mass_msun;

  /*! Maximal stellar mass considered by the IMF (in solar masses) */
  float imf_max_mass_msun;

  /*! Log 10 of the minimal stellar mass considered by the IMF (in solar masses)
   */
  float log10_imf_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered by the IMF (in solar masses)
   */
  float log10_imf_max_mass_msun;

  /* ------------ SNe feedback properties ------------ */

  /*! Log 10 of the minimal stellar mass considered for SNII feedback (in solar
   * masses) */
  float log10_SNII_min_mass_msun;

  /*! Log 10 of the maximal stellar mass considered for SNII feedback (in solar
   * masses) */
  float log10_SNII_max_mass_msun;

  /*! Number of type II supernovae per solar mass */
  float num_SNII_per_msun;

  /*! Wind delay time for SNII */
  float SNII_wind_delay;

  /*! Temperature increase induced by SNe feedback */
  float SNe_deltaT_desired;

  /*! Energy released by one supernova type II in cgs units */
  double E_SNII_cgs;

  /*! Energy released by one supernova type II in internal units */
  float E_SNII;

  /*! Minimal energy fraction for supernova type II feedback */
  double f_E_min;

  /*! Maximal energy fraction for supernova type II feedback */
  double f_E_max;

  /*! Pivot point for the metallicity dependance of the feedback energy fraction
   * model */
  double Z_0;

  /*! Pivot point for the density dependance of the feedback energy fraction
   * model */
  double n_0_cgs;

  /*! Slope of the density dependance of the feedback energy fraction model */
  double n_n;

  /*! Slope of the metallicity dependance of the feedback energy fraction model
   */
  double n_Z;
};

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);

#endif /* SWIFT_EAGLE_FEEDBACK_PROPERTIES_H */
