/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_FEEDBACK_ENRICHMENT_H
#define SWIFT_EAGLE_FEEDBACK_ENRICHMENT_H

#include "imf.h"
#include "interpolate.h"

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle assuming that all the SNII stars go off at once.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNII(const struct spart* sp,
                                     const struct feedback_props* props) {

  /* Note: For a Chabrier 2003 IMF and SNII going off
   * - between 6 and 100 M_sun, the first term is 0.017362 M_sun^-1 (EAGLE)
   * - between 8 and 100 M_sun, the first term is 0.011801 M_sun^-1 (EAGLE-XL)
   */
  return props->num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the number of supernovae of type II exploding for a given
 * star particle between two mass limits
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 * @param min_dying_mass_Msun Minimal star mass dying this step (in solar
 * masses).
 * @param max_dying_mass_Msun Maximal star mass dying this step (in solar
 * masses).
 */
double eagle_feedback_number_of_sampled_SNII(const struct spart* sp,
                                             const struct feedback_props* props,
                                             const double min_dying_mass_Msun,
                                             const double max_dying_mass_Msun) {

  /* The max dying star mass is below the SNII mass window
   * --> No SNII */
  if (max_dying_mass_Msun < props->SNII_min_mass_msun) return 0.;

  /* The min dying star mass is above the SNII mass window
   * --> No SNII */
  if (min_dying_mass_Msun > props->SNII_max_mass_msun) return 0.;

  /* Ok, we have some overlap with the SNII mass window. */

  double log10_min_mass_Msun = -10.;
  double log10_max_mass_Msun = -10.;

  /* The min dying star mass dies inside the SNII mass window */
  if (min_dying_mass_Msun <= props->SNII_max_mass_msun &&
      min_dying_mass_Msun > props->SNII_min_mass_msun) {

    /* Now, check the max dying star mass */

    /* The max dying star mass is also inside the SNII mass window */
    if (max_dying_mass_Msun <= props->SNII_max_mass_msun) {
      log10_min_mass_Msun = log10(min_dying_mass_Msun);
      log10_max_mass_Msun = log10(max_dying_mass_Msun);
    }

    /* The max dying star is above the SNII mass window */
    else {
      log10_min_mass_Msun = log10(min_dying_mass_Msun);
      log10_max_mass_Msun = props->log10_SNII_max_mass_msun;
    }

  }

  /* The min dying star mass dies below the SNII mass window */
  else if (min_dying_mass_Msun <= props->SNII_min_mass_msun) {

    /* Now, check the max dying star mass */

    /* The max dying star mass is inside the SNII mass window */
    if (max_dying_mass_Msun > props->SNII_min_mass_msun &&
        max_dying_mass_Msun <= props->SNII_max_mass_msun) {
      log10_min_mass_Msun = props->log10_SNII_min_mass_msun;
      log10_max_mass_Msun = log10(max_dying_mass_Msun);
    }

    /* The max dying star is above the SNII mass window */
    else if (max_dying_mass_Msun > props->SNII_max_mass_msun) {
      log10_min_mass_Msun = props->log10_SNII_min_mass_msun;
      log10_max_mass_Msun = props->log10_SNII_max_mass_msun;
    }

    /* The max dying star mass is also below the SNII mass window */
    else {

      /* We already excluded this at the star of the function */
#ifdef SWIFT_DEBUG_CHECKS
      error("Error in the logic");
#endif
    }
  }

  /* The min dying star mass dies above the SNII mass window */
  else {

    /* We already excluded this at the star of the function */
#ifdef SWIFT_DEBUG_CHECKS
    error("Error in the logic");
#endif
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (log10_min_mass_Msun == -10. || log10_max_mass_Msun == -10.)
    error("Something went wrong in the calculation of the number of SNII.");
#endif

  /* Calculate how many supernovae have exploded in this timestep
   * by integrating the IMF between the bounds we chose */
  const double num_SNII_per_msun =
      integrate_imf(log10_min_mass_Msun, log10_max_mass_Msun,
                    eagle_imf_integration_no_weight, NULL, props);

  return num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t0 and t1
 *
 * @param M_init The inital mass of the star particle in internal units.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNIa(const double M_init, const double t0,
                                     const double t1,
                                     const struct feedback_props* props) {

  double num_SNIa_per_Msun = 0.;

#ifdef SWIFT_DEBUG_CHECKS
  if (t1 < t0) error("Negative time range!");
  if (t0 < props->SNIa_DTD_delay_Gyr)
    error("Initial time smaller than the delay time!");
#endif

  switch (props->SNIa_DTD) {

    case eagle_feedback_SNIa_DTD_exponential: {

      /* We follow Foerster et al. 2006, MNRAS, 368 */

      /* The calculation is written as the integral between t0 and t1 of
       * eq. 3 of Schaye 2015 paper. */
      const double tau = props->SNIa_DTD_exp_timescale_Gyr_inv;
      const double nu = props->SNIa_DTD_exp_norm;
      num_SNIa_per_Msun = nu * (exp(-t0 * tau) - exp(-t1 * tau));
      break;
    }

    case eagle_feedback_SNIa_DTD_power_law: {

      /* We follow Graur et al. 2011, MNRAS, 417 */

      const double norm = props->SNIa_DTD_power_law_norm;
      num_SNIa_per_Msun = norm * log(t1 / t0);
      break;
    }

    default: {
      num_SNIa_per_Msun = 0.;
      error("Invalid choice of SNIa delay time distribution!");
    }
  }

  return num_SNIa_per_Msun * M_init * props->mass_to_solar_mass;
}

/**
 * @brief Find the bins and offset along the metallicity dimension of the
 * yields table.
 *
 * Note for Matthieu: In the normal EAGLE case, the arrays are of length 3 and
 * 5 for the AGB and SNII channels respectivly. Can we simplify this?
 *
 * @param index_Z_low (return) Lower index along the metallicity dimension.
 * @param index_Z_high (return) High index along the metallicity dimension.
 * @param dZ (return) Offset between the metallicity bin and Z.
 * @param log10_Z log10 of the star metallicity (metal mass fraction).
 * @param log_Z_bins bin of log10 of the metallicities in the table for this
 * enrichment channel.
 * @param N_bins The number of metallicity bins for this enrichment channel.
 */
INLINE static void determine_bin_yields(int* index_Z_low, int* index_Z_high,
                                        float* dZ, const float log10_Z,
                                        const double* const log_Z_bins,
                                        const int N_bins) {

  if (log10_Z > log10_min_metallicity) {

    /* Find metallicity bin which contains the star's metallicity */
    int j = 0;
    while (j < (N_bins - 1) && log10_Z > log_Z_bins[j + 1]) {
      j++;
    }

    /* Store the indices */
    *index_Z_low = j;
    *index_Z_high = j + 1;

    *index_Z_high = min(*index_Z_high, N_bins - 1);

    /* Compute offset from index_Z_low in the table*/
    if ((log10_Z >= log_Z_bins[0]) && (log10_Z <= log_Z_bins[N_bins - 1])) {
      *dZ = log10_Z - log_Z_bins[*index_Z_low];
    } else {
      *dZ = 0.f;
    }

    /* Normalize offset */
    const float delta_Z = log_Z_bins[*index_Z_high] - log_Z_bins[*index_Z_low];

    if (delta_Z > 0.f) {
      *dZ /= delta_Z;
    } else {
      *dZ = 0.f;
    }

  } else {
    *index_Z_low = 0;
    *index_Z_high = 0;
    *dZ = 0.f;
  }
}

/**
 * @brief compute enrichment and feedback due to SNIa. To do this compute the
 * number of SNIa that occur during the timestep, multiply by constants read
 * from tables.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param M_init The initial mass of the star particle (in internal units).
 * @param Z The total metallicity of the star (metal mass fraction).
 * @param props Properties of the feedback model.
 * @param star_age_Gyr age of star in Gyr
 * @param dt_Gyr timestep dt in Gyr
 * @param feedback_data (return) The #feedback_spart_data to fill with things to
 * distribute to the gas.
 */
INLINE static void evolve_SNIa(
    const double log10_min_mass, const double log10_max_mass,
    const double M_init, const double Z, const struct feedback_props* props,
    double star_age_Gyr, const double dt_Gyr,
    struct feedback_spart_data* const feedback_data) {

  const double star_age_end_step_Gyr = star_age_Gyr + dt_Gyr;

  /* Check if we're outside the mass range for SNIa */
  if (star_age_end_step_Gyr < props->SNIa_DTD_delay_Gyr) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (dt_Gyr < 0.) error("Negative time-step length!");
  if (star_age_Gyr < 0.) error("Negative age!");
#endif

  /* Only consider stars beyond the minimal age for SNIa */
  star_age_Gyr = max(star_age_Gyr, props->SNIa_DTD_delay_Gyr);

  /* Compute the number of SNIa */
  const float num_SNIa = eagle_feedback_number_of_SNIa(
      M_init, star_age_Gyr, star_age_end_step_Gyr, props);

  /* Compute mass of each metal */
  for (int i = 0; i < chemistry_element_count; i++) {
    feedback_data->to_distribute.metal_mass[i] +=
        num_SNIa * props->yield_SNIa_IMF_resampled[i] *
        props->solar_mass_to_mass;
  }

  /* Update the metallicity of the material released */
  feedback_data->to_distribute.metal_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Update the metal mass produced */
  feedback_data->to_distribute.total_metal_mass +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the mass produced by SNIa
   * Note: SNIa do not inject H or He so the mass injected is the same
   * as the metal mass injected. */
  feedback_data->to_distribute.mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the iron mass produced */
  feedback_data->to_distribute.Fe_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_IMF_resampled[chemistry_element_Fe] *
      props->solar_mass_to_mass;

  /* Compute the energy to be injected */
  if (props->with_SNIa_feedback) {
    feedback_data->to_distribute.energy += num_SNIa * props->E_SNIa;
  }
}

/**
 * @brief compute enrichment and feedback due to SNII. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param M_init The initial mass of the star particle (in internal units).
 * @param Z The total metallicity of the star (metal mass fraction).
 * @param abundances The individual metal abundances (mass fractions) of the
 * star.
 * @param props Properties of the feedback model.
 * @param feedback_data (return) The #feedback_spart_data to fill with things to
 * distribute to the gas.
 */
INLINE static void evolve_SNII(
    double log10_min_mass, double log10_max_mass, const double M_init,
    const double Z, const float* const abundances,
    const struct feedback_props* props,
    struct feedback_spart_data* const feedback_data) {

  /* Pull out common arrays */

  /* Metal mass produced by the star */
  const double* const total_yields =
      props->yield_SNII.total_metals_IMF_resampled;

  /* Individual elements produced by the star */
  const double* const metal_yields = props->yield_SNII.yield_IMF_resampled;

  /* Elements already in the stars that are ejected */
  const double* const ejecta = props->yield_SNII.ejecta_IMF_resampled;

  /* If mass at beginning of step is less than tabulated lower bound for IMF,
   * limit it.*/
  if (log10_min_mass < props->log10_SNII_min_mass_msun)
    log10_min_mass = props->log10_SNII_min_mass_msun;

  /* If mass at end of step is greater than tabulated upper bound for IMF, limit
   * it.*/
  if (log10_max_mass > props->log10_SNII_max_mass_msun)
    log10_max_mass = props->log10_SNII_max_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  int low_imf_mass_bin_index, high_imf_mass_bin_index;
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int index_Z_lo = 0, index_Z_hi = 0;
  float dZ = 0.;
  determine_bin_yields(&index_Z_lo, &index_Z_hi, &dZ, log10(Z),
                       props->yield_SNII.metallicity,
                       eagle_feedback_SNII_N_metals);

  /* Allocate temporary array for calculating imf weights */
  double stellar_yields[eagle_feedback_N_imf_bins];

  /*******************************
   * Compute metal mass produced *
   *******************************/
  double metal_mass_released[chemistry_element_count];

  /* Loop over all the elements */
  for (int elem = 0; elem < chemistry_element_count; elem++) {

    /* Loop over all the relevent IMF mass bins */
    for (int mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

      const int lo_index_3d = row_major_index_3d(
          index_Z_lo, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      const int hi_index_3d = row_major_index_3d(
          index_Z_hi, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);

      const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                                 eagle_feedback_SNII_N_metals,
                                                 eagle_feedback_N_imf_bins);
      const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                                 eagle_feedback_SNII_N_metals,
                                                 eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1.f - dZ) * (metal_yields[lo_index_3d] +
                        abundances[elem] * ejecta[lo_index_2d]) +
          (0.f + dZ) * (metal_yields[hi_index_3d] +
                        abundances[elem] * ejecta[hi_index_2d]);
    }
    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /*************************************
   * Compute total metal mass produced *
   *************************************/
  for (int mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

    const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                               eagle_feedback_SNII_N_metals,
                                               eagle_feedback_N_imf_bins);
    const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                               eagle_feedback_SNII_N_metals,
                                               eagle_feedback_N_imf_bins);

    stellar_yields[mass_bin_index] =
        (1.f - dZ) * (total_yields[lo_index_2d] + Z * ejecta[lo_index_2d]) +
        (0.f + dZ) * (total_yields[hi_index_2d] + Z * ejecta[hi_index_2d]);
  }
  double metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /************************************************
   * Compute the total mass ejected from the star *
   ************************************************/
  for (int mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

    const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                               eagle_feedback_SNII_N_metals,
                                               eagle_feedback_N_imf_bins);
    const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                               eagle_feedback_SNII_N_metals,
                                               eagle_feedback_N_imf_bins);

    stellar_yields[mass_bin_index] =
        (1 - dZ) * ejecta[lo_index_2d] + dZ * ejecta[hi_index_2d];
  }
  const double mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* Zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);
  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total mass released */
  const double mass_released = metal_mass_released_total +
                               metal_mass_released[chemistry_element_H] +
                               metal_mass_released[chemistry_element_He];

#ifdef SWIFT_DEBUG_CHECKS
  if (mass_released <= 0) {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
#endif

  /* Set normalisation factor. Note additional multiplication by the star
   * initial mass as tables are per initial mass */
  const double norm_factor = M_init * mass_ejected / mass_released;

  /* Store what we want to distribute */
  for (int i = 0; i < chemistry_element_count; i++) {
    feedback_data->to_distribute.metal_mass[i] +=
        metal_mass_released[i] * norm_factor;
    feedback_data->to_distribute.mass_from_SNII +=
        metal_mass_released[i] * norm_factor;
  }
  feedback_data->to_distribute.total_metal_mass +=
      metal_mass_released_total * norm_factor;
  feedback_data->to_distribute.metal_mass_from_SNII +=
      metal_mass_released_total * norm_factor;
}

/**
 * @brief compute enrichment and feedback due to AGB. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param M_init The initial mass of the star particle (in internal units).
 * @param Z The total metallicity of the star (metal mass fraction).
 * @param abundances The individual metal abundances (mass fractions) of the
 * star.
 * @param props Properties of the feedback model.
 * @param feedback_data (return) The #feedback_spart_data to fill with things to
 * distribute to the gas.
 */
INLINE static void evolve_AGB(const double log10_min_mass,
                              double log10_max_mass, const double M_init,
                              const double Z, const float* const abundances,
                              const struct feedback_props* props,
                              struct feedback_spart_data* const feedback_data) {

  /* Pull out common arrays */

  /* Metal mass produced by the star */
  const double* const total_yields =
      props->yield_AGB.total_metals_IMF_resampled;

  /* Individual elements produced by the star */
  const double* const metal_yields = props->yield_AGB.yield_IMF_resampled;

  /* Elements already in the stars that are ejected */
  const double* const ejecta = props->yield_AGB.ejecta_IMF_resampled;

  /* If mass at end of step is greater than tabulated lower bound for IMF, limit
   * it.*/
  if (log10_max_mass > props->log10_SNII_min_mass_msun)
    log10_max_mass = props->log10_SNII_min_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  int low_imf_mass_bin_index, high_imf_mass_bin_index;
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int index_Z_lo = 0, index_Z_hi = 0;
  float dZ = 0.f;
  determine_bin_yields(&index_Z_lo, &index_Z_hi, &dZ, log10(Z),
                       props->yield_AGB.metallicity,
                       eagle_feedback_AGB_N_metals);

  /* Allocate temporary array for calculating imf weights */
  double stellar_yields[eagle_feedback_N_imf_bins];

  /*******************************
   * Compute metal mass produced *
   *******************************/
  double metal_mass_released[chemistry_element_count];

  /* Loop over all the elements */
  for (int elem = 0; elem < chemistry_element_count; elem++) {

    /* Loop over all the relevent IMF mass bins */
    for (int mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

      const int lo_index_3d = row_major_index_3d(
          index_Z_lo, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      const int hi_index_3d = row_major_index_3d(
          index_Z_hi, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);

      const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                                 eagle_feedback_AGB_N_metals,
                                                 eagle_feedback_N_imf_bins);
      const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                                 eagle_feedback_AGB_N_metals,
                                                 eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1.f - dZ) * (metal_yields[lo_index_3d] +
                        abundances[elem] * ejecta[lo_index_2d]) +
          (0.f + dZ) * (metal_yields[hi_index_3d] +
                        abundances[elem] * ejecta[hi_index_2d]);
    }

    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /*************************************
   * Compute total metal mass produced *
   *************************************/
  for (int mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

    const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                               eagle_feedback_AGB_N_metals,
                                               eagle_feedback_N_imf_bins);
    const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                               eagle_feedback_AGB_N_metals,
                                               eagle_feedback_N_imf_bins);

    stellar_yields[mass_bin_index] =
        (1.f - dZ) * (total_yields[lo_index_2d] + Z * ejecta[lo_index_2d]) +
        (0.f + dZ) * (total_yields[hi_index_2d] + Z * ejecta[hi_index_2d]);
  }

  double metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /************************************************
   * Compute the total mass ejected from the star *
   ************************************************/

  for (int mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {

    const int lo_index_2d = row_major_index_2d(index_Z_lo, mass_bin_index,
                                               eagle_feedback_AGB_N_metals,
                                               eagle_feedback_N_imf_bins);
    const int hi_index_2d = row_major_index_2d(index_Z_hi, mass_bin_index,
                                               eagle_feedback_AGB_N_metals,
                                               eagle_feedback_N_imf_bins);

    stellar_yields[mass_bin_index] =
        (1.f - dZ) * ejecta[lo_index_2d] + dZ * ejecta[hi_index_2d];
  }
  const double mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* Zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);
  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total mass released */
  const double mass_released = metal_mass_released_total +
                               metal_mass_released[chemistry_element_H] +
                               metal_mass_released[chemistry_element_He];

#ifdef SWIFT_DEBUG_CHECKS
  if (mass_released <= 0) {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
#endif

  /* Set normalisation factor. Note additional multiplication by the stellar
   * initial mass as tables are per initial mass */
  const double norm_factor = M_init * mass_ejected / mass_released;

  for (int i = 0; i < chemistry_element_count; i++) {
    feedback_data->to_distribute.metal_mass[i] +=
        metal_mass_released[i] * norm_factor;
    feedback_data->to_distribute.mass_from_AGB +=
        metal_mass_released[i] * norm_factor;
  }
  feedback_data->to_distribute.total_metal_mass +=
      metal_mass_released_total * norm_factor;
  feedback_data->to_distribute.metal_mass_from_AGB +=
      metal_mass_released_total * norm_factor;
}

/**
 * @brief Zero pointers in yield_table structs
 *
 * @param table yield_table struct in which pointers to tables
 * set to NULL
 */
void zero_yield_table_pointers(struct yield_table* table) {

  table->mass = NULL;
  table->metallicity = NULL;
  table->yield_IMF_resampled = NULL;
  table->yield = NULL;
  table->ejecta_IMF_resampled = NULL;
  table->ejecta = NULL;
  table->total_metals_IMF_resampled = NULL;
  table->total_metals = NULL;
}

/**
 * @brief Restore feedback tables (if applicable) after
 * restart
 *
 * @param fp the #feedback_props structure
 */
void feedback_restore_tables(struct feedback_props* fp) {

  init_imf(fp);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);
}

#endif /* SWIFT_EAGLE_FEEDBACK_ENRICHMENT_H */
