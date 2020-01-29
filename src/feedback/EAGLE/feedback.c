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

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "hydro_properties.h"
#include "imf.h"
#include "inline.h"
#include "interpolate.h"
#include "timers.h"
#include "yield_tables.h"

/**
 * @brief Return the change in temperature (in internal units) to apply to a
 * gas particle affected by SNe feedback.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
double eagle_feedback_temperature_change(const struct spart* sp,
                                         const struct feedback_props* props) {

  /* In the EAGLE REF model, the change of temperature is constant */
  return props->SNe_deltaT_desired;
}

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
 * @brief Computes the fraction of the available super-novae energy to
 * inject for a given event.
 *
 * Note that the fraction can be > 1.
 *
 * We use equation 7 of Schaye et al. 2015.
 *
 * @param sp The #spart.
 * @param props The properties of the feedback model.
 */
double eagle_feedback_energy_fraction(const struct spart* sp,
                                      const struct feedback_props* props) {

  /* Model parameters */
  const double f_E_max = props->f_E_max;
  const double f_E_min = props->f_E_min;
  const double Z_0 = props->Z_0;
  const double n_0 = props->n_0_cgs;
  const double n_Z = props->n_Z;
  const double n_n = props->n_n;

  /* Star properties */

  /* Metallicity (metal mass fraction) at birth time of the star */
  const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(sp);

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->birth_density;
  const double n_birth = rho_birth * props->rho_to_n_cgs;

  /* Calculate f_E */
  const double Z_term = pow(max(Z, 1e-6) / Z_0, n_Z);
  const double n_term = pow(n_birth / n_0, -n_n);
  const double denonimator = 1. + Z_term * n_term;

  return f_E_min + (f_E_max - f_E_min) / denonimator;
}

/**
 * @brief Compute the properties of the SNII stochastic feedback energy
 * injection.
 *
 * Only does something if the particle reached the SNII age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel.
 * @param feedback_props The properties of the feedback model.
 * @param min_dying_mass_Msun Minimal star mass dying this step (in solar
 * masses).
 * @param max_dying_mass_Msun Maximal star mass dying this step (in solar
 * masses).
 */
INLINE static void compute_SNII_feedback(
    struct spart* sp, const double star_age, const double dt,
    const float ngb_gas_mass, const struct feedback_props* feedback_props,
    const double min_dying_mass_Msun, const double max_dying_mass_Msun) {

  /* Are we sampling the delay function or using a fixed delay? */
  const int SNII_sampled_delay = feedback_props->SNII_sampled_delay;

  /* Time after birth considered for SNII feedback (internal units)
   * when using a fixed delay */
  const double SNII_wind_delay = feedback_props->SNII_wind_delay;

  /* Are we doing feedback this step?
   * Note that since the ages are calculated using an interpolation table we
   * must allow some tolerance here*/
  if ((SNII_sampled_delay) || (star_age <= SNII_wind_delay &&
                               (star_age + 1.001 * dt) > SNII_wind_delay)) {

    /* Make sure a star does not do feedback twice
     * when using a fixed delay! */
    if (!SNII_sampled_delay && sp->f_E != -1.f) {
#ifdef SWIFT_DEBUG_CHECKS
      message("Star has already done feedback! sp->id=%lld age=%e d=%e", sp->id,
              star_age, dt);
#endif
      return;
    }

    /* Properties of the model (all in internal units) */
    const double delta_T =
        eagle_feedback_temperature_change(sp, feedback_props);
    const double E_SNe = feedback_props->E_SNII;
    const double f_E = eagle_feedback_energy_fraction(sp, feedback_props);

    /* Number of SNe at this time-step */
    double N_SNe;
    if (SNII_sampled_delay) {
      N_SNe = eagle_feedback_number_of_sampled_SNII(
          sp, feedback_props, min_dying_mass_Msun, max_dying_mass_Msun);
    } else {
      N_SNe = eagle_feedback_number_of_SNII(sp, feedback_props);
    }

    /* Abort if there are no SNe exploding this step */
    if (N_SNe == 0.) return;

    /* Conversion factor from T to internal energy */
    const double conv_factor = feedback_props->temp_to_u_factor;

    /* Calculate the default heating probability */
    double prob = f_E * E_SNe * N_SNe / (conv_factor * delta_T * ngb_gas_mass);

    /* Calculate the change in internal energy of the gas particles that get
     * heated */
    double delta_u;
    if (prob <= 1.) {

      /* Normal case */
      delta_u = delta_T * conv_factor;

    } else {

      /* Special case: we need to adjust the energy irrespective of the
         desired deltaT to ensure we inject all the available energy. */

      prob = 1.;
      delta_u = f_E * E_SNe * N_SNe / ngb_gas_mass;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (f_E < feedback_props->f_E_min || f_E > feedback_props->f_E_max)
      error("f_E is not in the valid range! f_E=%f sp->id=%lld", f_E, sp->id);
#endif

    /* Store all of this in the star for delivery onto the gas */
    sp->f_E = f_E;
    sp->feedback_data.to_distribute.SNII_heating_probability = prob;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
  }
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
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");

  if (sp->count_since_last_enrichment != 0)
    error("Computing feedback on a star that should not");
#endif

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_in_cgs = 1e9 * 365.25 * 24. * 3600.;
  const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const double conversion_factor = time_to_cgs / Gyr_in_cgs;
  const double dt_Gyr = dt * conversion_factor;
  const double star_age_Gyr = age * conversion_factor;

  /* Get the birth mass of the star */
  const double M_init = sp->mass_init;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const double Z = max(chemistry_get_total_metal_mass_fraction_for_feedback(sp),
                       exp10(log10_min_metallicity));

  /* Get the individual abundances (mass fractions at birth time) */
  const float* const abundances =
      chemistry_get_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const double max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const double min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Compute properties of the stochastic SNII feedback model. */
  if (feedback_props->with_SNII_feedback) {
    compute_SNII_feedback(sp, age, dt, ngb_gas_mass, feedback_props,
                          min_dying_mass_Msun, max_dying_mass_Msun);
  }

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const double log10_max_dying_mass_Msun = log10(max_dying_mass_Msun);
  const double log10_min_dying_mass_Msun = log10(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                feedback_props, star_age_Gyr, dt_Gyr, &sp->feedback_data);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                abundances, feedback_props, &sp->feedback_data);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
               abundances, feedback_props, &sp->feedback_data);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo) {

  const double Gyr_in_cgs = 1.0e9 * 365.25 * 24. * 3600.;

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNII_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_feedback");

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_enrichment");

  if (fp->with_SNIa_feedback && !fp->with_SNIa_enrichment) {
    error("Cannot run with SNIa feedback without SNIa enrichment.");
  }

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  /* Are we sampling the SNII lifetimes for feedback or using a fixed delay? */
  fp->SNII_sampled_delay =
      parser_get_param_int(params, "EAGLEFeedback:SNII_sampled_delay");

  if (!fp->SNII_sampled_delay) {

    /* Set the delay time before SNII occur */
    fp->SNII_wind_delay =
        parser_get_param_double(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
        Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  }

  /* Read the temperature change to use in stochastic heating */
  fp->SNe_deltaT_desired =
      parser_get_param_float(params, "EAGLEFeedback:SNII_delta_T_K");
  fp->SNe_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->SNII_min_mass_msun = SNII_min_mass_msun;
  fp->SNII_max_mass_msun = SNII_max_mass_msun;
  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Properties of the energy fraction model */
  fp->f_E_min =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_min");
  fp->f_E_max =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_max");
  fp->Z_0 =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_Z_0");
  fp->n_0_cgs = parser_get_param_double(
      params, "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3");
  fp->n_n =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_n_n");
  fp->n_Z =
      parser_get_param_double(params, "EAGLEFeedback:SNII_energy_fraction_n_Z");

  /* Check that it makes sense. */
  if (fp->f_E_max < fp->f_E_min) {
    error("Can't have the maximal energy fraction smaller than the minimal!");
  }

  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "EAGLEFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  fp->SNIa_DTD_delay_Gyr =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_DTD_delay_Gyr");

  char temp[32] = {0};
  parser_get_param_string(params, "EAGLEFeedback:SNIa_DTD", temp);

  if (strcmp(temp, "Exponential") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_exponential;

    /* Read SNIa exponential DTD model parameters */
    fp->SNIa_DTD_exp_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_norm_p_Msun");
    fp->SNIa_DTD_exp_timescale_Gyr = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_exp_timescale_Gyr");
    fp->SNIa_DTD_exp_timescale_Gyr_inv = 1.f / fp->SNIa_DTD_exp_timescale_Gyr;

  } else if (strcmp(temp, "PowerLaw") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_power_law;

    /* Read SNIa power-law DTD model parameters */
    fp->SNIa_DTD_power_law_norm = parser_get_param_float(
        params, "EAGLEFeedback:SNIa_DTD_power_law_norm_p_Msun");

    /* Renormalize everything such that the integral converges to
       'SNIa_DTD_power_law_norm' over 13.6 Gyr. */
    fp->SNIa_DTD_power_law_norm /= log(13.6 / fp->SNIa_DTD_delay_Gyr);

  } else {
    error("Invalid SNIa DTD model: '%s'", temp);
  }

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "EAGLEFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the enrichment down-sampling ----------------------------- */

  fp->stellar_evolution_age_cut =
      parser_get_param_double(params,
                              "EAGLEFeedback:stellar_evolution_age_cut_Gyr") *
      Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "EAGLEFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Check that we are not downsampling before reaching the SNII delay */
  if (!fp->SNII_sampled_delay &&
      fp->stellar_evolution_age_cut < fp->SNII_wind_delay)
    error(
        "Time at which the enrichment downsampling stars is lower than the "
        "SNII wind delay!");

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Calculate number of type II SN per unit solar mass based on our choice
   * of IMF and integration limits for type II SNe.
   * Note: No weighting by yields here. */
  fp->num_SNII_per_msun =
      integrate_imf(fp->log10_SNII_min_mass_msun, fp->log10_SNII_max_mass_msun,
                    eagle_imf_integration_no_weight,
                    /*(stellar_yields=)*/ NULL, fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "EAGLEFeedback:filename",
                          fp->yield_table_path);

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

  message("initialized stellar feedback");
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

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param feedback_props the feedback data structure.
 */
void feedback_clean(struct feedback_props* feedback_props) {

  swift_free("imf-tables", feedback_props->imf);
  swift_free("imf-tables", feedback_props->imf_mass_bin);
  swift_free("imf-tables", feedback_props->imf_mass_bin_log10);
  swift_free("feedback-tables", feedback_props->yields_SNIa);
  swift_free("feedback-tables", feedback_props->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_AGB.mass);
  swift_free("feedback-tables", feedback_props->yield_AGB.metallicity);
  swift_free("feedback-tables", feedback_props->yield_AGB.yield);
  swift_free("feedback-tables", feedback_props->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_AGB.ejecta);
  swift_free("feedback-tables", feedback_props->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_AGB.total_metals);
  swift_free("feedback-tables",
             feedback_props->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_SNII.mass);
  swift_free("feedback-tables", feedback_props->yield_SNII.metallicity);
  swift_free("feedback-tables", feedback_props->yield_SNII.yield);
  swift_free("feedback-tables", feedback_props->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_SNII.ejecta);
  swift_free("feedback-tables",
             feedback_props->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", feedback_props->yield_SNII.total_metals);
  swift_free("feedback-tables",
             feedback_props->yield_SNII.total_metals_IMF_resampled);
  swift_free("feedback-tables", feedback_props->lifetimes.mass);
  swift_free("feedback-tables", feedback_props->lifetimes.metallicity);
  swift_free("feedback-tables", feedback_props->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(feedback_props->lifetimes.dyingtime[i]);
  }
  free(feedback_props->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(feedback_props->SNIa_element_names[i]);
  }
  free(feedback_props->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(feedback_props->SNII_element_names[i]);
  }
  free(feedback_props->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(feedback_props->AGB_element_names[i]);
  }
  free(feedback_props->AGB_element_names);
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  feedback_restore_tables(feedback);
}
