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
 * star particle.
 *
 * @param sp The #spart.
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNII(const struct spart* sp,
                                     const struct feedback_props* props) {

  /* Note: For a Chabrier 2003 IMF and SNII going off between 6 and 100
   * M_sun, the first term is 0.017362 M_sun^-1 */
  return props->num_SNII_per_msun * sp->mass_init * props->mass_to_solar_mass;
}

/**
 * @brief Computes the number of supernovae of type Ia exploding for a given
 * star particle between time t and t+dt
 *
 * We follow Foerster et al. 2006, MNRAS, 368
 *
 * @param sp The #spart.
 * @param t0 The initial time (in Gyr).
 * @param t1 The final time (in Gyr).
 * @param props The properties of the stellar model.
 */
double eagle_feedback_number_of_SNIa(const struct spart* sp, const double t0,
                                     const double t1,
                                     const struct feedback_props* props) {

  /* The calculation is written as the integral between t0 and t1 of
   * eq. 3 of Schaye 2015 paper. */
  const double tau = props->SNIa_timescale_Gyr_inv;
  const double nu = props->SNIa_efficiency;
  const double num_SNe_per_Msun = nu * (exp(-t0 * tau) - exp(-t1 * tau));

  return num_SNe_per_Msun * sp->mass_init * props->mass_to_solar_mass;
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

  /* Smoothed metallicity (metal mass fraction) at birth time of the star */
  const double Z_smooth = sp->chemistry_data.smoothed_metal_mass_fraction_total;

  /* Physical density of the gas at the star's birth time */
  const double rho_birth = sp->birth_density;
  double n_birth = rho_birth * props->rho_to_n_cgs;

  /* Calculate f_E */
  const double Z_term = pow(max(Z_smooth, 1e-6) / Z_0, n_Z);
  const double n_term = pow(max(n_birth, 1e-6) / n_0, -n_n);
  const double denonimator = 1. + Z_term * n_term;

  return f_E_min + (f_E_max - f_E_min) / denonimator;
}

/**
 * @brief Compute the properties of the SNe feedback energy injection.
 *
 * Only does something if the particle reached the SNe age during this time
 * step.
 *
 * @param sp The star particle.
 * @param star_age Age of star at the beginning of the step in internal units.
 * @param dt Length of time-step in internal units.
 * @param ngb_gas_mass Total un-weighted mass in the star's kernel.
 * @param feedback_props The properties of the feedback model.
 */
INLINE static void compute_SNe_feedback(
    struct spart* sp, const double star_age, const double dt,
    const float ngb_gas_mass, const struct feedback_props* feedback_props) {

  /* Time after birth considered for SNII feedback (internal units) */
  const float SNII_wind_delay = feedback_props->SNII_wind_delay;

  /* Are we doing feedback this step? */
  if (star_age <= SNII_wind_delay && (star_age + dt) > SNII_wind_delay) {

#ifdef SWIFT_DEBUG_CHECKS
    if (sp->f_E != -1.f) error("Star has already done feedback!");
#endif

    /* Properties of the model (all in internal units) */
    const double delta_T =
        eagle_feedback_temperature_change(sp, feedback_props);
    const double N_SNe = eagle_feedback_number_of_SNII(sp, feedback_props);
    const double E_SNe = feedback_props->E_SNII;
    const double f_E = eagle_feedback_energy_fraction(sp, feedback_props);

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

    /* Store all of this in the star for delivery onto the gas */
    sp->f_E = f_E;
    sp->feedback_data.to_distribute.SNII_heating_probability = prob;
    sp->feedback_data.to_distribute.SNII_delta_u = delta_u;
  }
}

/**
 * @brief Find the bins and offset along the metallicity dimension of the
 * AGB yields table.
 *
 * @param iz_low (return) Lower index along the metallicity dimension.
 * @param iz_high (return) High index along the metallicity dimension.
 * @param dz (return) Offset between the metallicity bin and Z.
 * @param log10_Z log10 of the star metallicity (metal mass fraction).
 * @param props The properties of the feedback model.
 */
INLINE static void determine_bin_yield_AGB(int* iz_low, int* iz_high, float* dz,
                                           const float log10_Z,
                                           const struct feedback_props* props) {

  const double* AGB_Z = props->yield_AGB.metallicity;
  const int N_bins = eagle_feedback_AGB_N_metals;

  if (log10_Z > log10_min_metallicity) {

    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < (N_bins - 1) && log10_Z > AGB_Z[j + 1]; j++)
      ;

    /* Store the indices */
    *iz_low = j;
    *iz_high = *iz_low + 1;

    *iz_high = min(*iz_high, N_bins - 1);

    /* Compute offset */
    if ((log10_Z >= AGB_Z[0]) && (log10_Z <= AGB_Z[N_bins - 1])) {

      *dz = log10_Z - AGB_Z[*iz_low];
    } else {
      *dz = 0.f;
    }

    /* Normalize offset */
    const float delta_Z = AGB_Z[*iz_high] - AGB_Z[*iz_low];

    if (delta_Z > 0.f)
      *dz /= delta_Z;
    else
      *dz = 0.f;

  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0.f;
  }
}

/**
 * @brief Find the bins and offset along the metallicity dimension of the
 * SNII yields table.
 *
 * @param iz_low (return) Lower index along the metallicity dimension.
 * @param iz_high (return) High index along the metallicity dimension.
 * @param dz (return) Offset between the metallicity bin and Z.
 * @param log10_Z log10 of the star metallicity (metal mass fraction).
 * @param props The properties of the feedback model.
 */
INLINE static void determine_bin_yield_SNII(
    int* iz_low, int* iz_high, float* dz, const float log10_Z,
    const struct feedback_props* props) {

  const double* SNII_Z = props->yield_SNII.metallicity;
  const int N_bins = eagle_feedback_SNII_N_metals;

  if (log10_Z > log10_min_metallicity) {

    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < (N_bins - 1) && log10_Z > SNII_Z[j + 1]; j++)
      ;

    /* Store the indices */
    *iz_low = j;
    *iz_high = *iz_low + 1;

    *iz_high = min(*iz_high, N_bins - 1);

    /* Compute offset */
    if ((log10_Z >= SNII_Z[0]) && (log10_Z <= SNII_Z[N_bins - 1])) {

      *dz = log10_Z - SNII_Z[*iz_low];
    } else {
      *dz = 0.f;
    }

    /* Normalize offset */
    const float delta_Z = SNII_Z[*iz_high] - SNII_Z[*iz_low];

    if (delta_Z > 0.f)
      *dz = *dz / delta_Z;
    else
      *dz = 0.f;

  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0.f;
  }
}

/**
 * @brief compute enrichment and feedback due to SNIa. To do this compute the
 * number of SNIa that occur during the timestep, multiply by constants read
 * from tables.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param props properties of the feedback model
 * @param sp #spart we are computing feedback from
 * @param star_age_Gyr age of star in Gyr
 * @param dt_Gyr timestep dt in Gyr
 */
INLINE static void evolve_SNIa(const float log10_min_mass,
                               const float log10_max_mass,
                               const struct feedback_props* props,
                               struct spart* sp, float star_age_Gyr,
                               float dt_Gyr) {

  /* Check if we're outside the mass range for SNIa */
  if (log10_min_mass >= props->log10_SNIa_max_mass_msun) return;

  /* If the max mass is outside the mass range update it to be the maximum
   * and use updated values for the star's age and timestep in this function */
  if (log10_max_mass > props->log10_SNIa_max_mass_msun) {

    const float Z = sp->chemistry_data.metal_mass_fraction_total;
    const float max_mass = exp10f(props->log10_SNIa_max_mass_msun);
    const float lifetime_Gyr = lifetime_in_Gyr(max_mass, Z, props);

    dt_Gyr = star_age_Gyr + dt_Gyr - lifetime_Gyr;
    star_age_Gyr = lifetime_Gyr;
  }

  /* Compute the number of SNIa */
  const float num_SNIa = eagle_feedback_number_of_SNIa(
      sp, star_age_Gyr, star_age_Gyr + dt_Gyr, props);

  /* compute mass of each metal */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.to_distribute.metal_mass[i] +=
        num_SNIa * props->yield_SNIa_IMF_resampled[i] *
        props->solar_mass_to_mass;
  }

  /* Update the metallicity of the material released */
  sp->feedback_data.to_distribute.metal_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Update the metal mass produced */
  sp->feedback_data.to_distribute.total_metal_mass +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the mass produced by SNIa
   * Note: SNIa do not inject H or He so the mass injected is the same
   * as the metal mass injected. */
  sp->feedback_data.to_distribute.mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_total_metals_IMF_resampled *
      props->solar_mass_to_mass;

  /* Compute the iron mass produced */
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa +=
      num_SNIa * props->yield_SNIa_IMF_resampled[chemistry_element_Fe] *
      props->solar_mass_to_mass;

  /* Compute the energy to be injected */
  sp->feedback_data.to_distribute.energy += num_SNIa * props->E_SNIa;
}

/**
 * @brief compute enrichment and feedback due to SNII. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * Note for Matthieu: This function is poorly written and needs improving.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param props properties of the feedback model.
 * @param sp spart we are computing feedback from
 */
INLINE static void evolve_SNII(float log10_min_mass, float log10_max_mass,
                               float* stellar_yields,
                               const struct feedback_props* props,
                               struct spart* sp) {

  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

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
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d,
      high_index_2d;
  float dz = 0.;
  determine_bin_yield_SNII(&iz_low, &iz_high, &dz,
                           log10(sp->chemistry_data.metal_mass_fraction_total),
                           props);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_SNII_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      low_index_2d = row_major_index_2d(iz_low, mass_bin_index,
                                        eagle_feedback_SNII_N_metals,
                                        eagle_feedback_N_imf_bins);
      high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                         eagle_feedback_SNII_N_metals,
                                         eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (props->yield_SNII.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   props->yield_SNII.ejecta_IMF_resampled[low_index_2d]) +
          dz * (props->yield_SNII.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    props->yield_SNII.ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_SNII_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                       eagle_feedback_SNII_N_metals,
                                       eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (props->yield_SNII.total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        props->yield_SNII.ejecta_IMF_resampled[low_index_2d]) +
        dz * (props->yield_SNII.total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  props->yield_SNII.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star*/
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_SNII_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                       eagle_feedback_SNII_N_metals,
                                       eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * props->yield_SNII.ejecta_IMF_resampled[low_index_2d] +
        dz * props->yield_SNII.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* compute the total mass released */
  mass_released = metal_mass_released_total +
                  metal_mass_released[chemistry_element_H] +
                  metal_mass_released[chemistry_element_He];

  /* normalize the yields */
  if (mass_released > 0) {
    /* Set normalisation factor. Note additional multiplication by the star
     * initial mass as tables are per initial mass */
    const float norm_factor = sp->mass_init * mass_ejected / mass_released;

    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
    }
    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.mass_from_SNII +=
          sp->feedback_data.to_distribute.metal_mass[i];
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_SNII +=
        metal_mass_released_total * norm_factor;
  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
}

/**
 * @brief compute enrichment and feedback due to AGB. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * Note for Matthieu: This function is poorly written and needs improving.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param props Properties of the feedback model.
 * @param sp spart we are computing feedback for.
 */
INLINE static void evolve_AGB(const float log10_min_mass, float log10_max_mass,
                              float* stellar_yields,
                              const struct feedback_props* props,
                              struct spart* sp) {

  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* If mass at end of step is greater than tabulated lower bound for IMF, limit
   * it.*/
  if (log10_max_mass > props->log10_SNII_min_mass_msun)
    log10_max_mass = props->log10_SNII_min_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d,
      high_index_2d;
  float dz = 0.f;
  determine_bin_yield_AGB(&iz_low, &iz_high, &dz,
                          log10(sp->chemistry_data.metal_mass_fraction_total),
                          props);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = row_major_index_3d(
          iz_low, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      high_index_3d = row_major_index_3d(
          iz_high, elem, mass_bin_index, eagle_feedback_AGB_N_metals,
          chemistry_element_count, eagle_feedback_N_imf_bins);
      low_index_2d = row_major_index_2d(iz_low, mass_bin_index,
                                        eagle_feedback_AGB_N_metals,
                                        eagle_feedback_N_imf_bins);
      high_index_2d = row_major_index_2d(iz_high, mass_bin_index,
                                         eagle_feedback_AGB_N_metals,
                                         eagle_feedback_N_imf_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) * (props->yield_AGB.yield_IMF_resampled[low_index_3d] +
                      sp->chemistry_data.metal_mass_fraction[elem] *
                          props->yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
          dz * (props->yield_AGB.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    props->yield_AGB.ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] = integrate_imf(
        log10_min_mass, log10_max_mass, eagle_imf_integration_yield_weight,
        stellar_yields, props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d =
        row_major_index_2d(iz_high, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (props->yield_AGB.total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        props->yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
        dz * (props->yield_AGB.total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  props->yield_AGB.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d =
        row_major_index_2d(iz_low, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    high_index_2d =
        row_major_index_2d(iz_high, mass_bin_index, eagle_feedback_AGB_N_metals,
                           eagle_feedback_N_imf_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * props->yield_AGB.ejecta_IMF_resampled[low_index_2d] +
        dz * props->yield_AGB.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass,
                    eagle_imf_integration_yield_weight, stellar_yields, props);

  /* compute the total mass released */
  mass_released = metal_mass_released_total +
                  metal_mass_released[chemistry_element_H] +
                  metal_mass_released[chemistry_element_He];

  /* normalize the yields */
  if (mass_released > 0) {

    /* Set normalisation factor. Note additional multiplication by the stellar
     * initial mass as tables are per initial mass */
    const float norm_factor = sp->mass_init * mass_ejected / mass_released;

    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.metal_mass[i] +=
          metal_mass_released[i] * norm_factor;
      sp->feedback_data.to_distribute.mass_from_AGB +=
          metal_mass_released[i] * norm_factor;
    }
    sp->feedback_data.to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->feedback_data.to_distribute.metal_mass_from_AGB +=
        metal_mass_released_total * norm_factor;
  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
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
                               const struct unit_system* us, const float age,
                               const float dt) {

  TIMER_TIC;

  /* Allocate temporary array for calculating imf weights */
  float stellar_yields[eagle_feedback_N_imf_bins];

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_in_cgs = 1e9 * 365. * 24. * 3600.;
  const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  const float conversion_factor = time_to_cgs / Gyr_in_cgs;
  const float dt_Gyr = dt * conversion_factor;
  const float star_age_Gyr = age * conversion_factor;

  /* Get the metallicity */
  const float Z = sp->chemistry_data.metal_mass_fraction_total;

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

  /* Compute properties of the stochastic SNe feedback model. */
  if (feedback_props->with_SNe_feedback) {
    compute_SNe_feedback(sp, age, dt, ngb_gas_mass, feedback_props);
  }

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const float max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const float min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECK
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const float log10_max_dying_mass_Msun = log10f(max_dying_mass_Msun);
  const float log10_min_dying_mass_Msun = log10f(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                feedback_props, sp, star_age_Gyr, dt_Gyr);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
                stellar_yields, feedback_props, sp);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun,
               stellar_yields, feedback_props, sp);
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
 * By default, takes the values provided by the hydro.
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

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNe_feedback =
      parser_get_param_int(params, "EAGLEFeedback:use_SNe_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "EAGLEFeedback:use_SNIa_enrichment");

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:IMF_min_mass_Msun");

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);

  /* Properties of the SNII energy feedback model ------------------------- */

  /* Set the delay time before SNII occur */
  const float Gyr_in_cgs = 1e9 * 365 * 24 * 3600;
  fp->SNII_wind_delay =
      parser_get_param_float(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
      Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Read the temperature change to use in stochastic heating */
  fp->SNe_deltaT_desired =
      parser_get_param_float(params, "EAGLEFeedback:SNII_delta_T_K");
  fp->SNe_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNII_max_mass_Msun");

  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);

  /* Properties of the energy fraction model */
  fp->f_E_min =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_fraction_min");
  fp->f_E_max =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_fraction_max");
  fp->Z_0 =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_fraction_Z_0");
  fp->n_0_cgs = parser_get_param_double(
      params, "EAGLEFeedback:SNII_Energy_fraction_n_0_H_p_cm3");
  fp->n_n =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_fraction_n_n");
  fp->n_Z =
      parser_get_param_double(params, "EAGLEFeedback:SNII_Energy_fraction_n_Z");

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

  const double SNIa_max_mass_msun =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_max_mass_Msun");
  fp->log10_SNIa_max_mass_msun = log10(SNIa_max_mass_msun);

  /* Read SNIa timescale model parameters */
  fp->SNIa_efficiency =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_efficiency_p_Msun");
  fp->SNIa_timescale_Gyr =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_timescale_Gyr");
  fp->SNIa_timescale_Gyr_inv = 1.f / fp->SNIa_timescale_Gyr;

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "EAGLEFeedback:SNIa_Energy_erg");
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
