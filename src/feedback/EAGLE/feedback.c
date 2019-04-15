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


#include "feedback.h"

#include "hydro_properties.h"
#include "imf.h"
#include "interpolate.h"
#include "yield_tables.h"

/**
 * @brief determine which metallicity bin star belongs to for AGB, compute bin
 * indices and offsets
 *
 * @param iz_low Pointer to index of metallicity bin to which the star belongs
 * (to be calculated in this function)
 * @param iz_high Pointer to index of metallicity bin to above the star's
 * metallicity (to be calculated in this function)
 * @param dz metallicity bin offset
 * @param log10_metallicity log10 of star metallicity
 * @param star_properties stars_props data structure
 */
inline static void determine_bin_yield_AGB(
    int* iz_low, int* iz_high, float* dz, float log10_metallicity,
    const struct feedback_props* feedback_props) {

  // MATTHIEU
  
  /* if (log10_metallicity > log10_min_metallicity) { */
  /*   /\* Find metallicity bin which contains the star's metallicity *\/ */
  /*   int j; */
  /*   for (j = 0; j < star_properties->feedback.AGB_n_z - 1 && */
  /*               log10_metallicity > */
  /*                   star_properties->feedback.yield_AGB.metallicity[j + 1]; */
  /*        j++) */
  /*     ; */
  /*   *iz_low = j; */
  /*   *iz_high = *iz_low + 1; */

  /*   /\* Compute offset *\/ */
  /*   if (log10_metallicity >= */
  /*           star_properties->feedback.yield_AGB.metallicity[0] && */
  /*       log10_metallicity <= */
  /*           star_properties->feedback.yield_AGB */
  /*               .metallicity[star_properties->feedback.AGB_n_z - 1]) */
  /*     *dz = log10_metallicity - */
  /*           star_properties->feedback.yield_AGB.metallicity[*iz_low]; */
  /*   else */
  /*     *dz = 0; */

  /*   /\* Normalize offset *\/ */
  /*   float deltaz = star_properties->feedback.yield_AGB.metallicity[*iz_high] - */
  /*                  star_properties->feedback.yield_AGB.metallicity[*iz_low]; */

  /*   if (deltaz > 0) */
  /*     *dz /= deltaz; */
  /*   else */
  /*     dz = 0; */
  /* } else { */
  /*   *iz_low = 0; */
  /*   *iz_high = 0; */
  /*   *dz = 0; */
  /* } */
}

/**
 * @brief determine which metallicity bin star belongs to for SNII, compute bin
 * indices and offsets
 *
 * @param iz_low Pointer to index of metallicity bin to which the star belongs
 * (to be calculated in this function)
 * @param iz_high Pointer to index of metallicity bin to above the star's
 * metallicity (to be calculated in this function)
 * @param dz metallicity bin offset
 * @param log10_metallicity log10 of star metallicity
 * @param star_properties stars_props data structure
 */
inline static void determine_bin_yield_SNII(
    int* iz_low, int* iz_high, float* dz, float log10_metallicity,
    const struct feedback_props* feedback_props) {

  // MATTHIEU
  
  /* if (log10_metallicity > log10_min_metallicity) { */
  /*   /\* Find metallicity bin which contains the star's metallicity *\/ */
  /*   int j; */
  /*   for (j = 0; j < star_properties->feedback.SNII_n_z - 1 && */
  /*               log10_metallicity > */
  /*                   star_properties->feedback.yield_SNII.metallicity[j + 1]; */
  /*        j++) */
  /*     ; */
  /*   *iz_low = j; */
  /*   *iz_high = *iz_low + 1; */

  /*   /\* Compute offset *\/ */
  /*   if (log10_metallicity >= */
  /*           star_properties->feedback.yield_SNII.metallicity[0] && */
  /*       log10_metallicity <= */
  /*           star_properties->feedback.yield_SNII */
  /*               .metallicity[star_properties->feedback.SNII_n_z - 1]) */
  /*     *dz = log10_metallicity - */
  /*           star_properties->feedback.yield_SNII.metallicity[*iz_low]; */
  /*   else */
  /*     *dz = 0; */

  /*   /\* Normalize offset *\/ */
  /*   float deltaz = star_properties->feedback.yield_SNII.metallicity[*iz_high] - */
  /*                  star_properties->feedback.yield_SNII.metallicity[*iz_low]; */

  /*   if (deltaz > 0) */
  /*     *dz = *dz / deltaz; */
  /*   else */
  /*     dz = 0; */
  /* } else { */
  /*   *iz_low = 0; */
  /*   *iz_high = 0; */
  /*   *dz = 0; */
  /* } */
}


/**
 * @brief compute enrichment and feedback due to SNIa. To do this compute the
 * number of SNIa that occur during the timestep, multiply by constants read
 * from tables.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stars star properties data structure
 * @param sp spart we are computing feedback from
 * @param star_age_Gyr age of star in Gyr
 * @param dt_Gyr timestep dt in Gyr
 */
inline static void evolve_SNIa(float log10_min_mass, float log10_max_mass,
                               const struct feedback_props* feedback_props,
                               struct spart* restrict sp, float star_age_Gyr,
                               float dt_Gyr) {

  /* Check if we're outside the mass range for SNIa */
  if (log10_min_mass >= feedback_props->log10_SNIa_max_mass_msun) return;

  /* If the max mass is outside the mass range update it to be the maximum and
   * use updated values for the star's age and timestep in this function */
  if (log10_max_mass > feedback_props->log10_SNIa_max_mass_msun) {
    log10_max_mass = feedback_props->log10_SNIa_max_mass_msun;
    float lifetime_Gyr =
        lifetime_in_Gyr(exp(M_LN10 * feedback_props->log10_SNIa_max_mass_msun),
                        sp->chemistry_data.metal_mass_fraction_total, feedback_props);
    dt_Gyr = star_age_Gyr + dt_Gyr - lifetime_Gyr;
    star_age_Gyr = lifetime_Gyr;
  }

  /* compute the number of SNIa */
  /* Efolding (Forster 2006) */
  float num_SNIa_per_msun =
      feedback_props->SNIa_efficiency *
      (exp(-star_age_Gyr / feedback_props->SNIa_timescale_Gyr) -
       exp(-(star_age_Gyr + dt_Gyr) / feedback_props->SNIa_timescale_Gyr)) *
    sp->mass_init;

  sp->feedback_data.to_distribute.num_SNIa =
      num_SNIa_per_msun / feedback_props->const_solar_mass;

  /* compute mass fractions of each metal */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.to_distribute.metal_mass[i] +=
        num_SNIa_per_msun * feedback_props->yield_SNIa_IMF_resampled[i];
  }

  /* Update the metallicity of the material released */
  sp->feedback_data.to_distribute.metal_mass_from_SNIa +=
      num_SNIa_per_msun * feedback_props->yield_SNIa_total_metals_IMF_resampled;

  /* Update the metal mass produced */
  sp->feedback_data.to_distribute.total_metal_mass +=
      num_SNIa_per_msun * feedback_props->yield_SNIa_total_metals_IMF_resampled;

  /* Compute the mass produced by SNIa */
  sp->feedback_data.to_distribute.mass_from_SNIa +=
      num_SNIa_per_msun * feedback_props->yield_SNIa_total_metals_IMF_resampled;

  /* Compute the iron mass produced */
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa +=
      num_SNIa_per_msun *
      feedback_props->yield_SNIa_IMF_resampled[chemistry_element_Fe];
}

/**
 * @brief compute enrichment and feedback due to SNII. To do this, integrate the
 * IMF weighted by the yields read from tables for each of the quantities of
 * interest.
 *
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param stars star properties data structure
 * @param sp spart we are computing feedback from
 */
inline static void evolve_SNII(float log10_min_mass, float log10_max_mass,
                               float* stellar_yields,
                               const struct feedback_props* restrict feedback_props,
                               struct spart* restrict sp) {
  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* If mass at beginning of step is less than tabulated lower bound for IMF,
   * limit it.*/
  if (log10_min_mass < feedback_props->log10_SNII_min_mass_msun)
    log10_min_mass = feedback_props->log10_SNII_min_mass_msun;

  /* If mass at end of step is greater than tabulated upper bound for IMF, limit
   * it.*/
  if (log10_max_mass > feedback_props->log10_SNII_max_mass_msun)
    log10_max_mass = feedback_props->log10_SNII_max_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, feedback_props);

  /* Integrate IMF to determine number of SNII */
  sp->feedback_data.to_distribute.num_SNII =
      integrate_imf(log10_min_mass, log10_max_mass, 0, stellar_yields, feedback_props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d, high_index_2d;
  float dz = 0.;
  determine_bin_yield_SNII(&iz_low, &iz_high, &dz,
                           log10(sp->chemistry_data.metal_mass_fraction_total),
                           feedback_props);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = feedback_row_major_index_3d(
          iz_low, elem, mass_bin_index, feedback_props->SNII_n_z,
          chemistry_element_count, feedback_props->n_imf_mass_bins);
      high_index_3d = feedback_row_major_index_3d(
          iz_high, elem, mass_bin_index, feedback_props->SNII_n_z,
          chemistry_element_count, feedback_props->n_imf_mass_bins);
      low_index_2d = feedback_row_major_index_2d(
          iz_low, mass_bin_index, feedback_props->SNII_n_z,
          feedback_props->n_imf_mass_bins);
      high_index_2d = feedback_row_major_index_2d(
          iz_high, mass_bin_index, feedback_props->SNII_n_z,
          feedback_props->n_imf_mass_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (feedback_props->yield_SNII.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   feedback_props->yield_SNII
                       .ejecta_IMF_resampled[low_index_2d]) +
          dz * (feedback_props->yield_SNII.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    feedback_props->yield_SNII
                        .ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] =
        integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               feedback_props->SNII_n_z,
                                               feedback_props->n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, feedback_props->SNII_n_z,
        feedback_props->n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (feedback_props->yield_SNII
                        .total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        feedback_props->yield_SNII
                            .ejecta_IMF_resampled[low_index_2d]) +
        dz * (feedback_props->yield_SNII
                  .total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  feedback_props->yield_SNII
                      .ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star*/
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               feedback_props->SNII_n_z,
                                               feedback_props->n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, feedback_props->SNII_n_z,
        feedback_props->n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            feedback_props->yield_SNII.ejecta_IMF_resampled[low_index_2d] +
        dz * feedback_props->yield_SNII.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);

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
      sp->feedback_data.to_distribute.metal_mass[i] += metal_mass_released[i] * norm_factor;
    }
    for (int i = 0; i < chemistry_element_count; i++) {
      sp->feedback_data.to_distribute.mass_from_SNII += sp->feedback_data.to_distribute.metal_mass[i];
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
 * @param log10_min_mass log10 mass at the end of step
 * @param log10_max_mass log10 mass at the beginning of step
 * @param stellar_yields array to store calculated yields for passing to
 * integrate_imf
 * @param feedback_props star properties data structure
 * @param sp spart we are computing feedback from
 */
inline static void evolve_AGB(float log10_min_mass, float log10_max_mass,
                              float* stellar_yields,
                              const struct feedback_props* restrict feedback_props,
                              struct spart* restrict sp) {
  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* If mass at end of step is greater than tabulated lower bound for IMF, limit
   * it.*/
  if (log10_max_mass > feedback_props->log10_SNII_min_mass_msun)
    log10_max_mass = feedback_props->log10_SNII_min_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, feedback_props);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low = 0, iz_high = 0, low_index_3d, high_index_3d, low_index_2d, high_index_2d;
  float dz = 0.f;
  determine_bin_yield_AGB(&iz_low, &iz_high, &dz,
                          log10(sp->chemistry_data.metal_mass_fraction_total),
                          feedback_props);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = feedback_row_major_index_3d(
          iz_low, elem, mass_bin_index, feedback_props->AGB_n_z,
          chemistry_element_count, feedback_props->n_imf_mass_bins);
      high_index_3d = feedback_row_major_index_3d(
          iz_high, elem, mass_bin_index, feedback_props->AGB_n_z,
          chemistry_element_count, feedback_props->n_imf_mass_bins);
      low_index_2d = feedback_row_major_index_2d(
          iz_low, mass_bin_index, feedback_props->AGB_n_z,
          feedback_props->n_imf_mass_bins);
      high_index_2d = feedback_row_major_index_2d(
          iz_high, mass_bin_index, feedback_props->AGB_n_z,
          feedback_props->n_imf_mass_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (feedback_props->yield_AGB.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   feedback_props->yield_AGB
                       .ejecta_IMF_resampled[low_index_2d]) +
          dz * (feedback_props->yield_AGB.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    feedback_props->yield_AGB
                        .ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] =
        integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               feedback_props->AGB_n_z,
                                               feedback_props->n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, feedback_props->AGB_n_z,
        feedback_props->n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            (feedback_props->yield_AGB
                 .total_metals_IMF_resampled[low_index_2d] +
             sp->chemistry_data.metal_mass_fraction_total *
                 feedback_props->yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
        dz *
            (feedback_props->yield_AGB
                 .total_metals_IMF_resampled[high_index_2d] +
             sp->chemistry_data.metal_mass_fraction_total *
                 feedback_props->yield_AGB.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);

  /* yield normalization */
  float mass_ejected, mass_released;

  /* zero all negative values */
  for (int i = 0; i < chemistry_element_count; i++)
    metal_mass_released[i] = max(metal_mass_released[i], 0.f);

  metal_mass_released_total = max(metal_mass_released_total, 0.f);

  /* compute the total metal mass ejected from the star */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               feedback_props->AGB_n_z,
                                               feedback_props->n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, feedback_props->AGB_n_z,
        feedback_props->n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            feedback_props->yield_AGB.ejecta_IMF_resampled[low_index_2d] +
        dz * feedback_props->yield_AGB.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, feedback_props);

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
      sp->feedback_data.to_distribute.metal_mass[i] += metal_mass_released[i] * norm_factor;
      sp->feedback_data.to_distribute.mass_from_AGB += metal_mass_released[i] * norm_factor;
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
 * @param star_properties feedback_props_props data structure
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 */
void compute_stellar_evolution(
			       const struct feedback_props* feedback_props,
    struct spart* sp, const struct unit_system* us, float age,
    double dt) {

  /* Allocate temporary array for calculating imf weights */
  float* stellar_yields;
  stellar_yields =
      malloc(feedback_props->n_imf_mass_bins * sizeof(float));

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_in_cgs = 3.154e16;
  double dt_Gyr =
      dt * units_cgs_conversion_factor(us, UNIT_CONV_TIME) / Gyr_in_cgs;
  double star_age_Gyr =
      age * units_cgs_conversion_factor(us, UNIT_CONV_TIME) / Gyr_in_cgs;

  /* calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  double log10_max_dying_mass_msun = log10(dying_mass_msun(
      star_age_Gyr, sp->chemistry_data.metal_mass_fraction_total,
      feedback_props));
  double log10_min_dying_mass_msun = log10(dying_mass_msun(
      star_age_Gyr + dt_Gyr, sp->chemistry_data.metal_mass_fraction_total,
      feedback_props));

  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (log10_min_dying_mass_msun > log10_max_dying_mass_msun)
    error("min dying mass is greater than max dying mass");

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_msun. Return without doing any
   * feedback. */
  if (log10_min_dying_mass_msun == log10_max_dying_mass_msun) return;

  /* Evolve SNIa, SNII, AGB */
  evolve_SNIa(log10_min_dying_mass_msun, log10_max_dying_mass_msun,
              feedback_props, sp, star_age_Gyr, dt_Gyr);
  evolve_SNII(log10_min_dying_mass_msun, log10_max_dying_mass_msun,
              stellar_yields, feedback_props, sp);
  evolve_AGB(log10_min_dying_mass_msun, log10_max_dying_mass_msun,
             stellar_yields, feedback_props, sp);

  /* Compute the mass to distribute */
  sp->feedback_data.to_distribute.mass = sp->feedback_data.to_distribute.total_metal_mass +
                           sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
                           sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Clean up */
  free(stellar_yields);
}

/**
 * @brief Compute number of SN that should go off given the age of the spart
 *
 * @param sp spart we're evolving
 * @param stars_properties stars_props data structure
 * @param age age of star at the beginning of the step
 * @param dt length of step
 */
float compute_SNe(struct spart* sp,
                                const struct feedback_props* feedback_props,
                                float age, double dt) {
  
  if (age <= feedback_props->SNII_wind_delay &&
      age + dt > feedback_props->SNII_wind_delay) {

    return feedback_props->num_SNII_per_msun * sp->mass_init / 
           feedback_props->const_solar_mass;
  } else {
    return 0;
  }
}


/**
 * @brief Initializes constants related to stellar evolution, initializes imf,
 * reads and processes yield tables
 *
 * @param params swift_params parameters structure
 * @param stars stars_props data structure
 */
inline static void stars_evolve_init(struct swift_params* params,
                                     struct feedback_props* feedback_props) {

  /* Set number of elements found in yield tables */
  feedback_props->SNIa_n_elements = 42;
  feedback_props->SNII_n_mass = 11;
  feedback_props->SNII_n_elements = 11;
  feedback_props->SNII_n_z = 5;
  feedback_props->AGB_n_mass = 23;
  feedback_props->AGB_n_elements = 11;
  feedback_props->AGB_n_z = 3;
  feedback_props->lifetimes.n_mass = 30;
  feedback_props->lifetimes.n_z = 6;
  feedback_props->element_name_length = 15;

  /* Set bounds for imf  */
  feedback_props->log10_SNII_min_mass_msun = 0.77815125f;  // log10(6).
  feedback_props->log10_SNII_max_mass_msun = 2.f;          // log10(100).
  feedback_props->log10_SNIa_max_mass_msun = 0.90308999f;  // log10(8).

  /* Yield table filepath  */
  parser_get_param_string(params, "EAGLEFeedback:filename",
                          feedback_props->yield_table_path);
  parser_get_param_string(params, "EAGLEFeedback:imf_model",
                          feedback_props->IMF_Model);

  /* Initialise IMF */
  init_imf(feedback_props);

  /* Allocate yield tables  */
  allocate_yield_tables(feedback_props);

  /* Set factors for each element adjusting SNII yield */
  feedback_props->typeII_factor[0] = 1.f;
  feedback_props->typeII_factor[1] = 1.f;
  feedback_props->typeII_factor[2] = 0.5f;
  feedback_props->typeII_factor[3] = 1.f;
  feedback_props->typeII_factor[4] = 1.f;
  feedback_props->typeII_factor[5] = 1.f;
  feedback_props->typeII_factor[6] = 2.f;
  feedback_props->typeII_factor[7] = 1.f;
  feedback_props->typeII_factor[8] = 0.5f;

  /* Read the tables  */
  read_yield_tables(feedback_props);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (feedback_props->log10_imf_max_mass_msun -
       feedback_props->log10_imf_min_mass_msun) /
      (feedback_props->n_imf_mass_bins - 1);
  for (int i = 0; i < feedback_props->n_imf_mass_bins; i++)
    feedback_props->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + feedback_props->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(feedback_props);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(feedback_props);

  /* Calculate number of type II SN per solar mass
   * Note: since we are integrating the IMF without weighting it by the yields
   * pass NULL pointer for stellar_yields array */
  feedback_props->num_SNII_per_msun =
      integrate_imf(feedback_props->log10_SNII_min_mass_msun,
                    feedback_props->log10_SNII_max_mass_msun, 0,
                    /*(stellar_yields=)*/ NULL, feedback_props);

  message("initialized stellar feedback");
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * By default, takes the values provided by the hydro.
 *
 * @param sp The #feedback_properties.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The already read-in properties of the hydro scheme.
 */
void feedbakc_props_init(struct feedback_props *fp,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct hydro_props *p,
                                       const struct cosmology *cosmo) {

  /* Read SNIa timscale */
  fp->SNIa_timescale_Gyr =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_timescale_Gyr");

  /* Read the efficiency of producing SNIa */
  fp->SNIa_efficiency =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_efficiency");

  /* Are we doing continuous heating? */
  fp->continuous_heating =
      parser_get_param_int(params, "EAGLEFeedback:continuous_heating_switch");

  /* Set the delay time before SNII occur */
  const float Gyr_in_cgs = 1e9 * 365 * 24 * 3600;
  fp->SNII_wind_delay =
      parser_get_param_float(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
      Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Read the temperature change to use in stochastic heating */
  fp->SNe_deltaT_desired =
      parser_get_param_float(params, "EAGLEFeedback:SNe_heating_temperature_K");
  fp->SNe_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Set ejecta thermal energy */
  const float ejecta_velocity =
      1.0e6 / units_cgs_conversion_factor(
                  us, UNIT_CONV_SPEED);  // EAGLE parameter is 10 km/s
  fp->ejecta_specific_thermal_energy = 0.5 * ejecta_velocity * ejecta_velocity;

  /* Energy released by supernova */
  fp->total_energy_SNe =
      1.0e51 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Calculate temperature to internal energy conversion factor */
  fp->temp_to_u_factor =
      phys_const->const_boltzmann_k /
      (p->mu_ionised * (hydro_gamma_minus_one)*phys_const->const_proton_mass);

  /* Read birth time to set all stars in ICs to (defaults to -1 to indicate star
   * present in ICs) */
  fp->spart_first_init_birth_time = parser_get_opt_param_float(
      params, "EAGLEFeedback:birth_time_override", -1);

  /* Copy over solar mass */
  fp->const_solar_mass = phys_const->const_solar_mass;
}


