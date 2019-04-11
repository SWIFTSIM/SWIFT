/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_STARS_H
#define SWIFT_EAGLE_STARS_H

#include <float.h>
#include "imf.h"
#include "minmax.h"
#include "yield_tables.h"

/**
 * @brief Computes the gravity time-step of a given star particle.
 *
 * @param sp Pointer to the s-particle data.
 */
__attribute__((always_inline)) INLINE static float stars_compute_timestep(
    const struct spart* const sp) {

  return FLT_MAX;
}

/**
 * @brief Prepares a s-particle for its interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_init_spart(
    struct spart* sp) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    sp->ids_ngbs_density[i] = -1;
  sp->num_ngb_density = 0;
#endif

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
  sp->rho_gas = 0.f;

  sp->density_weighted_frac_normalisation_inv = 0.f;
  sp->ngb_mass = 0.f;
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp, const struct stars_props *stars_properties) {

  sp->time_bin = 0;
  sp->birth_density = -1.f;
  sp->birth_time = stars_properties->feedback.spart_first_init_birth_time;

  stars_init_spart(sp);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void stars_predict_extra(
    struct spart* restrict sp, float dt_drift) {

  // MATTHIEU
  /* const float h_inv = 1.f / sp->h; */

  /* /\* Predict smoothing length *\/ */
  /* const float w1 = sp->feedback.h_dt * h_inv * dt_drift; */
  /* if (fabsf(w1) < 0.2f) */
  /*   sp->h *= approx_expf(w1); /\* 4th order expansion of exp(w) *\/ */
  /* else */
  /*   sp->h *= expf(w1); */
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void stars_reset_predicted_values(
    struct spart* restrict sp) {}

/**
 * @brief Finishes the calculation of (non-gravity) forces acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_feedback(
    struct spart* sp) {

  sp->feedback.h_dt *= sp->h * hydro_dimension_inv;
}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void stars_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_density(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->rho_gas *= h_inv_dim;
  sp->density.wcount *= h_inv_dim;
  sp->density.wcount_dh *= h_inv_dim_plus_one;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_spart_has_no_neighbours(
    struct spart* restrict sp, const struct cosmology* cosmo) {

  /* Re-set problematic values */
  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
  sp->rho_gas = 0.f;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_acceleration(
    struct spart* restrict p) {
#ifdef DEBUG_INTERACTIONS_STARS
  p->num_ngb_force = 0;
#endif
}

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
    const struct stars_props* restrict star_properties) {

  if (log10_metallicity > log10_min_metallicity) {
    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < star_properties->feedback.AGB_n_z - 1 &&
                log10_metallicity >
                    star_properties->feedback.yield_AGB.metallicity[j + 1];
         j++)
      ;
    *iz_low = j;
    *iz_high = *iz_low + 1;

    /* Compute offset */
    if (log10_metallicity >=
            star_properties->feedback.yield_AGB.metallicity[0] &&
        log10_metallicity <=
            star_properties->feedback.yield_AGB
                .metallicity[star_properties->feedback.AGB_n_z - 1])
      *dz = log10_metallicity -
            star_properties->feedback.yield_AGB.metallicity[*iz_low];
    else
      *dz = 0;

    /* Normalize offset */
    float deltaz = star_properties->feedback.yield_AGB.metallicity[*iz_high] -
                   star_properties->feedback.yield_AGB.metallicity[*iz_low];

    if (deltaz > 0)
      *dz /= deltaz;
    else
      dz = 0;
  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0;
  }
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
    const struct stars_props* restrict star_properties) {

  if (log10_metallicity > log10_min_metallicity) {
    /* Find metallicity bin which contains the star's metallicity */
    int j;
    for (j = 0; j < star_properties->feedback.SNII_n_z - 1 &&
                log10_metallicity >
                    star_properties->feedback.yield_SNII.metallicity[j + 1];
         j++)
      ;
    *iz_low = j;
    *iz_high = *iz_low + 1;

    /* Compute offset */
    if (log10_metallicity >=
            star_properties->feedback.yield_SNII.metallicity[0] &&
        log10_metallicity <=
            star_properties->feedback.yield_SNII
                .metallicity[star_properties->feedback.SNII_n_z - 1])
      *dz = log10_metallicity -
            star_properties->feedback.yield_SNII.metallicity[*iz_low];
    else
      *dz = 0;

    /* Normalize offset */
    float deltaz = star_properties->feedback.yield_SNII.metallicity[*iz_high] -
                   star_properties->feedback.yield_SNII.metallicity[*iz_low];

    if (deltaz > 0)
      *dz = *dz / deltaz;
    else
      dz = 0;
  } else {
    *iz_low = 0;
    *iz_high = 0;
    *dz = 0;
  }
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
                               const struct stars_props* restrict stars,
                               struct spart* restrict sp, float star_age_Gyr,
                               float dt_Gyr) {

  /* Check if we're outside the mass range for SNIa */
  if (log10_min_mass >= stars->feedback.log10_SNIa_max_mass_msun) return;

  /* If the max mass is outside the mass range update it to be the maximum and
   * use updated values for the star's age and timestep in this function */
  if (log10_max_mass > stars->feedback.log10_SNIa_max_mass_msun) {
    log10_max_mass = stars->feedback.log10_SNIa_max_mass_msun;
    float lifetime_Gyr =
        lifetime_in_Gyr(exp(M_LN10 * stars->feedback.log10_SNIa_max_mass_msun),
                        sp->chemistry_data.metal_mass_fraction_total, stars);
    dt_Gyr = star_age_Gyr + dt_Gyr - lifetime_Gyr;
    star_age_Gyr = lifetime_Gyr;
  }

  /* compute the number of SNIa */
  /* Efolding (Forster 2006) */
  float num_SNIa_per_msun =
      stars->feedback.SNIa_efficiency *
      (exp(-star_age_Gyr / stars->feedback.SNIa_timescale_Gyr) -
       exp(-(star_age_Gyr + dt_Gyr) / stars->feedback.SNIa_timescale_Gyr)) *
      sp->mass_init;

  sp->to_distribute.num_SNIa =
      num_SNIa_per_msun / stars->feedback.const_solar_mass;

  /* compute mass fractions of each metal */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->to_distribute.metal_mass[i] +=
        num_SNIa_per_msun * stars->feedback.yield_SNIa_IMF_resampled[i];
  }

  /* Update the metallicity of the material released */
  sp->to_distribute.metal_mass_from_SNIa +=
      num_SNIa_per_msun * stars->feedback.yield_SNIa_total_metals_IMF_resampled;

  /* Update the metal mass produced */
  sp->to_distribute.total_metal_mass +=
      num_SNIa_per_msun * stars->feedback.yield_SNIa_total_metals_IMF_resampled;

  /* Compute the mass produced by SNIa */
  sp->to_distribute.mass_from_SNIa +=
      num_SNIa_per_msun * stars->feedback.yield_SNIa_total_metals_IMF_resampled;

  /* Compute the iron mass produced */
  sp->to_distribute.Fe_mass_from_SNIa +=
      num_SNIa_per_msun *
      stars->feedback.yield_SNIa_IMF_resampled[chemistry_element_Fe];
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
                               const struct stars_props* restrict stars,
                               struct spart* restrict sp) {
  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* If mass at beginning of step is less than tabulated lower bound for IMF,
   * limit it.*/
  if (log10_min_mass < stars->feedback.log10_SNII_min_mass_msun)
    log10_min_mass = stars->feedback.log10_SNII_min_mass_msun;

  /* If mass at end of step is greater than tabulated upper bound for IMF, limit
   * it.*/
  if (log10_max_mass > stars->feedback.log10_SNII_max_mass_msun)
    log10_max_mass = stars->feedback.log10_SNII_max_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, stars);

  /* Integrate IMF to determine number of SNII */
  sp->to_distribute.num_SNII =
      integrate_imf(log10_min_mass, log10_max_mass, 0, stellar_yields, stars);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low, iz_high, low_index_3d, high_index_3d, low_index_2d, high_index_2d;
  float dz;
  determine_bin_yield_SNII(&iz_low, &iz_high, &dz,
                           log10(sp->chemistry_data.metal_mass_fraction_total),
                           stars);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = feedback_row_major_index_3d(
          iz_low, elem, mass_bin_index, stars->feedback.SNII_n_z,
          chemistry_element_count, stars->feedback.n_imf_mass_bins);
      high_index_3d = feedback_row_major_index_3d(
          iz_high, elem, mass_bin_index, stars->feedback.SNII_n_z,
          chemistry_element_count, stars->feedback.n_imf_mass_bins);
      low_index_2d = feedback_row_major_index_2d(
          iz_low, mass_bin_index, stars->feedback.SNII_n_z,
          stars->feedback.n_imf_mass_bins);
      high_index_2d = feedback_row_major_index_2d(
          iz_high, mass_bin_index, stars->feedback.SNII_n_z,
          stars->feedback.n_imf_mass_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (stars->feedback.yield_SNII.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   stars->feedback.yield_SNII
                       .ejecta_IMF_resampled[low_index_2d]) +
          dz * (stars->feedback.yield_SNII.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    stars->feedback.yield_SNII
                        .ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] =
        integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               stars->feedback.SNII_n_z,
                                               stars->feedback.n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, stars->feedback.SNII_n_z,
        stars->feedback.n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) * (stars->feedback.yield_SNII
                        .total_metals_IMF_resampled[low_index_2d] +
                    sp->chemistry_data.metal_mass_fraction_total *
                        stars->feedback.yield_SNII
                            .ejecta_IMF_resampled[low_index_2d]) +
        dz * (stars->feedback.yield_SNII
                  .total_metals_IMF_resampled[high_index_2d] +
              sp->chemistry_data.metal_mass_fraction_total *
                  stars->feedback.yield_SNII
                      .ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);

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
                                               stars->feedback.SNII_n_z,
                                               stars->feedback.n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, stars->feedback.SNII_n_z,
        stars->feedback.n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            stars->feedback.yield_SNII.ejecta_IMF_resampled[low_index_2d] +
        dz * stars->feedback.yield_SNII.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);

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
      sp->to_distribute.metal_mass[i] += metal_mass_released[i] * norm_factor;
    }
    for (int i = 0; i < chemistry_element_count; i++) {
      sp->to_distribute.mass_from_SNII += sp->to_distribute.metal_mass[i];
    }
    sp->to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->to_distribute.metal_mass_from_SNII +=
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
 * @param stars star properties data structure
 * @param sp spart we are computing feedback from
 */
inline static void evolve_AGB(float log10_min_mass, float log10_max_mass,
                              float* stellar_yields,
                              const struct stars_props* restrict stars,
                              struct spart* restrict sp) {
  int low_imf_mass_bin_index, high_imf_mass_bin_index, mass_bin_index;

  /* If mass at end of step is greater than tabulated lower bound for IMF, limit
   * it.*/
  if (log10_max_mass > stars->feedback.log10_SNII_min_mass_msun)
    log10_max_mass = stars->feedback.log10_SNII_min_mass_msun;

  /* Don't do anything if the stellar mass hasn't decreased by the end of the
   * step */
  if (log10_min_mass >= log10_max_mass) return;

  /* determine which IMF mass bins contribute to the integral */
  determine_imf_bins(log10_min_mass, log10_max_mass, &low_imf_mass_bin_index,
                     &high_imf_mass_bin_index, stars);

  /* determine which metallicity bin and offset this star belongs to */
  int iz_low, iz_high, low_index_3d, high_index_3d, low_index_2d, high_index_2d;
  float dz;
  determine_bin_yield_AGB(&iz_low, &iz_high, &dz,
                          log10(sp->chemistry_data.metal_mass_fraction_total),
                          stars);

  /* compute metals produced */
  float metal_mass_released[chemistry_element_count], metal_mass_released_total;
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    for (mass_bin_index = low_imf_mass_bin_index;
         mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
      low_index_3d = feedback_row_major_index_3d(
          iz_low, elem, mass_bin_index, stars->feedback.AGB_n_z,
          chemistry_element_count, stars->feedback.n_imf_mass_bins);
      high_index_3d = feedback_row_major_index_3d(
          iz_high, elem, mass_bin_index, stars->feedback.AGB_n_z,
          chemistry_element_count, stars->feedback.n_imf_mass_bins);
      low_index_2d = feedback_row_major_index_2d(
          iz_low, mass_bin_index, stars->feedback.AGB_n_z,
          stars->feedback.n_imf_mass_bins);
      high_index_2d = feedback_row_major_index_2d(
          iz_high, mass_bin_index, stars->feedback.AGB_n_z,
          stars->feedback.n_imf_mass_bins);
      stellar_yields[mass_bin_index] =
          (1 - dz) *
              (stars->feedback.yield_AGB.yield_IMF_resampled[low_index_3d] +
               sp->chemistry_data.metal_mass_fraction[elem] *
                   stars->feedback.yield_AGB
                       .ejecta_IMF_resampled[low_index_2d]) +
          dz * (stars->feedback.yield_AGB.yield_IMF_resampled[high_index_3d] +
                sp->chemistry_data.metal_mass_fraction[elem] *
                    stars->feedback.yield_AGB
                        .ejecta_IMF_resampled[high_index_2d]);
    }

    metal_mass_released[elem] =
        integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);
  }

  /* Compute mass produced */
  for (mass_bin_index = low_imf_mass_bin_index;
       mass_bin_index < high_imf_mass_bin_index + 1; mass_bin_index++) {
    low_index_2d = feedback_row_major_index_2d(iz_low, mass_bin_index,
                                               stars->feedback.AGB_n_z,
                                               stars->feedback.n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, stars->feedback.AGB_n_z,
        stars->feedback.n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            (stars->feedback.yield_AGB
                 .total_metals_IMF_resampled[low_index_2d] +
             sp->chemistry_data.metal_mass_fraction_total *
                 stars->feedback.yield_AGB.ejecta_IMF_resampled[low_index_2d]) +
        dz *
            (stars->feedback.yield_AGB
                 .total_metals_IMF_resampled[high_index_2d] +
             sp->chemistry_data.metal_mass_fraction_total *
                 stars->feedback.yield_AGB.ejecta_IMF_resampled[high_index_2d]);
  }

  metal_mass_released_total =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);

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
                                               stars->feedback.AGB_n_z,
                                               stars->feedback.n_imf_mass_bins);
    high_index_2d = feedback_row_major_index_2d(
        iz_high, mass_bin_index, stars->feedback.AGB_n_z,
        stars->feedback.n_imf_mass_bins);
    stellar_yields[mass_bin_index] =
        (1 - dz) *
            stars->feedback.yield_AGB.ejecta_IMF_resampled[low_index_2d] +
        dz * stars->feedback.yield_AGB.ejecta_IMF_resampled[high_index_2d];
  }

  mass_ejected =
      integrate_imf(log10_min_mass, log10_max_mass, 2, stellar_yields, stars);

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
      sp->to_distribute.metal_mass[i] += metal_mass_released[i] * norm_factor;
      sp->to_distribute.mass_from_AGB += metal_mass_released[i] * norm_factor;
    }
    sp->to_distribute.total_metal_mass +=
        metal_mass_released_total * norm_factor;
    sp->to_distribute.metal_mass_from_AGB +=
        metal_mass_released_total * norm_factor;
  } else {
    error("wrong normalization!!!! mass_released = %e\n", mass_released);
  }
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param star_properties stars_props data structure
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 */
inline static void compute_stellar_evolution(
    const struct stars_props* restrict star_properties,
    struct spart* restrict sp, const struct unit_system* us, float age,
    double dt) {

  /* Allocate temporary array for calculating imf weights */
  float* stellar_yields;
  stellar_yields =
      malloc(star_properties->feedback.n_imf_mass_bins * sizeof(float));

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
      star_properties));
  double log10_min_dying_mass_msun = log10(dying_mass_msun(
      star_age_Gyr + dt_Gyr, sp->chemistry_data.metal_mass_fraction_total,
      star_properties));

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
              star_properties, sp, star_age_Gyr, dt_Gyr);
  evolve_SNII(log10_min_dying_mass_msun, log10_max_dying_mass_msun,
              stellar_yields, star_properties, sp);
  evolve_AGB(log10_min_dying_mass_msun, log10_max_dying_mass_msun,
             stellar_yields, star_properties, sp);

  /* Compute the mass to distribute */
  sp->to_distribute.mass = sp->to_distribute.total_metal_mass +
                           sp->to_distribute.metal_mass[chemistry_element_H] +
                           sp->to_distribute.metal_mass[chemistry_element_He];

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
inline static float compute_SNe(struct spart* sp,
                                const struct stars_props* stars_properties,
                                float age, double dt) {
  if (age <= stars_properties->feedback.SNII_wind_delay &&
      age + dt > stars_properties->feedback.SNII_wind_delay) {
    return stars_properties->feedback.num_SNII_per_msun * sp->mass_init /
           stars_properties->feedback.const_solar_mass;
  } else {
    return 0;
  }
}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 * @param stars_properties The #stars_props
 */
__attribute__((always_inline)) INLINE static void stars_evolve_spart(
    struct spart* restrict sp, const struct stars_props* stars_properties,
    const struct cosmology* cosmo, const struct unit_system* us,
    double star_age, double dt) {

  /* Zero the number of SN and amount of mass that is distributed */
  sp->to_distribute.num_SNIa = 0;
  sp->to_distribute.num_SNII = 0;
  sp->to_distribute.mass = 0;

  /* Zero the enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++)
    sp->to_distribute.metal_mass[i] = 0;
  sp->to_distribute.total_metal_mass = 0;
  sp->to_distribute.mass_from_AGB = 0;
  sp->to_distribute.metal_mass_from_AGB = 0;
  sp->to_distribute.mass_from_SNII = 0;
  sp->to_distribute.metal_mass_from_SNII = 0;
  sp->to_distribute.mass_from_SNIa = 0;
  sp->to_distribute.metal_mass_from_SNIa = 0;
  sp->to_distribute.Fe_mass_from_SNIa = 0;

  /* Compute amount of enrichment and feedback that needs to be done in this
   * step */
  compute_stellar_evolution(stars_properties, sp, us, star_age, dt);

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= sp->to_distribute.mass;

  /* Compute the number of type II SNe that went off */
  sp->to_distribute.num_SNe = compute_SNe(sp, stars_properties, star_age, dt);

  /* Compute energy change due to thermal and kinetic energy of ejecta */
  sp->to_distribute.d_energy = sp->to_distribute.mass *
      (stars_properties->feedback.ejecta_specific_thermal_energy +
       0.5 * (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
           cosmo->a2_inv);

  /* Compute probability of heating neighbouring particles */
  if (dt > 0 && sp->ngb_mass > 0)  sp->to_distribute.heating_probability = 
      stars_properties->feedback.total_energy_SNe *
      sp->to_distribute.num_SNe /
      (stars_properties->feedback.temp_to_u_factor *
      stars_properties->feedback.SNe_deltaT_desired * sp->ngb_mass);
}

/**
 * @brief Initializes constants related to stellar evolution, initializes imf,
 * reads and processes yield tables
 *
 * @param params swift_params parameters structure
 * @param stars stars_props data structure
 */
inline static void stars_evolve_init(struct swift_params* params,
                                     struct stars_props* restrict stars) {

  /* Set number of elements found in yield tables */
  stars->feedback.SNIa_n_elements = 42;
  stars->feedback.SNII_n_mass = 11;
  stars->feedback.SNII_n_elements = 11;
  stars->feedback.SNII_n_z = 5;
  stars->feedback.AGB_n_mass = 23;
  stars->feedback.AGB_n_elements = 11;
  stars->feedback.AGB_n_z = 3;
  stars->feedback.lifetimes.n_mass = 30;
  stars->feedback.lifetimes.n_z = 6;
  stars->feedback.element_name_length = 15;

  /* Set bounds for imf  */
  stars->feedback.log10_SNII_min_mass_msun = 0.77815125f;  // log10(6).
  stars->feedback.log10_SNII_max_mass_msun = 2.f;          // log10(100).
  stars->feedback.log10_SNIa_max_mass_msun = 0.90308999f;  // log10(8).

  /* Yield table filepath  */
  parser_get_param_string(params, "EAGLEFeedback:filename",
                          stars->feedback.yield_table_path);
  parser_get_param_string(params, "EAGLEFeedback:imf_model",
                          stars->feedback.IMF_Model);

  /* Initialise IMF */
  init_imf(stars);

  /* Allocate yield tables  */
  allocate_yield_tables(stars);

  /* Set factors for each element adjusting SNII yield */
  stars->feedback.typeII_factor[0] = 1.f;
  stars->feedback.typeII_factor[1] = 1.f;
  stars->feedback.typeII_factor[2] = 0.5f;
  stars->feedback.typeII_factor[3] = 1.f;
  stars->feedback.typeII_factor[4] = 1.f;
  stars->feedback.typeII_factor[5] = 1.f;
  stars->feedback.typeII_factor[6] = 2.f;
  stars->feedback.typeII_factor[7] = 1.f;
  stars->feedback.typeII_factor[8] = 0.5f;

  /* Read the tables  */
  read_yield_tables(stars);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (stars->feedback.log10_imf_max_mass_msun -
       stars->feedback.log10_imf_min_mass_msun) /
      (stars->feedback.n_imf_mass_bins - 1);
  for (int i = 0; i < stars->feedback.n_imf_mass_bins; i++)
    stars->feedback.yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + stars->feedback.log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(stars);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(stars);

  /* Calculate number of type II SN per solar mass
   * Note: since we are integrating the IMF without weighting it by the yields
   * pass NULL pointer for stellar_yields array */
  stars->feedback.num_SNII_per_msun =
      integrate_imf(stars->feedback.log10_SNII_min_mass_msun,
                    stars->feedback.log10_SNII_max_mass_msun, 0,
                    /*(stellar_yields=)*/ NULL, stars);

  message("initialized stellar feedback");
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_feedback(
    struct spart* restrict p) {

  /* Reset time derivative */
  p->feedback.h_dt = 0.f;

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    p->ids_ngbs_force[i] = -1;
  p->num_ngb_force = 0;
#endif
}

#endif /* SWIFT_EAGLE_STARS_H */
