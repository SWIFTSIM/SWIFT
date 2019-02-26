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
#ifndef SWIFT_CONST_STARS_H
#define SWIFT_CONST_STARS_H

#include <float.h>
#include "minmax.h"

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
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp) {

  sp->time_bin = 0;
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

  sp->omega_normalisation_inv = 0.f;
  sp->ngb_mass = 0.f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param p The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void stars_predict_extra(
    struct spart *restrict sp, float dt_drift) {

  const float h_inv = 1.f / sp->h;

  /* Predict smoothing length */
  const float w1 = sp->feedback.h_dt * h_inv * dt_drift;
  if (fabsf(w1) < 0.2f)
    sp->h *= approx_expf(w1); /* 4th order expansion of exp(w) */
  else
    sp->h *= expf(w1);

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
 * @brief Finishes the calculation of feedback acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_feedback(
    struct spart* sp) {}

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

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  sp->density.wcount = kernel_root * h_inv_dim;
  sp->density.wcount_dh = 0.f;
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
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 * @param stars_properties The #stars_props
 * @param dt Timestep over which the particle evolves.
 */
__attribute__((always_inline)) INLINE static void stars_evolve_spart(
    struct spart* restrict sp, const struct stars_props* stars_properties,
    const struct cosmology* cosmo, const struct unit_system* us, float current_time, double dt) {
  
  /* Proportion of quantities to be released each timestep */
  float feedback_factor = dt/stars_properties->feedback_timescale;

  /* Amount of mass to distribute in one step */
  sp->to_distribute.mass = sp->mass_init * feedback_factor;

  /* Set all enrichment quantities to constant values */
  for (int i = 0; i < chemistry_element_count; i++) sp->to_distribute.chemistry_data.metal_mass_fraction[i] = 1.f/chemistry_element_count;
  sp->to_distribute.chemistry_data.metal_mass_fraction_total = 1.f - 2.f/chemistry_element_count;
  sp->to_distribute.chemistry_data.mass_from_AGB = 1.0e-2 * sp->to_distribute.mass;
  sp->to_distribute.chemistry_data.metal_mass_fraction_from_AGB = 1.0e-2;
  sp->to_distribute.chemistry_data.mass_from_SNII = 1.0e-2 * sp->to_distribute.mass;
  sp->to_distribute.chemistry_data.metal_mass_fraction_from_SNII = 1.0e-2;
  sp->to_distribute.chemistry_data.mass_from_SNIa = 1.0e-2 * sp->to_distribute.mass;
  sp->to_distribute.chemistry_data.metal_mass_fraction_from_SNIa = 1.0e-2;
  sp->to_distribute.chemistry_data.iron_mass_fraction_from_SNIa = 1.0e-2;

  /* Set feedback to constant values */
  const float total_sn = sp->mass_init / stars_properties->const_solar_mass * stars_properties->sn_per_msun;

  /* Print total_sn and timescale to be read by test script for checking stochastic energy injection */
  if (dt == 0) {
    FILE *feedback_output = fopen("feedback_properties.dat","w");
    fprintf(feedback_output,"%.5e \n%.5e \n", total_sn, stars_properties->feedback_timescale);
    fclose(feedback_output);
  }

  /* Set number of supernovae to distribute */
  sp->to_distribute.num_SNIa = total_sn * feedback_factor;

  /* Set ejected thermal energy */
  sp->to_distribute.ejecta_specific_thermal_energy = 1.0e-3;

}

inline static void stars_evolve_init(struct swift_params *params, struct stars_props* restrict stars){}


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


#endif /* SWIFT_CONST_STARS_H */
