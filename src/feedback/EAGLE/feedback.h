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
#ifndef SWIFT_FEEDBACK_EAGLE_H
#define SWIFT_FEEDBACK_EAGLE_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const float age,
                               const float dt);

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 */
__attribute__((always_inline)) INLINE static int feedback_do_feedback(
    const struct spart* sp) {

  return (sp->birth_time != -1.);
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param time The current simulation time (Non-cosmological runs).
 * @param cosmo The cosmological model (cosmological runs).
 * @param with_cosmology Are we doing a cosmological run?
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const float time, const struct cosmology* cosmo,
    const int with_cosmology) {

  if (sp->birth_time == -1.) return 0;

  if (with_cosmology) {
    return ((float)cosmo->a) > sp->birth_scale_factor;
  } else {
    return time > sp->birth_time;
  }
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.to_collect.enrichment_weight_inv = 0.f;
  sp->feedback_data.to_collect.ngb_mass = 0.f;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the distribution weights */
  sp->feedback_data.to_distribute.enrichment_weight = 0.f;

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.to_distribute.mass = 0.f;

  /* Zero the metal enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.to_distribute.metal_mass[i] = 0.f;
  }
  sp->feedback_data.to_distribute.total_metal_mass = 0.f;
  sp->feedback_data.to_distribute.mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_AGB = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNII = 0.f;
  sp->feedback_data.to_distribute.mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.metal_mass_from_SNIa = 0.f;
  sp->feedback_data.to_distribute.Fe_mass_from_SNIa = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.to_distribute.energy = 0.f;

  /* Zero the SNII feedback probability */
  sp->feedback_data.to_distribute.SNII_heating_probability = 0.f;

  /* Zero the SNII feedback energy */
  sp->feedback_data.to_distribute.SNII_delta_u = 0.f;
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_first_init_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {

  feedback_init_spart(sp);
}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_spart(
    struct spart* sp, const struct feedback_props* feedback_props) {}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 */
__attribute__((always_inline)) INLINE static void feedback_evolve_spart(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const double star_age_beg_step, const double dt) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  /* Compute amount of enrichment and feedback that needs to be done in this
   * step */
  compute_stellar_evolution(feedback_props, cosmo, sp, us, star_age_beg_step,
                            dt);

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= sp->feedback_data.to_distribute.mass;
}

void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);

void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);

#endif /* SWIFT_FEEDBACK_EAGLE_H */
