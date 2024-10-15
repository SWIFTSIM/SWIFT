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
#ifndef SWIFT_FEEDBACK_EAGLE_THERMAL_H
#define SWIFT_FEEDBACK_EAGLE_THERMAL_H

#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "rays.h"
#include "units.h"

#include <strings.h>

void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin);

/**
 * @brief Update the properties of a particle fue to feedback effects after
 * the cooling was applied.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void feedback_update_part(
    struct part* p, struct xpart* xp, const struct engine* e) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here in the EAGLE-thermal model.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const struct engine* e) {

  return e->step <= 0 ||
         ((sp->birth_time != -1.) && (sp->count_since_last_enrichment == 0));
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.to_collect.enrichment_weight_inv = 0.f;
  sp->feedback_data.to_collect.ngb_N = 0;
  sp->feedback_data.to_collect.ngb_mass = 0.f;
  sp->feedback_data.to_collect.ngb_rho = 0.f;
  sp->feedback_data.to_collect.ngb_Z = 0.f;

  /* Reset all ray structs carried by this star particle */
  ray_init(sp->feedback_data.SNII_rays, eagle_SNII_feedback_num_of_rays);

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 0;
#endif
}

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * @param sp The #spart.
 * @param with_cosmology Are we running with cosmological time integration on?
 * @param cosmo The cosmological model.
 * @param time The current time (since the Big Bang / start of the run) in
 * internal units.
 * @param dt_star the length of this particle's time-step in internal units.
 * @return The length of the enrichment step in internal units.
 */
INLINE static double feedback_get_enrichment_timestep(
    const struct spart* sp, const int with_cosmology,
    const struct cosmology* cosmo, const double time, const double dt_star) {

  if (with_cosmology) {
    return cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->last_enrichment_time, cosmo->a);
  } else {
    return time - sp->last_enrichment_time;
  }
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

  /* Zero the SNII feedback energy */
  sp->feedback_data.to_distribute.SNII_delta_u = 0.f;

  /* Zero the SNII feedback properties */
  sp->feedback_data.to_distribute.SNII_num_of_thermal_energy_inj = 0;
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
 * @brief Prepare a #spart for the feedback task.
 *
 * In EAGLE, this function evolves the stellar properties of a #spart.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
__attribute__((always_inline)) INLINE static void feedback_prepare_feedback(
    struct spart* restrict sp, const struct feedback_props* feedback_props,
    const struct cosmology* cosmo, const struct unit_system* us,
    const struct phys_const* phys_const, const double star_age_beg_step,
    const double dt, const double time, const integertime_t ti_begin,
    const int with_cosmology) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1. && dt > 0.)
    error("Evolving a star particle that should not!");
#endif

  /* Start by finishing the loops over neighbours */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  sp->feedback_data.to_collect.ngb_rho *= h_inv_dim;
  const float rho = sp->feedback_data.to_collect.ngb_rho;
  const float rho_inv = rho != 0.f ? 1.f / rho : 0.f;
  sp->feedback_data.to_collect.ngb_Z *= h_inv_dim * rho_inv;

  /* Compute amount of enrichment and feedback that needs to be done in this
   * step */
  compute_stellar_evolution(feedback_props, phys_const, cosmo, sp, us,
                            star_age_beg_step, dt, ti_begin);

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= sp->feedback_data.to_distribute.mass;

  /* Mark this is the last time we did enrichment */
  if (with_cosmology)
    sp->last_enrichment_time = cosmo->a;
  else
    sp->last_enrichment_time = time;

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 1;
#endif
}

/**
 * @brief Will this star particle want to do feedback during the next time-step?
 *
 * This is called in the time step task and increases counters of time-steps
 * that have been performed.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param time The physical time in internal units.
 * @param with_cosmology Are we running with cosmology on?
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 */
__attribute__((always_inline)) INLINE static void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base) {

  /* Special case for new-born stars */
  if (with_cosmology) {
    if (sp->birth_scale_factor == (float)cosmo->a) {

      /* Set the counter to "let's do enrichment" */
      sp->count_since_last_enrichment = 0;

      /* Ok, we are done. */
      return;
    }
  } else {
    if (sp->birth_time == (float)time) {

      /* Set the counter to "let's do enrichment" */
      sp->count_since_last_enrichment = 0;

      /* Ok, we are done. */
      return;
    }
  }

  /* Calculate age of the star at current time */
  double age_of_star;
  if (with_cosmology) {
    age_of_star = cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->birth_scale_factor, cosmo->a);
  } else {
    age_of_star = time - (double)sp->birth_time;
  }

  /* Is the star still young? */
  if (age_of_star < feedback_props->stellar_evolution_age_cut) {

    /* Set the counter to "let's do enrichment" */
    sp->count_since_last_enrichment = 0;

  } else {

    /* Increment counter */
    sp->count_since_last_enrichment++;

    if ((sp->count_since_last_enrichment %
         feedback_props->stellar_evolution_sampling_rate) == 0) {

      /* Reset counter */
      sp->count_since_last_enrichment = 0;
    }
  }
}

void feedback_clean(struct feedback_props* fp);

void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);

void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);

#ifdef HAVE_HDF5
/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The properties of the feedback scheme.
 * @param h_grp The HDF5 group in which to write.
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "EAGLE (thermal)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_EAGLE_THERMAL_H */
