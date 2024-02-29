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
#ifndef SWIFT_FEEDBACK_SIMBA_H
#define SWIFT_FEEDBACK_SIMBA_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @brief Calculates speed particles will be kicked based on
 * host galaxy properties
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_kick_speed(struct spart* sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo) {

  /* Calculate circular velocity based on Baryonic Tully-Fisher relation
  const float v_circ = pow(sp->feedback_data.host_galaxy_mass /
                               feedback_props->simba_host_galaxy_mass_norm,
                           feedback_props->simba_v_circ_exp);*/
  const float sigma = sqrtf(sp->potential * cosmo->a_inv);

  /* Calculate wind speed  (orig Simba way)
  const float random_num = 1.;
  sp->feedback_data.to_distribute.v_kick =
      feedback_props->galsf_firevel *
      pow(v_circ * cosmo->a / feedback_props->scale_factor_norm,
          feedback_props->galsf_firevel_slope) *
      pow(feedback_props->scale_factor_norm,
          0.12 - feedback_props->galsf_firevel_slope) *
      (1. - feedback_props->vwvf_scatter -
       2. * feedback_props->vwvf_scatter * random_num) *
      v_circ; */

  // ALEXEI: temporarily set to arbitrary number for testing.
  sp->feedback_data.to_distribute.v_kick =
      feedback_props->scale_factor_norm *
      pow(sigma, feedback_props->galsf_firevel_slope);
}

/**
 * @brief Calculates speed particles will be kicked based on
 * host galaxy properties
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_mass_loading(struct spart* sp,
                                 const struct feedback_props* feedback_props){};

/**
 * @brief Calculates speed particles will be kicked based on
 * host galaxy properties
 *
 * @param sp The sparticle doing the feedback
 * @param feedback_props The properties of the feedback model
 */
inline void compute_heating(struct spart* sp,
                            const struct feedback_props* feedback_props){};

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
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param time The current simulation time (Non-cosmological runs).
 * @param cosmo The cosmological model (cosmological runs).
 * @param with_cosmology Are we doing a cosmological run?
 */
/*__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const float time, const struct cosmology* cosmo,
    const int with_cosmology) {*/
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const struct engine* e) {

  return (sp->birth_time != -1.);
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {
  // Temporarily set the particle's galaxy host mass artificially.
  sp->feedback_data.host_galaxy_mass = 1.;
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
    return time - (double)sp->last_enrichment_time;
  }
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {}

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
  sp->feedback_data.to_distribute.simba_delay_time =
      feedback_props->simba_delay_time;
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
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  /* Calculate the velocity to kick neighbouring particles with */
  compute_kick_speed(sp, feedback_props, cosmo);

  /* Compute wind mass loading */
  compute_mass_loading(sp, feedback_props);

  /* Compute residual heating */
  compute_heating(sp, feedback_props);

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
    const integertime_t ti_current, const double time_base) {}

static INLINE void feedback_clean(struct feedback_props* fp) {}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
static INLINE void feedback_struct_dump(const struct feedback_props* feedback,
                                        FILE* stream) {}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
static INLINE void feedback_struct_restore(struct feedback_props* feedback,
                                           FILE* stream) {}

#ifdef HAVE_HDF5
/**
 * @brief Writes the current model of feedback to the file
 *
 * @param feedback The properties of the feedback scheme.
 * @param h_grp The HDF5 group in which to write.
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "SIMBA (kinetic decoupled)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_SIMBA_H */
