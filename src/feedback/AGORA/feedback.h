/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.reavz@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_AGORA_H
#define SWIFT_FEEDBACK_AGORA_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @brief Update the properties of a particle fue to feedback effects after
 * the cooling was applied.
 *
 * Nothing to do here.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param cosmo The #cosmology.
 */
void feedback_update_part(struct part* p, struct xpart* xp,
                          const struct engine* e);

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart* sp);

void feedback_init_after_star_formation(
    struct spart* sp, const struct feedback_props* feedback_props);

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 */
__attribute__((always_inline)) INLINE static int feedback_do_feedback(
    const struct spart* sp) {

  return 0;
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * Note: Since this 'none' feedback mode is used for testing the neighbour
 * loops only, we want to always do feedback irrespective of the particle
 * or of the system's state.
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart* sp, const struct engine* e);

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * We just return the normal time-step here since particles do something every
 * regular time-step.
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

  /* Just return the regular step length */
  return dt_star;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
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
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props);

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props);

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 * This function allows for example to compute the SN rate before sending
 * this information to a different MPI rank.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
void feedback_prepare_feedback(struct spart* restrict sp,
                               const struct feedback_props* feedback_props,
                               const struct cosmology* cosmo,
                               const struct unit_system* us,
                               const struct phys_const* phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology);

/**
 * @brief Will this star particle want to do feedback during the next time-step?
 *
 * This is called in the time step task.
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
void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base);

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props* feedback);

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream);

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 * @param cosmo #cosmology structure
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream);

#ifdef HAVE_HDF5
/**
 * @brief Writes the current model of feedback to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void feedback_write_flavour(struct feedback_props* feedback,
                                          hid_t h_grp) {

  io_write_attribute_s(h_grp, "Feedback Model", "AGORA");
};
#endif

#endif /* SWIFT_FEEDBACK_AGORA_H */
