/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "feedback.h"

/* Local includes */
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "stellar_evolution.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part* restrict p, struct xpart* restrict xp,
                          const struct engine* restrict e) {

  /* Did the particle receive a supernovae */
  if (xp->feedback_data.delta_mass == 0) return;

  const struct cosmology* cosmo = e->cosmology;

  /* Turn off the cooling */
  xp->cooling_data.time_last_event = e->time;

  /* Update mass */
  const float old_mass = hydro_get_mass(p);
  const float new_mass = old_mass + xp->feedback_data.delta_mass;

  if (xp->feedback_data.delta_mass < 0.) {
    error("Delta mass smaller than 0");
  }

  hydro_set_mass(p, new_mass);

  xp->feedback_data.delta_mass = 0;

  /* Update the density */
  p->rho *= new_mass / old_mass;

  /* Update internal energy */
  const float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  const float u_new = u + xp->feedback_data.delta_u;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, u_new);

  xp->feedback_data.delta_u = 0.;

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = xp->feedback_data.delta_p[i] / new_mass;

    xp->v_full[i] += dv;
    p->v[i] += dv;

    xp->feedback_data.delta_p[i] = 0;
  }
}

/**
 * @brief Should we do feedback for this star?
 *
 * @param sp The star to consider.
 * @param feedback_props The #feedback_props.
 * @param with_cosmology Is the cosmology switch on?
 * @param cosmo The #cosmology.
 * @param time The current time.
 */
int feedback_will_do_feedback(const struct spart* sp,
                              const struct feedback_props* feedback_props,
                              const int with_cosmology,
                              const struct cosmology* cosmo,
                              const double time) {

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
int feedback_is_active(const struct spart* sp, const double time,
                       const struct cosmology* cosmo,
                       const int with_cosmology) {

  // TODO improve this with estimates for SNII and SNIa
  if (sp->birth_time == -1.) return 0;

  if (with_cosmology) {
    return ((double)cosmo->a) > sp->birth_scale_factor;
  } else {
    return time > sp->birth_time;
  }
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
double feedback_get_enrichment_timestep(const struct spart* sp,
                                        const int with_cosmology,
                                        const struct cosmology* cosmo,
                                        const double time,
                                        const double dt_star) {
  return dt_star;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart* sp) {

  sp->feedback_data.enrichment_weight = 0.f;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props) {

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;
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
void feedback_first_init_spart(struct spart* sp,
                               const struct feedback_props* feedback_props) {

  feedback_init_spart(sp);

  feedback_reset_feedback(sp, feedback_props);
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
void feedback_prepare_spart(struct spart* sp,
                            const struct feedback_props* feedback_props) {}

/**
 * @brief Evolve the stellar properties of a #spart.
 *
 * This function compute the SN rate and yields before sending
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
void feedback_evolve_spart(struct spart* restrict sp,
                           const struct feedback_props* feedback_props,
                           const struct cosmology* cosmo,
                           const struct unit_system* us,
                           const struct phys_const* phys_const,
                           const double star_age_beg_step, const double dt,
                           const double time, const integertime_t ti_begin,
                           const int with_cosmology) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  /* Reset the feedback */
  feedback_reset_feedback(sp, feedback_props);

  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */

  sp->feedback_data.enrichment_weight *= hi_inv_dim;

  /* Compute the stellar evolution */
  stellar_evolution_evolve_spart(sp, &feedback_props->stellar_model, cosmo, us,
                                 phys_const, ti_begin, star_age_beg_step, dt);

  /* Transform the number of SN to the energy */
  sp->feedback_data.energy_ejected =
      sp->feedback_data.number_sn * feedback_props->energy_per_supernovae;
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  restart_write_blocks((void*)feedback, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");

  stellar_evolution_dump(&feedback->stellar_model, stream);
}

/**
 * @brief Restore a feedback struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {

  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  stellar_evolution_restore(&feedback->stellar_model, stream);
}

/**
 * @brief Clean the allocated memory.
 *
 * @param feedback the #feedback_props.
 */
void feedback_clean(struct feedback_props* feedback) {

  stellar_evolution_clean(&feedback->stellar_model);
}
