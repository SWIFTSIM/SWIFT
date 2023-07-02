/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part* p, struct xpart* xp,
                          const struct engine* e) {

  /* Did the particle receive a supernovae */
  if (xp->feedback_data.delta_mass == 0) return;

  const struct cosmology* cosmo = e->cosmology;
  const struct pressure_floor_props* pressure_floor = e->pressure_floor_props;

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
  const float u =
      hydro_get_physical_internal_energy(p, xp, cosmo) * old_mass / new_mass;
  const float u_new = u + xp->feedback_data.delta_u;

  hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

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
 * @brief Compute the times for the stellar model.
 *
 * This function assumed to be called in the time step task.
 *
 * @param sp The #spart to act upon
 * @param with_cosmology Are we running with the cosmological expansion?
 * @param cosmo The current cosmological model.
 * @param star_age_beg_of_step (output) Age of the star at the beginning of the
 * step.
 * @param dt_enrichment (output) Time step for the stellar evolution.
 * @param ti_begin_star (output) Integer time at the beginning of the time step.
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 * @param time The current time (in double)
 */
void compute_time(struct spart* sp, const int with_cosmology,
                  const struct cosmology* cosmo, double* star_age_beg_of_step,
                  double* dt_enrichment, integertime_t* ti_begin_star,
                  const integertime_t ti_current, const double time_base,
                  const double time) {
  const integertime_t ti_step = get_integer_timestep(sp->time_bin);
  *ti_begin_star = get_integer_time_begin(ti_current, sp->time_bin);

  /* Get particle time-step */
  double dt_star;
  if (with_cosmology) {
    dt_star = cosmology_get_delta_time(cosmo, *ti_begin_star,
                                       *ti_begin_star + ti_step);
  } else {
    dt_star = get_timestep(sp->time_bin, time_base);
  }

  /* Calculate age of the star at current time */
  double star_age_end_of_step;
  if (with_cosmology) {
    if (cosmo->a > (double)sp->birth_scale_factor)
      star_age_end_of_step = cosmology_get_delta_time_from_scale_factors(
          cosmo, (double)sp->birth_scale_factor, cosmo->a);
    else
      star_age_end_of_step = 0.;
  } else {
    star_age_end_of_step = max(time - (double)sp->birth_time, 0.);
  }

  /* Get the length of the enrichment time-step */
  *dt_enrichment = feedback_get_enrichment_timestep(sp, with_cosmology, cosmo,
                                                    time, dt_star);

  *star_age_beg_of_step = star_age_end_of_step - *dt_enrichment;
}

/**
 * @brief Will this star particle want to do feedback during the next time-step?
 *
 * This is called in the time step task.
 *
 * In GEAR, we compute the full stellar evolution here.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param ti_current The current time (in integer)
 * @param time_base The time base.
 * @param time The physical time in internal units.
 */
void feedback_will_do_feedback(
    struct spart* sp, const struct feedback_props* feedback_props,
    const int with_cosmology, const struct cosmology* cosmo, const double time,
    const struct unit_system* us, const struct phys_const* phys_const,
    const integertime_t ti_current, const double time_base) {

  /* skip if the particle is idle for feedback (it already exploded) */
  if (sp->feedback_data.idle == 1) {
    sp->feedback_data.will_do_feedback = 0;
    return;
  }

  /* Compute the times */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;
  sp->feedback_data.will_do_feedback = 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
  if (star_age_beg_step + dt_enrichment < 0) {
    error("Negative age for a star");
  }
#endif

  /* Ensure that the age is positive (rounding errors) */
  const double star_age_beg_step_safe =
      star_age_beg_step < 0 ? 0 : star_age_beg_step;

  /* Wait until the star is older than the explosion time */
  if (star_age_beg_step_safe < feedback_props->supernovae_explosion_time)
    return;

  /* Number of exploding supernovae */
  double NSN = sp->mass / phys_const->const_solar_mass *
               feedback_props->ccsne_per_solar_mass;

  /* setup ejected energy */
  sp->feedback_data.energy_ejected = NSN * feedback_props->energy_per_CCSN;
  sp->feedback_data.energy_ejected *= feedback_props->supernovae_efficiency;

  /* setup ejected mass */
  sp->feedback_data.mass_ejected = NSN * feedback_props->ejected_mass_per_CCSN;

  /* setup ejected metals */
  sp->feedback_data.metal_mass_ejected[0] =
      NSN * feedback_props->ejected_Fe_mass_per_CCSN;
  sp->feedback_data.metal_mass_ejected[AGORA_CHEMISTRY_ELEMENT_COUNT - 1] =
      NSN * feedback_props->ejected_metal_mass_per_CCSN;

  /* Set the particle as doing some feedback */
  sp->feedback_data.will_do_feedback = sp->feedback_data.energy_ejected != 0.;

  /* set it as idle : the particle explode only once */
  sp->feedback_data.idle = 1;
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart* sp, const struct engine* e) {

  return sp->feedback_data.will_do_feedback;
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
 * @brief Reset the feedback field when the spart is not
 * in a correct state for feeedback_will_do_feedback.
 *
 * This function is called in the timestep task.
 */
void feedback_init_after_star_formation(
    struct spart* sp, const struct feedback_props* feedback_props) {
  feedback_init_spart(sp);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;

  /* The particle is not idle */
  sp->feedback_data.idle = 0;
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

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;

  /* The particle is not idle */
  sp->feedback_data.idle = 0;
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
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 *
 * In here, we only need to add the missing coefficients.
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
                               const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  sp->feedback_data.enrichment_weight *= hi_inv_dim;
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
}

/**
 * @brief Clean the allocated memory.
 *
 * @param feedback the #feedback_props.
 */
void feedback_clean(struct feedback_props* feedback) {}
