/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#include "feedback_common.h"

#include "cosmology.h"
#include "engine.h"
#include "hydro_properties.h"
#include "part.h"
#include "stellar_evolution.h"
#include "units.h"
#include "lifetime.h"

/**
 * @brief Computes the time-step length of a given star particle from feedback
 * physics
 *
 * @param sp Pointer to the s-particle data.
 * @param feedback_props Properties of the feedback model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 * @param ti_current The current time (in integer).
 * @param time  The current time (in double, used if running without cosmology).
 * @param time_base The time base.
 */
float feedback_compute_spart_timestep(
    const struct spart *const sp, const struct feedback_props *feedback_props,
    const struct phys_const *phys_const, const struct unit_system *us,
    const int with_cosmology, const struct cosmology *cosmo,
    const integertime_t ti_current, const double time, const double time_base) {


  /* Return FLT_MAX if in the initial fake step, as no feedback properties were calculated. Prevent timestep = 0.0*/
  if (sp->time_bin == 0) return FLT_MAX;

  /* ------------- CFL timestep for Stellar winds. -------------------- */

  float dt_cfl = 0.0;
  const float CFL_condition = feedback_props->CFL_condition;

  /* Convert relative velocity to physical (without the Hubble flow contribution, here we want the local relative velocity between gas and the star) */
  const double gas_v_phys[3] = {
      sp->feedback_data.relative_velocity_gas[0] * cosmo->a_inv,
      sp->feedback_data.relative_velocity_gas[1] * cosmo->a_inv,
      sp->feedback_data.relative_velocity_gas[2] * cosmo->a_inv};
  double gas_v_norm2 = gas_v_phys[0] * gas_v_phys[0] +
                       gas_v_phys[1] * gas_v_phys[1] +
                       gas_v_phys[2] * gas_v_phys[2];

  if (gas_v_norm2 == 0.0) {
    dt_cfl = FLT_MAX;
  } else {
    /* Compute the physical sound speed */
    const double gas_c_phys =
        sp->feedback_data.sound_speed_gas * cosmo->a_factor_sound_speed;
    const double gas_c_phys2 = gas_c_phys * gas_c_phys;
    const float denominator = sqrtf(gas_c_phys2 + gas_v_norm2);
    const float h_min =
        cosmo->a * kernel_gamma * min(sp->h, sp->feedback_data.minimal_h_gas);
    dt_cfl = 2.f * CFL_condition * h_min / denominator;
  

    //message("h_min = %e, CFL timestep = %e", h_min, dt_cfl);
  }

  /* ------------- Mass-Loss based timestep for Stellar winds. -------------------- */
  float dt_mass_loss = FLT_MAX;

  if (sp->feedback_data.preSN.mass_dot > 0) {
    const float mass_loss_condition = feedback_props->mass_loss_condition;
    dt_mass_loss = mass_loss_condition * min(sp->mass / sp->feedback_data.preSN.mass_dot, sp->feedback_data.total_gas_mass / sp->feedback_data.preSN.mass_dot);
    //message("init_mass = %e, mass_dot = %e, gas_mass = %e, dt_mass_loss = %e", sp->mass, sp->feedback_data.preSN.mass_dot, sp->feedback_data.total_gas_mass, dt_mass_loss);

  }


  /* ------------- Momentum based timestep for Stellar winds. -------------------- */

  float dt_momentum = FLT_MAX;
  float momentum_condition = feedback_props->momentum_condition;


  if (gas_v_norm2 > 0 && sp->feedback_data.preSN.mass_dot > 0) {
    const double wind_momentum_rate = sqrt(2. * sp->feedback_data.preSN.energy_dot * sp->feedback_data.preSN.mass_dot);

    if (wind_momentum_rate > 0 && sp->feedback_data.total_gas_mass > 0) {
      const double gas_c_phys =
        sp->feedback_data.sound_speed_gas * cosmo->a_factor_sound_speed;
      const double v_gas = min(sqrtf(gas_v_norm2), gas_c_phys);

      dt_momentum = momentum_condition * (sp->feedback_data.total_gas_mass * v_gas) / wind_momentum_rate;

      //message("wind_momentum_rate = %e, v_gas = %e, dt_momentum = %e", wind_momentum_rate, v_gas, dt_momentum);
    }
  }


  /* ------------- Energy based timestep for Stellar winds. -------------------- */

  float dt_energy = FLT_MAX;

  if (sp->feedback_data.preSN.energy_dot > 0) {
    const float energy_condition = feedback_props->energy_condition;
    dt_energy = energy_condition * (sp->feedback_data.total_internal_energy_gas + sp->feedback_data.total_kinetic_energy_gas) / sp->feedback_data.preSN.energy_dot;
    //message("total_internal_energy_gas = %e, total_kinetic_energy_gas = %e, energy_dot = %e, dt_energy = %e", sp->feedback_data.total_internal_energy_gas, sp->feedback_data.total_kinetic_energy_gas, sp->feedback_data.preSN.energy_dot, dt_energy);
  }

  /* If the star is dead, do not limit its timestep */
  if (sp->feedback_data.is_dead) {
    return FLT_MAX;
  } else {
    float dt_1 = min(dt_cfl, dt_mass_loss);
    float dt_2 = min(dt_energy, dt_momentum);
    float dt = min(dt_1, dt_2);
    // if (dt != FLT_MAX) {
    //   message("dt_CFL=%e; dt_mass=%e; dt_momentum=%e; dt_energy=%e; Returning timestep=%e", dt_cfl, dt_mass_loss, dt_momentum, dt_energy, dt);
    // }
    return dt;
  }
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
    struct spart *sp, const struct feedback_props *feedback_props,
    const int with_cosmology, const struct cosmology *cosmo, const double time,
    const struct unit_system *us, const struct phys_const *phys_const,
    const integertime_t ti_current, const double time_base) {

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;
  sp->feedback_data.preSN.energy_ejected = 0;
  sp->feedback_data.will_do_feedback = 0;

  /* Quit if the birth_scale_factor or birth_time is negative.
     No Feedback event for the initial fake step. */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0 || sp->time_bin == 0)
    return;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_star_total_iron_mass_fraction_for_feedback(sp);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then  sp is a first star particle. */
  const int is_first_star = metal < threshold;
  const struct stellar_model *model =
      is_first_star ? &feedback_props->stellar_model_first_stars
                    : &feedback_props->stellar_model;

  /* Compute the times */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);

  /* There is no feedback to do for newborn stars */
  const double star_age_end_step = star_age_beg_step + dt_enrichment;
  if (star_age_end_step == 0.0) {
    /* See the comment in feedback_init_after_star_formation(). But you will
       need to go through the feedback loops in the next timestep to compute
       all required quantitied for the stellar evolution. */
    sp->feedback_data.will_do_feedback = 1;
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
  if (star_age_beg_step + dt_enrichment < 0) {
    error("Negative age for a star");
  }
#endif

  /* Ensure that the age is positive (rounding errors) */
  const double star_age_beg_step_safe =
      star_age_beg_step < 0 ? 0 : star_age_beg_step;

  /* A single star */
  if (sp->star_type == single_star) {
    /* If the star has completely exploded, do not continue. This will also
       avoid NaN values in the liftetime if the mass is set to 0. Correction
       (28.04.2024): A bug fix in the mass of the star (see stellar_evolution.c
       in stellar_evolution_compute_X_feedback_properties, X=discrete,
       continuous) has changed the mass of the star from 0 to
       discrete_star_minimal_gravity_mass. Hence the fix is propagated here. */
    if (sp->mass <= model->discrete_star_minimal_gravity_mass) {
      return;
    }

    /* Now, compute the stellar evolution state for individual star particles.
     */
    stellar_evolution_evolve_individual_star(
        sp, model, cosmo, us, phys_const,
        feedback_props->with_stellar_wind_feedback, ti_begin,
        star_age_beg_step_safe, dt_enrichment);
  } else {
    /* Compute the stellar evolution including SNe energy. This function treats
       the case of particles representing the whole IMF (star_type =
       star_population) and the particles representing only the continuous part
       of the IMF (star_type = star_population_continuous_IMF) */
    stellar_evolution_evolve_spart(sp, model, cosmo, us, phys_const,
                                   feedback_props->with_stellar_wind_feedback,
                                   ti_begin, star_age_beg_step_safe,
                                   dt_enrichment);
  }

  /* Apply the energy efficiency factor */
  sp->feedback_data.energy_ejected *= feedback_props->supernovae_efficiency;

  /* Multiply pre-SN energy by the efficiency */
  sp->feedback_data.preSN.energy_ejected *= feedback_props->preSN_efficiency;

  /* Set the particle as doing some feedback */
  sp->feedback_data.will_do_feedback =
      sp->feedback_data.energy_ejected != 0. ||
      sp->feedback_data.preSN.energy_ejected != 0. ||
      !sp->feedback_data.is_dead;
}

/**
 * @brief Compute age of the star at the end of the current timestep.
 *
 * @param sp The #spart to act upon
 * @param with_cosmology Are we running with the cosmological expansion?
 * @param cosmo The current cosmological model.
 * @param time The current time (in double)
 */
double compute_star_age_end_of_step(const struct spart *sp,
                                    const int with_cosmology,
                                    const struct cosmology *cosmo,
                                    const double time) {
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
  return star_age_end_of_step;
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
void compute_time(const struct spart *sp, const int with_cosmology,
                  const struct cosmology *cosmo, double *star_age_beg_of_step,
                  double *dt_enrichment, integertime_t *ti_begin_star,
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
  const double star_age_end_of_step =
      compute_star_age_end_of_step(sp, with_cosmology, cosmo, time);

  /* Get the length of the enrichment time-step */
  *dt_enrichment = feedback_get_enrichment_timestep(sp, with_cosmology, cosmo,
                                                    time, dt_star);

  *star_age_beg_of_step = star_age_end_of_step - *dt_enrichment;
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
double feedback_get_enrichment_timestep(const struct spart *sp,
                                        const int with_cosmology,
                                        const struct cosmology *cosmo,
                                        const double time,
                                        const double dt_star) {
  return dt_star;
}

/**
 * @brief Prepare the feedback fields after a star is born.
 *
 * This function is called in the functions sink_copy_properties_to_star() and
 * star_formation_copy_properties().
 *
 * @param sp The #spart to act upon.
 * @param feedback_props The feedback perties to use.
 * @param star_type The stellar particle type.
 */
void feedback_init_after_star_formation(
    struct spart *sp, const struct feedback_props *feedback_props,
    const enum stellar_type star_type) {

  feedback_init_spart(sp);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* The star has nothing useful to do in this loop. Note that in GEAR, the
  order of operations are:
  1. Star formation: Form a star with age_beg_step < 0 and age_end_step = 0.
  2. Stars density, prep1-4, feedback apply: Nothing to do or distribute.
  sp->feedback_data.will_do_feedback = 0;
  3. Timestep: Call to feedback_will_do_feedback(), which calls the stellar
  evolution to be distributed in the next step. Since we have age_beg_step < 0
  and age_end_step = 0 now, there is nothing to compute or distribute for the
  next timestep. So we do not need to compute any feedback or stellar evolution
  now.
  */
  sp->feedback_data.will_do_feedback = 0;

  /* Give to the star its appropriate type: single star, continuous IMF star or
     single population star */
  sp->star_type = star_type;
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
void feedback_first_init_spart(struct spart *sp,
                               const struct feedback_props *feedback_props) {
  /* Initialize the feedback struct for the first time */
  feedback_init_spart(sp);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props *feedback, FILE *stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* Zero the stellar_evolution */
  stellar_evolution_zero_pointers(feedback_copy.stellar_model);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_zero_pointers(feedback_copy.stellar_model_first_stars);
  }

  restart_write_blocks((void *)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");

  /* Now dump the stellar evolution */
  stellar_evolution_dump(&feedback->stellar_model, stream);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_dump(&feedback->stellar_model_first_stars, stream);
  }
}

/**
 * @brief Restore a feedback struct from the given FILE as a stream of
 * bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props *feedback, FILE *stream) {

  restart_read_blocks((void *)feedback, sizeof(struct feedback_props), 1,
                      stream, NULL, "feedback function");

  stellar_evolution_restore(&feedback->stellar_model, stream,
                            feedback->with_stellar_wind_feedback);

  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_restore(&feedback->stellar_model_first_stars, stream,
                              feedback->with_stellar_wind_feedback);
  }
}

/**
 * @brief Clean the allocated memory.
 *
 * @param feedback the #feedback_props.
 */
void feedback_clean(struct feedback_props *feedback) {

  stellar_evolution_clean(&feedback->stellar_model);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_clean(&feedback->stellar_model_first_stars);
  }
}
