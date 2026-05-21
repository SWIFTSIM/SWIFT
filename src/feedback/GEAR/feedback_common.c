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

#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "hydro_properties.h"
#include "part.h"
#include "radiation.h"
#include "stellar_evolution.h"
#include "units.h"

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

  /* TODO: Compute timestep for feedback */
  float dt = FLT_MAX;

  /*----------------------------------------*/
  /* Timestep based on the evolutionnary stage */
  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metallicity =
      chemistry_get_star_total_iron_mass_fraction_for_feedback(sp);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then  sp is a first star particle. */
  const int is_first_star = metallicity < threshold;
  const struct stellar_model *sm =
      is_first_star ? &feedback_props->stellar_model_first_stars
                    : &feedback_props->stellar_model;

  /* Compute the times */
  double star_age_beg_step = 0;

  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);
  double star_age_end_step =
      compute_star_age_end_of_step(sp, with_cosmology, cosmo, time);

  /* Convert mass to M_sun. The lifetime function assumes solar masses */
  const float log_mass =
      (sp->star_type == single_star)
          ? log10(sp->sf_data.birth_mass / phys_const->const_solar_mass)
          : log10(1.0);

  const float lifetime_myr = pow(10, lifetime_get_log_lifetime_from_mass(
                                         &sm->lifetime, log_mass, metallicity));
  const float lifetime = lifetime_myr * 1e6 * phys_const->const_year;

  /* Adapt the factor depending on the star lifetime to provide adequate
     timesteps for different star liftimes. */
  float factor = 0.0;
  if (lifetime_myr >= 100) {
    factor = 1;
  } else if (lifetime_myr < 100 && lifetime_myr >= 50) {
    factor = 20;
  } else {
    factor = 300;
  }

  /* Ensure that the age is positive (rounding errors) */
  const double star_age_beg_step_safe =
      star_age_beg_step < 0 ? star_age_end_step : star_age_beg_step;

  /* To avoid very small timesteps for star_age_beg_step, take the mean */
  const double star_age = 0.5 * (star_age_beg_step_safe + star_age_end_step);

  float dt_evolution =
      (star_age_beg_step <= 0) ? FLT_MAX : star_age / (factor * lifetime);
  /*----------------------------------------*/

  /* If the star is dead, do not limit its timestep */
  if (sp->feedback_data.is_dead) {
    return FLT_MAX;
  } else {
    dt = min(dt, dt_evolution);
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
  sp->feedback_data.supernovae.energy_ejected = 0;
  sp->feedback_data.winds.energy_ejected = 0;
  sp->feedback_data.will_do_feedback = 0;
  sp->feedback_data.will_do_HII_ionization = 0;

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
    sp->feedback_data.will_do_HII_ionization = 1;
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
  sp->feedback_data.supernovae.energy_ejected *=
      feedback_props->supernovae_efficiency;

  /* Multiply pre-SN energy by the efficiency */
  sp->feedback_data.winds.energy_ejected *= feedback_props->winds_efficiency;

  /* Apply the radiation pressure efficiency factor */
  sp->feedback_data.radiation.L_bol *=
      feedback_props->radiation_pressure_efficiency;

  /* Set the particle as doing some feedback */
  sp->feedback_data.will_do_feedback =
      sp->feedback_data.supernovae.energy_ejected != 0. ||
      sp->feedback_data.winds.energy_ejected != 0. ||
      !sp->feedback_data.is_dead;

  /* Do we need to do HII ionization feedback? */
  const char do_photoionization =
      feedback_props->radiation_policy & radiation_policy_photoionization;
  const char has_enough_photons = feedback_get_star_ionization_rate(sp) > 0.0;
  sp->feedback_data.will_do_HII_ionization =
      do_photoionization && (!sp->feedback_data.is_dead && has_enough_photons);

  /* TODO: Add other conditions:
     1. If the star is dead
     2. Or it is too old to continue HII ionization
  This star stopped doing HII for the rest of its life. Hence, h_hii = 0. If
     we want to keep track of the last h_hii, we can store it in the tracers. */
  if (sp->feedback_data.is_dead && sp->h_hii != 0.0) {
    /* Store the value in the tracers */
    sp->tracers_data.final_HII_radius = sp->h_hii * kernel_gamma;

    /* Reset to 0. This prevents dead particles with large h_hii to force the
       tasks at higher levels than necessary. */
    sp->h_hii = 0.0;
  }
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
 * Get the #spart ionization photon emission rate.
 *
 * @param sp The star.
 * @return Ionizing photon rate.
 */
__attribute__((always_inline)) INLINE double feedback_get_star_ionization_rate(
    const struct spart *sp) {
  return sp->feedback_data.radiation.dot_N_ion;
}

/**
 * @brief Should this particle be doing any HII ionization feedback-related
 * operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_HII_ionization_active(const struct spart *sp,
                                      const struct engine *e) {

  /* the particle is inactive if its birth_scale_factor or birth_time is
   * negative */

  /* If the spart is dead, don't do anything */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_HII_ionization;
}

/**
 * Determines whether a gas #part can be ionized.
 *
 * @param phys_const Physical constants.
 * @param us Unit system.
 * @param hydro_properties The #hydro_props.
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param p The particle.
 * @param xp The extended data of the particle.
 * @return Is the particle ionized?
 */
__attribute__((always_inline)) INLINE char feedback_part_can_be_ionized(
    const struct part *p, const struct xpart *xp, const struct engine *e) {

  const struct phys_const *phys_const = e->physical_constants;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct feedback_props *feedback_props = e->feedback_props;
  const struct unit_system *us = e->internal_units;
  const struct cosmology *cosmo = e->cosmology;
  const struct cooling_function_data *cooling = e->cooling_func;

  /* Is T > 10^4 K ? */
  const float T = cooling_get_temperature(phys_const, hydro_props, us, cosmo,
                                          cooling, p, xp);
  const float ten_to_four_kelvin =
      1e4 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* The 1.1 factor is here for safety margin and numerical stability */
  const char is_cold = (T <= 1.01 * ten_to_four_kelvin);

  /* Density threshold criterion */
  const float rho = hydro_get_physical_density(p, cosmo);
  const float rho_threshold = feedback_props->HII_region_min_density;
  const char is_dense = rho >= rho_threshold;

  /* Can the particle be ionized? */
  return (is_cold && is_dense && !radiation_is_part_tagged_as_ionized(p, xp));
}

/**
 * @brief Perform the ionization of a gas particle by a star.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being ionized.
 * @param xpj The #xpart (gas) being ionized.
 * @param r2 Squared distance between star and gas.
 * @param phys_const Physics constants.
 * @param hydro_props Hydrodynamics properties.
 * @param us Internal unit system.
 * @param cosmo Cosmology.
 * @param cooling Cooling function data.
 * @param ti_begin Integer time at the start of the step (for RNG).
 */
__attribute__((always_inline)) INLINE void feedback_iact_HII_ionization(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, const struct phys_const *phys_const,
    const struct hydro_props *hydro_props, const struct unit_system *us,
    const struct cosmology *cosmo, const struct cooling_function_data *cooling,
    const integertime_t ti_begin) {

  /* If already ionized (by another thread), just move on. */
  if (radiation_is_part_tagged_as_ionized(pj, xpj)) return;

  const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(
      phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

  /* Case 1: Ionization is guaranteed */
  if (Delta_dot_N_ion <= feedback_get_star_ionization_rate(si)) {
    if (atomic_cas(&xpj->tracers_data.HII_region.is_ionized, 0, 1) == 0) {

      /* Flag the particle to be synchronized on the timeline.
         We still use atomic_or here because other stars might be
         tripping the limiter in feedback loop. */
      atomic_or(&pj->limiter_data.to_be_synchronized, 1);

      /* Add the star ID */
      xpj->tracers_data.HII_region.star_id = si->id;

      /* Consume photons from the star */
      radiation_consume_ionizing_photons(si, Delta_dot_N_ion);

      /* Update HII region properties */
      si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
      si->h_hii = max(si->h_hii, sqrtf(r2) * kernel_gamma_inv);
    }
    /* If CAS failed, someone else grabbed it; just continue. */
  } else {

    /* If we cannot fully ionize, compute a probability to determine if
       we fully ionize pj or not and draw the random number.  */
    const double dot_N_ion = feedback_get_star_ionization_rate(si);
    const float proba = dot_N_ion / Delta_dot_N_ion;
    const float random_number =
        random_unit_interval(si->id, ti_begin, random_number_HII_regions);

    /* If we are lucky, do the ionization */
    if (random_number <= proba) {
      /* We won the roll! Now try to claim the particle. */
      if (atomic_cas(&xpj->tracers_data.HII_region.is_ionized, 0, 1) == 0) {
        atomic_or(&pj->limiter_data.to_be_synchronized, 1);

        xpj->tracers_data.HII_region.star_id = si->id;

        /* The star is now empty of photons */
        radiation_consume_ionizing_photons(si, Delta_dot_N_ion);
        si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
        si->h_hii = max(si->h_hii, sqrtf(r2) * kernel_gamma_inv);
      }
    } else {
      /* We lost the roll. We still consume the remaining photons. */
      radiation_consume_ionizing_photons(si, Delta_dot_N_ion);
    }
  } /* End of probability handling */
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
  sp->feedback_data.supernovae.energy_ejected = 0;

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
  sp->feedback_data.will_do_HII_ionization = 0;

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

  /* Zero energies and masses */
  sp->feedback_data.supernovae.energy_ejected = 0.0;
  sp->feedback_data.supernovae.mass_ejected = 0.0;
  sp->feedback_data.winds.energy_ejected = 0.0;
  sp->feedback_data.winds.mass_ejected = 0.0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;
  sp->feedback_data.will_do_HII_ionization = 1;
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
