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
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "part.h"
#include "physical_constants.h"
#include "radiation.h"
#include "stellar_evolution.h"
#include "stromgren_sphere.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Computes the time-step length of a given star particle from feedback
 * physics
 *
 * @param sp Pointer to the s-particle data.
 * @param feedback_props Properties of the feedback model.
 * @param phys_const The #phys_const.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 */
float feedback_compute_spart_timestep(
    const struct spart* const sp, const struct feedback_props* feedback_props,
    const struct phys_const* phys_const, const struct unit_system* us,
    const int with_cosmology, const struct cosmology* cosmo,
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
  const struct stellar_model* sm =
      is_first_star ? &feedback_props->stellar_model_first_stars
                    : &feedback_props->stellar_model;

  /* Compute the times */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);

  const float log_mass =
      (sp->star_type == single_star)
          ? log10(sp->sf_data.birth_mass / phys_const->const_solar_mass)
          : log10(sm->imf.mass_min / phys_const->const_solar_mass);

  const float lifetime_myr = pow(10, lifetime_get_log_lifetime_from_mass(
                                         &sm->lifetime, log_mass, metallicity));
  const float lifetime = lifetime_myr * 1e6 * phys_const->const_year;

  float factor = 0.0;

  if (lifetime_myr >= 100) {
    factor = 1;
  } else if (lifetime_myr < 100 && lifetime_myr >= 50) {
    factor = 20;
  } else {
    factor = 300;
  }

  float dt_evolution = (star_age_beg_step == 0)
                           ? FLT_MAX
                           : star_age_beg_step / (factor * lifetime);

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
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 */
void feedback_update_part(struct part* p, struct xpart* xp,
                          const struct engine* e) {
  const struct cosmology* cosmo = e->cosmology;
  const struct pressure_floor_props* pressure_floor = e->pressure_floor_props;

  /* TODO: Treat the pre-SN case
     WARNING: Do not comment out this line, because it will mess-up with
     SF/sinks. (I think it injects something that it should not...) */
  /* Did the particle receive a supernovae */
  if (xp->feedback_data.delta_mass != 0) {

    /* Turn off the cooling */
    cooling_set_part_time_cooling_off(p, xp, e->time);

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

  /*----------------------------------------*/
  /* Treat ionisation */
  if (radiation_is_part_tagged_as_ionized(p, xp)) {
    const struct unit_system* us = e->internal_units;
    const struct phys_const* phys_const = e->physical_constants;
    const struct hydro_props* hydro_props = e->hydro_properties;
    const struct cooling_function_data* cooling = e->cooling_func;
    const double m_p = phys_const->const_proton_mass;
    const double k_B = phys_const->const_boltzmann_k;

    /* Get the current internal energy */
    const float u = hydro_get_physical_internal_energy(p, xp, cosmo);

    /* Get the internal energy increase to ionization */
    const double N_H = radiation_get_part_number_hydrogen_atoms(
        phys_const, hydro_props, us, cosmo, cooling, p, xp);
    const double E_ion =
        2.17872e-11 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
    const double Delta_u_ionized = N_H * E_ion / hydro_get_mass(p);

    /* Get internal energy due to collisions */
    const double Z = chemistry_get_total_metal_mass_fraction_for_feedback(p);
    const double Z_sun = 0.02;
    const double mu = cooling_get_mean_molecular_weight(
        phys_const, us, cosmo, hydro_props, cooling, p, xp);

    /* Here we need to treat the cases Z << Z_sun otherwise we have T<0*/
    const double ten_to_four_K =
        1e4 * units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
    const double tmp = 0.86 / (1 + 0.22 * log(Z / Z_sun));
    const double T_collisional = ten_to_four_K * min(6.62, tmp);
    const double u_collisional =
        cooling_internal_energy_from_T(T_collisional, mu, k_B, m_p);

    const float u_new = min(Delta_u_ionized, u_collisional);

    message("u = %e, Delta_u_ionized = %e, u_coll = %e, T_col = %e, tmp = %e",
            u, Delta_u_ionized, u_collisional, T_collisional, tmp);

    hydro_set_physical_internal_energy(p, xp, cosmo, u_new);
    hydro_set_drifted_physical_internal_energy(p, cosmo, pressure_floor, u_new);

#if COOLING_GRACKLE_MODE > 0
    /* If we have the non-equilibrium cooling, we can also set these values and
       let grackle solve the network. */
    xp->cooling_data.HI_frac = 0;
    xp->cooling_data.HII_frac = 1;
#endif

    /* Reset the ionization tag */
    radiation_reset_part_ionized_tag(p, xp);
  }

  /*----------------------------------------*/
  /* Radiation pressure */

  if (e->feedback_props->radiation_pressure_efficiency != 0.0) {
    for (int i = 0; i < 3; i++) {
      const float dv =
          xp->feedback_data.radiation.delta_p[i] / hydro_get_mass(p);
      xp->v_full[i] += dv;
      p->v[i] += dv;

      /* Reset */
      xp->feedback_data.radiation.delta_p[i] = 0;
    }
  }
}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here in the GEAR model.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
void feedback_reset_part(struct part* p, struct xpart* xp) {}

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
void compute_time(const struct spart* sp, const int with_cosmology,
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

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;
  sp->feedback_data.will_do_feedback = 0;

  /* quit if the birth_scale_factor or birth_time is negative */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return;

  /* Pick the correct table. (if only one table, threshold is < 0) */
  const float metal =
      chemistry_get_star_total_iron_mass_fraction_for_feedback(sp);
  const float threshold = feedback_props->metallicity_max_first_stars;

  /* If metal < threshold, then  sp is a first star particle. */
  const int is_first_star = metal < threshold;
  const struct stellar_model* model =
      is_first_star ? &feedback_props->stellar_model_first_stars
                    : &feedback_props->stellar_model;

  /* Compute the times */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin = 0;
  compute_time(sp, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin, ti_current, time_base, time);

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
    stellar_evolution_evolve_individual_star(sp, model, cosmo, us, phys_const,
                                             ti_begin, star_age_beg_step_safe,
                                             dt_enrichment);
  } else {
    /* Compute the stellar evolution including SNe energy. This function treats
       the case of particles representing the whole IMF (star_type =
       star_population) and the particles representing only the continuous part
       of the IMF (star_type = star_population_continuous_IMF) */
    stellar_evolution_evolve_spart(sp, model, cosmo, us, phys_const, ti_begin,
                                   star_age_beg_step_safe, dt_enrichment);
  }

  /* Apply the energy efficiency factor */
  sp->feedback_data.energy_ejected *= feedback_props->supernovae_efficiency;

  /* TODO: See if we need to add something about pre-SN */
  /* Set the particle as doing some feedback */
  sp->feedback_data.will_do_feedback =
      sp->feedback_data.energy_ejected != 0. || !sp->feedback_data.is_dead;

  /* TODO: Do we want to multiply pre-SN energy bu the efficiency? */
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart* sp, const struct engine* e) {

  /* the particle is inactive if its birth_scale_factor or birth_time is
   * negative */

  /* If the spart is dead, don't do anything */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_feedback;
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
 * Note: In GEAR, this function must not reset the data as the are computed
 * at the end of the tasks by the stellar_evolution functions.
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart* sp) {

  sp->feedback_data.enrichment_weight = 0.f;

  sp->feedback_data.rho_star = 0.0;
  sp->feedback_data.grad_rho_star[0] = 0.0;
  sp->feedback_data.grad_rho_star[1] = 0.0;
  sp->feedback_data.grad_rho_star[2] = 0.0;

  sp->feedback_data.Z_star = 0.0;

  /* Radiation fields */
  sp->feedback_data.num_ngbs = 0;
  stromgren_shell_init(sp->feedback_data.radiation.stromgren_sphere,
                       GEAR_STROMGREN_NUMBER_NEIGHBOURS);
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
    struct spart* sp, const struct feedback_props* feedback_props,
    const enum stellar_type star_type) {

  feedback_init_spart(sp);

  /* Zero the energy of supernovae */
  sp->feedback_data.energy_ejected = 0;

  /* Activate the feedback loop for the first step */
  sp->feedback_data.will_do_feedback = 1;

  /* Give to the star its appropriate type: single star, continuous IMF star or
     single population star */
  sp->star_type = star_type;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 */
void feedback_reset_feedback(struct spart* sp,
                             const struct feedback_props* feedback_props) {}

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

  /* The spart is dead if its birth_time is negative */
  sp->feedback_data.is_dead =
      (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0);
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
  const float hi_inv_dim = pow_dimension(hi_inv);        /* 1/h^d */
  const float hi_inv_dim_plus_one = hi_inv_dim * hi_inv; /* 1/h^(d+1) */
  sp->feedback_data.enrichment_weight *= hi_inv_dim;

  sp->feedback_data.rho_star *= hi_inv;
  sp->feedback_data.grad_rho_star[0] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[1] *= hi_inv_dim_plus_one;
  sp->feedback_data.grad_rho_star[2] *= hi_inv_dim_plus_one;

  sp->feedback_data.Z_star *= hi_inv / sp->feedback_data.rho_star;

  const float Sigma_gas = radiation_get_star_gas_column_density(sp);
  const float kappa_IR = radiation_get_IR_opacity(sp, us, phys_const);
  const float tau_IR = radiation_get_IR_optical_depth(sp, us, phys_const);

  message(
      "rho_star = %e, Grad rho = (%e %e %e), Z = %e, Sigma_gas = %e, kappa_IR "
      "= %e, tau_IR = %e",
      sp->feedback_data.rho_star, sp->feedback_data.grad_rho_star[0],
      sp->feedback_data.grad_rho_star[1], sp->feedback_data.grad_rho_star[2],
      sp->feedback_data.Z_star, Sigma_gas, kappa_IR, tau_IR);

  /*----------------------------------------*/
  /* Do the HII ionization */

  if (feedback_props->do_photoionization) {
    const struct stromgren_shell_data* stromgren =
        sp->feedback_data.radiation.stromgren_sphere;
    const int num_ngb =
        min(sp->feedback_data.num_ngbs, GEAR_STROMGREN_NUMBER_NEIGHBOURS);

    /* Loop over the sorted gas neighbours */
    for (int i = 0; i < num_ngb; i++) {
      /* This is recomputed at each iteration */
      const double dot_N_ion = radiation_get_star_ionization_rate(sp);

      if (dot_N_ion <= 0.0) {
        sp->feedback_data.radiation.dot_N_ion = 0.0;
        break;
      }

      const double Delta_dot_N_ion = stromgren[i].Delta_N_dot;

      if (Delta_dot_N_ion <= dot_N_ion) {
        /* We can fully ionize this particle */
        /* Update the Stromgren sphere radius */
        sp->feedback_data.radiation.R_stromgren = stromgren[i].distance;

        /* Consume the photons */
        radiation_consume_ionizing_photons(sp, Delta_dot_N_ion);
      } else {
        /* If we cannot fully ionize, compute a probability to determine if we
           fully ionize pj or not and draw the random number.  */
        const float proba = dot_N_ion / Delta_dot_N_ion;
        const float random_number =
            random_unit_interval(sp->id, ti_begin, random_number_HII_regions);

        /* If we are lucky or we are the first particle, do the ionization */
        if (random_number <= proba || i == 0) {
          /* Update the Stromgren sphere radius */
          sp->feedback_data.radiation.R_stromgren = stromgren[i].distance;
        }

        /* Consume the photons in all cases */
        radiation_consume_ionizing_photons(sp, Delta_dot_N_ion);
      }
    }
  }
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
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {

  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  stellar_evolution_restore(&feedback->stellar_model, stream);

  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_restore(&feedback->stellar_model_first_stars, stream);
  }
}

/**
 * @brief Clean the allocated memory.
 *
 * @param feedback the #feedback_props.
 */
void feedback_clean(struct feedback_props* feedback) {

  stellar_evolution_clean(&feedback->stellar_model);
  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_clean(&feedback->stellar_model_first_stars);
  }
}
