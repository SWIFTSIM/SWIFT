/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_FEEDBACK_KIARA_H
#define SWIFT_FEEDBACK_KIARA_H

#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"
#include "timers.h"
#include "timestep_sync.h"
#include "timestep_sync_part.h"

#include <strings.h>

double feedback_wind_probability(struct part* p, struct xpart* xp, const struct engine* e, 
                                 const struct cosmology* cosmo,
                                 const struct feedback_props* fb_props, 
                                 const integertime_t ti_current, 
                                 const double dt_part,
                                 double *rand_for_sf_wind,
                                 double *wind_mass);
void feedback_kick_and_decouple_part(struct part* p, struct xpart* xp, 
                                     const struct engine* e, 
                                     const struct cosmology* cosmo,
                                     const struct feedback_props* fb_props, 
                                     const integertime_t ti_current,
                                     const int with_cosmology,
                                     const double dt_part,
                                     const double wind_mass);
double feedback_get_lum_from_star_particle(const struct spart *sp, 
				           double age,
                                           const struct feedback_props* fb_props);
void feedback_get_ejecta_from_star_particle(const struct spart* sp,
                                            double age,
                                            const struct feedback_props* fb_props,
                                            double dt,
                                            float *N_SNe,
                                            float *ejecta_energy,
                                            float *ejecta_mass,
                                            float *ejecta_unprocessed,
                                            float ejecta_metal_mass[chem5_element_count]);
float feedback_life_time(const struct feedback_props* fb_props,
                         const float m, 
                         const float z);
float feedback_imf(const struct feedback_props* fb_props, 
                   const float m);
void feedback_set_turnover_mass(const struct feedback_props* fb_props, 
                                const float z, double* LFLT2);
float feedback_get_turnover_mass(const struct feedback_props* fb_props, 
                                 const float t, const float z);
void feedback_prepare_interpolation_tables(const struct feedback_props* fb_props);

/**
 * @brief Recouple wind particles.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_recouple_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology) {

  /* No reason to do this is the decoupling time is zero */
  if (p->feedback_data.decoupling_delay_time > 0.f) {
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, p->time_bin);

    /* Get particle time-step */
    double dt_part;
    if (with_cosmology) {
      dt_part =
          cosmology_get_delta_time(e->cosmology, ti_begin, ti_begin + ti_step);
    } else {
      dt_part = get_timestep(p->time_bin, e->time_base);
    }

    p->feedback_data.decoupling_delay_time -= dt_part;
    if (p->feedback_data.decoupling_delay_time < 0.f) {
      p->feedback_data.decoupling_delay_time = 0.f;

      /* Make sure to sync the newly coupled part on the timeline */
      timestep_sync_part(p);
    }
  } else {
    /* Because we are using floats, always make sure to set exactly zero */
    p->feedback_data.decoupling_delay_time = 0.f;
  }
}

/**
 * @brief Determine if particles that ignore cooling should start cooling again.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_ready_to_cool(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology) {

  /* No reason to do this is the decoupling time is zero */
  if (p->feedback_data.cooling_shutoff_delay_time > 0.f) {
    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const integertime_t ti_begin =
        get_integer_time_begin(e->ti_current - 1, p->time_bin);

    /* Get particle time-step */
    double dt_part;
    if (with_cosmology) {
      dt_part =
          cosmology_get_delta_time(e->cosmology, ti_begin, ti_begin + ti_step);
    } else {
      dt_part = get_timestep(p->time_bin, e->time_base);
    }

    p->feedback_data.cooling_shutoff_delay_time -= dt_part;
    if (p->feedback_data.cooling_shutoff_delay_time < 0.f) {
      p->feedback_data.cooling_shutoff_delay_time = 0.f;

      /* Make sure to sync the newly coupled part on the timeline */
      timestep_sync_part(p);
    }
  } else {
    /* Because we are using floats, always make sure to set exactly zero */
    p->feedback_data.cooling_shutoff_delay_time = 0.f;
  }
}

/**
 * @brief Update the properties of a particle due to feedback effects after
 * the cooling was applied.
 *
 * Nothing to do here in the SIMBA model.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param with_cosmology Is this a cosmological simulation?
 */
__attribute__((always_inline)) INLINE static void feedback_update_part(
    struct part* p, struct xpart* xp, const struct engine* e,
    const int with_cosmology) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * @param p The particle.
 * @param xp The extended data of the particle.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_part(
    struct part* p, struct xpart* xp) {
#if COOLING_GRACKLE_MODE >= 2
    p->feedback_data.G0 = 0.;
    p->feedback_data.SNe_ThisTimeStep = 0.;
#endif
}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int feedback_is_active(
    const struct spart* sp, const struct engine* e) {

  return 1; /* always */
}

/**
 * @brief Should this particle be doing any DM looping?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int stars_dm_loop_is_active(
    const struct spart* sp, const struct engine* e) {
  /* Active stars always do the DM loop for the SIMBA model */
  return 0;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void feedback_init_spart(
    struct spart* sp) {

  sp->feedback_data.enrichment_weight_inv = 0.f;
  sp->feedback_data.ngb_N = 0;
  sp->feedback_data.ngb_mass = 0.f;
  sp->feedback_data.ngb_rho = 0.f;
  sp->feedback_data.ngb_Z = 0.f;
#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 0;
#endif
}

/**
 * @brief Returns the length of time since the particle last did
 * enrichment/feedback.
 *
 * @param sp The #spart.
 * @param dm_ngb_N the integer number of neighbours from the previous loop
 * @param dm_mean_velocity the mass-weighted (unnormalized) three components of
 * velocity
 */
INLINE static void feedback_intermediate_density_normalize(
    struct spart* sp, const int dm_ngb_N, float dm_mean_velocity[3]) { }

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
  /* D. Rennehan: We always use the timestep of the star particle */
  return dt_star;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 */
__attribute__((always_inline)) INLINE static void feedback_reset_feedback(
    struct spart* sp, const struct feedback_props* feedback_props) {

  /* Zero the distribution weights */
  sp->feedback_data.enrichment_weight = 0.f;

  /* Zero the amount of mass that is distributed */
  sp->feedback_data.mass = 0.f;

  /* Zero the metal enrichment quantities */
  for (int i = 0; i < chemistry_element_count; i++) {
    sp->feedback_data.metal_mass[i] = 0.f;
  }
  sp->feedback_data.total_metal_mass = 0.f;

  /* Zero the energy to inject */
  sp->feedback_data.energy = 0.f;

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
 * In SIMBA, this function evolves the stellar properties of a #spart.
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

  if (sp->feedback_data.ngb_rho <= 0.) {
    warning("Star %lld has zero neighbor gas density.", sp->id);
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->birth_time == -1.) error("Evolving a star particle that should not!");
#endif

  /* Start by finishing the loops over neighbours */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  sp->feedback_data.ngb_rho *= h_inv_dim;

  const float rho_inv = 1.f / sp->feedback_data.ngb_rho;
  sp->feedback_data.ngb_Z *= h_inv_dim * rho_inv;

  TIMER_TIC;

#ifdef SIMBA_DEBUG_CHECKS
  if (sp->feedback_data.ngb_rho <= 0) {
    warning("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
    return;
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  if (star_age_beg_step < 0.f) error("Negative age for a star.");
  if (sp->feedback_data.ngb_rho <= 0)
    error("Star %lld with mass %g has no neighbors!",
            sp->id, sp->mass);
  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.enrichment_weight < 0.)
    error("Negative weight!");
#endif

#if COOLING_GRACKLE_MODE >= 2
  /* Compute Habing luminosity of star for use in ISRF (only with Grackle subgrid ISM model) */
  sp->feedback_data.lum_habing = feedback_get_lum_from_star_particle(sp, star_age_beg_step, feedback_props);
#endif

  /* Do chem5 chemical evolution model */
  int elem;
  float N_SNe = 0.f;
  float ejecta_energy = 0.f;
  float ejecta_mass = 0.f;
  float ejecta_unprocessed = 0.f;
  float ejecta_metal_mass[chem5_element_count];
  for (elem = 0; elem < chem5_element_count; elem++) ejecta_metal_mass[elem] = 0.f;

  feedback_get_ejecta_from_star_particle(sp, star_age_beg_step, feedback_props, dt, 
		  			 &N_SNe,
                                         &ejecta_energy, 
                                         &ejecta_mass, 
                                         &ejecta_unprocessed, 
                                         ejecta_metal_mass);

  if (isnan(ejecta_mass)) {
    for (elem = 0; elem < chem5_element_count; elem++) {
      message("ejecta_metal_mass[%d]=%g", elem, ejecta_metal_mass[elem]);
    }

    message("[Fe/H] = %g", sp->chemistry_data.metal_mass_fraction[chemistry_element_Fe] / sp->chemistry_data.metal_mass_fraction[chemistry_element_H]);
    message("Z = %g", sp->chemistry_data.metal_mass_fraction_total);

    error("Star particle %lld with mass %g (init_mass %g) is trying to give away NaN mass (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, ejecta_unprocessed);
  }

  if (ejecta_energy < 0.f) {
    warning("Star particle %lld with mass %g (init_mass %g) is trying to give away negative energy (Mejecta=%g, Energy=%g, Unprocessed=%g)!",
          sp->id, sp->mass, sp->mass_init, ejecta_mass, ejecta_energy, ejecta_unprocessed);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  if (sp->mass-ejecta_mass < 0.2 * sp->mass_init) {
    warning("Star particle %lld with mass %g is trying to lower its mass past 0.2 of initial (Mejecta=%g)!",
          sp->id, sp->mass, ejecta_mass);
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* D. Rennehan: Do some magic that I still don't understand 
   * https://www.youtube.com/watch?v=cY2xBNWrBZ4
   */
  float dum = 0.f;
  int flag_negative = 0;
  /* Here we can loop over Swift metals because metal_mass_fraction 
   * would be zero for the unique Chem5 metals anyway, and would
   * not activate the condition.
   */
  for (elem = 0; elem < chemistry_element_count; elem++) {
    dum = ejecta_unprocessed * sp->chemistry_data.metal_mass_fraction[elem];
    ejecta_metal_mass[feedback_props->element_index_conversions[elem]] += dum;
    if (ejecta_metal_mass[feedback_props->element_index_conversions[elem]] < 0.f) {
      ejecta_metal_mass[feedback_props->element_index_conversions[elem]] = 0.f;
      flag_negative = 1;
      /* Do not break here, we need the zeroed elements where negative */
    }
  }
  
  /* Check for any remaining that have negative mass after adding unprocessed */
  for (elem = 0; elem < chem5_element_count; elem++) {
    if (ejecta_metal_mass[elem] < 0.f) {
      ejecta_metal_mass[elem] = 0.f;
      flag_negative = 1;
    }
  }

  /* If ANY element ended up negative we recompute everything */
  if (flag_negative) {
    ejecta_mass = 0.f;
    ejecta_metal_mass[chem5_element_Z] = 0.f;
    for (elem = chem5_element_H; elem < chem5_element_Zn; elem++) ejecta_mass += ejecta_metal_mass[elem];
    for (elem = chem5_element_C; elem < chem5_element_Zn; elem++) ejecta_metal_mass[chem5_element_Z] += ejecta_metal_mass[elem];
  }

  /* Now we loop over the Swift metals and set the proper values using the conversion map */
  sp->feedback_data.total_metal_mass = 0.f;
  for (elem = 0; elem < chemistry_element_count; elem++) {
    sp->feedback_data.metal_mass[elem] = ejecta_metal_mass[feedback_props->element_index_conversions[elem]];

    /* Only count real metals in total_metal_mass */
    if (elem != chemistry_element_H && elem != chemistry_element_He) {
      sp->feedback_data.total_metal_mass += ejecta_metal_mass[feedback_props->element_index_conversions[elem]];
    }
  }

#ifdef SIMBA_DEBUG_CHECKS
    if (sp->mass/sp->mass_init<0.2) message("Star particle %lld with mass %g (init %g) is giving away %g Msun and %g erg (%g Msun metals).",
          sp->id, 
          sp->mass, 
          sp->mass_init, 
          ejecta_mass * feedback_props->mass_to_solar_mass, 
          ejecta_energy * feedback_props->energy_to_cgs,
          sp->feedback_data.total_metal_mass * feedback_props->mass_to_solar_mass);
#endif

  /* Compute the total mass to distribute */
  sp->feedback_data.mass = ejecta_mass;

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.energy = ejecta_energy;

  /* Decrease star mass by amount of mass distributed to gas neighbours */
  sp->mass -= ejecta_mass;

#if COOLING_GRACKLE_MODE >= 2
  sp->feedback_data.SNe_ThisTimeStep = N_SNe;
#endif

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->has_done_feedback = 1;
#endif

  TIMER_TOC(timer_do_star_evol);

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
    const integertime_t ti_current, const double time_base) { }

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

  io_write_attribute_s(h_grp, "Feedback Model", "KIARA (decoupled kinetic + chem5 enrichment)");
}
#endif  // HAVE_HDF5

#endif /* SWIFT_FEEDBACK_KIARA_H */
