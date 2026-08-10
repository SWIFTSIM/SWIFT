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
#include "timestep_sync_part.h"
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
     timesteps for different star lifetimes. */
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

  const float dt_evolution =
      (star_age_beg_step <= 0) ? FLT_MAX : star_age / (factor * lifetime);
  /*----------------------------------------*/

  /* HII region constraint to rebuild the HII region every
   * feedback_props->HII_rebuild_time. */
  const float HII_region_rebuild_dt = feedback_props->HII_rebuild_time;
  const double HII_max_age = feedback_props->HII_max_age;

  float dt_HII_safe = FLT_MAX;

  if (HII_region_rebuild_dt < 0.0 && star_age_beg_step < HII_max_age) {
    /* Negative rebuild time means "rebuild every step"
       (feedback_will_do_feedback()'s need_HII_region_rebuild gate already
       treats this as unconditional). Force a small step so the star is
       reliably active every step; stars_compute_timestep()'s final floor
       guards against this ever reaching 0.0. */
    dt_HII_safe = (float)feedback_props->HII_rebuild_floor_Myr;
  } else if (HII_region_rebuild_dt > 0.0 && star_age_beg_step < HII_max_age) {
    /* HII_region_last_rebuild zero-inits with the rest of the spart
       struct and is only ever written a non-negative age (see the stamp
       in runner_dosub_stars_hii_ionization_feedback()), so treating it
       as "0.0 = never rebuilt yet" below is correct for a brand-new
       star too -- no separate first-time case is needed. */
    const double HII_region_last_rebuild =
        sp->feedback_data.radiation.HII_region_last_rebuild;

    /* Calculate the absolute target age for the next rebuild */
    const double target = HII_region_last_rebuild + HII_region_rebuild_dt;
    double next_rebuild_age;
    if (star_age_beg_step > target) {
      /* Already past due (e.g. after a lag) -- jump to the next integer
         interval instead of the immediately-overshot one. */
      const double intervals =
          ceil((star_age_beg_step - HII_region_last_rebuild) /
               HII_region_rebuild_dt);
      next_rebuild_age =
          HII_region_last_rebuild + intervals * HII_region_rebuild_dt;
    } else {
      next_rebuild_age = target;
    }

    /* Limit next rebuild age by the maximum possible HII life stage */
    if (next_rebuild_age > HII_max_age) {
      next_rebuild_age = HII_max_age;
    }

    /* The maximum allowable time-step size to avoid overshooting the
       rebuild point. Floored uniformly with every other criterion in
       stars_compute_timestep(). */
    dt_HII_safe = (float)(next_rebuild_age - star_age_beg_step);
  }

  /*----------------------------------------*/
  /* If the star is dead, do not limit its timestep */
  if (sp->feedback_data.is_dead) {
    return FLT_MAX;
  } else {
    /* Not floored here: dt_evolution and dt_HII_safe are combined with
       the age-based criterion and floored uniformly in
       stars_compute_timestep() (src/stars/GEAR/stars.h), the
       model-agnostic layer where the sibling age-based bounds
       (max_time_step_young/old) already live. */
    dt = min3(dt, dt_evolution, dt_HII_safe);
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

  /* Zero the radiation pressure luminosity. Unlike the energies above, this
     is not otherwise reset when a star is dead (mass floored, or age past
     lifetime): every early-return path below skips the fresh assignment in
     stellar_evolution_compute_preSN_feedback_*, so without this a dead star
     would keep exerting its last living L_bol as a radiation-pressure kick
     on its gas neighbours indefinitely. */
  sp->feedback_data.radiation.L_bol = 0;

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
    /* A fully-exploded discrete star's mass floors at
       discrete_star_minimal_gravity_mass (not 0, see
       stellar_evolution_compute_discrete_feedback_properties), so compare
       against that floor rather than 0 to detect it and avoid NaNs in the
       lifetime calculation below. */
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

  /* TODO: Move all these into a dedicated function...? */
  /* Do we need to do HII ionization feedback? */
  const char do_photoionization =
      feedback_props->radiation_policy & radiation_policy_photoionization;
  /* A rate question, asked before the pass opens its photon budget: is this
     star emitting at all? */
  const char has_enough_photons =
      feedback_get_star_ionization_rate(sp, 0) > 0.0;
  const char is_HII_eligible = star_age_beg_step <= feedback_props->HII_max_age;

  /* Only rebuild every HII_rebuild_time (internal units) and at the
     first timestep the star is born. A negative rebuild time means "every
     step". */
  const double HII_region_last_rebuild =
      sp->feedback_data.radiation.HII_region_last_rebuild;
  const float HII_rebuild_time = feedback_props->HII_rebuild_time;

  char need_HII_region_rebuild = 0;
  if (HII_rebuild_time < 0.0) {
    need_HII_region_rebuild = 1;
  } else {
    /* HII_region_last_rebuild zero-inits with the rest of the spart
       struct (never seeded negative), so "0.0 = never rebuilt yet" below
       already covers a brand-new star -- its very first rebuild is
       triggered separately, by the star_age_end_step == 0.0 special case
       above. */
    const double next_rebuild_target =
        HII_region_last_rebuild + HII_rebuild_time;
    const double eps = 1e-4 * HII_rebuild_time;

    if (star_age_end_step >= next_rebuild_target - eps) {
      need_HII_region_rebuild = 1;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message(
      "HII_region_last_rebuild = %e, age_beg_step = %e, age_end_step = %e, "
      "next_rebuild_time "
      "= %e, need_HII_region_rebuild = %i",
      HII_region_last_rebuild, star_age_beg_step, star_age_end_step,
      HII_region_last_rebuild + HII_rebuild_time, need_HII_region_rebuild);
#endif

  sp->feedback_data.will_do_HII_ionization =
      do_photoionization && (!sp->feedback_data.is_dead && has_enough_photons &&
                             is_HII_eligible && need_HII_region_rebuild);

  /* Reset the accumulated HII-region mass exactly when we're about to redo
     the search this step, so the value read by anything outside the
     rebuild step (e.g. a snapshot dump) reflects the last completed
     rebuild instead of always reading 0 (see feedback_init_spart(), which
     used to reset this unconditionally every step). */
  if (sp->feedback_data.will_do_HII_ionization) {
    sp->feedback_data.radiation.mass_HII_region = 0.0;
  }

  /* This star stopped doing HII for the rest of its life. Its region is no
     longer maintained and will lapse, so retire both of its extent measures
     to the tracers rather than leave snapshots reporting a region that no
     longer exists. */
  if ((sp->feedback_data.is_dead || !is_HII_eligible) && sp->h_hii != 0.0) {
    /* Store the values in the tracers */
    sp->tracers_data.final_HII_radius = sp->h_hii * kernel_gamma;
    sp->tracers_data.final_HII_mass =
        sp->feedback_data.radiation.mass_HII_region;

    /* Reset to 0. This prevents dead particles with large h_hii to force the
       tasks at higher levels than necessary. */
    sp->h_hii = 0.0;
    sp->feedback_data.radiation.mass_HII_region = 0.f;
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
 * Get the #spart ionization photon emission rate for a given angular pixel.
 *
 * A pure rate, never debited during a rebuild pass: use it for physical
 * quantities derived from the star's emission (the rate-coupled flux, the
 * bootstrap Stromgren radius), not for the pass's spending decisions --
 * those go through #feedback_get_star_ionization_budget.
 *
 * @param sp The star.
 * @param pixel The angular pixel.
 * @return Ionizing photon emission rate of that pixel.
 */
__attribute__((always_inline)) INLINE double feedback_get_star_ionization_rate(
    const struct spart *sp, int pixel) {
  return sp->feedback_data.radiation.dot_N_ion_pix[pixel];
}

/**
 * Get the #spart's remaining ionizing photon *count* for a given angular
 * pixel this rebuild pass.
 *
 * @param sp The star.
 * @param pixel The angular pixel.
 * @return Ionizing photon count remaining in that pixel.
 */
__attribute__((always_inline)) INLINE double
feedback_get_star_ionization_budget(const struct spart *sp, int pixel) {
  return sp->feedback_data.radiation.N_ion_budget_pix[pixel];
}

/**
 * Get the largest remaining ionizing photon count across all of the
 * #spart's active angular pixels. Used only for loop-termination/retry
 * decisions -- one exhausted pixel doesn't mean the star is done.
 *
 * @param sp The star.
 * @return Largest remaining ionizing photon count over all active pixels.
 */
__attribute__((always_inline)) INLINE double
feedback_get_star_ionization_budget_max(const struct spart *sp) {
  double budget_max = 0.0;
  for (int p = 0; p < sp->feedback_data.radiation.n_HII_pixels; p++) {
    budget_max =
        max(budget_max, sp->feedback_data.radiation.N_ion_budget_pix[p]);
  }
  return budget_max;
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
  const float rho_threshold = feedback_props->HII_min_density;
  const char is_dense = rho >= rho_threshold;

  /* Can the particle be ionized? */
  return (is_cold && is_dense && !radiation_is_part_tagged_as_ionized(p, xp));
}

/**
 * @brief Photoionization rate coefficient this star delivers at a gas
 * particle's location, frozen at tag time.
 *
 * The cooling task cannot recompute it later (it has no neighbour search),
 * so it is stored on the particle. The intermediate photon flux would
 * overflow float32 in this unit system, so only the final coefficient --
 * computed in double up to that point -- is returned. Zero unless
 * GEARFeedback:HII_couple_ionization_rate is on.
 *
 * @param si The #spart (star) providing photons.
 * @param r2 Squared comoving distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param us Internal unit system.
 * @param cosmo The current cosmological model.
 * @param cooling Cooling function data.
 */
__attribute__((always_inline)) INLINE static float
feedback_hii_photoionization_rate_HI(
    const struct spart *restrict si, float r2, int pixel,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling) {

  if (!cooling->HII_couple_ionization_rate) return 0.f;

  const float Omega_pixel =
      4.0f * (float)M_PI / si->feedback_data.radiation.n_HII_pixels;

  /* r2 is comoving, so the physical 1/r^2 dilution carries an a^-2. Floored at
     a fraction of the star's own smoothing length: a gas particle sitting on
     top of the star would otherwise divide by zero and hand Grackle an
     infinite heating rate. */
  const float r2_min = 1e-6f * si->h * si->h;
  const double r2_phys = max(r2, r2_min) * cosmo->a * cosmo->a;
  const double ionizing_flux_HI =
      feedback_get_star_ionization_rate(si, pixel) / (Omega_pixel * r2_phys);
  return (float)radiation_get_photoionization_rate_coefficient_from_flux_HI(
      us, ionizing_flux_HI);
}

/**
 * @brief Absolute time until which an ionization tag stamped now stays valid.
 *
 * Sized on the interval the pass just covered, so a tag outlives the gap to
 * its star's next rebuild and cooling cannot expire it first.
 *
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 */
__attribute__((always_inline)) INLINE static double feedback_hii_tag_end_time(
    const double time, const double dt_back) {
  return time + RADIATION_TAG_LIFETIME_INTERVALS * dt_back;
}

/**
 * @brief Tag a gas particle as ionized by this star and charge its photons.
 *
 * Shared by the three acceptance branches of #feedback_iact_HII_ionization.
 * No atomics: task_type_stars_hii_ionization_feedback's cell locking
 * (src/task.c) already serializes every writer of these fields.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being ionized.
 * @param xpj The #xpart (gas) being ionized.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param us Internal unit system.
 * @param cosmo The current cosmological model.
 * @param cooling Cooling function data.
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 * @param cost Ionizing photon count to charge for this particle.
 */
__attribute__((always_inline)) INLINE static void feedback_hii_claim_part(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const double time,
    const double dt_back, const double cost) {

  if (xpj->tracers_data.HII_region.is_ionized != 0) return;

  radiation_tag_part_as_ionized(
      pj, xpj, si->id, feedback_hii_tag_end_time(time, dt_back),
      si->feedback_data.radiation.mean_excess_photon_energy_HI,
      feedback_hii_photoionization_rate_HI(si, r2, pixel, us, cosmo, cooling));
  timestep_sync_part(pj);

  radiation_consume_ionizing_photons(si, pixel, cost);

  /* Update HII region properties */
  si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
  si->h_hii = max(si->h_hii, sqrtf(r2) * kernel_gamma_inv);
}

/**
 * @brief Charge an already-ionized gas particle's recombination losses and
 * renew its tag.
 *
 * Recombinations must be replaced continuously to hold gas ionized, so a
 * particle already held ionized costs photons every pass -- charging only
 * newly-claimed ones is what let the ionized volume grow without bound as the
 * rebuild cadence was refined. There is no one-off N_H term here: those
 * electrons are already stripped.
 *
 * Renewal is budget-gated: once a pixel is spent, its remaining gas is not
 * renewed and cooling expires those tags at end_time. Note this is charged in
 * traversal order (own cell in array order, then the radiation_in links), NOT
 * in radius order as Phase 2's distance-sorted buffer is, so the set that
 * lapses when a star cannot afford its whole region is spatially arbitrary
 * rather than an outer shell. At equilibrium the shortfall is only the
 * marginal budget deficit, so that set is thin.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being maintained.
 * @param xpj The #xpart (gas) being maintained.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param phys_const Physics constants.
 * @param hydro_props Hydrodynamics properties.
 * @param us Internal unit system.
 * @param cosmo Cosmology.
 * @param cooling Cooling function data.
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 * @return Photons charged for this particle, 0 if the budget was exhausted.
 */
__attribute__((always_inline)) INLINE double
feedback_iact_HII_maintain_ionized_part(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const double time,
    const double dt_back) {

  /* Counted whether or not its photons can be afforded: this field is the
     mass the star HOLDS ionized, matching h_hii's extent (see stars_io.h). */
  si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);

  if (feedback_get_star_ionization_budget(si, pixel) <= 0.0) return 0.0;

  const double cost =
      dt_back * radiation_get_part_rate_to_fully_ionize(
                    phys_const, hydro_props, us, cosmo, cooling, pj, xpj);
  radiation_consume_ionizing_photons(si, pixel, cost);

  /* Renew in place: the tag's other fields would otherwise keep values
     frozen at first claim, growing staler every pass the particle survives.
     Deliberately no timestep_sync_part() here, unlike a fresh claim: waking
     every held particle every pass would drag the whole region down to the
     shortest time bin. Expiry is therefore checked on the gas particle's own
     next cooling task, so recession can lag by up to one of its timesteps. */
  radiation_tag_part_as_ionized(
      pj, xpj, si->id, feedback_hii_tag_end_time(time, dt_back),
      si->feedback_data.radiation.mean_excess_photon_energy_HI,
      feedback_hii_photoionization_rate_HI(si, r2, pixel, us, cosmo, cooling));

  return cost;
}

/**
 * @brief Perform the ionization of a gas particle by a star.
 *
 * @param si The #spart (star) providing photons.
 * @param pj The #part (gas) being ionized.
 * @param xpj The #xpart (gas) being ionized.
 * @param r2 Squared distance between star and gas.
 * @param pixel The angular pixel this gas particle was assigned to.
 * @param phys_const Physics constants.
 * @param hydro_props Hydrodynamics properties.
 * @param us Internal unit system.
 * @param cosmo Cosmology.
 * @param cooling Cooling function data.
 * @param feedback_props The #feedback_props.
 * @param ti_begin Integer time at the start of the step (for RNG).
 * @param time The current simulation time.
 * @param dt_back Time elapsed since this star's previous HII rebuild pass.
 */
__attribute__((always_inline)) INLINE void feedback_iact_HII_ionization(
    struct spart *restrict si, struct part *restrict pj,
    struct xpart *restrict xpj, float r2, int pixel,
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling,
    const struct feedback_props *feedback_props, const integertime_t ti_begin,
    const double time, const double dt_back) {

  const int deterministic_boundary =
      feedback_props->HII_deterministic_boundary_ionization;

  /* If already ionized (by another thread, or by an earlier pass), just move
     on: recombination losses for gas already held ionized are charged by
     feedback_iact_HII_maintain_ionized_part, not here. */
  if (radiation_is_part_tagged_as_ionized(pj, xpj)) return;

  /* Photons this candidate costs: a one-off payment to strip its remaining
     NEUTRAL hydrogen (not its total hydrogen content -- a particle whose
     tag lapsed on a marginal budget shortfall, or one pre-ionized by a UV
     background, is already partway or fully stripped, and re-paying full
     N_H on reclaim would be a cadence-coupled photon sink), plus a
     maintenance reserve sized by the *elapsed* interval (the next
     interval's recombinations are charged by
     feedback_iact_HII_maintain_ionized_part, not here). That reserve is what
     makes region growth implicit -- dS/dt = (Q-S)/t_rec integrates as
     dS = (Q-S)*dt/(t_rec+dt) -- and so unconditionally stable at the
     dt >> t_rec the default HII_rebuild_time_Myr produces. Dropping it would
     give explicit Euler, which overshoots and then churns. */
  const double N_HI = radiation_get_part_number_neutral_hydrogen_atoms(
      phys_const, hydro_props, us, cosmo, cooling, pj, xpj);
  const double Delta_dot_N_ion_maintenance_rate =
      radiation_get_part_rate_to_fully_ionize(phys_const, hydro_props, us,
                                              cosmo, cooling, pj, xpj);
  const double cost = N_HI + dt_back * Delta_dot_N_ion_maintenance_rate;

  const double budget = feedback_get_star_ionization_budget(si, pixel);

  /* Case 1: Ionization is guaranteed */
  if (cost <= budget) {
    feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                            dt_back, cost);
  } else if (deterministic_boundary) {
    /* Deterministic mode: always ionize the boundary particle, letting the
       pixel's budget go slightly negative (bounded by one particle's
       ionization cost). */
    feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                            dt_back, cost);
  } else {
    /* Probabilistic mode: a weighted coin flip decides whether to fully
       ionize pj. On a win, the full cost is consumed (more than what was
       left, going slightly negative); on a loss, nothing is consumed, so
       the budget survives intact to try the next, farther candidate. This
       keeps the expected photon consumption equal to what's actually
       available: proba*cost + (1-proba)*0 = budget. */
    const float proba = budget / cost;
    /* Keyed on the star and the pixel, deliberately *not* on pj: this must be
       one trial per pixel per pass for the identity below to hold. The same
       number is drawn for every candidate offered to this pixel, so a loss
       rejects the whole remaining shell -- which is exactly the (1 - proba)
       branch. Rolling per candidate instead would give each of the up-to
       HII_max_retry_full_buffer x max_ngbs candidates an independent shot in a
       loop that stops at the first win, so a pass would claim one particle
       almost surely and spend `cost` rather than `budget`. */
    const float random_number = random_unit_interval_part_ID_and_index(
        si->id, pixel, ti_begin, random_number_HII_regions);

    if (random_number <= proba) {
      /* We won the roll! Claim the particle. */
      feedback_hii_claim_part(si, pj, xpj, r2, pixel, us, cosmo, cooling, time,
                              dt_back, cost);
    }
    /* Lost the roll: consume nothing, budget carries over. */
  } /* End of probability handling */
}

/**
 * @brief Number of active angular pixels this star's ionizing photon budget
 * is currently split across (1 = spherical, HEALPix disabled).
 *
 * Exists so callers outside this feedback model (e.g. the generic
 * runner_radiation_feedback.c) never need to reach directly into
 * sp->feedback_data.radiation, a GEAR-specific struct.
 *
 * @param sp The #spart to query.
 */
int feedback_get_star_HII_pixel_count(const struct spart *sp) {
  return sp->feedback_data.radiation.n_HII_pixels;
}

/**
 * @brief Record the (star-relative) age at which this star's HII region was
 * last rebuilt.
 *
 * @param sp The #spart to update.
 * @param star_age_beg_step The star's age at the start of the current step.
 */
void feedback_set_star_HII_last_rebuild(struct spart *sp,
                                        double star_age_beg_step) {
  sp->feedback_data.radiation.HII_region_last_rebuild = star_age_beg_step;
}

/**
 * @brief The (star-relative) age at which this star's HII region was last
 * rebuilt.
 *
 * @param sp The #spart to query.
 */
double feedback_get_star_HII_last_rebuild(const struct spart *sp) {
  return sp->feedback_data.radiation.HII_region_last_rebuild;
}

/**
 * @brief Record the (star-relative) age at which this star was last given a
 * chance to rebuild its HII region, whether or not that chance found gas.
 *
 * @param sp The #spart to update.
 * @param star_age_beg_step The star's age at the start of the current step.
 */
void feedback_set_star_HII_last_attempt(struct spart *sp,
                                        double star_age_beg_step) {
  sp->feedback_data.radiation.HII_region_last_attempt = star_age_beg_step;
}

/**
 * @brief The (star-relative) age at which this star was last given a chance
 * to rebuild its HII region, whether or not that chance found gas.
 *
 * @param sp The #spart to query.
 */
double feedback_get_star_HII_last_attempt(const struct spart *sp) {
  return sp->feedback_data.radiation.HII_region_last_attempt;
}

/**
 * @brief The HII rebuild cadence currently in force for this star, i.e. the
 * interval one scheduled pass is expected to cover.
 *
 * Used to bound the elapsed interval the photon budget is integrated over
 * (#HII_DT_BACK_MAX_INTERVALS), so passes skipped by a gas-free cell cannot
 * bank photons that escaped instead of being stored.
 *
 * @param feedback_props The #feedback_props.
 * @param dt_enrichment This star's current enrichment timestep, the stand-in
 * cadence in "rebuild every step" mode.
 */
double feedback_get_star_HII_nominal_interval(
    const struct feedback_props *feedback_props, const double dt_enrichment) {

  const double floor_interval = (double)feedback_props->HII_rebuild_floor_Myr;

  if (feedback_props->HII_rebuild_time > 0.f)
    return (double)feedback_props->HII_rebuild_time;

  /* Negative rebuild time means "every step", so the star's own timestep is
     the cadence. */
  return max(dt_enrichment, floor_interval);
}

/**
 * @brief Dispatch wrapper exposing radiation_open_ionizing_photon_budget()
 * (radiation.c) to runner_radiation_feedback.c, same reasoning as
 * #feedback_get_star_HII_pixel_count.
 *
 * @param sp The star.
 * @param dt_back Time elapsed since this star's last HII rebuild pass.
 */
void feedback_open_star_ionizing_photon_budget(struct spart *sp,
                                               double dt_back) {
  radiation_open_ionizing_photon_budget(sp, dt_back);
}

/**
 * @brief Is this gas particle currently tagged as HII-ionized?
 *
 * Thin dispatch wrapper so callers outside this feedback model (e.g.
 * star_formation/GEAR, sink/GEAR, both of which are selectable
 * independently of the feedback model) can query ionization state without
 * depending on this model being the one actually compiled in -- every
 * feedback model provides this function, matching #radiation_is_part_
 * tagged_as_ionized() here for GEAR and unconditionally returning false
 * everywhere else.
 *
 * @param p The #part to query.
 * @param xp The #part's extended data.
 */
char feedback_is_part_tagged_as_ionized(const struct part *p,
                                        const struct xpart *xp) {
  return radiation_is_part_tagged_as_ionized(p, xp);
}

/**
 * @brief Current ionized mass of this star's HII region.
 *
 * Dispatch wrapper so callers outside this feedback model (e.g. the GEAR
 * stars I/O converter, which is compiled whenever stars=GEAR regardless of
 * the feedback model) can read this without depending on GEAR feedback
 * being the one actually compiled in.
 *
 * @param sp The #spart to query.
 */
float feedback_get_star_HII_mass(const struct spart *sp) {
  return sp->feedback_data.radiation.mass_HII_region;
}

/**
 * @brief Seed value for the stars.h_hii_max cell tracker (radiation task
 * placement only -- see radiation_get_h_hii_max_seed()).
 *
 * Age/eligibility-gated on top of that helper: outside the HII-active
 * window (dead, or aged past HII_max_age) the seed collapses to sp->h_hii
 * like the real region does once nothing is held any more, instead of
 * pinning the tracker at the seeded reach forever after the star stops
 * doing any ionization. Deliberately NOT gated on will_do_HII_ionization
 * (feedback_is_HII_ionization_active): that flag is true only on the
 * specific step a rebuild pass actually runs, and gating on it would make
 * the seed flicker between rebuild-cadence steps instead of tracking the
 * star's age smoothly.
 *
 * Dispatch wrapper so callers outside this feedback model (the generic
 * cell/space rebuild bookkeeping that collects stars.h_hii_max, guarded
 * there by #IONIZATION_FEEDBACK_LOOP) can reach this without including
 * GEAR-internal headers.
 *
 * @param sp The #spart to query.
 * @param e The #engine (for physical constants, units, cooling, current
 * time/cosmology, and GEARFeedback:HII_max_age_Myr).
 */
float feedback_get_star_h_hii_max_seed(const struct spart *sp,
                                       const struct engine *e) {

  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return sp->h_hii;

  const int with_cosmology = e->policy & engine_policy_cosmology;
  const double star_age =
      compute_star_age_end_of_step(sp, with_cosmology, e->cosmology, e->time);
  if (star_age > e->feedback_props->HII_max_age) return sp->h_hii;

  return radiation_get_h_hii_max_seed(sp, e->physical_constants,
                                      e->internal_units, e->cooling_func);
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

  /* The star has nothing useful to do in this loop yet. Order of operations
     for a newly-formed star:
     1. Star formation: age_beg_step < 0, age_end_step = 0.
     2. Stars density/prep/feedback-apply tasks this step: nothing to
        distribute (will_do_feedback = 0 below).
     3. Timestep task: feedback_will_do_feedback() runs the stellar
        evolution, to be distributed starting next step. */
  sp->feedback_data.will_do_feedback = 0;
  sp->feedback_data.will_do_HII_ionization = 0;

  /* 0.0 reads as "never attempted", so the first attempt's dt_elapsed
     measures from the star's actual birth age. */
  sp->feedback_data.radiation.HII_region_last_attempt = 0.0;

  /* A newborn star owes no photons: radiation_open_ionizing_photon_budget()
     carries a pixel's overdraft into the next pass, so this must start clean.
     Over the whole array, not n_HII_pixels, which can still grow. */
  for (int p = 0; p < HII_MAX_ANGULAR_PIXELS; p++) {
    sp->feedback_data.radiation.N_ion_budget_pix[p] = 0.0;
  }

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

  /* n_HII_pixels bounds the loop in feedback_get_star_ionization_budget_max();
     set it before the first stellar_evolution call so a population star
     that returns early (already past the IMF's alive range) never reads
     it uninitialized. */
  sp->feedback_data.radiation.n_HII_pixels = 1;

  /* Same reasoning as the identical seed in
     feedback_init_after_star_formation(). */
  sp->feedback_data.radiation.HII_region_last_attempt = 0.0;

  /* No photon debt carried in from before the run, same reasoning as the
     identical seed in feedback_init_after_star_formation(). */
  for (int p = 0; p < HII_MAX_ANGULAR_PIXELS; p++) {
    sp->feedback_data.radiation.N_ion_budget_pix[p] = 0.0;
  }

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
 * @param us The unit system.
 * @param phys_const The physical constants in internal units.
 */
void feedback_struct_restore(struct feedback_props *feedback, FILE *stream,
                             const struct unit_system *us,
                             const struct phys_const *phys_const) {

  restart_read_blocks((void *)feedback, sizeof(struct feedback_props), 1,
                      stream, NULL, "feedback function");

  stellar_evolution_restore(&feedback->stellar_model, stream,
                            feedback->with_stellar_wind_feedback, us,
                            phys_const);

  if (feedback->metallicity_max_first_stars != -1) {
    stellar_evolution_restore(&feedback->stellar_model_first_stars, stream,
                              feedback->with_stellar_wind_feedback, us,
                              phys_const);
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
