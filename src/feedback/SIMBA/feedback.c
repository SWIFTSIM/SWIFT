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

/* This file's header */
#include "feedback.h"

/* Local includes. */
#include "../EAGLE/enrichment.h"
#include "../EAGLE/yield_tables.h"
#include "hydro_properties.h"
#include "inline.h"
#include "random.h"
#include "timers.h"


/**
 * @brief Determine if the star forming gas should kick a particle,
 *        and then kick it.
 *
 * @param p The #part to consider.
 * @param xp The #xpart to consider.
 * @param e The #engine.
 * @param fb_props The feedback properties.
 * @param ti_current The current timestep.
 */
double feedback_wind_probability(struct part* p, struct xpart* xp, 
                                 const struct engine* e, 
                                 const struct cosmology* cosmo,
                                 const struct feedback_props* fb_props, 
                                 const integertime_t ti_current, 
                                 const double dt_part,
                                 double *rand_for_sf_wind) {

  /* First thing we will do is generate a random number */
  *rand_for_sf_wind = random_unit_interval(p->id, ti_current,
                                           random_number_stellar_feedback_1);

  double galaxy_stellar_mass = p->gpart->fof_data.group_stellar_mass;
  if (galaxy_stellar_mass <= 0.) return 0.;

  const double stellar_mass_this_step = xp->sf_data.SFR * dt_part;
  if (stellar_mass_this_step <= 0.) return 0.;

  /* If M* is non-zero, make sure it is at least resolved in the
   * following calculations.
   */
  if (galaxy_stellar_mass < fb_props->minimum_galaxy_stellar_mass) {
    galaxy_stellar_mass = fb_props->minimum_galaxy_stellar_mass;
  }

  /* When early wind suppression is enabled, we alter the minimum
   * stellar mass to be safe.
   */
  if (fb_props->early_wind_suppression_enabled) {
    const double early_minimum_stellar_mass =
        fb_props->early_stellar_mass_norm *
        exp(
          -1. *
          (
            (cosmo->a / fb_props->early_wind_suppression_scale_factor) *
            (cosmo->a / fb_props->early_wind_suppression_scale_factor)
          )
        );
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      galaxy_stellar_mass = early_minimum_stellar_mass;
    }
  }

  double wind_mass = 
      fb_props->FIRE_eta_normalization * stellar_mass_this_step;
  if (galaxy_stellar_mass < fb_props->FIRE_eta_break) {
    wind_mass *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_lower_slope /*-0.317*/
    );
  } else {
    wind_mass *= pow(
      galaxy_stellar_mass / fb_props->FIRE_eta_break, 
      fb_props->FIRE_eta_upper_slope /*-0.761*/
    );
  }

  /* Suppress stellar feedback in the early universe when galaxies are
   * too small. Star formation can destroy unresolved galaxies, so
   * we must suppress the stellar feedback.
   */
  if (fb_props->early_wind_suppression_enabled) {
    if (cosmo->a < fb_props->early_wind_suppression_scale_factor) {
      wind_mass *= pow(cosmo->a / fb_props->early_wind_suppression_scale_factor, 
                       fb_props->early_wind_suppression_slope);
    }
  }

  const float probability_to_kick = 1. - exp(-wind_mass / hydro_get_mass(p));
  return probability_to_kick;

}

void feedback_kick_and_decouple_part(struct part* p, struct xpart* xp, 
                                     const struct engine* e, 
                                     const struct cosmology* cosmo,
                                     const struct feedback_props* fb_props, 
                                     const integertime_t ti_current,
                                     const int with_cosmology,
                                     const double dt_part) {

  const double galaxy_stellar_mass = 
      p->gpart->fof_data.group_stellar_mass;
  const double galaxy_gas_stellar_mass_Msun = 
      p->gpart->fof_data.group_mass * fb_props->mass_to_solar_mass;
  if (galaxy_gas_stellar_mass_Msun <= 0. || galaxy_stellar_mass <= 0.) return;

  /* Physical circular velocity km/s */
  const double v_circ_km_s = 
      pow(galaxy_gas_stellar_mass_Msun / 102.329, 0.26178) *
      pow(cosmo->H / cosmo->H0, 1. / 3.);
  const double rand_for_scatter = random_unit_interval(p->id, ti_current,
                                      random_number_stellar_feedback_2);

  /* physical km/s */
  const double wind_velocity =
      fb_props->FIRE_velocity_normalization *
      pow(v_circ_km_s / 200., fb_props->FIRE_velocity_slope) *
      (
        1. - fb_props->kick_velocity_scatter + 
        2. * fb_props->kick_velocity_scatter * rand_for_scatter
      ) *
      v_circ_km_s *
      fb_props->kms_to_internal;

  const double dir[3] = {
    p->gpart->a_grav[1] * p->gpart->v_full[2] - 
        p->gpart->a_grav[2] * p->gpart->v_full[1],
    p->gpart->a_grav[2] * p->gpart->v_full[0] - 
        p->gpart->a_grav[0] * p->gpart->v_full[2],
    p->gpart->a_grav[0] * p->gpart->v_full[1] - 
        p->gpart->a_grav[1] * p->gpart->v_full[0]
  };
  const double norm = sqrt(
    dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
  );
  const double prefactor = cosmo->a * wind_velocity / norm;

  xp->v_full[0] += prefactor * dir[0];
  xp->v_full[1] += prefactor * dir[1];
  xp->v_full[2] += prefactor * dir[2];

  /* Decouple the particles from the hydrodynamics */
  p->feedback_data.decoupling_delay_time = 
      fb_props->wind_decouple_time_factor * 
      cosmology_get_time_since_big_bang(cosmo, cosmo->a);

  p->feedback_data.number_of_times_decoupled += 1;

  /* Wind cannot be star forming */
  if (xp->sf_data.SFR > 0.f) {

    /* Record the current time as an indicator of when this particle was last
       star-forming. */
    if (with_cosmology) {
      xp->sf_data.SFR = -e->cosmology->a;
    } else {
      xp->sf_data.SFR = -e->time;
    }

  }

  /**
   * z pid dt M* Mb vkick vkx vky vkz h x y z vx vy vz T rho v_sig decoupletime 
   * Ndecouple
   */
  const float length_convert = cosmo->a * fb_props->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
  const float u_convert = 
      cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;
  printf("WIND_LOG %.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",
          cosmo->z,
          p->id, 
          dt_part * fb_props->time_to_Myr,
          galaxy_stellar_mass * fb_props->mass_to_solar_mass,
          galaxy_gas_stellar_mass_Msun,
          wind_velocity / fb_props->kms_to_internal,
          prefactor * dir[0] * velocity_convert,
          prefactor * dir[1] * velocity_convert,
          prefactor * dir[2] * velocity_convert,
          p->h * cosmo->a * fb_props->length_to_kpc,
          p->x[0] * length_convert, 
          p->x[1] * length_convert, 
          p->x[2] * length_convert,
          p->gpart->v_full[0] * velocity_convert, 
          p->gpart->v_full[1] * velocity_convert, 
          p->gpart->v_full[2] * velocity_convert,
          p->u * u_convert, 
          p->rho * rho_convert, 
          p->viscosity.v_sig * velocity_convert,
          p->feedback_data.decoupling_delay_time * fb_props->time_to_Myr, 
          p->feedback_data.number_of_times_decoupled);
}

/**
 * @brief calculates stellar mass in spart that died over the timestep, calls
 * functions to calculate feedback due to SNIa, SNII and AGB
 *
 * @param feedback_props feedback_props data structure
 * @param phys_const The physical constants in internal units.
 * @param cosmo The cosmological model.
 * @param sp spart that we're evolving
 * @param us unit_system data structure
 * @param age age of spart at beginning of step
 * @param dt length of current timestep
 * @param ti_begin The current integer time (for random number hashing).
 */
void compute_stellar_evolution(const struct feedback_props* feedback_props,
                               const struct phys_const* phys_const,
                               const struct cosmology* cosmo, struct spart* sp,
                               const struct unit_system* us, const double age,
                               const double dt, const integertime_t ti_begin) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (age < 0.f) error("Negative age for a star.");

  if (sp->count_since_last_enrichment != 0 && engine_current_step > 0)
    error("Computing feedback on a star that should not");
#endif

  /* Convert dt and stellar age from internal units to Gyr. */
  const double Gyr_inv = 1. / (phys_const->const_year * 1e9);
  const double dt_Gyr = dt * Gyr_inv;
  const double star_age_Gyr = age * Gyr_inv;

  /* Get the birth mass of the star */
  const double M_init = sp->mass_init;

  /* Get the total metallicity (metal mass fraction) at birth time and impose a
   * minimum */
  const double Z =
      max(chemistry_get_star_total_metal_mass_fraction_for_feedback(sp),
          exp10(log10_min_metallicity));

  /* Get the individual abundances (mass fractions at birth time) */
  const float* const abundances =
      chemistry_get_star_metal_mass_fraction_for_feedback(sp);

  /* Properties collected in the stellar density loop. */
  const float ngb_gas_mass = sp->feedback_data.to_collect.ngb_mass;

  /* Check if there are neighbours, otherwise exit */
  if (ngb_gas_mass == 0.f || sp->density.wcount * pow_dimension(sp->h) < 1e-4) {
    feedback_reset_feedback(sp, feedback_props);
    return;
  }

  /* Update the enrichment weights */
  const float enrichment_weight_inv =
      sp->feedback_data.to_collect.enrichment_weight_inv;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_collect.enrichment_weight_inv < 0.)
    error("Negative inverse weight!");
#endif

  /* Now we start filling the data structure for information to apply to the
   * particles. Do _NOT_ read from the to_collect substructure any more. */

  /* Zero all the output fields */
  feedback_reset_feedback(sp, feedback_props);

  /* Update the weights used for distribution */
  const float enrichment_weight =
      (enrichment_weight_inv != 0.f) ? 1.f / enrichment_weight_inv : 0.f;
  sp->feedback_data.to_distribute.enrichment_weight = enrichment_weight;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.enrichment_weight < 0.)
    error("Negative weight!");
#endif

  /* Calculate mass of stars that has died from the star's birth up to the
   * beginning and end of timestep */
  const double max_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr, Z, feedback_props);
  const double min_dying_mass_Msun =
      dying_mass_msun(star_age_Gyr + dt_Gyr, Z, feedback_props);

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity check. Worth investigating if necessary as functions for evaluating
   * mass of stars dying might be strictly decreasing.  */
  if (min_dying_mass_Msun > max_dying_mass_Msun)
    error("min dying mass is greater than max dying mass");
#endif

  /* Integration interval is zero - this can happen if minimum and maximum
   * dying masses are above imf_max_mass_Msun. Return without doing any
   * enrichment. */
  if (min_dying_mass_Msun == max_dying_mass_Msun) return;

  /* Life is better in log */
  const double log10_max_dying_mass_Msun = log10(max_dying_mass_Msun);
  const double log10_min_dying_mass_Msun = log10(min_dying_mass_Msun);

  /* Compute elements, energy and momentum to distribute from the
   *  three channels SNIa, SNII, AGB */
  if (feedback_props->with_SNIa_enrichment) {
    evolve_SNIa(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                feedback_props, star_age_Gyr, dt_Gyr, &sp->feedback_data);
  }
  if (feedback_props->with_SNII_enrichment) {
    evolve_SNII(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
                abundances, feedback_props, &sp->feedback_data);
  }
  if (feedback_props->with_AGB_enrichment) {
    evolve_AGB(log10_min_dying_mass_Msun, log10_max_dying_mass_Msun, M_init, Z,
               abundances, feedback_props, &sp->feedback_data);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->feedback_data.to_distribute.mass != 0.f)
    error("Injected mass will be lost");
#endif

  /* Compute the total mass to distribute (H + He  metals) */
  sp->feedback_data.to_distribute.mass =
      sp->feedback_data.to_distribute.total_metal_mass +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_H] +
      sp->feedback_data.to_distribute.metal_mass[chemistry_element_He];

  /* Compute energy change due to kinetic energy of ejectas */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass *
      feedback_props->AGB_ejecta_specific_kinetic_energy;

  /* Compute energy change due to kinetic energy of the star */
  sp->feedback_data.to_distribute.energy +=
      sp->feedback_data.to_distribute.mass * 0.5f *
      (sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2]) *
      cosmo->a2_inv;

  TIMER_TOC(timer_do_star_evol);
}

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
void feedback_props_init(struct feedback_props* fp,
                         const struct phys_const* phys_const,
                         const struct unit_system* us,
                         struct swift_params* params,
                         const struct hydro_props* hydro_props,
                         const struct cosmology* cosmo) {

  /* Gather common conversion factors --------------------------------------- */

  /* Calculate internal mass to solar mass conversion factor */
  const double Msun_cgs = phys_const->const_solar_mass *
                          units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  const double unit_mass_cgs = units_cgs_conversion_factor(us, UNIT_CONV_MASS);
  fp->mass_to_solar_mass = unit_mass_cgs / Msun_cgs;
  fp->solar_mass_to_mass = 1. / fp->mass_to_solar_mass;

  /* Calculate temperature to internal energy conversion factor (all internal
   * units) */
  const double k_B = phys_const->const_boltzmann_k;
  const double m_p = phys_const->const_proton_mass;
  const double mu = hydro_props->mu_ionised;
  fp->temp_to_u_factor = k_B / (mu * hydro_gamma_minus_one * m_p);

  /* Calculate conversion factor from rho to n_H
   * Note this assumes primoridal abundance */
  const double X_H = hydro_props->hydrogen_mass_fraction;
  fp->rho_to_n_cgs =
      (X_H / m_p) * units_cgs_conversion_factor(us, UNIT_CONV_NUMBER_DENSITY);

  fp->kms_to_internal = 1.0e5f / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  fp->time_to_Myr = units_cgs_conversion_factor(us, UNIT_CONV_TIME) /
      (1.e6f * 365.25f * 24.f * 60.f * 60.f);

  fp->length_to_kpc = 
      units_cgs_conversion_factor(us, UNIT_CONV_LENGTH) / 3.08567758e21f;

  /* Main operation modes ------------------------------------------------- */

  fp->with_SNII_feedback =
      parser_get_param_int(params, "SIMBAFeedback:use_SNII_feedback");

  fp->with_SNIa_feedback =
      parser_get_param_int(params, "SIMBAFeedback:use_SNIa_feedback");

  fp->with_AGB_enrichment =
      parser_get_param_int(params, "SIMBAFeedback:use_AGB_enrichment");

  fp->with_SNII_enrichment =
      parser_get_param_int(params, "SIMBAFeedback:use_SNII_enrichment");

  fp->with_SNIa_enrichment =
      parser_get_param_int(params, "SIMBAFeedback:use_SNIa_enrichment");

  if (fp->with_SNIa_feedback && !fp->with_SNIa_enrichment) {
    error("Cannot run with SNIa feedback without SNIa enrichment.");
  }

  /* Properties of the IMF model ------------------------------------------ */

  /* Minimal and maximal mass considered */
  fp->imf_max_mass_msun =
      parser_get_param_double(params, "SIMBAFeedback:IMF_max_mass_Msun");
  fp->imf_min_mass_msun =
      parser_get_param_double(params, "SIMBAFeedback:IMF_min_mass_Msun");

  /* Check that it makes sense. */
  if (fp->imf_max_mass_msun < fp->imf_min_mass_msun) {
    error("Can't have the max IMF mass smaller than the min IMF mass!");
  }

  fp->log10_imf_max_mass_msun = log10(fp->imf_max_mass_msun);
  fp->log10_imf_min_mass_msun = log10(fp->imf_min_mass_msun);


  /* Properties of the SNII model */

  /* Energy released by supernova type II */
  fp->E_SNII_cgs =
      parser_get_param_double(params, "SIMBAFeedback:SNII_energy_erg");
  fp->E_SNII =
      fp->E_SNII_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Stellar mass limits for SNII feedback */
  const double SNII_min_mass_msun =
      parser_get_param_double(params, "SIMBAFeedback:SNII_min_mass_Msun");
  const double SNII_max_mass_msun =
      parser_get_param_double(params, "SIMBAFeedback:SNII_max_mass_Msun");

  /* Check that it makes sense. */
  if (SNII_max_mass_msun < SNII_min_mass_msun) {
    error("Can't have the max SNII mass smaller than the min SNII mass!");
  }

  fp->SNII_min_mass_msun = SNII_min_mass_msun;
  fp->SNII_max_mass_msun = SNII_max_mass_msun;
  fp->log10_SNII_min_mass_msun = log10(SNII_min_mass_msun);
  fp->log10_SNII_max_mass_msun = log10(SNII_max_mass_msun);


  /* Properties of the SNII enrichment model -------------------------------- */

  /* Set factors for each element adjusting SNII yield */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "SIMBAFeedback:SNII_yield_factor_%s",
            chemistry_get_element_name((enum chemistry_element)elem));

    fp->SNII_yield_factor[elem] =
        parser_get_opt_param_float(params, buffer, 1.f);
  }

  /* Properties of the SNIa enrichment model -------------------------------- */

  fp->SNIa_DTD_delay_Gyr =
      parser_get_param_double(params, "SIMBAFeedback:SNIa_DTD_delay_Gyr");

  char temp[32] = {0};
  parser_get_param_string(params, "SIMBAFeedback:SNIa_DTD", temp);

  if (strcmp(temp, "Exponential") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_exponential;

    /* Read SNIa exponential DTD model parameters */
    fp->SNIa_DTD_exp_norm = parser_get_param_float(
        params, "SIMBAFeedback:SNIa_DTD_exp_norm_p_Msun");
    fp->SNIa_DTD_exp_timescale_Gyr = parser_get_param_float(
        params, "SIMBAFeedback:SNIa_DTD_exp_timescale_Gyr");
    fp->SNIa_DTD_exp_timescale_Gyr_inv = 1.f / fp->SNIa_DTD_exp_timescale_Gyr;

  } else if (strcmp(temp, "PowerLaw") == 0) {

    fp->SNIa_DTD = eagle_feedback_SNIa_DTD_power_law;

    /* Read SNIa power-law DTD model parameters */
    fp->SNIa_DTD_power_law_norm = parser_get_param_float(
        params, "SIMBAFeedback:SNIa_DTD_power_law_norm_p_Msun");

    /* Renormalize everything such that the integral converges to
       'SNIa_DTD_power_law_norm' over 13.6 Gyr. */
    fp->SNIa_DTD_power_law_norm /= log(13.6 / fp->SNIa_DTD_delay_Gyr);

  } else {
    error("Invalid SNIa DTD model: '%s'", temp);
  }

  /* Energy released by supernova type Ia */
  fp->E_SNIa_cgs =
      parser_get_param_double(params, "SIMBAFeedback:SNIa_energy_erg");
  fp->E_SNIa =
      fp->E_SNIa_cgs / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Properties of the SNIa enrichment model -------------------------------- */

  /* Read AGB ejecta velocity */
  const float ejecta_velocity_km_p_s = parser_get_param_float(
      params, "SIMBAFeedback:AGB_ejecta_velocity_km_p_s");

  /* Convert to internal units */
  const float ejecta_velocity_cgs = ejecta_velocity_km_p_s * 1e5;
  const float ejecta_velocity =
      ejecta_velocity_cgs / units_cgs_conversion_factor(us, UNIT_CONV_SPEED);

  /* Convert to specific thermal energy */
  fp->AGB_ejecta_specific_kinetic_energy =
      0.5f * ejecta_velocity * ejecta_velocity;

  /* Properties of the enrichment down-sampling ----------------------------- */

  fp->stellar_evolution_age_cut =
      parser_get_param_double(params,
                              "SIMBAFeedback:stellar_evolution_age_cut_Gyr") *
      phys_const->const_year * 1e9;

  fp->stellar_evolution_sampling_rate = parser_get_param_double(
      params, "SIMBAFeedback:stellar_evolution_sampling_rate");

  if (fp->stellar_evolution_sampling_rate < 1 ||
      fp->stellar_evolution_sampling_rate >= (1 << (8 * sizeof(char) - 1)))
    error("Stellar evolution sampling rate too large. Must be >0 and <%d",
          (1 << (8 * sizeof(char) - 1)));

  /* Properties of Simba kinetic winds -------------------------------------- */

  fp->FIRE_velocity_normalization =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_velocity_normalization");
  fp->FIRE_velocity_slope =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_velocity_slope");
  fp->FIRE_eta_normalization =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_eta_normalization");
  fp->FIRE_eta_break =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_eta_break_Msun");
  fp->FIRE_eta_break *= fp->solar_mass_to_mass;
  fp->FIRE_eta_lower_slope =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_eta_lower_slope");
  fp->FIRE_eta_upper_slope =
      parser_get_param_double(params, "SIMBAFeedback:FIRE_eta_upper_slope");

  fp->early_stellar_mass_norm =
      parser_get_param_double(params, "SIMBAFeedback:early_stellar_mass_norm_Msun");
  fp->early_stellar_mass_norm *= fp->solar_mass_to_mass;
  fp->early_wind_suppression_enabled = 
      parser_get_param_int(params, "SIMBAFeedback:early_wind_suppression_enabled");
  fp->early_wind_suppression_scale_factor =
      parser_get_param_double(params, 
          "SIMBAFeedback:early_wind_suppression_scale_factor");
  fp->early_wind_suppression_slope =
      parser_get_param_double(params, "SIMBAFeedback:early_wind_suppression_slope");

  fp->minimum_galaxy_stellar_mass =
      parser_get_param_double(params, "SIMBAFeedback:minimum_galaxy_stellar_mass_Msun");
  fp->minimum_galaxy_stellar_mass *= fp->solar_mass_to_mass;

  fp->kick_velocity_scatter =
      parser_get_param_double(params, "SIMBAFeedback:kick_velocity_scatter");

  fp->wind_decouple_time_factor = parser_get_param_double(
      params, "SIMBAFeedback:wind_decouple_time_factor");

  /* Initialise the IMF ------------------------------------------------- */

  init_imf(fp);

  /* Calculate number of type II SN per unit solar mass based on our choice
   * of IMF and integration limits for type II SNe.
   * Note: No weighting by yields here. */
  fp->num_SNII_per_msun =
      integrate_imf(fp->log10_SNII_min_mass_msun, fp->log10_SNII_max_mass_msun,
                    eagle_imf_integration_no_weight,
                    /*(stellar_yields=)*/ NULL, fp);

  /* Initialise the yields ---------------------------------------------- */

  /* Read yield table filepath  */
  parser_get_param_string(params, "SIMBAFeedback:filename",
                          fp->yield_table_path);

  /* Allocate yield tables  */
  allocate_yield_tables(fp);

  /* Read the tables  */
  read_yield_tables(fp);

  /* Set yield_mass_bins array */
  const float imf_log10_mass_bin_size =
      (fp->log10_imf_max_mass_msun - fp->log10_imf_min_mass_msun) /
      (eagle_feedback_N_imf_bins - 1);

  for (int i = 0; i < eagle_feedback_N_imf_bins; i++)
    fp->yield_mass_bins[i] =
        imf_log10_mass_bin_size * i + fp->log10_imf_min_mass_msun;

  /* Resample yields from mass bins used in tables to mass bins used in IMF  */
  compute_yields(fp);

  /* Resample ejecta contribution to enrichment from mass bins used in tables to
   * mass bins used in IMF  */
  compute_ejecta(fp);

  if (engine_rank == 0) {
    message("Feedback model is SIMBA");
    message("Feedback FIRE velocity normalization: %g", 
            fp->FIRE_velocity_normalization);
    message("Feedback FIRE velocity slope: %g", fp->FIRE_velocity_slope);
    message("Feedback velocity scatter: %g", fp->kick_velocity_scatter);
    message("Feedback FIRE eta normalization: %g", fp->FIRE_eta_normalization);
    message("Feedback FIRE eta break: %g", fp->FIRE_eta_break);
    message("Feedback FIRE eta upper slope: %g", fp->FIRE_eta_upper_slope);
    message("Feedback FIRE eta lower slope: %g", fp->FIRE_eta_lower_slope);

    message("Feedback early suppression enabled: %d", 
            fp->early_wind_suppression_enabled);
    message("Feedback early stellar mass norm: %g Msun", 
            fp->early_stellar_mass_norm);
    message("Feedback early suppression scale factor: %g", 
            fp->early_wind_suppression_scale_factor);
    message("Feedback early suppression slope: %g", fp->early_wind_suppression_slope);

  }
}

/**
 * @brief Clean-up the memory allocated for the feedback routines
 *
 * We simply free all the arrays.
 *
 * @param fp the feedback data structure.
 */
void feedback_clean(struct feedback_props* fp) {

  swift_free("imf-tables", fp->imf);
  swift_free("imf-tables", fp->imf_mass_bin);
  swift_free("imf-tables", fp->imf_mass_bin_log10);
  swift_free("feedback-tables", fp->yields_SNIa);
  swift_free("feedback-tables", fp->yield_SNIa_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.mass);
  swift_free("feedback-tables", fp->yield_AGB.metallicity);
  swift_free("feedback-tables", fp->yield_AGB.yield);
  swift_free("feedback-tables", fp->yield_AGB.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.ejecta);
  swift_free("feedback-tables", fp->yield_AGB.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_AGB.total_metals);
  swift_free("feedback-tables", fp->yield_AGB.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.mass);
  swift_free("feedback-tables", fp->yield_SNII.metallicity);
  swift_free("feedback-tables", fp->yield_SNII.yield);
  swift_free("feedback-tables", fp->yield_SNII.yield_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.ejecta);
  swift_free("feedback-tables", fp->yield_SNII.ejecta_IMF_resampled);
  swift_free("feedback-tables", fp->yield_SNII.total_metals);
  swift_free("feedback-tables", fp->yield_SNII.total_metals_IMF_resampled);
  swift_free("feedback-tables", fp->lifetimes.mass);
  swift_free("feedback-tables", fp->lifetimes.metallicity);
  swift_free("feedback-tables", fp->yield_mass_bins);
  for (int i = 0; i < eagle_feedback_lifetime_N_metals; i++) {
    free(fp->lifetimes.dyingtime[i]);
  }
  free(fp->lifetimes.dyingtime);
  for (int i = 0; i < eagle_feedback_SNIa_N_elements; i++) {
    free(fp->SNIa_element_names[i]);
  }
  free(fp->SNIa_element_names);
  for (int i = 0; i < eagle_feedback_SNII_N_elements; i++) {
    free(fp->SNII_element_names[i]);
  }
  free(fp->SNII_element_names);
  for (int i = 0; i < eagle_feedback_AGB_N_elements; i++) {
    free(fp->AGB_element_names[i]);
  }
  free(fp->AGB_element_names);
}

/**
 * @brief Write a feedback struct to the given FILE as a stream of bytes.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_dump(const struct feedback_props* feedback, FILE* stream) {

  /* To make sure everything is restored correctly, we zero all the pointers to
     tables. If they are not restored correctly, we would crash after restart on
     the first call to the feedback routines. Helps debugging. */
  struct feedback_props feedback_copy = *feedback;

  /* zero AGB and SNII table pointers */
  zero_yield_table_pointers(&feedback_copy.yield_AGB);
  zero_yield_table_pointers(&feedback_copy.yield_SNII);

  /* zero SNIa table pointers */
  feedback_copy.yield_SNIa_IMF_resampled = NULL;
  feedback_copy.yields_SNIa = NULL;
  feedback_copy.yield_SNIa_total_metals_IMF_resampled = 0;

  /* zero element name tables */
  feedback_copy.SNIa_element_names = NULL;
  feedback_copy.SNII_element_names = NULL;
  feedback_copy.AGB_element_names = NULL;

  /* zero mass bins table */
  feedback_copy.yield_mass_bins = NULL;

  /* zero lifetime tracks */
  feedback_copy.lifetimes.mass = NULL;
  feedback_copy.lifetimes.metallicity = NULL;
  feedback_copy.lifetimes.dyingtime = NULL;

  /* zero IMF tables */
  feedback_copy.imf = NULL;
  feedback_copy.imf_mass_bin = NULL;
  feedback_copy.imf_mass_bin_log10 = NULL;

  restart_write_blocks((void*)&feedback_copy, sizeof(struct feedback_props), 1,
                       stream, "feedback", "feedback function");
}

/**
 * @brief Restore a hydro_props struct from the given FILE as a stream of
 * bytes.
 *
 * Read the structure from the stream and restore the feedback tables by
 * re-reading them.
 *
 * @param feedback the struct
 * @param stream the file stream
 */
void feedback_struct_restore(struct feedback_props* feedback, FILE* stream) {
  restart_read_blocks((void*)feedback, sizeof(struct feedback_props), 1, stream,
                      NULL, "feedback function");

  if (strlen(feedback->yield_table_path) != 0)
    feedback_restore_tables(feedback);
}
