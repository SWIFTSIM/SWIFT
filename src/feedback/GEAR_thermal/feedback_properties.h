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
#ifndef SWIFT_GEAR_FEEDBACK_PROPERTIES_H
#define SWIFT_GEAR_FEEDBACK_PROPERTIES_H

#include "../GEAR/stellar_evolution.h"
#include "../GEAR/stellar_evolution_struct.h"
#include "chemistry.h"
#include "hydro_properties.h"

#define default_HII_min_density_Hpcm3 1.0
#define default_HII_max_age_Myr 50.0
#define default_HII_rebuild_time_Myr 0.5
#define default_HII_deterministic_boundary_ionization 0
#define default_HII_rebuild_floor_Myr 1e-4
#define default_dt_evolution_factor_max 300.0
#define default_event_dt_floor_Myr 1e-4

/**
 * @brief The different subgrid radiation feedback processes GEAR models.
 */
enum radiation_policy {
  radiation_policy_none = 0,
  /*! Do we want the ionization effect (Strömgren sphere)? */
  radiation_policy_photoionization = (1 << 0),
  /*! Radiation pressure from the stars' bolometric luminosity */
  radiation_policy_radiation_pressure = (1 << 1),
  /* Photoelectric (PE) heating by FUV radiation on dust */
  radiation_policy_photoelectric_heating = (1 << 2),
};

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {

  /*! Supernovae energy effectively deposited */
  float supernovae_efficiency;

  /* ------------- Stellar model properties ------------- */

  /*! The stellar model */
  struct stellar_model stellar_model;

  /*! The stellar model for first stars */
  struct stellar_model stellar_model_first_stars;

  /*! Metallicity limits for the first stars */
  float metallicity_max_first_stars;

  /*! Metallicity [Fe/H] transition for the first stars */
  float imf_transition_metallicity;

  /* ------------- Star evolution timestep properties ------------- */

  /*! Timestep refinement factor as lifetime_myr -> 0, used by
   * feedback_compute_spart_timestep()'s logistic transition (SSP particles
   * only; single_star particles use the exact, zero-tuning-parameter
   * dt_event instead -- see event_dt_floor_Myr below). The logistic's
   * midpoint and steepness are fixed internal constants
   * (GEAR_dt_evolution_lifetime_myr_0/GEAR_dt_evolution_steepness in
   * feedback_common.c), not parsed from params.yml. */
  float dt_evolution_factor_max;

  /*! Floors the event-anchored timestep terms (single_star's dt_event,
   * and dt_HII_safe for every star type) in internal units, before they are
   * combined with the coarser min_star_timestep-floored non-event terms
   * (stars_compute_timestep(), src/stars/GEAR/stars.h). Runs for every
   * star regardless of radiation, so parsed unconditionally -- unlike
   * HII_rebuild_floor_Myr, which stays 0.0 (no floor at all) whenever
   * with_photoionization is off, this constant must never be zero, since
   * dt_event applies to every single_star particle unconditionally. */
  float event_dt_floor_Myr;

  /* ------------- Subgrid Radiation properties ------------- */

  /* The radiation processes enabled */
  int radiation_policy;

  /*! Radiation pressure momentum effectively injected */
  float radiation_pressure_efficiency;

  /*! Minimal density to consider a particle eligible for HII ionization */
  float HII_min_density;

  /*! HII region rebuild frequency */
  float HII_rebuild_time;

  /*! Maximun age of star particle to trigger the HII region algorithm */
  float HII_max_age;

  /*! How to treat the boundary gas particle a star's remaining photon
   * budget cannot fully ionize: 0 = probabilistic (weighted coin flip,
   * unbiased in expectation), 1 = deterministic (always ionize it,
   * letting the budget go slightly negative). */
  char HII_deterministic_boundary_ionization;

  /*! Floors the elapsed interval the per-pass ionizing photon budget is
   * integrated over, in every cadence mode. */
  float HII_rebuild_floor_Myr;

  /* ------------- Stellar winds properties ------------- */

  /*! Pre-supernova feedback energy effectively deposited */
  float winds_efficiency;

  /*! Do stellar wind feedback? */
  char with_stellar_wind_feedback;
};

/**
 * @brief Print the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void feedback_props_print(
    const struct feedback_props *feedback_props) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  /* Print the name of the elements */
  char txt[GEAR_CHEMISTRY_ELEMENT_COUNT * (GEAR_LABELS_SIZE + 2)] = "";
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (i != 0) {
      strcat(txt, ", ");
    }
    strcat(txt, stellar_evolution_get_element_name(
                    &feedback_props->stellar_model, i));
  }

  if (engine_rank == 0) {
    message("Chemistry elements: %s", txt);
  }

  /* Print the feedback properties */
  message("Supernovae efficiency                                      = %.2g",
          feedback_props->supernovae_efficiency);
  message("Stellar wind feedback                                      = %s",
          feedback_props->with_stellar_wind_feedback ? "ON" : "OFF");
  message("Stellar winds efficiency                                   = %.2g",
          feedback_props->winds_efficiency);
  message("dt_evolution factor_max                                    = %g",
          feedback_props->dt_evolution_factor_max);
  message("event_dt_floor (internal units)                            = %g",
          feedback_props->event_dt_floor_Myr);

  const char do_photoionization =
      feedback_props->radiation_policy & radiation_policy_photoionization;
  message("Photoionization                                            = %i",
          do_photoionization);

  if (do_photoionization) {
    message("HII region minimal gas density (internal units)            = %g",
            feedback_props->HII_min_density);
    message("HII boundary ionization mode                               = %s",
            feedback_props->HII_deterministic_boundary_ionization
                ? "deterministic"
                : "probabilistic");
    message("HII max age (internal units)                               = %g",
            feedback_props->HII_max_age);
    message("HII rebuild time (internal units)                          = %g",
            feedback_props->HII_rebuild_time);
    message("HII rebuild floor (internal units)                         = %g",
            feedback_props->HII_rebuild_floor_Myr);
  }

  message(
      "Radiation pressure                                         = %i",
      feedback_props->radiation_policy & radiation_policy_radiation_pressure);
  message("Radiation pressure efficiency                              = %.2g",
          feedback_props->radiation_pressure_efficiency);
  message("Photo-electric heating                                     = %i",
          feedback_props->radiation_policy &
              radiation_policy_photoelectric_heating);

  message("Yields table                                               = %s",
          feedback_props->stellar_model.yields_table);

  /* Print the stellar model */
  stellar_model_print(&feedback_props->stellar_model);

  /* Print the first stars */
  if (feedback_props->metallicity_max_first_stars != -1) {
    message("Yields table first stars                                 = %s",
            feedback_props->stellar_model_first_stars.yields_table);
    stellar_model_print(&feedback_props->stellar_model_first_stars);
    message("Metallicity max for the first stars (in abundance)       = %g",
            feedback_props->imf_transition_metallicity);
    message("Metallicity max for the first stars (in mass fraction)   = %g",
            feedback_props->metallicity_max_first_stars);
  }
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
__attribute__((always_inline)) INLINE static void feedback_props_init(
    struct feedback_props *fp, const struct phys_const *phys_const,
    const struct unit_system *us, struct swift_params *params,
    const struct hydro_props *hydro_props, const struct cosmology *cosmo) {

  /* Supernovae energy efficiency */
  double e_efficiency =
      parser_get_param_double(params, "GEARFeedback:supernovae_efficiency");
  fp->supernovae_efficiency = e_efficiency;

  /* Activate the stellar wind feedback */
  char with_stellar_wind_feedback = (char)parser_get_param_int(
      params, "GEARFeedback:with_stellar_wind_feedback");
  fp->with_stellar_wind_feedback = with_stellar_wind_feedback;

  /* Are we running with photoionization? Read early, like
   * with_stellar_wind_feedback above, so stellar_evolution_props_init() can
   * skip opening the radiation table for non-radiation runs. */
  const char with_photoionization = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_photoionization", 0);

  /* Radiation pressure efficiency, read early for the same reason: its sign
   * decides whether the radiation table is needed too, independently of
   * with_photoionization. */
  const float radiation_pressure_efficiency = parser_get_opt_param_float(
      params, "GEARFeedback:radiation_pressure_efficiency", 0.0);

  /* The radiation table backs both the HII photoionization band and the
   * bolometric radiation-pressure band (see
   * stellar_evolution_compute_preSN_feedback_individual_star()/_spart());
   * photoelectric heating has no downstream consumer yet, so it does not
   * gate this. */
  const char with_radiation =
      with_photoionization || (radiation_pressure_efficiency > 0.0f);

  /* Pre-Supernovae energy efficiency */
  double w_efficiency = 0.0;
  if (with_stellar_wind_feedback) {
    w_efficiency = parser_get_param_double(
        params, "GEARFeedback:stellar_winds_efficiency");
  }

  fp->winds_efficiency = w_efficiency;

  /* filename of the chemistry tables. */
  parser_get_param_string(params, "GEARFeedback:yields_table",
                          fp->stellar_model.yields_table);

  /* Initialize the stellar models. */
  stellar_evolution_props_init(&fp->stellar_model, phys_const, us, params,
                               cosmo, fp->with_stellar_wind_feedback,
                               with_radiation);

  /* Read the metallicity threashold */
  fp->imf_transition_metallicity = parser_get_opt_param_float(
      params, "GEARFeedback:imf_transition_metallicity", 0);

  /* Read and get the solar abundances */
  struct chemistry_global_data data;
  bzero(&data, sizeof(struct chemistry_global_data));
  chemistry_read_solar_abundances(params, &data);

  const int iFe = stellar_evolution_get_element_index(&fp->stellar_model, "Fe");
  const float XFe = data.solar_abundances[iFe];

  if (fp->imf_transition_metallicity == 0)
    fp->metallicity_max_first_stars = -1;
  else
    fp->metallicity_max_first_stars =
        exp10(fp->imf_transition_metallicity) * XFe;

  /* Now initialize the first stars. */
  if (fp->metallicity_max_first_stars == -1) {
    message("First stars are disabled.");
  } else {
    if (fp->metallicity_max_first_stars < 0) {
      error(
          "The metallicity threshold for the first stars is in mass fraction. "
          "It cannot be lower than 0.");
    }
    if (engine_rank == 0) {
      message("Reading the stellar model for the first stars");
    }
    parser_get_param_string(params, "GEARFeedback:yields_table_first_stars",
                            fp->stellar_model_first_stars.yields_table);
    stellar_evolution_props_init(&fp->stellar_model_first_stars, phys_const, us,
                                 params, cosmo, fp->with_stellar_wind_feedback,
                                 with_radiation);
  }

  /* ------------- Star evolution timestep properties ------------- */
  /* Runs for every star regardless of radiation, so parsed unconditionally
   * (feedback_compute_spart_timestep() uses this factor for all stars). */

  /* Needed unconditionally below (event_dt_floor_Myr's conversion), unlike
   * HII_max_age/HII_rebuild_time/HII_rebuild_floor_Myr further down, which
   * stay gated behind with_photoionization. */
  const double Myr_internal_units = 1e6 * phys_const->const_year;

  fp->dt_evolution_factor_max =
      parser_get_opt_param_float(params, "GEARFeedback:dt_evolution_factor_max",
                                 default_dt_evolution_factor_max);

  if (fp->dt_evolution_factor_max < 1.f)
    error("GEARFeedback:dt_evolution_factor_max must be >= 1 (got %g).",
          fp->dt_evolution_factor_max);

  fp->event_dt_floor_Myr = parser_get_opt_param_float(
      params, "GEARFeedback:event_dt_floor_Myr", default_event_dt_floor_Myr);

  if (fp->event_dt_floor_Myr <= 0.f)
    error(
        "GEARFeedback:event_dt_floor_Myr must be > 0 (got %g): it floors "
        "single_star's event-anchored death timestep (dt_event) in every "
        "run, not just radiation ones -- <= 0 reopens the get_spart_timestep "
        "crash this parameter exists to prevent.",
        fp->event_dt_floor_Myr);

  fp->event_dt_floor_Myr *= Myr_internal_units;

  /* Startup cross-check: event_dt_floor_Myr is the sole crash guard for
   * every single_star death (get_spart_timestep()'s dt_min error()), so it
   * must itself clear dt_min. TimeIntegration:dt_min is a mandatory
   * parameter already fully parsed into `params` at this point (engine_init
   * has not run yet, but that only matters for e->dt_min the struct field --
   * the parsed value is available directly from `params` regardless of
   * init order). */
  const double dt_min =
      parser_get_param_double(params, "TimeIntegration:dt_min");
  if (fp->event_dt_floor_Myr <= dt_min)
    error(
        "GEARFeedback:event_dt_floor_Myr (%g, internal units) must exceed "
        "TimeIntegration:dt_min (%g, internal units): otherwise the floor "
        "itself can still trip get_spart_timestep()'s dt_min error() at a "
        "single_star's death.",
        fp->event_dt_floor_Myr, dt_min);

  /* ------------- Subgrid Radiation properties ------------- */
  fp->radiation_policy = 0;

  /* TODO: For the future, enforce these to have a non-zero value */

  /* Radiation pressure */
  fp->radiation_pressure_efficiency = radiation_pressure_efficiency;

  if (fp->radiation_pressure_efficiency > 0.0) {
    fp->radiation_policy |= radiation_policy_radiation_pressure;
  }

  const int with_photoelectric_heating = parser_get_opt_param_int(
      params, "GEARFeedback:with_photoelectric_heating", 0);

  if (with_photoelectric_heating) {
    fp->radiation_policy |= radiation_policy_photoelectric_heating;
  }

  if (with_photoionization) {
    fp->radiation_policy |= radiation_policy_photoionization;

    /* Read the minimal density */
    fp->HII_min_density =
        parser_get_opt_param_float(params, "GEARFeedback:HII_min_density_Hpcm3",
                                   default_HII_min_density_Hpcm3);

    /* Read the HII region maximal age */
    fp->HII_max_age = parser_get_opt_param_float(
        params, "GEARFeedback:HII_max_age_Myr", default_HII_max_age_Myr);

    /* Read the HII region rebuild frequency */
    fp->HII_rebuild_time =
        parser_get_opt_param_float(params, "GEARFeedback:HII_rebuild_time_Myr",
                                   default_HII_rebuild_time_Myr);

    /* Read the boundary-particle ionization mode */
    fp->HII_deterministic_boundary_ionization = parser_get_opt_param_int(
        params, "GEARFeedback:HII_deterministic_boundary_ionization",
        default_HII_deterministic_boundary_ionization);

    fp->HII_rebuild_floor_Myr =
        parser_get_opt_param_float(params, "GEARFeedback:HII_rebuild_floor_Myr",
                                   default_HII_rebuild_floor_Myr);
    if (fp->HII_rebuild_floor_Myr <= 0.f)
      error(
          "GEARFeedback:HII_rebuild_floor_Myr must be > 0 (got %g): it "
          "floors the interval every per-pass ionizing photon budget is "
          "integrated over, in every cadence mode -- <= 0 silently zeroes "
          "every star's first-pass budget.",
          fp->HII_rebuild_floor_Myr);

    /* Convert to internal units */
    const double m_p_cgs = phys_const->const_proton_mass *
                           units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    fp->HII_min_density *=
        m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

    /* Myr_internal_units already computed above, unconditionally. */
    fp->HII_max_age *= Myr_internal_units;
    fp->HII_rebuild_time *= Myr_internal_units;
    fp->HII_rebuild_floor_Myr *= Myr_internal_units;
  }

  /* -------------------------------------------- */
  /* Print the stellar properties */
  feedback_props_print(fp);

  /* Print a final message. */
  if (engine_rank == 0) {
    message("Stellar feedback initialized");
  }
}

#endif /* SWIFT_GEAR_FEEDBACK_PROPERTIES_H */
