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

#define default_HII_region_min_density_Hpcm3 1.0
#define default_HII_region_max_age_Myr 50.0
#define default_HII_region_rebuild_time_Myr 0.5

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

  /* ------------- Subgrid Radiation properties ------------- */

  /* The radiation processes enabled */
  int radiation_policy;

  /*! Radiation pressure momentum effectively injected */
  float radiation_pressure_efficiency;

  /*! Minimal density to consider a particle eligible for HII ionization */
  float HII_region_min_density;

  /*! HII region rebuild frequency */
  float HII_region_rebuild_time;

  /*! Maximun age of star particle to trigger the HII region algorithm */
  float HII_region_max_age;

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

  const char do_photoionization =
      feedback_props->radiation_policy & radiation_policy_photoionization;
  message("Photoionization                                            = %i",
          do_photoionization);

  if (do_photoionization)
    message("HII region minimal gas density (internal units)            = %g",
            feedback_props->HII_region_min_density);

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
                               cosmo, fp->with_stellar_wind_feedback);

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
                                 params, cosmo, fp->with_stellar_wind_feedback);
  }

  /* ------------- Subgrid Radiation properties ------------- */
  fp->radiation_policy = 0;

  /* TODO: For the future, enforce these to have a non-zero value */

  /* Radiation pressure */
  fp->radiation_pressure_efficiency = parser_get_opt_param_float(
      params, "GEARFeedback:radiation_pressure_efficiency", 0.0);

  if (fp->radiation_pressure_efficiency > 0.0) {
    fp->radiation_policy |= radiation_policy_radiation_pressure;
  }

  const int with_photoelectric_heating = parser_get_opt_param_int(
      params, "GEARFeedback:with_photoelectric_heating", 0);

  if (with_photoelectric_heating) {
    fp->radiation_policy |= radiation_policy_photoelectric_heating;
  }

  /* Are we running with photoionization? */
  const int with_photoionization =
      parser_get_opt_param_int(params, "GEARFeedback:with_photoionization", 0);

  if (with_photoionization) {
    fp->radiation_policy |= radiation_policy_photoionization;
    
    /* Read the minimal density */
    fp->HII_region_min_density = parser_get_opt_param_float(
        params, "GEARFeedback:HII_region_min_density_Hpcm3",
        default_HII_region_min_density_Hpcm3);

    /* Read the HII region maximal age */
    fp->HII_region_max_age = parser_get_opt_param_float(
        params, "GEARFeedback:HII_region_max_age_Myr",
        default_HII_region_max_age_Myr);

    /* Read the HII region rebuild frequency */
    fp->HII_region_rebuild_time = parser_get_opt_param_float(
        params, "GEARFeedback:HII_region_rebuild_time_Myr",
        default_HII_region_rebuild_time_Myr);

    /* Convert to internal units */
    const double m_p_cgs = phys_const->const_proton_mass *
                           units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    fp->HII_region_min_density *=
        m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

    const double Myr_internal_units = 1e6 * phys_const->const_year;
    fp->HII_region_max_age *= Myr_internal_units;
    fp->HII_region_rebuild_time *= Myr_internal_units;

    if (fp->HII_region_rebuild_time <= 0.0) {
      /* TODO: What do we do? We rebuild at every step the star is active */
    }
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
