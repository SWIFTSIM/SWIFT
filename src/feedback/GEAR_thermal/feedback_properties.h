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

#include "../GEAR/radiation_isrf.h"
#include "../GEAR/stellar_evolution.h"
#include "../GEAR/stellar_evolution_struct.h"
#include "chemistry.h"
#include "hydro_properties.h"
#include "minmax.h"

#include <math.h>
#include <string.h>

#define default_HII_min_density_Hpcm3 1.0
#define default_HII_max_age_Myr 50.0
#define default_HII_rebuild_time_Myr 0.5
#define default_HII_deterministic_boundary_ionization 0
#define default_HII_rebuild_floor_Myr 1e-4
#define default_dt_evolution_factor_max 300.0
#define default_event_dt_floor_Myr 1e-4
/* One kernel support radius: the largest path with the star inside the
 * receiver's kernel. */
#define default_ISRF_extinction_path_in_kernel_radii 1.0

/* Temperature cap in K of the temperature_capped_jeans extinction path:
 * Safranek-Shrader et al. (2017), MNRAS 465, 885, Section 3.4. */
#define default_ISRF_extinction_jeans_temperature_cap_K 40.0

/**
 * @brief The subgrid radiation feedback processes.
 */
enum radiation_policy {
  radiation_policy_none = 0,
  /*! Photoionization (Strömgren sphere). */
  radiation_policy_photoionization = (1 << 0),
  /*! Radiation pressure from the stars' bolometric luminosity */
  radiation_policy_radiation_pressure = (1 << 1),

  /*! Interstellar radiation field: photoelectric heating and H2
   * photodissociation. */
  radiation_policy_isrf = (1 << 2),
};

/**
 * @brief Scheme that sets the ISRF propagation speed `c_hyp_i`.
 *
 * The parameter file uses the names #isrf_c_hyp_scheme_name_fixed_fraction and
 * #isrf_c_hyp_scheme_name_kernel_local. The integer values are internal and
 * go to the restart file.
 */
enum isrf_c_hyp_scheme {
  /*! `c_hyp_i = ISRF_c_hyp_fixed_fraction_of_c * c` for every particle. A
   * timestep term enforces the CFL condition. Parameter value
   * `fixed_fraction`. */
  isrf_c_hyp_scheme_fixed_fraction = 2,
  /*! `c_hyp_i = min(C_hyp*h_i/dt_max(i), c)`, `dt_max(i)` the longest
   * timestep in the kernel of i. Default. Parameter value `kernel_local`. */
  isrf_c_hyp_scheme_kernel_local_reduced_flux = 4,
};

/*! Value of GEARFeedback:ISRF_c_hyp_scheme that selects
 * #isrf_c_hyp_scheme_fixed_fraction. */
#define isrf_c_hyp_scheme_name_fixed_fraction "fixed_fraction"

/*! Value of GEARFeedback:ISRF_c_hyp_scheme that selects
 * #isrf_c_hyp_scheme_kernel_local_reduced_flux. */
#define isrf_c_hyp_scheme_name_kernel_local "kernel_local"

/**
 * @brief Convert a GEARFeedback:ISRF_c_hyp_scheme value to #isrf_c_hyp_scheme.
 *
 * @param name The parameter value.
 *
 * @return The scheme, or -1 if @p name is neither
 * #isrf_c_hyp_scheme_name_fixed_fraction nor
 * #isrf_c_hyp_scheme_name_kernel_local.
 */
__attribute__((always_inline)) INLINE static int
feedback_props_c_hyp_scheme_from_name(const char *name) {
  if (strcmp(name, isrf_c_hyp_scheme_name_fixed_fraction) == 0)
    return isrf_c_hyp_scheme_fixed_fraction;
  if (strcmp(name, isrf_c_hyp_scheme_name_kernel_local) == 0)
    return isrf_c_hyp_scheme_kernel_local_reduced_flux;
  return -1;
}

/**
 * @brief Convert a #isrf_c_hyp_scheme to its GEARFeedback:ISRF_c_hyp_scheme
 * value.
 *
 * @param scheme The scheme.
 *
 * @return The parameter value, or "invalid" if @p scheme is not a valid
 * #isrf_c_hyp_scheme.
 */
__attribute__((always_inline)) INLINE static const char *
feedback_props_c_hyp_scheme_name(int scheme) {
  if (scheme == isrf_c_hyp_scheme_fixed_fraction)
    return isrf_c_hyp_scheme_name_fixed_fraction;
  if (scheme == isrf_c_hyp_scheme_kernel_local_reduced_flux)
    return isrf_c_hyp_scheme_name_kernel_local;
  return "invalid";
}

/**
 * @brief Dust extinction path of the LW/PE bands
 * (GEARFeedback:ISRF_extinction_path).
 */
enum isrf_extinction_path_mechanism {
  /*! `l = R * kernel_gamma * h_j`, with `R` =
   * #feedback_props.ISRF_extinction_path_in_kernel_radii. */
  isrf_extinction_path_constant_kernel_path = 0,
  /*! `l = r`, the star-to-particle separation of the pair. */
  isrf_extinction_path_pair_separation = 1,
  /*! `l = min(lambda_J(min(T, T_cap)), kernel_gamma * h_j)`, with `T_cap` =
   * #feedback_props.ISRF_extinction_jeans_temperature_cap_K. */
  isrf_extinction_path_temperature_capped_jeans = 2,
};

/**
 * @brief Properties of the GEAR feedback model.
 */
struct feedback_props {

  /*! Whether sinks are configured; set by sink_props_init(), 0 otherwise. */
  int with_sinks;

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

  /*! Timestep refinement factor of SSP stars as lifetime_myr -> 0. */
  float dt_evolution_factor_max;

  /*! Floor on the event-anchored star timestep terms (dt_event,
   * dt_HII_safe), in internal units after init. Never zero. */
  float event_dt_floor_Myr;

  /* ------------- Subgrid Radiation properties ------------- */

  /* The radiation processes enabled */
  int radiation_policy;

  /*! Radiation pressure momentum effectively injected */
  float radiation_pressure_efficiency;

  /*! Run the hyperbolic M1 propagation update? Only meaningful when
   * radiation_policy_isrf is set. */
  char ISRF_propagation;

  /*! Band-edge weights lambda_E(PE), lambda_E(LW), lambda_N(LW). Never 0 or
   * 1. */
  double band_edge_weight_pe;
  double band_edge_weight_lw;
  double band_edge_photon_weight_lw;

  /*! Active #isrf_extinction_path_mechanism. */
  char ISRF_extinction_path_mechanism;

  /*! Path R of #isrf_extinction_path_constant_kernel_path, in kernel support
   * radii. Only parsed for that mechanism. */
  float ISRF_extinction_path_in_kernel_radii;

  /*! Temperature cap in K of #isrf_extinction_path_temperature_capped_jeans. */
  float ISRF_extinction_jeans_temperature_cap_K;

  /*! Stability margin C_hyp in `c_hyp_i = C_hyp*h_i/dt_max(i)`, see
   * feedback_props_init() for its range. */
  float ISRF_c_hyp_margin;

  /*! Active #isrf_c_hyp_scheme, parsed from a string and stored as an
   * integer. 0 (not a valid scheme) when the ISRF is off. */
  int ISRF_c_hyp_scheme;

  /*! Fraction f with `c_hyp_i = f*c`, fixed_fraction scheme only. 0 when
   * unused. */
  float ISRF_c_hyp_fixed_fraction_of_c;

  /*! Debug only: skip the timestep term `C_hyp*h_i/(f*c)` of the
   * fixed_fraction scheme. */
  char ISRF_c_hyp_timestep_term_off_for_debugging;

  /*! Ceiling of the negativity-triggered dissipation coefficient. 0 disables
   * it. The allowed range depends on #ISRF_c_hyp_margin. */
  float ISRF_dissipation_alpha_max;

  /*! Relative undershoot of `rho_prev*u` below the neighbours' kernel mean
   * at which the trigger reaches #ISRF_dissipation_alpha_max. */
  float ISRF_dissipation_negativity_threshold;

  /*! Floor `alpha_floor/(1+(h*kappa/eps_lambda)^4)` under the trigger,
   * combined with it by a max. 0 disables it. */
  float ISRF_dissipation_alpha_floor;

  /*! Screening-length budget `eps_lambda` of the floor roll-off. */
  float ISRF_dissipation_floor_h_over_lambda;

  /*! Threshold `eps_R` in [0, 1] of the flux-relaxation residual gate on
   * the floor. 0 disables it. */
  float ISRF_dissipation_floor_relaxation_residual;

  /*! Debug only: when positive, pin every dissipation coefficient to this
   * value. */
  float ISRF_dissipation_alpha_pin_for_debugging;

  /*! Minimal density to consider a particle eligible for HII ionization */
  float HII_min_density;

  /*! HII region rebuild frequency */
  float HII_rebuild_time;

  /*! Maximun age of star particle to trigger the HII region algorithm */
  float HII_max_age;

  /*! Boundary particle the photon budget cannot fully ionize: 0 =
   * probabilistic, 1 = always ionize it. */
  char HII_deterministic_boundary_ionization;

  /*! Floor on the interval the photon budget is integrated over. */
  float HII_rebuild_floor_Myr;

  /* ------------- Stellar winds properties ------------- */

  /*! Pre-supernova feedback energy effectively deposited */
  float winds_efficiency;

  /*! Do stellar wind feedback? */
  char with_stellar_wind_feedback;
};

/**
 * @brief Does this run need cooling_init() to have run, so that the ISRF
 * reads a real local_dust_to_gas_ratio?
 *
 * @param feedback_props The #feedback_props.
 * @return True if the ISRF is enabled.
 */
__attribute__((always_inline)) INLINE static int
feedback_props_needs_cooling_initialized(
    const struct feedback_props *feedback_props) {
  return (feedback_props->radiation_policy & radiation_policy_isrf) != 0;
}

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

  /* Grouped like GEARFeedback in parameter_example.yml. */

  /* Stellar evolution */
  message("Yields table                                               = %s",
          feedback_props->stellar_model.yields_table);
  message("dt_evolution factor_max                                    = %g",
          feedback_props->dt_evolution_factor_max);
  message("event_dt_floor (internal units)                            = %g",
          feedback_props->event_dt_floor_Myr);

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

  /* Supernovae */
  message("Supernovae efficiency                                      = %.2g",
          feedback_props->supernovae_efficiency);

  /* Stellar winds */
  message("Stellar wind feedback                                      = %s",
          feedback_props->with_stellar_wind_feedback ? "ON" : "OFF");
  message("Stellar winds efficiency                                   = %.2g",
          feedback_props->winds_efficiency);

  /* Radiation pressure */
  message(
      "Radiation pressure                                         = %i",
      feedback_props->radiation_policy & radiation_policy_radiation_pressure);
  message("Radiation pressure efficiency                              = %.2g",
          feedback_props->radiation_pressure_efficiency);

  /* HII regions */
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

  /* ISRF */
  const char do_photoelectric_heating =
      feedback_props->radiation_policy & radiation_policy_isrf;
  message("Photo-electric heating / H2 photodissociation (ISRF)       = %i",
          do_photoelectric_heating);
  if (do_photoelectric_heating) {
    message("ISRF propagation                                           = %s",
            feedback_props->ISRF_propagation ? "ON" : "OFF (injection only)");
    if (feedback_props->ISRF_propagation) {
      message("ISRF propagation speed margin (C_hyp)                      = %g",
              feedback_props->ISRF_c_hyp_margin);
      message(
          "ISRF c_hyp scheme                                          = %s",
          feedback_props_c_hyp_scheme_name(feedback_props->ISRF_c_hyp_scheme));
      if (feedback_props->ISRF_c_hyp_fixed_fraction_of_c > 0.f) {
        message(
            "ISRF c_hyp fixed fraction of c                             = %g",
            feedback_props->ISRF_c_hyp_fixed_fraction_of_c);
        message(
            "ISRF c_hyp fixed fraction timestep term                    = %s",
            feedback_props->ISRF_c_hyp_timestep_term_off_for_debugging
                ? "OFF for debugging (diagnostic only, not a supported "
                  "configuration)"
                : "ON");
      }
      message("ISRF dissipation alpha_max                                 = %g",
              feedback_props->ISRF_dissipation_alpha_max);
      message("ISRF dissipation negativity threshold                      = %g",
              feedback_props->ISRF_dissipation_negativity_threshold);
      message("ISRF dissipation alpha_floor                               = %g",
              feedback_props->ISRF_dissipation_alpha_floor);
      message("ISRF dissipation floor h/lambda budget (eps_lambda)        = %g",
              feedback_props->ISRF_dissipation_floor_h_over_lambda);
      message("ISRF dissipation floor relaxation-residual gate (eps_R)    = %g",
              feedback_props->ISRF_dissipation_floor_relaxation_residual);
      if (feedback_props->ISRF_dissipation_alpha_pin_for_debugging > 0.f)
        message(
            "ISRF dissipation alpha pinned for debugging                = %g",
            feedback_props->ISRF_dissipation_alpha_pin_for_debugging);
    }
  }
}

/**
 * @brief Enforce that #feedback_props.ISRF_c_hyp_scheme and
 * #feedback_props.ISRF_c_hyp_fixed_fraction_of_c are a matched pair.
 *
 * The two schemes are alternatives. Standalone so a unit test can call it.
 *
 * @param scheme #feedback_props.ISRF_c_hyp_scheme's parsed value.
 * @param fixed_fraction #feedback_props.ISRF_c_hyp_fixed_fraction_of_c's
 * parsed value.
 */
__attribute__((always_inline)) INLINE static void
feedback_props_check_c_hyp_scheme(int scheme, float fixed_fraction) {
  if (scheme == isrf_c_hyp_scheme_fixed_fraction && fixed_fraction <= 0.f)
    error(
        "GEARFeedback:ISRF_c_hyp_scheme is set to fixed_fraction "
        "but GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c is 0: there is "
        "nothing to select. Set a positive fraction.");
  if (scheme != isrf_c_hyp_scheme_fixed_fraction && fixed_fraction > 0.f)
    error(
        "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c is set (%g) but "
        "GEARFeedback:ISRF_c_hyp_scheme is %s, not fixed_fraction: "
        "the two speed schemes are alternatives, not layers. Set "
        "ISRF_c_hyp_scheme to fixed_fraction to use this fraction, or leave "
        "it at 0 if it was set by mistake.",
        fixed_fraction, feedback_props_c_hyp_scheme_name(scheme));
}

/**
 * @brief Enforce that the operator -> owning-moment map is the first-match
 * reverse of the moment -> operator map.
 *
 * @p owner[o] must be the lowest-index entry of @p forward equal to @p o.
 * Every operator must own a moment. The maps are arguments so a unit test can
 * pass synthetic ones.
 *
 * @param forward Moment -> operator map (#radiation_isrf_moment_to_operator
 * or a test stand-in).
 * @param n_moments Number of entries in @p forward.
 * @param owner Operator -> owning-moment map (#radiation_isrf_operator_owner
 * or a test stand-in).
 * @param n_operators Number of entries in @p owner.
 */
__attribute__((always_inline)) INLINE static void
feedback_check_isrf_operator_owner_map(
    const enum radiation_isrf_operator *forward, int n_moments,
    const enum radiation_isrf_moment *owner, int n_operators) {

  for (int m = 0; m < n_moments; m++)
    if ((int)forward[m] < 0 || (int)forward[m] >= n_operators)
      error(
          "the moment->operator map's entry %d is %d, out of range "
          "[0, %d): would index feedback_part_data.isrf_operator out of "
          "bounds.",
          m, (int)forward[m], n_operators);

  for (int o = 0; o < n_operators; o++)
    if ((int)owner[o] < 0 || (int)owner[o] >= n_moments)
      error(
          "the operator->owner map's entry %d is %d, out of range "
          "[0, %d): would index moment-keyed state out of bounds.",
          o, (int)owner[o], n_moments);

  for (int o = 0; o < n_operators; o++) {
    int first = -1;
    for (int m = 0; m < n_moments; m++)
      if ((int)forward[m] == o) {
        first = m;
        break;
      }
    if (first == -1)
      error(
          "no moment maps to operator %d through the moment->operator "
          "map: every operator must own at least one moment.",
          o);
    if ((int)owner[o] != first)
      error(
          "the operator->owner map's entry %d is %d, but the first moment "
          "mapping to operator %d is %d: the owner must be the first "
          "moment that maps to the operator, not merely any moment that "
          "does.",
          o, (int)owner[o], o, first);
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

  /* Several fields are only assigned in conditional branches but read
   * unconditionally downstream. sink_props_init() writes with_sinks after
   * this function returns. */
  bzero(fp, sizeof(struct feedback_props));

  /* Every operator-state writer trusts the owner map. Mirrored on restart by
   * feedback_struct_restore(). */
  feedback_check_isrf_operator_owner_map(
      radiation_isrf_moment_to_operator, ISRF_MOMENT_COUNT,
      radiation_isrf_operator_owner, ISRF_OPERATOR_COUNT);

  /* Supernovae energy efficiency */
  double e_efficiency =
      parser_get_param_double(params, "GEARFeedback:supernovae_efficiency");

  /* 0 is legal: it runs the enrichment channel without the thermal one. */
  if (e_efficiency < 0.0)
    error(
        "GEARFeedback:supernovae_efficiency is %g (< 0): a negative "
        "efficiency makes a supernova take thermal energy out of its gas "
        "neighbours instead of injecting it. Use 0 to inject no energy, or "
        "a positive value.",
        e_efficiency);

  fp->supernovae_efficiency = e_efficiency;

  /* Activate the stellar wind feedback */
  char with_stellar_wind_feedback = (char)parser_get_param_int(
      params, "GEARFeedback:with_stellar_wind_feedback");
  fp->with_stellar_wind_feedback = with_stellar_wind_feedback;

  /* Are we running with photoionization? */
  const char with_photoionization = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_photoionization", 0);

  /* Radiation pressure. The efficiency is read only when this is on. */
  const char with_radiation_pressure = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_radiation_pressure", 0);

  /* Stays 0 when radiation pressure is off. */
  float radiation_pressure_efficiency = 0.f;
  if (with_radiation_pressure) {
    radiation_pressure_efficiency = parser_get_opt_param_float(
        params, "GEARFeedback:radiation_pressure_efficiency", 0.0);
    if (radiation_pressure_efficiency <= 0.0f)
      error(
          "GEARFeedback:with_radiation_pressure is on but "
          "GEARFeedback:radiation_pressure_efficiency is %g (<= 0): there is "
          "nothing to inject. Set radiation_pressure_efficiency to a "
          "positive value (1 reproduces the table's own unboosted L_bol).",
          radiation_pressure_efficiency);
  }

  /* Are we running with the local Lyman-Werner/PE feedback (photoelectric
   * heating + H2 photodissociation)? */
  const char with_interstellar_radiation_field = (char)parser_get_opt_param_int(
      params, "GEARFeedback:with_interstellar_radiation_field", 0);

  /* Every channel switch must appear here. Keep in sync with
   * feedback_struct_restore(). */
  const char with_radiation = with_photoionization || with_radiation_pressure ||
                              with_interstellar_radiation_field;

  /* Pre-Supernovae energy efficiency */
  double w_efficiency = 0.0;
  if (with_stellar_wind_feedback) {
    w_efficiency = parser_get_param_double(
        params, "GEARFeedback:stellar_winds_efficiency");
  }

  /* A negative efficiency gives a NaN wind momentum, which -ffast-math cannot
   * catch downstream. */
  if (w_efficiency < 0.0)
    error(
        "GEARFeedback:stellar_winds_efficiency is %g (< 0): the wind "
        "momentum is a square root of the ejected energy, so a negative "
        "efficiency produces NaN momentum and internal energy. Use 0 to "
        "inject no wind energy, or a positive value.",
        w_efficiency);

  fp->winds_efficiency = w_efficiency;

  /* filename of the chemistry tables. */
  parser_get_param_string(params, "GEARFeedback:yields_table",
                          fp->stellar_model.yields_table);

  /* Initialize the stellar models. */
  stellar_evolution_props_init(&fp->stellar_model, phys_const, us, params,
                               cosmo, fp->with_stellar_wind_feedback,
                               with_radiation);

  /* Announce the H2 photodissociation coefficient, once per run (not for the
     first-stars model). */
  radiation_set_lw_photon_energy_cgs(&fp->stellar_model.rad,
                                     &fp->stellar_model);

  /* Announce the band-edge weights. */
  radiation_set_band_edge_coefficients(fp, &fp->stellar_model.rad,
                                       &fp->stellar_model);

  /* Read the metallicity threshold */
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
  /* Parsed unconditionally: applies to every star. */
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
        "run, not just radiation ones: <= 0 reopens the get_spart_timestep "
        "crash this parameter exists to prevent.",
        fp->event_dt_floor_Myr);

  fp->event_dt_floor_Myr *= Myr_internal_units;

  /* The floor must exceed dt_min to avoid get_spart_timestep()'s error. */
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

  /* Radiation pressure. The policy bit gates the injection. */
  fp->radiation_pressure_efficiency = radiation_pressure_efficiency;

  if (with_radiation_pressure) {
    fp->radiation_policy |= radiation_policy_radiation_pressure;
  }

  /* The magnitude keys are parsed only by the mechanism that reads them. */
  fp->ISRF_extinction_path_in_kernel_radii = 0.f;
  fp->ISRF_extinction_jeans_temperature_cap_K = 0.f;

  char extinction_path[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "GEARFeedback:ISRF_extinction_path",
                              extinction_path, "pair_separation");

  if (strcmp(extinction_path, "constant_kernel_path") == 0) {
    fp->ISRF_extinction_path_mechanism =
        (char)isrf_extinction_path_constant_kernel_path;
    fp->ISRF_extinction_path_in_kernel_radii = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_extinction_path_in_kernel_radii",
        default_ISRF_extinction_path_in_kernel_radii);
    if (fp->ISRF_extinction_path_in_kernel_radii <= 0.f)
      error(
          "GEARFeedback:ISRF_extinction_path_in_kernel_radii must be "
          "positive, got %g.",
          fp->ISRF_extinction_path_in_kernel_radii);
  } else if (strcmp(extinction_path, "pair_separation") == 0) {
    fp->ISRF_extinction_path_mechanism =
        (char)isrf_extinction_path_pair_separation;
  } else if (strcmp(extinction_path, "temperature_capped_jeans") == 0) {
    fp->ISRF_extinction_path_mechanism =
        (char)isrf_extinction_path_temperature_capped_jeans;
    fp->ISRF_extinction_jeans_temperature_cap_K = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_extinction_jeans_temperature_cap_K",
        default_ISRF_extinction_jeans_temperature_cap_K);
    if (fp->ISRF_extinction_jeans_temperature_cap_K <= 0.f)
      error(
          "GEARFeedback:ISRF_extinction_jeans_temperature_cap_K must be "
          "positive, got %g.",
          fp->ISRF_extinction_jeans_temperature_cap_K);
  } else {
    /* Stop the run on an unknown value rather than fall back to a default. */
    error(
        "GEARFeedback:ISRF_extinction_path must be one of "
        "constant_kernel_path, pair_separation or temperature_capped_jeans, "
        "got '%s'. A run archived before this parameter existed recorded no "
        "value and ran at two kernel support radii: reproduce it with "
        "constant_kernel_path and "
        "ISRF_extinction_path_in_kernel_radii: 2.0.",
        extinction_path);
  }

  if (with_interstellar_radiation_field) {
    fp->radiation_policy |= radiation_policy_isrf;

    fp->ISRF_propagation = (char)parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_propagation", 0);

    /* Parsed even with ISRF_propagation off, so validation runs can set it.
     * Default is kernel_local. */
    char c_hyp_scheme[PARSER_MAX_LINE_SIZE];
    parser_get_opt_param_string(params, "GEARFeedback:ISRF_c_hyp_scheme",
                                c_hyp_scheme,
                                isrf_c_hyp_scheme_name_kernel_local);
    fp->ISRF_c_hyp_scheme = feedback_props_c_hyp_scheme_from_name(c_hyp_scheme);
    if (fp->ISRF_c_hyp_scheme < 0)
      error(
          "GEARFeedback:ISRF_c_hyp_scheme must be one of %s (kernel-local "
          "speed) or %s (fixed fraction of c), got '%s'. The integer values "
          "are no longer accepted.",
          isrf_c_hyp_scheme_name_kernel_local,
          isrf_c_hyp_scheme_name_fixed_fraction, c_hyp_scheme);

    /* 0 disables it. Only valid with the fixed_fraction scheme, enforced
     * below. */
    fp->ISRF_c_hyp_fixed_fraction_of_c = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c", 0.0f);
    fp->ISRF_c_hyp_timestep_term_off_for_debugging =
        (char)parser_get_opt_param_int(
            params, "GEARFeedback:ISRF_c_hyp_timestep_term_off_for_debugging",
            0);

    if (fp->ISRF_c_hyp_fixed_fraction_of_c < 0.f ||
        fp->ISRF_c_hyp_fixed_fraction_of_c > 1.f)
      error(
          "GEARFeedback:ISRF_c_hyp_fixed_fraction_of_c must lie in [0, 1] "
          "(got %g): it is a fraction of the true speed of light c. 0 "
          "disables it.",
          fp->ISRF_c_hyp_fixed_fraction_of_c);

    /* The two speed schemes are alternatives. */
    feedback_props_check_c_hyp_scheme(fp->ISRF_c_hyp_scheme,
                                      fp->ISRF_c_hyp_fixed_fraction_of_c);

    /* Parsed even with ISRF_propagation off, like the dissipation
     * coefficients below. */
    fp->ISRF_c_hyp_margin = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_c_hyp_margin", 0.5f);

    /* Loosest case of the joint bound below, at alpha = 0:
     * C_hyp <= sqrt(2/0.70) ~ 1.6903. */
    const float ISRF_c_hyp_absolute_bound = sqrtf(2.f / 0.70f);

    if (fp->ISRF_c_hyp_margin <= 0.f ||
        fp->ISRF_c_hyp_margin > ISRF_c_hyp_absolute_bound)
      error(
          "GEARFeedback:ISRF_c_hyp_margin must lie in "
          "(0, %g] (got %g): above this bound the staggered "
          "exact-relaxation scheme is no longer stable for every "
          "lambda/h at this project's kernel/eta_neighbours choice, even "
          "with dissipation (ISRF_dissipation_alpha_max/alpha_floor) "
          "fully disabled (6.2*alpha*C_hyp + 0.70*C_hyp^2 <= 2 at "
          "alpha = 0).",
          ISRF_c_hyp_absolute_bound, fp->ISRF_c_hyp_margin);

    /* Negativity-triggered artificial dissipation. */
    fp->ISRF_dissipation_alpha_max = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_max", 0.5f);
    fp->ISRF_dissipation_negativity_threshold = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_negativity_threshold", 0.01f);
    fp->ISRF_dissipation_alpha_pin_for_debugging = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging", 0.0f);

    if (fp->ISRF_dissipation_alpha_max < 0.f)
      error("GEARFeedback:ISRF_dissipation_alpha_max must be >= 0 (got %g).",
            fp->ISRF_dissipation_alpha_max);

    if (fp->ISRF_dissipation_negativity_threshold <= 0.f ||
        fp->ISRF_dissipation_negativity_threshold > 1.f)
      error(
          "GEARFeedback:ISRF_dissipation_negativity_threshold must lie "
          "in (0, 1] (got %g).",
          fp->ISRF_dissipation_negativity_threshold);

    /* Diffuse-phase floor under the trigger. */
    fp->ISRF_dissipation_alpha_floor = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_alpha_floor", 0.5f);
    fp->ISRF_dissipation_floor_h_over_lambda = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_floor_h_over_lambda", 0.5f);

    if (fp->ISRF_dissipation_alpha_floor < 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_alpha_floor must be >= 0 (got "
          "%g).",
          fp->ISRF_dissipation_alpha_floor);

    if (fp->ISRF_dissipation_floor_h_over_lambda <= 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_floor_h_over_lambda must be > 0 "
          "(got %g).",
          fp->ISRF_dissipation_floor_h_over_lambda);

    /* Flux-relaxation residual gate: 0 disables it. */
    fp->ISRF_dissipation_floor_relaxation_residual = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_dissipation_floor_relaxation_residual",
        0.40f);

    if (fp->ISRF_dissipation_floor_relaxation_residual < 0.f ||
        fp->ISRF_dissipation_floor_relaxation_residual > 1.f)
      error(
          "GEARFeedback:ISRF_dissipation_floor_relaxation_residual must "
          "lie in [0, 1] (got %g): the residual R it gates is itself bounded "
          "to [0, 1] by the triangle inequality, so a value above 1 would "
          "scale the floor down on every particle, fronts included. 0 "
          "disables the gate.",
          fp->ISRF_dissipation_floor_relaxation_residual);

    /* Joint stability bound: 6.2 = 2*I_W and 0.70 = nu_max_coeff^2/2, the
     * Wendland-C2 lattice constants (I_W = 3.10, nu_max_coeff = 1.18). The
     * floor is not gated on negativity, so check max(alpha_max,
     * alpha_floor). */
    const float C_hyp = fp->ISRF_c_hyp_margin;
    const float alpha_bound = (2.f - 0.70f * C_hyp * C_hyp) / (6.2f * C_hyp);
    const float alpha_joint_ceiling =
        max(fp->ISRF_dissipation_alpha_max, fp->ISRF_dissipation_alpha_floor);
    if (alpha_joint_ceiling > 0.f && alpha_joint_ceiling > alpha_bound)
      error(
          "max(GEARFeedback:ISRF_dissipation_alpha_max, "
          "GEARFeedback:ISRF_dissipation_alpha_floor) = %g exceeds the "
          "stability bound %g at GEARFeedback:ISRF_c_hyp_margin = %g "
          "(6.2*alpha*C_hyp + 0.70*C_hyp^2 <= 2). alpha_max = %g, "
          "alpha_floor = %g. Lower GEARFeedback:ISRF_c_hyp_margin to "
          "raise the bound, or lower alpha_max/alpha_floor to fit it.",
          alpha_joint_ceiling, alpha_bound, C_hyp,
          fp->ISRF_dissipation_alpha_max, fp->ISRF_dissipation_alpha_floor);

    if (fp->ISRF_dissipation_alpha_pin_for_debugging < 0.f)
      error(
          "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging must be "
          ">= 0 (got %g).",
          fp->ISRF_dissipation_alpha_pin_for_debugging);

    if (fp->ISRF_dissipation_alpha_pin_for_debugging > 0.f &&
        fp->ISRF_dissipation_alpha_pin_for_debugging > alpha_bound)
      warning(
          "GEARFeedback:ISRF_dissipation_alpha_pin_for_debugging (%g) "
          "exceeds the stability bound %g at the run's C_hyp: not itself "
          "clamped (debug/test only). Never use this in a production run.",
          fp->ISRF_dissipation_alpha_pin_for_debugging, alpha_bound);

    if (fp->ISRF_propagation) {
      /* The IC metallicity is not visible here, so warn unconditionally. */
      warning(
          "GEARFeedback:ISRF_propagation is on together with "
          "GEARFeedback:with_interstellar_radiation_field. The propagation's "
          "only loss channel is dust absorption, whose rate is proportional "
          "to the gas metallicity. Gas at or near zero metallicity has no "
          "loss channel. Check GEARChemistry:initial_metallicity and any "
          "MetalMassFraction field in the initial conditions before a long "
          "run.");
    }
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
          "integrated over, in every cadence mode: <= 0 silently zeroes "
          "every star's first-pass budget.",
          fp->HII_rebuild_floor_Myr);

    /* Convert to internal units */
    const double m_p_cgs = phys_const->const_proton_mass *
                           units_cgs_conversion_factor(us, UNIT_CONV_MASS);
    fp->HII_min_density *=
        m_p_cgs / units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);

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
