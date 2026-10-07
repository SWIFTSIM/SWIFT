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
#ifndef SWIFT_FEEDBACK_GEAR_RADIATION_PROPERTIES_H
#define SWIFT_FEEDBACK_GEAR_RADIATION_PROPERTIES_H

/**
 * @file src/feedback/GEAR/radiation_properties.h
 * @brief Parameters of the GEAR subgrid radiation, shared by the GEAR
 * feedback modules.
 *
 * The functions read and write the radiation fields of #feedback_props by
 * name. Include this file after the module's definition of #feedback_props.
 */

#include "error.h"
#include "minmax.h"
#include "parser.h"
#include "physical_constants.h"
#include "radiation_selection.h"
#include "radiation_struct.h"
#include "units.h"

#include <math.h>
#include <string.h>

#define default_HII_min_density_Hpcm3 1.0
#define default_HII_max_age_Myr 50.0
#define default_HII_rebuild_time_Myr 0.5
#define default_HII_deterministic_boundary_ionization 0
#define default_HII_rebuild_floor_Myr 1e-4
/* One kernel support radius: the largest path with the star inside the
 * receiver's kernel. */
#define default_ISRF_extinction_path_in_kernel_radii 1.0

/* Temperature cap in K of the temperature_capped_jeans extinction path:
 * Safranek-Shrader et al. (2017), MNRAS 465, 885, Section 3.4. */
#define default_ISRF_extinction_jeans_temperature_cap_K 40.0

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
 * @brief Print the subgrid radiation part of the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void
feedback_props_print_radiation(const struct feedback_props *feedback_props) {

  message("Subgrid radiation parts compiled (--with-subgrid-radiation) =%s%s%s",
          RADIATION_COMPILED_PRESSURE ? " rp" : "",
          RADIATION_COMPILED_HII ? " hii" : "",
          RADIATION_COMPILED_ISRF ? " isrf" : "");

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
 * @brief Check and announce the subgrid radiation part of a #feedback_props
 * just read from a restart file.
 *
 * Stops with error() if the restored ISRF c_hyp scheme is not a valid one.
 *
 * @param feedback The restored #feedback_props.
 */
__attribute__((always_inline)) INLINE static void
feedback_props_restore_radiation(const struct feedback_props *feedback) {

  /* The two ISRF operator/moment maps are compile-time constants, so
   * feedback_props_init()'s own check of them applies unchanged here; it
   * does not run on a restart, so its call is mirrored explicitly. */
  feedback_check_isrf_operator_owner_map(
      radiation_isrf_moment_to_operator, ISRF_MOMENT_COUNT,
      radiation_isrf_operator_owner, ISRF_OPERATOR_COUNT);

  /* The flat block read of feedback_struct_restore() bypasses
   * feedback_props_init()'s parse-time check of the scheme, so a restart
   * written by a run that used a removed scheme would otherwise resume with a
   * speed rule nothing sets. The field is only meaningful, and only validated
   * at parse time, when the interstellar radiation field is on. */
  if ((feedback->radiation_policy & radiation_policy_isrf) &&
      feedback->ISRF_c_hyp_scheme != isrf_c_hyp_scheme_fixed_fraction &&
      feedback->ISRF_c_hyp_scheme !=
          isrf_c_hyp_scheme_kernel_local_reduced_flux)
    error(
        "The restart file holds GEARFeedback:ISRF_c_hyp_scheme = %d (internal "
        "value), which is neither fixed_fraction (2) nor kernel_local (4). A "
        "restart written with a removed scheme cannot be resumed. Rerun the "
        "simulation from its initial conditions with fixed_fraction or "
        "kernel_local.",
        feedback->ISRF_c_hyp_scheme);

  /* feedback->band_edge_weight_pe/lw/photon_weight_lw need no re-derivation,
     unlike radiation_lw_photon_energy_cgs: they are plain fields of
     *feedback, already restored verbatim by the flat restart_read_blocks()
     call of feedback_struct_restore(). Announcing the
     restored value (not re-deriving it) still lets a restarted run's log be
     checked against its own start-up announcement, confirming the restart
     path preserves this value across a change to the radiation sub-struct. */
  if (engine_rank == 0 && feedback->radiation_policy != 0)
    message(
        "Band-edge weights restored from the restart file: lambda_E(PE)=%.5g, "
        "lambda_E(LW)=%.5g, lambda_N(LW)=%.5g",
        feedback->band_edge_weight_pe, feedback->band_edge_weight_lw,
        feedback->band_edge_photon_weight_lw);
}

/**
 * @brief Read the switches of the subgrid radiation channels and the
 * radiation pressure efficiency, and set #feedback_props.radiation_policy.
 *
 * Called before the stellar models are read: they load the radiation tables
 * only when a channel is on.
 *
 * @param fp The #feedback_props.
 * @param params The parsed parameters.
 * @return True if any subgrid radiation channel is on.
 */
__attribute__((always_inline)) INLINE static char
feedback_props_init_radiation_switches(struct feedback_props *fp,
                                       struct swift_params *params) {

  /* The keys of a part this build lacks are ignored, with one warning. */
  radiation_selection_warn_absent_parts(params);

  /* Are we running with photoionization? */
  const char with_photoionization = (char)radiation_selection_get_switch(
      params, "GEARFeedback:with_photoionization", RADIATION_COMPILED_HII);

  /* Radiation pressure. The efficiency is read only when this is on. */
  const char with_radiation_pressure = (char)radiation_selection_get_switch(
      params, "GEARFeedback:with_radiation_pressure",
      RADIATION_COMPILED_PRESSURE);

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
  const char with_interstellar_radiation_field =
      (char)radiation_selection_get_switch(
          params, "GEARFeedback:with_interstellar_radiation_field",
          RADIATION_COMPILED_ISRF);

  fp->radiation_policy = 0;

  /* TODO: For the future, enforce these to have a non-zero value */

  /* Radiation pressure. The policy bit gates the injection. */
  fp->radiation_pressure_efficiency = radiation_pressure_efficiency;
  if (with_radiation_pressure)
    fp->radiation_policy |= radiation_policy_radiation_pressure;
  if (with_interstellar_radiation_field)
    fp->radiation_policy |= radiation_policy_isrf;
  if (with_photoionization)
    fp->radiation_policy |= radiation_policy_photoionization;

  /* Every channel switch must appear here. Keep in sync with
   * feedback_struct_restore(). */
  return with_photoionization || with_radiation_pressure ||
         with_interstellar_radiation_field;
}

/**
 * @brief Read the parameters of the subgrid radiation channels selected by
 * feedback_props_init_radiation_switches().
 *
 * @param fp The #feedback_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 */
__attribute__((always_inline)) INLINE static void feedback_props_init_radiation(
    struct feedback_props *fp, const struct phys_const *phys_const,
    const struct unit_system *us, struct swift_params *params) {

  const double Myr_internal_units = 1e6 * phys_const->const_year;

  /* The magnitude keys are parsed only by the mechanism that reads them. */
  fp->ISRF_extinction_path_in_kernel_radii = 0.f;
  fp->ISRF_extinction_jeans_temperature_cap_K = 0.f;

  /* A build without the ISRF keeps the default and reads no ISRF key. */
  char extinction_path[PARSER_MAX_LINE_SIZE] = "pair_separation";
  if (RADIATION_COMPILED_ISRF)
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

  if (fp->radiation_policy & radiation_policy_isrf) {

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

  if (fp->radiation_policy & radiation_policy_photoionization) {

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
}
#endif /* SWIFT_FEEDBACK_GEAR_RADIATION_PROPERTIES_H */
