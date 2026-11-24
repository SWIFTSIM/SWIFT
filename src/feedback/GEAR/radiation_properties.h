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
#include "radiation_isrf.h"
#include "radiation_selection.h"
#include "radiation_struct.h"
#include "units.h"

#include <math.h>
#include <string.h>

#define default_HII_min_density_Hpcm3 1.0
#define default_HII_max_age_Myr 50.0
#define default_HII_rebuild_time_Myr 0.1
#define default_HII_deterministic_boundary_ionization 0
#define default_HII_random_cone_frame_per_star 1
#define default_HII_cone_frame_redraw_period 0
#define default_HII_rebuild_floor_Myr 1e-4
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
 * @brief Speed limit of the light of the packet transport
 * (GEARFeedback:ISRF_packet_light_speed). The integer values go to the restart
 * file.
 */
enum isrf_packet_light_speed {
  /*! A share arrives at the next step of its receiver, whatever the distance.
   * The light can move faster than c. For debugging only. */
  isrf_packet_light_speed_unlimited = 0,
  /*! A share arrives no earlier than its path divided by the limit speed,
   * in the frame of the gas: gas that moves away from the sender faster than
   * 0.99 of the limit speed gets its light early. The step is not bounded:
   * the light moves at most one kernel per gas step, so the gas step sets its
   * speed, typically far below the limit. Default. */
  isrf_packet_light_speed_cap = 1,
  /*! As the cap, and the time step is bounded so that the light can move at
   * the limit speed, at the cost of many more steps. */
  isrf_packet_light_speed_fixed = 2,
};

/*! Values of GEARFeedback:ISRF_packet_light_speed. */
#define isrf_packet_light_speed_name_unlimited "unlimited"
#define isrf_packet_light_speed_name_cap "cap"
#define isrf_packet_light_speed_name_fixed "fixed"

/*! Default limit speed of the packet transport, as a fraction of the speed of
 * light. A share to a receiver that moves away from the sender at the speed
 * `u` arrives after the separation over `v - u`, at most 100 times the
 * separation over `v`. The limit holds in the frame of the gas: a sender that
 * moves towards the receiver between its stamp and its booking, or a receiver
 * that speeds up after it, can make the front move at up to `v + u` in the
 * frame of the box, which is harmless for `v` near `c`. */
#define default_ISRF_packet_light_speed_fraction_of_c 0.99f

/*! Default power of the soft transmission cone. */
#define default_ISRF_packet_cone_power 4.f

/*! Default time constant of the running mean of the packet field, in Myr. */
#define default_ISRF_packet_read_mean_time_Myr 0.f

/*! Default angles in degrees of the frame of the packet reception slots: a
 * generic rotation, not aligned with the box. */
#define default_ISRF_packet_octant_frame_angle_x_deg 24.f
#define default_ISRF_packet_octant_frame_angle_y_deg 36.f
#define default_ISRF_packet_octant_frame_angle_z_deg 66.f

/*! Default seed offset of the packet octant frame. 0 is the frame of the
 * angles alone. */
#define default_ISRF_packet_octant_frame_seed_offset 0

/*! Default distance of the origin of the test beams behind the particle, in
 * smoothing lengths. */
#define default_ISRF_packet_test_beam_origin_distance_h 1000.f

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
 * @brief Do the stars build HII regions (GEARFeedback:with_photoionization)?
 *
 * @param feedback_props The #feedback_props.
 * @return True if the HII photoionization is on.
 */
__attribute__((always_inline)) INLINE static int
feedback_props_do_photoionization(const struct feedback_props *feedback_props) {
  return (feedback_props->radiation_policy &
          radiation_policy_photoionization) != 0;
}

/**
 * @brief Print the subgrid radiation part of the feedback model.
 *
 * @param feedback_props The #feedback_props
 */
__attribute__((always_inline)) INLINE static void
feedback_props_print_radiation(const struct feedback_props *feedback_props) {

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
    message("HII random cone frame per star                             = %s",
            feedback_props->HII_random_cone_frame_per_star ? "on" : "off");
    message("HII cone frame redraw period (star passes)                 = %d",
            feedback_props->HII_cone_frame_redraw_period);
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
    message("ISRF LW population dissociation factor (tables)            = %s",
            feedback_props->ISRF_packet_constants.LW_population_track ? "ON"
                                                                      : "OFF");
    message("ISRF LW dissociation factor bound                          = %.6g",
            feedback_props->ISRF_packet_constants.LW_dissociating_factor_max);
    if (feedback_props->ISRF_propagation) {
      message(
          "ISRF transport                                             = packet "
          "hop (the only scheme)");
      message("ISRF packet transmission cone power                        = %g",
              feedback_props->ISRF_packet_constants.cone_power);
      message("ISRF packet light speed limit                              = %s",
              feedback_props->ISRF_packet_light_speed ==
                      isrf_packet_light_speed_fixed
                  ? isrf_packet_light_speed_name_fixed
              : feedback_props->ISRF_packet_light_speed ==
                      isrf_packet_light_speed_cap
                  ? isrf_packet_light_speed_name_cap
                  : isrf_packet_light_speed_name_unlimited);
      if (feedback_props->ISRF_packet_light_speed !=
          isrf_packet_light_speed_unlimited)
        message(
            "ISRF packet limit speed (fraction of c)                    = %g",
            feedback_props->ISRF_packet_light_speed_fraction_of_c);
      if (feedback_props->ISRF_packet_light_speed ==
          isrf_packet_light_speed_fixed)
        message(
            "ISRF packet light speed time step margin                   = %g",
            feedback_props->ISRF_packet_light_speed_time_step_margin);
      message("ISRF packet field read mean time (internal units, 0 = off) = %g",
              feedback_props->ISRF_packet_read_mean_time);
      message("ISRF H2 line absorption along the packet path              = %s",
              feedback_props->ISRF_H2_line_absorption ? "ON" : "OFF");
      message(
          "ISRF packet octant frame rotation (row-major)              = "
          "%g %g %g, %g %g %g, %g %g %g",
          feedback_props->ISRF_packet_constants.octant_frame[0],
          feedback_props->ISRF_packet_constants.octant_frame[1],
          feedback_props->ISRF_packet_constants.octant_frame[2],
          feedback_props->ISRF_packet_constants.octant_frame[3],
          feedback_props->ISRF_packet_constants.octant_frame[4],
          feedback_props->ISRF_packet_constants.octant_frame[5],
          feedback_props->ISRF_packet_constants.octant_frame[6],
          feedback_props->ISRF_packet_constants.octant_frame[7],
          feedback_props->ISRF_packet_constants.octant_frame[8]);
    }
  }
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
 * @brief Read and check the speed limit of the light of the packet transport
 * (GEARFeedback:ISRF_packet_light_speed and
 * GEARFeedback:ISRF_packet_light_speed_fraction_of_c).
 *
 * A key left out gives the mode cap with ISRF_propagation, and unlimited
 * otherwise. Writes a warning for the mode unlimited, for which light can move
 * faster than c, and for a fraction above 0.99, for which gas faster than `(1
 * - fraction) c` can push the front above c. A limited mode also prints a
 * message that the limit applies in the gas frame, so that gas which moves
 * away faster than 0.99 of the limit speed gets its light early. Stops the run
 * on an unknown mode, on a fraction outside (0, 1] and on a limit without
 * ISRF_propagation, which must already be set in #fp.
 *
 * @param params The parameter file.
 * @param fp The #feedback_props, with ISRF_propagation set.
 */
INLINE static void feedback_props_read_packet_light_speed(
    struct swift_params *params, struct feedback_props *fp) {

  const int packet_hop = fp->ISRF_propagation;
  char light_speed[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(
      params, "GEARFeedback:ISRF_packet_light_speed", light_speed,
      packet_hop ? isrf_packet_light_speed_name_cap
                 : isrf_packet_light_speed_name_unlimited);
  if (strcmp(light_speed, isrf_packet_light_speed_name_unlimited) == 0)
    fp->ISRF_packet_light_speed = isrf_packet_light_speed_unlimited;
  else if (strcmp(light_speed, isrf_packet_light_speed_name_cap) == 0)
    fp->ISRF_packet_light_speed = isrf_packet_light_speed_cap;
  else if (strcmp(light_speed, isrf_packet_light_speed_name_fixed) == 0)
    fp->ISRF_packet_light_speed = isrf_packet_light_speed_fixed;
  else
    error(
        "GEARFeedback:ISRF_packet_light_speed must be %s, %s or %s, got '%s'.",
        isrf_packet_light_speed_name_unlimited,
        isrf_packet_light_speed_name_cap, isrf_packet_light_speed_name_fixed,
        light_speed);
  fp->ISRF_packet_light_speed_fraction_of_c = parser_get_opt_param_float(
      params, "GEARFeedback:ISRF_packet_light_speed_fraction_of_c",
      default_ISRF_packet_light_speed_fraction_of_c);
  if (fp->ISRF_packet_light_speed != isrf_packet_light_speed_unlimited &&
      !(fp->ISRF_packet_light_speed_fraction_of_c > 0.f &&
        fp->ISRF_packet_light_speed_fraction_of_c <= 1.f))
    error(
        "GEARFeedback:ISRF_packet_light_speed_fraction_of_c must be in "
        "(0, 1] (got %g).",
        fp->ISRF_packet_light_speed_fraction_of_c);
  fp->ISRF_packet_constants.light_speed_limited =
      fp->ISRF_packet_light_speed != isrf_packet_light_speed_unlimited;
  if (fp->ISRF_packet_light_speed != isrf_packet_light_speed_unlimited &&
      !fp->ISRF_propagation)
    error(
        "GEARFeedback:ISRF_packet_light_speed: %s acts only with "
        "GEARFeedback:ISRF_propagation: 1.",
        light_speed);
  if (packet_hop &&
      fp->ISRF_packet_light_speed == isrf_packet_light_speed_unlimited)
    warning(
        "GEARFeedback:ISRF_packet_light_speed: unlimited lets the light of "
        "the packet transport move faster than c. Use it for debugging only.");
  if (packet_hop &&
      fp->ISRF_packet_light_speed != isrf_packet_light_speed_unlimited)
    message(
        "GEARFeedback:ISRF_packet_light_speed: %s, light limited to v = f c "
        "in the gas frame.",
        light_speed);
  if (packet_hop &&
      fp->ISRF_packet_light_speed != isrf_packet_light_speed_unlimited &&
      fp->ISRF_packet_light_speed_fraction_of_c > 0.99f)
    warning(
        "GEARFeedback:ISRF_packet_light_speed_fraction_of_c = %g: gas that "
        "moves faster than (1 - fraction) c can make the light front move "
        "faster than c in the frame of the box.",
        fp->ISRF_packet_light_speed_fraction_of_c);
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

  struct radiation_packet_constants *packet_constants =
      &fp->ISRF_packet_constants;
  packet_constants->LW_population_track = 0;
  packet_constants->LW_dissociating_factor_max = 1.f;

  if (fp->radiation_policy & radiation_policy_isrf) {
    fp->ISRF_propagation = (char)parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_propagation", 0);

    /* The tables were read above. A model without the datasets keeps the
       factor 1 for its stars. */
    const struct radiation *rad_main = &fp->stellar_model.rad;
    if (rad_main->has_l_lw_dissociating) {
      packet_constants->LW_population_track = 1;
      packet_constants->LW_dissociating_factor_max =
          rad_main->lw_dissociating_factor_max;
    }
    if (fp->metallicity_max_first_stars != -1) {
      const struct radiation *rad_first = &fp->stellar_model_first_stars.rad;
      if (rad_first->has_l_lw_dissociating) {
        packet_constants->LW_population_track = 1;
        packet_constants->LW_dissociating_factor_max =
            max(packet_constants->LW_dissociating_factor_max,
                rad_first->lw_dissociating_factor_max);
      }
    }
    if (!(packet_constants->LW_dissociating_factor_max >= 1.f))
      error("The LW dissociation factor bound is %g, below 1.",
            packet_constants->LW_dissociating_factor_max);
    if (!packet_constants->LW_population_track &&
        packet_constants->LW_dissociating_factor_max != 1.f)
      error("The LW dissociation factor bound is %g without the factor.",
            packet_constants->LW_dissociating_factor_max);

    const int H2_line_absorption = parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_H2_line_absorption", 0);
    if (H2_line_absorption != 0 && H2_line_absorption != 1)
      error("GEARFeedback:ISRF_H2_line_absorption must be 0 or 1 (got %d).",
            H2_line_absorption);
    fp->ISRF_H2_line_absorption = (char)H2_line_absorption;
    if (fp->ISRF_H2_line_absorption && !fp->ISRF_propagation)
      error(
          "GEARFeedback:ISRF_H2_line_absorption: 1 acts on the transported "
          "light, which needs GEARFeedback:ISRF_propagation: 1.");
#ifndef GEAR_COOLING
    if (fp->ISRF_H2_line_absorption)
      error(
          "GEARFeedback:ISRF_H2_line_absorption: 1 needs the geajay cooling "
          "with an H2 network.");
#endif
#ifdef GEAR_CHEMISTRY_ELEMENT_COUNT
    const int n_elements = GEAR_CHEMISTRY_ELEMENT_COUNT;
#else
    const int n_elements = 0;
#endif
    if (fp->ISRF_H2_line_absorption && sizeof(struct part) == 4096)
      error(
          "sizeof(struct part) is 4096 B in this build (%d chemistry "
          "elements). A particle of exactly 4096 bytes costs about 7 per cent "
          "more CPU time in a test without H2 line absorption; the cost with "
          "it on was not measured. Change the size of "
          "feedback_part_data.ISRF_packet_stride_padding or the number of "
          "chemistry elements.",
          n_elements);

    /* Parsed even with ISRF_propagation off. */
    fp->ISRF_packet_light_speed_time_step_margin = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_packet_light_speed_time_step_margin", 0.5f);

    /* Largest accepted margin. */
    const float ISRF_time_step_margin_bound = sqrtf(2.f / 0.70f);

    if (fp->ISRF_packet_light_speed_time_step_margin <= 0.f ||
        fp->ISRF_packet_light_speed_time_step_margin >
            ISRF_time_step_margin_bound)
      error(
          "GEARFeedback:ISRF_packet_light_speed_time_step_margin must lie in "
          "(0, %g] (got %g).",
          ISRF_time_step_margin_bound,
          fp->ISRF_packet_light_speed_time_step_margin);

    feedback_props_read_packet_light_speed(params, fp);

    fp->ISRF_packet_constants.cone_power = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_packet_cone_power",
        default_ISRF_packet_cone_power);
    fp->ISRF_packet_read_mean_time = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_packet_read_mean_time_Myr",
        default_ISRF_packet_read_mean_time_Myr);
    if (!(fp->ISRF_packet_read_mean_time >= 0.f))
      error(
          "GEARFeedback:ISRF_packet_read_mean_time_Myr must be >= 0 (got %g).",
          fp->ISRF_packet_read_mean_time);
    fp->ISRF_packet_read_mean_time *= Myr_internal_units;

    float frame_angles_deg[3] = {default_ISRF_packet_octant_frame_angle_x_deg,
                                 default_ISRF_packet_octant_frame_angle_y_deg,
                                 default_ISRF_packet_octant_frame_angle_z_deg};
    parser_get_opt_param_float_array(
        params, "GEARFeedback:ISRF_packet_octant_frame_angles_deg", 3,
        frame_angles_deg);
    const int frame_seed_offset = parser_get_opt_param_int(
        params, "GEARFeedback:ISRF_packet_octant_frame_seed_offset",
        default_ISRF_packet_octant_frame_seed_offset);
    if (frame_seed_offset < 0)
      error(
          "GEARFeedback:ISRF_packet_octant_frame_seed_offset must be >= 0 "
          "(got %d).",
          frame_seed_offset);
    for (int k = 0; k < 3; k++)
      if (!(fabsf(frame_angles_deg[k]) < 1e4f))
        error(
            "GEARFeedback:ISRF_packet_octant_frame_angles_deg must be finite "
            "angles in degrees (got %g).",
            frame_angles_deg[k]);
    radiation_packet_build_octant_frame(frame_angles_deg, frame_seed_offset,
                                        fp->ISRF_packet_constants.octant_frame);

    const int n_test_beams =
        radiation_packet_read_test_beams(params, fp->ISRF_packet_test_beam);
    fp->ISRF_packet_test_beam_origin_distance_h = parser_get_opt_param_float(
        params, "GEARFeedback:ISRF_packet_test_beam_origin_distance_h",
        default_ISRF_packet_test_beam_origin_distance_h);
    uint32_t origin_bits;
    memcpy(&origin_bits, &fp->ISRF_packet_test_beam_origin_distance_h,
           sizeof(origin_bits));
    volatile uint32_t origin_opaque = origin_bits;
    if (((origin_opaque >> 23) & 0xFFu) == 0xFFu ||
        !(fp->ISRF_packet_test_beam_origin_distance_h > 0.f &&
          fp->ISRF_packet_test_beam_origin_distance_h <= 1e4f))
      error(
          "GEARFeedback:ISRF_packet_test_beam_origin_distance_h must be in "
          "(0, 1e4] (got %g).",
          fp->ISRF_packet_test_beam_origin_distance_h);
    if (n_test_beams > 0 && !fp->ISRF_propagation)
      error(
          "GEARFeedback:ISRF_packet_test_beam_* is set, but it acts only with "
          "GEARFeedback:ISRF_propagation: 1.");

    if (fp->ISRF_propagation) {
#ifndef HAVE_CHEALPIX
      error(
          "GEARFeedback:ISRF_propagation: 1 needs the HEALPix C library for "
          "the emission cones of the stars. Configure with --with-chealpix.");
#endif
      if (!(fp->ISRF_packet_constants.cone_power >= 1.f))
        error("GEARFeedback:ISRF_packet_cone_power must be >= 1 (got %g).",
              fp->ISRF_packet_constants.cone_power);
      /* Stars emit energy and photon number from one spectrum per band. */
      if (!fp->stellar_model.rad.has_mean_photon_energy_pe ||
          (fp->metallicity_max_first_stars != -1 &&
           !fp->stellar_model_first_stars.rad.has_mean_photon_energy_pe))
        error(
            "GEARFeedback:ISRF_propagation: 1 needs the datasets "
            "MeanPhotonEnergyPE and Integrated_MeanPhotonEnergyPE in the "
            "Data/Radiation group of every yields table.");
    }

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

    const int random_cone_frame = parser_get_opt_param_int(
        params, "GEARFeedback:HII_random_cone_frame_per_star",
        default_HII_random_cone_frame_per_star);
    if (random_cone_frame != 0 && random_cone_frame != 1)
      error(
          "GEARFeedback:HII_random_cone_frame_per_star must be 0 or 1, got %d.",
          random_cone_frame);
    fp->HII_random_cone_frame_per_star = (char)random_cone_frame;

    fp->HII_cone_frame_redraw_period = parser_get_opt_param_int(
        params, "GEARFeedback:HII_cone_frame_redraw_period",
        default_HII_cone_frame_redraw_period);
    if (fp->HII_cone_frame_redraw_period < 0)
      error("GEARFeedback:HII_cone_frame_redraw_period must be >= 0, got %d.",
            fp->HII_cone_frame_redraw_period);

    fp->HII_rebuild_floor_Myr =
        parser_get_opt_param_float(params, "GEARFeedback:HII_rebuild_floor_Myr",
                                   default_HII_rebuild_floor_Myr);
    if (fp->HII_rebuild_floor_Myr <= 0.f)
      error(
          "GEARFeedback:HII_rebuild_floor_Myr must be > 0 (got %g): it "
          "sets the interval a star's first HII pass issues ahead: <= 0 "
          "silently zeroes every star's first-pass budget.",
          fp->HII_rebuild_floor_Myr);
    if (fp->HII_rebuild_time > 0.f &&
        fp->HII_rebuild_floor_Myr > fp->HII_rebuild_time)
      error(
          "GEARFeedback:HII_rebuild_floor_Myr (%g) must not exceed "
          "GEARFeedback:HII_rebuild_time_Myr (%g): the first pass issues the "
          "floor interval ahead, so it must fit inside one rebuild interval "
          "(above twice the rebuild time, photons are lost).",
          fp->HII_rebuild_floor_Myr, fp->HII_rebuild_time);

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
