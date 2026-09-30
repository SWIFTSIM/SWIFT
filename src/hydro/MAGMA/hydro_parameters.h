/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifndef SWIFT_MAGMA_HYDRO_PARAMETERS_H
#define SWIFT_MAGMA_HYDRO_PARAMETERS_H

/* Configuration file */
#include <config.h>

/* Global headers */
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif
#include <math.h>

/* Local headers */
#include "common_io.h"
#include "error.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "parser.h"

// #define TRADITIONAL_SPH_ACCELERATION_TERM
// #define USE_ZEROTH_ORDER_VELOCITIES
// #define USE_STANDARD_KERNEL_GRADIENTS
#define GRAVITY_DIFF_VELOCITY

/**
 * @file MAGMA/hydro_parameters.h
 * @brief MAGMA-2 implementation of SPH (Rosswog. 2020)
 *
 *        This file defines a number of things that are used in
 *        hydro_properties.c as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 *
 *        All the constants of the scheme can be changed at run time through
 *        the (optional) parameters listed in viscosity_init() and
 *        diffusion_init(). The interaction functions have no access to the
 *        #hydro_props, so the values are also kept in the global variables
 *        magma_viscosity and magma_diffusion (defined in hydro.c), following
 *        the pattern of the equation-of-state parameters.
 */

/* Default values of the run-time parameters */

/*! Linear viscosity coefficient alpha (Rosswog 2020, eq. 14) */
#define hydro_props_default_viscosity_alpha 1.0f

/*! Quadratic viscosity coefficient beta (Rosswog 2020, eq. 14) */
#define hydro_props_default_viscosity_beta 2.0f

/*! Softening of the viscosity velocity jump, in units of h (eq. 15). Also
 * softens the approach velocity of the Courant condition (eq. 36). */
#define hydro_props_default_viscosity_epsilon 0.1f

/*! Distance (in units of h) below which the slope limiter switches the
 * reconstruction off (eq. 21-23). A negative value means the mean
 * inter-particle separation, 1 / eta_neighbours, as in the paper (whose
 * eta_crit = (32 pi / 3 N_ngb)^(1/3) is that separation in units of half the
 * kernel support). */
#define hydro_props_default_limiter_eta_crit -1.f

/*! Width (in units of h) of the exponential cut-off of the limiter (eq. 21).
 * A negative value means the paper's 0.2 in units of half the kernel support,
 * i.e. 0.1 * kernel_gamma. */
#define hydro_props_default_limiter_width -1.f

/*! Maximal condition number of the C-matrix before falling back to standard
 * SPH gradients for the particle */
#define hydro_props_default_gradient_max_condition_number 60.f

/*! Maximal angle (in radians) between the gradient function G and the
 * separation vector of a pair before falling back to kernel gradients for the
 * pair */
#define hydro_props_default_gradient_angle_limit 0.5f

/*! Thermal conduction coefficient alpha_u (Rosswog 2020, eq. 24) */
#define hydro_props_default_diffusion_alpha 0.05f

/* Structs that store the relevant variables */

/*! Artificial viscosity (and gradient-function) parameters */
struct viscosity_global_data {

  /*! Linear viscosity coefficient (eq. 14) */
  float alpha;

  /*! Quadratic viscosity coefficient (eq. 14) */
  float beta;

  /*! Softening of the velocity jump, in units of h (eq. 15) */
  float epsilon;

  /*! Limiter cut-off distance, in units of h (eq. 21-23) */
  float eta_crit;

  /*! Limiter cut-off width, in units of h (eq. 21) */
  float limiter_width;

  /*! Maximal condition number of the C-matrix */
  float max_condition_number;

  /*! Maximal angle between G and the pair separation (radians) */
  float angle_limit;

  /*! cos(angle_limit), pre-computed */
  float cos_angle_limit;
};

/*! Thermal diffusion parameters */
struct diffusion_global_data {

  /*! Conduction coefficient (eq. 24) */
  float alpha;
};

/*! Copies of the parameters accessible to the interaction functions */
extern struct viscosity_global_data magma_viscosity;
extern struct diffusion_global_data magma_diffusion;

/* Functions for reading from parameter file */

/* Forward declarations */
struct swift_params;
struct phys_const;
struct unit_system;

/* Viscosity */

/**
 * @brief Sets the viscosity parameters to their default values.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 * @param eta_neighbours: resolution parameter (h in units of the mean
 * inter-particle separation) used to resolve the automatic limiter settings.
 **/
static INLINE void viscosity_set_defaults(
    struct viscosity_global_data *viscosity, const float eta_neighbours) {

  viscosity->alpha = hydro_props_default_viscosity_alpha;
  viscosity->beta = hydro_props_default_viscosity_beta;
  viscosity->epsilon = hydro_props_default_viscosity_epsilon;
  viscosity->eta_crit = hydro_props_default_limiter_eta_crit;
  viscosity->limiter_width = hydro_props_default_limiter_width;
  viscosity->max_condition_number =
      hydro_props_default_gradient_max_condition_number;
  viscosity->angle_limit = hydro_props_default_gradient_angle_limit;

  if (viscosity->eta_crit < 0.f) viscosity->eta_crit = 1.f / eta_neighbours;
  if (viscosity->limiter_width < 0.f)
    viscosity->limiter_width = 0.1f * kernel_gamma;
  viscosity->cos_angle_limit = cosf(viscosity->angle_limit);
}

/**
 * @brief Initialises the viscosity parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init(struct swift_params *params,
                                  const struct unit_system *us,
                                  const struct phys_const *phys_const,
                                  struct viscosity_global_data *viscosity) {

  viscosity->alpha = parser_get_opt_param_float(
      params, "SPH:viscosity_alpha", hydro_props_default_viscosity_alpha);
  viscosity->beta = parser_get_opt_param_float(
      params, "SPH:viscosity_beta", hydro_props_default_viscosity_beta);
  viscosity->epsilon = parser_get_opt_param_float(
      params, "SPH:viscosity_epsilon", hydro_props_default_viscosity_epsilon);
  viscosity->eta_crit = parser_get_opt_param_float(
      params, "SPH:limiter_eta_crit", hydro_props_default_limiter_eta_crit);
  viscosity->limiter_width = parser_get_opt_param_float(
      params, "SPH:limiter_width", hydro_props_default_limiter_width);
  viscosity->max_condition_number = parser_get_opt_param_float(
      params, "SPH:gradient_max_condition_number",
      hydro_props_default_gradient_max_condition_number);
  viscosity->angle_limit =
      parser_get_opt_param_float(params, "SPH:gradient_angle_limit",
                                 hydro_props_default_gradient_angle_limit);

  /* Automatic limiter settings: the paper's prescription in SWIFT's h
   * convention (kernel support = kernel_gamma * h, the paper uses 2 h) */
  if (viscosity->eta_crit < 0.f) {
    const float eta_neighbours =
        parser_get_param_float(params, "SPH:resolution_eta");
    viscosity->eta_crit = 1.f / eta_neighbours;
  }
  if (viscosity->limiter_width < 0.f)
    viscosity->limiter_width = 0.1f * kernel_gamma;

  if (viscosity->alpha < 0.f || viscosity->beta < 0.f ||
      viscosity->epsilon <= 0.f || viscosity->limiter_width <= 0.f ||
      viscosity->max_condition_number < 1.f || viscosity->angle_limit < 0.f ||
      viscosity->angle_limit > M_PI_2)
    error("Invalid MAGMA viscosity / gradient parameters.");

  viscosity->cos_angle_limit = cosf(viscosity->angle_limit);

  /* Make the values available to the interaction functions */
  magma_viscosity = *viscosity;
}

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init_no_hydro(
    struct viscosity_global_data *viscosity) {

  viscosity_set_defaults(viscosity, /*eta_neighbours=*/1.2348f);
  magma_viscosity = *viscosity;
}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data *viscosity) {
  message(
      "Artificial viscosity parameters set to alpha: %.3f, beta: %.3f, "
      "epsilon: %.3f.",
      viscosity->alpha, viscosity->beta, viscosity->epsilon);
  message("Slope limiter parameters set to eta_crit: %.4f h, width: %.4f h.",
          viscosity->eta_crit, viscosity->limiter_width);
  message(
      "Gradient-function fallbacks: condition number > %.1f (particle), "
      "angle > %.3f rad (pair).",
      viscosity->max_condition_number, viscosity->angle_limit);
}

#if defined(HAVE_HDF5)
/**
 * @brief Prints the viscosity information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param viscosity: pointer to the viscosity_global_data struct.
 **/
static INLINE void viscosity_print_snapshot(
    hid_t h_grpsph, const struct viscosity_global_data *viscosity) {

  io_write_attribute_f(h_grpsph, "Alpha viscosity", viscosity->alpha);
  io_write_attribute_f(h_grpsph, "Beta viscosity", viscosity->beta);
  io_write_attribute_f(h_grpsph, "Epsilon viscosity", viscosity->epsilon);
  io_write_attribute_f(h_grpsph, "Limiter eta_crit", viscosity->eta_crit);
  io_write_attribute_f(h_grpsph, "Limiter width", viscosity->limiter_width);
  io_write_attribute_f(h_grpsph, "Gradient max condition number",
                       viscosity->max_condition_number);
  io_write_attribute_f(h_grpsph, "Gradient angle limit",
                       viscosity->angle_limit);
}
#endif

/* Diffusion */

/**
 * @brief Initialises the diffusion parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param diffusion: pointer to the diffusion struct to be filled.
 **/
static INLINE void diffusion_init(struct swift_params *params,
                                  const struct unit_system *us,
                                  const struct phys_const *phys_const,
                                  struct diffusion_global_data *diffusion) {

  diffusion->alpha = parser_get_opt_param_float(
      params, "SPH:diffusion_alpha", hydro_props_default_diffusion_alpha);

  if (diffusion->alpha < 0.f) error("Invalid MAGMA diffusion parameters.");

  /* Make the values available to the interaction functions */
  magma_diffusion = *diffusion;
}

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data *diffusion) {

  diffusion->alpha = hydro_props_default_diffusion_alpha;
  magma_diffusion = *diffusion;
}

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data *diffusion) {
  message("Artificial conduction parameters set to alpha_u: %.3f.",
          diffusion->alpha);
}

#ifdef HAVE_HDF5
/**
 * @brief Prints the diffusion information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param diffusion: pointer to the diffusion_global_data struct.
 **/
static INLINE void diffusion_print_snapshot(
    hid_t h_grpsph, const struct diffusion_global_data *diffusion) {
  io_write_attribute_f(h_grpsph, "Alpha diffusion", diffusion->alpha);
}
#endif

#endif /* SWIFT_MAGMA_HYDRO_PARAMETERS_H */
