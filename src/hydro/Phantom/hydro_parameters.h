/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
 *
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

#ifndef SWIFT_PHANTOM_HYDRO_PARAMETERS_H
#define SWIFT_PHANTOM_HYDRO_PARAMETERS_H

/* Configuration file */
#include <config.h>

/* Global headers */
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local headers */
#include "common_io.h"
#include "error.h"
#include "inline.h"

/**
 * @file Phantom/hydro_parameters.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added diffusive physics (Cullen & Denhen 2011 AV,
 *        Price 2017 (PHANTOM) diffusion) (default compile-time
 *        parameters).
 *
 *        This is a base reference implementation
 *        similar to the one presented in Price 2018.
 *
 *        This file defines a number of things that are used in
 *        hydro_properties.c as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */

/*! Viscosity parameters -- FIXED -- MUST BE DEFINED AT COMPILE-TIME */

/*! Cosmology default beta=3.0.
 * Alpha can be set in the parameter file.
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta 3.0f

/*! The viscosity that the particles are reset to after being hit by a
 * feedback event. This should be set to the same value as the
 * hydro_props_default_viscosity_alpha in fixed schemes, and likely
 * to hydro_props_default_viscosity_alpha_max in variable schemes. */
#define hydro_props_default_viscosity_alpha_feedback_reset 2.0f

/* Viscosity paramaters -- Defaults; can be changed at run-time */

/*! The "initial" hydro viscosity, or the fixed value for non-variable
 * schemes. This usually takes the value 0.8. */
#define hydro_props_default_viscosity_alpha 0.1f

/*! Minimal value for the viscosity alpha in variable schemes. */
#define hydro_props_default_viscosity_alpha_min 0.0f

/*! Maximal value for the viscosity alpha in variable schemes. */
#define hydro_props_default_viscosity_alpha_max 2.0f

/*! Decay length for the viscosity scheme. This is scheme dependent. In
 * non-variable schemes this must be defined but is not used. This also
 * sets the decay length for the diffusion. */
#define hydro_props_default_viscosity_length 0.25f

/* Diffusion parameters -- Defaults; can be changed at run-time */

/*! The "initial" diffusion, or the fixed value for non-variable
 * schemes. For this fixed scheme, this just takes the value of unity. */
#define hydro_props_default_diffusion_alpha 1.0f

/* Structs that store the relevant variables */

/*! Artificial viscosity parameters */
struct viscosity_global_data {
  /*! For the fixed, simple case. Also used to set the initial AV
      coefficient for variable schemes. */
  float alpha;

  /*! Artificial viscosity (max) for the variable case (e.g. M&M) */
  float alpha_max;

  /*! Artificial viscosity (min) for the variable case (e.g. M&M) */
  float alpha_min;

  /*! The decay length of the artificial viscosity (used in M&M, etc.) */
  float length;
};

/*! Thermal diffusion parameters */
struct diffusion_global_data {

  /*! Initialisation value, or the case for constant thermal diffusion coeffs
   */
  float alpha;
};

/* Functions for reading from parameter file */

/* Forward declartions */
struct swift_params;
struct phys_const;
struct unit_system;

/* Viscosity */

/**
 * @brief Initialises the viscosity parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param unit_system: pointer to the unit system
 * @param phys_const: pointer to the physical constants system
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init(struct swift_params* params,
                                  const struct unit_system* us,
                                  const struct phys_const* phys_const,
                                  struct viscosity_global_data* viscosity) {

  /* Read the artificial viscosity parameters from the file, if they exist,
   * otherwise set them to the defaults defined above. */

  viscosity->alpha = parser_get_opt_param_float(
      params, "SPH:viscosity_alpha", hydro_props_default_viscosity_alpha);

  viscosity->alpha_max =
      parser_get_opt_param_float(params, "SPH:viscosity_alpha_max",
                                 hydro_props_default_viscosity_alpha_max);

  viscosity->alpha_min =
      parser_get_opt_param_float(params, "SPH:viscosity_alpha_min",
                                 hydro_props_default_viscosity_alpha_min);

  viscosity->length = parser_get_opt_param_float(
      params, "SPH:viscosity_length", hydro_props_default_viscosity_length);
}

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init_no_hydro(
    struct viscosity_global_data* viscosity) {
  viscosity->alpha = hydro_props_default_viscosity_alpha;
  viscosity->alpha_max = hydro_props_default_viscosity_alpha_max;
  viscosity->alpha_min = hydro_props_default_viscosity_alpha_min;
  viscosity->length = hydro_props_default_viscosity_length;
}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data* viscosity) {
  message(
      "Artificial viscosity parameters set to alpha: %.3f, max: %.3f, "
      "min: %.3f, length: %.3f.",
      viscosity->alpha, viscosity->alpha_max, viscosity->alpha_min,
      viscosity->length);
}

#if defined(HAVE_HDF5)
/**
 * @brief Prints the viscosity information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param viscosity: pointer to the viscosity_global_data struct.
 **/
static INLINE void viscosity_print_snapshot(
    hid_t h_grpsph, const struct viscosity_global_data* viscosity) {

  io_write_attribute_f(h_grpsph, "Alpha viscosity", viscosity->alpha);
  io_write_attribute_f(h_grpsph, "Alpha viscosity (max)", viscosity->alpha_max);
  io_write_attribute_f(h_grpsph, "Alpha viscosity (min)", viscosity->alpha_min);
  io_write_attribute_f(h_grpsph, "Viscosity decay length [internal units]",
                       viscosity->length);
  io_write_attribute_f(h_grpsph, "Beta viscosity", const_viscosity_beta);
}
#endif

/* Diffusion */

/**
 * @brief Initialises the diffusion parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param unit_system: pointer to the unit system
 * @param phys_const: pointer to the physical constants system
 * @param diffusion_global_data: pointer to the diffusion struct to be filled.
 **/
static INLINE void diffusion_init(struct swift_params* params,
                                  const struct unit_system* us,
                                  const struct phys_const* phys_const,
                                  struct diffusion_global_data* diffusion) {

  diffusion->alpha = parser_get_opt_param_float(
      params, "SPH:diffusion_alpha", hydro_props_default_diffusion_alpha);
}

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data* diffusion) {
  diffusion->alpha = hydro_props_default_diffusion_alpha;
}

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data* diffusion) {
  message("Artificial diffusion parameters set to alpha: %.3f (fixed).",
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
    hid_t h_grpsph, const struct diffusion_global_data* diffusion) {
  io_write_attribute_f(h_grpsph, "Diffusion alpha", diffusion->alpha);
}
#endif

#endif /* SWIFT_PHANTOM_HYDRO_PARAMETERS_H */
