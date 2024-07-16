/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_HYDRO_PARAMETERS_H
#define SWIFT_PLANETARY_HYDRO_PARAMETERS_H

/* Configuration file */
#include "config.h"

/* Global headers */
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local headers */
#include "common_io.h"
#include "error.h"
#include "inline.h"

/**
 * @file Planetary/hydro_parameters.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2024) (default parameters)
 *
 * This file defines a number of things that are used in
 * hydro_properties.c as defaults for run-time parameters
 * as well as a number of compile-time parameters.
 */

 /* Viscosity paramaters -- Defaults; can be changed at run-time */

 /*! Default REMIX artificial viscosity parameters */
 #define hydro_props_default_viscosity_alpha 1.5f
 #define hydro_props_default_viscosity_beta 3.f
 #define hydro_props_default_viscosity_epsilon 0.1f
 #define hydro_props_default_remix_a_visc 2.0f / 3.0f
 #define hydro_props_default_remix_b_visc 1.0f / 3.0f
 #define hydro_props_default_remix_a_difn_u 0.05f
 #define hydro_props_default_remix_b_difn_u 0.95f
 #define hydro_props_default_remix_a_difn_rho 0.05f
 #define hydro_props_default_remix_b_difn_rho 0.95f
 #define hydro_props_default_remix_alpha_norm 1.0f
 #define hydro_props_default_remix_eta 1.0f

 #define hydro_slope_limiter_exp_denom 0.04f

/* The viscosity that the particles are reset to after being hit by a
 * feedback event. This should be set to the same value as the
 * hydro_props_default_viscosity_alpha in fixed schemes, and likely
 * to hydro_props_default_viscosity_alpha_max in variable schemes. */
#define hydro_props_default_viscosity_alpha_feedback_reset 1.5f

/* Structs that store the relevant variables */

/*! Artificial viscosity parameters */
struct viscosity_global_data {
  /*! For the fixed, simple case. Also used to set the initial AV
      coefficient for variable schemes. */
  float alpha;
  float beta;
  float epsilon;
  float a_visc;
  float b_visc;
  float eta_crit;
};

/*! Artificial diffusion parameters */
struct diffusion_global_data {
  float a_difn_u;
  float b_difn_u;
  float a_difn_rho;
  float b_difn_rho;
  float alpha_norm;
};

/* Functions for reading from parameter file */

/* Forward declartions */
struct swift_params;
struct phys_const;
struct unit_system;

extern struct viscosity_global_data viscosity_global;
extern struct diffusion_global_data diffusion_global;

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
  viscosity->beta = parser_get_opt_param_float(
      params, "SPH:viscosity_beta", hydro_props_default_viscosity_beta);
  viscosity->epsilon = parser_get_opt_param_float(
      params, "SPH:viscosity_epsilon", hydro_props_default_viscosity_epsilon);
  viscosity->a_visc = parser_get_opt_param_float(
      params, "SPH:viscosity_a", hydro_props_default_remix_a_visc);
  viscosity->b_visc = parser_get_opt_param_float(
      params, "SPH:viscosity_b", hydro_props_default_remix_b_visc);
  viscosity->eta_crit = 1.0f / parser_get_opt_param_float(
      params, "SPH:resolution_eta", hydro_props_default_remix_eta);

  viscosity_global = *viscosity;
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
  viscosity->beta = 0.f;
  viscosity->epsilon = 1.f;
  viscosity->a_visc = 0.f;
  viscosity->b_visc = 0.f;
  viscosity->eta_crit = 0.f;
}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data* viscosity) {
  message("Artificial viscosity parameters set to alpha=%.3f, beta=%.3f, "
          "epsilon=%.3f, a_visc=%.3f, b_visc=%.3f, eta_crit=%.3f",
          viscosity->alpha, viscosity->beta, viscosity->epsilon, viscosity->a_visc,
          viscosity->b_visc, viscosity->eta_crit);
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

  io_write_attribute_f(h_grpsph, "Viscosity alpha", viscosity->alpha);
  io_write_attribute_f(h_grpsph, "Viscosity beta", viscosity->beta);
  io_write_attribute_f(h_grpsph, "Viscosity epsilon", viscosity->epsilon);
  io_write_attribute_f(h_grpsph, "Viscosity a_visc", viscosity->a_visc);
  io_write_attribute_f(h_grpsph, "Viscosity b_visc", viscosity->b_visc);
  io_write_attribute_f(h_grpsph, "Slope limiter eta_crit", viscosity->eta_crit);
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

    /* Read the artificial diffusion parameters from the file, if they exist,
   * otherwise set them to the defaults defined above. */

  diffusion->a_difn_u = parser_get_opt_param_float(
      params, "SPH:diffusion_a_u", hydro_props_default_remix_a_difn_u);
  diffusion->b_difn_u = parser_get_opt_param_float(
      params, "SPH:diffusion_b_u", hydro_props_default_remix_b_difn_u);
  diffusion->a_difn_rho = parser_get_opt_param_float(
      params, "SPH:diffusion_a_rho", hydro_props_default_remix_a_difn_rho);
  diffusion->b_difn_rho = parser_get_opt_param_float(
      params, "SPH:vdiffusion_b_rho", hydro_props_default_remix_b_difn_rho);
  diffusion->alpha_norm = parser_get_opt_param_float(
      params, "SPH:alpha_norm", hydro_props_default_remix_alpha_norm);

  diffusion_global = *diffusion;
}

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data* diffusion) {
  diffusion->a_difn_u = 0.f;
  diffusion->b_difn_u = 0.f;
  diffusion->a_difn_rho = 0.f;
  diffusion->b_difn_rho = 0.f;
  diffusion->alpha_norm = 0.f;
}

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data* diffusion) {
  message("Artificial diffusion parameters set to a_difn_u=%.3f, b_difn_u=%.3f, "
          "a_difn_rho=%.3f, b_difn_rho=%.3f, alpha_norm=%.3f",
          diffusion->a_difn_u, diffusion->b_difn_u, diffusion->a_difn_rho,
          diffusion->b_difn_rho, diffusion->alpha_norm);
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
  io_write_attribute_f(h_grpsph, "Diffusion a_difn_u", diffusion->a_difn_u);
  io_write_attribute_f(h_grpsph, "Diffusion b_difn_u", diffusion->b_difn_u);
  io_write_attribute_f(h_grpsph, "Diffusion a_difn_rho", diffusion->a_difn_rho);
  io_write_attribute_f(h_grpsph, "Diffusion b_difn_rho", diffusion->b_difn_rho);
  io_write_attribute_f(h_grpsph, "Normalising alpha_norm", diffusion->alpha_norm);
}
#endif

#endif /* SWIFT_PLANETARY_HYDRO_PARAMETERS_H */
