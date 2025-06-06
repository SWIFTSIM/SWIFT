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
#ifndef SWIFT_REMIX_HYDRO_PARAMETERS_H
#define SWIFT_REMIX_HYDRO_PARAMETERS_H

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
 * @file REMIX/hydro_parameters.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2025) (default parameters)
 *
 * This file defines a number of things that are used in
 * hydro_properties.c as defaults for run-time parameters
 * as well as a number of compile-time parameters.
 */

/* Viscosity paramaters -- Defaults; can be changed at run-time */

/*! Default REMIX artificial viscosity parameters */
#define const_remix_visc_alpha 1.5f
#define const_remix_visc_beta 3.f
#define const_remix_visc_epsilon 0.1f
#define const_remix_visc_a 2.0f / 3.0f
#define const_remix_visc_b 1.0f / 3.0f
#define const_remix_difn_a_u 0.05f
#define const_remix_difn_b_u 0.95f
#define const_remix_difn_a_rho 0.05f
#define const_remix_difn_b_rho 0.95f
#define const_remix_norm_alpha 1.0f
#define const_remix_slope_limiter_exp_denom 0.04f

/* The viscosity that the particles are reset to after being hit by a
 * feedback event. This should be set to the same value as the
 * const_viscosity_alpha in fixed schemes, and likely
 * to const_viscosity_alpha_max in variable schemes. */
#define hydro_props_default_viscosity_alpha_feedback_reset 1.5f

/* Structs that store the relevant variables */

/*! Artificial viscosity parameters */
struct viscosity_global_data {};

/*! Artificial diffusion parameters */
struct diffusion_global_data {};

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
                                  struct viscosity_global_data* viscosity) {}

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init_no_hydro(
    struct viscosity_global_data* viscosity) {}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data* viscosity) {}

#if defined(HAVE_HDF5)
/**
 * @brief Prints the viscosity information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param viscosity: pointer to the viscosity_global_data struct.
 **/
static INLINE void viscosity_print_snapshot(
    hid_t h_grpsph, const struct viscosity_global_data* viscosity) {

  io_write_attribute_f(h_grpsph, "Viscosity alpha", const_remix_visc_alpha);
  io_write_attribute_f(h_grpsph, "Viscosity beta", const_remix_visc_beta);
  io_write_attribute_f(h_grpsph, "Viscosity epsilon", const_remix_visc_epsilon);
  io_write_attribute_f(h_grpsph, "Viscosity a_visc", const_remix_visc_a);
  io_write_attribute_f(h_grpsph, "Viscosity b_visc", const_remix_visc_b);
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
                                  struct diffusion_global_data* diffusion) {}

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data* diffusion) {}

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data* diffusion) {}

#ifdef HAVE_HDF5
/**
 * @brief Prints the diffusion information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param diffusion: pointer to the diffusion_global_data struct.
 **/
static INLINE void diffusion_print_snapshot(
    hid_t h_grpsph, const struct diffusion_global_data* diffusion) {
  io_write_attribute_f(h_grpsph, "Diffusion a_difn_u", const_remix_difn_a_u);
  io_write_attribute_f(h_grpsph, "Diffusion b_difn_u", const_remix_difn_b_u);
  io_write_attribute_f(h_grpsph, "Diffusion a_difn_rho",
                       const_remix_difn_a_rho);
  io_write_attribute_f(h_grpsph, "Diffusion b_difn_rho",
                       const_remix_difn_b_rho);
}
#endif

#endif /* SWIFT_REMIX_HYDRO_PARAMETERS_H */
