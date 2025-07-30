/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
 *               2025 Doug Rennehan (douglas.rennehan@gmail.com)
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

#ifndef SWIFT_MAGMA2_HYDRO_PARAMETERS_H
#define SWIFT_MAGMA2_HYDRO_PARAMETERS_H

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
 * @file MAGMA2/hydro_parameters.h
 * @brief Density-Energy non-conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (default compile-time
 *        parameters).
 *
 *        This file defines a number of things that are used in
 *        hydro_properties.c as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */


/* ---------- Viscosity & Conductivitiy parameters ---------- */


/*! Alpha viscosity, usually ~2. For lower neighbor number, should be higher */
#define const_viscosity_alpha 1.0

/*! Artificial conductivity alpha */
#define const_conductivity_alpha 0.05

/*! Desired number of neighbours -- CRITICAL that this matches hydro props */
#if defined(HYDRO_DIMENSION_1D)
#define const_kernel_target_neighbours 4.0
#elif defined(HYDRO_DIMENSION_2D)
#define const_kernel_target_neighbours 17.0
#else
#define const_kernel_target_neighbours 128.//57.0
#endif

/*! Use the alternative viscosity weighting (each particle has mu_i and mu_j) */
#define hydro_props_use_asymmetric_viscosity_mu

/*! Use the second-order velocities in v_ij * G_ij  */
#define hydro_props_use_second_order_velocities_in_divergence

/*! Use double precision for all matrix/vector operations */
#define hydro_props_use_double_precision


/* ---------- These parameters should not be changed ---------- */


/*! Consider matrix inversion to be ill-conditioned above this limit */
#ifdef hydro_props_use_double_precision
#define const_condition_number_upper_limit 99999.
#else
#define const_condition_number_upper_limit 999.
#endif

/*! Mean interparticle spacing for this kernel and neighbour number */
#define const_kernel_mean_spacing (kernel_gamma*(4. * M_PI / (3. * \
    powf((float)const_kernel_target_neighbours, 1. / 3.))))

/*! eta_crit Rosswog 2020 Eq 23. Of order the mean interparticle spacing. */
#define const_slope_limiter_eta_crit (4.f * const_kernel_mean_spacing)

/*! eta_fold from Frontiere+'17 Equation 51 */
#define const_slope_limiter_eta_fold 0.2

/*! Softening squared (epsilon^2) in Eq. 15 Rosswog 2020 */
#define const_viscosity_epsilon2 0.01

/*! Cosmology default const_viscosity_beta=2*const_viscosity_alpha
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta (2.0 * const_viscosity_alpha)


/* ---------- Structures for below ---------- */


/*! Artificial viscosity parameters */
struct viscosity_global_data { };

/*! Thermal diffusion parameters */
struct diffusion_global_data { };

/* Functions for reading from parameter file */

/* Forward declartions */
struct swift_params;
struct phys_const;
struct unit_system;

/* Define float or double depending on hydro_props_use_double_precision */
#if defined(hydro_props_use_double_precision)
typedef double hydro_real_t;
#else
typedef float hydro_real_t;
#endif

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
                                  struct viscosity_global_data* viscosity) { }

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void viscosity_init_no_hydro(
    struct viscosity_global_data* viscosity) { }

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void viscosity_print(
    const struct viscosity_global_data* viscosity) {
  message("Artificial viscosity alpha set to %.3f", const_viscosity_alpha);
  message("Artificial viscosity beta set to %.3f", const_viscosity_beta);
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

  io_write_attribute_f(h_grpsph, "Alpha viscosity", const_viscosity_alpha);
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
                                  struct diffusion_global_data* diffusion) { }

/**
 * @brief Initialises a diffusion struct to sensible numbers for mocking
 *        purposes.
 *
 * @param diffusion: pointer to the diffusion_global_data struct to be filled.
 **/
static INLINE void diffusion_init_no_hydro(
    struct diffusion_global_data* diffusion) { }

/**
 * @brief Prints out the diffusion parameters at the start of a run.
 *
 * @param diffusion: pointer to the diffusion_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void diffusion_print(
    const struct diffusion_global_data* diffusion) {
  message("Artificial conductivity alpha set to %.3f",
          const_conductivity_alpha);
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
  io_write_attribute_f(h_grpsph, "Conductivity alpha",
                       const_conductivity_alpha);
}
#endif

#endif /* SWIFT_MAGMA2_HYDRO_PARAMETERS_H */
