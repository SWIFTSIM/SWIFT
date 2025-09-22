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


/*! Alpha viscosity, usually =1.0. For lower N_ngb, should be higher */
#define const_viscosity_alpha 2.0

/*! Alpha conductivity, usually =0.05. At lower N_ngb, should be higher */
#define const_conductivity_alpha 0.075

/*! Desired number of neighbours -- CRITICAL that this matches hydro props */
#if defined(HYDRO_DIMENSION_1D)
#define const_kernel_target_neighbours 8.0
#elif defined(HYDRO_DIMENSION_2D)
#define const_kernel_target_neighbours 34.0
#else
#define const_kernel_target_neighbours 114.0
#endif


/* ---------- These parameters should not be changed ---------- */

/*! Use a Swift-like estimator for dh/dt rather than the correct formula
 * 0 = Simple mass flow estimator
 * 1 = Correct formula based on number density constraint
 * 2 = Using v_ij dot G_ij with simple mass flow estimator
 */
#define hydro_props_dh_dt_estimator_type 0

/*! Flag to use Balsara limiter */
#define hydro_props_use_balsara_limiter

/*! Flag to use additional slope limiting procedures */
//#define hydro_props_use_extra_slope_limiter

/* Flag to disallow sign flip in reconstructed quantities */
//#define hydro_props_use_strict_minmod_limiter

/* Slope limiter length, fraction of max. distance in kernel */
#ifdef hydro_props_use_extra_slope_limiter
#define const_grad_overshoot_length 0.25

/*! Slope limiter tolerance */
#define const_grad_overshoot_tolerance 0.0
#endif

/* Viscosity floor when G_ij is extremely misaligned with dx_ij */
#define const_viscosity_cosine_limit 0.1736

/* Viscosity weighting scheme: 
 *    0 = (rho_i * q_i + rho_j * q_j) / (rho_i * rho_j)
 *    1 = (rho_i * q_ij + rho_j * q_ij) / (rho_i * rho_j)
 *    2 = 2.0 * q_ij / (rho_i + rho_j) */
#define hydro_props_viscosity_weighting_type 2

/* Flag to use radial gradients for viscosity and conductivity */
//#define hydro_props_use_radial_artificial_terms

/*! Use the correction terms to make the internal energy match the mass flux */
//#define hydro_props_use_adiabatic_correction

/* Kernel gradient weighting scheme:
 *    0 = 0.5 * (G_i + G_j)
 *    1 = 0.5 * (f_i * G_i + f_j * G_j)
 *    2 = 0.5 * f_ij * (G_i + G_j)
 *        with f_ij = 0.5 * (f_i + f_j)
 *    3 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = 2 * f_i * f_j / (f_i + f_j)
 *    4 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = sqrt(f_i * f_j)
 *    5 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = (f_i * rho_i + f_j * rho_j) / (rho_i + rho_j)
 *    6 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = 2 * f_i * f_j / (f_i * rho_i + f_j * rho_j)
 *    7 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = (f_i * P_i + f_j * P_j) / (P_i + P_j)
 *    8 = 0.5 * f_ij * (G_i + G_j) 
 *        with f_ij = 2 * f_i * f_j / (f_i * P_i + f_j * P_j)
 *   9 = 0.5 * f_ij * (G_i + G_j)
 *       with f_ij = (f_i * V_i + f_j * V_j) / (V_i + V_j)
 */
#define hydro_props_kernel_gradient_weighting 0

/*! Use double precision for all matrix/vector operations */
//#define hydro_props_use_double_precision

#ifdef hydro_props_use_double_precision
/*! Consider matrix inversion to be ill-conditioned above this limit */
#define const_condition_number_upper_limit 300.
/*! Mean interparticle spacing for this kernel and neighbour number */
#define const_kernel_mean_spacing (kernel_gamma*pow(4. * M_PI / (3. * \
    (double)const_kernel_target_neighbours), 1. / 3.))
#else
/*! Consider matrix inversion to be ill-conditioned above this limit */
#define const_condition_number_upper_limit 60.
/*! Mean interparticle spacing for this kernel and neighbour number */
#define const_kernel_mean_spacing (kernel_gamma*powf(4. * M_PI / (3. * \
    (float)const_kernel_target_neighbours), 1. / 3.))
#endif

/*! eta_crit Rosswog 2020 Eq 23. Of order the mean interparticle spacing. */
#define const_slope_limiter_eta_crit (const_kernel_mean_spacing)

/*! eta_fold from Frontiere+'17 Equation 51 */
#define const_slope_limiter_eta_fold 0.2

/*! Softening squared (epsilon^2) in Eq. 15 Rosswog 2020 */
#define const_viscosity_epsilon2 0.01

/*! Cosmology default const_viscosity_beta=2*const_viscosity_alpha
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta (2.0*const_viscosity_alpha)

/*! Prefactor for alpha term in signal velocity */
#define const_viscosity_alpha_prefactor (1.25 * (1. + \
  0.75 * const_viscosity_alpha))

/*! Prefactor for beta term in signal velocity */
#define const_viscosity_beta_prefactor (1.25 * 0.75 * const_viscosity_beta)

/*! Fallback multiplier for alpha/beta terms to reduce spread */
#define const_fallback_reduction_factor 0.25

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
