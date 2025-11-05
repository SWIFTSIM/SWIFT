/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_MHD_PARAMETERS_H
#define SWIFT_NONE_MHD_PARAMETERS_H

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
 * @file None/mhd_parameters.h
 * @brief NO MHD but default parameters for other schemes
 *
 *        This file defines a number of things that are used in
 *        mhd schemes as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */

/* Dedner cleaning -- FIXED -- MUST BE DEFINED AT COMPILE-TIME */

/* if set to 0 NO dedner cleaning
 * hyperbolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_hyperbolic 0.0f

/*
 * parabolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_parabolic 0.0f

/* Magnetic Diffusion parameters -- Defaults can be changed in RunTime */

/* Magnetic Diffusion, if set to 0 IDEAL mhd
 *  */
#define mhd_propos_default_difussion_eta 0.0f

/*! MHD parameters */
struct mhd_global_data {};

/* Functions for reading from parameter file */

/**
 * @brief Initialises the mhd parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param mhd: pointer to the mhd_global_data struct to be filled.
 **/
static INLINE void mhd_init(struct swift_params* params,
                            const struct unit_system* us,
                            const struct phys_const* phys_const,
                            struct mhd_global_data* mhd) {}

/**
 * @brief Prints out the mhd parameters at the start of a run.
 *
 * @param mhd: pointer to the mhd_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void mhd_print(const struct mhd_global_data* mhd) {}

#endif /* SWIFT_NONE_MHD_PARAMETERS_H */
