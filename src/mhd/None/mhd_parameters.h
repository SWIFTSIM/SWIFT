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
 */

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

#if defined(HAVE_HDF5)
/**
 * @brief Prints the MHD information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param mhd_data: pointer to the mhd_global_data struct.
 **/
static INLINE void mhd_print_snapshot(hid_t h_grpsph,
                                      const struct mhd_global_data* mhd_data) {}
#endif

#endif /* SWIFT_NONE_MHD_PARAMETERS_H */
