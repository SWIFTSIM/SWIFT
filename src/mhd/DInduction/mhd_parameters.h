/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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

#ifndef SWIFT_DI_MHD_PARAMETERS_H
#define SWIFT_DI_MHD_PARAMETERS_H

/* Configuration file */
#include "config.h"

#include <math.h>

/* Global headers */
#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local headers */
#include "adiabatic_index.h"
#include "common_io.h"
#include "error.h"
#include "inline.h"

/**
 * @file DInduction/mhd_parameters.h
 * @brief NO MHD but default parameters for other schemes
 *
 *        This file defines a number of things that are used in
 *        mhd schemes as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */

/* Freedom to choose the way the Comoving Bfield behaves
 * the comoving conversion goes like:
 * B_phi = a^MHD_COMOVING_FACTOR * B_co
 */
#define mhd_comoving_factor -2.f
//#define mhd_comoving_factor -3.f/2.f*(hydro_gamma-1.f)

/* if set to 0 NO dedner cleaning
 * hyperbolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_hyperbolic 1.0f

/*
 * parabolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_parabolic 2.f

/* Magnetic Diffusion parameters -- Defaults can be changed in RunTime */
/* Magnetic Diffusion, if set to 0 IDEAL mhd
 *  */
#define mhd_propos_default_resistive_eta 0.0f

/* Structs that store the relevant variables */

/*! MHD parameters */
struct mhd_global_data {
  /*! For the fixed, simple case of direct induction. */
  float hyp_dedner;
  float par_dedner;
  float mhd_eta;
};

/* Functions for reading from parameter file */

/* Forward declartions */
// struct swift_params;
// struct phys_const;
// struct unit_system;

/* Viscosity */

/**
 * @brief Initialises the viscosity parameters in the struct from
 *        the parameter file, or sets them to defaults.
 *
 * @param params: the pointer to the swift_params file
 * @param us: pointer to the internal unit system
 * @param phys_const: pointer to the physical constants system
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
static INLINE void mhd_init(struct swift_params* params,
                            const struct unit_system* us,
                            const struct phys_const* phys_const,
                            struct mhd_global_data* mhd) {

  /* Read the mhd parameters from the file, if they exist,
   * otherwise set them to the defaults defined above. */

  mhd->hyp_dedner = parser_get_param_float(params, "MHD:hyperbolic_dedner");
  mhd->par_dedner = parser_get_param_float(params, "MHD:parabolic_dedner");
  mhd->mhd_eta = parser_get_param_float(params, "MHD:resistive_eta");
}

/**
 * @brief Initialises a viscosity struct to sensible numbers for mocking
 *        purposes.
 *
 * @param viscosity: pointer to the viscosity_global_data struct to be filled.
 **/
// static INLINE void viscosity_init_no_hydro(
//    struct viscosity_global_data* viscosity) {
//  viscosity->alpha = hydro_props_default_viscosity_alpha;
//}

/**
 * @brief Prints out the viscosity parameters at the start of a run.
 *
 * @param viscosity: pointer to the viscosity_global_data struct found in
 *                   hydro_properties
 **/
static INLINE void mhd_print(const struct mhd_global_data* mhd) {

  message("Dedner Hyperbolic/Parabolic: %.3f, %.3f ", mhd->hyp_dedner,
          mhd->par_dedner);
  message("MHD global dissipation Eta: %.3f", mhd->mhd_eta);
}

#if defined(HAVE_HDF5)
/**
 * @brief Prints the MHD information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param mhd_data: pointer to the mhd_global_data struct.
 **/
static INLINE void mhd_print_snapshot(hid_t h_grpsph,
                                      const struct mhd_global_data* mhd_data) {

  io_write_attribute_f(h_grpsph, "Dedner Hyperbolic Constant",
                       mhd_data->hyp_dedner);
  io_write_attribute_f(h_grpsph, "Dedner Parabolic Constant",
                       mhd_data->par_dedner);
  io_write_attribute_f(h_grpsph, "Resistive Eta", mhd_data->mhd_eta);
  io_write_attribute_f(h_grpsph, "Comoving exponent", mhd_comoving_factor);
}
#endif

#endif /* SWIFT_DI_MHD_PARAMETERS_H */
