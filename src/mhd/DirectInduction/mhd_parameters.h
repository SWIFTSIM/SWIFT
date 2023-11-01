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

#ifndef SWIFT_DIRECT_INDUCTION_MHD_PARAMETERS_H
#define SWIFT_DIRECT_INDUCTION_MHD_PARAMETERS_H

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
#include "math.h"

/**
 * @file DirectInduction/mhd_parameters.h
 * @brief NO MHD but default parameters for other schemes
 *
 *        This file defines a number of things that are used in
 *        mhd schemes as defaults for run-time parameters
 *        as well as a number of compile-time parameters.
 */

/* Tensile instability correction factor */

#define mhd_props_tensile_instability_correction_prefactor 1.0f

/* Dedner cleaning -- FIXED -- MUST BE DEFINED AT COMPILE-TIME */

/* Standard Hyperbolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_hyperbolic 1.0f

/* Additional hyperbolic term of Dender Scalar field evolution, propotrional to
 * div(v) */
#define mhd_props_dedner_hyperbolic_divv 0.5f

/* Parabolic term of Dender Scalar field evolution */
#define mhd_propos_dedner_parabolic 1.0f

/* Magnetic Diffusion parameters -- Defaults can be changed in RunTime */
/* Magnetic Diffusion, if set to 0 IDEAL mhd
 *  */
#define mhd_propos_default_resistive_eta 0.0f

/* Artificial diffussion term */
#define mhd_props_artificial_diffusion_beta 1.0f

//#define monopole_beta 1.0f
//#define resistivity_beta 1.0f
//#define dedner_beta 1.0f
//#define dedner_gamma 0.5f

/* Structs that store the relevant variables */

/*! MHD parameters */
struct mhd_global_data {
  /*! For the fixed, simple case of direct induction. */
  float monopole_subtraction;
  float art_diffusion;
  float hyp_dedner;
  float hyp_dedner_divv;
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
  mhd->monopole_subtraction =
      parser_get_param_float(params, "MHD:monopole_subtraction");
  mhd->art_diffusion =
      parser_get_param_float(params, "MHD:artificial_diffusion");
  mhd->hyp_dedner_divv =
      parser_get_param_float(params, "MHD:hyperbolic_dedner_divv");
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

  message("MHD tensile instability correction prefactor: %.3f ",
          mhd->monopole_subtraction);
  message("Artificial diffusion: %.3f ", mhd->art_diffusion);
  message("Dedner Hyperbolic/Hyperbolic div(v)/Parabolic: %.3f, %.3f, %.3f ",
          mhd->hyp_dedner, mhd->hyp_dedner_divv, mhd->par_dedner);
  message("MHD global Resistive Eta: %.3f", mhd->mhd_eta);
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
  io_write_attribute_f(h_grpsph, "MHD Tensile Instability Correction Prefactor",
                       mhd_data->monopole_subtraction);
  io_write_attribute_f(h_grpsph, "Artificial Diffusion Constant",
                       mhd_data->art_diffusion);
  io_write_attribute_f(h_grpsph, "Dedner Hyperbolic Constant",
                       mhd_data->hyp_dedner);
  io_write_attribute_f(h_grpsph, "Dedner Hyperbolic div(v) Constant",
                       mhd_data->hyp_dedner_divv);
  io_write_attribute_f(h_grpsph, "Dedner Parabolic Constant",
                       mhd_data->par_dedner);
  io_write_attribute_f(h_grpsph, "Resistive Eta", mhd_data->mhd_eta);
  // io_write_attribute_f(h_grpsph, "Comoving exponent", mhd_comoving_factor);
}
#endif

#if defined(HAVE_HDF5)
/** XXX TO BE IMPLEMENTED
 * @brief Prints the viscosity information to the snapshot when writing.
 *
 * @param h_grpsph: the SPH group in the ICs to write attributes to.
 * @param viscosity: pointer to the viscosity_global_data struct.
 **/
// static INLINE void viscosity_print_snapshot(
//    hid_t h_grpsph, const struct viscosity_global_data* viscosity) {
//
//  io_write_attribute_f(h_grpsph, "Alpha viscosity", viscosity->alpha);
//  io_write_attribute_f(h_grpsph, "Beta viscosity", const_viscosity_beta);
//}
#endif

#endif /* SWIFT_NONEDIRECT_INDUCTION_PARAMETERS_H */
