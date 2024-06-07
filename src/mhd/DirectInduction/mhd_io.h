/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leideuniv.nl)
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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_IO_H
#define SWIFT_DIRECT_INDUCTION_MHD_IO_H

#include "adiabatic_index.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @return number of fields readed
 */
INLINE static int mhd_read_particles(struct part* parts,
                                     struct io_props* list) {

  list[0] =
      io_make_input_field("MagneticFluxDensities", FLOAT, 3, COMPULSORY,
                          UNIT_CONV_MAGNETIC_FIELD, parts, mhd_data.B_over_rho);

  return 1;
}

INLINE static void convert_B(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {

  ret[0] = xp->mhd_data.B_over_rho_full[0] * p->rho;
  ret[1] = xp->mhd_data.B_over_rho_full[1] * p->rho;
  ret[2] = xp->mhd_data.B_over_rho_full[2] * p->rho;
}

/**
 * @brief Compute the R0 error.
 *
 * This is div(B) * h / B and only triggers if above a fixed signal-to-noise
 * ratio. The noise is estimated from the SPH approximation to grad(1).
 */
INLINE static void calculate_R0(const struct engine* e, const struct part* p,
                                const struct xpart* xp, float* ret) {

  /* Calculate R0 error metric */
  const float B[3] = {xp->mhd_data.B_over_rho_full[0] * p->rho,
                      xp->mhd_data.B_over_rho_full[1] * p->rho,
                      xp->mhd_data.B_over_rho_full[2] * p->rho};
  const float B_abs = sqrtf(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  const float divB_abs = fabsf(p->mhd_data.divB);

  ret[0] = divB_abs * p->h / (B_abs + FLT_MIN);

  /* Now do some filtering based on the noise level */

  /* Strength of noise masking.
   * Default 10.
   * 1 - no mask, 100 - strong masking*/
  const float signal_to_noise = 10;

  /* SPH approximation of grad 1 */
  const float SPH_gr_1[3] = {p->mhd_data.mean_grad_SPH_err[0],
                             p->mhd_data.mean_grad_SPH_err[1],
                             p->mhd_data.mean_grad_SPH_err[2]};

  /* divB_err = |(B * <grad*1>)| */
  const float divB_err_abs =
      fabsf(B[0] * SPH_gr_1[0] + B[1] * SPH_gr_1[1] + B[2] * SPH_gr_1[2]);

  /* Zero the output if less than the signal to noise ratio */
  if (divB_abs < signal_to_noise * divB_err_abs) ret[0] = 0.f;
}

INLINE static void calculate_R1(const struct engine* e, const struct part* p,
                                const struct xpart* xp, float* ret) {

  /* Calculate R1 error metric */
  const float B[3] = {xp->mhd_data.B_over_rho_full[0] * p->rho,
                      xp->mhd_data.B_over_rho_full[1] * p->rho,
                      xp->mhd_data.B_over_rho_full[2] * p->rho};
  const float B_abs = sqrtf(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

  const float Fmag[3] = {p->mhd_data.tot_mag_F[0], p->mhd_data.tot_mag_F[1],
                         p->mhd_data.tot_mag_F[2]};
  const float Fmag_abs =
      sqrtf(Fmag[0] * Fmag[0] + Fmag[1] * Fmag[1] + Fmag[2] * Fmag[2]);

  const float dot_Fmag_B_abs =
      fabsf(B[0] * Fmag[0] + B[1] * Fmag[1] + B[2] * Fmag[2]);

  ret[0] = dot_Fmag_B_abs / (B_abs * Fmag_abs + FLT_MIN);

  /* Now do some filtering based on the noise level */

  /* Strength of noise masking.
   * Default 10.
   * 1 - no mask, 100 - strong masking*/

  const float signal_to_noise = 10;

  /* Difference between exact 1 and SPH approximation of 1*/
  const float SPH_1_diff = 1 - p->mhd_data.mean_SPH_err;

  /* SPH approximation of grad 1 */
  const float SPH_gr_1[3] = {p->mhd_data.mean_grad_SPH_err[0],
                             p->mhd_data.mean_grad_SPH_err[1],
                             p->mhd_data.mean_grad_SPH_err[2]};

  /* Get total force acting on particles*/
  const float Ftot[3] = {p->a_hydro[0], p->a_hydro[1], p->a_hydro[2]};
  const float Ftot_abs =
      sqrtf(Ftot[0] * Ftot[0] + Ftot[1] * Ftot[1] + Ftot[2] * Ftot[2]);

  /* Relative contribution of magnetic force to the total force */
  const float Fmag_fraction = Fmag_abs / (Ftot_abs + Fmag_abs + FLT_MIN);

  /* Get vacuum permeability */
  const float mu0 = e->physical_constants->const_vacuum_permeability;

  /* Estimate noise level in Fmag from SPH aproximation of gradients */
  const float two_Pmag_over_rho = B_abs * B_abs / (p->rho * mu0 + FLT_MIN);
  const float Fmag_SPH_gr_err[3] = {fabsf(two_Pmag_over_rho * SPH_gr_1[0]),
                                    fabsf(two_Pmag_over_rho * SPH_gr_1[1]),
                                    fabsf(two_Pmag_over_rho * SPH_gr_1[2])};

  /* Estimate noise level in Fmag from SPH sums */
  const float Fmag_SPH_1_err[3] = {fabsf(SPH_1_diff * Fmag[0]),
                                   fabsf(SPH_1_diff * Fmag[1]),
                                   fabsf(SPH_1_diff * Fmag[2])};

  /* Total Fmag error estimate*/
  const float Fmag_err[3] = {Fmag_SPH_gr_err[0] + Fmag_SPH_1_err[0],
                             Fmag_SPH_gr_err[1] + Fmag_SPH_1_err[1],
                             Fmag_SPH_gr_err[2] + Fmag_SPH_1_err[2]};
  const float Fmag_err_abs =
      sqrtf(Fmag_err[0] * Fmag_err[0] + Fmag_err[1] * Fmag_err[1] +
            Fmag_err[2] * Fmag_err[2]);

  /* Zero the output if less than the signal to noise ratio or if magnetic force
   * contribution is less than 10%  */
  if (Fmag_abs < signal_to_noise * Fmag_err_abs || Fmag_fraction < 0.1) {
    ret[0] = 0.f;
  }
}

INLINE static void calculate_R2(const struct engine* e, const struct part* p,
                                const struct xpart* xp, float* ret) {

  /* Calculate R2 error metric */
  const float curlB[3] = {p->mhd_data.curl_B[0], p->mhd_data.curl_B[1],
                          p->mhd_data.curl_B[2]};
  const float curlB_abs =
      sqrtf(curlB[0] * curlB[0] + curlB[1] * curlB[1] + curlB[2] * curlB[2]);

  const float divB_abs = fabsf(p->mhd_data.divB);

  ret[0] = divB_abs / (curlB_abs + FLT_MIN);

  /* Now do some filtering based on the noise level */

  /* Strength of noise masking.
   * Default 10.
   * 1 - no mask, 100 - strong masking*/
  const float signal_to_noise = 10;

  const float B[3] = {xp->mhd_data.B_over_rho_full[0] * p->rho,
                      xp->mhd_data.B_over_rho_full[1] * p->rho,
                      xp->mhd_data.B_over_rho_full[2] * p->rho};

  /* SPH approximation of grad 1 */
  const float SPH_gr_1[3] = {p->mhd_data.mean_grad_SPH_err[0],
                             p->mhd_data.mean_grad_SPH_err[1],
                             p->mhd_data.mean_grad_SPH_err[2]};

  /* divB_err = |(B * <grad*1>)| */
  const float divB_err_abs =
      fabsf(B[0] * SPH_gr_1[0] + B[1] * SPH_gr_1[1] + B[2] * SPH_gr_1[2]);

  /* curlB_err =|[B x <grad*1>]| */
  const float curlB_err[3] = {B[1] * SPH_gr_1[2] - B[2] * SPH_gr_1[1],
                              B[2] * SPH_gr_1[0] - B[0] * SPH_gr_1[2],
                              B[0] * SPH_gr_1[1] - B[1] * SPH_gr_1[0]};
  const float curlB_err_abs =
      sqrtf(curlB_err[0] * curlB_err[0] + curlB_err[1] * curlB_err[1] +
            curlB_err[2] * curlB_err[2]);

  /* Zero the output if less than the signal to noise ratio */
  if (divB_abs < signal_to_noise * divB_err_abs ||
      curlB_abs < signal_to_noise * curlB_err_abs) {
    ret[0] = 0.f;
  }
}

INLINE static void calculate_R3(const struct engine* e, const struct part* p,
                                const struct xpart* xp, float* ret) {

  /* Calculate R3 error metric */

  const float B[3] = {xp->mhd_data.B_over_rho_full[0] * p->rho,
                      xp->mhd_data.B_over_rho_full[1] * p->rho,
                      xp->mhd_data.B_over_rho_full[2] * p->rho};
  const float B_abs = sqrtf(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

  const float curlB[3] = {p->mhd_data.curl_B[0], p->mhd_data.curl_B[1],
                          p->mhd_data.curl_B[2]};
  const float curlB_abs =
      sqrtf(curlB[0] * curlB[0] + curlB[1] * curlB[1] + curlB[2] * curlB[2]);

  ret[0] = curlB_abs * p->h / (B_abs + FLT_MIN);

  /* Now do some filtering based on the noise level */

  /* Strength of noise masking.
   * Default 10.
   * 1 - no mask, 100 - strong masking*/
  const float signal_to_noise = 10;

  /* SPH approximation of grad 1 */
  const float SPH_gr_1[3] = {p->mhd_data.mean_grad_SPH_err[0],
                             p->mhd_data.mean_grad_SPH_err[1],
                             p->mhd_data.mean_grad_SPH_err[2]};

  /* curlB_err =|[B x <grad*1>]| */
  const float curlB_err[3] = {B[1] * SPH_gr_1[2] - B[2] * SPH_gr_1[1],
                              B[2] * SPH_gr_1[0] - B[0] * SPH_gr_1[2],
                              B[0] * SPH_gr_1[1] - B[1] * SPH_gr_1[0]};
  const float curlB_err_abs =
      sqrtf(curlB_err[0] * curlB_err[0] + curlB_err[1] * curlB_err[1] +
            curlB_err[2] * curlB_err[2]);

  /* Zero the output if less than the signal to noise ratio */
  if (curlB_abs < signal_to_noise * curlB_err_abs) {
    ret[0] = 0.f;
  }
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @return num_fields The number of i/o fields to write.
 */
INLINE static int mhd_write_particles(const struct part* parts,
                                      const struct xpart* xparts,
                                      struct io_props* list) {

  list[0] = io_make_output_field_convert_part(
      "MagneticFluxDensities", FLOAT, 3, UNIT_CONV_MAGNETIC_FIELD,
      -1.5f * hydro_gamma, parts, xparts, convert_B,
      "Magnetic flux densities of the particles");

  list[1] = io_make_output_field(
      "MagneticDivergences", FLOAT, 1, UNIT_CONV_MAGNETIC_DIVERGENCE,
      -1.5f * hydro_gamma - 1.f, parts, mhd_data.divB,
      "co-moving DivB  of the particle");

  list[2] = io_make_output_field(
      "DednerScalars", FLOAT, 1, UNIT_CONV_ELECTRIC_CHARGE_FIELD_STRENGTH,
      -1.5f * hydro_gamma - 1.f, parts, mhd_data.psi_over_ch,
      "Dedner scalar associated to the particle");

  list[3] = io_make_output_field(
      "DednerScalarsdt", FLOAT, 1,
      UNIT_CONV_ELECTRIC_CHARGE_FIELD_STRENGTH_RATE, 1.f, parts,
      mhd_data.psi_over_ch_dt,
      "Time derivative of Dedner scalar associated to particle");

  list[4] = io_make_output_field(
      "MagneticFluxDensitiesdt", FLOAT, 3, UNIT_CONV_MAGNETIC_FIELD_PER_TIME,
      1.f, parts, mhd_data.B_over_rho_dt,
      "Time derivative of Magnetic flux densities of the particles");

  list[5] = io_make_output_field(
      "MagneticFluxCurl", FLOAT, 3, UNIT_CONV_MAGNETIC_CURL,
      -1.5f * hydro_gamma - 1.f, parts, mhd_data.curl_B,
      "The curl of Magnetic flux densities of the particles");

  list[6] = io_make_output_field(
      "AlphaAR", FLOAT, 1, UNIT_CONV_NO_UNITS, 1.f, parts, mhd_data.alpha_AR,
      "Artificial resistivity switch of the particles");

  list[7] = io_make_output_field(
      "MagneticFluxDensitiesdtAR", FLOAT, 3, UNIT_CONV_MAGNETIC_FIELD_PER_TIME,
      0.f, parts, mhd_data.B_over_rho_dt_AR,
      "AR contribution to time derivative of Magnetic flux densities of the particles");

  list[8] = io_make_output_field(
      "ThermalEnergiesdtAR", FLOAT, 3, UNIT_CONV_ENERGY_PER_UNIT_MASS_PER_TIME,
      0.f, parts, mhd_data.u_dt_AR,
      "AR contribution to time derivative of thermal energies of the particles");
  
  /* Error metrics */
  list[9] = io_make_output_field_convert_part(
      "R0", FLOAT, 1, UNIT_CONV_NO_UNITS, 0, parts, xparts, calculate_R0,
      "Classical error metric, indicates places with large divergence. "
      "Sensetivity to particle noise depends on signal_to_noise parameter, "
      "default is 10 (if 1 - weak noise filtering, if 100 - strong noise "
      "filtering)");
  list[10] = io_make_output_field_convert_part(
      "R1", FLOAT, 1, UNIT_CONV_NO_UNITS, 0, parts, xparts, calculate_R1,
      "Error metric, angle between B field and total Fmag. Indicates unpysical "
      "magnetic force. Sensetivity to particle noise depends on "
      "signal_to_noise parameter, default is 10 (if 1 - weak noise filtering, "
      "if 100 - strong noise filtering)");
  list[11] = io_make_output_field_convert_part(
      "R2", FLOAT, 1, UNIT_CONV_NO_UNITS, 0, parts, xparts, calculate_R2,
      "Error metric, ratio of divB and |curlB|. Estimates upper limit on "
      "B_monopole/B_physical. Sensetivity to particle noise depends on "
      "signal_to_noise parameter, default is 10 (if 1 - weak noise filtering, "
      "if 100 - strong noise filtering)");
  list[12] = io_make_output_field_convert_part(
      "R3", FLOAT, 1, UNIT_CONV_NO_UNITS, 0, parts, xparts, calculate_R3,
      "Error metric, shows relation of smoothing length to characteristic B "
      "gradient scale. Sensetivity to particle noise depends on "
      "signal_to_noise parameter, default is 10 (if 1 - weak noise filtering, "
      "if 100 - strong noise filtering)");

  return 13;
}

/**
 * @brief Writes the current model of MHD to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void mhd_write_flavour(hid_t h_grpsph) {

  io_write_attribute_s(
      h_grpsph, "MHD Flavour",
      "Orestis - Direct Induction, divB Subtraction, "
      "Artificial Resistivity & Dedner Cleaning. Price et al. (2018).");
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_IO_H */
