/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_STRUCT_H
#define SWIFT_DIRECT_INDUCTION_MHD_STRUCT_H

/**
 * @brief Particle-carried fields for the MHD scheme.
 */
struct mhd_part_data {

  float B_over_rho[3];

  float divB;

  /*! dB Direct Induction */
  float B_over_rho_dt[3];

  float v_fm;

  float curl_B[3];

  /* Resistive Eta */
  float resistive_eta;

  float grad_B_tensor[3][3];

  float alpha_AR;

  float psi_over_ch;

  float psi_over_ch_dt;

  /*! Monopole subtraction in Lorentz Force*/
  float monopole_beta;

  /*! Artifical Diffusion */
  float art_diff_beta;

  float B_over_rho_dt_AR[3];

  float u_dt_AR;

  float v_sig_AR;
  
  /* SPH <1> error */
  float mean_SPH_err;

  /* SPH <grad1> error */
  float mean_grad_SPH_err[3];

  /* Magnetic force */
  float tot_mag_F[3];
};

/**
 * @brief Particle-carried extra fields for the MHD scheme.
 */
struct mhd_xpart_data {

  /*! Full Step Magnetic field */
  float B_over_rho_full[3];
};

#endif /* SWIFT_DIRECT_INDUCTION_MHD_STRUCT_H */
