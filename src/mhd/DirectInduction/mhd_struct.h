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

  /* Predicted magnetic field over density */
  float B_over_rho[3];

  /* Time derivative of magnetic field over density */
  float B_over_rho_dt[3];

  /* Alfven speed (=sqrt(B2/(mu_0 * rho))) of the particle drifted to the
   * current time */
  float Alfven_speed;

  /* Divergence of the magnetic field */
  float divB;

  /* Curl of the magnetic field */
  float curl_B[3];

  /* Tensile instability correction multiplicative prefactor */
  float monopole_beta;

  /* Local plasma beta mean square and its normalising factor */
  float plasma_beta_mean_square;
  float plasma_beta_mean_square_norm;
  
  /* Artifical resistivity multiplicative prefactor */
  float art_diff_beta;

  /* Spatial gradient tensor of the magnetic field */
  float grad_B_tensor[3][3];

  /* Artificial resistivity gradient based switch */
  float alpha_AR;

  /* Artificial resistivity contribution to the time derivative of magnetic
   * field over density */
  float B_over_rho_dt_AR[3];

  /* Artificial resistivity contribution to the time derivative of thermal
   * energy */
  float u_dt_AR;

  /* Physical resistive parameter */
  float resistive_eta;

  /* Predicted Dedner scalar over divergence cleaning speed */
  float psi_over_ch;

  /* Time derivative of Dedner scalar over divergence cleaning speed */
  float psi_over_ch_dt;

  /* SPH <1> error */
  float mean_SPH_err;

  /* SPH <grad1> error */
  float mean_grad_SPH_err[3];

  /* Magnetic force */
  float tot_mag_F[3];

  /* B advection source */
  float Adv_B_source[3];

  /* B total diffusion source */
  float Diff_B_source[3];

  /* Laplacian B */
  float Delta_B[3];
};

/**
 * @brief Particle-carried extra fields for the MHD scheme.
 */
struct mhd_xpart_data {

  /* Full step magnetic field over density */
  float B_over_rho_full[3];

  /* Full step dedner scalar over divergence cleaning speed */
  float psi_over_ch_full;
};

#endif /* SWIFT_DIRECT_INDUCTION_MHD_STRUCT_H */
