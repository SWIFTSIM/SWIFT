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
#ifndef SWIFT_DI_MHD_STRUCT_H
#define SWIFT_DI_MHD_STRUCT_H

/**
 * @brief Particle-carried fields for the MHD scheme.
 */
struct mhd_part_data {

  /*! Predicted Bfield */
  float BPred[3];
  /*! Full step Divergence of B */
  float divB;
  /*! limiters Induction and force */
  float Q1, Q0;
  /*! dB Direct Induction */
  float dBdt[3];
  /* Predicted Dedner Cleaning Scalar */
  float phi;
  /* Resistive Eta */
  float resistive_eta;
  /* Curl B */
  float curlB[3];
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

  /*! Full Step Magnetic field */
  float Bfld_full[3];
};

#endif /* SWIFT_DI_MHD_STRUCT_H */