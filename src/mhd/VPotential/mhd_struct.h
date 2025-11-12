/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
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
#ifndef SWIFT_VECTOR_POTENTIAL_MHD_STRUCT_H
#define SWIFT_VECTOR_POTENTIAL_MHD_STRUCT_H

/**
 * @brief Particle-carried fields for the MHD scheme.
 */
#include "symmetric_matrix.h"

struct mhd_part_data {

  /*! Predicted Bfield */
  float BPred[3];

  /*! Predicted BSmooth */
  float BSmooth[3];

  /* Alfven speed (=sqrt(B2/(mu_0 * rho))) of the particle drifted to the
   * current time */
  float Alfven_speed;

  /*! Full step Divergence of B */
  float divB;

  /* Curl of the magnetic field */
  float curl_B[3];

  /* predicted VPotencial */
  float APred[3];

  /* predicted step Gauge, divA */
  float Gau, divA;

  /* Time derivative of Gauge */
  float Gau_dt;

  /* Spatial gradient tensor of the magnetic field */
  float grad_B_tensor[3][3];

  /* VP evolution */
  float dAdt[3];

  /* Artificial resistivity gradient based switch */
  float alpha_AR;

  /* mute variable */
  float Q0;

  /* Resistive Eta */
  float resistive_eta;

  /* SPH <1> error */
  float mean_SPH_err;

  /* SPH <grad1> error */
  float mean_grad_SPH_err[3];

  /* Magnetic force */
  float tot_mag_F[3];

  /* A advection source */
  float Adv_A_source[3];

  /* B total diffusion source */
  float Diff_A_source[3];

  /* Laplacian A */
  float Delta_A[3];

  struct {

    /*! The inverse of 'correction matrix' (e.q. 6) - It's symmetric */
    struct sym_matrix c_matrix_inv;

    /*! Gradient per component of the Afield means Bfield*/
    float Mat_b[3][3];
    /*! Gradient per component of the dAdt*/
    float Mat_da[3][3];

  } grad;

  struct {

    /*! The 'correction matrix' (e.q. 6) - It's symmetric */
    struct sym_matrix c_matrix;

    /*! Gradient per component of the Afield means Bfield*/
    float Mat_b[3][3];

  } force;
};

/**
 * @brief Particle-carried extra fields for the MHD scheme.
 */
struct mhd_xpart_data {
  /* Full step magnetic vector potential */
  float Afull[3];

  /* Full step dedner scalar over divergence cleaning speed */
  float Gaufull;
};

#endif /* SWIFT_VECTOR_POTENTIAL_MHD_STRUCT_H */
