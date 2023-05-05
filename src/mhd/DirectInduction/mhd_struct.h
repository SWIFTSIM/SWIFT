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

  float B_over_rho_dt[3];

  float v_fm;

  float B_mon;

  float curl_B[3];

  /*union{
     float B_mon;
     float divB;
  };
  */

  float psi_over_ch;

  float psi_over_ch_dt;
};

/**
 * @brief Particle-carried extra fields for the MHD scheme.
 */
struct mhd_xpart_data {

  float B_over_rho_full[3];
};

#endif /* SWIFT_DIRECT_INDUCTION_MHD_STRUCT_H */
