/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_STRUCT_COMPTON_H
#define SWIFT_COOLING_STRUCT_COMPTON_H

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Compton rate in cgs [g cm^2 s^-3 K^-1] */
  double const_Compton_rate_cgs;

  /*! Temperature of the CMB at redshift 0 in cgs [K] */
  double const_T_CMB_0;

  /*! Conversion factor from internal units to cgs for density */
  double conv_factor_density_to_cgs;

  /*! Conversion factor from internal units to cgs for internal energy */
  double conv_factor_energy_to_cgs;

  /*! Conversion factor from internal units from cgs for internal energy
   * derivative */
  double conv_factor_energy_rate_from_cgs;

  /*! Inverse of the proton mass in cgs units [g^-1] */
  double proton_mass_cgs_inv;
};

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct cooling_part_data {};

/**
 * @brief Properties of the cooling stored in the particle data.
 */
struct cooling_xpart_data {

  /*! Energy radiated away by this particle since the start of the run */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_COMPTON_H */
