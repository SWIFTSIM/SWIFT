/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Stefan Arridge  (stefan.arridge@durham.ac.uk)
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
#ifndef SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H
#define SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H

/**
 * @file src/cooling/const_lambda/cooling_properties.h
 * @brief Structures related to the "constant lambda" cooling function.
 *
 * This model assumes a constant cooling rate Lambda irrespective of redshift
 * or density.
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling rate / nH^2 in physical cgs units [erg * s^-1 * cm^3] */
  double lambda_nH2_cgs;

  /*! Conversion factor from internal units to cgs for density */
  double conv_factor_density_to_cgs;

  /*! Conversion factor from internal units from cgs for internal energy
   * derivative */
  double conv_factor_energy_rate_from_cgs;

  /*! Inverse of the proton mass in cgs units [g^-1] */
  double proton_mass_cgs_inv;

  /*! Constant multiplication factor for time-step criterion */
  float cooling_tstep_mult;
};

#endif /* SWIFT_COOLING_PROPERTIES_CONST_LAMBDA_H */
