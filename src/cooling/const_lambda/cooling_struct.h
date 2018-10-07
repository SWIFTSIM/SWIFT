/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Richard Bower (r.g.bower@durham.ac.uk)
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

#ifndef SWIFT_COOLING_STRUCT_CONST_LAMBDA_H
#define SWIFT_COOLING_STRUCT_CONST_LAMBDA_H

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling rate in cgs units */
  double lambda_cgs;

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

/**
 * @brief Properties of the cooling stored in the particle data.
 */
struct cooling_xpart_data {

  /*! Energy radiated away by this particle since the start of the run */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_CONST_LAMBDA_H */
