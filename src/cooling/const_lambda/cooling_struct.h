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

  /*! Cooling rate in cgs units. Defined by 'rho * du/dt = -lambda * n_H^2'*/
  float lambda;

  /*! Minimum temperature (in Kelvin) for all gas particles*/
  float min_temperature;

  /*! Fraction of gas mass that is Hydrogen. Used to calculate n_H*/
  float hydrogen_mass_abundance;

  /* 'mu', used to convert min_temperature to min_internal energy*/
  float mean_molecular_weight;

  /*! Minimally allowed internal energy of the particles */
  float min_energy;
  float min_energy_cgs;

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
