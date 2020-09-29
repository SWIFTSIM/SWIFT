/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_PROPERTIES_CONST_DU_H
#define SWIFT_COOLING_PROPERTIES_CONST_DU_H

/**
 * @file src/cooling/const_du/cooling_properties.h
 * @brief Structures related to the "constant cooling" cooling function.
 *
 * This is the simplest possible cooling function. A constant cooling rate
 * (du/dt) with a minimal energy floor is applied. Should be used as a template
 * for more realistic functions.
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling rate in internal units. du_dt = -cooling_rate */
  float cooling_rate;

  /*! Minimally allowed internal energy of the particles */
  float min_energy;

  /*! Constant multiplication factor for time-step criterion */
  float cooling_tstep_mult;
};

#endif /* SWIFT_COOLING_PROPERTIES_CONST_DU_H */
