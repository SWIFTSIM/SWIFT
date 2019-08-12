/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_FLOOR_PART_GEAR_H
#define SWIFT_PRESSURE_FLOOR_PART_GEAR_H

/**
 * Structure containing the required variables for the pressure
 * floor in the density loop.
 */
struct pressure_floor_part_data {
  /*! Estimation of local turbulence (squared) */
  float sigma2;
};


#endif // SWIFT_PRESSURE_FLOOR_PART_GEAR_H
