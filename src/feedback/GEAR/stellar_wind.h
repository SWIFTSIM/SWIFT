/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Zachary Rey (zachary.rey@epfl.ch)
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
#ifndef SWIFT_GEAR_STELLAR_WIND_H
#define SWIFT_GEAR_STELLAR_WIND_H

#include "hdf5_functions.h"
#include "interpolation.h"
#include "stellar_evolution_struct.h"

/* The cut off parameter for the ejected mass function.
 * In ratio of log_10(M/M_sol)
 */
static const float stellar_wind_x0 = 2;

float calculate_b_parameter(const float log_Z, const float a[]);

float stellar_wind_get_ejected_mass(const float log_Z,const float log_m);

float stellar_wind_get_wind_velocity(const float log_Z,const float log_m);

float stellar_wind_get_energy_dot(const float mass_loss, const float v_infinity);

#endif 