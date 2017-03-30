/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#ifndef SWIFT_VORONOI_CELL_H
#define SWIFT_VORONOI_CELL_H

#if defined(HYDRO_DIMENSION_1D)
#include "voronoi1d_cell.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "voronoi2d_cell.h"
#elif defined(HYDRO_DIMENSION_3D)
#include "voronoi3d_cell.h"
#else
#error "You have to select a dimension for the hydro!"
#endif

#endif  // SWIFT_VORONOI_CELL_H
