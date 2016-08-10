/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#ifndef SWIFT_HYDRO_GRADIENTS_H
#define SWIFT_HYDRO_GRADIENTS_H

#define SPH_GRADIENTS

#if defined(SPH_GRADIENTS)
#include "hydro_gradients_sph.h"
#elif defined(GIZMO_GRADIENTS)
#include "hydro_gradients_gizmo.h"
#else
/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#include "hydro_gradients_none.h"
#endif

#endif  // SWIFT_HYDRO_GRADIENTS_H
