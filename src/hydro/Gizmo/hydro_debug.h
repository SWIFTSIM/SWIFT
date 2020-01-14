/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_DEBUG_H
#define SWIFT_GIZMO_HYDRO_DEBUG_H

/* Import the right definition */
#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro_debug.h"
#elif defined(GIZMO_MFM_SPH)
#include "MFM/hydro_debug.h"
#endif

#endif /* SWIFT_GIZMO_HYDRO_DEBUG_H */
