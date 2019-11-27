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
#ifndef SWIFT_GIZMO_HYDRO_PARAMETERS_H
#define SWIFT_GIZMO_HYDRO_PARAMETERS_H

/**
 * @file src/hydro_parameters.h
 * @brief Contains all the parameters of the hydro schemes, included from
 *        their own local header.
 */

/* Import the right hydro header */
#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro_parameters.h"
#elif defined(GIZMO_MFM_SPH)
#include "MFM/hydro_parameters.h"
#endif

#endif /* SWIFT_HYDRO_PARAMETERS_H */
