/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_FVPM_GEOMETRY_STRUCT_H
#define SWIFT_FVPM_GEOMETRY_STRUCT_H

/* Config parameters. */
#include <config.h>

/* Import the right geometry struct definition */
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH) || defined(RT_GEAR)
#include "./fvpm_geometry/Gizmo/fvpm_geometry_struct.h"
#else
#include "./fvpm_geometry/None/fvpm_geometry_struct.h"
#endif

#endif /* SWIFT_FVPM_GEOMETRY_STRUCT_H */
