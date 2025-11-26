/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025   Jacob Kegerreis (j.kegerreis@imperial.ac.uk).
 *               2025   Thomas Sandnes (thomas.d.sandnes@durham.ac.uk).
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
#ifndef SWIFT_HYDRO_PART_STRENGTH_H
#define SWIFT_HYDRO_PART_STRENGTH_H

/* Config parameters */
#include <config.h>

/* Hydro functions */
#ifdef MATERIAL_STRENGTH

#if defined(PLANETARY_SPH)
#include "./hydro/Planetary/strength/default/hydro_part_strength.h"
#elif defined(REMIX_SPH)
#include "./hydro/REMIX/strength/default/hydro_part_strength.h"
#else
#error "Choice of SPH variant not valid for material strength"
#endif

#else /* !MATERIAL_STRENGTH */

#if defined(PLANETARY_SPH)
#include "./hydro/Planetary/strength/none/hydro_part_strength.h"
#elif defined(REMIX_SPH)
#include "./hydro/REMIX/strength/none/hydro_part_strength.h"
#endif

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_HYDRO_PART_STRENGTH_H */
