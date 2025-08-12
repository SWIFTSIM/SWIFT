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
#ifndef SWIFT_STRENGTH_H
#define SWIFT_STRENGTH_H

/* Config parameters */
#include <config.h>

/* Strength functions */
#ifdef MATERIAL_STRENGTH
// Utilities
#include "./strength/strength_utilities.h"

// Strength models
#include "./strength/strength_artificial_stress.h"
#include "./strength/strength_damage.h"
#include "./strength/strength_yield_stress_softening.h"
#include "./strength/strength_yield_stress.h"

// Stress tensor, combining strength models
#include "./strength/strength_stress_tensor.h"
#endif /* MATERIAL_STRENGTH */

/* Hydro functions */
#ifdef MATERIAL_STRENGTH

#if defined(PLANETARY_SPH)
#include "./hydro/Planetary/strength/default/hydro_strength.h"
#include "./hydro/Planetary/strength/default/hydro_iact_strength.h"
#include "./hydro/Planetary/strength/default/hydro_io_strength.h"
#elif defined(REMIX_SPH)
#include "./hydro/REMIX/strength/default/hydro_strength.h"
#include "./hydro/REMIX/strength/default/hydro_iact_strength.h"
#include "./hydro/REMIX/strength/default/hydro_io_strength.h"
#else
#error "Choice of SPH variant not valid for material strength"
#endif

#else /* !MATERIAL_STRENGTH */

#if defined(PLANETARY_SPH)
#include "./hydro/Planetary/strength/none/hydro_strength.h"
#include "./hydro/Planetary/strength/none/hydro_iact_strength.h"
#include "./hydro/Planetary/strength/none/hydro_io_strength.h"
#elif defined(REMIX_SPH)
#include "./hydro/REMIX/strength/none/hydro_strength.h"
#include "./hydro/REMIX/strength/none/hydro_iact_strength.h"
#include "./hydro/REMIX/strength/none/hydro_io_strength.h"
#endif

#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_STRENGTH_H */
