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
#ifndef SWIFT_ENTROPY_FLOOR_H
#define SWIFT_ENTROPY_FLOOR_H

/**
 * @file src/entropy_floor.h
 * @brief Branches between the different entropy floor models
 */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "common_io.h"
#include "error.h"
#include "inline.h"

/* Pre-declarations to avoid cyclic inclusions */
static INLINE float hydro_get_physical_density(const struct part *restrict p,
                                               const struct cosmology *cosmo);
static INLINE float hydro_get_comoving_density(const struct part *restrict p);

/* Import the right entropy floor definition */
#if defined(ENTROPY_FLOOR_NONE)
#include "./entropy_floor/none/entropy_floor.h"
#elif defined(ENTROPY_FLOOR_QLA)
#include "./entropy_floor/QLA/entropy_floor.h"
#elif defined(ENTROPY_FLOOR_EAGLE)
#include "./entropy_floor/EAGLE/entropy_floor.h"
#else
#error "Invalid choice of entropy floor"
#endif

#endif /* SWIFT_ENTROPY_FLOOR_H */
