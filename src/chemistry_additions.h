/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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
#ifndef SWIFT_CHEMISTRY_ADDITIONS_H
#define SWIFT_CHEMISTRY_ADDITIONS_H

/**
 * @file src/chemistry_additions.h
 * @brief Branches between the different additional functions required outside
 * of the chemistry files (e.g. in hydro loops);
 * Required to avoid circular inclusions.
 **/

/* Config parameters. */
#include <config.h>

/* Import the right chemistry definition */
#if defined(CHEMISTRY_NONE)
#include "./chemistry/none/chemistry_additions.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/GEAR/chemistry_additions.h"
#elif defined(CHEMISTRY_GEAR_DIFFUSION)
#include "./chemistry/GEAR_DIFFUSION/chemistry_additions.h"
#elif defined(CHEMISTRY_AGORA)
#include "./chemistry/AGORA/chemistry_additions.h"
#elif defined(CHEMISTRY_QLA)
#include "./chemistry/QLA/chemistry_additions.h"
#elif defined(CHEMISTRY_EAGLE)
#include "./chemistry/EAGLE/chemistry_additions.h"
#else
#error "Invalid choice of chemistry function."
#endif

#endif  // SWIFT_CHEMISTRY_ADDITIONS_H
