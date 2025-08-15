/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_SIDM_IO_H
#define SWIFT_SIDM_IO_H

#include <config.h>

/* Local includes */
#include "engine.h"

/* Load the correct SIDM model */
#if defined(SIDM_NONE)
#include "./sidm/None/sidm_io.h"
#elif defined(SIDM_BASIC)
#include "./sidm/Basic/sidm_io.h"
#else
#error "Invalid choice of SIDM model"
#endif

#endif /* SWIFT_SIDM_IO_H */