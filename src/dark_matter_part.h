/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Camila Correa (camila.correa@cea.fr)
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
#ifndef SWIFT_DARK_MATTER_PART_H
#define SWIFT_DARK_MATTER_PART_H

/* Config parameters. */
#include <config.h>

/* Select the correct sidm model */
#if defined(SIDM_NONE)
#include "./dark_matter/Default/dark_matter_part.h"
#elif defined(SIDM_MODEL)
#include "./dark_matter/TangoSIDM/dark_matter_part.h"
#else
#error "Invalid choice of sidm model"
#endif

#endif /* SWIFT_DARK_MATTER_PART_H */
