/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_LOGGER_CHEMISTRY_H
#define SWIFT_LOGGER_CHEMISTRY_H

/* Config parameters. */
#include "../config.h"

/* Import the right functions */
#if defined(CHEMISTRY_NONE)
#include "./chemistry/none/logger_chemistry.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/GEAR/logger_chemistry.h"
#elif defined(CHEMISTRY_GEAR_DIFFUSION)
#error TODO
#elif defined(CHEMISTRY_QLA)
#error TODO
#elif defined(CHEMISTRY_EAGLE)
#error TODO
#else
#error "Invalid choice of chemistry function."
#endif

#endif /* SWIFT_CHEMISTRY_H */
