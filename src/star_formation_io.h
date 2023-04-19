/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_STAR_FORMATION_IO_H
#define SWIFT_STAR_FORMATION_IO_H

/**
 * @file src/star_formation_io.h
 * @brief Branches between the i/o routines for the SF code.
 */

/* Config parameters. */
#include <config.h>

/* Import the right cooling definition */
#if defined(STAR_FORMATION_NONE)
#include "./star_formation/none/star_formation_io.h"
#elif defined(STAR_FORMATION_QLA)
#include "./star_formation/QLA/star_formation_io.h"
#elif defined(STAR_FORMATION_EAGLE)
#include "./star_formation/EAGLE/star_formation_io.h"
#elif defined(STAR_FORMATION_GEAR)
#include "./star_formation/GEAR/star_formation_io.h"
#else
#error "Invalid choice of star formation model."
#endif

#endif /* SWIFT_STAR_FORMATION_IO_H */
