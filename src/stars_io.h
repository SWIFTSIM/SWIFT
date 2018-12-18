/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_STARS_IO_H
#define SWIFT_STARS_IO_H

#include "../config.h"
#include "./const.h"

/* Load the correct star type */
#if defined(STARS_NONE)
#include "./stars/Default/stars_io.h"
#elif defined(STARS_EAGLE)
#include "./stars/EAGLE/stars_io.h"
#elif defined(STARS_GEAR)
#include "./stars/GEAR/stars_io.h"
#else
#error "Invalid choice of star model"
#endif

#endif /* SWIFT_STARS_IO_H */
