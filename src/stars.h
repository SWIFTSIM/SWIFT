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
#ifndef SWIFT_STARS_H
#define SWIFT_STARS_H

/* Config parameters. */
#include "../config.h"

/* Select the correct star model */
#if defined(STARS_NONE)
#include "./stars/Default/stars.h"
#include "./stars/Default/stars_iact.h"
#elif defined(STARS_EAGLE)
#include "./stars/EAGLE/stars.h"
#include "./stars/EAGLE/stars_iact.h"
#else
#error "Invalid choice of star model"
#endif

#endif /* SWIFT_STARS_H */
