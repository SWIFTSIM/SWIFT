/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_GRAVITY_IO_H
#define SWIFT_GRAVITY_IO_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "./const.h"

/* Import the right functions */
#if defined(DEFAULT_GRAVITY)
#include "./gravity/Default/gravity_io.h"
#elif defined(MULTI_SOFTENING_GRAVITY)
#include "./gravity/MultiSoftening/gravity_io.h"
#else
#error "Invalid choice of gravity variant"
#endif

#endif /* SWIFT_GRAVITY_IO_H */
