/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_INLINE_H
#define SWIFT_INLINE_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief Defines inline
 */
#ifndef INLINE
#ifdef __cplusplus
#define INLINE inline
#else
#if __GNUC__ && !__GNUC_STDC_INLINE__
#define INLINE extern inline
#else
#define INLINE inline
#endif
#endif
#endif

#endif /* SWIFT_INLINE_H */
