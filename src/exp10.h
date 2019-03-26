/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EXP10_H
#define SWIFT_EXP10_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "inline.h"

#if !defined(HAVE_EXP10) && !defined(HAVE___EXP10)

/**
 * @brief Raises 10 to the power of the argument.
 *
 * This function is only used as a replacement for compilers that do
 * not implement GNU extensions to the C language.
 *
 * @param x The input value.
 */
__attribute__((always_inline, const)) INLINE static double exp10(
    const double x) {

  return exp(x * M_LN10);
}

#endif

#if !defined(HAVE_EXP10F) && !defined(HAVE___EXP10F)

/**
 * @brief Raises 10 to the power of the argument.
 *
 * This function is only used as a replacement for compilers that do
 * not implement GNU extensions to the C language.
 *
 * @param x The input value.
 */
__attribute__((always_inline, const)) INLINE static float exp10f(
    const float x) {

  return expf(x * (float)M_LN10);
}

#endif

/* Use the __exp10 and __exp10f versions if needed. */
#if !defined(HAVE_EXP10) && defined(HAVE___EXP10)
#define exp10(x) __exp10(x)
#endif

#if !defined(HAVE_EXP10F) && defined(HAVE___EXP10F)
#define exp10f(x) __exp10f(x)
#endif

#endif /* SWIFT_EXP10_H */
