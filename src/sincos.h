/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SINCOS_H
#define SWIFT_SINCOS_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "inline.h"

#if !defined(HAVE_SINCOS) && !defined(HAVE___SINCOS)

/**
 * @brief Compute both the sin() and cos() of a number.
 *
 * This function is only used as a replacement for compilers that do
 * not implement GNU extensions to the C language.
 *
 * @param x The input value.
 * @param s (return) The sine of x.
 * @param c (return) The cosine of x.
 */
__attribute__((always_inline)) INLINE static void sincos(const double x,
                                                         double *restrict s,
                                                         double *restrict c) {

  *s = sin(x);
  *c = cos(x);
}

#endif

#if !defined(HAVE_SINCOSF) && !defined(HAVE___SINCOSF)

/**
 * @brief Compute both the sin() and cos() of a number.
 *
 * This function is only used as a replacement for compilers that do
 * not implement GNU extensions to the C language.
 *
 * @param x The input value.
 * @param s (return) The sine of x.
 * @param c (return) The cosine of x.
 */
__attribute__((always_inline)) INLINE static void sincosf(const float x,
                                                          float *restrict s,
                                                          float *restrict c) {

  *s = sinf(x);
  *c = cosf(x);
}

#endif

/* Use the __sincos and __sincosf versions if needed. */
#if !defined(HAVE_SINCOS) && defined(HAVE___SINCOS)
#define sincos(x, s, c) __sincos(x, s, c)
#endif

#if !defined(HAVE_SINCOSF) && defined(HAVE___SINCOSF)
#define sincosf(x, s, c) __sincosf(x, s, c)
#endif

#endif /* SWIFT_SINCOS_H */
