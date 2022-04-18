/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_SIGN_H
#define SWIFT_SIGN_H

/**
 * @brief Return the sign of a floating point number
 *
 * @param x The number of interest.
 * @return >0 if positive, 0 if negative.
 */
__attribute__((always_inline, const)) INLINE static int signf(float x) {
#ifdef __GNUC__
  return !signbit(x);
#else
  return (0.f < val) - (val < 0.f);
#endif
}

/**
 * @brief Return 1 if two numbers have the same sign, 0 otherwise
 *
 * @param x The first number
 * @param y The second number
 */
__attribute__((always_inline, const)) INLINE static int same_signf(float x,
                                                                   float y) {
  return signf(x) == signf(y);
}

#endif
