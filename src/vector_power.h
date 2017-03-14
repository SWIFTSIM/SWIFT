/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_VECTOR_POWER_H
#define SWIFT_VECTOR_POWER_H

/**
 * @file vector_power.h
 * @brief Powers of 3D vectors to a multi-index with factorial.
 *
 * These expressions are to be used in 3D Taylor series.
 *
 * We use the notation of Dehnen, Computational Astrophysics and Cosmology,
 * 1, 1, pp. 24 (2014), arXiv:1405.2255.
 *
 * We compute \f$ \frac{1}{\vec{m}!}\vec{v}^{\vec{m}} \f$ for all relevant m.
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/***************************/
/* 0th order vector powers */
/***************************/

/**
 * @brief \f$ \frac{1}{(0,0,0)!}\vec{v}^{(0,0,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_000(const double v[3]) {

  return 1.;
}

/***************************/
/* 1st order vector powers */
/***************************/

/**
 * @brief \f$ \frac{1}{(1,0,0)!}\vec{v}^{(1,0,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_100(const double v[3]) {

  return v[0];
}

/**
 * @brief \f$ \frac{1}{(0,1,0)!}\vec{v}^{(0,1,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_010(const double v[3]) {

  return v[1];
}

/**
 * @brief \f$ \frac{1}{(0,0,1)!}\vec{v}^{(0,0,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_001(const double v[3]) {

  return v[2];
}

/***************************/
/* 2nd order vector powers */
/***************************/

/**
 * @brief \f$ \frac{1}{(2,0,0)!}\vec{v}^{(2,0,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_200(const double v[3]) {

  return 0.5 * v[0] * v[0];
}

/**
 * @brief \f$ \frac{1}{(0,2,0)!}\vec{v}^{(0,2,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_020(const double v[3]) {

  return 0.5 * v[1] * v[1];
}

/**
 * @brief \f$ \frac{1}{(0,0,2)!}\vec{v}^{(0,0,2)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_002(const double v[3]) {

  return 0.5 * v[2] * v[2];
}

/**
 * @brief \f$ \frac{1}{(1,1,0)!}\vec{v}^{(1,1,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_110(const double v[3]) {

  return v[0] * v[1];
}

/**
 * @brief \f$ \frac{1}{(1,0,1)!}\vec{v}^{(1,0,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_101(const double v[3]) {

  return v[0] * v[2];
}

/**
 * @brief \f$ \frac{1}{(0,1,1)!}\vec{v}^{(0,1,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_011(const double v[3]) {

  return v[1] * v[2];
}

/***************************/
/* 3rd order vector powers */
/***************************/

/**
 * @brief \f$ \frac{1}{(3,0,0)!}\vec{v}^{(3,0,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_300(const double v[3]) {

  return 0.1666666666666667 * v[0] * v[0] * v[0];
}

/**
 * @brief \f$ \frac{1}{(0,3,0)!}\vec{v}^{(0,3,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_030(const double v[3]) {

  return 0.1666666666666667 * v[1] * v[1] * v[1];
}

/**
 * @brief \f$ \frac{1}{(0,0,3)!}\vec{v}^{(0,0,3)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_003(const double v[3]) {

  return 0.1666666666666667 * v[2] * v[2] * v[2];
}

/**
 * @brief \f$ \frac{1}{(2,1,0)!}\vec{v}^{(2,1,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_210(const double v[3]) {

  return 0.5 * v[0] * v[0] * v[1];
}

/**
 * @brief \f$ \frac{1}{(2,0,1)!}\vec{v}^{(2,0,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_201(const double v[3]) {

  return 0.5 * v[0] * v[0] * v[2];
}

/**
 * @brief \f$ \frac{1}{(1,2,0)!}\vec{v}^{(1,2,0)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_120(const double v[3]) {

  return 0.5 * v[0] * v[1] * v[1];
}

/**
 * @brief \f$ \frac{1}{(0,2,1)!}\vec{v}^{(0,2,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_021(const double v[3]) {

  return 0.5 * v[1] * v[1] * v[2];
}

/**
 * @brief \f$ \frac{1}{(1,0,2)!}\vec{v}^{(1,0,2)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_102(const double v[3]) {

  return 0.5 * v[0] * v[2] * v[2];
}

/**
 * @brief \f$ \frac{1}{(0,1,2)!}\vec{v}^{(0,1,2)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_012(const double v[3]) {

  return 0.5 * v[1] * v[2] * v[2];
}

/**
 * @brief \f$ \frac{1}{(1,1,1)!}\vec{v}^{(1,1,1)} \f$.
 *
 * @param v vector (\f$ v \f$).
 */
__attribute__((always_inline)) INLINE static double X_111(const double v[3]) {

  return v[0] * v[1] * v[2];
}

#endif /* SWIFT_VECTOR_POWER_H */
