/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_DIMENSION_H
#define SWIFT_DIMENSION_H

/**
 * @file dimension.h
 * @brief Defines the dimensionality \f$d\f$ of the problem and (fast)
 * mathematical functions involving it
 */

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "inline.h"
#include "vector.h"

/* First define some constants */
#if defined(HYDRO_DIMENSION_3D)

#define hydro_dimension 3.f
#define hydro_dimension_inv 0.3333333333f
#define hydro_dimention_unit_sphere ((float)(4. * M_PI / 3.))

#elif defined(HYDRO_DIMENSION_2D)

#define hydro_dimension 2.f
#define hydro_dimension_inv 0.5f
#define hydro_dimention_unit_sphere ((float)M_PI)

#elif defined(HYDRO_DIMENSION_1D)

#define hydro_dimension 1.f
#define hydro_dimension_inv 1.f
#define hydro_dimention_unit_sphere 2.f

#else

#error "A problem dimensionality must be chosen in const.h !"

#endif

/**
 * @brief Returns the argument to the power given by the dimension
 *
 * Computes \f$x^d\f$.
 */
__attribute__((always_inline)) INLINE static float pow_dimension(float x) {

#if defined(HYDRO_DIMENSION_3D)

  return x * x * x;

#elif defined(HYDRO_DIMENSION_2D)

  return x * x;

#elif defined(HYDRO_DIMENSION_1D)

  return x;

#else

  error("The dimension is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns the argument to the power given by the dimension plus one
 *
 * Computes \f$x^{d+1}\f$.
 */
__attribute__((always_inline)) INLINE static float pow_dimension_plus_one(
    float x) {

#if defined(HYDRO_DIMENSION_3D)

  const float x2 = x * x;
  return x2 * x2;

#elif defined(HYDRO_DIMENSION_2D)

  return x * x * x;

#elif defined(HYDRO_DIMENSION_1D)

  return x * x;

#else

  error("The dimension is not defined !");
  return 0.f;

#endif
}

/* ------------------------------------------------------------------------- */
#ifdef WITH_VECTORIZATION

/**
 * @brief Returns the argument to the power given by the dimension (vector
 * version)
 *
 * Computes \f$x^d\f$.
 */
__attribute__((always_inline)) INLINE static vector pow_dimension_vec(
    vector x) {

#if defined(HYDRO_DIMENSION_3D)

  return (vector)(x.v * x.v * x.v);

#elif defined(HYDRO_DIMENSION_2D)

  return (vector)(x.v * x.v);

#elif defined(HYDRO_DIMENSION_1D)

  return x;

#else

  error("The dimension is not defined !");
  return vec_set(0.f);

#endif
}

/**
 * @brief Returns the argument to the power given by the dimension plus one
 * (vector version)
 *
 * Computes \f$x^{d+1}\f$.
 */
__attribute__((always_inline)) INLINE static vector pow_dimension_plus_one_vec(
    vector x) {

#if defined(HYDRO_DIMENSION_3D)

  const vector x2 = (vector) (x.v * x.v);
  return (vector)(x2.v * x2.v);

#elif defined(HYDRO_DIMENSION_2D)

  return (vector)(x.v * x.v * x.v);

#elif defined(HYDRO_DIMENSION_1D)

  return (vector)(x.v * x.v);

#else

  error("The dimension is not defined !");
  return vec_set(0.f);

#endif
}
#endif

#endif /* SWIFT_DIMENSION_H */
