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
#include "inline.h"
#include "vector.h"

#include <math.h>

/* First define some constants */
#if defined(HYDRO_DIMENSION_3D)

#define hydro_dimension 3.f
#define hydro_dimension_inv 0.3333333333f
#define hydro_dimension_unit_sphere ((float)(4. * M_PI / 3.))
#define hydro_dimension_unit_sphere_inv ((float)(3. * M_1_PI / 4.))

#elif defined(HYDRO_DIMENSION_2D)

#define hydro_dimension 2.f
#define hydro_dimension_inv 0.5f
#define hydro_dimension_unit_sphere ((float)M_PI)
#define hydro_dimension_unit_sphere_inv ((float)M_1_PI)

#elif defined(HYDRO_DIMENSION_1D)

#define hydro_dimension 1.f
#define hydro_dimension_inv 1.f
#define hydro_dimension_unit_sphere 2.f
#define hydro_dimension_unit_sphere_inv 0.5f

#else

#error "A problem dimensionality must be chosen in config.h !"

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
 * @brief Returns the argument to the power given by the inverse of the
 * dimension
 *
 * Computes \f$x^{1/d}\f$.
 */
__attribute__((always_inline)) INLINE static float pow_inv_dimension(float x) {

#if defined(HYDRO_DIMENSION_3D)

  return cbrtf(x);

#elif defined(HYDRO_DIMENSION_2D)

  return sqrtf(x);

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

/**
 * @brief Returns the argument to the power given by the dimension minus one
 *
 * Computes \f$x^{d-1}\f$.
 */
__attribute__((always_inline)) INLINE static float pow_dimension_minus_one(
    float x) {

#if defined(HYDRO_DIMENSION_3D)

  return x * x;

#elif defined(HYDRO_DIMENSION_2D)

  return x;

#elif defined(HYDRO_DIMENSION_1D)

  return 1.f;

#else

  error("The dimension is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Inverts the given dimension by dimension matrix (in place)
 *
 * @param A A 3x3 matrix of which we want to invert the top left dxd part
 */
__attribute__((always_inline)) INLINE static void
invert_dimension_by_dimension_matrix(float A[3][3]) {

#if defined(HYDRO_DIMENSION_3D)

  float detA, Ainv[3][3];

  detA = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
         A[0][2] * A[1][0] * A[2][1] - A[0][2] * A[1][1] * A[2][0] -
         A[0][1] * A[1][0] * A[2][2] - A[0][0] * A[1][2] * A[2][1];

  if (detA && !isnan(detA)) {
    Ainv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / detA;
    Ainv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / detA;
    Ainv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / detA;
    Ainv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / detA;
    Ainv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / detA;
    Ainv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / detA;
    Ainv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detA;
    Ainv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / detA;
    Ainv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / detA;
  } else {
    Ainv[0][0] = 0.0f;
    Ainv[0][1] = 0.0f;
    Ainv[0][2] = 0.0f;
    Ainv[1][0] = 0.0f;
    Ainv[1][1] = 0.0f;
    Ainv[1][2] = 0.0f;
    Ainv[2][0] = 0.0f;
    Ainv[2][1] = 0.0f;
    Ainv[2][2] = 0.0f;
  }

  A[0][0] = Ainv[0][0];
  A[0][1] = Ainv[0][1];
  A[0][2] = Ainv[0][2];
  A[1][0] = Ainv[1][0];
  A[1][1] = Ainv[1][1];
  A[1][2] = Ainv[1][2];
  A[2][0] = Ainv[2][0];
  A[2][1] = Ainv[2][1];
  A[2][2] = Ainv[2][2];

#elif defined(HYDRO_DIMENSION_2D)

  float detA, Ainv[2][2];

  detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

  if (detA && !isnan(detA)) {
    Ainv[0][0] = A[1][1] / detA;
    Ainv[0][1] = -A[0][1] / detA;
    Ainv[1][0] = -A[1][0] / detA;
    Ainv[1][1] = A[0][0] / detA;
  } else {
    Ainv[0][0] = 0.0f;
    Ainv[0][1] = 0.0f;
    Ainv[1][0] = 0.0f;
    Ainv[1][1] = 0.0f;
  }

  A[0][0] = Ainv[0][0];
  A[0][1] = Ainv[0][1];
  A[1][0] = Ainv[1][0];
  A[1][1] = Ainv[1][1];

#elif defined(HYDRO_DIMENSION_1D)

  if (A[0][0] && !isnan(A[0][0])) {
    A[0][0] = 1.0f / A[0][0];
  } else {
    A[0][0] = 0.0f;
  }

#else

  error("The dimension is not defined !");

#endif
}

/**
 * @brief Get the radius of a dimension sphere with the given volume
 *
 * @param volume Volume of the dimension sphere
 * @return Radius of the dimension sphere
 */
__attribute__((always_inline)) INLINE static float get_radius_dimension_sphere(
    float volume) {

#if defined(HYDRO_DIMENSION_3D)

  return cbrtf(volume * hydro_dimension_unit_sphere_inv);

#elif defined(HYDRO_DIMENSION_2D)

  return sqrtf(volume * hydro_dimension_unit_sphere_inv);

#elif defined(HYDRO_DIMENSION_1D)

  return volume * hydro_dimension_unit_sphere_inv;

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

  return (vector)(vec_mul(vec_mul(x.v, x.v), x.v));

#elif defined(HYDRO_DIMENSION_2D)

  return (vector)(vec_mul(x.v, x.v));

#elif defined(HYDRO_DIMENSION_1D)

  return x;

#else

  error("The dimension is not defined !");
  return vec_set1(0.f);

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

  const vector x2 = (vector)(vec_mul(x.v, x.v));
  return (vector)(vec_mul(x2.v, x2.v));

#elif defined(HYDRO_DIMENSION_2D)

  return (vector)(vec_mul(x.v, vec_mul(x.v, x.v)));

#elif defined(HYDRO_DIMENSION_1D)

  return (vector)(vec_mul(x.v, x.v));

#else

  error("The dimension is not defined !");
  return vec_set1(0.f);

#endif
}
#endif

#endif /* SWIFT_DIMENSION_H */
