/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include <config.h>

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
#define hydro_dimension_integer 3

#elif defined(HYDRO_DIMENSION_2D)

#define hydro_dimension 2.f
#define hydro_dimension_inv 0.5f
#define hydro_dimension_unit_sphere ((float)M_PI)
#define hydro_dimension_unit_sphere_inv ((float)M_1_PI)
#define hydro_dimension_integer 2

#elif defined(HYDRO_DIMENSION_1D)

#define hydro_dimension 1.f
#define hydro_dimension_inv 1.f
#define hydro_dimension_unit_sphere 2.f
#define hydro_dimension_unit_sphere_inv 0.5f
#define hydro_dimension_integer 1

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
 * @param A 3x3 matrix of which we want to invert the top left dxd part
 * @param min_cond_num Minimal condition number to attempt an inversion. Smaller
 * values will trigger the singular matrix case.
 * @return Exit code: 0 for success, 1 if a singular matrix was detected.
 */
__attribute__((always_inline)) INLINE static int
invert_dimension_by_dimension_matrix(
    float A[hydro_dimension_integer][hydro_dimension_integer],
    const float min_cond_num) {

#if defined(HYDRO_DIMENSION_3D)

  /* Calculate root mean square of elements. */
  float rms = 0.0f;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rms += A[i][j] * A[i][j];
    }
  }
  rms = sqrtf(rms / 9.0f);

  /* Early abort if the matrix is all zeros. */
  if (rms == 0.0f) {
    return 1;
  }

  /* Scale to avoid any issues if values are tiny. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A[i][j] /= rms;
    }
  }

  int pivot[3];
  for (int i = 0; i < 3; i++) {
    int imax = i;
    float Smax = fabsf(A[imax][i]);
    for (int j = i + 1; j < 3; j++) {
      const float this_Smax = fabsf(A[j][i]);
      if (this_Smax > Smax) {
        Smax = this_Smax;
        imax = j;
      }
    }

    if (Smax < min_cond_num) {
      /* singular matrix. Early abort */
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          A[j][k] = 0.0f;
        }
      }
      return 1;
    }

    pivot[i] = imax;
    if (i != imax) {
      for (int j = 0; j < 3; j++) {
        const float temp = A[i][j];
        A[i][j] = A[imax][j];
        A[imax][j] = temp;
      }
    }

    const float Aii_inv = 1.0f / A[i][i];
    for (int j = i + 1; j < 3; j++) {
      A[j][i] *= Aii_inv;
    }

    for (int j = i + 1; j < 3; j++) {
      for (int k = i + 1; k < 3; k++) {
        A[j][k] -= A[j][i] * A[i][k];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    A[i][i] = 1.0f / A[i][i];
    for (int j = i + 1; j < 3; j++) {
      float Aij = 0.0f;
      for (int k = i; k < j; k++) {
        Aij -= A[i][k] * A[k][j];
      }
      A[i][j] = Aij / A[j][j];
    }
  }

  float work[3];
  for (int jp1 = 3; jp1 > 0; jp1--) {
    const int j = jp1 - 1;
    for (int i = 0; i < jp1; i++) {
      work[i] = A[i][j];
    }
    for (int i = jp1; i < 3; i++) {
      work[i] = 0.0f;
    }
    for (int k = jp1; k < 3; k++) {
      for (int i = 0; i < 3; i++) {
        work[i] -= A[i][k] * A[k][j];
      }
    }
    for (int i = 0; i < 3; i++) {
      A[i][j] = work[i];
    }
  }

  for (int jp1 = 3; jp1 > 0; jp1--) {
    const int j = jp1 - 1;
    const int jp = pivot[j];
    if (jp != j) {
      for (int i = 0; i < 3; i++) {
        const float temp = A[i][j];
        A[i][j] = A[i][jp];
        A[i][jp] = temp;
      }
    }
  }

  /* Scale the inverted matrix back. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A[i][j] /= rms;
    }
  }

  return 0;

#elif defined(HYDRO_DIMENSION_2D)

  /* Calculate root mean square of elements. */
  float rms = 0.0f;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      rms += A[i][j] * A[i][j];
    }
  }
  rms = sqrtf(rms / 4.0f);

  /* Early abort if the matrix is all zeros. */
  if (rms == 0.0f) {
    return 1;
  }

  /* Scale to avoid any issues if values are tiny. */
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      A[i][j] /= rms;
    }
  }

  float Ainv[2][2];

  const float detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

  if (fabsf(detA) < min_cond_num) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        A[j][k] = 0.0f;
      }
    }
    return 1;
  }

  Ainv[0][0] = A[1][1] / detA;
  Ainv[0][1] = -A[0][1] / detA;
  Ainv[1][0] = -A[1][0] / detA;
  Ainv[1][1] = A[0][0] / detA;

  A[0][0] = Ainv[0][0];
  A[0][1] = Ainv[0][1];
  A[1][0] = Ainv[1][0];
  A[1][1] = Ainv[1][1];

  /* Scale the inverted matrix back. */
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      A[i][j] /= rms;
    }
  }

  return 0;

#elif defined(HYDRO_DIMENSION_1D)

  if (A[0][0] == 0.0f || isnan(A[0][0]) || isinf(A[0][0])) {
    A[0][0] = 0.0f;
    return 1;
  }
  A[0][0] = 1.0f / A[0][0];

  return 0;

#else

  error("The dimension is not defined !");

#endif
}

/**
 * @brief Computes the 2-norm condition number of a 3x3
 * row-major matrix.
 *
 * Textbook implementation matchin the result of a GSL call to
 * gsl_linalg_SV_decomp() and taking the max/min ratio of the sigma values.
 *
 * @param m The matrix.
 */
__attribute__((always_inline)) INLINE static double
matrix_3x3_2norm_condition_number(const double m[3][3]) {

  /* Form the symmetric matrix S = m^T * m */
  const double s0 = m[0][0] * m[0][0] + m[1][0] * m[1][0] + m[2][0] * m[2][0];
  const double s1 = m[0][0] * m[0][1] + m[1][0] * m[1][1] + m[2][0] * m[2][1];
  const double s2 = m[0][0] * m[0][2] + m[1][0] * m[1][2] + m[2][0] * m[2][2];
  const double s4 = m[0][1] * m[0][1] + m[1][1] * m[1][1] + m[2][1] * m[2][1];
  const double s5 = m[0][1] * m[0][2] + m[1][1] * m[1][2] + m[2][1] * m[2][2];
  const double s8 = m[0][2] * m[0][2] + m[1][2] * m[1][2] + m[2][2] * m[2][2];

  /* Compute invariants of S (coefficients of characteristic polynomial) */
  const double c2 = s0 + s4 + s8;
  const double c1 =
      (s0 * s4 - s1 * s1) + (s0 * s8 - s2 * s2) + (s4 * s8 - s5 * s5);
  const double c0 = s0 * (s4 * s8 - s5 * s5) - s1 * (s1 * s8 - s2 * s5) +
                    s2 * (s1 * s5 - s2 * s4);

  /* Solve the cubic equation analytically using Cardano's formulation */
  const double p = c1 - (c2 * c2) / 3.0;
  const double q = c0 - (c2 * c1) / 3.0 + (2.0 * c2 * c2 * c2) / 27.0;

  double ev_max, ev_min;

  if (p >= 0.0) {
    ev_max = c2 / 3.0;
    ev_min = c2 / 3.0;
  } else {
    const double r = sqrt(-p / 3.0);
    double val = q / (2.0 * r * r * r);

    /* Clamp to prevent out-of-bounds inputs to acos due to numerical drift */
    if (val > 1.0) val = 1.0;
    if (val < -1.0) val = -1.0;

    const double phi = acos(val);

    /* Find the three roots (eigenvalues of m^T*m) */
    const double r1 = c2 / 3.0 + 2.0 * r * cos(phi / 3.0);
    const double r2 = c2 / 3.0 + 2.0 * r * cos((phi + 2.0 * M_PI) / 3.0);
    const double r3 = c2 / 3.0 + 2.0 * r * cos((phi + 4.0 * M_PI) / 3.0);

    /* Sort to isolate max and min eigenvalues */
    ev_max = r1;
    if (r2 > ev_max) ev_max = r2;
    if (r3 > ev_max) ev_max = r3;
    ev_min = r1;
    if (r2 < ev_min) ev_min = r2;
    if (r3 < ev_min) ev_min = r3;
  }

  /* Return condition number (sigma_max / sigma_min) */
  if (ev_min <= 1e-15) return INFINITY;

  return sqrt(ev_max / ev_min);
}

/**
 * @brief Compute the inverse of 3x3 matrix using LU decomposition.
 *
 * Textbook implementation matching a call to gsl_linalg_LU_decomp()
 * and gsl_linalg_LU_invert().
 *
 * @param A the matrix to invert.
 * @param inv (return) The inverse of the matrix.
 * @return 1 if the inversion failed, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int invert3x3_matrix_LU(
    const double A[3][3], double inv[3][3]) {
  double mat[3][3];
  double scale_factors[3];

  /* Initialize identity matrix and compute row scale factors */
  for (int i = 0; i < 3; ++i) {
    double max_val = 0.0;
    for (int j = 0; j < 3; ++j) {
      mat[i][j] = A[i][j];
      inv[i][j] = (i == j) ? 1.0 : 0.0;

      double abs_val = fabs(A[i][j]);
      if (abs_val > max_val) max_val = abs_val;
    }

    /* If an entire row is 0, the matrix is singular */
    if (max_val < 1e-15) return 1;
    scale_factors[i] = 1.0 / max_val;
  }

  /* Gaussian elimination with scaled partial pivoting */
  for (int i = 0; i < 3; ++i) {

    /* Find pivot row using scaled values to eliminate magnitude bias */
    int pivot = i;
    double max_scaled_val = fabs(mat[i][i]) * scale_factors[i];

    for (int r = i + 1; r < 3; ++r) {
      const double scaled_val = fabs(mat[r][i]) * scale_factors[r];
      if (scaled_val > max_scaled_val) {
        max_scaled_val = scaled_val;
        pivot = r;
      }
    }

    /* Singularity check based on machine epsilon threshold */
    if (max_scaled_val < 1e-15) return 1;

    /* Swap rows if a better pivot was found */
    if (pivot != i) {

      /* Swap scale factors */
      const double ts = scale_factors[i];
      scale_factors[i] = scale_factors[pivot];
      scale_factors[pivot] = ts;

      /* Swap working matrix rows */
      for (int j = 0; j < 3; ++j) {
        double t = mat[i][j];
        mat[i][j] = mat[pivot][j];
        mat[pivot][j] = t;
        t = inv[i][j];
        inv[i][j] = inv[pivot][j];
        inv[pivot][j] = t;
      }
    }

    /* Eliminate column elements below the pivot */
    for (int r = i + 1; r < 3; ++r) {
      double factor = mat[r][i] / mat[i][i];
      for (int j = 0; j < 3; ++j) {
        mat[r][j] -= factor * mat[i][j];
        inv[r][j] -= factor * inv[i][j];
      }
    }
  }

  /* Back-substitution to finalize the inverse */
  for (int i = 2; i >= 0; --i) {
    for (int r = i - 1; r >= 0; --r) {
      double factor = mat[r][i] / mat[i][i];
      for (int j = 0; j < 3; ++j) {
        inv[r][j] -= factor * inv[i][j];
      }
    }
    double d = 1.0 / mat[i][i];
    for (int j = 0; j < 3; ++j) inv[i][j] *= d;
  }

  return 0;
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
