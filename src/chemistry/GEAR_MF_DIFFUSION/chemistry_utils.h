/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@ealumni.pfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UTILS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UTILS_H

#include "inline.h"
#include "minmax.h"

#include <math.h>

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_utils.h
 * @brief File containing functions to diagonalize real 3x3 symmetric
 * matrices. Functions in this file were adapted from the C++ implementation of
 * David Eberly (Geometric Tools) non iterative eigen-solver
 * (https://www.geometrictools.com/GTE/Mathematics/SymmetricEigensolver3x3.h).
 * Mathematical and numerical details are provided in
 * https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf.
 * */

/**
 * @brief Computes the cross product of two 3D vectors.
 *
 * @param u A 3D vector (array of 3 doubles).
 * @param v A 3D vector (array of 3 doubles).
 * @param result (return) The resulting 3D vector (array of 3 doubles) from the
 * cross product of u and v.
 */
__attribute__((always_inline)) INLINE static void chemistry_utils_cross_product(
    const double u[3], const double v[3], double result[3]) {
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}

/**
 * @brief Computes the dot product of two 3D vectors.
 *
 * @param u A 3D vector (array of 3 doubles).
 * @param v A 3D vector (array of 3 doubles).
 * @return The dot product (scalar) of vectors u and v.
 */
__attribute__((always_inline)) INLINE static double chemistry_utils_dot_product(
    const double u[3], const double v[3]) {
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

/**
 * @brief Normalizes a 3D vector to unit length.
 *
 * @param v (return) A 3D vector (array of 3 doubles) to be normalized. The
 * vector is modified in place.
 *
 * @note If the vector's magnitude is zero, the behavior is undefined.
 */
__attribute__((always_inline)) INLINE static void chemistry_utils_normalize(
    double v[3]) {
  const double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  for (int i = 0; i < 3; i++) {
    v[i] /= norm;
  }
}

/**
 * @brief Computes an orthogonal complement to a given 3D unit vector.
 *
 * @param W A 3D unit vector (array of 3 doubles). Assumes W is normalized.
 * @param U (return) A 3D vector (array of 3 doubles) to store one orthogonal
 * complement.
 * @param V (return) A 3D vector (array of 3 doubles) to store the second
 * orthogonal complement.
 *
 * @details Generates a robust right-handed orthonormal basis { U, V, W },
 *          where U and V are orthogonal to W and to each other. The vector W
 *          is assumed to have unit length, ensuring stability in the
 * computations.
 *
 * @note The orthonormal basis is computed such that:
 *       - `chemistry_utils_dot_product(U, W) == 0`
 *       - `chemistry_utils_dot_product(V, W) == 0`
 *       - `chemistry_utils_dot_product(U, V) == 0`
 */
__attribute__((always_inline)) INLINE static void
chemistry_utils_compute_orthogonal_complement(const double W[3], double U[3],
                                              double V[3]) {
  /* Robustly compute a right-handed orthonormal set { U, V, W }. The vector W
     is guaranteed to be unit-length, in which case there is no need to worry
     about a division by zero when computing invLength. */
  double invLength = 0.0;
  if (fabs(W[0]) > fabs(W[1])) {
    /* The component of maximum absolute value is either W[0] or W[2]. */
    invLength = 1.0 / sqrt(W[0] * W[0] + W[2] * W[2]);
    U[0] = -W[2] * invLength;
    U[1] = 0.0;
    U[2] = +W[0] * invLength;
  } else {
    /* The component of maximum absolute value is either W[1] or W[2]. */
    invLength = 1.0 / sqrt(W[1] * W[1] + W[2] * W[2]);
    U[0] = 0.0;
    U[1] = +W[2] * invLength;
    U[2] = -W[1] * invLength;
  }
  chemistry_utils_cross_product(W, U, V);
}

/**
 * @brief Computes the eigenvector corresponding to the smallest eigenvalue of
 * a symmetric 3x3 matrix.
 *
 * This function calculates a unit-length eigenvector corresponding to the
 * given eigenvalue for a symmetric 3x3 matrix. It uses a robust method to
 * ensure numerical stability by selecting two rows of the matrix whose cross
 * product has the largest magnitude, ensuring linear independence.
 *
 * @param a00 Diagonal element of the matrix (row 0, column 0).
 * @param a01 Off-diagonal element of the matrix (row 0, column 1).
 * @param a02 Off-diagonal element of the matrix (row 0, column 2).
 * @param a11 Diagonal element of the matrix (row 1, column 1).
 * @param a12 Off-diagonal element of the matrix (row 1, column 2).
 * @param a22 Diagonal element of the matrix (row 2, column 2).
 * @param eigenvalue0 Smallest eigenvalue of the matrix.
 * @param (return) eigenvector0 Eigenvector corresponding to the smallest
 * eigenvalue.
 */
__attribute__((always_inline)) INLINE static void compute_eigenvector_0(
    const double a00, const double a01, const double a02, const double a11,
    const double a12, const double a22, const double eigenvalue0,
    double eigenvector0[3]) {
  /* Compute a unit-length eigenvector for eigenvalue[i0]. The  matrix is rank
     2, so two of the rows are linearly independent. For a robust computation
     of the eigenvector, select the two rows whose cross product has largest
     length of all pairs of rows. */
  const double row0[3] = {a00 - eigenvalue0, a01, a02};
  const double row1[3] = {a01, a11 - eigenvalue0, a12};
  const double row2[3] = {a02, a12, a22 - eigenvalue0};

  double r0xr1[3], r0xr2[3], r1xr2[3];
  chemistry_utils_cross_product(row0, row1, r0xr1);
  chemistry_utils_cross_product(row0, row2, r0xr2);
  chemistry_utils_cross_product(row1, row2, r1xr2);

  const double d0 = chemistry_utils_dot_product(r0xr1, r0xr1);
  const double d1 = chemistry_utils_dot_product(r0xr2, r0xr2);
  const double d2 = chemistry_utils_dot_product(r1xr2, r1xr2);

  double dmax = d0;
  int imax = 0;
  if (d1 > dmax) {
    dmax = d1;
    imax = 1;
  }
  if (d2 > dmax) {
    imax = 2;
  }

  if (imax == 0) {
    eigenvector0[0] = r0xr1[0] / sqrt(d0);
    eigenvector0[1] = r0xr1[1] / sqrt(d0);
    eigenvector0[2] = r0xr1[2] / sqrt(d0);
  } else if (imax == 1) {
    eigenvector0[0] = r0xr2[0] / sqrt(d1);
    eigenvector0[1] = r0xr2[1] / sqrt(d1);
    eigenvector0[2] = r0xr2[2] / sqrt(d1);
  } else {
    eigenvector0[0] = r1xr2[0] / sqrt(d2);
    eigenvector0[1] = r1xr2[1] / sqrt(d2);
    eigenvector0[2] = r1xr2[2] / sqrt(d2);
  }
}

/**
 * @brief Computes the second eigenvector for a symmetric 3x3 matrix.
 *
 * This function calculates the second eigenvector of a symmetric 3x3 matrix
 * using a right-handed orthonormal basis and a robust numerical method to solve
 * the reduced 2x2 linear system.
 *
 * @param a00 Diagonal element of the matrix (row 0, column 0).
 * @param a01 Off-diagonal element of the matrix (row 0, column 1).
 * @param a02 Off-diagonal element of the matrix (row 0, column 2).
 * @param a11 Diagonal element of the matrix (row 1, column 1).
 * @param a12 Off-diagonal element of the matrix (row 1, column 2).
 * @param a22 Diagonal element of the matrix (row 2, column 2).
 * @param eigenvector0 First eigenvector of the matrix.
 * @param eigenvalue1 Second eigenvalue of the matrix.
 * @param (return) eigenvector1 Second eigenvector of the matrix.
 */
__attribute__((always_inline)) INLINE static void compute_eigenvector_1(
    const double a00, const double a01, const double a02, const double a11,
    const double a12, const double a22, const double eigenvector0[3],
    const double eigenvalue1, double eigenvector1[3]) {

  /* Robustly compute a right-handed orthonormal set { U, V, eigenvector0 }. */
  double U[3], V[3];
  chemistry_utils_compute_orthogonal_complement(eigenvector0, U, V);

  /* Let e be eval1 and let E be a corresponding eigenvector which is a
   solution to the linear system (A - e*I)*E = 0.  The matrix (A - e*I) is 3x3,
   not invertible (so infinitely many solutions), and has rank 2 when eval1 and
   eval are different. It has rank 1 when eval1 and eval2 are equal.
   Numerically, it is difficult to compute robustly the rank of a
   matrix. Instead,the 3x3 linear system is reduced to a 2x2 system as follows.
   Define the 3x2 matrix J = [U V] whose columns are the U and V computed
   previously. Define the 2x1 vector X = J*E.  The 2x2 system is 0 = M * X =
   (J^T * (A - e*I) * J) * X where J^T is the transpose of J and M = J^T * (A -
   e*I) * J is a 2x2 matrix. The system may be written as
       +-                        -++-  -+       +-  -+
       | U^T*A*U - e  U^T*A*V     || x0 | = e * | x0 |
       | V^T*A*U      V^T*A*V - e || x1 |       | x1 |
       +-                        -++   -+       +-  -+
       where X has row entries x0 and x1. */

  const double AU[3] = {a00 * U[0] + a01 * U[1] + a02 * U[2],
                        a01 * U[0] + a11 * U[1] + a12 * U[2],
                        a02 * U[0] + a12 * U[1] + a22 * U[2]};

  const double AV[3] = {a00 * V[0] + a01 * V[1] + a02 * V[2],
                        a01 * V[0] + a11 * V[1] + a12 * V[2],
                        a02 * V[0] + a12 * V[1] + a22 * V[2]};

  double m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eigenvalue1;
  double m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
  double m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eigenvalue1;

  /* For robustness, choose the largest-length row of M to compute the
   eigenvector.  The 2-tuple of coefficients of U and V in the assignments to
   eigenvector[1] lies on a circle, and U and V are unit length and
   perpendicular, so eigenvector[1] is unit length (within numerical
   tolerance). */
  const double absM00 = fabs(m00);
  const double absM01 = fabs(m01);
  const double absM11 = fabs(m11);
  double maxAbsComp;
  if (absM00 >= absM11) {
    maxAbsComp = max(absM00, absM01);
    if (maxAbsComp > 0.0) {
      if (absM00 >= absM01) {
        m01 /= m00;
        m00 = 1.0 / sqrt(1.0 + m01 * m01);
        m01 *= m00;
      } else {
        m00 /= m01;
        m01 = 1.0 / sqrt(1.0 + m00 * m00);
        m00 *= m01;
      }
      /* eigenvector1 = Subtract(Multiply(m01, U), Multiply(m00, V))
         Intermediate results: m01 * U and m00 * V */
      double m01U[3], m00V[3];

      /* Multiply m01 with U and m00 with V */
      for (int i = 0; i < 3; i++) {
        m01U[i] = m01 * U[i];
        m00V[i] = m00 * V[i];
      }
      /* Subtract m00 * V from m01 * U to get eigenvector1 */
      for (int i = 0; i < 3; i++) {
        eigenvector1[i] = m01U[i] - m00V[i];
      }
    } else {
      eigenvector1 = U;
    }
  } else {
    maxAbsComp = max(absM11, absM01);
    if (maxAbsComp > 0.0) {
      if (absM11 >= absM01) {
        m01 /= m11;
        m11 = 1.0 / sqrt(1.0 + m01 * m01);
        m01 *= m11;
      } else {
        m11 /= m01;
        m01 = 1.0 / sqrt(1.0 + m11 * m11);
        m11 *= m01;
      }
      /* eigenvector1 = Subtract(Multiply(m11, U), Multiply(m01, V))
         Intermediate results: m11 * U and m01 * V */
      double m11U[3], m01V[3];

      /* Multiply m11 with U and m01 with V */
      for (int i = 0; i < 3; i++) {
        m11U[i] = m11 * U[i];
        m01V[i] = m01 * V[i];
      }

      /* Subtract m01 * V from m11 * U to get eigenvector1 */
      for (int i = 0; i < 3; i++) {
        eigenvector1[i] = m11U[i] - m01V[i];
      }
    } else {
      eigenvector1 = U;
    }
  }
}

/**
 * @brief Diagonalizes a real symmetric 3x3 matrix.
 *
 * This function computes the eigenvalues and eigenvectors of a real symmetric
 * 3x3 matrix. It preconditions the matrix to prevent numerical instability,
 * calculates eigenvalues using an analytical method, and computes eigenvectors
 * to ensure they form an orthonormal and right-handed set.
 *
 * @param A The symmetric 3x3 matrix to diagonalize (input).
 * @param (return) eigenvalues Array to store the three eigenvalues of the
 * matrix.
 * @param (return) eigenvector0 Array to store the first eigenvector.
 * @param (return) eigenvector1 Array to store the second eigenvector.
 * @param (return) eigenvector2 Array to store the third eigenvector.
 */
__attribute__((always_inline)) INLINE static void
chemistry_utils_diagonalize_3x3(const double A[3][3], double eigenvalues[3],
                                double eigenvector0[3], double eigenvector1[3],
                                double eigenvector2[3]) {
  /* Precondition the matrix by factoring out the maximum absolute value of the
     components. This guards against floating-point overflow when computing
     the eigenvalues. */
  const double max0 = max(fabs(A[0][0]), fabs(A[0][1]));
  const double max1 = max(fabs(A[0][2]), fabs(A[1][1]));
  const double max2 = max(fabs(A[1][2]), fabs(A[2][2]));
  const double max01 = max(max0, max1);
  const double maxAbsElement = max(max01, max2);
  if (maxAbsElement == 0.0) {
    /* A is the zero matrix. */
    eigenvalues[0] = 0.0;
    eigenvalues[1] = 0.0;
    eigenvalues[2] = 0.0;
    eigenvector0[0] = 1.0;
    eigenvector0[1] = 0.0;
    eigenvector0[2] = 0.0;
    eigenvector1[0] = 0.0;
    eigenvector1[1] = 1.0;
    eigenvector1[2] = 0.0;
    eigenvector2[0] = 0.0;
    eigenvector2[1] = 0.0;
    eigenvector2[2] = 1.0;
    return;
  }

  const double invMaxAbsElement = 1.0 / maxAbsElement;
  double A_p[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A_p[i][j] = invMaxAbsElement * A[i][j];
    }
  }

  const double norm =
      A_p[0][1] * A_p[0][1] + A_p[0][2] * A_p[0][2] + A_p[1][2] * A_p[1][2];

  if (norm > 0.0) {
    /* Compute the eigenvalues of A. */
    const double a00 = A_p[0][0], a01 = A_p[0][1], a02 = A_p[0][2];
    const double a11 = A_p[1][1], a12 = A_p[1][2], a22 = A_p[2][2];

    /* Compute the trace q = tr(A)/3 */
    const double q = (a00 + a11 + a22) / 3.0;

    /* Matrix A - q*I, where b00, b11, b22 are the diagonal entries */
    const double b00 = a00 - q;
    const double b11 = a11 - q;
    const double b22 = a22 - q;

    /* Compute p = sqrt(tr((A - q*I)^2)/6) */
    const double p = sqrt((b00 * b00 + b11 * b11 + b22 * b22 + 2 * norm) / 6.0);

    /* Cofactor expansion to compute det(A - q*I) */
    const double c00 = b11 * b22 - a12 * a12;
    const double c01 = a01 * b22 - a12 * a02;
    const double c02 = a01 * a12 - b11 * a02;
    const double det = (b00 * c00 - a01 * c01 + a02 * c02) / (p * p * p);

    /* Compute halfDet and clamp it to [-1, 1] */
    double halfDet = det * 0.5;
    const double tmp = max(halfDet, -1.0);
    halfDet = min(tmp, 1.0);

    /* Compute the angle for the eigenvalues */
    const double angle = acos(halfDet) / 3.0;
    const double twoThirdsPi = 2.09439510239319549;

    /* Compute the eigenvalues of B (beta0, beta1, beta2) */
    const double beta2 = cos(angle) * 2.0;
    const double beta0 = cos(angle + twoThirdsPi) * 2.0;
    const double beta1 = -(beta0 + beta2);

    /* The eigenvalues of A are then calculated as */
    eigenvalues[0] = q + p * beta0;
    eigenvalues[1] = q + p * beta1;
    eigenvalues[2] = q + p * beta2;

    /* Compute the eigenvectors so that the set {eigenvector0, eigenvector1,
       eigenvector2} is right handed and orthonormal. */
    if (halfDet >= 0.0) {
      compute_eigenvector_0(a00, a01, a02, a11, a12, a22, eigenvalues[2],
                            eigenvector2);
      compute_eigenvector_1(a00, a01, a02, a11, a12, a22, eigenvector2,
                            eigenvalues[1], eigenvector1);
      chemistry_utils_cross_product(eigenvector1, eigenvector2, eigenvector0);
    } else {
      compute_eigenvector_0(a00, a01, a02, a11, a12, a22, eigenvalues[0],
                            eigenvector0);
      compute_eigenvector_1(a00, a01, a02, a11, a12, a22, eigenvector0,
                            eigenvalues[1], eigenvector1);
      chemistry_utils_cross_product(eigenvector0, eigenvector1, eigenvector2);
    }
  } else /* The matrix is diagonal */ {
    eigenvalues[0] = A_p[0][0];
    eigenvalues[1] = A_p[1][1];
    eigenvalues[2] = A_p[2][2];
    eigenvector0[0] = 1.0;
    eigenvector0[1] = 0.0;
    eigenvector0[2] = 0.0;
    eigenvector1[0] = 0.0;
    eigenvector1[1] = 1.0;
    eigenvector1[2] = 0.0;
    eigenvector2[0] = 0.0;
    eigenvector2[1] = 0.0;
    eigenvector2[2] = 1.0;
  }

  /* The preconditioning scaled the matrix A, which scales the
     eigenvalues. Revert the scaling. */
  eigenvalues[0] *= maxAbsElement;
  eigenvalues[1] *= maxAbsElement;
  eigenvalues[2] *= maxAbsElement;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_UTILS_H */
