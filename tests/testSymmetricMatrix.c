/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/**
 * @file testSymmetricMatrix.c
 * @brief Verify the matrix functions of dimension.h and symmetric_matrix.h
 * against the equivalent GSL calls.
 *
 * Matrices are built with a prescribed spectrum (M = Q diag(lambda) Q^T, or
 * U diag(sigma) V^T for general matrices) such that the condition number is
 * controlled, and are then scaled over many orders of magnitude to verify that
 * the results do not depend on the unit system. Inputs of the float functions
 * are rounded to float before being handed to either SWIFT or the reference
 * so that both see exactly the same matrix.
 *
 * Every inverse is verified via its residual |A^-1 A - I|. When GSL is
 * available, all results are additionally compared to the equivalent GSL
 * calls (SVD, LU, BLAS). Without GSL, condition numbers are compared to the
 * value the matrices were constructed with and products to plain
 * double-precision loops.
 *
 * Floating-point exceptions (division by zero, invalid, overflow) are trapped
 * when supported, except around the one check feeding a NaN on purpose.
 *
 * Usage: testSymmetricMatrix [seed]
 */
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local headers */
#include "dimension.h"
#include "error.h"
#include "symmetric_matrix.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#endif

/* Dimension of the symmetric matrices */
#define DIM hydro_dimension_integer

/* Number of random matrices per (scale, condition number) combination */
#define NUM_TRIALS 50

/* Maximal number of failures reported in detail */
#define MAX_REPORTS 20

static int num_checks = 0;
static int num_failures = 0;

/**
 * @brief Record the outcome of a check and report the first few failures.
 */
static void check(const int ok, const char *test, const char *what,
                  const double value, const double tolerance) {
  num_checks++;
  if (ok) return;
  num_failures++;
  if (num_failures <= MAX_REPORTS)
    message("FAIL [%s] %s: value=%.6e tolerance=%.6e", test, what, value,
            tolerance);
}

/**
 * @brief Switch the trapping of floating-point exceptions on or off.
 */
static void fpe_traps(const int on) {
#ifdef HAVE_FE_ENABLE_EXCEPT
  const int excepts = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
  feclearexcept(FE_ALL_EXCEPT);
  if (on)
    feenableexcept(excepts);
  else
    fedisableexcept(excepts);
#endif
}

/* ---------------------------------------------------------------------------
 * Random numbers and matrix construction
 * ------------------------------------------------------------------------ */

static double rand_uniform(const double a, const double b) {
  return a + (b - a) * drand48();
}

static double rand_gaussian(void) {
  const double u1 = 1. - drand48();
  const double u2 = drand48();
  return sqrt(-2. * log(u1)) * cos(2. * M_PI * u2);
}

/**
 * @brief Random orthogonal n x n matrix (Gram-Schmidt on Gaussian vectors).
 *
 * Two orthogonalisation passes ("twice is enough") make Q orthogonal to
 * machine precision, so that the constructed spectrum is exact to ~eps.
 */
static void random_orthogonal(const int n, double Q[3][3]) {
  for (int i = 0; i < n; ++i) {
    double norm;
    do {
      for (int k = 0; k < n; ++k) Q[k][i] = rand_gaussian();
      for (int pass = 0; pass < 2; ++pass) {
        for (int j = 0; j < i; ++j) {
          double dot = 0.;
          for (int k = 0; k < n; ++k) dot += Q[k][i] * Q[k][j];
          for (int k = 0; k < n; ++k) Q[k][i] -= dot * Q[k][j];
        }
      }
      norm = 0.;
      for (int k = 0; k < n; ++k) norm += Q[k][i] * Q[k][i];
      norm = sqrt(norm);
    } while (norm < 1e-3);
    for (int k = 0; k < n; ++k) Q[k][i] /= norm;
  }
}

/**
 * @brief Spectrum with max/min ratio exactly cond and the others log-uniform
 * in between. Signs are random if indefinite is set.
 */
static void random_spectrum(const int n, const double cond,
                            const int indefinite, double lambda[3]) {
  for (int i = 0; i < n; ++i) lambda[i] = pow(cond, drand48());
  lambda[0] = 1.;
  if (n > 1) lambda[n - 1] = cond;
  for (int i = 0; i < n; ++i)
    if (indefinite && drand48() < 0.5) lambda[i] = -lambda[i];
}

/**
 * @brief General n x n matrix U diag(sigma) V^T * scale, optionally rounded
 * to float.
 */
static void random_general_matrix(const int n, const double cond,
                                  const double scale, const int round,
                                  double M[3][3]) {
  double U[3][3], V[3][3], sigma[3];
  random_orthogonal(n, U);
  random_orthogonal(n, V);
  random_spectrum(n, cond, /*indefinite=*/0, sigma);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m = 0.;
      for (int k = 0; k < n; ++k) m += U[i][k] * sigma[k] * V[j][k];
      M[i][j] = round ? (float)(scale * m) : scale * m;
    }
  }
}

/**
 * @brief Symmetric n x n matrix Q diag(lambda) Q^T * scale, optionally rounded
 * to float.
 */
static void random_symmetric_matrix(const int n, const double cond,
                                    const int indefinite, const double scale,
                                    const int round, double M[3][3]) {
  double Q[3][3], lambda[3];
  random_orthogonal(n, Q);
  random_spectrum(n, cond, indefinite, lambda);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      double m = 0.;
      for (int k = 0; k < n; ++k) m += Q[i][k] * lambda[k] * Q[j][k];
      M[i][j] = M[j][i] = round ? (float)(scale * m) : scale * m;
    }
  }
}

/* ---------------------------------------------------------------------------
 * Conversions and GSL reference implementations
 * ------------------------------------------------------------------------ */

static void to_sym_matrix(const double M[3][3], struct sym_matrix *S) {
  float A[DIM][DIM];
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j) A[i][j] = (float)M[i][j];
  get_sym_matrix_from_matrix(S, A);
}

static void from_sym_matrix(const struct sym_matrix *S, double M[3][3]) {
  float A[DIM][DIM];
  get_matrix_from_sym_matrix(A, S);
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j) M[i][j] = A[i][j];
}

#ifdef HAVE_LIBGSL

/**
 * @brief max_ij |A - B| / max_ij |B|
 */
static double rel_max_diff(const int n, const double A[3][3],
                           const double B[3][3]) {
  double diff = 0., norm = 0.;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      diff = fmax(diff, fabs(A[i][j] - B[i][j]));
      norm = fmax(norm, fabs(B[i][j]));
    }
  }
  return norm > 0. ? diff / norm : diff;
}

static gsl_matrix *to_gsl(const int n, const double M[3][3]) {
  gsl_matrix *G = gsl_matrix_alloc(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) gsl_matrix_set(G, i, j, M[i][j]);
  return G;
}

/**
 * @brief 2-norm condition number from the GSL SVD.
 */
static double gsl_condition_number(const int n, const double M[3][3]) {
  gsl_matrix *A = to_gsl(n, M);
  gsl_matrix *V = gsl_matrix_alloc(n, n);
  gsl_vector *S = gsl_vector_alloc(n);
  gsl_vector *work = gsl_vector_alloc(n);
  gsl_linalg_SV_decomp(A, V, S, work);
  const double s_max = gsl_vector_get(S, 0);
  const double s_min = gsl_vector_get(S, n - 1);
  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_vector_free(S);
  gsl_vector_free(work);
  return s_min > 0. ? s_max / s_min : INFINITY;
}

/**
 * @brief Inverse via GSL LU decomposition. Returns non-zero on failure.
 */
static int gsl_inverse(const int n, const double M[3][3], double inv[3][3]) {
  gsl_matrix *A = to_gsl(n, M);
  gsl_matrix *Ainv = gsl_matrix_alloc(n, n);
  gsl_permutation *p = gsl_permutation_alloc(n);
  int signum;
  int res = gsl_linalg_LU_decomp(A, p, &signum);
  if (res == GSL_SUCCESS) res = gsl_linalg_LU_invert(A, p, Ainv);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) inv[i][j] = gsl_matrix_get(Ainv, i, j);
  gsl_matrix_free(A);
  gsl_matrix_free(Ainv);
  gsl_permutation_free(p);
  return res;
}

#endif /* HAVE_LIBGSL */

/**
 * @brief Reference condition number: the GSL SVD if available, otherwise the
 * value the matrix was constructed with (always 1 for a 1x1 matrix, see
 * random_spectrum()).
 */
static double ref_condition_number(const int n, const double M[3][3],
                                   const double constructed_cond) {
#ifdef HAVE_LIBGSL
  return gsl_condition_number(n, M);
#else
  return n == 1 ? 1. : constructed_cond;
#endif
}

#if defined(HYDRO_DIMENSION_3D)
/**
 * @brief Reference product C = A B (GSL dgemm if available).
 */
static void ref_matmul(const int n, const double A[3][3], const double B[3][3],
                       double C[3][3]) {
#ifdef HAVE_LIBGSL
  gsl_matrix *GA = to_gsl(n, A);
  gsl_matrix *GB = to_gsl(n, B);
  gsl_matrix *GC = gsl_matrix_alloc(n, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., GA, GB, 0., GC);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) C[i][j] = gsl_matrix_get(GC, i, j);
  gsl_matrix_free(GA);
  gsl_matrix_free(GB);
  gsl_matrix_free(GC);
#else
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      C[i][j] = 0.;
      for (int k = 0; k < n; ++k) C[i][j] += A[i][k] * B[k][j];
    }
  }
#endif
}
#endif /* HYDRO_DIMENSION_3D */

/**
 * @brief Reference product y = M x for symmetric M (GSL dsymv if available).
 */
static void ref_symv(const int n, const double M[3][3], const double x[3],
                     double y[3]) {
#ifdef HAVE_LIBGSL
  gsl_matrix *G = to_gsl(n, M);
  gsl_vector *gx = gsl_vector_alloc(n);
  gsl_vector *gy = gsl_vector_alloc(n);
  for (int i = 0; i < n; ++i) gsl_vector_set(gx, i, x[i]);
  gsl_blas_dsymv(CblasUpper, 1., G, gx, 0., gy);
  for (int i = 0; i < n; ++i) y[i] = gsl_vector_get(gy, i);
  gsl_matrix_free(G);
  gsl_vector_free(gx);
  gsl_vector_free(gy);
#else
  for (int i = 0; i < n; ++i) {
    y[i] = 0.;
    for (int j = 0; j < n; ++j) y[i] += M[i][j] * x[j];
  }
#endif
}

/**
 * @brief Verify an inverse: residual max_ij |(inv M - I)_ij| <= res_tol, and,
 * when GSL is available, max_ij |inv - inv_GSL| / max_ij |inv_GSL| <= gsl_tol.
 */
static void check_inverse(const char *name, const int n, const double M[3][3],
                          const double inv[3][3], const double res_tol,
                          const double gsl_tol) {
  double residual = 0.;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double r = (i == j) ? -1. : 0.;
      for (int k = 0; k < n; ++k) r += inv[i][k] * M[k][j];
      residual = fmax(residual, fabs(r));
    }
  }
  check(residual <= res_tol, name, "residual |inv M - I|", residual, res_tol);

#ifdef HAVE_LIBGSL
  double ref[3][3];
  const int res = gsl_inverse(n, M, ref);
  check(res == GSL_SUCCESS, name, "GSL inversion failed", res, 0.);
  const double err = rel_max_diff(n, inv, ref);
  check(err <= gsl_tol, name, "relative error vs. GSL LU", err, gsl_tol);
#endif
}

/* Scales used to verify the unit-independence of the results. Chosen such
 * that the matrices and their inverses remain normal floats. */
static const double scales[] = {1e-15, 1e-10, 1e-5, 1., 1e5, 1e10, 1e15};
static const int num_scales = sizeof(scales) / sizeof(scales[0]);

static const double conds[] = {1., 2., 10., 1e2, 1e3, 1e4, 1e5};
static const int num_conds = sizeof(conds) / sizeof(conds[0]);

/* ---------------------------------------------------------------------------
 * Tests of the double-precision 3x3 helpers of dimension.h
 * ------------------------------------------------------------------------ */

/**
 * @brief matrix_3x3_2norm_condition_number() vs. gsl_linalg_SV_decomp().
 */
static void test_condition_number(void) {

  const char *name = "cond_3x3";

  for (int s = 0; s < num_scales; ++s) {
    for (int c = 0; c < num_conds; ++c) {
      for (int t = 0; t < NUM_TRIALS; ++t) {
        double M[3][3];
        if (t % 2)
          random_general_matrix(3, conds[c], scales[s], /*round=*/0, M);
        else
          random_symmetric_matrix(3, conds[c], t % 4 == 0, scales[s],
                                  /*round=*/0, M);

        const double ref = ref_condition_number(3, M, conds[c]);
        const double val = matrix_3x3_2norm_condition_number(M);

        /* Working with M^T M implies a relative error ~eps * cond^2 */
        const double tol = 1e-13 + 1e-14 * ref * ref;
        const double err = fabs(val - ref) / ref;
        check(err <= tol, name, "relative error vs. reference", err, tol);
      }
    }
  }

  /* Singular matrices must be reported as (at least) very badly
   * conditioned, whatever their scale. */
  for (int s = 0; s < num_scales; ++s) {
    const double a = scales[s];
    const double zero[3][3] = {{0.}};
    const double rank1[3][3] = {{a, 2. * a, 3. * a},
                                {2. * a, 4. * a, 6. * a},
                                {3. * a, 6. * a, 9. * a}};
    const double rank2[3][3] = {
        {a, 2. * a, 3. * a}, {2. * a, 4. * a, 6. * a}, {a, 0., a}};
    /* Not isinf(), which is folded away under -ffast-math */
    check(matrix_3x3_2norm_condition_number(zero) > 1e300, name,
          "zero matrix is not reported as singular", a, 0.);
    check(matrix_3x3_2norm_condition_number(rank1) > 1e7, name,
          "rank-1 matrix is not reported as singular", a, 1e7);
    check(matrix_3x3_2norm_condition_number(rank2) > 1e7, name,
          "rank-2 matrix is not reported as singular", a, 1e7);
  }
}

/**
 * @brief matrix_3x3_symmetric_2norm_condition_number() vs.
 * gsl_linalg_SV_decomp().
 */
static void test_symmetric_condition_number(void) {

  const char *name = "sym_cond_3x3";
  const double high_conds[] = {1e6, 1e7};

  for (int s = 0; s < num_scales; ++s) {
    for (int c = 0; c < num_conds + 2; ++c) {
      const double cond = c < num_conds ? conds[c] : high_conds[c - num_conds];
      for (int t = 0; t < NUM_TRIALS; ++t) {
        double M[3][3];
        random_symmetric_matrix(3, cond, t % 2, scales[s], /*round=*/0, M);

        const double ref = ref_condition_number(3, M, cond);
        const double val = matrix_3x3_symmetric_2norm_condition_number(M);

        /* Eigenvalues of M itself: relative error ~eps * cond */
        const double tol = 1e-13 + 1e-14 * ref;
        const double err = fabs(val - ref) / ref;
        check(err <= tol, name, "relative error vs. reference", err, tol);
      }
    }
  }

  /* Exactly degenerate spectra (worst case of closed-form solutions) */
  for (int s = 0; s < num_scales; ++s) {
    const double a = scales[s];
    const double Id[3][3] = {{a, 0., 0.}, {0., a, 0.}, {0., 0., a}};
    const double deg2[3][3] = {{2. * a, a, 0.}, {a, 2. * a, 0.}, {0., 0., a}};
    const double v1 = matrix_3x3_symmetric_2norm_condition_number(Id);
    const double v2 = matrix_3x3_symmetric_2norm_condition_number(deg2);
    const double v3 = matrix_3x3_2norm_condition_number(deg2);
    check(fabs(v1 - 1.) <= 1e-15, name, "identity", v1 - 1., 1e-15);
    check(fabs(v2 - 3.) <= 3e-15, name, "double eigenvalue", v2 - 3., 3e-15);
    check(fabs(v3 - 3.) <= 3e-15, "cond_3x3", "double eigenvalue", v3 - 3.,
          3e-15);
  }

  /* Singular matrices must be reported as such, whatever their scale. */
  for (int s = 0; s < num_scales; ++s) {
    const double a = scales[s];
    const double zero[3][3] = {{0.}};
    const double rank1[3][3] = {{a, 2. * a, 3. * a},
                                {2. * a, 4. * a, 6. * a},
                                {3. * a, 6. * a, 9. * a}};
    const double rank2[3][3] = {{a, a, 0.}, {a, a, 0.}, {0., 0., a}};
    check(matrix_3x3_symmetric_2norm_condition_number(zero) > 1e300, name,
          "zero matrix is not reported as singular", a, 0.);
    check(matrix_3x3_symmetric_2norm_condition_number(rank1) > 1e14, name,
          "rank-1 matrix is not reported as singular", a, 1e14);
    check(matrix_3x3_symmetric_2norm_condition_number(rank2) > 1e14, name,
          "rank-2 matrix is not reported as singular", a, 1e14);
  }
}

/**
 * @brief invert3x3_matrix_LU() vs. gsl_linalg_LU_decomp/invert().
 */
static void test_LU_inverse(void) {

  const char *name = "LU_3x3";

  for (int s = 0; s < num_scales; ++s) {
    for (int c = 0; c < num_conds; ++c) {
      for (int t = 0; t < NUM_TRIALS; ++t) {
        double M[3][3], inv[3][3];
        if (t % 2)
          random_general_matrix(3, conds[c], scales[s], /*round=*/0, M);
        else
          random_symmetric_matrix(3, conds[c], t % 4 == 0, scales[s],
                                  /*round=*/0, M);

        const int res = invert3x3_matrix_LU(M, inv);
        check(res == 0, name, "inversion failed", res, 0.);
        if (res == 0) {
          const double tol = 1e-13 * conds[c] + 1e-14;
          check_inverse(name, 3, M, inv, tol, tol);
        }
      }
    }
  }

  /* Matrices that require pivoting */
  const double perm[3][3] = {{0., 1., 0.}, {1., 0., 0.}, {0., 0., 1.}};
  const double tiny_pivot[3][3] = {{1e-20, 1., 0.}, {1., 1., 0.}, {0., 0., 1.}};
  const double *pivots[2] = {&perm[0][0], &tiny_pivot[0][0]};
  for (int k = 0; k < 2; ++k) {
    double M[3][3], inv[3][3];
    for (int i = 0; i < 9; ++i) M[i / 3][i % 3] = pivots[k][i];
    const int res = invert3x3_matrix_LU(M, inv);
    check(res == 0, name, "pivoting case failed", res, 0.);
    if (res == 0) check_inverse(name, 3, M, inv, 1e-14, 1e-14);
  }

  /* Singular matrices must be rejected, whatever their scale. */
  for (int s = 0; s < num_scales; ++s) {
    const double a = scales[s];
    double inv[3][3];
    const double zero[3][3] = {{0.}};
    const double zero_row[3][3] = {{a, 2. * a, 0.}, {0., 0., 0.}, {a, 0., a}};
    const double dep_rows[3][3] = {
        {a, 2. * a, 3. * a}, {2. * a, 4. * a, 6. * a}, {a, 0., a}};
    check(invert3x3_matrix_LU(zero, inv) == 1, name, "zero matrix not rejected",
          a, 0.);
    check(invert3x3_matrix_LU(zero_row, inv) == 1, name,
          "zero-row matrix not rejected", a, 0.);
    check(invert3x3_matrix_LU(dep_rows, inv) == 1, name,
          "rank-deficient matrix not rejected", a, 0.);
  }
}

/* ---------------------------------------------------------------------------
 * Tests of the (float) sym_matrix functions in the compiled dimension
 * ------------------------------------------------------------------------ */

/**
 * @brief sym_matrix_invert() vs. GSL SVD condition number and LU inverse.
 *
 * Matrices well inside the max_cond_num limit must be inverted accurately,
 * matrices well outside must be rejected and return the null matrix.
 */
static void test_sym_matrix_invert(void) {

  const char *name = "sym_matrix_invert";
  const double max_conds[] = {10., 60., 1e3, 1e6};

  for (int m = 0; m < 4; ++m) {
    const double max_cond = max_conds[m];
    for (int s = 0; s < num_scales; ++s) {
      for (int c = 0; c < num_conds; ++c) {
        for (int t = 0; t < NUM_TRIALS; ++t) {
          double M[3][3], inv[3][3];
          random_symmetric_matrix(DIM, conds[c], t % 2, scales[s],
                                  /*round=*/1, M);

          struct sym_matrix S, S_inv;
          to_sym_matrix(M, &S);
          for (int i = 0; i < sym_matrix_num_elements; ++i)
            S_inv.elements[i] = 42.f;

          const int res = sym_matrix_invert(&S_inv, &S, max_cond);

          /* Without GSL, the constructed value is off by ~eps_float * cond
           * because of the rounding, which the factors 2 below absorb. */
          const double ref_cond = ref_condition_number(DIM, M, conds[c]);

          if (ref_cond < 0.5 * max_cond) {
            check(res == 0, name, "well-conditioned matrix rejected", ref_cond,
                  max_cond);
            if (res == 0) {
              from_sym_matrix(&S_inv, inv);
              /* Output is stored in float */
              check_inverse(name, DIM, M, inv, 1e-6 * ref_cond,
                            1e-6 + 1e-13 * ref_cond);
            }
          } else if (ref_cond > 2. * max_cond) {
            check(res == 1, name, "ill-conditioned matrix accepted", ref_cond,
                  max_cond);
            check(sym_matrix_is_null(&S_inv), name,
                  "M_inv not null after failure", ref_cond, max_cond);
          }
        }
      }
    }
  }

  /* Degenerate inputs must fail and return the null matrix. */
  struct sym_matrix S, S_inv;
  zero_sym_matrix(&S);
  S_inv.xx = 42.f;
  check(sym_matrix_invert(&S_inv, &S, 1e6) == 1 && sym_matrix_is_null(&S_inv),
        name, "zero matrix", 0., 0.);

  /* The NaN is deliberate here: no traps */
  fpe_traps(0);
  sym_matrix_identity(&S);
  S.xx = NAN;
  S_inv.xx = 42.f;
  check(sym_matrix_invert(&S_inv, &S, 1e6) == 1 && sym_matrix_is_null(&S_inv),
        name, "NaN matrix", 0., 0.);
  fpe_traps(1);

#if DIM > 1
  /* Exactly singular (rank-deficient) matrices at all scales */
  for (int s = 0; s < num_scales; ++s) {
    double M[3][3] = {{0.}};
    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j) M[i][j] = scales[s] * (i + 1) * (j + 1);
    to_sym_matrix(M, &S);
    S_inv.xx = 42.f;
    check(sym_matrix_invert(&S_inv, &S, 1e6) == 1 && sym_matrix_is_null(&S_inv),
          name, "rank-1 matrix", scales[s], 0.);
  }
#endif
}

/**
 * @brief sym_matrix_multiply_by_vector() vs. gsl_blas_dsymv().
 */
static void test_sym_matrix_multiply_by_vector(void) {

  const char *name = "sym_matrix_mult_vec";

  for (int s = 0; s < num_scales; ++s) {
    for (int t = 0; t < NUM_TRIALS; ++t) {
      double M[3][3];
      random_symmetric_matrix(DIM, 1e3, t % 2, scales[s], /*round=*/1, M);
      struct sym_matrix S;
      to_sym_matrix(M, &S);

      float v[DIM], out[DIM];
      double x[3], y[3];
      for (int i = 0; i < DIM; ++i) {
        v[i] = (float)rand_uniform(-1., 1.);
        x[i] = v[i];
      }

      sym_matrix_multiply_by_vector(out, &S, v);
      ref_symv(DIM, M, x, y);

      /* Float arithmetic: error relative to sum_j |M_ij v_j| */
      for (int i = 0; i < DIM; ++i) {
        double norm = 0.;
        for (int j = 0; j < DIM; ++j) norm += fabs(M[i][j] * v[j]);
        const double err = fabs(out[i] - y[i]) / norm;
        check(err <= 1e-6, name, "relative error vs. reference", err, 1e-6);
      }
    }
  }
}

/**
 * @brief sym_matrix_multiplication_ABA() vs. two gsl_blas_dgemm() calls.
 */
static void test_sym_matrix_multiplication_ABA(void) {
#if defined(HYDRO_DIMENSION_3D)

  const char *name = "sym_matrix_ABA";

  for (int s = 0; s < num_scales; ++s) {
    for (int t = 0; t < NUM_TRIALS; ++t) {
      double A[3][3], B[3][3], out[3][3], BA[3][3], ABA[3][3];
      random_symmetric_matrix(3, 1e2, t % 2, 1., /*round=*/1, A);
      random_symmetric_matrix(3, 1e2, (t / 2) % 2, scales[s], /*round=*/1, B);

      struct sym_matrix SA, SB, S_out;
      to_sym_matrix(A, &SA);
      to_sym_matrix(B, &SB);
      sym_matrix_multiplication_ABA(&S_out, &SA, &SB);
      from_sym_matrix(&S_out, out);

      ref_matmul(3, B, A, BA);
      ref_matmul(3, A, BA, ABA);

      /* Float arithmetic: error relative to |A| |B| |A| */
      double norm_A = 0., norm_B = 0., err = 0.;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          norm_A = fmax(norm_A, fabs(A[i][j]));
          norm_B = fmax(norm_B, fabs(B[i][j]));
          err = fmax(err, fabs(out[i][j] - ABA[i][j]));
        }
      }
      err /= norm_A * norm_B * norm_A;
      check(err <= 1e-5, name, "relative error vs. reference", err, 1e-5);
    }
  }
#endif
}

/**
 * @brief invert_dimension_by_dimension_matrix() vs. GSL LU inverse.
 */
static void test_invert_dimension_by_dimension(void) {

  const char *name = "invert_dim_by_dim";

  for (int s = 0; s < num_scales; ++s) {
    for (int c = 0; c < 5; ++c) { /* cond <= 1e3: the inversion is in float */
      for (int t = 0; t < NUM_TRIALS; ++t) {
        double M[3][3] = {{0.}}, inv[3][3] = {{0.}};
        random_general_matrix(DIM, conds[c], scales[s], /*round=*/1, M);

        /* The top-left DIM x DIM part of a 3x3 array is inverted */
        float A[3][3] = {{0.f}};
        for (int i = 0; i < DIM; ++i)
          for (int j = 0; j < DIM; ++j) A[i][j] = (float)M[i][j];

        const int res = invert_dimension_by_dimension_matrix(A, 1e-8f);
        check(res == 0, name, "inversion failed", conds[c], 0.);

        if (res == 0) {
          for (int i = 0; i < DIM; ++i)
            for (int j = 0; j < DIM; ++j) inv[i][j] = A[i][j];
          check_inverse(name, DIM, M, inv, 1e-5 * conds[c], 1e-5 * conds[c]);
        }
      }
    }
  }
}

/**
 * @brief Conversions, identity, null test and scalar multiplication.
 */
static void test_sym_matrix_basics(void) {

  const char *name = "sym_matrix_basics";

  for (int t = 0; t < NUM_TRIALS; ++t) {
    double M[3][3], back[3][3];
    random_symmetric_matrix(DIM, 10., t % 2, 1., /*round=*/1, M);

    /* Round-trip must be exact */
    struct sym_matrix S;
    to_sym_matrix(M, &S);
    from_sym_matrix(&S, back);
    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        check(back[i][j] == M[i][j], name, "round-trip not exact",
              back[i][j] - M[i][j], 0.);

    /* Scalar multiplication vs. the double-precision product */
    const float alpha = (float)rand_uniform(-10., 10.);
    sym_matrix_multiply_by_scalar(&S, alpha);
    from_sym_matrix(&S, back);
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < DIM; ++j) {
        const double ref = alpha * M[i][j];
        const double err = fabs(back[i][j] - ref) / fabs(ref);
        check(err <= 1e-6, name, "scalar multiplication", err, 1e-6);
      }
    }
  }

  struct sym_matrix I;
  sym_matrix_identity(&I);
  double Id[3][3];
  from_sym_matrix(&I, Id);
  for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
      check(Id[i][j] == (i == j ? 1. : 0.), name, "identity", Id[i][j], 0.);
  check(!sym_matrix_is_null(&I), name, "identity reported as null", 0., 0.);

  zero_sym_matrix(&I);
  check(sym_matrix_is_null(&I), name, "zero matrix not reported as null", 0.,
        0.);
}

int main(int argc, char *argv[]) {

  const long seed = argc > 1 ? atol(argv[1]) : 1234;
  srand48(seed);

#ifdef HAVE_LIBGSL
  /* We check the error codes of GSL ourselves */
  gsl_set_error_handler_off();
  const char *reference = "GSL";
#else
  const char *reference = "residuals only (no GSL)";
#endif

#ifdef HAVE_FE_ENABLE_EXCEPT
  const char *traps = "on";
#else
  const char *traps = "unsupported";
#endif

  message("Testing matrix functions in %dD against %s (seed=%ld, FPE traps %s)",
          DIM, reference, seed, traps);

  /* Choke on FPEs */
  fpe_traps(1);

  test_condition_number();
  test_symmetric_condition_number();
  test_LU_inverse();
  test_sym_matrix_invert();
  test_sym_matrix_multiply_by_vector();
  test_sym_matrix_multiplication_ABA();
  test_invert_dimension_by_dimension();
  test_sym_matrix_basics();

  message("%d checks performed, %d failures.", num_checks, num_failures);
  return num_failures ? 1 : 0;
}
