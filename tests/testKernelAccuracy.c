/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
 * @file testKernelAccuracy.c
 * @brief Compares the tabulated kernels of kernel_hydro.h with the exact
 * analytic expressions of Dehnen & Aly (2012), including the edge of the
 * support, the double-precision version and the hand-vectorised versions.
 *
 * The reference uses the (float) kernel_gamma of the header as the support,
 * which is the convention of the tables.
 */

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "error.h"
#include "kernel_hydro.h"
#include "vector.h"

/* System includes. */
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Dimension as an integer */
#if defined(HYDRO_DIMENSION_3D)
#define DIM 3
#elif defined(HYDRO_DIMENSION_2D)
#define DIM 2
#elif defined(HYDRO_DIMENSION_1D)
#define DIM 1
#endif

/* Tolerances */
static const double tol_abs_W = 1e-6;    /* |W - W_ref| in units of W(0) */
static const double tol_abs_dW = 1e-6;   /* |dW - dW_ref| in units of max|dW| */
static const double tol_rel_edge = 1e-4; /* relative, for u >= 0.9 gamma */
static const double tol_dbl = 1e-13;     /* kernel_eval_double, units of W(0) */
static const double tol_dbl_edge = 1e-10; /* relative, for u >= 0.9 gamma */
static const float tol_vec_ulps = 8.f;    /* between implementations, ulps */
static const double tol_norm = 1e-5;      /* integral of W over the support */
static const double tiny = 1e-30; /* smallest value checked relatively */

/* ------------------------------------------------------------------------- */
/* Dehnen & Aly (2012), table 1: normalisation C and shape f(q), q = r / H.  */
/* Written in terms of q and s = 1 - q (passed exactly to avoid cancellation */
/* near the edge of the support).                                            */

static double pos(const double x) { return x > 0. ? x : 0.; }

#if defined(CUBIC_SPLINE_KERNEL)
#if DIM == 3
static const double C_ref = 16. * M_1_PI;
#elif DIM == 2
static const double C_ref = 80. * M_1_PI / 7.;
#else
static const double C_ref = 8. / 3.;
#endif
static double f_ref(const double q, const double s) {
  return pow(pos(s), 3) - 4. * pow(pos(0.5 - q), 3);
}
static double fp_ref(const double q, const double s) {
  return -3. * pow(pos(s), 2) + 12. * pow(pos(0.5 - q), 2);
}

#elif defined(QUARTIC_SPLINE_KERNEL)
#if DIM == 3
static const double C_ref = 15625. * M_1_PI / 512.;
#elif DIM == 2
static const double C_ref = 46875. * M_1_PI / 2398.;
#else
static const double C_ref = 3125. / 768.;
#endif
static double f_ref(const double q, const double s) {
  return pow(pos(s), 4) - 5. * pow(pos(0.6 - q), 4) +
         10. * pow(pos(0.2 - q), 4);
}
static double fp_ref(const double q, const double s) {
  return -4. * pow(pos(s), 3) + 20. * pow(pos(0.6 - q), 3) -
         40. * pow(pos(0.2 - q), 3);
}

#elif defined(QUINTIC_SPLINE_KERNEL)
#if DIM == 3
static const double C_ref = 2187. * M_1_PI / 40.;
#elif DIM == 2
static const double C_ref = 15309. * M_1_PI / 478.;
#else
static const double C_ref = 243. / 40.;
#endif
static double f_ref(const double q, const double s) {
  return pow(pos(s), 5) - 6. * pow(pos(2. / 3. - q), 5) +
         15. * pow(pos(1. / 3. - q), 5);
}
static double fp_ref(const double q, const double s) {
  return -5. * pow(pos(s), 4) + 30. * pow(pos(2. / 3. - q), 4) -
         75. * pow(pos(1. / 3. - q), 4);
}

#elif defined(WENDLAND_C2_KERNEL)
#if DIM == 3
static const double C_ref = 21. * M_1_PI / 2.;
#elif DIM == 2
static const double C_ref = 7. * M_1_PI;
#else
static const double C_ref = 5. / 4.;
#endif
static double f_ref(const double q, const double s) {
#if DIM == 1
  return pow(pos(s), 3) * (1. + 3. * q);
#else
  return pow(pos(s), 4) * (1. + 4. * q);
#endif
}
static double fp_ref(const double q, const double s) {
#if DIM == 1
  return -12. * q * pow(pos(s), 2);
#else
  return -20. * q * pow(pos(s), 3);
#endif
}

#elif defined(WENDLAND_C4_KERNEL)
#if DIM == 3
static const double C_ref = 495. * M_1_PI / 32.;
#elif DIM == 2
static const double C_ref = 9. * M_1_PI;
#else
#error "Wendland C4 kernel not defined in 1D."
#endif
static double f_ref(const double q, const double s) {
  return pow(pos(s), 6) * (1. + 6. * q + 35. / 3. * q * q);
}
static double fp_ref(const double q, const double s) {
  return -56. / 3. * q * pow(pos(s), 5) * (1. + 5. * q);
}

#elif defined(WENDLAND_C6_KERNEL)
#if DIM == 3
static const double C_ref = 1365. * M_1_PI / 64.;
#elif DIM == 2
static const double C_ref = 78. * M_1_PI / 7.;
#else
#error "Wendland C6 kernel not defined in 1D."
#endif
static double f_ref(const double q, const double s) {
  return pow(pos(s), 8) * (1. + 8. * q + 25. * q * q + 32. * q * q * q);
}
static double fp_ref(const double q, const double s) {
  return -22. * q * pow(pos(s), 7) * (1. + 7. * q + 16. * q * q);
}
#endif

/* The exact kernel and its derivative for the float gamma of the header. */
static const double gamma_d = (double)kernel_gamma;

static double W_ref(const double u) {
  const double q = u / gamma_d;
  const double s = (gamma_d - u) / gamma_d;
  return C_ref / pow(gamma_d, DIM) * f_ref(q, s);
}

static double dW_ref(const double u) {
  const double q = u / gamma_d;
  const double s = (gamma_d - u) / gamma_d;
  return C_ref / pow(gamma_d, DIM + 1) * fp_ref(q, s);
}

/* ------------------------------------------------------------------------- */

/* Number of grid points on [0, 1.05 gamma] */
#define N_GRID 200001

/* Edge-approach points u = gamma * (1 - 2^-k) */
#define N_EDGE 24

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  message("Kernel: %s, %dD, gamma=%.9g, %d sub-intervals of degree %d",
          kernel_name, DIM, kernel_gamma, kernel_poly_ivals,
          kernel_poly_degree);

  /* The header's normalisation constant must be the one of the reference */
  if (kernel_constant != (float)C_ref)
    error("kernel_constant=%.9g differs from the reference %.9g",
          kernel_constant, (float)C_ref);

  /* Scales */
  const double W0 = W_ref(0.);
  double dWmax = 0.;
  for (int i = 0; i < 1000; ++i) {
    const double u = gamma_d * i / 1000.;
    dWmax = fmax(dWmax, fabs(dW_ref(u)));
  }

  /* ------------------------------------------------------------------- */
  /* Build the list of evaluation points */
  const int N = N_GRID + N_EDGE + 6;
  float *u = (float *)malloc(N * sizeof(float));
  if (u == NULL) error("Error allocating u");
  for (int i = 0; i < N_GRID; ++i)
    u[i] = (float)(1.05 * gamma_d * i / (N_GRID - 1));
  for (int k = 0; k < N_EDGE; ++k)
    u[N_GRID + k] = (float)(gamma_d * (1. - pow(2., -(k + 1))));
  u[N_GRID + N_EDGE + 0] = nextafterf(kernel_gamma, 0.f);
  u[N_GRID + N_EDGE + 1] = kernel_gamma;
  u[N_GRID + N_EDGE + 2] = nextafterf(kernel_gamma, 2.f * kernel_gamma);
  u[N_GRID + N_EDGE + 3] = kernel_gamma * (1.f + 1e-6f);
  u[N_GRID + N_EDGE + 4] = 2.f * kernel_gamma;
  u[N_GRID + N_EDGE + 5] = 0.f;

  /* ------------------------------------------------------------------- */
  /* Scalar functions against the exact kernel */
  double max_err_W = 0., max_err_dW = 0., max_rel_W = 0., max_rel_dW = 0.;
  double max_err_Wd = 0., max_rel_Wd = 0.;
  float W_prev = FLT_MAX;

  for (int i = 0; i < N; ++i) {

    const float ui = u[i];
    const double Wr = W_ref(ui);
    const double dWr = dW_ref(ui);
    const int in_edge_zone = (ui >= 0.9f * kernel_gamma) && (ui < kernel_gamma);

    /* All scalar evaluations */
    float W, dW, W2, dW2;
    double Wd;
    kernel_deval(ui, &W, &dW);
    kernel_eval(ui, &W2);
    kernel_eval_dWdx(ui, &dW2);
    kernel_eval_double((double)ui, &Wd);

    /* Sign conditions */
    if (W < 0.f) error("Kernel is negative u=%e W=%e", ui, W);
    if (dW > 0.f) error("Kernel derivative is positive u=%e dW=%e", ui, dW);

    /* Zero outside the support (exactly) */
    if (ui >= kernel_gamma && (W != 0.f || dW != 0.f || Wd != 0.))
      error("Kernel not zero outside the support u=%e W=%e dW=%e Wd=%e", ui, W,
            dW, Wd);

    /* Absolute errors */
    const double err_W = fabs((double)W - Wr) / W0;
    const double err_dW = fabs((double)dW - dWr) / dWmax;
    const double err_Wd = fabs(Wd - Wr) / W0;
    max_err_W = fmax(max_err_W, err_W);
    max_err_dW = fmax(max_err_dW, err_dW);
    max_err_Wd = fmax(max_err_Wd, err_Wd);
    if (err_W > tol_abs_W)
      error("W inaccurate: u=%.9g W=%.9g W_ref=%.17g err/W0=%.3e", ui, W, Wr,
            err_W);
    if (err_dW > tol_abs_dW)
      error("dW inaccurate: u=%.9g dW=%.9g dW_ref=%.17g err/dWmax=%.3e", ui, dW,
            dWr, err_dW);
    if (err_Wd > tol_dbl)
      error("W (double) inaccurate: u=%.9g W=%.17g W_ref=%.17g err/W0=%.3e", ui,
            Wd, Wr, err_Wd);

    /* Relative errors near the edge of the support */
    if (in_edge_zone && Wr > tiny) {
      const double rel = fabs((double)W - Wr) / Wr;
      const double rel_d = fabs(Wd - Wr) / Wr;
      max_rel_W = fmax(max_rel_W, rel);
      max_rel_Wd = fmax(max_rel_Wd, rel_d);
      if (rel > tol_rel_edge)
        error(
            "W inaccurate near the edge: u=%.9g (gamma-u=%.3e) W=%.9g "
            "W_ref=%.17g rel=%.3e",
            ui, gamma_d - ui, W, Wr, rel);
      if (rel_d > tol_dbl_edge)
        error(
            "W (double) inaccurate near the edge: u=%.9g W=%.17g "
            "W_ref=%.17g rel=%.3e",
            ui, Wd, Wr, rel_d);
    }
    if (in_edge_zone && fabs(dWr) > tiny) {
      const double rel = fabs((double)dW - dWr) / fabs(dWr);
      max_rel_dW = fmax(max_rel_dW, rel);
      if (rel > tol_rel_edge)
        error(
            "dW inaccurate near the edge: u=%.9g (gamma-u=%.3e) dW=%.9g "
            "dW_ref=%.17g rel=%.3e",
            ui, gamma_d - ui, dW, dWr, rel);
    }

    /* Consistency between the scalar functions (same arithmetic; they can
     * only differ by the compiler's choice of contraction / association) */
    if (fabsf(W2 - W) > tol_vec_ulps * FLT_EPSILON * fmaxf(W, kernel_root))
      error("kernel_eval() and kernel_deval() disagree: u=%.9g %.9g %.9g", ui,
            W2, W);
    if (fabsf(dW2 - dW) >
        tol_vec_ulps * FLT_EPSILON * fmaxf(fabsf(dW), kernel_root))
      error("kernel_eval_dWdx() and kernel_deval() disagree: u=%.9g %.9g %.9g",
            ui, dW2, dW);

    /* Monotonicity along the grid */
    if (i < N_GRID) {
      if (W > W_prev + 2.f * FLT_EPSILON * (float)W0)
        error("Kernel not monotonic: u=%.9g W=%.9g W_prev=%.9g", ui, W, W_prev);
      W_prev = W;
    }
  }

  /* ------------------------------------------------------------------- */
  /* Special points */
  float W, dW;
  kernel_deval(0.f, &W, &dW);
  if (W != kernel_root) error("W(0)=%.9g != kernel_root=%.9g", W, kernel_root);
  if (dW != 0.f) error("dW(0)=%.9g != 0", dW);
  if (kernel_poly_index(0.f) != 0) error("Wrong sub-interval for u=0");
  if (kernel_poly_index(nextafterf(kernel_gamma, 0.f)) != kernel_poly_ivals - 1)
    error("Wrong sub-interval just below the edge");
  if (kernel_poly_index(kernel_gamma) != kernel_poly_ivals)
    error("Wrong sub-interval at the edge");
  if (kernel_poly_index(2.f * kernel_gamma) != kernel_poly_ivals)
    error("Wrong sub-interval beyond the edge");
  if (kernel_poly_index(1e30f) != kernel_poly_ivals)
    error("Wrong sub-interval for a huge u");

  /* ------------------------------------------------------------------- */
  /* Normalisation: integral of the tabulated W over its support */
  double integral = 0.;
  {
    const int M = 400000;
    const double du = gamma_d / M;
    for (int i = 0; i <= M; ++i) {
      const float ui = (float)(i * du);
      float Wi;
      kernel_eval(ui, &Wi);
      const double weight = (i == 0 || i == M) ? 0.5 : 1.;
#if DIM == 3
      integral += weight * 4. * M_PI * ui * ui * Wi;
#elif DIM == 2
      integral += weight * 2. * M_PI * ui * Wi;
#else
      integral += weight * 2. * Wi;
#endif
    }
    integral *= du;
  }
  if (fabs(integral - 1.) > tol_norm)
    error("Kernel not normalised: integral=%.9g", integral);

  message("Scalar: max |W-W_ref|/W0=%.2e, max |dW-dW_ref|/dWmax=%.2e",
          max_err_W, max_err_dW);
  message("Scalar: edge (u > 0.9 gamma) max relative error W=%.2e dW=%.2e",
          max_rel_W, max_rel_dW);
  message("Double: max |W-W_ref|/W0=%.2e, edge max relative error=%.2e",
          max_err_Wd, max_rel_Wd);
  message("Integral of W over the support: 1 %+.2e", integral - 1.);

  /* ------------------------------------------------------------------- */
  /* Hand-vectorised functions against the scalar ones */
#ifdef WITH_VECTORIZATION

  float max_diff_vec = 0.f;

  for (int i = 0; i + VEC_SIZE <= N; i += VEC_SIZE) {

    vector vu, vu2, vW, vdW, vW2, vdW2, vWb, vdWb, vdWc, vdWd, vdWe;
    for (int j = 0; j < VEC_SIZE; j++) {
      vu.f[j] = u[i + j];
      vu2.f[j] = u[N - 1 - i - j]; /* second vector: reversed order */
    }

    kernel_deval_1_vec(&vu, &vW, &vdW);
    kernel_deval_2_vec(&vu, &vW2, &vdW2, &vu2, &vWb, &vdWb);
    kernel_eval_W_vec(&vu, &vWb);
    kernel_eval_dWdx_vec(&vu, &vdWc);
    kernel_eval_dWdx_force_vec(&vu, &vdWd);
    kernel_eval_dWdx_force_2_vec(&vu, &vdWe, &vu2, &vdWb);

    for (int j = 0; j < VEC_SIZE; j++) {

      float Ws, dWs, Ws2, dWs2;
      kernel_deval(vu.f[j], &Ws, &dWs);
      kernel_deval(vu2.f[j], &Ws2, &dWs2);

      const float tol_W = tol_vec_ulps * FLT_EPSILON * fmaxf(Ws, kernel_root);
      const float tol_dW =
          tol_vec_ulps * FLT_EPSILON * fmaxf(fabsf(dWs), kernel_root);
      const float tol_dW2 =
          tol_vec_ulps * FLT_EPSILON * fmaxf(fabsf(dWs2), kernel_root);

      const float dvals[7] = {fabsf(vW.f[j] - Ws),   fabsf(vdW.f[j] - dWs),
                              fabsf(vW2.f[j] - Ws),  fabsf(vdW2.f[j] - dWs),
                              fabsf(vWb.f[j] - Ws),  fabsf(vdWc.f[j] - dWs),
                              fabsf(vdWd.f[j] - dWs)};
      for (int k = 0; k < 7; ++k) max_diff_vec = fmaxf(max_diff_vec, dvals[k]);
      max_diff_vec = fmaxf(max_diff_vec, fabsf(vdWe.f[j] - dWs));
      max_diff_vec = fmaxf(max_diff_vec, fabsf(vdWb.f[j] - dWs2));

      if (dvals[0] > tol_W || dvals[2] > tol_W || dvals[4] > tol_W)
        error(
            "Vector W differs from scalar: u=%.9g scalar=%.9g vector=%.9g "
            "%.9g %.9g",
            vu.f[j], Ws, vW.f[j], vW2.f[j], vWb.f[j]);
      if (dvals[1] > tol_dW || dvals[3] > tol_dW || dvals[5] > tol_dW ||
          dvals[6] > tol_dW || fabsf(vdWe.f[j] - dWs) > tol_dW)
        error(
            "Vector dW differs from scalar: u=%.9g scalar=%.9g vector=%.9g "
            "%.9g %.9g %.9g %.9g",
            vu.f[j], dWs, vdW.f[j], vdW2.f[j], vdWc.f[j], vdWd.f[j], vdWe.f[j]);
      if (fabsf(vdWb.f[j] - dWs2) > tol_dW2)
        error(
            "Vector dW (2nd vector) differs from scalar: u=%.9g scalar=%.9g "
            "vector=%.9g",
            vu2.f[j], dWs2, vdWb.f[j]);

      /* Exactly zero outside the support */
      if (vu.f[j] >= kernel_gamma &&
          (vW.f[j] != 0.f || vdW.f[j] != 0.f || vWb.f[j] != 0.f ||
           vdWc.f[j] != 0.f || vdWd.f[j] != 0.f || vdWe.f[j] != 0.f))
        error("Vector kernel not zero outside the support: u=%.9g", vu.f[j]);
    }
  }

  message("Vector (VEC_SIZE=%d): max |vector - scalar| = %.2e", VEC_SIZE,
          max_diff_vec);
#endif

  free(u);
  message("All kernel accuracy checks passed");
  return 0;
}
