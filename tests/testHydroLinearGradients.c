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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local headers. */
#include "swift.h"

/**
 * @file testHydroLinearGradients.c
 * @brief Checks that the matrix-inversion gradients of the MAGMA scheme are
 * exact for linear fields.
 *
 * Particles on a lattice and on a jittered lattice carry a linear velocity
 * field v = A x + b and a linear internal energy u = c . x + d. After the
 * density and gradient loops (brute force, every pair within the kernel
 * support) and the gradient finishing, the matrix gradients of v and u must
 * reproduce A and c, the divergence the trace of A, and the SPH Laplacian of u
 * with its linear correction must vanish (from the second gradient pass on,
 * since the correction uses the gradient of the previous pass). With SPHENIX
 * the same set-up reports the errors of its kernel-gradient divergence and
 * uncorrected Laplacian for comparison, without an exactness requirement
 * (they are exact on the lattice by symmetry only). The fields are not
 * periodic, so only
 * the particles whose kernel support stays inside the box are checked. The
 * test runs in the configured dimension and does nothing for the other hydro
 * schemes.
 */

#if !defined(MAGMA_SPH) && !defined(SPHENIX_SPH)

int main(int argc, char *argv[]) {
  message("Neither a MAGMA nor a SPHENIX build: nothing to test.");
  return 0;
}

#else

/* A deterministic pseudo-random number in [-1, 1] */
static double rand_sym(void) { return 2.0 * rand() / (double)RAND_MAX - 1.0; }

/**
 * @brief Runs one configuration and returns the number of failed checks.
 *
 * @param n Particles per dimension.
 * @param jitter Random displacement of the lattice, in units of the spacing.
 * @param eta Resolution parameter (h in units of the spacing).
 * @param hp The hydro properties.
 * @param cosmo The (non-cosmological) cosmology.
 * @param floor The pressure floor properties.
 * @param label Name of the configuration for the messages.
 */
static int run_configuration(const int n, const double jitter, const double eta,
                             const struct hydro_props *hp,
                             const struct cosmology *cosmo,
                             const struct pressure_floor_props *floor,
                             const char *label) {

  const double L = 1.0;
  const double dx = L / n;
  const double rho0 = 1.0;
  const int N = (int)pow(n, hydro_dimension);
  const float h = eta * dx;
  const float H = kernel_gamma * h;
  const double mass = rho0 * pow(dx, hydro_dimension);

  /* The linear fields: v = A x + b, u = c . x + d (u stays positive) */
  const double A[3][3] = {
      {0.7, -0.3, 0.2}, {0.4, 1.1, -0.5}, {-0.6, 0.25, 0.9}};
  const double b[3] = {0.1, -0.2, 0.05};
  const double c[3] = {0.8, -0.35, 0.5};
  const double d = 3.0;

  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  if (posix_memalign((void **)&parts, part_align, N * sizeof(struct part)) != 0)
    error("Failed to allocate the particles");
  if (posix_memalign((void **)&xparts, xpart_align, N * sizeof(struct xpart)) !=
      0)
    error("Failed to allocate the extended particles");
  bzero(parts, N * sizeof(struct part));
  bzero(xparts, N * sizeof(struct xpart));

  /* Place the particles and set the fields */
  srand(42);
  for (int i = 0; i < N; ++i) {
    struct part *p = &parts[i];
    int idx = i;
    double x[3] = {0.5 * dx, 0.5 * dx, 0.5 * dx};
    for (int k = 0; k < (int)hydro_dimension; ++k) {
      x[k] = (idx % n + 0.5) * dx + jitter * dx * rand_sym();
      idx /= n;
    }
    for (int k = (int)hydro_dimension; k < 3; ++k) x[k] = 0.;
    for (int k = 0; k < 3; ++k) p->x[k] = x[k];
    for (int k = 0; k < 3; ++k) {
      double v = b[k];
      for (int l = 0; l < (int)hydro_dimension; ++l) v += A[k][l] * x[l];
      p->v[k] = (k < (int)hydro_dimension) ? v : 0.f;
    }
    double u = d;
    for (int l = 0; l < (int)hydro_dimension; ++l) u += c[l] * x[l];
    p->u = u;
    p->mass = mass;
    p->h = h;
    p->id = i + 1;
    p->time_bin = 1;
    hydro_first_init_part(p, &xparts[i]);
    p->time_bin = 1;
    hydro_init_part(p, NULL);
  }

  const float hig2 = H * H;

  /* Density loop: every ordered pair within the kernel support of i */
  for (int i = 0; i < N; ++i) {
    struct part *pi = &parts[i];
    for (int j = 0; j < N; ++j) {
      if (i == j) continue;
      struct part *pj = &parts[j];
      const float dxv[3] = {(float)(pi->x[0] - pj->x[0]),
                            (float)(pi->x[1] - pj->x[1]),
                            (float)(pi->x[2] - pj->x[2])};
      const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
      if (r2 < hig2)
        runner_iact_nonsym_density(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
    }
  }
  for (int i = 0; i < N; ++i) {
    hydro_end_density(&parts[i], cosmo);
    hydro_prepare_gradient(&parts[i], &xparts[i], cosmo, hp, floor);
    hydro_reset_gradient(&parts[i]);
  }

  /* Gradient loop and finishing. With MAGMA twice: its SPH Laplacian of u
   * removes the linear part of the field with the gradient of the previous
   * evaluation (hydro_prepare_force() stores it), so it is only exact from
   * the second pass on; the matrix gradients are exact from the first.
   * SPHENIX resets its Laplacian in the density phase only, so one pass. */
#if defined(MAGMA_SPH)
  const int n_passes = 2;
#else
  const int n_passes = 1;
#endif
  for (int pass = 0; pass < n_passes; ++pass) {
    if (pass > 0)
      for (int i = 0; i < N; ++i) hydro_reset_gradient(&parts[i]);
    for (int i = 0; i < N; ++i) {
      struct part *pi = &parts[i];
      for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        struct part *pj = &parts[j];
        const float dxv[3] = {(float)(pi->x[0] - pj->x[0]),
                              (float)(pi->x[1] - pj->x[1]),
                              (float)(pi->x[2] - pj->x[2])};
        const float r2 = dxv[0] * dxv[0] + dxv[1] * dxv[1] + dxv[2] * dxv[2];
        if (r2 < hig2)
          runner_iact_nonsym_gradient(r2, dxv, pi->h, pj->h, pi, pj, 1.f, 0.f);
      }
    }
    for (int i = 0; i < N; ++i) {
      hydro_end_gradient(&parts[i]);
      hydro_prepare_force(&parts[i], &xparts[i], cosmo, hp, floor,
                          /*dt_alpha=*/0.f, /*dt_therm=*/0.f);
    }
  }

  /* Check the interior particles */
  const double margin = 2.0 * H;
#if defined(MAGMA_SPH)
  const double tol_grad = 1e-4; /* relative to the gradient norm */
  const double tol_lap = 1e-3;  /* |laplace u| h / |c| */
#else
  const double tol_grad = 1.; /* a kernel estimate: sanity only */
  const double tol_lap = 1.;
#endif
  double norm_A = 0., norm_c = 0.;
  for (int k = 0; k < (int)hydro_dimension; ++k) {
    norm_c += c[k] * c[k];
    for (int l = 0; l < (int)hydro_dimension; ++l) norm_A += A[k][l] * A[k][l];
  }
  norm_A = sqrt(norm_A);
  norm_c = sqrt(norm_c);
  double trace_A = 0.;
  for (int k = 0; k < (int)hydro_dimension; ++k) trace_A += A[k][k];

  int checked = 0, failed = 0, fallbacks = 0;
  double worst_v = 0., worst_u = 0., worst_div = 0., worst_lap = 0.;
  for (int i = 0; i < N; ++i) {
    const struct part *p = &parts[i];
    int interior = 1;
    for (int k = 0; k < (int)hydro_dimension; ++k)
      if (p->x[k] < margin || p->x[k] > L - margin) interior = 0;
    if (!interior) continue;
    ++checked;
#if defined(MAGMA_SPH)
    if (p->fallback_flags & magma_fallback_base_sph_mask) {
      ++fallbacks;
      continue;
    }

    /* Velocity gradient: row k of A is the gradient of v_k */
    const float *grad_v[3] = {p->force.gradient_vx, p->force.gradient_vy,
                              p->force.gradient_vz};
    double err_v = 0.;
    for (int k = 0; k < (int)hydro_dimension; ++k)
      for (int l = 0; l < (int)hydro_dimension; ++l)
        err_v = fmax(err_v, fabs(grad_v[k][l] - A[k][l]));
    err_v /= norm_A;

    double err_u = 0.;
    for (int l = 0; l < (int)hydro_dimension; ++l)
      err_u = fmax(err_u, fabs(p->force.gradient_u[l] - c[l]));
    err_u /= norm_c;

    const double err_div = fabs(p->div_v - trace_A) / norm_A;
    const double err_lap = fabs(p->laplace_u) * h / norm_c;
#else
    /* SPHENIX: no gradient tensor nor gradient of u; the kernel-gradient
     * divergence and the uncorrected Laplacian are reported, not required
     * to be exact. */
    const double err_v = 0., err_u = 0.;
    const double err_div = fabs(p->viscosity.div_v - trace_A) / norm_A;
    const double err_lap = fabs(p->diffusion.laplace_u) * h / norm_c;
#endif

    worst_v = fmax(worst_v, err_v);
    worst_u = fmax(worst_u, err_u);
    worst_div = fmax(worst_div, err_div);
    worst_lap = fmax(worst_lap, err_lap);

    if (err_v > tol_grad || err_u > tol_grad || err_div > tol_grad ||
        err_lap > tol_lap) {
      if (failed < 5)
        message(
            "%s: particle %lld at (%.3f, %.3f, %.3f): grad v err %.2e, "
            "grad u err %.2e, div v err %.2e, laplace u err %.2e",
            label, p->id, p->x[0], p->x[1], p->x[2], err_v, err_u, err_div,
            err_lap);
      ++failed;
    }
  }

  message(
      "%s: %d particles, %d checked, %d fallbacks, %d failed; worst errors: "
      "grad v %.2e, grad u %.2e, div v %.2e, laplace u %.2e",
      label, N, checked, fallbacks, failed, worst_v, worst_u, worst_div,
      worst_lap);

  if (checked == 0) error("%s: no interior particle to check", label);
  if (fallbacks > 0) {
    message("%s: %d interior particles fell back to standard SPH", label,
            fallbacks);
    failed += fallbacks;
  }

  free(parts);
  free(xparts);
  return failed;
}

int main(int argc, char *argv[]) {

  /* Initialise the hydro properties, cosmology and floors */
  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);
  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);
  struct pressure_floor_props floor;
  bzero(&floor, sizeof(struct pressure_floor_props));

  /* Particles per dimension: enough interior particles at eta 1.5 */
#if defined(HYDRO_DIMENSION_1D)
  const int n = 512;
#elif defined(HYDRO_DIMENSION_2D)
  const int n = 64;
#else
  const int n = 24;
#endif
  const double eta = 1.5;

  int failed = 0;
  failed += run_configuration(n, 0.0, eta, &hp, &cosmo, &floor, "lattice");
  failed +=
      run_configuration(n, 0.3, eta, &hp, &cosmo, &floor, "jittered lattice");

  if (failed > 0) {
    message("%d checks failed", failed);
    return 1;
  }
  message("All checks passed.");
  return 0;
}

#endif /* MAGMA_SPH */
