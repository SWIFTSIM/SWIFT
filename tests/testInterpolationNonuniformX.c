/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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
#include "feedback/GEAR/interpolation.h"
#include "swift.h"

#include <config.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NZ 13
#define NM 8

/* The shipped PopII Data/SW/MetallicityDependent metallicity axis. Its
   |z0| / dz of 134 is what makes recovering dz by differencing two adjacent
   float32 nodes lose about 8e-06 of relative precision. */
#define NZ_UNIFORM 110
static const float uniform_z0 = -9.847711655616944f;
static const float uniform_dz = 0.07339449541284404f;

/* Relative bar, with an absolute floor for values near zero. Two orders of
   magnitude below the 2.6e-05 a step recovered by differencing costs. */
#define TOL_REL 1e-6

/* The shipped PopII Data/Radiation metallicity rungs: strictly increasing,
   and spaced so no uniform log grid reproduces them. */
static const float z_rungs[NZ] = {1e-11f,  1e-6f, 1e-4f, 1e-3f, 2e-3f,
                                  4e-3f,   6e-3f, 8e-3f, 1e-2f, 1.4e-2f,
                                  1.7e-2f, 2e-2f, 3e-2f};

/* Linear in log_z, not in the rung index, so an interpolation that got the
   node spacing wrong lands on a different value. */
static double rung_value(double log_z, int j) {
  return 1.0 + 3.0 * log_z + 0.25 * j;
}

/* Linear in the x index, so the uniform leg's expected value is the same
   closed form the plain (log_x - xmin) / dx arithmetic produces. */
static double uniform_value(double x_index, int j) {
  return 1.0 + x_index + 0.25 * j;
}

static void check_finite(double v, const char *what) {
  if (!(v > -DBL_MAX && v < DBL_MAX)) error("%s is not finite (%g)", what, v);
}

static void check_close(double got, double want, const char *what) {
  check_finite(got, what);
  const double bar = TOL_REL * (fabs(want) > 1. ? fabs(want) : 1.);
  if (fabs(got - want) > bar)
    error("%s: got %.9g, expected %.9g (difference %.3e, bar %.3e)", what, got,
          want, fabs(got - want), bar);
}

static double *allocate_data(int nx, int ny, const char *what) {
  double *d = (double *)malloc(sizeof(double) * nx * ny);
  if (d == NULL) error("Failed to allocate the %s test data", what);
  return d;
}

int main(int argc, char *argv[]) {

  const float log_m_min = 0.f, log_m_max = 1.f;
  const float dm = (log_m_max - log_m_min) / (NM - 1.f);

  /**************************************************************************/
  /* Non-uniform axis: the table's own metallicity rungs.                   */
  /**************************************************************************/

  float log_z[NZ];
  for (int i = 0; i < NZ; i++) log_z[i] = log10f(z_rungs[i]);

  double *data = allocate_data(NZ, NM, "non-uniform");
  for (int i = 0; i < NZ; i++)
    for (int j = 0; j < NM; j++) data[i * NM + j] = rung_value(log_z[i], j);

  /* Node lists identical on input and output: no resample. */
  struct interpolation_2d interp;
  interpolate_2d_init(&interp, log_z, NZ, log_z, NZ, log_m_min, log_m_max, NM,
                      log_m_min, dm, NM, data, boundary_condition_const,
                      boundary_condition_const);

  /* A query on a rung must return that rung's own row. */
  for (int i = 0; i < NZ; i++) {
    for (int j = 0; j < NM; j++) {
      const float log_m = log_m_min + j * dm;
      check_close(interpolate_2d(&interp, log_z[i], log_m),
                  rung_value(log_z[i], j), "rung query");
    }
  }

  /* Between two rungs the result must follow log_z, not the rung index: a
     quarter of the way along in log_z is a quarter of the way along in
     value only if the interval's own width was used. */
  const float fractions[3] = {0.25f, 0.5f, 0.75f};
  for (int i = 0; i < NZ - 1; i++) {
    for (int k = 0; k < 3; k++) {
      const float log_z_q = log_z[i] + fractions[k] * (log_z[i + 1] - log_z[i]);
      check_close(interpolate_2d(&interp, log_z_q, log_m_min),
                  rung_value(log_z_q, 0), "between-rung query");
    }
  }

  /* Below the first rung. Pristine gas floors to 1e-300, so this is the
     common production path, and the constant boundary returns row 0. */
  for (int j = 0; j < NM; j++)
    check_close(interpolate_2d(&interp, -300.f, log_m_min + j * dm),
                rung_value(log_z[0], j), "below-first-rung query");

  /* Above the last rung, same policy, last row. */
  for (int j = 0; j < NM; j++)
    check_close(
        interpolate_2d(&interp, log_z[NZ - 1] + 1.f, log_m_min + j * dm),
        rung_value(log_z[NZ - 1], j), "above-last-rung query");

  interpolate_2d_free(&interp);
  free(data);

  /**************************************************************************/
  /* Uniform axis: the stellar-wind metallicity axis, which must keep       */
  /* dividing by its own step rather than one recovered from the nodes.     */
  /**************************************************************************/

  const float ux_min = uniform_z0;
  const float ux_max = ux_min + uniform_dz * (NZ_UNIFORM - 1.f);
  const float udx = (ux_max - ux_min) / (NZ_UNIFORM - 1.f);

  double *udata = allocate_data(NZ_UNIFORM, NM, "uniform");
  for (int i = 0; i < NZ_UNIFORM; i++)
    for (int j = 0; j < NM; j++) udata[i * NM + j] = uniform_value(i, j);

  struct interpolation_2d uni;
  interpolate_2d_init_uniform_x(
      &uni, ux_min, ux_max, NZ_UNIFORM, log_m_min, log_m_max, NM, ux_min,
      log_m_min, uniform_dz, dm, NZ_UNIFORM, NM, udata,
      boundary_condition_const, boundary_condition_const);

  /* Each stored row sits at the source index the build-time division puts
     it at, and each query blends two stored rows at the fraction the
     query-time division puts it at. Both divisions here are written out in
     full rather than read back from the table, so the check fails if the
     table ever recovers a step from its node list instead. */
  for (int k = 0; k <= 400; k++) {
    const float log_x = ux_min + k * (ux_max - ux_min) / 400.f;
    const float idxf = (log_x - ux_min) / udx;

    int i0 = (int)idxf;
    if (i0 > NZ_UNIFORM - 2) i0 = NZ_UNIFORM - 2;
    const double f = idxf - i0;

    const float node_lo = ux_min + i0 * udx;
    const float node_hi = ux_min + (i0 + 1) * udx;
    const double row_lo = uniform_value((node_lo - ux_min) / uniform_dz, 0);
    const double row_hi = uniform_value((node_hi - ux_min) / uniform_dz, 0);
    const double want = row_lo * (1. - f) + row_hi * f;

    check_close(interpolate_2d(&uni, log_x, log_m_min), want, "uniform query");
  }

  interpolate_2d_free(&uni);
  free(udata);

  message("All non-uniform x interpolation tests passed.");
  return 0;
}
