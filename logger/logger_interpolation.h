/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#ifndef LOGGER_LOGGER_INTERPOLATION_H
#define LOGGER_LOGGER_INTERPOLATION_H

#include "logger_tools.h"

/**
 * @brief Compute the quintic hermite spline interpolation.
 * See https://en.wikipedia.org/wiki/Hermite_interpolation for more information.
 * @f$ p(x) = f(x_0) + f'(x_0) (x - x_0) + \frac{1}{2}f''(x_0) (x - x_0)^2 +
 \frac{f(x_1) - f(x_0) - f'(x_0) (x_1 - x_0) - \frac{1}{2} f''(x_0) (x_1 -
 x_0)^2}{(x_1 - x_0)^3} (x - x_0)^3
  + \frac{3 f(x_0) - 3 f(x_1) + \left( 2 f'(x_0) + f'(x_1) \right) (x_1 - x_0) +
 \frac{1}{2} f''(x_0) (x_1 - x_0)^2}{(x_1 - x_0)^4} (x - x_0)^3 (x - x_1)
  + \frac{6 f(x_1) - 6 f(x_0) - 3 \left( f'(x_0) + f'(x_1) \right) (x_1 - x_0) +
 \frac{1}{2}\left( f''(x_1) - f''(x_0) \right) (x_1 - x_0)^2}{(x_1 - x_0)^5} (x
 - x_0)^3 (x - x_1)^2 @f$
 *
 * @param t0 The time at the left of the interval.
 * @param x0 The function at the left of the interval.
 * @param v0 The first derivative at the left of the interval.
 * @param a0 The second derivative at the left of the interval.
 * @param t1 The time at the right of the interval.
 * @param x1 The function at the right of the interval.
 * @param v1 The first derivative at the right of the interval.
 * @param a1 The second derivative at the right of the interval.
 * @param t The time of the interpolation.
 *
 * @return The function evaluated at t.
 */
__attribute__((always_inline)) INLINE static double
interpolate_quintic_hermite_spline(const double t0, const double x0,
                                   const float v0, const float a0,
                                   const double t1, const double x1,
                                   const float v1, const float a1,
                                   const double t) {

  /* Generates recurring variables  */
  /* Time differences */
  const double dt = t1 - t0;
  const double dt2 = dt * dt;
  const double dt3 = dt2 * dt;
  const double dt4 = dt3 * dt;
  const double dt5 = dt4 * dt;

  const double t_t0 = t - t0;
  const double t_t0_2 = t_t0 * t_t0;
  const double t_t0_3 = t_t0_2 * t_t0;
  const double t_t1 = t - t1;
  const double t_t1_2 = t_t1 * t_t1;

  /* Derivatives */
  const double v0_dt = v0 * dt;
  const double a0_dt2 = 0.5 * a0 * dt2;
  const double v1_dt = v1 * dt;
  const double a1_dt2 = 0.5 * a1 * dt2;

  /* Do the first 3 terms of the hermite spline */
  double x = x0 + v0 * t_t0 + 0.5 * a0 * t_t0_2;

  /* Cubic term */
  x += (x1 - x0 - v0_dt - a0_dt2) * t_t0_3 / dt3;

  /* Quartic term */
  x += (3. * x0 - 3. * x1 + v1_dt + 2. * v0_dt + a0_dt2) * t_t0_3 * t_t1 / dt4;

  /* Quintic term */
  x += (6. * x1 - 6. * x0 - 3. * v0_dt - 3. * v1_dt + a1_dt2 - a0_dt2) *
       t_t0_3 * t_t1_2 / dt5;

  return x;
}

/**
 * @brief Compute the cubic hermite spline interpolation.
 * See https://en.wikipedia.org/wiki/Hermite_interpolation for more information.
 * @f$ p(x) = f(x_0) + f'(x_0) (x - x_0) + \frac{f(x_1) - f(x_0) - f'(x_0) (x_1
 - x_0)}{(x_1 - x_0)^2} (x - x_0)^2
  + \frac{2 f(x_0) - 2 f(x_1) + f'(x_0) (x_1 - x_0) + f'(x_1) (x_1 - x_0)}{(x_1
 - x_0)^3} (x - x_0)^2 (x - x_1) @f$
 *
 * @param t0 The time at the left of the interval.
 * @param v0 The first derivative at the left of the interval.
 * @param a0 The second derivative at the left of the interval.
 * @param t1 The time at the right of the interval.
 * @param v1 The first derivative at the right of the interval.
 * @param a1 The second derivative at the right of the interval.
 * @param t The time of the interpolation.
 *
 * @return The function evaluated at t.
 */
__attribute__((always_inline)) INLINE static float
interpolate_cubic_hermite_spline(const double t0, const float v0,
                                 const float a0, const double t1,
                                 const float v1, const float a1,
                                 const double t) {

  /* Generates recurring variables  */
  /* Time differences */
  const float dt = t1 - t0;
  const float dt2 = dt * dt;
  const float dt3 = dt2 * dt;

  const float t_t0 = t - t0;
  const float t_t0_2 = t_t0 * t_t0;
  const float t_t1 = t - t1;

  /* Derivatives */
  const float a0_dt = a0 * dt;
  const float a1_dt = a1 * dt;

  /* Do the first 2 terms of the hermite spline */
  float x = v0 + a0 * t_t0;

  /* Square term */
  x += (v1 - v0 - a0_dt) * t_t0_2 / dt2;

  /* Cubic term */
  x += (2. * v0 - 2. * v1 + a1_dt + a0_dt) * t_t0_2 * t_t1 / dt3;

  return x;
}

/**
 * @brief Interpolates a field in N dimension using the first and second
 * derivatives if possible.
 *
 * The field is in double and both derivatives in float.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 */
__attribute__((always_inline)) INLINE static void
interpolate_quintic_double_float_ND(const double t_before,
                                    const struct logger_field *restrict before,
                                    const double t_after,
                                    const struct logger_field *restrict after,
                                    void *restrict output, const double t,
                                    const int dimension) {
  /* Compute the interpolation scaling. */
  const double wa = (t - t_before) / (t_after - t_before);
  const double wb = 1. - wa;

  for (int i = 0; i < dimension; i++) {
    double *x = (double *)output;
    const double *x_bef = (double *)before->field;
    const float *v_bef = (float *)before->first_deriv;
    const float *a_bef = (float *)before->second_deriv;

    const double *x_aft = (double *)after->field;
    const float *v_aft = (float *)after->first_deriv;
    const float *a_aft = (float *)after->second_deriv;

    /* Use quintic hermite spline. */
    if (v_bef && v_aft && a_bef && a_aft) {
      x[i] = interpolate_quintic_hermite_spline(t_before, x_bef[i], v_bef[i],
                                                a_bef[i], t_after, x_aft[i],
                                                v_aft[i], a_aft[i], t);
    }
    /* Use cubic hermite spline. */
    else if (v_bef && v_aft) {
      x[i] = interpolate_cubic_hermite_spline(t_before, x_bef[i], v_bef[i],
                                              t_after, x_aft[i], v_aft[i], t);
    }
    /* Use linear interpolation. */
    else {
      x[i] = wa * x_aft[i] + wb * x_bef[i];
    }
  }
}

/**
 * @brief Interpolates a field in N dimension using the first derivative if
 * possible.
 *
 * The field and the first derivatives are floats.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 */
__attribute__((always_inline)) INLINE static void interpolate_cubic_float_ND(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int dimension) {

  /* Compute the interpolation scaling. */
  const float wa = (t - t_before) / (t_after - t_before);
  const float wb = 1. - wa;

  for (int i = 0; i < dimension; i++) {
    float *v = (float *)output;
    const float *v_bef = (float *)before->field;
    const float *a_bef = (float *)before->first_deriv;

    const float *v_aft = (float *)after->field;
    const float *a_aft = (float *)after->first_deriv;

    /* Use a cubic hermite spline. */
    if (a_bef && a_aft) {
      v[i] = interpolate_cubic_hermite_spline(t_before, v_bef[i], a_bef[i],
                                              t_after, v_aft[i], a_aft[i], t);
    }
    /* Use linear interpolation. */
    else {
      v[i] = wa * v_aft[i] + wb * v_bef[i];
    }
  }
}

/**
 * @brief Interpolates a field in N dimension.
 *
 * The field is in float.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 */
__attribute__((always_inline)) INLINE static void interpolate_linear_float_ND(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t, const int dimension) {

  /* Compute the interpolation scaling. */
  const float wa = (t - t_before) / (t_after - t_before);
  const float wb = 1. - wa;

  /* interpolate vectors. */
  for (int i = 0; i < dimension; i++) {
    float *a = (float *)output;
    const float *a_bef = (float *)before->field;
    const float *a_aft = (float *)after->field;
    a[i] = wa * a_aft[i] + wb * a_bef[i];
  }
}

/**
 * @brief Interpolates a field stored as a float.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 */
__attribute__((always_inline)) INLINE static void interpolate_linear_float(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t) {

  /* Compute the interpolation scaling. */
  const float wa = (t - t_before) / (t_after - t_before);
  const float wb = 1. - wa;
  ((float *)output)[0] =
      wa * ((float *)after->field)[0] + wb * ((float *)before->field)[0];
}

/**
 * @brief Interpolation function for the ids.
 *
 * As the IDs should not change during the simulation, we ensure
 * that the IDs are the same between two records.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param t_before Time of field_before (< t).
 * @param t_after Time of field_after (> t).
 * @param t Requested time.
 */
__attribute__((always_inline)) INLINE static void interpolate_ids(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t) {
  if (*(long long *)after->field != *(long long *)before->field) {
    error("Interpolating different particles");
  }
  *(long long *)output = *(long long *)after->field;
}
#endif  // LOGGER_LOGGER_INTERPOLATION_H
