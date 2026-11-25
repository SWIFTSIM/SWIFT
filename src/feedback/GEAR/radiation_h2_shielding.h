/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_H2_SHIELDING_GEAR_H
#define SWIFT_RADIATION_H2_SHIELDING_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_h2_shielding.h
 * @brief Self-shielding of the H2 lines by a column, its inverse and the hop
 * of the line transmission of a packet share.
 *
 * The factor is the fit of Wolcott-Green & Haiman (2019) in the form the
 * geajay library uses, with a floor on the density. All quantities are double
 * and in cgs, except the velocities, which share one unit.
 */

#include "inline.h"

#include <math.h>

/*! Lowest gas density of the fit, in cm^-3. Below it the exponent of the fit
    turns negative and the factor rises with the column. */
#define RADIATION_H2_SHIELDING_DENSITY_FLOOR_CM3 1e-8

/*! Highest gas density of the fit, in cm^-3. */
#define RADIATION_H2_SHIELDING_DENSITY_CAP_CM3 1e7

/*! Temperature range of the fit, in K. */
#define RADIATION_H2_SHIELDING_TEMPERATURE_FLOOR_K 100.
#define RADIATION_H2_SHIELDING_TEMPERATURE_CAP_K 8000.

/*! Column of the unit of the dimensionless column `x` of the fit, in cm^-2. */
#define RADIATION_H2_SHIELDING_COLUMN_UNIT_CM2 5e14

/*! Largest column that a share carries, in cm^-2. A transmission that is not
    positive maps to it. */
#define RADIATION_H2_SHIELDING_COLUMN_CAP_CM2 1e25

/*! Smallest column that the inverse resolves, in cm^-2. A transmission above
    the factor at this column maps to 0. */
#define RADIATION_H2_SHIELDING_COLUMN_FLOOR_CM2 1e2

/*! Relative velocity of the share and the gas, in thermal line widths, beyond
    which the column no longer shields. */
#define RADIATION_H2_SHIELDING_VELOCITY_RESET_WIDTHS 3.

/*! Tolerance of the inverse in the natural logarithm of the column: 2^-26, a
    quarter of the float rounding of a stored column. */
#define RADIATION_H2_SHIELDING_INVERSE_TOLERANCE 1.4901161193847656e-8

/*! Largest number of iterations of the inverse. Bisection alone needs 32 to
    reach the tolerance from the full range of columns. */
#define RADIATION_H2_SHIELDING_INVERSE_MAX_ITERATIONS 64

/**
 * @brief Exponent of the fit, a function of the gas temperature and density.
 *
 * Both arguments are clamped to the range of the fit, so any real value,
 * zero or negative included, gives a finite exponent, about 0.139 at the
 * lowest.
 *
 * @param T Gas temperature in K.
 * @param n Gas particle density in cm^-3.
 * @return The exponent of the fit.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_exponent(const double T, const double n) {
  double Tc = T;
  if (Tc < RADIATION_H2_SHIELDING_TEMPERATURE_FLOOR_K)
    Tc = RADIATION_H2_SHIELDING_TEMPERATURE_FLOOR_K;
  if (Tc > RADIATION_H2_SHIELDING_TEMPERATURE_CAP_K)
    Tc = RADIATION_H2_SHIELDING_TEMPERATURE_CAP_K;
  double nc = n;
  if (nc < RADIATION_H2_SHIELDING_DENSITY_FLOOR_CM3)
    nc = RADIATION_H2_SHIELDING_DENSITY_FLOOR_CM3;
  if (nc > RADIATION_H2_SHIELDING_DENSITY_CAP_CM3)
    nc = RADIATION_H2_SHIELDING_DENSITY_CAP_CM3;
  const double lgT = log10(Tc);
  return (0.8711 * lgT - 1.928) * exp(-0.2856 * log10(nc)) +
         (-0.9639 * lgT + 3.892);
}

/**
 * @brief The two terms of the self-shielding factor of an H2 column
 * (#radiation_h2_shielding_factor): the power term of the line cores and the
 * damped term of the line wings.
 *
 * The power is evaluated as an exponential of a logarithm, so that a large
 * exponent and column cannot overflow. The column is clamped to
 * [0, #RADIATION_H2_SHIELDING_COLUMN_CAP_CM2].
 *
 * @param N H2 column in cm^-2.
 * @param alpha Exponent of the fit (#radiation_h2_shielding_exponent).
 * @param b5 Thermal line width of H2 in km/s. Positive.
 * @param x (return) The clamped column in units of
 * #RADIATION_H2_SHIELDING_COLUMN_UNIT_CM2.
 * @param power (return) The power term.
 * @param damped (return) The damped term.
 */
__attribute__((always_inline)) INLINE static void
radiation_h2_shielding_factor_terms(const double N, const double alpha,
                                    const double b5, double *x, double *power,
                                    double *damped) {
  double Nc = N;
  if (!(Nc > 0.)) Nc = 0.;
  if (Nc > RADIATION_H2_SHIELDING_COLUMN_CAP_CM2)
    Nc = RADIATION_H2_SHIELDING_COLUMN_CAP_CM2;
  *x = Nc / RADIATION_H2_SHIELDING_COLUMN_UNIT_CM2;
  const double root = sqrt(1. + *x);
  *power = 0.965 * exp(-alpha * log1p(*x / b5));
  *damped = 0.035 * exp(-8.5e-4 * root) / root;
}

/**
 * @brief Self-shielding factor of an H2 column: the dissociation rate at the
 * column over the rate at column 0.
 *
 * @param N H2 column in cm^-2, clamped as in
 * #radiation_h2_shielding_factor_terms.
 * @param alpha Exponent of the fit (#radiation_h2_shielding_exponent).
 * @param b5 Thermal line width of H2 in km/s. Positive.
 * @return The factor, at most 1.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_factor(const double N, const double alpha,
                              const double b5) {
  double x, power, damped;
  radiation_h2_shielding_factor_terms(N, alpha, b5, &x, &power, &damped);
  const double f = power + damped;
  return f < 1. ? f : 1.;
}

/**
 * @brief Logarithm of the self-shielding factor of an H2 column, and its slope
 * in the logarithm of the column.
 *
 * @param N H2 column in cm^-2, positive and at most
 * #RADIATION_H2_SHIELDING_COLUMN_CAP_CM2.
 * @param alpha Positive exponent of the fit (#radiation_h2_shielding_exponent).
 * @param b5 Thermal line width of H2 in km/s. Positive.
 * @param slope (return) d ln S / d ln N, negative.
 * @return ln S, without the cap of the factor at 1.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_log_factor_slope(const double N, const double alpha,
                                        const double b5, double *slope) {
  double x, power, damped;
  radiation_h2_shielding_factor_terms(N, alpha, b5, &x, &power, &damped);
  const double f = power + damped;
  const double wing_slope = 8.5e-4 * sqrt(1. + x) + 1.;
  *slope = -(power * alpha * x / (b5 + x) +
             0.5 * damped * wing_slope * x / (1. + x)) /
           f;
  return log(f);
}

/**
 * @brief Self-shielding factor at column 0, the same for every fit
 * (#radiation_h2_shielding_factor).
 *
 * @return The factor at column 0.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_factor_at_zero(void) {
  return 0.965 + 0.035 * exp(-8.5e-4);
}

/**
 * @brief Column whose shielding factor is a given transmission, searched from
 * a first guess.
 *
 * The factor falls with the column, so the root is bracketed in the logarithm
 * u of the column and found by Newton steps on ln S(u) - ln t. A step that
 * leaves the bracket, or that does not halve the step before it, is replaced
 * by a bisection of the bracket. The search ends when a step is below
 * #RADIATION_H2_SHIELDING_INVERSE_TOLERANCE, where the next error of a Newton
 * step is of the order of its square. A transmission at or above the factor
 * at #RADIATION_H2_SHIELDING_COLUMN_FLOOR_CM2 gives 0, one at or below the
 * factor at #RADIATION_H2_SHIELDING_COLUMN_CAP_CM2 (a zero or negative one
 * included) gives that cap.
 *
 * @param t Transmission.
 * @param alpha Positive exponent of the fit (#radiation_h2_shielding_exponent).
 * @param b5 Thermal line width of H2 in km/s. Positive.
 * @param N_guess First guess of the column in cm^-2. A guess outside the
 * bracket, 0 included, is replaced by the column of the power term alone.
 * @return The column in cm^-2.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_inverse_from(const double t, const double alpha,
                                    const double b5, const double N_guess) {
  if (!(t > 0.)) return RADIATION_H2_SHIELDING_COLUMN_CAP_CM2;
  if (t >= radiation_h2_shielding_factor(
               RADIATION_H2_SHIELDING_COLUMN_FLOOR_CM2, alpha, b5))
    return 0.;
  if (t <= radiation_h2_shielding_factor(RADIATION_H2_SHIELDING_COLUMN_CAP_CM2,
                                         alpha, b5))
    return RADIATION_H2_SHIELDING_COLUMN_CAP_CM2;

  double lo = log(RADIATION_H2_SHIELDING_COLUMN_FLOOR_CM2);
  double hi = log(RADIATION_H2_SHIELDING_COLUMN_CAP_CM2);
  const double log_t = log(t);

  double u = 0.5 * (lo + hi);
  if (N_guess > RADIATION_H2_SHIELDING_COLUMN_FLOOR_CM2 &&
      N_guess < RADIATION_H2_SHIELDING_COLUMN_CAP_CM2) {
    u = log(N_guess);
  } else if (t < 0.965) {
    /* The column at which the power term alone is t. It is below the cap,
       since t exceeds the factor at the cap. */
    const double x = b5 * expm1(-log(t / 0.965) / alpha);
    if (x > 0.) {
      const double u_power = log(x * RADIATION_H2_SHIELDING_COLUMN_UNIT_CM2);
      if (u_power > lo && u_power < hi) u = u_power;
    }
  }

  double step_before = hi - lo;
  for (int k = 0; k < RADIATION_H2_SHIELDING_INVERSE_MAX_ITERATIONS; k++) {
    double slope;
    const double residual =
        radiation_h2_shielding_log_factor_slope(exp(u), alpha, b5, &slope) -
        log_t;
    if (residual > 0.)
      lo = u;
    else
      hi = u;

    const double step = slope < 0. ? -residual / slope : hi - lo;
    if (fabs(step) <= RADIATION_H2_SHIELDING_INVERSE_TOLERANCE)
      return exp(u + step);

    double u_next = u + step;
    if (!(u_next > lo && u_next < hi) || 2. * fabs(step) > fabs(step_before))
      u_next = 0.5 * (lo + hi);
    step_before = u_next - u;
    u = u_next;
    if (hi - lo <= RADIATION_H2_SHIELDING_INVERSE_TOLERANCE) break;
  }
  return exp(u);
}

/**
 * @brief Column whose shielding factor is a given transmission
 * (#radiation_h2_shielding_inverse_from without a first guess).
 *
 * @param t Transmission.
 * @param alpha Positive exponent of the fit (#radiation_h2_shielding_exponent).
 * @param b5 Thermal line width of H2 in km/s. Positive.
 * @return The column in cm^-2.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_inverse(const double t, const double alpha,
                               const double b5) {
  return radiation_h2_shielding_inverse_from(t, alpha, b5, 0.);
}

/**
 * @brief Add a column to a share, lower its transmission by the ratio of the
 * factor of the new column to the factor of the old one, and give the factor
 * of the new column in the fit of the receiver.
 *
 * The ratio never exceeds 1. A column step that is not positive leaves the
 * share as it is. Without a line width the share is left as it is and counts
 * as clear light: the factor is the one at column 0. The factor is never below
 * about 1e-58, so its inverse is finite.
 *
 * @param t (in, out) Transmission of the share.
 * @param N (in, out) Column of the share in cm^-2.
 * @param dN Column of the gas on the path of the hop in cm^-2.
 * @param alpha Exponent of the fit of the receiver.
 * @param b5 Thermal line width of H2 of the receiver in km/s.
 * @return The factor of the column after the hop.
 */
__attribute__((always_inline)) INLINE static double
radiation_h2_shielding_hop_factor(double *t, double *N, const double dN,
                                  const double alpha, const double b5) {
  if (!(b5 > 0.) || !(*N > 0. || dN > 0.))
    return radiation_h2_shielding_factor_at_zero();
  if (!(dN > 0.)) return radiation_h2_shielding_factor(*N, alpha, b5);
  const double f_old = radiation_h2_shielding_factor(*N, alpha, b5);
  const double f_new = radiation_h2_shielding_factor(*N + dN, alpha, b5);
  double N_new = *N + dN;
  if (N_new > RADIATION_H2_SHIELDING_COLUMN_CAP_CM2)
    N_new = RADIATION_H2_SHIELDING_COLUMN_CAP_CM2;
  *N = N_new;
  double ratio = f_new / f_old;
  if (ratio > 1.) ratio = 1.;
  *t *= ratio;
  return f_new;
}

/**
 * @brief The hop of #radiation_h2_shielding_hop_factor without the factor.
 *
 * @param t (in, out) Transmission of the share.
 * @param N (in, out) Column of the share in cm^-2.
 * @param dN Column of the gas on the path of the hop in cm^-2.
 * @param alpha Exponent of the fit of the receiver.
 * @param b5 Thermal line width of H2 of the receiver in km/s.
 */
__attribute__((always_inline)) INLINE static void radiation_h2_shielding_hop(
    double *t, double *N, const double dN, const double alpha,
    const double b5) {
  if (!(dN > 0.) || !(b5 > 0.)) return;
  radiation_h2_shielding_hop_factor(t, N, dN, alpha, b5);
}

/**
 * @brief Add the velocity shift of a hop to the shift of a share, and give the
 * share a fresh start if the sum exceeds the limit.
 *
 * With a fresh start the column is 0, the transmission 1 and the shift 0. A
 * limit that is not positive, for example T <= 0, or k_B or m_H not positive,
 * leaves the share as it is.
 *
 * @param t (in, out) Transmission of the share.
 * @param N (in, out) Column of the share in cm^-2.
 * @param dv (in, out) Velocity shift of the share.
 * @param dv_hop Shift added by the hop, signed, in the unit of @p limit.
 * @param limit Largest absolute shift that keeps the column.
 * @return 1 if the share started again, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int radiation_h2_shielding_shift(
    double *t, double *N, double *dv, const double dv_hop, const double limit) {
  if (!(limit > 0.)) return 0;
  const double shifted = *dv + dv_hop;
  if (fabs(shifted) <= limit) {
    *dv = shifted;
    return 0;
  }
  *t = 1.;
  *N = 0.;
  *dv = 0.;
  return 1;
}

/**
 * @brief #radiation_h2_shielding_shift for light whose transmission also
 * carries the dissociation factor of its stars: a fresh start drops only the
 * shielding of the column, `t S(0) / S(N)` in the fit of the receiver.
 *
 * The factor is exact for a column built in one fit. A factor that is not
 * positive, as for light whose transmission flushed to 0, is replaced by 1.
 *
 * @param t (in, out) Transmission of the share.
 * @param N (in, out) Column of the share in cm^-2.
 * @param dv (in, out) Velocity shift of the share.
 * @param dv_hop Shift added by the hop, signed, in the unit of @p limit.
 * @param limit Largest absolute shift that keeps the column.
 * @param alpha Exponent of the fit of the receiver.
 * @param b5 Thermal line width of H2 of the receiver in km/s.
 * @param t_max Bound of the dissociation factor.
 * @return 1 if the share started again, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
radiation_h2_shielding_shift_keep_factor(double *t, double *N, double *dv,
                                         const double dv_hop,
                                         const double limit, const double alpha,
                                         const double b5, const double t_max) {
  if (!(limit > 0.)) return 0;
  const double shifted = *dv + dv_hop;
  if (fabs(shifted) <= limit) {
    *dv = shifted;
    return 0;
  }
  double factor = *t;
  if (*N > 0. && b5 > 0.)
    factor *= radiation_h2_shielding_factor_at_zero() /
              radiation_h2_shielding_factor(*N, alpha, b5);
  if (!(factor > 0.)) factor = 1.;
  if (factor > t_max) factor = t_max;
  *t = factor;
  *N = 0.;
  *dv = 0.;
  return 1;
}

/**
 * @brief Quantities of the gas of a receiver that the hop of a share needs.
 *
 * A mass fraction below 0, a density, a temperature or a mean molecular
 * weight that is not positive gives no column and no line width, so a hop
 * through such a particle leaves a share as it is.
 *
 * @param X_H2 H2 mass fraction.
 * @param rho Gas density in g cm^-3.
 * @param mu Mean molecular weight.
 * @param T Gas temperature in K.
 * @param m_H Proton mass in g.
 * @param k_B Boltzmann constant in erg K^-1.
 * @param length_cm Centimetres in one unit of internal length.
 * @param velocity_cm_s Centimetres per second in one unit of internal
 * velocity.
 * @param column_per_length (return) H2 column per unit of internal length of
 * path, in cm^-2.
 * @param exponent (return) Exponent of the fit.
 * @param b5 (return) Thermal line width of H2 in km/s, 0 if T is not positive.
 * @param velocity_limit (return) Velocity shift in internal units beyond which
 * the column no longer shields, 0 if T is not positive.
 */
__attribute__((always_inline)) INLINE static void
radiation_h2_shielding_receiver_cache(
    const double X_H2, const double rho, const double mu, const double T,
    const double m_H, const double k_B, const double length_cm,
    const double velocity_cm_s, double *column_per_length, double *exponent,
    double *b5, double *velocity_limit) {

  double n_H2 = 0., n_gas = 0.;
  if (rho > 0. && m_H > 0.) {
    if (X_H2 > 0.) n_H2 = X_H2 * rho / (2. * m_H);
    n_gas = (mu > 0. ? rho / (mu * m_H) : rho / m_H);
  }
  *column_per_length = n_H2 * length_cm;
  *exponent = radiation_h2_shielding_exponent(T, n_gas);
  *b5 = 0.;
  *velocity_limit = 0.;
  if (T > 0. && k_B > 0. && m_H > 0.) {
    const double b_cgs = sqrt(k_B * T / m_H);
    *b5 = 1e-5 * b_cgs;
    *velocity_limit =
        RADIATION_H2_SHIELDING_VELOCITY_RESET_WIDTHS * b_cgs / velocity_cm_s;
  }
}

#endif /* SWIFT_RADIATION_H2_SHIELDING_GEAR_H */
