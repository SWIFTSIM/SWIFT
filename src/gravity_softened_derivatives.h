/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_GRAVITY_SOFTENED_DERIVATIVE_H
#define SWIFT_GRAVITY_SOFTENED_DERIVATIVE_H

/**
 * @file gravity_softened_derivatives.h
 * @brief Derivatives of the softened gravitational potential.
 *
 * We use the notation of Dehnen, Computational Astrophysics and Cosmology,
 * 1, 1, pp. 24 (2014), arXiv:1405.2255
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "inline.h"

/*************************/
/* 0th order derivatives */
/*************************/

__attribute__((always_inline)) INLINE static double D_soft_0(double u) {

  /* phi(u) = 3u^7 - 15u^6 + 28u^5 - 21u^4 + 7u^2 - 3 */
  double phi = 3. * u - 15.;
  phi = phi * u + 28.;
  phi = phi * u - 21.;
  phi = phi * u;
  phi = phi * u + 7.;
  phi = phi * u;
  phi = phi * u - 3.;

  return phi;
}

/**
 * @brief \f$ \phi(r_x, r_y, r_z, h) \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_000(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return eps_inv * D_soft_0(u);
}

/*************************/
/* 1st order derivatives */
/*************************/

__attribute__((always_inline)) INLINE static double D_soft_1(double u) {

  /* phi(u) = 21u^6 - 90u^5 + 140u^4 - 84u^3 + 14u */
  double phi = 21. * u - 90.;
  phi = phi * u + 140.;
  phi = phi * u - 84.;
  phi = phi * u;
  phi = phi * u + 14.;
  phi = phi * u;

  return phi;
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z, h)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_100(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return -r_x * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z, h)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_010(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return -r_y * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z, h)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_001(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return -r_z * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/*************************/
/* 2nd order derivatives */
/*************************/

__attribute__((always_inline)) INLINE static double D_soft_2(double u) {

  /* phi(u) = 126u^5 - 450u^4 + 560u^3 - 252u^2 + 14 */
  double phi = 126. * u - 450.;
  phi = phi * u + 560.;
  phi = phi * u - 252.;
  phi = phi * u;
  phi = phi * u + 14.;

  return phi;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_x^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_200(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_x * r_x * eps_inv * eps_inv * D_soft_2(u) +
         (r_y * r_y + r_z * r_z) * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_y^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_020(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_y * r_y * eps_inv * eps_inv * D_soft_2(u) +
         (r_x * r_x + r_z * r_z) * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_z^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_002(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_z * r_z * eps_inv * eps_inv * D_soft_2(u) +
         (r_x * r_x + r_y * r_y) * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_x\partial r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_110(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_x * r_y * eps_inv * eps_inv * D_soft_2(u) -
         r_x * r_y * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_x\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_101(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_x * r_z * eps_inv * eps_inv * D_soft_2(u) -
         r_x * r_z * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z, h)}{\partial r_y\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_011(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  return r_y * r_z * eps_inv * eps_inv * D_soft_2(u) -
         r_y * r_z * eps_inv * eps_inv * eps_inv * D_soft_1(u);
}

#endif /* SWIFT_GRAVITY_SOFTENED_DERIVATIVE_H */
