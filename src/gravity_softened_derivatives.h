/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "inline.h"
#include "kernel_gravity.h"

#if 0

/*************************/
/* 0th order derivatives */
/*************************/

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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv3 = eps_inv * eps_inv2;
  const double eps_inv5 = eps_inv3 * eps_inv2;
  return r_x * r_x * eps_inv5 * D_soft_2(u) - eps_inv3 * D_soft_1(u);
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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv3 = eps_inv * eps_inv2;
  const double eps_inv5 = eps_inv3 * eps_inv2;
  return r_y * r_y * eps_inv5 * D_soft_2(u) - eps_inv3 * D_soft_1(u);
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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv3 = eps_inv * eps_inv2;
  const double eps_inv5 = eps_inv3 * eps_inv2;
  return r_z * r_z * eps_inv5 * D_soft_2(u) - eps_inv3 * D_soft_1(u);
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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  return r_x * r_y * eps_inv5 * D_soft_2(u);
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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  return r_x * r_z * eps_inv5 * D_soft_2(u);
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
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  return r_y * r_z * eps_inv5 * D_soft_2(u);
}

/*************************/
/* 3rd order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_x^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_300(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_x * r_x * r_x * eps_inv7 * D_soft_3(u) +
         3. * r_x * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_y^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_030(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_y * r_y * r_y * eps_inv7 * D_soft_3(u) +
         3. * r_y * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_z^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_003(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_z * r_z * r_z * eps_inv7 * D_soft_3(u) +
         3. * r_z * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_x^2\partial
 * r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_210(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_x * r_x * r_y * eps_inv7 * D_soft_3(u) +
         r_y * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_x^2\partial
 * r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_201(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_x * r_x * r_z * eps_inv7 * D_soft_3(u) +
         r_z * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_x\partial
 * r_y^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_120(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_x * r_y * r_y * eps_inv7 * D_soft_3(u) +
         r_x * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_y^2\partial
 * r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_021(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_y * r_y * r_z * eps_inv7 * D_soft_3(u) +
         r_z * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_x\partial
 * r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_102(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_x * r_z * r_z * eps_inv7 * D_soft_3(u) +
         r_x * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_y\partial
 * r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_012(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv5 = eps_inv2 * eps_inv2 * eps_inv;
  const double eps_inv7 = eps_inv5 * eps_inv2;
  return -r_y * r_z * r_z * eps_inv7 * D_soft_3(u) +
         r_y * eps_inv5 * D_soft_2(u);
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z, h)}{\partial r_z\partial
 * r_y\partial r_z} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r Norm of the distance vector (\f$ |r| \f$).
 * @param eps_inv Inverse of the softening length (\f$ 1/h \f$).
 */
__attribute__((always_inline)) INLINE static double D_soft_111(
    double r_x, double r_y, double r_z, double r, double eps_inv) {

  const double u = r * eps_inv;
  const double eps_inv2 = eps_inv * eps_inv;
  const double eps_inv4 = eps_inv2 * eps_inv2;
  const double eps_inv7 = eps_inv4 * eps_inv2 * eps_inv;
  return -r_x * r_y * r_z * eps_inv7 * D_soft_3(u);
}

#endif

#endif /* SWIFT_GRAVITY_SOFTENED_DERIVATIVE_H */
