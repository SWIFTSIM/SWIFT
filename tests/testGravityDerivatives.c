/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

/*************************/
/* 0th order derivatives */
/*************************/

/**
 * @brief \f$ \phi(r_x, r_y, r_z) \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_000(double r_x, double r_y, double r_z, double r_inv) {

  return r_inv;
}

/*************************/
/* 1st order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_100(double r_x, double r_y, double r_z, double r_inv) {

  return -r_x * r_inv * r_inv * r_inv;
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_010(double r_x, double r_y, double r_z, double r_inv) {

  return -r_y * r_inv * r_inv * r_inv;
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_001(double r_x, double r_y, double r_z, double r_inv) {

  return -r_z * r_inv * r_inv * r_inv;
}

/*************************/
/* 2nd order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_200(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_x * r_x * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_y^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_020(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_y * r_y * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_z^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_002(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_z * r_z * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x\partial r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_110(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_x * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_101(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_x * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_y\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_011(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_y * r_z * r_inv5;
}

/*************************/
/* 3rd order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_300(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_x * r_inv7 + 9. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_030(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_y * r_y * r_y * r_inv7 + 9. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_z^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_003(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_z * r_z * r_z * r_inv7 + 9. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^2\partial r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_210(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_y * r_inv7 + 3. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^2\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_201(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_z * r_inv7 + 3. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x\partial r_y^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_120(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_y * r_y * r_inv7 + 3. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y^2\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_021(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_z * r_y * r_y * r_inv7 + 3. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x\partial r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_102(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_z * r_z * r_inv7 + 3. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y\partial r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_012(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_y * r_z * r_z * r_inv7 + 3. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_z\partial
 * r_y\partial r_z} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
INLINE static double D_111(double r_x, double r_y, double r_z, double r_inv) {
  const double r_inv3 = r_inv * r_inv * r_inv;
  const double r_inv7 = r_inv3 * r_inv3 * r_inv;
  return -15. * r_x * r_y * r_z * r_inv7;
}

/*********************************/
/* 4th order gravity derivatives */
/*********************************/

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_z^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_004(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^1 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_013(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^2 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_022(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^3 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_031(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_040(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_y * r_y) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_103(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^1 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_112(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_y);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^2 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_121(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_z);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_130(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_y);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_202(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_y^1 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_211(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_z);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_y^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_220(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_y) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^3 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_301(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^3 \partial_y^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_310(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_y);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_400(double r_x, double r_y, double r_z, double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_x * r_x) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/*********************************/
/* 5th order gravity derivatives */
/*********************************/

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_z^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_005(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_z * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_z);
  /* 26 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^1 \partial_z^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_014(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^2 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_023(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^3 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_032(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^4 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_041(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_050(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_y);
  /* 26 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_z^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_104(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^1 \partial_z^3
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_113(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^2 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_122(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_x);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^3 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_131(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_140(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_203(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^1 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_212(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_y);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^2 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_221(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_230(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_302(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_y^1 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_311(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_y^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_320(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^4 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_401(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_x * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^4 \partial_y^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_410(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_x * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
INLINE static double D_500(double r_x, double r_y, double r_z, double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_x * r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_x);
  /* 26 zero-valued terms not written out */
}

void test(double x, double y, double tol, double min, const char* name) {

  double diff = fabs(x - y);
  double norm = 0.5 * fabs(x + y);
  if (diff > norm * tol && norm > min)
    error(
        "Relative difference (%e) for '%s' (swift=%e) and (exact=%e) exceeds "
        "tolerance (%e)",
        diff / norm, name, x, y, tol);
  /* else */
  /*   message("'%s' (%e -- %e) OK!", name, x, y); */
}

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Relative tolerance */
  double tol = 1e-4;

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  for (int i = 0; i < 100; ++i) {

    const double dx = 100. * ((double)rand() / (RAND_MAX));
    const double dy = 100. * ((double)rand() / (RAND_MAX));
    const double dz = 100. * ((double)rand() / (RAND_MAX));

    message("Testing gravity for r=(%e %e %e)", dx, dy, dz);

    /* Compute distance */
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double r_inv = 1. / sqrt(r2);
    const double r = r2 * r_inv;
    const double eps = r / 10.;
    const double eps_inv = 1. / eps;

    /* Compute all derivatives */
    struct potential_derivatives_M2L pot;
    compute_potential_derivatives_M2L(dx, dy, dz, r2, r_inv, eps, eps_inv,
                                      &pot);

    /* Minimal value we care about */
    const double min = 1e-9;

    /* Now check everything... */

    /* 0th order terms */
    test(pot.D_000, D_000(dx, dy, dz, r_inv), tol, min, "D_000");

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

    /* 1st order terms */
    test(pot.D_100, D_100(dx, dy, dz, r_inv), tol, min, "D_100");
    test(pot.D_010, D_010(dx, dy, dz, r_inv), tol, min, "D_010");
    test(pot.D_001, D_001(dx, dy, dz, r_inv), tol, min, "D_001");
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

    /* 2nd order terms */
    test(pot.D_200, D_200(dx, dy, dz, r_inv), tol, min, "D_200");
    test(pot.D_020, D_020(dx, dy, dz, r_inv), tol, min, "D_020");
    test(pot.D_002, D_002(dx, dy, dz, r_inv), tol, min, "D_002");
    test(pot.D_110, D_110(dx, dy, dz, r_inv), tol, min, "D_110");
    test(pot.D_101, D_101(dx, dy, dz, r_inv), tol, min, "D_101");
    test(pot.D_011, D_011(dx, dy, dz, r_inv), tol, min, "D_011");
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

    tol *= 2.;

    /* 3rd order terms */
    test(pot.D_300, D_300(dx, dy, dz, r_inv), tol, min, "D_300");
    test(pot.D_030, D_030(dx, dy, dz, r_inv), tol, min, "D_030");
    test(pot.D_003, D_003(dx, dy, dz, r_inv), tol, min, "D_003");
    test(pot.D_210, D_210(dx, dy, dz, r_inv), tol, min, "D_210");
    test(pot.D_201, D_201(dx, dy, dz, r_inv), tol, min, "D_201");
    test(pot.D_120, D_120(dx, dy, dz, r_inv), tol, min, "D_120");
    test(pot.D_021, D_021(dx, dy, dz, r_inv), tol, min, "D_021");
    test(pot.D_102, D_102(dx, dy, dz, r_inv), tol, min, "D_102");
    test(pot.D_012, D_012(dx, dy, dz, r_inv), tol, min, "D_012");
    test(pot.D_111, D_111(dx, dy, dz, r_inv), tol, min, "D_111");
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

    /* 4th order terms */
    test(pot.D_400, D_400(dx, dy, dz, r_inv), tol, min, "D_400");
    test(pot.D_040, D_040(dx, dy, dz, r_inv), tol, min, "D_040");
    test(pot.D_004, D_004(dx, dy, dz, r_inv), tol, min, "D_004");
    test(pot.D_310, D_310(dx, dy, dz, r_inv), tol, min, "D_310");
    test(pot.D_301, D_301(dx, dy, dz, r_inv), tol, min, "D_301");
    test(pot.D_130, D_130(dx, dy, dz, r_inv), tol, min, "D_130");
    test(pot.D_031, D_031(dx, dy, dz, r_inv), tol, min, "D_031");
    test(pot.D_103, D_103(dx, dy, dz, r_inv), tol, min, "D_103");
    test(pot.D_013, D_013(dx, dy, dz, r_inv), tol, min, "D_013");
    test(pot.D_220, D_220(dx, dy, dz, r_inv), tol, min, "D_220");
    test(pot.D_202, D_202(dx, dy, dz, r_inv), tol, min, "D_202");
    test(pot.D_022, D_022(dx, dy, dz, r_inv), tol, min, "D_022");
    test(pot.D_211, D_211(dx, dy, dz, r_inv), tol, min, "D_211");
    test(pot.D_121, D_121(dx, dy, dz, r_inv), tol, min, "D_121");
    test(pot.D_112, D_112(dx, dy, dz, r_inv), tol, min, "D_112");
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

    tol *= 2.;

    /* 5th order terms */
    test(pot.D_500, D_500(dx, dy, dz, r_inv), tol, min, "D_500");
    test(pot.D_050, D_050(dx, dy, dz, r_inv), tol, min, "D_050");
    test(pot.D_005, D_005(dx, dy, dz, r_inv), tol, min, "D_005");
    test(pot.D_410, D_410(dx, dy, dz, r_inv), tol, min, "D_410");
    test(pot.D_401, D_401(dx, dy, dz, r_inv), tol, min, "D_401");
    test(pot.D_140, D_140(dx, dy, dz, r_inv), tol, min, "D_140");
    test(pot.D_041, D_041(dx, dy, dz, r_inv), tol, min, "D_041");
    test(pot.D_104, D_104(dx, dy, dz, r_inv), tol, min, "D_104");
    test(pot.D_014, D_014(dx, dy, dz, r_inv), tol, min, "D_014");
    test(pot.D_320, D_320(dx, dy, dz, r_inv), tol, min, "D_320");
    test(pot.D_302, D_302(dx, dy, dz, r_inv), tol, min, "D_302");
    test(pot.D_230, D_230(dx, dy, dz, r_inv), tol, min, "D_230");
    test(pot.D_032, D_032(dx, dy, dz, r_inv), tol, min, "D_032");
    test(pot.D_203, D_203(dx, dy, dz, r_inv), tol, min, "D_203");
    test(pot.D_023, D_023(dx, dy, dz, r_inv), tol, min, "D_023");
    test(pot.D_311, D_311(dx, dy, dz, r_inv), tol, min, "D_311");
    test(pot.D_131, D_131(dx, dy, dz, r_inv), tol, min, "D_131");
    test(pot.D_113, D_113(dx, dy, dz, r_inv), tol, min, "D_113");
    test(pot.D_122, D_122(dx, dy, dz, r_inv), tol, min, "D_122");
    test(pot.D_212, D_212(dx, dy, dz, r_inv), tol, min, "D_212");
    test(pot.D_221, D_221(dx, dy, dz, r_inv), tol, min, "D_221");

#endif
    message("All good!");
  }
  return 0;
}
