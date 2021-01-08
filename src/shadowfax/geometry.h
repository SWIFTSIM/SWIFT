/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file geometry.h
 *
 * @brief Arbitrary exact and non-exact geometrical tests.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#ifndef SWIFT_GEOMETRY_H
#define SWIFT_GEOMETRY_H

#include <gmp.h>

/**
 * @brief Auxiliary variables used by the arbirary exact tests. Since allocating
 * and deallocating these variables poses a significant overhead, they are best
 * reused.
 */
struct geometry {
  /*! @brief Arbitrary exact vertex coordinates */
  mpz_t aix, aiy, bix, biy, cix, ciy, dix, diy;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s2x, s2y, s3x, s3y;

  /*! @brief Temporary variables used to store intermediate results. */
  mpz_t tmp1, tmp2;

  /*! @brief Temporary variable used to store final exact results, before their
   *  sign is evaluated and returned as a finite precision integer. */
  mpz_t result;
};

/**
 * @brief Initialize the geometry object.
 *
 * This allocates and initialises the auxiliary arbitrary precision variables.
 *
 * @param g Geometry object.
 */
inline static void geometry_init(struct geometry* restrict g) {
  mpz_inits(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
            g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y, g->tmp1, g->tmp2,
            g->result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry object.
 *
 * @param g Geometry object.
 */
inline static void geometry_destroy(struct geometry* restrict g) {
  mpz_clears(g->aix, g->aiy, g->bix, g->biy, g->cix, g->ciy, g->dix, g->diy,
             g->s1x, g->s1y, g->s2x, g->s2y, g->s3x, g->s3y, g->tmp1, g->tmp2,
             g->result, NULL);
}

/**
 * @brief Non-exact 2D orientation test.
 *
 * This test returns a positive value if the points a, b and c are ordered
 * counterclockwise, a negative value if they are ordered clockwise, and 0 if
 * they are colinear. The return value corresponds to twice the signed area
 * of the triangle formed by a, b and c, as computed from (Springel, 2010):
 * @f[
 *   \begin{vmatrix}
 *     1 & a_x & a_y \\
 *     1 & b_x & b_y \\
 *     1 & c_x & c_y
 *   \end{vmatrix}
 * @f]
 *
 * Due to roundoff error, this test can be inconsistent with other geometric
 * tests within the same floating point coordinate basis.
 *
 * @param ax, ay First point.
 * @param bx, by Second point.
 * @param cx, cy Third point.
 * @return Signed double area of the triangle formed by a, b and c.
 */
inline static double geometry_orient2d(double ax, double ay, double bx,
                                       double by, double cx, double cy) {

  /* the code below stays as close as possible to the implementation of the
     exact test below */
  double s1x, s1y, s2x, s2y, result;

  /* compute the relative positions of a and b w.r.t. */
  s1x = ax - cx;
  s1y = ay - cy;

  s2x = bx - cx;
  s2y = by - cy;

  /* compute the result in two steps */
  result = 0.;
  result += s1x * s2y;
  result -= s1y * s2x;

  return result;
}

/**
 * @brief Arbitrary exact alternative for geometry_orient2d().
 *
 * This function calculates exactly the same thing as the non-exact version, but
 * does so in an integer coordinate basis, using arbitrarily large integers.
 * Since we are not interested in the value of the final result, but only in its
 * sign, this function only returns -1, 0, or 1, depending on the sign of the
 * signed integer triangle area.
 *
 * @param g Geometry object (containing temporary variables that will be used).
 * @param ax, ay, bx, by, cx, cy Integer coordinates of the three points.
 */
inline static int geometry_orient2d_exact(
    struct geometry* restrict g, unsigned long int ax, unsigned long int ay,
    unsigned long int bx, unsigned long int by, unsigned long int cx,
    unsigned long int cy) {

  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->cix);
  mpz_sub(g->s1y, g->aiy, g->ciy);

  mpz_sub(g->s2x, g->bix, g->cix);
  mpz_sub(g->s2y, g->biy, g->ciy);

  /* now compute the result using the same 2 steps as the non-exact test */
  mpz_set_ui(g->result, 0);
  mpz_addmul(g->result, g->s1x, g->s2y);
  mpz_submul(g->result, g->s1y, g->s2x);

  /* evaluate the sign of result and return */
  return mpz_sgn(g->result);
}

/**
 * @brief Non-exact 2D in-circle test.
 *
 * This test determines whether the point d is inside the circle through the
 * points a, b and c. Assuming a, b and c are positively oriented (positive
 * return value of geometry_orient2d()), a negative return value means d is
 * outside the circle and a positive value it is inside. A return value of 0
 * means d lies on the circle.
 *
 * This test computes the following (Springel, 2010):
 * @f[
 *   \begin{vmatrix}
 *     1 & a_x & a_y & a_x^2 + a_y^2 \\
 *     1 & b_x & b_y & b_x^2 + b_y^2 \\
 *     1 & c_x & c_y & c_x^2 + c_y^2 \\
 *     1 & d_x & d_y & d_x^2 + d_y^2
 *   \end{vmatrix}
 * @f]
 *
 * Due to roundoff error, this test can be inconsistent with other geometric
 * tests within the same floating point coordinate basis.
 *
 * @param ax, ay First point.
 * @param bx, by Second point.
 * @param cx, cy Third point.
 * @param dx, dy Fourth point.
 */
inline static double geometry_in_circle(double ax, double ay, double bx,
                                        double by, double cx, double cy,
                                        double dx, double dy) {

  /* the code below stays as close as possible to the implementation of the
     exact test below */
  double s1x, s1y, s2x, s2y, s3x, s3y, tmp1, tmp2, result;

  /* compute the relative coordinates of a, b and c w.r.t. d */
  s1x = ax - dx;
  s1y = ay - dy;

  s2x = bx - dx;
  s2y = by - dy;

  s3x = cx - dx;
  s3y = cy - dy;

  /* compute the result in 3 steps */
  result = 0.;

  /* accumulate some terms in tmp1 and tmp2 and then update the result */
  tmp1 = s2x * s3y;
  tmp1 -= s3x * s2y;
  tmp2 = s1x * s1x;
  tmp2 += s1y * s1y;
  result += tmp1 * tmp2;

  tmp1 = s3x * s1y;
  tmp1 -= s1x * s3y;
  tmp2 = s2x * s2x;
  tmp2 += s2y * s2y;
  result += tmp1 * tmp2;

  tmp1 = s1x * s2y;
  tmp1 -= s2x * s1y;
  tmp2 = s3x * s3x;
  tmp2 += s3y * s3y;
  result += tmp1 * tmp2;

  return result;
}

/**
 * @brief Arbitrary exact alternative for geometry_in_circle().
 *
 * This function calculates exactly the same thing as the non-exact version, but
 * does so in an integer coordinate basis, using arbitrarily large integers.
 * Since we are not interested in the value of the final result, but only in its
 * sign, this function only returns -1, 0, or 1, depending on the sign of the
 * signed integer triangle area.
 *
 * @param g Geometry object (containing temporary variables that will be used).
 * @param ax, ay, bx, by, cx, cy, dx, dy Integer coordinates of the four points.
 */
inline static int geometry_in_circle_exact(
    struct geometry* restrict g, unsigned long int ax, unsigned long int ay,
    unsigned long int bx, unsigned long int by, unsigned long int cx,
    unsigned long int cy, unsigned long int dx, unsigned long int dy) {

  /* copy the coordinate values into the large integer temporary variables */
  mpz_set_ui(g->aix, ax);
  mpz_set_ui(g->aiy, ay);

  mpz_set_ui(g->bix, bx);
  mpz_set_ui(g->biy, by);

  mpz_set_ui(g->cix, cx);
  mpz_set_ui(g->ciy, cy);

  mpz_set_ui(g->dix, dx);
  mpz_set_ui(g->diy, dy);

  /* compute the relative coordinates using large integers */
  mpz_sub(g->s1x, g->aix, g->dix);
  mpz_sub(g->s1y, g->aiy, g->diy);

  mpz_sub(g->s2x, g->bix, g->dix);
  mpz_sub(g->s2y, g->biy, g->diy);

  mpz_sub(g->s3x, g->cix, g->dix);
  mpz_sub(g->s3y, g->ciy, g->diy);

  /* compute the result using the same 3 steps as in the non-exact version */
  mpz_set_ui(g->result, 0);

  /* accumulate temporary terms in tmp1 and tmp2 and update result */
  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_mul(g->tmp2, g->s1x, g->s1x);
  mpz_addmul(g->tmp2, g->s1y, g->s1y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s3x, g->s1y);
  mpz_submul(g->tmp1, g->s1x, g->s3y);
  mpz_mul(g->tmp2, g->s2x, g->s2x);
  mpz_addmul(g->tmp2, g->s2y, g->s2y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s1x, g->s2y);
  mpz_submul(g->tmp1, g->s2x, g->s1y);
  mpz_mul(g->tmp2, g->s3x, g->s3x);
  mpz_addmul(g->tmp2, g->s3y, g->s3y);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  /* evaluate the sign of the result and return */
  return mpz_sgn(g->result);
}

#endif /* SWIFT_GEOMETRY_H */
