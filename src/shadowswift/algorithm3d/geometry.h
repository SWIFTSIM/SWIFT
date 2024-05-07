//
// Created by yuyttenh on 10/05/22.
//

#ifndef SWIFTSIM_GITLAB_GEOMETRY_3D_H
#define SWIFTSIM_GITLAB_GEOMETRY_3D_H

#include "error.h"
#include "shadowswift/algorithm3d/delaunay_vertex.h"
#include "shadowswift/algorithm3d/shadowswift_ray.h"
#include "shadowswift/delaunay.h"
#include "shadowswift/utils.h"

#include <float.h>
#include <gmp.h>
#include <math.h>
#ifdef DELAUNAY_3D_HAND_VEC
#include <immintrin.h>
#endif

/**
 * @brief Auxiliary variables used by the arbirary exact tests. Since allocating
 * and deallocating these variables poses a significant overhead, they are best
 * reused.
 */
struct geometry3d {
  /*! @brief Arbitrary exact vertex coordinates */
  mpz_t aix, aiy, aiz, bix, biy, biz, cix, ciy, ciz, dix, diy, diz, eix, eiy,
      eiz;

  /*! @brief Temporary variables used to store relative vertex coordinates. */
  mpz_t s1x, s1y, s1z, s2x, s2y, s2z, s3x, s3y, s3z, s4x, s4y, s4z;

  /*! @brief Temporary variables used to store intermediate results. */
  mpz_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

  /*! @brief Temporary variable used to store final exact results, before their
   *  sign is evaluated and returned as a finite precision integer. */
  mpz_t result;

  /*! @brief Temporary variable to store floating point values */
  mpf_t frac_n, frac_d, frac_result;
};

/**
 * @brief Initialize the geometry3d object.
 *
 * This allocates and initialises the auxiliary arbitrary precision variables.
 *
 * @param g Geometry object.
 */
inline static void geometry3d_init(struct geometry3d* restrict g) {
  mpz_inits(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
            g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
            g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
            g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->tmp3, g->tmp4, g->tmp5,
            g->tmp6, g->tmp7, g->tmp8, g->result, NULL);
  mpf_inits(g->frac_n, g->frac_d, g->frac_result, NULL);
}

/**
 * @brief Deallocate all memory occupied by the geometry3d object.
 *
 * @param g Geometry object.
 */
inline static void geometry3d_destroy(struct geometry3d* restrict g) {
  mpz_clears(g->aix, g->aiy, g->aiz, g->bix, g->biy, g->biz, g->cix, g->ciy,
             g->ciz, g->dix, g->diy, g->diz, g->eix, g->eiy, g->eiz, g->s1x,
             g->s1y, g->s1z, g->s2x, g->s2y, g->s2z, g->s3x, g->s3y, g->s3z,
             g->s4x, g->s4y, g->s4z, g->tmp1, g->tmp2, g->tmp3, g->tmp4,
             g->tmp5, g->tmp6, g->tmp7, g->tmp8, g->result, NULL);
  mpf_clears(g->frac_n, g->frac_d, g->frac_result, NULL);
}

/**
 * @brief Perform 4 orientation tests simultaneously using AVX2 intrinsics
 */
#ifdef DELAUNAY_3D_HAND_VEC
inline static void geometry3d_orient_simd(__m256d ax, __m256d ay, __m256d az,
                                          __m256d bx, __m256d by, __m256d bz,
                                          __m256d cx, __m256d cy, __m256d cz,
                                          __m256d dx, __m256d dy, __m256d dz,
                                          __m128i* tests) {
  /* Convert to relative coordinates */
  ax = _mm256_sub_pd(ax, dx);
  ay = _mm256_sub_pd(ay, dy);
  az = _mm256_sub_pd(az, dz);

  bx = _mm256_sub_pd(bx, dx);
  by = _mm256_sub_pd(by, dy);
  bz = _mm256_sub_pd(bz, dz);

  cx = _mm256_sub_pd(cx, dx);
  cy = _mm256_sub_pd(cy, dy);
  cz = _mm256_sub_pd(cz, dz);

  /* Compute product terms */
  __m256d bxcy = _mm256_mul_pd(bx, cy);
  __m256d cxby = _mm256_mul_pd(cx, by);

  __m256d cxay = _mm256_mul_pd(cx, ay);
  __m256d axcy = _mm256_mul_pd(ax, cy);

  __m256d axby = _mm256_mul_pd(ax, by);
  __m256d bxay = _mm256_mul_pd(bx, ay);

  /* Compute the result */
  __m256d result = _mm256_mul_pd(_mm256_sub_pd(bxcy, cxby), az);
  result = _mm256_add_pd(result, _mm256_mul_pd(_mm256_sub_pd(cxay, axcy), bz));
  result = _mm256_add_pd(result, _mm256_mul_pd(_mm256_sub_pd(axby, bxay), cz));

  /* Compute the error bounds */
  __m256d bound = _mm256_mul_pd(
      _mm256_add_pd(mm256_abs_pd(bxcy), mm256_abs_pd(cxby)), mm256_abs_pd(az));
  bound = _mm256_add_pd(bound, _mm256_mul_pd(_mm256_add_pd(mm256_abs_pd(cxay),
                                                           mm256_abs_pd(axcy)),
                                             mm256_abs_pd(bz)));
  bound = _mm256_add_pd(bound, _mm256_mul_pd(_mm256_add_pd(mm256_abs_pd(axby),
                                                           mm256_abs_pd(bxay)),
                                             mm256_abs_pd(cz)));
  static const double dbl_fac = 4. * DBL_EPSILON;
  static const __m256d fac = {dbl_fac, dbl_fac, dbl_fac, dbl_fac};
  bound = _mm256_mul_pd(bound, fac);

  /* Compute tests */
  __m256d sign = mm256_sgn_pd(result);
  __m256d abs_result = mm256_abs_pd(result);
  __m256d test = _mm256_cmp_pd(abs_result, bound, _CMP_GT_OS);
  *tests = _mm256_cvtpd_epi32(_mm256_blendv_pd(test, sign, test));
}
#endif

/**
 * @brief Inexact, but fast orientation test.
 *
 * This function calculates a maximal error bound on the result due to numerical
 * roundoff error. If the result is smaller than this error bound, the function
 * returns 0. Else this functions behaves similar to it's exact counterpart.
 *
 * */

inline static int geometry3d_orient(const double* restrict a,
                                    const double* restrict b,
                                    const double* restrict c,
                                    const double* restrict d) {
  /* Compute relative coordinates */
  const double adx = a[0] - d[0];
  const double ady = a[1] - d[1];
  const double adz = a[2] - d[2];

  const double bdx = b[0] - d[0];
  const double bdy = b[1] - d[1];
  const double bdz = b[2] - d[2];

  const double cdx = c[0] - d[0];
  const double cdy = c[1] - d[1];
  const double cdz = c[2] - d[2];

  /* Compute intermediate terms */
  const double bdxcdy = bdx * cdy;
  const double cdxbdy = cdx * bdy;

  const double cdxady = cdx * ady;
  const double adxcdy = adx * cdy;

  const double adxbdy = adx * bdy;
  const double bdxady = bdx * ady;

  /* Compute error bounds */
  double errbound = (fabs(bdxcdy) + fabs(cdxbdy)) * fabs(adz) +
                    (fabs(cdxady) + fabs(adxcdy)) * fabs(bdz) +
                    (fabs(adxbdy) + fabs(bdxady)) * fabs(cdz);
  // not really the right factor (probably too large), but this will do
  //  errbound *= 1.e-10;
  errbound *= DBL_EPSILON * 4;

  /* Compute result */
  double result = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
                  cdz * (adxbdy - bdxady);
  if (result < -errbound || result > errbound) {
    return sgn(result);
  }

  return 0;
}

/**
 * @brief Test the orientation of the tetrahedron that has the four given
 * points as vertex_indices.
 *
 * The test returns a positive result if the fourth vertex is below the plane
 * through the three other vertex_indices, with above the direction from which
 * the three points are ordered counterclockwise.
 *
 * E.g. if the four points are (0, 0, 0), (0, 0, 1), (0, 1, 0), and (1, 0, 0),
 * then this function returns 1.
 *
 * If the four points are exactly coplanar, then this function returns 0.
 *
 * @param a First vertex.
 * @param b Second vertex.
 * @param c Third vertex.
 * @param d Fourth vertex.
 * @return -1, 0, or 1, depending on the orientation of the tetrahedron.
 */
inline static int geometry3d_orient_exact(struct geometry3d* g,
                                          const unsigned long* a,
                                          const unsigned long* b,
                                          const unsigned long* c,
                                          const unsigned long* d) {

  /* store the input coordinates into the temporary large integer variables */
  /* also: explicitely set the f64's exponent bits to zero */
  mpz_set_ui(g->aix, a[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiy, a[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiz, a[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->bix, b[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biy, b[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biz, b[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->cix, c[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciy, c[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciz, c[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->dix, d[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diy, d[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diz, d[2] & 0xFFFFFFFFFFFFFllu);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->dix);
  mpz_sub(g->s1y, g->aiy, g->diy);
  mpz_sub(g->s1z, g->aiz, g->diz);

  mpz_sub(g->s2x, g->bix, g->dix);
  mpz_sub(g->s2y, g->biy, g->diy);
  mpz_sub(g->s2z, g->biz, g->diz);

  mpz_sub(g->s3x, g->cix, g->dix);
  mpz_sub(g->s3y, g->ciy, g->diy);
  mpz_sub(g->s3z, g->ciz, g->diz);

  /* Compute the result in 3 steps */
  mpz_set_ui(g->result, 0);

  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_addmul(g->result, g->s1z, g->tmp1);

  mpz_mul(g->tmp1, g->s3x, g->s1y);
  mpz_submul(g->tmp1, g->s1x, g->s3y);
  mpz_addmul(g->result, g->s2z, g->tmp1);

  mpz_mul(g->tmp1, g->s1x, g->s2y);
  mpz_submul(g->tmp1, g->s2x, g->s1y);
  mpz_addmul(g->result, g->s3z, g->tmp1);

  return mpz_sgn(g->result);
}

/**
 * @brief Adaptive 3D orientation test
 *
 * Returns +1 if the fourth point is below the plane through the
 * three other points, -1 if it is above and 0 if the four points
 * are coplanar.
 *
 * Above is defined as the direction from which the three points in the plane
 * are seen in counterclockwise order, e.g. if the points are A (0,0,0), B
 * (0,0,1), C (0,1,0) and D (1,0,0), then this function returns 1.
 *
 * @param a Coordinates of the first point, have to be in the interval [1,2]
 * @param b Coordinates of the second point, have to be in the interval [1,2]
 * @param c Coordinates of the third point, have to be in the interval [1,2]
 * @param d Coordinates of the fourth point, have to be in the interval [1,2]
 * @return A positive, negative or zero value, depending on the outcome of the
 * test
 */
inline static int geometry3d_orient_adaptive(
    struct geometry3d* restrict g, const delaunay_vertex_t* restrict a,
    const delaunay_vertex_t* restrict b, const delaunay_vertex_t* restrict c,
    const delaunay_vertex_t* restrict d) {

  int result = geometry3d_orient(a->x_f64, b->x_f64, c->x_f64, d->x_f64);

  if (result == 0) {
    result = geometry3d_orient_exact(g, a->x_u64, b->x_u64, c->x_u64, d->x_u64);
  }

  return result;
}

inline static void geometry3d_orient_4(struct geometry3d* restrict g,
                                       const delaunay_vertex_t* restrict a,
                                       const delaunay_vertex_t* restrict b,
                                       const delaunay_vertex_t* restrict c,
                                       const delaunay_vertex_t* restrict d,
                                       const delaunay_vertex_t* restrict e,
                                       int* tests) {
#ifndef DELAUNAY_3D_HAND_VEC
  error("Should not be calling this function!");
#else
  /* Get the data in the right format */
  __m256d v0x = {b->x_f64[0], a->x_f64[0], a->x_f64[0], a->x_f64[0]};
  __m256d v0y = {b->x_f64[1], a->x_f64[1], a->x_f64[1], a->x_f64[1]};
  __m256d v0z = {b->x_f64[2], a->x_f64[2], a->x_f64[2], a->x_f64[2]};

  __m256d v1x = {d->x_f64[0], c->x_f64[0], d->x_f64[0], b->x_f64[0]};
  __m256d v1y = {d->x_f64[1], c->x_f64[1], d->x_f64[1], b->x_f64[1]};
  __m256d v1z = {d->x_f64[2], c->x_f64[2], d->x_f64[2], b->x_f64[2]};

  __m256d v2x = {c->x_f64[0], d->x_f64[0], b->x_f64[0], c->x_f64[0]};
  __m256d v2y = {c->x_f64[1], d->x_f64[1], b->x_f64[1], c->x_f64[1]};
  __m256d v2z = {c->x_f64[2], d->x_f64[2], b->x_f64[2], c->x_f64[2]};

  __m256d v3x = {e->x_f64[0], e->x_f64[0], e->x_f64[0], e->x_f64[0]};
  __m256d v3y = {e->x_f64[1], e->x_f64[1], e->x_f64[1], e->x_f64[1]};
  __m256d v3z = {e->x_f64[2], e->x_f64[2], e->x_f64[2], e->x_f64[2]};

  /* Compute the tests */
  geometry3d_orient_simd(v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y,
                         v3z, (__m128i*)tests);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify result */
  int test = geometry3d_orient(b->x_f64, d->x_f64, c->x_f64, e->x_f64);
  if (test != tests[0]) error("AVX result doesn't match regular calculation!");
  test = geometry3d_orient(a->x_f64, c->x_f64, d->x_f64, e->x_f64);
  if (test != tests[1]) error("AVX result doesn't match regular calculation!");
  test = geometry3d_orient(a->x_f64, d->x_f64, b->x_f64, e->x_f64);
  if (test != tests[2]) error("AVX result doesn't match regular calculation!");
  test = geometry3d_orient(a->x_f64, b->x_f64, c->x_f64, e->x_f64);
  if (test != tests[3]) error("AVX result doesn't match regular calculation!");
#endif

  /* Do we need to run some exact tests? */
  if (tests[0] == 0) {
    tests[0] =
        geometry3d_orient_exact(g, b->x_u64, d->x_u64, c->x_u64, e->x_u64);
  }
  if (tests[1] == 0) {
    tests[1] =
        geometry3d_orient_exact(g, a->x_u64, c->x_u64, d->x_u64, e->x_u64);
  }
  if (tests[2] == 0) {
    tests[2] =
        geometry3d_orient_exact(g, a->x_u64, d->x_u64, b->x_u64, e->x_u64);
  }
  if (tests[3] == 0) {
    tests[3] =
        geometry3d_orient_exact(g, a->x_u64, b->x_u64, c->x_u64, e->x_u64);
  }
#endif
}

inline static int geometry3d_in_sphere_simd(const double* a, const double* b,
                                            const double* c, const double* d,
                                            const double* e) {
#ifndef DELAUNAY_3D_HAND_VEC
  error("Should not be calling this function!");
  return 0;
#else
  __m256d x = {a[0], b[0], c[0], d[0]};
  __m256d e_x = {e[0], e[0], e[0], e[0]};
  x = _mm256_sub_pd(x, e_x);

  __m256d y = {a[1], b[1], c[1], d[1]};
  __m256d e_y = {e[1], e[1], e[1], e[1]};
  y = _mm256_sub_pd(y, e_y);

  __m256d z = {a[2], b[2], c[2], d[2]};
  __m256d e_z = {e[2], e[2], e[2], e[2]};
  z = _mm256_sub_pd(z, e_z);

  // x^2 + y^2 + z^2
  __m256d n =
      _mm256_add_pd(_mm256_mul_pd(x, x),
                    _mm256_add_pd(_mm256_mul_pd(y, y), _mm256_mul_pd(z, z)));

  double err_bound;
  double result = determinant_4x4_bounded_simd(x, y, z, n, &err_bound);
  if (fabs(result) > err_bound) return sgn(result);

  return 0;
#endif
}

inline static int geometry3d_in_sphere(const double* restrict a,
                                       const double* restrict b,
                                       const double* restrict c,
                                       const double* restrict d,
                                       const double* restrict e) {

  /* Compute relative coordinates */
  const double aex = a[0] - e[0];
  const double aey = a[1] - e[1];
  const double aez = a[2] - e[2];

  const double bex = b[0] - e[0];
  const double bey = b[1] - e[1];
  const double bez = b[2] - e[2];

  const double cex = c[0] - e[0];
  const double cey = c[1] - e[1];
  const double cez = c[2] - e[2];

  const double dex = d[0] - e[0];
  const double dey = d[1] - e[1];
  const double dez = d[2] - e[2];

  /* Compute intermediate values */
  const double aexbey = aex * bey;
  const double bexaey = bex * aey;
  const double ab = aexbey - bexaey;
  const double bexcey = bex * cey;
  const double cexbey = cex * bey;
  const double bc = bexcey - cexbey;
  const double cexdey = cex * dey;
  const double dexcey = dex * cey;
  const double cd = cexdey - dexcey;
  const double dexaey = dex * aey;
  const double aexdey = aex * dey;
  const double da = dexaey - aexdey;
  const double aexcey = aex * cey;
  const double cexaey = cex * aey;
  const double ac = aexcey - cexaey;
  const double bexdey = bex * dey;
  const double dexbey = dex * bey;
  const double bd = bexdey - dexbey;

  const double abc = aez * bc - bez * ac + cez * ab;
  const double bcd = bez * cd - cez * bd + dez * bc;
  const double cda = cez * da + dez * ac + aez * cd;
  const double dab = dez * ab + aez * bd + bez * da;

  const double aenrm2 = aex * aex + aey * aey + aez * aez;
  const double benrm2 = bex * bex + bey * bey + bez * bez;
  const double cenrm2 = cex * cex + cey * cey + cez * cez;
  const double denrm2 = dex * dex + dey * dey + dez * dez;

  /* Compute errorbound */
  const double aezplus = fabs(aez);
  const double bezplus = fabs(bez);
  const double cezplus = fabs(cez);
  const double dezplus = fabs(dez);
  const double aexbeyplus = fabs(aexbey);
  const double bexaeyplus = fabs(bexaey);
  const double bexceyplus = fabs(bexcey);
  const double cexbeyplus = fabs(cexbey);
  const double cexdeyplus = fabs(cexdey);
  const double dexceyplus = fabs(dexcey);
  const double dexaeyplus = fabs(dexaey);
  const double aexdeyplus = fabs(aexdey);
  const double aexceyplus = fabs(aexcey);
  const double cexaeyplus = fabs(cexaey);
  const double bexdeyplus = fabs(bexdey);
  const double dexbeyplus = fabs(dexbey);

  double errbound = ((cexdeyplus + dexceyplus) * bezplus +
                     (dexbeyplus + bexdeyplus) * cezplus +
                     (bexceyplus + cexbeyplus) * dezplus) *
                        aenrm2 +
                    ((dexaeyplus + aexdeyplus) * cezplus +
                     (aexceyplus + cexaeyplus) * dezplus +
                     (cexdeyplus + dexceyplus) * aezplus) *
                        benrm2 +
                    ((aexbeyplus + bexaeyplus) * dezplus +
                     (bexdeyplus + dexbeyplus) * aezplus +
                     (dexaeyplus + aexdeyplus) * bezplus) *
                        cenrm2 +
                    ((bexceyplus + cexbeyplus) * aezplus +
                     (cexaeyplus + aexceyplus) * bezplus +
                     (aexbeyplus + bexaeyplus) * cezplus) *
                        denrm2;
  // not really the right factor (which is smaller), but this will do
  //  errbound *= 1.e-10;
  errbound *= DBL_EPSILON * 5.;

  /* Compute result */
  const double result =
      (denrm2 * abc - cenrm2 * dab) + (benrm2 * cda - aenrm2 * bcd);

  if (result < -errbound || result > errbound) {
    return sgn(result);
  }

  return 0;
}

/**
 * @brief Check if the fifth given point is inside (-1) the circumsphere of the
 * tetrahedron formed by the other four given points.
 *
 * It is assumed that the first four points are the vertex_indices of a
 * positively oriented tetrahedron, as defined by a negative return value of
 * orient3d().
 *
 * If the fifth point is exactly on the circumsphere of the tetrahedron, this
 * functions returns 0.
 *
 * @param g Geometry struct
 * @param a, b, c, d, e Coordinates of the vertices of the tetrahedron and the
 * test point (e).
 * @return -1, 0, or 1, depending on the outcome of the geometric test
 */
inline static int geometry3d_in_sphere_exact(struct geometry3d* restrict g,
                                             const unsigned long* restrict a,
                                             const unsigned long* restrict b,
                                             const unsigned long* restrict c,
                                             const unsigned long* restrict d,
                                             const unsigned long* restrict e) {
  /* store the input coordinates into the temporary large integer variables */
  /* We also need to erase the leading bits (due to rescaling) */
  mpz_set_ui(g->aix, a[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiy, a[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiz, a[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->bix, b[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biy, b[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biz, b[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->cix, c[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciy, c[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciz, c[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->dix, d[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diy, d[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diz, d[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->eix, e[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->eiy, e[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->eiz, e[2] & 0xFFFFFFFFFFFFFllu);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->aix, g->eix);
  mpz_sub(g->s1y, g->aiy, g->eiy);
  mpz_sub(g->s1z, g->aiz, g->eiz);

  mpz_sub(g->s2x, g->bix, g->eix);
  mpz_sub(g->s2y, g->biy, g->eiy);
  mpz_sub(g->s2z, g->biz, g->eiz);

  mpz_sub(g->s3x, g->cix, g->eix);
  mpz_sub(g->s3y, g->ciy, g->eiy);
  mpz_sub(g->s3z, g->ciz, g->eiz);

  mpz_sub(g->s4x, g->dix, g->eix);
  mpz_sub(g->s4y, g->diy, g->eiy);
  mpz_sub(g->s4z, g->diz, g->eiz);

  /* compute intermediate values */
  mpz_mul(g->tmp3, g->s1x, g->s2y);
  mpz_submul(g->tmp3, g->s2x, g->s1y);

  mpz_mul(g->tmp4, g->s2x, g->s3y);
  mpz_submul(g->tmp4, g->s3x, g->s2y);

  mpz_mul(g->tmp5, g->s3x, g->s4y);
  mpz_submul(g->tmp5, g->s4x, g->s3y);

  mpz_mul(g->tmp6, g->s4x, g->s1y);
  mpz_submul(g->tmp6, g->s1x, g->s4y);

  mpz_mul(g->tmp7, g->s1x, g->s3y);
  mpz_submul(g->tmp7, g->s3x, g->s1y);

  mpz_mul(g->tmp8, g->s2x, g->s4y);
  mpz_submul(g->tmp8, g->s4x, g->s2y);

  /* compute the result in 4 steps */
  mpz_set_ui(g->result, 0);

  mpz_mul(g->tmp1, g->s4x, g->s4x);
  mpz_addmul(g->tmp1, g->s4y, g->s4y);
  mpz_addmul(g->tmp1, g->s4z, g->s4z);
  mpz_mul(g->tmp2, g->s1z, g->tmp4);
  mpz_submul(g->tmp2, g->s2z, g->tmp7);
  mpz_addmul(g->tmp2, g->s3z, g->tmp3);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s3x, g->s3x);
  mpz_addmul(g->tmp1, g->s3y, g->s3y);
  mpz_addmul(g->tmp1, g->s3z, g->s3z);
  mpz_mul(g->tmp2, g->s4z, g->tmp3);
  mpz_addmul(g->tmp2, g->s1z, g->tmp8);
  mpz_addmul(g->tmp2, g->s2z, g->tmp6);
  mpz_submul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s2x, g->s2x);
  mpz_addmul(g->tmp1, g->s2y, g->s2y);
  mpz_addmul(g->tmp1, g->s2z, g->s2z);
  mpz_mul(g->tmp2, g->s3z, g->tmp6);
  mpz_addmul(g->tmp2, g->s4z, g->tmp7);
  mpz_addmul(g->tmp2, g->s1z, g->tmp5);
  mpz_addmul(g->result, g->tmp1, g->tmp2);

  mpz_mul(g->tmp1, g->s1x, g->s1x);
  mpz_addmul(g->tmp1, g->s1y, g->s1y);
  mpz_addmul(g->tmp1, g->s1z, g->s1z);
  mpz_mul(g->tmp2, g->s2z, g->tmp5);
  mpz_submul(g->tmp2, g->s3z, g->tmp8);
  mpz_addmul(g->tmp2, g->s4z, g->tmp4);
  mpz_submul(g->result, g->tmp1, g->tmp2);

  return mpz_sgn(g->result);
}

inline static int geometry3d_in_sphere_adaptive(struct geometry3d* restrict g,
                                                const delaunay_vertex_t* a,
                                                const delaunay_vertex_t* b,
                                                const delaunay_vertex_t* c,
                                                const delaunay_vertex_t* d,
                                                const delaunay_vertex_t* e) {

#ifdef DELAUNAY_3D_HAND_VEC
  int result = geometry3d_in_sphere_simd(a->x_f64, b->x_f64, c->x_f64, d->x_f64,
                                         e->x_f64);
#ifdef SWIFT_DEBUG_CHECKS
  int result_normal =
      geometry3d_in_sphere(a->x_f64, b->x_f64, c->x_f64, d->x_f64, e->x_f64);
  assert(result == result_normal);
#endif
#else
  int result =
      geometry3d_in_sphere(a->x_f64, b->x_f64, c->x_f64, d->x_f64, e->x_f64);
#endif
  if (result == 0) {
    result = geometry3d_in_sphere_exact(g, a->x_u64, b->x_u64, c->x_u64,
                                        d->x_u64, e->x_u64);
  }

  return result;
}

static inline int geometry3d_compute_circumcenter_relative_non_exact_no_errb(
    const double* restrict v0, const double* restrict v1,
    const double* restrict v2, const double* restrict v3,
    double* restrict circumcenter) {
  /* Compute relative coordinates */
  const double r1x = v1[0] - v0[0];
  const double r1y = v1[1] - v0[1];
  const double r1z = v1[2] - v0[2];
  const double r2x = v2[0] - v0[0];
  const double r2y = v2[1] - v0[1];
  const double r2z = v2[2] - v0[2];
  const double r3x = v3[0] - v0[0];
  const double r3y = v3[1] - v0[1];
  const double r3z = v3[2] - v0[2];

  /* Compute squared norm of relative coordinates */
  const double r1_sqrd = r1x * r1x + r1y * r1y + r1z * r1z;
  const double r2_sqrd = r2x * r2x + r2y * r2y + r2z * r2z;
  const double r3_sqrd = r3x * r3x + r3y * r3y + r3z * r3z;

  /* Compute a */
  const double r2yr3z = r2y * r3z;
  const double r3yr2z = r3y * r2z;
  const double r1yr3z = r1y * r3z;
  const double r3yr1z = r3y * r1z;
  const double r1yr2z = r1y * r2z;
  const double r2yr1z = r2y * r1z;
  const double a = r1x * (r2yr3z - r3yr2z) - r2x * (r1yr3z - r3yr1z) +
                   r3x * (r1yr2z - r2yr1z);

  /* Compute Dx */
  const double Dx = r1_sqrd * (r2yr3z - r3yr2z) - r2_sqrd * (r1yr3z - r3yr1z) +
                    r3_sqrd * (r1yr2z - r2yr1z);

  /* Compute Dy */
  const double r2xr3z = r2x * r3z;
  const double r3xr2z = r3x * r2z;
  const double r1xr3z = r1x * r3z;
  const double r3xr1z = r3x * r1z;
  const double r1xr2z = r1x * r2z;
  const double r2xr1z = r2x * r1z;
  const double Dy = -r1_sqrd * (r2xr3z - r3xr2z) + r2_sqrd * (r1xr3z - r3xr1z) -
                    r3_sqrd * (r1xr2z - r2xr1z);

  /* Compute Dz */
  const double r2xr3y = r2x * r3y;
  const double r3xr2y = r3x * r2y;
  const double r1xr3y = r1x * r3y;
  const double r3xr1y = r3x * r1y;
  const double r1xr2y = r1x * r2y;
  const double r2xr1y = r2x * r1y;
  const double Dz = r1_sqrd * (r2xr3y - r3xr2y) - r2_sqrd * (r1xr3y - r3xr1y) +
                    r3_sqrd * (r1xr2y - r2xr1y);

  /* Compute circumcenter */
  const double denominator = 2. * a;
  circumcenter[0] = Dx / denominator;
  circumcenter[1] = Dy / denominator;
  circumcenter[2] = Dz / denominator;
  return 1;
}

static inline int geometry3d_compute_circumcenter_relative_non_exact(
    double v0x, double v0y, double v0z, double v1x, double v1y, double v1z,
    double v2x, double v2y, double v2z, double v3x, double v3y, double v3z,
    double* circumcenter) {

  double errbound_factor = 1.e-10;
  // double errbound_factor = DBL_EPSILON * 4;

  /* Compute relative coordinates */
  const double r1x = v1x - v0x;
  const double r1y = v1y - v0y;
  const double r1z = v1z - v0z;
  const double r2x = v2x - v0x;
  const double r2y = v2y - v0y;
  const double r2z = v2z - v0z;
  const double r3x = v3x - v0x;
  const double r3y = v3y - v0y;
  const double r3z = v3z - v0z;

  /* Compute squared norm of relative coordinates */
  const double r1_sqrd = r1x * r1x + r1y * r1y + r1z * r1z;
  const double r2_sqrd = r2x * r2x + r2y * r2y + r2z * r2z;
  const double r3_sqrd = r3x * r3x + r3y * r3y + r3z * r3z;

  const double r2yr3z = r2y * r3z;
  const double r3yr2z = r3y * r2z;
  const double r1yr3z = r1y * r3z;
  const double r3yr1z = r3y * r1z;
  const double r1yr2z = r1y * r2z;
  const double r2yr1z = r2y * r1z;
  const double a = r1x * (r2yr3z - r3yr2z) - r2x * (r1yr3z - r3yr1z) +
                   r3x * (r1yr2z - r2yr1z);
  double errbound = fabs(r1x) * (fabs(r2yr3z) + fabs(r3yr2z)) +
                    fabs(r2x) * (fabs(r1yr3z) + fabs(r3yr1z)) +
                    fabs(r3x) * (fabs(r1yr2z) + fabs(r2yr1z));
  errbound *= errbound_factor;
  if (a >= -errbound && a <= errbound) return 0;

  /* Compute Dx */
  const double Dx = r1_sqrd * (r2yr3z - r3yr2z) - r2_sqrd * (r1yr3z - r3yr1z) +
                    r3_sqrd * (r1yr2z - r2yr1z);
  errbound = fabs(r1_sqrd) * (fabs(r2yr3z) - fabs(r3yr2z)) -
             fabs(r2_sqrd) * (fabs(r1yr3z) - fabs(r3yr1z)) +
             fabs(r3_sqrd) * (fabs(r1yr2z) - fabs(r2yr1z));
  errbound *= errbound_factor;
  if (Dx >= -errbound && Dx <= errbound) return 0;

  /* Compute Dy */
  const double r2xr3z = r2x * r3z;
  const double r3xr2z = r3x * r2z;
  const double r1xr3z = r1x * r3z;
  const double r3xr1z = r3x * r1z;
  const double r1xr2z = r1x * r2z;
  const double r2xr1z = r2x * r1z;
  const double Dy = -r1_sqrd * (r2xr3z - r3xr2z) + r2_sqrd * (r1xr3z - r3xr1z) -
                    r3_sqrd * (r1xr2z - r2xr1z);
  errbound = fabs(r1_sqrd) * (fabs(r2xr3z) - fabs(r3xr2z)) +
             fabs(r2_sqrd) * (fabs(r1xr3z) - fabs(r3xr1z)) -
             fabs(r3_sqrd) * (fabs(r1xr2z) - fabs(r2xr1z));
  errbound *= errbound_factor;
  if (Dy >= -errbound && Dy <= errbound) return 0;

  /* Compute Dz */
  const double r2xr3y = r2x * r3y;
  const double r3xr2y = r3x * r2y;
  const double r1xr3y = r1x * r3y;
  const double r3xr1y = r3x * r1y;
  const double r1xr2y = r1x * r2y;
  const double r2xr1y = r2x * r1y;
  const double Dz = r1_sqrd * (r2xr3y - r3xr2y) - r2_sqrd * (r1xr3y - r3xr1y) +
                    r3_sqrd * (r1xr2y - r2xr1y);
  errbound = fabs(r1_sqrd) * (fabs(r2xr3y) - fabs(r3xr2y)) -
             fabs(r2_sqrd) * (fabs(r1xr3y) - fabs(r3xr1y)) +
             fabs(r3_sqrd) * (fabs(r1xr2y) - fabs(r2xr1y));
  errbound *= errbound_factor;
  if (Dz >= -errbound && Dz <= errbound) return 0;

  const double denominator = 2. * a;
  circumcenter[0] = Dx / denominator;
  circumcenter[1] = Dy / denominator;
  circumcenter[2] = Dz / denominator;
  return 1;
}

static inline void geometry3d_compute_circumcenter_relative_exact(
    struct geometry3d* g, unsigned long ax, unsigned long ay, unsigned long az,
    unsigned long bx, unsigned long by, unsigned long bz, unsigned long cx,
    unsigned long cy, unsigned long cz, unsigned long dx, unsigned long dy,
    unsigned long dz, double* circumcenter) {
  /* store the input coordinates into the temporary large integer variables */
  mpz_set_ui(g->aix, ax & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiy, ay & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiz, az & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->bix, bx & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biy, by & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biz, bz & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->cix, cx & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciy, cy & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciz, cz & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->dix, dx & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diy, dy & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->diz, dz & 0xFFFFFFFFFFFFFllu);

  /* compute large integer relative coordinates */
  mpz_sub(g->s1x, g->bix, g->aix);
  mpz_sub(g->s1y, g->biy, g->aiy);
  mpz_sub(g->s1z, g->biz, g->aiz);

  mpz_sub(g->s2x, g->cix, g->aix);
  mpz_sub(g->s2y, g->ciy, g->aiy);
  mpz_sub(g->s2z, g->ciz, g->aiz);

  mpz_sub(g->s3x, g->dix, g->aix);
  mpz_sub(g->s3y, g->diy, g->aiy);
  mpz_sub(g->s3z, g->diz, g->aiz);

  /* Calculate denominator (->frac_d) */
  mpz_mul(g->tmp1, g->s2y, g->s3z);
  mpz_submul(g->tmp1, g->s3y, g->s2z);
  mpz_mul(g->tmp2, g->s1y, g->s3z);
  mpz_submul(g->tmp2, g->s3y, g->s1z);
  mpz_mul(g->tmp3, g->s1y, g->s2z);
  mpz_submul(g->tmp3, g->s2y, g->s1z);

  mpz_mul(g->tmp7, g->s1x, g->tmp1);
  mpz_submul(g->tmp7, g->s2x, g->tmp2);
  mpz_addmul(g->tmp7, g->s3x, g->tmp3);
  mpz_mul_ui(g->tmp8, g->tmp7, 2);
  assert(mpz_sgn(g->tmp8) != 0);
  mpf_set_z(g->frac_d, g->tmp8);

  /* Compute squared norm of relative coordinates */
  mpz_mul(g->tmp4, g->s1x, g->s1x);
  mpz_addmul(g->tmp4, g->s1y, g->s1y);
  mpz_addmul(g->tmp4, g->s1z, g->s1z);
  mpz_mul(g->tmp5, g->s2x, g->s2x);
  mpz_addmul(g->tmp5, g->s2y, g->s2y);
  mpz_addmul(g->tmp5, g->s2z, g->s2z);
  mpz_mul(g->tmp6, g->s3x, g->s3x);
  mpz_addmul(g->tmp6, g->s3y, g->s3y);
  mpz_addmul(g->tmp6, g->s3z, g->s3z);

  /* Calculate Dx (->frac_n) */
  mpz_mul(g->tmp7, g->tmp4, g->tmp1);
  mpz_submul(g->tmp7, g->tmp5, g->tmp2);
  mpz_addmul(g->tmp7, g->tmp6, g->tmp3);
  mpf_set_z(g->frac_n, g->tmp7);

  /* Calculate circumcenter[0] */
  mpf_div(g->frac_result, g->frac_n, g->frac_d);
  circumcenter[0] = mpf_get_d(g->frac_result) / 0x10000000000000llu;

  /* Calculate Dy (->frac_n) */
  mpz_mul(g->tmp1, g->s2x, g->s3z);
  mpz_submul(g->tmp1, g->s3x, g->s2z);
  mpz_mul(g->tmp2, g->s1x, g->s3z);
  mpz_submul(g->tmp2, g->s3x, g->s1z);
  mpz_mul(g->tmp3, g->s1x, g->s2z);
  mpz_submul(g->tmp3, g->s2x, g->s1z);

  mpz_mul(g->tmp7, g->tmp5, g->tmp2);
  mpz_submul(g->tmp7, g->tmp4, g->tmp1);
  mpz_submul(g->tmp7, g->tmp6, g->tmp3);
  mpf_set_z(g->frac_n, g->tmp7);

  /* Calculate circumcenter[1] */
  mpf_div(g->frac_result, g->frac_n, g->frac_d);
  circumcenter[1] = mpf_get_d(g->frac_result) / 0x10000000000000llu;

  /* Calculate Dz (->frac_n) */
  mpz_mul(g->tmp1, g->s2x, g->s3y);
  mpz_submul(g->tmp1, g->s3x, g->s2y);
  mpz_mul(g->tmp2, g->s1x, g->s3y);
  mpz_submul(g->tmp2, g->s3x, g->s1y);
  mpz_mul(g->tmp3, g->s1x, g->s2y);
  mpz_submul(g->tmp3, g->s2x, g->s1y);

  mpz_mul(g->tmp7, g->tmp4, g->tmp1);
  mpz_submul(g->tmp7, g->tmp5, g->tmp2);
  mpz_addmul(g->tmp7, g->tmp6, g->tmp3);
  mpf_set_z(g->frac_n, g->tmp7);

  /* Calculate circumcenter[2] */
  mpf_div(g->frac_result, g->frac_n, g->frac_d);
  circumcenter[2] = mpf_get_d(g->frac_result) / 0x10000000000000llu;
}

/*! @brief Compute circumcenter relative to v0 */
static inline void geometry3d_compute_circumcenter_relative_adaptive(
    struct geometry3d* restrict g, const delaunay_vertex_t* restrict v0,
    const delaunay_vertex_t* restrict v1, const delaunay_vertex_t* restrict v2,
    const delaunay_vertex_t* restrict v3, double* restrict circumcenter) {

#ifdef DELAUNAY_3D_ADAPTIVE_CIRCUMCENTER
  int result_non_exact = geometry3d_compute_circumcenter_relative_non_exact(
      v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0],
      v3[1], v3[2], circumcenter);

  if (!result_non_exact) {
    geometry3d_compute_circumcenter_relative_exact(
        g, v0ul[0], v0ul[1], v0ul[2], v1ul[0], v1ul[1], v1ul[2], v2ul[0],
        v2ul[1], v2ul[2], v3ul[0], v3ul[1], v3ul[2], circumcenter);
  }
#else
  geometry3d_compute_circumcenter_relative_non_exact_no_errb(
      v0->x_f64, v1->x_f64, v2->x_f64, v3->x_f64, circumcenter);
#endif
}

/**
 * @brief Compute the coordinates of the circumcenter of the tetrahedron
 * (v0, v1, v2, v3).
 *
 * See https://mathworld.wolfram.com/Circumsphere.html
 *
 * @param v0, v1, v2, v3 Rescaled coordinates of the corners of the tetrahedron.
 * @param v0ul, v1ul, v2ul, v3ul Integer coordinates of the corners of the
 * tetrahedron.
 * @param circumcenter (Returned) coordinates of center of circumsphere
 * @param box_side Side of box used to rescale coordinates
 * @param box_anchor Anchor of box used to rescale coordinates
 */
static inline void geometry3d_compute_circumcenter_adaptive(
    struct geometry3d* restrict g, const delaunay_vertex_t* restrict v0,
    const delaunay_vertex_t* restrict v1, const delaunay_vertex_t* restrict v2,
    const delaunay_vertex_t* restrict v3, double* restrict circumcenter,
    double box_side, const double* restrict box_anchor) {

  /* Calculate relative circumcenter (relative to v0, rescaled coordinates) */
  geometry3d_compute_circumcenter_relative_adaptive(g, v0, v1, v2, v3,
                                                    circumcenter);

  /* Translate and rescale circumcenter */
  circumcenter[0] =
      (circumcenter[0] + v0->x_f64[0] - 1.) * box_side + box_anchor[0];
  circumcenter[1] =
      (circumcenter[1] + v0->x_f64[1] - 1.) * box_side + box_anchor[1];
  circumcenter[2] =
      (circumcenter[2] + v0->x_f64[2] - 1.) * box_side + box_anchor[2];
}

static inline double geometry3d_compute_circumradius2_adaptive(
    struct geometry3d* restrict g, const delaunay_vertex_t* restrict v0,
    const delaunay_vertex_t* restrict v1, const delaunay_vertex_t* restrict v2,
    const delaunay_vertex_t* restrict v3, double box_side) {

  /* Calculate relative circumcenter (relative to v0, rescaled coordinates) */
  double circumcenter[3];
  geometry3d_compute_circumcenter_relative_adaptive(g, v0, v1, v2, v3,
                                                    circumcenter);

  /* Calculate and rescale radius */
  double radius = circumcenter[0] * circumcenter[0] +
                  circumcenter[1] * circumcenter[1] +
                  circumcenter[2] * circumcenter[2];
  return radius * box_side * box_side;
}

inline static double geometry3d_compute_area_triangle(double ax, double ay,
                                                      double az, double bx,
                                                      double by, double bz,
                                                      double cx, double cy,
                                                      double cz) {
  const double abx = bx - ax;
  const double aby = by - ay;
  const double abz = bz - az;
  const double acx = cx - ax;
  const double acy = cy - ay;
  const double acz = cz - az;

  const double Dx = aby * acz - abz * acy;
  const double Dy = abz * acx - abx * acz;
  const double Dz = abx * acy - aby * acx;

  return sqrt(Dx * Dx + Dy * Dy + Dz * Dz) / 2.;
}

inline static void geometry3d_compute_centroid_triangle(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double* result) {
  result[0] = (ax + bx + cx) / 3.;
  result[1] = (ay + by + cy) / 3.;
  result[2] = (az + bz + cz) / 3.;
}

inline static double geometry3d_compute_volume_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz) {
  /* Compute relative coordinates */
  const double dax = ax - dx;
  const double day = ay - dy;
  const double daz = az - dz;
  const double dbx = bx - dx;
  const double dby = by - dy;
  const double dbz = bz - dz;
  const double dcx = cx - dx;
  const double dcy = cy - dy;
  const double dcz = cz - dz;

  /* compute (b - d) x (c - d) */
  const double cross_x = dby * dcz - dcy * dbz;
  const double cross_y = dbx * dcz - dcx * dbz;
  const double cross_z = dbx * dcy - dcx * dby;

  return fabs(dax * cross_x - day * cross_y + daz * cross_z) / 6.;
}

inline static void geometry3d_compute_centroid_tetrahedron_exact(
    unsigned long ax, unsigned long ay, unsigned long az, unsigned long bx,
    unsigned long by, unsigned long bz, unsigned long cx, unsigned long cy,
    unsigned long cz, unsigned long dx, unsigned long dy, unsigned long dz,
    unsigned long* result) {
  /* x coordinate */
  unsigned long a_rem = ax % 4;
  unsigned long b_rem = bx % 4;
  unsigned long c_rem = cx % 4;
  unsigned long d_rem = dx % 4;
  unsigned long rem_sum = a_rem + b_rem + c_rem + d_rem;
  unsigned long average = ax / 4 + bx / 4 + cx / 4 + dx / 4 + rem_sum / 4;
  if (rem_sum % 4 > 1) average++;
  result[0] = average;

  /* y coordinate */
  a_rem = ay % 4;
  b_rem = by % 4;
  c_rem = cy % 4;
  d_rem = dy % 4;
  rem_sum = a_rem + b_rem + c_rem + d_rem;
  average = ay / 4 + by / 4 + cy / 4 + dy / 4 + rem_sum / 4;
  if (rem_sum % 4 > 1) average++;
  result[1] = average;

  /* z coordinate */
  a_rem = az % 4;
  b_rem = bz % 4;
  c_rem = cz % 4;
  d_rem = dz % 4;
  rem_sum = a_rem + b_rem + c_rem + d_rem;
  average = az / 4 + bz / 4 + cz / 4 + dz / 4 + rem_sum / 4;
  if (rem_sum % 4 > 1) average++;
  result[2] = average;
}

inline static void geometry3d_compute_centroid_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz, double* result) {
  result[0] = (ax + bx + cx + dx) / 4.;
  result[1] = (ay + by + cy + dy) / 4.;
  result[2] = (az + bz + cz + dz) / 4.;
}

inline static double geometry3d_compute_centroid_volume_tetrahedron(
    double ax, double ay, double az, double bx, double by, double bz, double cx,
    double cy, double cz, double dx, double dy, double dz, double* result) {
  geometry3d_compute_centroid_tetrahedron(ax, ay, az, bx, by, bz, cx, cy, cz,
                                          dx, dy, dz, result);
  return geometry3d_compute_volume_tetrahedron(ax, ay, az, bx, by, bz, cx, cy,
                                               cz, dx, dy, dz);
}

inline static double geometry3d_compute_centroid_area(
    const double* restrict points, int n_points, double* result) {

  if (n_points < 2) {
    error("Must pass at least 3 points!");
  }

  /* Calculate area and centroid from triangles (more robust) */
  double area = 0.;
  result[0] = 0.;
  result[1] = 0.;
  result[2] = 0.;

  const double v0x = points[0];
  const double v0y = points[1];
  const double v0z = points[2];

  for (int i = 2; i < n_points; i++) {
    double area_triangle = geometry3d_compute_area_triangle(
        v0x, v0y, v0z, points[3 * i - 3], points[3 * i - 2], points[3 * i - 1],
        points[3 * i], points[3 * i + 1], points[3 * i + 2]);
    area += area_triangle;

    double centroid_triangle[3];
    geometry3d_compute_centroid_triangle(
        v0x, v0y, v0z, points[3 * i - 3], points[3 * i - 2], points[3 * i - 1],
        points[3 * i], points[3 * i + 1], points[3 * i + 2], centroid_triangle);
    result[0] += area_triangle * centroid_triangle[0];
    result[1] += area_triangle * centroid_triangle[1];
    result[2] += area_triangle * centroid_triangle[2];
  }
  result[0] /= area;
  result[1] /= area;
  result[2] /= area;
  return area;
}

inline static void geometry3d_cross(const double* v1, const double* v2,
                                    double* restrict out_cross) {
  out_cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  out_cross[1] = v2[0] * v1[2] - v2[2] * v1[0];
  out_cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

inline static double geometry3d_dot(const double* v1, const double* v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/*! @returns the signed distance from the ray origin along the ray
 * direction to the intersection with the plane given by p1, p2, p3. */
inline static double geometry3d_ray_plane_intersect(
    const struct shadowswift_ray* r, const delaunay_vertex_t* restrict p1,
    const delaunay_vertex_t* restrict p2,
    const delaunay_vertex_t* restrict p3) {

  /* Setup useful variables */
  const double EPSILON = 1e-13;
  /* Vectors determining plane */
  const double v1[3] = {p1->x_f64[0] - p3->x_f64[0],
                        p1->x_f64[1] - p3->x_f64[1],
                        p1->x_f64[2] - p3->x_f64[2]};
  const double v2[3] = {p2->x_f64[0] - p3->x_f64[0],
                        p2->x_f64[1] - p3->x_f64[1],
                        p2->x_f64[2] - p3->x_f64[2]};
  /* Normal vector to plane */
  double n[3];
  geometry3d_cross(v1, v2, n);

  /* Compute result (see Camps 2013) */
  double numerator = n[0] * (p3->x_f64[0] - r->origin[0]) +
                     n[1] * (p3->x_f64[1] - r->origin[1]) +
                     n[2] * (p3->x_f64[2] - r->origin[2]);
  if (fabs(numerator) < EPSILON) {
    /* Point lies on the plane... */
    return 0.;
  }

  double denominator = geometry3d_dot(n, r->direction);
  if (fabs(denominator) < EPSILON) {
    /* Ray parallel to plane... */
    return DBL_MAX;
  }
  return numerator / denominator;
}

/*! @brief
 * returns -1 when exact test is needed */
inline static int geometry3d_ray_triangle_intersect_non_exact(
    const struct shadowswift_ray* restrict r,
    const delaunay_vertex_t* restrict p1, const delaunay_vertex_t* restrict p2,
    const delaunay_vertex_t* restrict p3, double* restrict out_distance) {

  const double errbound_factor = 1e-10;

  /* Setup useful variables */
  /* edges of triangle */
  const double e1[3] = {p2->x_f64[0] - p1->x_f64[0],
                        p2->x_f64[1] - p1->x_f64[1],
                        p2->x_f64[2] - p1->x_f64[2]};
  const double e2[3] = {p3->x_f64[0] - p1->x_f64[0],
                        p3->x_f64[1] - p1->x_f64[1],
                        p3->x_f64[2] - p1->x_f64[2]};

  const double dxe2y = r->direction[0] * e2[1];
  const double dye2x = r->direction[1] * e2[0];
  const double dxe2z = r->direction[0] * e2[2];
  const double dze2x = r->direction[2] * e2[0];
  const double dye2z = r->direction[1] * e2[2];
  const double dze2y = r->direction[2] * e2[1];

  const double h[3] = {dye2z - dze2y, dze2x - dxe2z, dxe2y - dye2x};

  const double a = geometry3d_dot(e1, h);
  double errbound = fabs(e1[0]) * (fabs(dxe2y) + fabs(dye2x)) +
                    fabs(e1[1]) * (fabs(dze2x) + fabs(dxe2z)) +
                    fabs(e1[2]) * (fabs(dye2z) + fabs(dze2y));
  errbound *= errbound_factor;
  if (-errbound < a && a < errbound) {
    /* Ray approximately parallel to triangle, try exact method */
    return -1;
  }

  double f = 1.0 / a;
  double s[3] = {r->origin[0] - p1->x_f64[0], r->origin[1] - p1->x_f64[1],
                 r->origin[2] - p1->x_f64[2]};
  double u = f * geometry3d_dot(s, h);
  errbound = fabs(f) * (fabs(e1[0]) * (fabs(dxe2y) + fabs(dye2x)) +
                        fabs(e1[1]) * (fabs(dze2x) + fabs(dxe2z)) +
                        fabs(e1[2]) * (fabs(dye2z) + fabs(dze2y)));
  errbound *= errbound_factor;
  if (-errbound < u && u < errbound) return -1;

  const double sxe1y = s[0] * e1[1];
  const double sye1x = s[1] * e1[0];
  const double sxe1z = s[0] * e1[2];
  const double sze1x = s[2] * e1[0];
  const double sye1z = s[1] * e1[2];
  const double sze1y = s[2] * e1[1];

  double q[3] = {sye1z - sze1y, sze1x - sxe1z, sxe1y - sye1x};
  double v = f * geometry3d_dot(r->direction, q);
  errbound = fabs(f) * (fabs(r->direction[0]) * (fabs(sye1z) + fabs(sze1y)) +
                        fabs(r->direction[1]) * (fabs(sze1x) + fabs(sxe1z)) +
                        fabs(r->direction[2]) * (fabs(sxe1y) + fabs(sye1x)));
  errbound *= errbound_factor;
  if (-errbound < v && v < errbound) return -1;

  double u_plus_v_minus_1 = u + v - 1;
  errbound = errbound_factor * (fabs(u) + fabs(v) + 1.);
  if (-errbound < u_plus_v_minus_1 && u_plus_v_minus_1 < errbound) {
    /* exact test needed */
    return -1;
  }

  double d = f * geometry3d_dot(e2, q);
  errbound = fabs(f) * (fabs(e2[0]) * (fabs(sye1z) + fabs(sze1y)) +
                        fabs(e2[1]) * (fabs(sze1x) + fabs(sxe1z)) +
                        fabs(e2[2]) * (fabs(sxe1y) + fabs(sye1x)));
  errbound *= errbound_factor;
  if (-errbound < d && d < errbound) return -1;

  *out_distance = d;
  return (u >= 0.0 && v >= 0.0 && u_plus_v_minus_1 <= 0.0);
}

inline static int geometry3d_ray_triangle_intersect_exact(
    struct geometry3d* g, const struct shadowswift_ray* restrict r,
    const delaunay_vertex_t* restrict p1, const delaunay_vertex_t* restrict p2,
    const delaunay_vertex_t* restrict p3, double* restrict out_distance) {
  /* Set exponent and sign bits of f64 to zero */
  mpz_set_ui(g->aix, p1->x_u64[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiy, p1->x_u64[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->aiz, p1->x_u64[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->bix, p2->x_u64[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biy, p2->x_u64[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->biz, p2->x_u64[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->cix, p3->x_u64[0] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciy, p3->x_u64[1] & 0xFFFFFFFFFFFFFllu);
  mpz_set_ui(g->ciz, p3->x_u64[2] & 0xFFFFFFFFFFFFFllu);

  mpz_set_ui(g->dix, r->origin_ul[0]);
  mpz_set_ui(g->diy, r->origin_ul[1]);
  mpz_set_ui(g->diz, r->origin_ul[2]);

  mpz_set_ui(g->eix, r->end_ul[0]);
  mpz_set_ui(g->eiy, r->end_ul[1]);
  mpz_set_ui(g->eiz, r->end_ul[2]);

  /* edge1 */
  mpz_sub(g->s1x, g->bix, g->aix);
  mpz_sub(g->s1y, g->biy, g->aiy);
  mpz_sub(g->s1z, g->biz, g->aiz);

  /* edge2 */
  mpz_sub(g->s2x, g->cix, g->aix);
  mpz_sub(g->s2y, g->ciy, g->aiy);
  mpz_sub(g->s2z, g->ciz, g->aiz);

  /* direction */
  mpz_sub(g->s3x, g->eix, g->dix);
  mpz_sub(g->s3y, g->eiy, g->diy);
  mpz_sub(g->s3z, g->eiz, g->diz);

  /* calculate h = direction x edge2 */
  mpz_mul(g->s4x, g->s3y, g->s2z);
  mpz_submul(g->s4x, g->s3z, g->s2y);
  mpz_mul(g->s4y, g->s3z, g->s2x);
  mpz_submul(g->s4y, g->s3x, g->s2z);
  mpz_mul(g->s4z, g->s3x, g->s2y);
  mpz_submul(g->s4z, g->s3y, g->s2x);

  /* calculate a = edge1 . h */
  mpz_mul(g->result, g->s1x, g->s4x);
  mpz_addmul(g->result, g->s1y, g->s4y);
  mpz_addmul(g->result, g->s1z, g->s4z);
  if (mpz_sgn(g->result) == 0) {
    /* Ray parallel to plane */
    *out_distance = INFINITY;
    return 0;
  }

  /* calculate s = ray_origin - p1 */
  mpz_sub(g->tmp1, g->dix, g->aix);
  mpz_sub(g->tmp2, g->diy, g->aiy);
  mpz_sub(g->tmp3, g->diz, g->aiz);

  /* calculate u * a = s . h */
  mpz_mul(g->tmp7, g->tmp1, g->s4x);
  mpz_addmul(g->tmp7, g->tmp2, g->s4y);
  mpz_addmul(g->tmp7, g->tmp3, g->s4z);

  /* calculate q = s x edge1 */
  mpz_mul(g->tmp4, g->tmp2, g->s1z);
  mpz_submul(g->tmp4, g->tmp3, g->s1y);
  mpz_mul(g->tmp5, g->tmp3, g->s1x);
  mpz_submul(g->tmp5, g->tmp1, g->s1z);
  mpz_mul(g->tmp6, g->tmp1, g->s1y);
  mpz_submul(g->tmp6, g->tmp2, g->s1x);

  /* calculate v * a = direction . q */
  mpz_mul(g->tmp8, g->s3x, g->tmp4);
  mpz_addmul(g->tmp8, g->s3y, g->tmp5);
  mpz_addmul(g->tmp8, g->s3z, g->tmp6);

  /* calculate a * (u + v) */
  mpz_add(g->tmp1, g->tmp7, g->tmp8);

  /* calculate a * distance = edge2 . q */
  mpz_mul(g->tmp2, g->s2x, g->tmp4);
  mpz_addmul(g->tmp2, g->s2y, g->tmp5);
  mpz_addmul(g->tmp2, g->s2z, g->tmp6);

  /* Calculate squared norm of direction */
  mpz_mul(g->tmp3, g->s3x, g->s3x);
  mpz_addmul(g->tmp3, g->s3y, g->s3y);
  mpz_addmul(g->tmp3, g->s3z, g->s3z);

  /* calculate distance */
  mpf_set_z(g->frac_n, g->tmp2);
  mpf_set_z(g->frac_d, g->result);
  mpf_div(g->frac_result, g->frac_n, g->frac_d);
  /* Multiply result by norm to compensate for un-normalized direction in a */
  mpf_set_z(g->frac_d, g->tmp3);
  mpf_sqrt(g->frac_n, g->frac_d);
  mpf_mul(g->frac_result, g->frac_result, g->frac_n);
  *out_distance = mpf_get_d(g->frac_result) / 0x10000000000000llu;

  /* intersects or not? (u >= 0, v >= 0, u + v <= 1.) */
  int sgn_a = mpz_sgn(g->result);
  if ((mpz_sgn(g->tmp7) == sgn_a) && (mpz_sgn(g->tmp8) == sgn_a)) {
    if (sgn_a > 0) {
      return mpz_cmp(g->tmp1, g->result) <= 0;
    } else {
      return mpz_cmp(g->tmp1, g->result) >= 0;
    }
  }
  return 0;
}

inline static int geometry3d_ray_triangle_intersect(
    struct geometry3d* restrict g, const struct shadowswift_ray* restrict r,
    const delaunay_vertex_t* restrict p1, const delaunay_vertex_t* restrict p2,
    const delaunay_vertex_t* restrict p3, double* restrict out_distance) {
  int intersection_result =
      geometry3d_ray_triangle_intersect_non_exact(r, p1, p2, p3, out_distance);

  if (intersection_result == -1) {
    return geometry3d_ray_triangle_intersect_exact(g, r, p1, p2, p3,
                                                   out_distance);
  }
  return intersection_result;
}

#endif  // SWIFTSIM_GITLAB_GEOMETRY_3D_H
