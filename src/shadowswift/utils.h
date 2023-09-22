//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_UTILS_H
#define SWIFTSIM_SHADOWSWIFT_UTILS_H

#include <float.h>
#include <math.h>

static inline int sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

/**
 * @brief Check whether two doubles are equal up to the given precision.
 *
 * @param double1
 * @param double2
 * @param precision
 * @return 1 for equality 0 else.
 */
inline static int double_cmp(double double1, double double2,
                             unsigned long precision) {
  long long1, long2;
  if (double1 > 0)
    long1 = (long)(double1 * precision + .5);
  else
    long1 = (long)(double1 * precision - .5);
  if (double2 > 0)
    long2 = (long)(double2 * precision + .5);
  else
    long2 = (long)(double2 * precision - .5);
  return (long1 == long2);
}

inline static int approx_equals(double a, double b, double epsilon) {
  if (a == b) return 1;

  double abs_a = fabs(a);
  double abs_b = fabs(b);
  double diff = fabs(a - b);

  if (a == 0. || b == 0. || abs_a + abs_b < DBL_EPSILON) {
    return diff < epsilon * DBL_EPSILON;
  } else {
    return diff / fmin(abs_a + abs_b, DBL_MAX) < epsilon;
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_UTILS_H
