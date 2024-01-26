//
// Created by yuyttenh on 1/26/24.
//

#ifndef SWIFTSIM_SHADOWSWIFT_QSELECT_H
#define SWIFTSIM_SHADOWSWIFT_QSELECT_H

#include "memswap.h"

/** @brief A modified version of the Hoare partitioning scheme to be used in the
 * QuickSelect algorithm.
 *
 * In particular this version makes sure the original pivot value is at the
 * returned position, which is assumed in the quickselect algorithm.
 *
 * NOTE: This function assumes that there is at least one v <= pivot and
 * one v >= pivot in the remaining values.
 * (this allows us to skip the i <= j check in the while loops)
 *
 * @param pivot The index of the pivot element. NOTE: Must be selected so that
 * the pivot value is *not* strictly the largest or smallest of all values.
 * @param left The start index of the slice values to partition.
 * @param right The end index (inclusive) of the slice of values to partition.
 * @param values
 * @param size The size in bytes of a single value
 * @param compar Function used to compare two values. compar(a, b, arg) should
 * return -1 for a < b, 0 for a == b and 1 for a > b.
 * @param arg Extra information required by the compar function.
 * @returns The final index of the pivot. The slice of values will also be
 * reordered such that all values before the pivot are smaller than or equal to
 * the pivot value and all values after it larger than or equal.
 */
inline static int partition_r(int pivot, int left, int right, void *values,
                              size_t size,
                              int (*compar)(const void *, const void *, void *),
                              void *arg) {
  // First swap pivot to the left
  void *p = values + left * size;
  memswap(p, values + pivot * size, size);

  // Now partition the left + 1:right subarray using our pivot (using the hoare
  // partition scheme).
  int i = left;
  int j = right + 1;
  void *a, *b;
  while (1) {
    // find i and j such that they must be inverted
    a = values + i * size;
    b = values + j * size;
    do {
      i++;
      a += size;
    } while (/* i <= j &&  */ compar(a, p, arg) < 0);
    do {
      j--;
      b -= size;
    } while (/* i <= j &&  */ compar(b, p, arg) > 0);
    if (i >= j) {
      // Don't return j just jet, we must still insert the pivot.
      break;
    }
    // swap the values located at index i and j
    memswap(a, b, size);
  }

  // Now the array looks as follows: [pivot:A:B], with pivot >= a, for all a in
  // A and pivot <= b for all b in B Swap the pivot to the end of A to obtain
  // the desired partition.
  memswap(p, b, size);
  return j;
}

/** @brief Computes a pivot as the median value of the start, end and middle
 * value of the slice.
 *
 * This ensures that our pivot satisfies the requirements of
 * our partition_r function and also protects against poor scaling of the
 * QuickSelect algorithm for (partially) sorted data.
 *
 * @param left The start index of the slice values to partition.
 * @param right The end index (inclusive) of the slice of values to partition.
 * We assume that right >= left + 2.
 * @param values
 * @param size The size in bytes of a single value
 * @param compar Function used to compare two values. compar(a, b, arg) should
 * return -1 for a < b, 0 for a == b and 1 for a > b.
 * @param arg Extra information required by the compar function.
 * @returns The index of the desired pivot.
 */
inline static int median_3_pivot_r(
    int left, int right, const void *values, size_t size,
    int (*compar)(const void *, const void *, void *), void *arg) {
  int mid = left + (right - left) / 2;
  const void *l = values + left * size;
  const void *m = values + mid * size;
  const void *r = values + right * size;
  const int lm = compar(l, m, arg);
  const int lr = compar(l, r, arg);
  const int mr = compar(m, r, arg);
  if (lm <= 0 && mr <= 0) return mid;
  if (lm >= 0 && mr >= 0) return mid;
  if (lm >= 0 && lr <= 0) return left;
  if (lm <= 0 && lr >= 0) return left;
  // (mr <= 0 && lr >= 0) || (mr >= 0 && lr <= 0)
  return right;
}

/** @brief An implementation of the QuickSelect algorithm.
 *
 * This function reorders a values array such that for all i < k,
 * values[i] <= values[k] and for all i > k values[i] >= k. I.e. the (k+1)-th
 * element of the array is located at index k.
 * This can be used to efficiently (in linear time) compute the median of a
 * given array.
 *
 * @param k The index around which to order the elements.
 * @param values The array to order.
 * @param count The length of #values
 * @param size The size in bytes of the values
 * @param compar Function used to compare two values. compar(a, b, arg) should
 * return -1 for a < b, 0 for a == b and 1 for a > b.
 * @param arg Extra information required by the #compar function.
 * @returns VOID. The array will be reordered as described above once this
 * function completes.
 */
inline static void qselect_r(int k, void *values, size_t count, size_t size,
                             int (*compar)(const void *, const void *, void *),
                             void *arg) {
  int n_iter = 0;
  int left = 0;
  int right = count - 1;
  // QuickSelect should have linear time complexity
  const int max_iter = 10 * count * count;
  while (n_iter++ < max_iter) {
    // Base cases
    if (left == right) {
      return;
    } else if (left + 1 == right) {
      void *l = values + left * size;
      void *r = values + right * size;
      if (compar(l, r, arg) > 0) {
        memswap(l, r, size);
      }
      return;
    }
#ifdef SWIFT_DEBUG_CHECKS
    if (left > right)
      error("Trying to run QuickSelect step with left > right!");
#endif
    int pivot = median_3_pivot_r(left, right, values, size, compar, arg);
    pivot = partition_r(pivot, left, right, values, size, compar, arg);
    if (k == pivot) {
      return;
    } else if (k < pivot) {
      right = pivot - 1;
    } else {  // pivot < k
      left = pivot + 1;
    }
  }
  error("QuickSelect failed to converge after %d iterations!", max_iter);
}
#endif  // SWIFTSIM_SHADOWSWIFT_QSELECT_H
