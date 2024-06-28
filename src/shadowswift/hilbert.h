//
// Created by yuyttenh on 22/04/22.
//

#ifndef SWIFTSIM_HILBERT_H
#define SWIFTSIM_HILBERT_H

/**
 * @brief Comparison function for two unsigned long values.
 *
 * @param a First value.
 * @param b Second value.
 * @return -1 if a < b, 0 if a == b, +1 if a > b.
 */
inline static int compare_unsigned_long(const unsigned long a,
                                        const unsigned long b) {
  if (a < b) {
    return -1;
  } else {
    if (a > b) {
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 * @brief Sorting function used to sort vertex_indices on their Hilbert key.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Hilbert key array to sort.
 * @return Return value of compare_unsigned_long() for the hilbert key of
 * vertex_indices a and b.
 */
inline static int sort_h_comp(const void *a, const void *b, void *x) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  unsigned long *keys = (unsigned long *)x;
  unsigned long ah = keys[ai];
  unsigned long bh = keys[bi];
  return compare_unsigned_long(ah, bh);
}

static const unsigned int t2d[8][4][2] = {
    {{7, 0}, {0, 1}, {6, 3}, {0, 2}}, {{1, 2}, {7, 3}, {1, 1}, {6, 0}},
    {{2, 1}, {2, 2}, {4, 0}, {5, 3}}, {{4, 3}, {5, 0}, {3, 2}, {3, 1}},
    {{3, 3}, {4, 2}, {2, 0}, {4, 1}}, {{5, 1}, {3, 0}, {5, 2}, {2, 3}},
    {{6, 2}, {6, 1}, {0, 3}, {1, 0}}, {{0, 0}, {1, 3}, {7, 1}, {7, 2}}};

inline static unsigned long hilbert_get_key_2d(unsigned long *bits,
                                               unsigned int nbits) {
  unsigned long key = 0;
  unsigned long mask = 1;
  mask <<= (nbits - 1);
  unsigned int x[2];
  unsigned int ci;
  unsigned int si = 7;
  for (unsigned int i = nbits; i--;) {
    key <<= 2;
    x[0] = (bits[0] & mask) > 0;
    x[1] = (bits[1] & mask) > 0;
    ci = (x[0] << 1) | x[1];
    key += t2d[si][ci][1];
    si = t2d[si][ci][0];
    mask >>= 1;
  }
  return key;
}

static const unsigned int t3d[12][8][2] = {
    {{5, 0}, {1, 7}, {4, 1}, {2, 6}, {3, 3}, {3, 4}, {4, 2}, {2, 5}},
    {{6, 4}, {2, 7}, {6, 3}, {8, 0}, {0, 5}, {0, 6}, {7, 2}, {7, 1}},
    {{1, 6}, {0, 7}, {1, 5}, {9, 4}, {10, 1}, {11, 0}, {10, 2}, {9, 3}},
    {{9, 2}, {8, 5}, {0, 3}, {0, 4}, {9, 1}, {8, 6}, {6, 0}, {10, 7}},
    {{0, 0}, {5, 1}, {8, 3}, {5, 2}, {11, 7}, {6, 6}, {8, 4}, {6, 5}},
    {{4, 0}, {10, 3}, {9, 7}, {10, 4}, {0, 1}, {0, 2}, {7, 6}, {7, 5}},
    {{11, 6}, {11, 5}, {3, 1}, {3, 2}, {4, 7}, {1, 4}, {9, 0}, {1, 3}},
    {{9, 6}, {8, 1}, {5, 7}, {1, 0}, {9, 5}, {8, 2}, {11, 4}, {11, 3}},
    {{1, 2}, {4, 3}, {1, 1}, {7, 0}, {10, 5}, {4, 4}, {10, 6}, {3, 7}},
    {{2, 4}, {5, 5}, {7, 7}, {5, 6}, {2, 3}, {6, 2}, {3, 0}, {6, 1}},
    {{11, 2}, {11, 1}, {3, 5}, {3, 6}, {5, 3}, {2, 0}, {5, 4}, {8, 7}},
    {{7, 4}, {7, 3}, {4, 5}, {2, 2}, {6, 7}, {10, 0}, {4, 6}, {2, 1}}};

inline static unsigned long hilbert_get_key_3d(unsigned long *bits,
                                               unsigned int nbits) {
  unsigned long key = 0;
  unsigned long mask = 1;
  mask <<= nbits - 1;
  unsigned int x[3];
  unsigned int ci;
  unsigned int si = 4;
  for (unsigned int i = nbits; i--;) {
    key <<= 3;
    x[0] = (bits[0] & mask) > 0;
    x[1] = (bits[1] & mask) > 0;
    x[2] = (bits[2] & mask) > 0;
    ci = (x[0] << 2) | (x[1] << 1) | x[2];
    key += t3d[si][ci][1];
    si = t3d[si][ci][0];
    mask >>= 1;
  }
  return key;
}

#endif  // SWIFTSIM_HILBERT_H
