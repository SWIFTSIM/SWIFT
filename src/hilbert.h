/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Will Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_HILBERT_H
#define SWIFT_HILBERT_H

/* This header is adapted from cVoroni by Bert Vandenbroucke. */

/* Config parameters. */
#include <config.h>

/**
 * @brief 3D Hilbert key type able to represent up to 64 bits per dimension.
 *
 * We need 3 * 64 = 192 bits for a depth-64 Hilbert curve in 3D. We represent
 * this as three 64-bit words, with w[0] the most-significant word and w[2]
 * the least-significant word.
 */
typedef struct hilbert_key_3d {
  unsigned long w[3];
} hilbert_key_3d;

/**
 * @brief Zero-initialize a 3D Hilbert key.
 *
 * @param key Pointer to the key to clear.
 */
inline static void hilbert_key_3d_clear(hilbert_key_3d *key) {
  key->w[0] = 0UL;
  key->w[1] = 0UL;
  key->w[2] = 0UL;
}

/**
 * @brief Logical left shift of a 3D Hilbert key by s bits.
 *
 * Only values 0 < s < word size are valid. Internally we only use s == 3.
 *
 * @param key Pointer to the key to shift.
 * @param s   Number of bits to shift left.
 */
inline static void hilbert_key_3d_shl(hilbert_key_3d *key,
                                      const unsigned int s) {
  const unsigned int word_bits = (unsigned int)(sizeof(unsigned long) * 8U);

  /* Nothing to do if we shift by zero bits. */
  if (s == 0U) return;

  /* We never call this with s >= word_bits. */
  unsigned long carry = 0UL;

  /* Shift from least-significant word to most-significant, propagating carry.
   */
  for (int i = 2; i >= 0; --i) {
    const unsigned long v = key->w[i];
    const unsigned long new_carry = v >> (word_bits - s);
    key->w[i] = (v << s) | carry;
    carry = new_carry;
  }
}

/**
 * @brief key = (key << 3) + digit, where digit is an octant in [0, 7].
 *
 * This is the basic operation we need for building a Hilbert key one 3-bit
 * digit at a time.
 *
 * @param key   Pointer to the key to update.
 * @param digit Octant value to append (0 <= digit < 8).
 */
inline static void hilbert_key_3d_push_digit(hilbert_key_3d *key,
                                             const unsigned int digit) {
  /* Shift left by three bits to make room for the new digit. */
  hilbert_key_3d_shl(key, 3U);

  /* Add the new digit to the least-significant word. */
  const unsigned long d = (unsigned long)digit;
  unsigned long tmp = key->w[2] + d;
  key->w[2] = tmp;

  /* Propagate carry if we overflowed the least-significant word. */
  if (tmp < d) {
    ++key->w[1];
    if (key->w[1] == 0UL) {
      ++key->w[0];
    }
  }
}

/**
 * @brief Lexicographic comparison of two 192-bit 3D Hilbert keys.
 *
 * @param a First key.
 * @param b Second key.
 * @return -1 if a < b, 0 if a == b, +1 if a > b.
 */
inline static int compare_hilbert_key_3d(const hilbert_key_3d *a,
                                         const hilbert_key_3d *b) {
  if (a->w[0] < b->w[0]) return -1;
  if (a->w[0] > b->w[0]) return 1;

  if (a->w[1] < b->w[1]) return -1;
  if (a->w[1] > b->w[1]) return 1;

  if (a->w[2] < b->w[2]) return -1;
  if (a->w[2] > b->w[2]) return 1;

  return 0;
}

/**
 * @brief Sorting function used to sort vertex_indices on their Hilbert key.
 *
 * @param a First index.
 * @param b Second index.
 * @param x Hilbert key array to sort.
 * @return -1 if key[a] < key[b], 0 if equal, +1 if key[a] > key[b].
 */
inline static int sort_h_comp(const void *a, const void *b, void *x) {
  const int ai = *(const int *)a;
  const int bi = *(const int *)b;
  const hilbert_key_3d *keys = (const hilbert_key_3d *)x;
  return compare_hilbert_key_3d(&keys[ai], &keys[bi]);
}

/* 3D Hilbert state transition table. */
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

/**
 * @brief Compute a 3D Hilbert key for coordinates up to 64 bits per axis.
 *
 * The input coordinates are given as bitfields in @p bits[0..2].
 * The @p nbits parameter sets how many of the most-significant bits of each
 * coordinate are used to build the Hilbert index.
 *
 * The returned key is 192 bits wide and can represent up to nbits == 64
 * without collisions.
 *
 * This implementation builds the key using an explicit 192-bit accumulator
 * in three 64-bit words and performs a branchless 3-bit left shift at each
 * step to append the next octant digit.
 *
 * @param bits  Array of length 3 with the x, y, z coordinates as unsigned long.
 * @param nbits Number of bits per coordinate to use (1 <= nbits <= 64).
 *
 * @return 3D Hilbert key as a 192-bit value (hilbert_key_3d).
 */
inline static hilbert_key_3d hilbert_get_key_3d(const unsigned long *bits,
                                                unsigned int nbits) {

  hilbert_key_3d key;
  hilbert_key_3d_clear(&key);

  const unsigned int word_bits = (unsigned int)(sizeof(unsigned long) * 8U);

  /* Guard against bogus input. */
  if (nbits == 0U) {
    return key;
  }
  if (nbits > word_bits) {
    nbits = word_bits;
  }

  /* Local 192-bit accumulator in three words (MSW w0, then w1, LSW w2). */
  unsigned long w0 = 0UL;
  unsigned long w1 = 0UL;
  unsigned long w2 = 0UL;

  /* We always shift by three bits per Hilbert digit. */
  const unsigned int s = 3U;
  const unsigned int r = word_bits - s;

  /* Start with a mask on the most-significant bit we actually use. */
  unsigned long mask = 1UL;
  mask <<= (nbits - 1U);

  /* Initial state, as in the original implementation. */
  unsigned int si = 4U;

  for (unsigned int i = nbits; i--;) {

    /* Extract the next bit of each coordinate. */
    const unsigned int x0 = (bits[0] & mask) != 0UL;
    const unsigned int x1 = (bits[1] & mask) != 0UL;
    const unsigned int x2 = (bits[2] & mask) != 0UL;

    const unsigned int ci = (x0 << 2) | (x1 << 1) | x2;
    const unsigned int digit = t3d[si][ci][1];

    /* Branchless 192-bit left shift by 3 bits and append new digit. */
    const unsigned long new_w0 = (w0 << s) | (w1 >> r);
    const unsigned long new_w1 = (w1 << s) | (w2 >> r);
    const unsigned long new_w2 = (w2 << s) | (unsigned long)digit;

    w0 = new_w0;
    w1 = new_w1;
    w2 = new_w2;

    /* Update state for the next level. */
    si = t3d[si][ci][0];

    /* Move on to the next less-significant bit. */
    mask >>= 1U;
  }

  key.w[0] = w0;
  key.w[1] = w1;
  key.w[2] = w2;

  return key;
}

#endif /* SWIFT_HILBERT_H */
