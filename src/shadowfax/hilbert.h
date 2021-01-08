/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file hilbert.h
 *
 * @brief Hilbert space-filling curve.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#ifndef SWIFT_HILBERT_H
#define SWIFT_HILBERT_H

static const unsigned int t2d[8][4][2] = {
    {{7, 0}, {0, 1}, {6, 3}, {0, 2}}, {{1, 2}, {7, 3}, {1, 1}, {6, 0}},
    {{2, 1}, {2, 2}, {4, 0}, {5, 3}}, {{4, 3}, {5, 0}, {3, 2}, {3, 1}},
    {{3, 3}, {4, 2}, {2, 0}, {4, 1}}, {{5, 1}, {3, 0}, {5, 2}, {2, 3}},
    {{6, 2}, {6, 1}, {0, 3}, {1, 0}}, {{0, 0}, {1, 3}, {7, 1}, {7, 2}}};

unsigned long hilbert_get_key_2d(unsigned long* bits, unsigned int nbits) {
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

unsigned long hilbert_get_key_3d(unsigned long* bits, unsigned int nbits) {
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

#endif  // SWIFT_HILBERT_H
