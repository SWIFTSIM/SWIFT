/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MATRIX_H
#define SWIFT_MATRIX_H

/* Local includes */
#include "error.h"

#if !defined(HYDRO_DIMENSION_3D)
#error "Not yet defined!"
#endif

struct sym_matrix {

  union {
    struct {
      float elements[6];
    };
    struct {
      float xx;
      float yy;
      float zz;
      float xy;
      float xz;
      float yz;
    };
  };
};

__attribute__((always_inline)) INLINE static void zero_sym_matrix(
    struct sym_matrix *M) {
  for (int i = 0; i < 6; ++i) M->elements[i] = 0;
}

__attribute__((always_inline)) INLINE static void get_matrix_from_sym_matrix(
    float out[3][3], const struct sym_matrix *in) {

  out[0][0] = in->xx;
  out[0][1] = in->xy;
  out[0][2] = in->xz;
  out[1][0] = in->xy;
  out[1][1] = in->yy;
  out[1][2] = in->yz;
  out[2][0] = in->xz;
  out[2][1] = in->yz;
  out[2][2] = in->zz;
}

__attribute__((always_inline)) INLINE static void get_sym_matrix_from_matrix(
    struct sym_matrix *out, const float in[3][3]) {
  out->xx = in[0][0];
  out->yy = in[1][1];
  out->zz = in[2][2];
  out->xy = in[0][1];
  out->xz = in[0][2];
  out->yz = in[1][2];
}

__attribute__((always_inline)) INLINE static void sym_matrix_multiply_by_vector(
    float out[3], const struct sym_matrix *M, const float v[3]) {

  out[0] = M->xx * v[0] + M->xy * v[1] + M->xz * v[2];
  out[1] = M->xy * v[0] + M->yy * v[1] + M->yz * v[2];
  out[2] = M->xz * v[0] + M->yz * v[1] + M->zz * v[2];
}

__attribute__((always_inline)) INLINE static void sym_matrix_print(
    const struct sym_matrix *M) {
  message("|%.4f %.4f %.4f|", M->xx, M->xy, M->xz);
  message("|%.4f %.4f %.4f|", M->xy, M->yy, M->yz);
  message("|%.4f %.4f %.4f|", M->xz, M->yz, M->zz);
}

#endif /* SWIFT_MATRIX_H */
