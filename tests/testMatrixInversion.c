/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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
#include <stdlib.h>
#include <string.h>

/* Local headers */
#include "const.h"
#include "dimension.h"
#include "error.h"
#include "tools.h"

void setup_matrix(float A[3][3]) {
  A[0][0] = random_uniform(-1.0, 1.0);
  A[0][1] = random_uniform(-1.0, 1.0);
  A[0][2] = random_uniform(-1.0, 1.0);
  A[1][0] = random_uniform(-1.0, 1.0);
  A[1][1] = random_uniform(-1.0, 1.0);
  A[1][2] = random_uniform(-1.0, 1.0);
  A[2][0] = random_uniform(-1.0, 1.0);
  A[2][1] = random_uniform(-1.0, 1.0);
  A[2][2] = random_uniform(-1.0, 1.0);
}

int is_unit_matrix(float A[3][3]) {
  int check = 1;

  check &= (fabsf(A[0][0] - 1.0f) < 1.e-6f);

#if defined(HYDRO_DIMENSION_2D) && defined(HYDRO_DIMENSION_3D)
  check &= (fabsf(A[0][1]) < 1.e-6f);
  check &= (fabsf(A[1][0]) < 1.e-6f);
  check &= (fabsf(A[1][1] - 1.0f) < 1.e-6f);
#if defined(HYDRO_DIMENSION_3D)
  check &= (fabsf(A[0][2]) < 1.e-6f);
  check &= (fabsf(A[1][2]) < 1.e-6f);
  check &= (fabsf(A[2][0]) < 1.e-6f);
  check &= (fabsf(A[2][1]) < 1.e-6f);
  check &= (fabsf(A[2][2] - 1.0f) < 1.e-6f);
#endif  // 3D
#endif  // 2D and 3D

  return check;
}

void print_matrix(float A[3][3], const char* s) {
  message("Matrix %s:", s);
#if defined(HYDRO_DIMENSION_1D)
  message("[%.3e]", A[0][0]);
#elif defined(HYDRO_DIMENSION_2D)
  message("[%.3e, %.3e]", A[0][0], A[0][1]);
  message("[%.3e, %.3e]", A[1][0], A[1][1]);
#elif defined(HYDRO_DIMENSION_3D)
  message("[%.8e, %.8e, %.8e]", A[0][0], A[0][1], A[0][2]);
  message("[%.8e, %.8e, %.8e]", A[1][0], A[1][1], A[1][2]);
  message("[%.8e, %.8e, %.8e]", A[2][0], A[2][1], A[2][2]);
#endif
}

void multiply_matrices(float A[3][3], float B[3][3], float C[3][3]) {
#if defined(HYDRO_DIMENSION_1D)
  C[0][0] = A[0][0] * B[0][0];
#elif defined(HYDRO_DIMENSION_2D)
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      C[i][j] = 0.0f;
      for (int k = 0; k < 2; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
#elif defined(HYDRO_DIMENSION_3D)
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      C[i][j] = 0.0f;
      for (int k = 0; k < 3; ++k) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
#endif
}

int main(int argc, char* argv[]) {

  float A[3][3], B[3][3], C[3][3];
  setup_matrix(A);

  memcpy(B, A, 9 * sizeof(float));

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (A[i][j] != B[i][j]) {
        error("Matrices not equal after copy!");
      }
    }
  }

  invert_dimension_by_dimension_matrix(A);

  multiply_matrices(A, B, C);

  if (!is_unit_matrix(C)) {
    print_matrix(A, "A");
    print_matrix(B, "B");
    print_matrix(C, "C");
    error("Inverted matrix is wrong!");
  }

  return 0;
}
