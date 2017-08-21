/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    James Willis (james.s.willis@durham.ac.uk)
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

#include "kernel_hydro.h"
#include "vector.h"

#include <stdlib.h>
#include <strings.h>

#define numPoints (1 << 4)

int main() {

  const float h = 1.2348f;

  float u[numPoints] = {0.f};
  float W[numPoints] = {0.f};
  float dW[numPoints] = {0.f};

  printf("\nSerial Output\n");
  printf("-------------\n");
  const float numPoints_inv = 1. / numPoints;

  for (int i = 0; i < numPoints; ++i) {
    u[i] = i * 2.25f * numPoints_inv / h;
  }

  for (int i = 0; i < numPoints; ++i) {

    kernel_deval(u[i], &W[i], &dW[i]);

    printf("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i, h,
           h * kernel_gamma, u[i] * h, W[i], dW[i]);
  }

  printf("\nVector Output for VEC_SIZE=%d\n", VEC_SIZE);
  printf("-------------\n");

#ifdef WITH_VECTORIZATION

  printf("\nVector Output for kernel_deval_1_vec\n");
  printf("-------------\n");

  /* Test vectorised kernel that uses one vector. */
  for (int i = 0; i < numPoints; i += VEC_SIZE) {

    vector vx, vx_h;
    vector W_vec, dW_vec;

    for (int j = 0; j < VEC_SIZE; j++) {
      vx.f[j] = (i + j) * 2.25f / numPoints;
    }

    vx_h.v = vec_mul(vx.v, vec_set1(1.f / h));

    kernel_deval_1_vec(&vx_h, &W_vec, &dW_vec);

    for (int j = 0; j < VEC_SIZE; j++) {
      printf("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h,
             h * kernel_gamma, vx.f[j], W_vec.f[j], dW_vec.f[j]);

      if (fabsf(W_vec.f[j] - W[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", W[i + j],
               W_vec.f[j]);
        return 1;
      }
      if (fabsf(dW_vec.f[j] - dW[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", dW[i + j],
               dW_vec.f[j]);
        return 1;
      }
    }
  }

  printf("\nVector Output for kernel_deval_2_vec\n");
  printf("-------------\n");

  /* Test vectorised kernel that uses two vectors. */
  for (int i = 0; i < numPoints; i += VEC_SIZE) {

    vector vx, vx_h;
    vector W_vec, dW_vec;

    vector vx_2, vx_h_2;
    vector W_vec_2, dW_vec_2;

    for (int j = 0; j < VEC_SIZE; j++) {
      vx.f[j] = (i + j) * 2.25f / numPoints;
      vx_2.f[j] = (i + j) * 2.25f / numPoints;
    }

    vx_h.v = vec_mul(vx.v, vec_set1(1.f / h));
    vx_h_2.v = vec_mul(vx_2.v, vec_set1(1.f / h));

    kernel_deval_2_vec(&vx_h, &W_vec, &dW_vec, &vx_h_2, &W_vec_2, &dW_vec_2);

    /* Check first vector results. */
    for (int j = 0; j < VEC_SIZE; j++) {
      printf("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h,
             h * kernel_gamma, vx.f[j], W_vec.f[j], dW_vec.f[j]);

      if (fabsf(W_vec.f[j] - W[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", W[i + j],
               W_vec.f[j]);
        return 1;
      }
      if (fabsf(dW_vec.f[j] - dW[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", dW[i + j],
               dW_vec.f[j]);
        return 1;
      }
    }

    /* Check second vector results. */
    for (int j = 0; j < VEC_SIZE; j++) {
      printf("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h,
             h * kernel_gamma, vx_2.f[j], W_vec_2.f[j], dW_vec_2.f[j]);

      if (fabsf(W_vec_2.f[j] - W[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", W[i + j],
               W_vec_2.f[j]);
        return 1;
      }
      if (fabsf(dW_vec_2.f[j] - dW[i + j]) > 2e-7) {
        printf("Invalid value ! scalar= %e, vector= %e\n", dW[i + j],
               dW_vec_2.f[j]);
        return 1;
      }
    }
  }

  printf("\nAll values are consistent\n");

#endif
  return 0;
}
