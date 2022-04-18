/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "../config.h"
#include "align.h"
#include "kernel_hydro.h"
#include "vector.h"

#include <fenv.h>
#include <stdlib.h>
#include <strings.h>

const int numPoints = (1 << 28);

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const float h = 1.2348f;

  float *u, *W, *dW;
  if (posix_memalign((void **)&u, SWIFT_CACHE_ALIGNMENT,
                     numPoints * sizeof(float)) != 0)
    error("Error allocating u");
  if (posix_memalign((void **)&W, SWIFT_CACHE_ALIGNMENT,
                     numPoints * sizeof(float)) != 0)
    error("Error allocating W");
  if (posix_memalign((void **)&dW, SWIFT_CACHE_ALIGNMENT,
                     numPoints * sizeof(float)) != 0)
    error("Error allocating dW");

  message("Serial Output");
  message("-------------");
  const float numPoints_inv = 1. / numPoints;

  for (int i = 0; i < numPoints; ++i)
    u[i] = i * 1.2f * kernel_gamma * numPoints_inv / h;

  for (int i = 0; i < numPoints; ++i) {

    kernel_deval(u[i], &W[i], &dW[i]);

    if (W[i] < 0.f) error("Kernel is negative u=%e W=%e", u[i], W[i]);
    if (dW[i] > 0.f)
      error("Kernel derivatibe is positive u=%e dW=%e", u[i], dW[i]);
  }

  /* Test some additional special cases */
  float Wtest, dWtest;
  kernel_deval(1.930290, &Wtest, &dWtest);
  if (Wtest < 0.f) error("Kernel is negative u=%e W=%e", 1.930290, Wtest);
  if (dWtest > 0.f)
    error("Kernel derivative is positive u=%e dW=%e", 1.930290, dWtest);

#ifdef WITH_VECTORIZATION

  message("Vector Output for VEC_SIZE=%d", VEC_SIZE);
  message("-------------");

  message("Vector Output for kernel_deval_1_vec");
  message("-------------");

  /* Test vectorised kernel that uses one vector. */
  for (int i = 0; i < numPoints; i += VEC_SIZE) {

    vector vx, vx_h;
    vector W_vec, dW_vec;

    for (int j = 0; j < VEC_SIZE; j++)
      vx.f[j] = (i + j) * 1.2f * kernel_gamma / numPoints;

    vx_h.v = vec_mul(vx.v, vec_set1(1.f / h));

    kernel_deval_1_vec(&vx_h, &W_vec, &dW_vec);

    for (int j = 0; j < VEC_SIZE; j++) {

      /* message("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h, */
      /*        h * kernel_gamma, vx.f[j], W_vec.f[j], dW_vec.f[j]); */

      if (W_vec.f[j] < 0.f)
        error("Kernel is negative u=%e W=%e", u[i + j], W_vec.f[j]);
      if (dW_vec.f[j] > 0.f)
        error("Kernel derivative is positive u=%e dW=%e", u[i + j],
              dW_vec.f[j]);

      if (fabsf(W_vec.f[j] - W[i + j]) > 2e-6)
        error("Invalid Wvalue ! scalar= %e, vector= %e\n", W[i + j],
              W_vec.f[j]);
      if (fabsf(dW_vec.f[j] - dW[i + j]) > 2e-6)
        error("Invalid dW value ! scalar= %e, vector= %e %e %e\n", dW[i + j],
              dW_vec.f[j], fabsf(dW_vec.f[j] - dW[i + j]), fabsf(dW[i + j]));
    }
  }

  message("Vector Output for kernel_deval_2_vec");
  message("-------------");

  /* Test vectorised kernel that uses two vectors. */
  for (int i = 0; i < numPoints; i += VEC_SIZE) {

    vector vx, vx_h;
    vector W_vec, dW_vec;

    vector vx_2, vx_h_2;
    vector W_vec_2, dW_vec_2;

    for (int j = 0; j < VEC_SIZE; j++) {
      vx.f[j] = (i + j) * 1.2f * kernel_gamma / numPoints;
      vx_2.f[j] = (i + j) * 1.2f * kernel_gamma / numPoints;
    }

    vx_h.v = vec_mul(vx.v, vec_set1(1.f / h));
    vx_h_2.v = vec_mul(vx_2.v, vec_set1(1.f / h));

    kernel_deval_2_vec(&vx_h, &W_vec, &dW_vec, &vx_h_2, &W_vec_2, &dW_vec_2);

    /* Check first vector results. */
    for (int j = 0; j < VEC_SIZE; j++) {

      /* message("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h, */
      /*        h * kernel_gamma, vx.f[j], W_vec.f[j], dW_vec.f[j]); */

      if (W_vec.f[j] < 0.f)
        error("Kernel is negative u=%e W=%e", u[i + j], W_vec.f[j]);
      if (dW_vec.f[j] > 0.f)
        error("Kernel derivative is positive u=%e dW=%e", u[i + j],
              dW_vec.f[j]);

      if (fabsf(W_vec.f[j] - W[i + j]) > 2e-6)
        error("Invalid value ! scalar= %e, vector= %e\n", W[i + j], W_vec.f[j]);
      if (fabsf(dW_vec.f[j] - dW[i + j]) > 2e-6)
        error("Invalid value ! scalar= %e, vector= %e\n", dW[i + j],
              dW_vec.f[j]);
    }

    /* Check second vector results. */
    for (int j = 0; j < VEC_SIZE; j++) {

      /* message("%2d: h= %f H= %f x=%f W(x,h)=%f dW(x,h)=%f\n", i + j, h, */
      /* h * kernel_gamma, vx_2.f[j], W_vec_2.f[j], dW_vec_2.f[j]); */

      if (W_vec_2.f[j] < 0.f)
        error("Kernel is negative u=%e W=%e", u[i + j], W_vec_2.f[j]);
      if (dW_vec_2.f[j] > 0.f)
        error("Kernel derivative is positive u=%e dW=%e", u[i + j],
              dW_vec_2.f[j]);

      if (fabsf(W_vec_2.f[j] - W[i + j]) > 2e-6)
        error("Invalid value ! scalar= %e, vector= %e\n", W[i + j],
              W_vec_2.f[j]);
      if (fabsf(dW_vec_2.f[j] - dW[i + j]) > 2e-6)
        error("Invalid value ! scalar= %e, vector= %e\n", dW[i + j],
              dW_vec_2.f[j]);
    }
  }

  message("All values are consistent");

#endif

  free(u);
  free(W);
  free(dW);
  return 0;
}
