/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_IACT_GRAV_H
#define SWIFT_RUNNER_IACT_GRAV_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "vector.h"

/**
 * @brief Gravity potential
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav(
    float r2, float *dx, struct gpart *pi, struct gpart *pj) {

  float ir, r;
  float w, acc;
  float mi = pi->mass, mj = pj->mass;
  int k;

  /* Get the absolute distance. */
  ir = 1.0f / sqrtf(r2);
  r = r2 * ir;

  /* Evaluate the gravity kernel. */
  kernel_grav_eval(r, &acc);

  /* Scale the acceleration. */
  acc *= const_G * ir * ir * ir;

  /* Aggregate the accelerations. */
  for (k = 0; k < 3; k++) {
    w = acc * dx[k];
    pi->a_grav[k] -= w * mj;
    pj->a_grav[k] += w * mi;
  }
}

/**
 * @brief Gravity potential (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_grav(
    float *R2, float *Dx, struct gpart **pi, struct gpart **pj) {

#ifdef VECTORIZE

  vector ir, r, r2, dx[3];
  vector w, acc, ai, aj;
  vector mi, mj;
  int j, k;

#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ir.v = vec_rsqrt(r2.v);
  ir.v = ir.v - vec_set1(0.5f) * ir.v * (r2.v * ir.v * ir.v - vec_set1(1.0f));
  r.v = r2.v * ir.v;

  /* Evaluate the gravity kernel. */
  blender_eval_vec(&r, &acc);

  /* Scale the acceleration. */
  acc.v *= vec_set1(const_G) * ir.v * ir.v * ir.v;

  /* Aggregate the accelerations. */
  for (k = 0; k < 3; k++) {
    w.v = acc.v * dx[k].v;
    ai.v = w.v * mj.v;
    aj.v = w.v * mi.v;
    for (j = 0; j < VEC_SIZE; j++) {
      pi[j]->a_grav[k] -= ai.f[j];
      pj[j]->a_grav[k] += aj.f[j];
    }
  }

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_grav(R2[k], &Dx[3 * k], pi[k], pj[k]);

#endif
}

#endif /* SWIFT_RUNNER_IACT_GRAV_H */
