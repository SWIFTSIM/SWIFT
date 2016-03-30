/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_MULTIPOLE_H
#define SWIFT_MULTIPOLE_H

/* Some standard headers. */
#include <math.h>

/* Includes. */
#include "const.h"
#include "inline.h"
#include "kernel.h"
#include "part.h"

/* Some constants. */
#define multipole_order 1

/* Multipole struct. */
struct multipole {

  /* Multipole location. */
  double x[3];

  /* Acceleration on this multipole. */
  float a[3];

  /* Multipole coefficients. */
  float coeffs[multipole_order * multipole_order];
};

/* Multipole function prototypes. */
static void multipole_iact_mm(struct multipole *ma, struct multipole *mb,
                              double *shift);
void multipole_merge(struct multipole *ma, struct multipole *mb);
void multipole_addpart(struct multipole *m, struct gpart *p);
void multipole_addparts(struct multipole *m, struct gpart *p, int N);
void multipole_init(struct multipole *m, struct gpart *parts, int N);
void multipole_reset(struct multipole *m);

/**
 * @brief Compute the pairwise interaction between two multipoles.
 *
 * @param ma The first #multipole.
 * @param mb The second #multipole.
 * @param shift The periodicity correction.
 */

__attribute__((always_inline)) INLINE static void multipole_iact_mm(
    struct multipole *ma, struct multipole *mb, double *shift) {

  float dx[3], ir, r, r2 = 0.0f, acc;
  int k;

  /* Compute the multipole distance. */
  for (k = 0; k < 3; k++) {
    dx[k] = ma->x[k] - mb->x[k] - shift[k];
    r2 += dx[k] * dx[k];
  }

  /* Compute the normalized distance vector. */
  ir = 1.0f / sqrtf(r2);
  r = r2 * ir;

  /* Evaluate the gravity kernel. */
  kernel_grav_eval(r, &acc);

  /* Scale the acceleration. */
  acc *= const_G * ir * ir * ir;

/* Compute the forces on both multipoles. */
#if multipole_order == 1
  float mma = ma->coeffs[0], mmb = mb->coeffs[0];
  for (k = 0; k < 3; k++) {
    ma->a[k] -= dx[k] * acc * mmb;
    mb->a[k] += dx[k] * acc * mma;
  }
#else
#error( "Multipoles of order %i not yet implemented." , multipole_order )
#endif
}

/**
 * @brief Compute the interaction of a multipole on a particle.
 *
 * @param m The #multipole.
 * @param p The #gpart.
 * @param shift The periodicity correction.
 */

__attribute__((always_inline)) INLINE static void multipole_iact_mp(
    struct multipole *m, struct gpart *p, double *shift) {

  float dx[3], ir, r, r2 = 0.0f, acc;
  int k;

  /* Compute the multipole distance. */
  for (k = 0; k < 3; k++) {
    dx[k] = m->x[k] - p->x[k] - shift[k];
    r2 += dx[k] * dx[k];
  }

  /* Compute the normalized distance vector. */
  ir = 1.0f / sqrtf(r2);
  r = r2 * ir;

  /* Evaluate the gravity kernel. */
  kernel_grav_eval(r, &acc);

  /* Scale the acceleration. */
  acc *= const_G * ir * ir * ir * m->coeffs[0];

/* Compute the forces on both multipoles. */
#if multipole_order == 1
  for (k = 0; k < 3; k++) p->a_grav[k] += dx[k] * acc;
#else
#error( "Multipoles of order %i not yet implemented." , multipole_order )
#endif
}

#endif /* SWIFT_MULTIPOLE_H */
