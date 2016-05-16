/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "multipole.h"

/**
 * @brief Add the second multipole to the first one.
 *
 * @param m_sum The #multipole which will contain the sum.
 * @param m_other The #multipole to add.
 */

void multipole_add(struct multipole *m_sum, const struct multipole *m_other) {

#if multipole_order != 1
#error "Multipoles of order >1 not yet implemented."
#endif

  /* Correct the position. */
  const float M_tot = m_sum->mass + m_other->mass;
  const float M_tot_inv = 1.f / M_tot;
  for (int k = 0; k < 3; k++)
    m_sum->CoM[k] =
        (m_sum->CoM[k] * m_sum->mass + m_other->CoM[k] * m_other->mass) *
        M_tot_inv;

  /* Add the particle to the moments. */
  m_sum->mass = M_tot;
}

/**
 * @brief Add a particle to the given multipole.
 *
 * @param m The #multipole.
 * @param p The #gpart.
 */

void multipole_addpart(struct multipole *m, struct gpart *p) {

  /* #if multipole_order == 1 */

  /*   /\* Correct the position. *\/ */
  /*   float mm = m->coeffs[0], mp = p->mass; */
  /*   float w = 1.0f / (mm + mp); */
  /*   for (int k = 0; k < 3; k++) m->x[k] = (m->x[k] * mm + p->x[k] * mp) * w;
   */

  /*   /\* Add the particle to the moments. *\/ */
  /*   m->coeffs[0] = mm + mp; */

  /* #else */
  /* #error( "Multipoles of order %i not yet implemented." , multipole_order )
   */
  /* #endif */
}

/**
 * @brief Add a group of particles to the given multipole.
 *
 * @param m The #multipole.
 * @param p The #gpart array.
 * @param N Number of parts to add.
 */

void multipole_addparts(struct multipole *m, struct gpart *p, int N) {

  /* #if multipole_order == 1 */

  /*   /\* Get the combined mass and positions. *\/ */
  /*   double xp[3] = {0.0, 0.0, 0.0}; */
  /*   float mp = 0.0f, w; */
  /*   for (int k = 0; k < N; k++) { */
  /*     w = p[k].mass; */
  /*     mp += w; */
  /*     xp[0] += p[k].x[0] * w; */
  /*     xp[1] += p[k].x[1] * w; */
  /*     xp[2] += p[k].x[2] * w; */
  /*   } */

  /*   /\* Correct the position. *\/ */
  /*   float mm = m->coeffs[0]; */
  /*   w = 1.0f / (mm + mp); */
  /*   for (int k = 0; k < 3; k++) m->x[k] = (m->x[k] * mm + xp[k]) * w; */

  /*   /\* Add the particle to the moments. *\/ */
  /*   m->coeffs[0] = mm + mp; */

  /* #else */
  /* #error( "Multipoles of order %i not yet implemented." , multipole_order )
   */
  /* #endif */
}

/**
* @brief Reset the data of a #multipole.
*
* @param m The #multipole.
*/
void multipole_reset(struct multipole *m) {

  /* Just bzero the struct. */
  bzero(m, sizeof(struct multipole));
}

/**
* @brief Init a multipole from a set of particles.
*
* @param m The #multipole.
* @param parts The #gpart.
* @param N The number of particles.
*/
void multipole_init(struct multipole *m, const struct gpart *gparts,
                    int gcount) {

#if multipole_order != 1
#error "Multipoles of order >1 not yet implemented."
#endif

  /* Zero everything */
  multipole_reset(m);

  float mass = 0.0f;
  double x[3] = {0.0, 0.0, 0.0};

  /* Collect the particle data. */
  for (int k = 0; k < gcount; k++) {
    const float w = gparts[k].mass;
    mass += w;
    x[0] += gparts[k].x[0] * w;
    x[1] += gparts[k].x[1] * w;
    x[2] += gparts[k].x[2] * w;
  }

  /* Store the data on the multipole. */
  m->mass = mass;
  m->CoM[0] = x[0] / mass;
  m->CoM[1] = x[1] / mass;
  m->CoM[2] = x[2] / mass;
}
