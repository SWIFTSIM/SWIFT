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
#include <strings.h>

/* This object's header. */
#include "multipole.h"

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
* @param gparts The #gpart.
* @param gcount The number of particles.
*/
void multipole_init(struct multipole *m, const struct gpart *gparts,
                    int gcount) {

#if multipole_order > 2
#error "Multipoles of order >2 not yet implemented."
#endif

  /* Zero everything */
  multipole_reset(m);

  /* Temporary variables */
  double mass = 0.0;
  double com[3] = {0.0, 0.0, 0.0};

#if multipole_order >= 2
  double I_xx = 0.0, I_yy = 0.0, I_zz = 0.0;
  double I_xy = 0.0, I_xz = 0.0, I_yz = 0.0;
#endif

  /* Collect the particle data. */
  for (int k = 0; k < gcount; k++) {
    const float w = gparts[k].mass;

    mass += w;
    com[0] += gparts[k].x[0] * w;
    com[1] += gparts[k].x[1] * w;
    com[2] += gparts[k].x[2] * w;

#if multipole_order >= 2
    I_xx += gparts[k].x[0] * gparts[k].x[0] * w;
    I_yy += gparts[k].x[1] * gparts[k].x[1] * w;
    I_zz += gparts[k].x[2] * gparts[k].x[2] * w;
    I_xy += gparts[k].x[0] * gparts[k].x[1] * w;
    I_xz += gparts[k].x[0] * gparts[k].x[2] * w;
    I_yz += gparts[k].x[1] * gparts[k].x[2] * w;
#endif
  }

  const double imass = 1.0 / mass;

  /* Store the data on the multipole. */
  m->mass = mass;
  m->CoM[0] = com[0] * imass;
  m->CoM[1] = com[1] * imass;
  m->CoM[2] = com[2] * imass;

#if multipole_order >= 2
  m->I_xx = I_xx - imass * com[0] * com[0];
  m->I_yy = I_yy - imass * com[1] * com[1];
  m->I_zz = I_zz - imass * com[2] * com[2];
  m->I_xy = I_xy - imass * com[0] * com[1];
  m->I_xz = I_xz - imass * com[0] * com[2];
  m->I_yz = I_yz - imass * com[1] * com[2];
#endif
}

/**
 * @brief Add the second multipole to the first one.
 *
 * @param ma The #multipole which will contain the sum.
 * @param mb The #multipole to add.
 */

void multipole_add(struct multipole *ma, const struct multipole *mb) {

#if multipole_order > 2
#error "Multipoles of order >2 not yet implemented."
#endif

  /* Correct the position. */
  const double ma_mass = ma->mass;
  const double mb_mass = mb->mass;
  const double M_tot = ma_mass + mb_mass;
  const double M_tot_inv = 1.0 / M_tot;

  const double ma_CoM[3] = {ma->CoM[0], ma->CoM[1], ma->CoM[2]};
  const double mb_CoM[3] = {mb->CoM[0], mb->CoM[1], mb->CoM[2]};

#if multipole_order >= 2
  const double ma_I_xx = (double)ma->I_xx + ma_mass * ma_CoM[0] * ma_CoM[0];
  const double ma_I_yy = (double)ma->I_yy + ma_mass * ma_CoM[1] * ma_CoM[1];
  const double ma_I_zz = (double)ma->I_zz + ma_mass * ma_CoM[2] * ma_CoM[2];
  const double ma_I_xy = (double)ma->I_xy + ma_mass * ma_CoM[0] * ma_CoM[1];
  const double ma_I_xz = (double)ma->I_xz + ma_mass * ma_CoM[0] * ma_CoM[2];
  const double ma_I_yz = (double)ma->I_yz + ma_mass * ma_CoM[1] * ma_CoM[2];

  const double mb_I_xx = (double)mb->I_xx + mb_mass * mb_CoM[0] * mb_CoM[0];
  const double mb_I_yy = (double)mb->I_yy + mb_mass * mb_CoM[1] * mb_CoM[1];
  const double mb_I_zz = (double)mb->I_zz + mb_mass * mb_CoM[2] * mb_CoM[2];
  const double mb_I_xy = (double)mb->I_xy + mb_mass * mb_CoM[0] * mb_CoM[1];
  const double mb_I_xz = (double)mb->I_xz + mb_mass * mb_CoM[0] * mb_CoM[2];
  const double mb_I_yz = (double)mb->I_yz + mb_mass * mb_CoM[1] * mb_CoM[2];
#endif

  /* New mass */
  ma->mass = M_tot;

  /* New CoM */
  ma->CoM[0] = (ma_CoM[0] * ma_mass + mb_CoM[0] * mb_mass) * M_tot_inv;
  ma->CoM[1] = (ma_CoM[1] * ma_mass + mb_CoM[1] * mb_mass) * M_tot_inv;
  ma->CoM[2] = (ma_CoM[2] * ma_mass + mb_CoM[2] * mb_mass) * M_tot_inv;

/* New quadrupole */
#if multipole_order >= 2
  ma->I_xx = (ma_I_xx + mb_I_xx) - M_tot * ma->CoM[0] * ma->CoM[0];
  ma->I_yy = (ma_I_yy + mb_I_yy) - M_tot * ma->CoM[1] * ma->CoM[1];
  ma->I_zz = (ma_I_zz + mb_I_zz) - M_tot * ma->CoM[2] * ma->CoM[2];
  ma->I_xy = (ma_I_xy + mb_I_xy) - M_tot * ma->CoM[0] * ma->CoM[1];
  ma->I_xz = (ma_I_xz + mb_I_xz) - M_tot * ma->CoM[0] * ma->CoM[2];
  ma->I_yz = (ma_I_yz + mb_I_yz) - M_tot * ma->CoM[1] * ma->CoM[2];
#endif
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
