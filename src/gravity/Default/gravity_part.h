/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
/* Some standard headers. */
#include <stdlib.h>

/* Gravity particle. */
struct gpart {

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v_full[3];

  /* Particle acceleration. */
  float a_grav[3];

  /* Particle mass. */
  float mass;

  /* Particle time of beginning of time-step. */
  int ti_begin;

  /* Particle time of end of time-step. */
  int ti_end;

  /* Anonymous union for id/part. */
  union {

    /* Particle ID. */
    long long id;

    /* Pointer to corresponding SPH part. */
    struct part* part;
  };

} __attribute__((aligned(gpart_align)));
