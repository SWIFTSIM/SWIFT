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
#ifndef SWIFT_GPART_H
#define SWIFT_GPART_H

/* Some standard headers. */
#include <stdlib.h>


/* properties of external potential */
static struct ExternalPointMass
{
  const float Mass;
  const float Position[3];
} PointMass = {.Mass = 1, .Position={0.,0.,0.}};


/* Gravity particle. */
struct gpart {

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

  /* Particle external gravity acceleration */
  float a_grav_external[3];

  /* Particle mass. */
  float mass;

  /* Particle time of beginning of time-step. */
  float t_begin;

  /* Particle time of end of time-step. */
  float t_end;

  /* Anonymous union for id/part. */
  union {

    /* Particle ID. */
    size_t id;

    /* Pointer to corresponding SPH part. */
    struct part* part;
  };

} __attribute__((aligned(part_align)));
#endif /* SWIFT_GPART_H */
