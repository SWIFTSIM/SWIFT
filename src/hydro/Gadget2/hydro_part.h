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

/* Extra particle data not needed during the computation. */
struct xpart {

  /* Old position, at last tree rebuild. */
  double x_old[3];

  /* Velocity at the last full step. */
  float v_full[3];

} __attribute__((aligned(xpart_align)));

/* Data of a single particle. */
struct part {

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle time of beginning of time-step. */
  float t_begin;

  /* Particle time of end of time-step. */
  float t_end;

  struct {

    /* Number of neighbours */
    float wcount;

    /* Number of neighbours spatial derivative */
    float wcount_dh;

  } density;

  struct {

    /* Time derivative of the smoothing length */
    float h_dt;

  } force;

  /* Particle entropy. */
  float entropy;

  /* Particle density. */
  float rho;

  /* Derivative of the density with respect to this particle's smoothing length.
   */
  float rho_dh;

  /* Particle mass. */
  float mass;

  /* Particle pressure */
  float pressure;

  /* Entropy time derivative */
  float entropy_dt;

  /* Velocity divergence */
  float div_v;

  /* Velocity curl */
  float curl_v;
  float rot_v[3];

  /* Signal velocity */
  float v_sig;

  /* Particle ID. */
  unsigned long long id;

  /* Pointer to corresponding gravity part. */
  struct gpart* gpart;

} __attribute__((aligned(part_align)));
