/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

#define GFLOAT float

/* Extra particle data not needed during the computation. */
struct xpart {

  /* Old position, at last tree rebuild. */
  double x_old[3];

  /* Velocity at the last full step. */
  float v_full[3];

  /* Entropy at the half-step. */
  float u_hdt;

  /* Old density. */
  float omega;

} __attribute__((aligned(xpart_align)));

/* Data of a single particle. */
struct part {

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  float mass;  // MATTHIEU
  float h_dt;
  float rho;
  float rho_dh;

  /* Particle cutoff radius. */
  float h;

  /* Particle time of beginning of time-step. */
  int ti_begin;

  /* Particle time of end of time-step. */
  int ti_end;

  /* The primitive hydrodynamical variables */
  struct {

    /* fluid velocity */
    GFLOAT v[3];

    /* density */
    GFLOAT rho;

    /* pressure */
    GFLOAT P;

    struct {

      GFLOAT rho[3];

      GFLOAT v[3][3];

      GFLOAT P[3];

    } gradients;

    struct {

      /* extreme values among the neighbours */
      GFLOAT rho[2];

      GFLOAT v[3][2];

      GFLOAT P[2];

      /* maximal distance to all neighbouring faces */
      float maxr;

    } limiter;

  } primitives;

  /* The conserved hydrodynamical variables */
  struct {

    /* fluid momentum */
    GFLOAT momentum[3];

    /* fluid mass */
    GFLOAT mass;

    /* fluid energy */
    GFLOAT energy;

  } conserved;

  /* Geometrical quantities used for hydro */
  struct {

    /* volume of the particle */
    float volume;

    /* gradient matrix */
    float matrix_E[3][3];

  } geometry;

  struct {

    float vmax;

  } timestepvars;

  /* Quantities used during the density loop */
  struct {

    /* Particle velocity divergence. */
    float div_v;

    /* Derivative of particle number density. */
    float wcount_dh;

    /* Particle velocity curl. */
    float curl_v[3];

    /* Particle number density. */
    float wcount;

  } density;

  /* Particle ID. */
  unsigned long long id;

  /* Associated gravitas. */
  struct gpart *gpart;

} __attribute__((aligned(part_align)));
