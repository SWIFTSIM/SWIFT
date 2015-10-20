/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durhm.ac.uk)
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
#ifndef SWIFT_PART_H
#define SWIFT_PART_H

/* Some constants. */
#define part_maxwait 3
#define part_maxunlock 39
#define part_dtmax 10
#define part_align 64

/* Extra particle data not needed during the computation. */
struct xpart {

  /* Old position, at last tree rebuild. */
  double x_old[3];

  /* Velocity at the half-step. */
  float v_hdt[3];

  /* Entropy at the half-step. */
  float u_hdt;

  /* Old density. */
  float omega;

  /* particle's current time-step. */
  float dt_curr;

} __attribute__((aligned(32)));

/* Gravity particle. */
struct gpart {

  /* Particle position. */
  double x[3];

  /* Particle velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

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
    struct part *part;
  };

} __attribute__((aligned(part_align)));

/* Data of a single particle. */
struct part {

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle time of beginning of time-step. */
  float t_begin;

  /* Particle time of end of time-step. */
  float t_end;

  /* Particle internal energy. */
  float u;

  /* Particle density. */
  float rho;

  /* Derivative of the density with respect to this particle's smoothing length.
   */
  float rho_dh;

#ifndef LEGACY_GADGET2_SPH
  /* Particle viscosity parameter */
  float alpha;
#endif

  /* Store density/force specific stuff. */
  union {

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

    struct {

      /* Balsara switch */
      float balsara;

      /* Aggregate quantities. */
      float POrho2;

      /* Change in particle energy over time. */
      float u_dt;

      /* Change in smoothing length over time. */
      float h_dt;

      /* Signal velocity */
      float v_sig;

      /* Sound speed */
      float c;

    } force;
  };

  /* Particle pressure. */
  // float P;

  /* Particle mass. */
  float mass;

  /* Particle ID. */
  unsigned long long id;

  /* Associated gravitas. */
  struct gpart *gpart;

} __attribute__((aligned(part_align)));

#endif /* SWIFT_PART_H */
