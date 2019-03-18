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
#ifndef SWIFT_DEFAULT_HYDRO_PART_H
#define SWIFT_DEFAULT_HYDRO_PART_H

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "star_formation_struct.h"
#include "tracers_struct.h"

/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /* Velocity at the last full step. */
  float v_full[3];

  /* Gravitational acceleration at the last full step. */
  float a_grav[3];

  /* Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

  /* Additional data used by the tracers */
  struct tracers_xpart_data tracers_data;

  /* Additional data used by the star formation */
  struct star_formation_xpart_data sf_data;

  float u_full;

  /* Old density. */
  float omega;

} SWIFT_STRUCT_ALIGN;

/* Data of a single particle. */
struct part {

  /* Particle ID. */
  long long id;

  /* Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle internal energy. */
  float u;

  /* Particle density. */
  float rho;

  /* Derivative of the density with respect to this particle's smoothing length.
   */
  float rho_dh;

  /* Particle viscosity parameter */
  float alpha;

  /* Store density/force specific stuff. */
  union {

    struct {

      /* Particle velocity divergence. */
      float div_v;

      /* Derivative of particle number density. */
      float wcount_dh;

      /* Particle velocity curl. */
      float rot_v[3];

      /* Particle number density. */
      float wcount;

    } density;

    struct {

      /* Balsara switch */
      float balsara;

      /* Aggregate quantities. */
      float P_over_rho2;

      /* Change in particle energy over time. */
      float u_dt;

      /* Signal velocity */
      float v_sig;

      /* Sound speed */
      float soundspeed;

      /* Change in smoothing length over time. */
      float h_dt;

    } force;
  };

  /* Particle mass. */
  float mass;

  /* Chemistry information */
  struct chemistry_part_data chemistry_data;

  /* Particle time-bin */
  timebin_t time_bin;

  /* Need waking-up ? */
  timebin_t wakeup;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DEFAULT_HYDRO_PART_H */
