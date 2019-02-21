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
#ifndef SWIFT_GADGET2_HYDRO_PART_H
#define SWIFT_GADGET2_HYDRO_PART_H

/**
 * @file Gadget2/hydro_part.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 * Springel, V., MNRAS, Volume 364, Issue 4, pp. 1105-1134.
 * We use the same numerical coefficients as the Gadget-2 code. When used with
 * the Spline-3 kernel, the results should be equivalent to the ones obtained
 * with Gadget-2 up to the rounding errors and interactions missed by the
 * Gadget-2 tree-code neighbours search.
 */

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "logger.h"
#include "star_formation_struct.h"
#include "tracers_struct.h"

/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /* Velocity at the last full step. */
  float v_full[3];

  /* Gravitational acceleration at the last full step. */
  float a_grav[3];

  /* Entropy at the last full step. */
  float entropy_full;

  /* Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

  /* Additional data used by the tracers */
  struct tracers_xpart_data tracers_data;

  /* Additional data used by the star formation */
  struct star_formation_xpart_data sf_data;

#ifdef WITH_LOGGER
  /* Additional data for the particle logger */
  struct logger_part_data logger_data;
#endif

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

  /* Particle mass. */
  float mass;

  /* Particle density. */
  float rho;

  /* Particle entropy. */
  float entropy;

  /* Entropy time derivative */
  float entropy_dt;

  union {

    struct {

      /* Number of neighbours. */
      float wcount;

      /* Number of neighbours spatial derivative. */
      float wcount_dh;

      /* Derivative of the density with respect to h. */
      float rho_dh;

      /* Particle velocity curl. */
      float rot_v[3];

      /* Particle velocity divergence. */
      float div_v;

    } density;

    struct {

      /* Balsara switch */
      float balsara;

      /*! "Grad h" term */
      float f;

      /* Pressure over density squared  */
      float P_over_rho2;

      /* Particle sound speed. */
      float soundspeed;

      /* Signal velocity. */
      float v_sig;

      /* Time derivative of the smoothing length */
      float h_dt;

    } force;
  };

  /* Chemistry information */
  struct chemistry_part_data chemistry_data;

  /* Time-step length */
  timebin_t time_bin;

  /* Need waking-up ? */
  timebin_t wakeup;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef DEBUG_INTERACTIONS_SPH
  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[MAX_NUM_OF_NEIGHBOURS];

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[MAX_NUM_OF_NEIGHBOURS];

  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;
#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_GADGET2_HYDRO_PART_H */
