/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEBUG_INTERACTIONS_HYDRO_PART_H
#define SWIFT_DEBUG_INTERACTIONS_HYDRO_PART_H

/**
 * @file DebugInteractions/hydro_part.h
 * @brief Empty SPH implementation used solely to test the SELF/PAIR routines.
 */

#include "cooling_struct.h"

/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /* Velocity at the last full step. */
  float v_full[3];

  /* Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

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

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density;

  struct {

    float h_dt;

  } force;

  /* Time-step length */
  timebin_t time_bin;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

  /*! List of interacting particles in the density SELF and PAIR */
  long long ids_ngbs_density[256];

  /*! List of interacting particles in the force SELF and PAIR */
  long long ids_ngbs_force[256];

  /*! Number of interactions in the density SELF and PAIR */
  int num_ngb_density;

  /*! Number of interactions in the force SELF and PAIR */
  int num_ngb_force;

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_DEBUG_INTERACTIONS_HYDRO_PART_H */
