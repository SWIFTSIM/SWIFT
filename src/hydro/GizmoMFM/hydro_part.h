/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2014 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#ifndef SWIFT_GIZMO_MFM_HYDRO_PART_H
#define SWIFT_GIZMO_MFM_HYDRO_PART_H

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "star_formation_struct.h"
#include "tracers_struct.h"

/* Extra particle data not needed during the computation. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between the current position and position at the last sort. */
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

} SWIFT_STRUCT_ALIGN;

/* Data of a single particle. */
struct part {

  /* Particle ID. */
  long long id;

  /* Associated gravitas. */
  struct gpart *gpart;

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle smoothing length. */
  float h;

  /* Density. */
  float rho;

  /* Pressure. */
  float P;

  union {
    /* Quantities used during the volume (=density) loop. */
    struct {

      /* Derivative of particle number density. */
      float wcount_dh;

      /* Particle number density. */
      float wcount;

    } density;

    /* Quantities needed by the slope limiter. */
    struct {

      /* Extreme values of the density among the neighbours. */
      float rho[2];

      /* Extreme values of the fluid velocity among the neighbours. */
      float v[3][2];

      /* Extreme values of the pressure among the neighbours. */
      float P[2];

      /* Maximal distance to all neighbouring faces. */
      float maxr;

    } limiter;

    struct {
      /* Fluxes. */
      struct {

        /* No mass flux, since it is always zero. */

        /* Momentum flux. */
        float momentum[3];

        /* Energy flux. */
        float energy;

      } flux;

      /* Variables used for timestep calculation. */
      struct {

        /* Maximum signal velocity among all the neighbours of the particle. The
         * signal velocity encodes information about the relative fluid
         * velocities
         * AND particle velocities of the neighbour and this particle, as well
         * as
         * the sound speed of both particles. */
        float vmax;

      } timestepvars;

      /* Quantities used during the force loop. */
      struct {

        /* Needed to drift the primitive variables. */
        float h_dt;

      } force;
    };
  };

  /* Gradients of the primitive variables. */
  struct {

    /* Density gradients. */
    float rho[3];

    /* Fluid velocity gradients. */
    float v[3][3];

    /* Pressure gradients. */
    float P[3];

  } gradients;

  /* The conserved hydrodynamical variables. */
  struct {

    /* Fluid mass */
    float mass;

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid thermal energy (not per unit mass!). */
    float energy;

  } conserved;

  /* Geometrical quantities used for hydro. */
  struct {

    /* Volume of the particle. */
    float volume;

    /* Geometrical shear matrix used to calculate second order accurate
       gradients */
    float matrix_E[3][3];

    /* Centroid of the "cell". */
    float centroid[3];

    /* Correction factor for wcount. */
    float wcorr;

  } geometry;

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

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_GIZMO_MFM_HYDRO_PART_H */
