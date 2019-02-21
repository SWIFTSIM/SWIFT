/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "tracers_struct.h"
#include "voronoi_cell.h"

/* Extra particle data not needed during the computation. */
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

  /* Particle cutoff radius. */
  float h;

  /* The primitive hydrodynamical variables. */
  struct {

    /* Fluid velocity. */
    float v[3];

    /* Density. */
    float rho;

    /* Pressure. */
    float P;

    /* Gradients of the primitive variables. */
    struct {

      /* Density gradients. */
      float rho[3];

      /* Fluid velocity gradients. */
      float v[3][3];

      /* Pressure gradients. */
      float P[3];

    } gradients;

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

  } primitives;

  /* The conserved hydrodynamical variables. */
  struct {

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid mass (this field already exists outside of this struct as well). */
    float mass;

    /* Fluid thermal energy (not per unit mass!). */
    float energy;

    /* Fluxes. */
    struct {

      /* Mass flux. */
      float mass;

      /* Momentum flux. */
      float momentum[3];

      /* Energy flux. */
      float energy;

    } flux;

  } conserved;

  /* Variables used for timestep calculation (currently not used). */
  struct {

    /* Maximum fluid velocity among all neighbours. */
    float vmax;

  } timestepvars;

  /* Quantities used during the volume (=density) loop. */
  struct {

    /* Derivative of particle number density. */
    float wcount_dh;

    /* Particle number density. */
    float wcount;

  } density;

  /* Quantities used during the force loop. */
  struct {

    /* Needed to drift the primitive variables. */
    float h_dt;

    /* Physical time step of the particle. */
    float dt;

    /* Active flag. */
    char active;

    /* Actual velocity of the particle. */
    float v_full[3];

  } force;

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

  /* Voronoi cell. */
  struct voronoi_cell cell;

} SWIFT_STRUCT_ALIGN;
