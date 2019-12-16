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
#ifndef SWIFT_EULER_HYDRO_PART_H
#define SWIFT_EULER_HYDRO_PART_H

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "star_formation_struct.h"
#include "timestep_limiter_struct.h"
#include "tracers_struct.h"
#include "black_holes_struct.h"
#include "density_diffusion_part_data.h"
#include "particle_shifting_part_data.h"


/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /* Velocity at the last full step. */
  float v_full[3];

  float u_full;


} SWIFT_STRUCT_ALIGN;

/* Data of a single particle. */
struct part {

  /* Particle ID. */
  long long id;

  struct gpart *gpart;

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];
  float v_minus1[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle constant acceleration (e.g. input acceleration) */
  float a_constant[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle internal energy. */
  float u;

  /* Particle density. */
  double /*float*/ rho;
 
  float rho_base;

  /* Particle density at the previous timestep. */
  double /*float*/ rho_t_minus1;

  /* Derivative of the density with respect to time */
  /*float*/ double drho_dt;

  /* Kimenatic Viscosity of the fluid particle */
  float viscosity;

  /* Particle pressure. */
  float pressure;

  /* Particle soundspeed. */
  float soundspeed;

  /* Particle mass. */
  float mass;

  /* Tau stress tensor values */
  float tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;

  /* Velocity gradient values */
  float grad_v_xx, grad_v_xy, grad_v_xz, grad_v_yy, grad_v_yz, grad_v_zz;

  /* Is the particle a boundary particle. */
  int is_boundary;

  /* Timesteps since last euler step.*/
  int since_euler;

  /* Max viscosity value needed for CFL condition. */
  float max_visc;

  struct density_diffusion_part_data density_diffusive_data;
  /* Extra data only required for the particle shifting. */
  struct particle_shifting_part_data particle_shifting_data;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

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
