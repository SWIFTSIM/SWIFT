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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_PART_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_PART_H

/**
 * @file PressureEnergy/hydro_part.h
 * @brief Pressure-Energy implementation of SPH (Particle definition)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (16), (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

#include "chemistry_struct.h"
#include "cooling_struct.h"

/* Extra particle data not needed during the SPH loops over neighbours. */
struct xpart {

  /*! Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /*! Velocity at the last full step. */
  float v_full[3];

  /*! Gravitational acceleration */
  float a_grav[3];

  /*! Full internal energy */
  float u_full;

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

} SWIFT_STRUCT_ALIGN;

/* Data of a single particle. */
struct part {

  /*! Particle ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /*! Particle predicted velocity. */
  float v[3];

  /*! Particle acceleration. */
  float a_hydro[3];

  /*! Particle mass. */
  float mass;

  /*! SPH Density */
  float rho;

  /*! Smoothed particle pressure. */
  float pressure_bar;

  /*! Smoothed particle pressure's spatial derivative */
  float pressure_bar_dh;

  /*! Particle internal energy */
  float u;

  /*! Differential of the internal energy with respect to time */
  float u_dt;

  /*! Entropy (for if people want it) */ 
  float entropy;

  /*! Particle cutoff radius. */
  float h;


    struct {

      /*! Number of neighbours. */
      float wcount;

      /*! Number of neighbours spatial derivative. */
      float wcount_dh;

      /*! Velocity curl */
      float rot_v[3];

      /*! Velocity divergence */
      float div_v;

      /*! d\rho/dh */
      float rho_dh;

      /*! dP/dh */
      float pressure_dh;

    } density;

    struct {

      /*! Signal velocity. */
      float v_sig;

      /*! Time derivative of the smoothing length */
      float h_dt;

      /*! Sound speed */
      float soundspeed;
      
      /*! Balsara switch */
      float balsara;

      /*! F_{ij} -- not actually possible to set this. */
      float f;

      /*! P/\rho^2 -- not actually required for Pressure Energy but this is
      needed for cross-compatibility with the unit tests */
      float P_over_rho2;

    } force;

  /* Chemistry information */
  struct chemistry_part_data chemistry_data;

  /* Time-step length */
  timebin_t time_bin;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_PART_H */
