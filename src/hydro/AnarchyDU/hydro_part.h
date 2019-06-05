/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_ANARCHY_DU_HYDRO_PART_H
#define SWIFT_ANARCHY_DU_HYDRO_PART_H

/**
 * @file AnarchyDU/hydro_part.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added ANARCHY physics (Cullen & Denhen 2011 AV,
 *        Price 2008 thermal diffusion (particle definition)
 */

#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "star_formation_struct.h"
#include "tracers_struct.h"

/**
 * @brief Particle fields not needed during the SPH loops over neighbours.
 *
 * This structure contains the particle fields that are not used in the
 * density or force loops. Quantities should be used in the kick, drift and
 * potentially ghost tasks only.
 */
struct xpart {

  /*! Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /*! Offset between the current position and position at the last sort. */
  float x_diff_sort[3];

  /*! Velocity at the last full step. */
  float v_full[3];

  /*! Gravitational acceleration at the last full step. */
  float a_grav[3];

  /*! Internal energy at the last full step. */
  float u_full;

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

  /* Additional data used by the tracers */
  struct tracers_xpart_data tracers_data;

  /* Additional data used by the tracers */
  struct star_formation_xpart_data sf_data;

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Particle fields for the SPH particles
 *
 * The density and force substructures are used to contain variables only used
 * within the density and force loops over neighbours. All more permanent
 * variables should be declared in the main part of the part structure,
 */
struct part {

  /*! Particle unique ID. */
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

  /*! Particle smoothing length. */
  float h;

  /*! Particle internal energy. */
  float u;

  /*! Time derivative of the internal energy. */
  float u_dt;

  /*! Particle density. */
  float rho;

  /* Store viscosity information in a separate struct. */
  struct {

    /*! Particle velocity divergence */
    float div_v;

    /*! Particle velocity divergence from previous step */
    float div_v_previous_step;

    /*! Artificial viscosity parameter */
    float alpha;

    /*! Signal velocity */
    float v_sig;

  } viscosity;

  /* Store thermal diffusion information in a separate struct. */
  struct {

    /*! del^2 u, a smoothed quantity */
    float laplace_u;

    /*! Thermal diffusion coefficient */
    float alpha;

  } diffusion;

  /* Store density/force specific stuff. */
  union {

    /**
     * @brief Structure for the variables only used in the density loop over
     * neighbours.
     *
     * Quantities in this sub-structure should only be accessed in the density
     * loop over neighbours and the ghost task.
     */
    struct {

      /*! Neighbour number count. */
      float wcount;

      /*! Derivative of the neighbour number with respect to h. */
      float wcount_dh;

      /*! Derivative of density with respect to h */
      float rho_dh;

      /*! Particle velocity curl. */
      float rot_v[3];

    } density;

    /**
     * @brief Structure for the variables only used in the force loop over
     * neighbours.
     *
     * Quantities in this sub-structure should only be accessed in the force
     * loop over neighbours and the ghost, drift and kick tasks.
     */
    struct {

      /*! "Grad h" term -- only partial in P-U */
      float f;

      /*! Particle pressure. */
      float pressure;

      /*! Particle soundspeed. */
      float soundspeed;

      /*! Time derivative of smoothing length  */
      float h_dt;

      /*! Balsara switch */
      float balsara;

    } force;
  };

  /* Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Time-step length */
  timebin_t time_bin;

  /* Need waking up ? */
  timebin_t wakeup;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_ANARCHY_DU_HYDRO_PART_H */
