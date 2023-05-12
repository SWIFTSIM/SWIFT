/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SHADOWSWIFT_HYDRO_PART_H
#define SWIFT_SHADOWSWIFT_HYDRO_PART_H

/**
 * @file Shadowswift/hydro_part.h
 * @brief Moving mesh hydrodynamics implementation
 */

#include "black_holes_struct.h"
#include "chemistry_struct.h"
#include "const.h"
#include "cooling_struct.h"
#include "feedback_struct.h"
#include "particle_splitting_struct.h"
#include "rt_struct.h"
#include "sink_struct.h"
#include "star_formation_struct.h"
#include "timestep_limiter_struct.h"
#include "tracers_struct.h"

enum kick_type {
  KICK1,
  KICK2,
  ROLLBACK,
  RESTORE_AFTER_ROLLBACK,
};

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

  /*! Gravitational acceleration at the end of the last step */
  float a_grav[3];

  /*! Internal energy at the last full step. */
  float u_full;

  /*! Additional data used to record particle splits */
  struct particle_splitting_data split_data;

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

  /* Additional data used by the tracers */
  struct tracers_xpart_data tracers_data;

  /* Additional data used by the tracers */
  struct star_formation_xpart_data sf_data;

  /* Additional data used by the feedback */
  struct feedback_xpart_data feedback_data;

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Particle fields for the Moving mesh particles
 *
 */
struct part {

  /*! Particle unique ID. */
  long long id;

  /*! Pointer to corresponding gravity part. */
  struct gpart* gpart;

  /*! Particle position. */
  double x[3];

  /*! Particle acceleration. */
  float a_hydro[3];

  /*! Particle smoothing length. We use this as the search radius for completing
   * the Voronoi mesh. */
  float h;

  /* Density. */
  float rho;

  /* Fluid velocity. */
  float v[3];

  /* Particle velocity */
  float v_full[3];

  /* Pressure. */
  float P;

  /* Entropic function */
  float A;

  /* Fluid thermal energy (not per unit mass!). */
  float thermal_energy;

  /* Gradients of the primitive variables. */
  struct {

    /* Density gradients. */
    float rho[3];

    /* Fluid velocity gradients. */
    float v[3][3];

    /* Pressure gradients. */
    float P[3];

    /* Entropic function gradients */
    float A[3];

#ifdef SHADOWSWIFT_GRADIENTS_WLS
    /* Matrix to invert during gradient calculation */
    float matrix_wls[3][3];
#endif

  } gradients;

  /* Quantities needed by the slope limiter. */
  struct {

    /* Extreme values of the density among the neighbours. */
    float rho[2];

    /* Extreme values of the fluid velocity among the neighbours. */
    float v[3][2];

    /* Extreme values of the pressure among the neighbours. */
    float P[2];

    /* Extreme values of the entropic function among the neighbours */
    float A[2];

    /* Maximal value kinetic energy among the neighbours */
    float Ekin;

    struct {

      /* Extreme values of the extrapolated density towards the neighbours. */
      float rho[2];

      /* Extreme values of the extrapolated fluid velocity towards the
       * neighbours. */
      float v[3][2];

      /* Extreme values of the extrapolated pressure towards the neighbours. */
      float P[2];

      /* Extreme values of the extrapolated entropic function towards the
       * neighbours.*/
      float A[2];

    } extrapolations;

  } limiter;

#ifdef SHADOWSWIFT_EXTRAPOLATE_TIME
  /* Time extrapolations of primitive variables (cumulative over timestep) */
  float dW_time[6];
#endif

  /* The conserved hydrodynamical variables. */
  struct {

    /* Fluid mass */
    float mass;

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid total energy. */
    float energy;

    /* This is related to the thermodynamical entropy (see Springel 2010 eq.
     * 49).
     * Note that this is in general *not* conserved, but conservation of
     * entropy may be given precedence in smooth cold flows */
    float entropy;

  } conserved;

  /* Flux counter, should be conserved */
  long flux_count;

  /* Fluxes. */
  struct {

    /* Mass flux. */
    float mass;

    /* Momentum flux. */
    float momentum[3];

    /* Energy flux. */
    float energy;

    /* Entropy flux. */
    float entropy;

    /* Timestep for flux calculation. */
    float dt;

  } flux;

  struct {

    /*! Voronoi cell volume. */
    float volume;

    /*! Voronoi cell centroid, relative to this particles position. */
    float centroid[3];

    /*! Number of faces of this voronoi cell. */
    int nface;

    /*! cell_pair_connections offset in voronoi tesselation */
    int pair_connections_offset;

    /*! Estimate of the minimal distance from centroid to a face */
    float min_face_dist;

    /*! Flags indicating to which neighbouring cells this particle has already
     * been added. */
    int delaunay_flags;

  } geometry;

  struct {

    /* Signal velocity */
    float vmax;

    /* Last kick applied to this particle. */
    enum kick_type last_kick;

  } timestepvars;

  /* Unused in the ShadowSWIFT scheme */
  struct {

    /* Derivative of particle number density. */
    float wcount_dh;

    /* Particle number density. */
    float wcount;

  } density;

  /* Unused in the ShadowSWIFT scheme. */
  struct {

    /* Needed to drift the primitive variables. */
    float h_dt;

  } force;

  /* Specific stuff for the gravity-hydro coupling. */
  struct {

    /* Current value of the mass flux vector. */
    float mflux[3];

  } gravity;

  /*! Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Cooling information */
  struct cooling_part_data cooling_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data feedback_data;

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /*! Sink information (e.g. swallowing ID) */
  struct sink_part_data sink_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;

  /*! RT sub-cycling time stepping data */
  struct rt_timestepping_data rt_time_data;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_SHADOWSWIFT_HYDRO_PART_H */
