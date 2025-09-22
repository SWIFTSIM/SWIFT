/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_MAGMA2_HYDRO_PART_H
#define SWIFT_MAGMA2_HYDRO_PART_H

/**
 * @file MAGMA2/hydro_part.h
 * @brief Density-Energy non-conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (particle definition)
 */

#include "adaptive_softening_struct.h"
#include "black_holes_struct.h"
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "csds.h"
#include "feedback_struct.h"
#include "hydro_parameters.h"
#include "mhd_struct.h"
#include "particle_splitting_struct.h"
#include "pressure_floor_struct.h"
#include "rt_struct.h"
#include "sink_struct.h"
#include "star_formation_struct.h"
#include "timestep_limiter_struct.h"
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

  /*! Additional data used by the MHD scheme */
  struct mhd_xpart_data mhd_data;

#ifdef WITH_CSDS
  /* Additional data for the particle csds */
  struct csds_part_data csds_data;
#endif

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

  /*! Minimum smoothing length in the kernel */
  float h_min;

  /*! Minimum time-step amongst neighbours */
  float dt_min;
 
  /*! Conduction du/dt */
  float u_dt_cond;
 
  #ifdef MAGMA2_DEBUG_CHECKS
  struct {
    /*! Correction matrix at the last time it was ill-conditioned */
    hydro_real_t correction_matrix[3][3];

    /*! Velocity tensor at ill-condition time */
    hydro_real_t velocity_tensor_aux[3][3];

    /*! Velocity tensor norm ill-conditioned */
    hydro_real_t velocity_tensor_aux_norm[3][3];
    
    /*! u_aux tensor at ill-condition time */
    hydro_real_t u_aux[3];

    /*! u_aux_norm tensor ill-conditioned */
    hydro_real_t u_aux_norm[3];

    /*! Number of times correction_matrix was ill-conditioned */
    int C_ill_conditioned_count;

    /*! Number of times velocity_tensor_aux_norm was ill-conditioned */
    int D_ill_conditioned_count;

    /*! Number of times u_aux_norm was ill-conditioned */
    int u_ill_conditioned_count;

    /*! How many low-order SPH gradients in force interactions */
    int N_force_low_order_grad;

    /*! How many high-order SPH gradients in force interactions */
    int N_force_high_order_grad;

    /*! Number of neighbors in the kernel */
    int num_ngb;

    /*! The maximum viscous signal speed this step */
    hydro_real_t v_sig_visc_max;

    /*! The maximum conductive signal speed this step */
    hydro_real_t v_sig_cond_max;
  } debug;
#endif

  /* Store gradients in a separate struct */
  struct {
#ifdef hydro_props_use_adiabatic_correction
    /*! Adiabatic kernel correction factor numerator */
    hydro_real_t adiabatic_f_numerator;

    /*! Adiabatic kernel correction factor denominator */
    hydro_real_t adiabatic_f_denominator;
#endif

    /*! Sum of the kernel weights */
    hydro_real_t wcount;

    /*! Correction matrix (C^ki in Rosswog 2020) */
    hydro_real_t correction_matrix[3][3];

    /*! Flag for whether C is ill-conditioned */
    char C_well_conditioned;

    /*! Full velocity gradient tensor */
    hydro_real_t velocity_tensor[3][3];

    /*! Auxiliary full velocity gradient tensor */
    hydro_real_t velocity_tensor_aux[3][3];

    /*! Flag for whether D (velocity_tensor_aux) is ill-conditioned */
    char D_well_conditioned;

    /*! Normalization for computing velocity_tensor_aux */
    hydro_real_t velocity_tensor_aux_norm[3][3];

    /*! Hessian tensor */
    hydro_real_t velocity_hessian[3][3][3];

    /*! Internal energy gradient */
    hydro_real_t u[3];

    /*! Auxiliary internal energy gradient */
    hydro_real_t u_aux[3];

    /*! Normalization for computing u_aux */
    hydro_real_t u_aux_norm[3];
    
    /*! Flag for whether u_aux_norm is ill-conditioned */
    char u_well_conditioned;

    /*! Internal energy Hessian */
    hydro_real_t u_hessian[3][3];

    /*! Minimum delta u across kernel */
    hydro_real_t du_min;

    /*! Maximum delta u across kernel */
    hydro_real_t du_max;

    /*! Minimum dv across kernel */
    hydro_real_t dv_min[3];

    /*! Maximum dv across kernel */
    hydro_real_t dv_max[3];

    /*! Kernel size for slope limiting */
    hydro_real_t kernel_size;

    /*! Balsara limiter for divergences */
    hydro_real_t balsara;

  } gradients;

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

    } force;
  };

  /*! Additional data used for adaptive softening */
  struct adaptive_softening_part_data adaptive_softening_data;

  /*! Additional data used by the MHD scheme */
  struct mhd_part_data mhd_data;

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

  /*! Additional data used by the pressure floor */
  struct pressure_floor_part_data pressure_floor_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;

  /*! RT sub-cycling time stepping data */
  struct rt_timestepping_data rt_time_data;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char depth_h;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_MAGMA2_HYDRO_PART_H */
