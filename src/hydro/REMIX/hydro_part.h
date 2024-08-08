/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_HYDRO_PART_H
#define SWIFT_PLANETARY_HYDRO_PART_H

/**
 * @file Planetary/hydro_part.h
 * @brief REMIX implementation of SPH (Sandnes et al. 2024)
 */

#include "black_holes_struct.h"
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "equation_of_state.h"  // For enum material_id
#include "feedback_struct.h"
#include "mhd_struct.h"
#include "particle_splitting_struct.h"
#include "rt_struct.h"
#include "sink_struct.h"
#include "star_formation_struct.h"
#include "symmetric_matrix.h"
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

  /*! Evolved density at the last full step. */
  float rho_evol_full;

  /*! Additional data used to record particle splits */
  struct particle_splitting_data split_data;

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data cooling_data;

  /*! Additional data used by the tracers */
  struct tracers_xpart_data tracers_data;

  /*! Additional data used by the star formation */
  struct star_formation_xpart_data sf_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data feedback_data;

  /*! Additional data used by the MHD scheme */
  struct mhd_xpart_data mhd_data;

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

  /*! Particle density (standard SPH estimate). */
  float rho;

  /*! Particle evolved density (primary density for REMIX equations). */
  float rho_evol;

  /*! Time derivative of the evolved density. */
  float drho_dt;

  /*! Gradient of velocity, calculated using normalised kernel. */
  float dv_norm_kernel[3][3];

  /*! Gradient of internal energy, calculated using normalised kernel. */
  float du_norm_kernel[3];

  /*! Gradient of density, calculated using normalised kernel. */
  float drho_norm_kernel[3];

  /*! Gradient of smoothing length, calculated using normalised kernel. */
  float dh_norm_kernel[3];

  /*! Geometric moment m_0. */
  float m0;

  /*! Gradient of geometric moment m_0. */
  float grad_m0[3];

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
     * @brief Structure for the variables only used in the gradient loop over
     * neighbours.
     *
     * Quantities should only be accessed in the gradient loop over neighbours
     * and the extra ghost task.
     */
    struct {

      /*! Symmetrised kernel geometric moment m0. */
      float m0_bar;

      /*! Gradient of symmetrised kernel geometric moment m0. */
      float grad_m0_bar[3];

      /*! Contributing term for grad-h component of grad_m0_bar. */
      float grad_m0_bar_gradhterm;

      /*! Symmetrised kernel geometric moment m1. */
      float m1_bar[3];

      /*! Gradient of symmetrised kernel geometric moment m1. */
      float grad_m1_bar[3][3];

      /*! Contributing term for grad-h component of grad_m1_bar. */
      float grad_m1_bar_gradhterm[3];

      /*! Symmetrised kernel geometric moment m2. */
      struct sym_matrix m2_bar;

      /*! Gradient of symmetrised kernel geometric moment m2 (x, y, z
       * components). */
      struct sym_matrix grad_m2_bar[3];

      /*! Contributing term for grad-h component of grad_m2_bar. */
      struct sym_matrix grad_m2_bar_gradhterm;

    } gradient;

    /**
     * @brief Structure for the variables only used in the force loop over
     * neighbours.
     *
     * Quantities in this sub-structure should only be accessed in the force
     * loop over neighbours and the ghost, drift and kick tasks.
     */
    struct {

      /*! Particle pressure. */
      float pressure;

      /*! Particle soundspeed. */
      float soundspeed;

      /*! Particle signal velocity */
      float v_sig;

      /*! Time derivative of smoothing length  */
      float h_dt;

      /*! Balsara switch */
      float balsara;

      /*! Linear-order reproducing kernel coefficient */
      float A;

      /*! Linear-order reproducing kernel coefficient */
      float B[3];

      /*! Gradient of linear-order reproducing kernel coefficient */
      float grad_A[3];

      /*! Gradient of linear-order reproducing kernel coefficient */
      float grad_B[3][3];

      /*! Variable switch to identify proximity to vacuum boundaries. */
      float vac_switch;

    } force;
  };

  /*! Additional data used by the MHD scheme */
  struct mhd_part_data mhd_data;

  /*! Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Cooling information */
  struct cooling_part_data cooling_data;

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /*! Sink information (e.g. swallowing ID) */
  struct sink_part_data sink_data;

  /*! Material identifier flag */
  enum eos_planetary_material_id mat_id;

  /*! Phase state flag */
  enum eos_phase_state phase_state;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data rt_data;

  /*! RT sub-cycling time stepping data */
  struct rt_timestepping_data rt_time_data;

  /*! Time-step length */
  timebin_t time_bin;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

  /* Whether or not the particle has h=h_max ('1' or '0') */
  char is_h_max;

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

#ifdef PLANETARY_FIXED_ENTROPY
  /* Fixed specific entropy */
  float s_fixed;
#endif

#ifdef MATERIAL_STRENGTH
  // Stress tensor
  struct sym_matrix stress_tensor;

  // Principal stresses (Eigenvalues of stress_tensor)
  float principal_stress_eigen[3];

  // Deviatoric stress tensor
  struct sym_matrix deviatoric_stress_tensor;

  // Time derivative of deviatoric stress tensor
  struct sym_matrix dS_dt;

  // Gradient of velocity, calculated using linear-order reproducing kernel.
  float dv_lin_repr_kernel[3][3];
#endif /* MATERIAL_STRENGTH */

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_PLANETARY_HYDRO_PART_H */
