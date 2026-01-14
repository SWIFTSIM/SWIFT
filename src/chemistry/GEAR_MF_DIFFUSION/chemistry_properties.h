/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_PROPERTIES_GEAR_MF_DIFFUSION_H
#define SWIFT_CHEMISTRY_PROPERTIES_GEAR_MF_DIFFUSION_H

#define GEAR_LABELS_SIZE 10  // redumndant with the one defined in

/* Cell slope limiter beta factors. The default values come from Hopkins (2015)
 */
#define GIZMO_SLOPE_LIMITER_BETA_MIN 1.0
#define GIZMO_SLOPE_LIMITER_BETA_MAX 2.0

/* Tolerance for the positivity preserving cell limiter. We allow a small
   overshoot for diffusion */
#define GEAR_FVPM_DIFFUSION_CELL_LIMITER_SHOOT_TOLERANGE 0.2

/* Tolerance for the Riemann solver wavespeed estimate to avoid large numbers
   when the wavespeeds are close */
#define GEAR_FVMP_DIFFUSION_WAVESPEED_ESTIMATE_DIFFERENCE_TOLERANCE 1e-8

/* Enable the flux limiters. You can turn off the ones you don't want by
   commenting these lines*/
#define GEAR_FVPM_DIFFUSION_FLUX_LIMITER_AGRESSIVE_RESCALING

/* Verbosity for the flux limiters */
#define GEAR_FVPM_DIFFUSION_FLUX_LIMITER_VERBOSITY 0

/**
 * @brief The diffusion mode
 */
enum chemistry_diffusion_mode {
  isotropic_constant,    /* Constant isotropic diffusion */
  isotropic_smagorinsky, /* Smagorinsky turbulent diffusion \propto |S| */
  anisotropic_gradient   /* Rennehan (2021) gradient model \propto S */
};

/**
 * @brief The relaxation time mode
 */
enum chemistry_relaxation_time_mode {
  constant_mode,  /* Constant */
  turbulent_mode, /* Based on the gas turbulences: tau \propto 1/ ||Shear|| */
};

/**
 * @brief Global chemical abundance information.
 */
struct chemistry_global_data {

  /* Initial mass fraction */
  double initial_metallicities[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* Solar mass abundances read from the chemistry table */
  float solar_abundances[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Name of the different elements */
  char elements_name[GEAR_CHEMISTRY_ELEMENT_COUNT * GEAR_LABELS_SIZE];

  /***************************************************************************/
  /* Parameter related to diffusion model */

  /*! Diffusion normalisation constant: \kappa \propto C */
  float diffusion_coefficient;

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  /*! Relaxation time for the constant relaxation time case. In physical units.
   */
  double tau;

  /*! Relaxation time mode. 0: constant, 1: . */
  enum chemistry_relaxation_time_mode relaxation_time_mode;
#endif

  /*! Diffusion mode. 0: isotropic with constant coefficient, 1: Smagorinsky
      isotrpoic diffusion, 2: anistropic diffusion with the shear tensor. */
  enum chemistry_diffusion_mode diffusion_mode;

  /***************************************************************************/
  /* HLL Riemann solver parameters
   * See Hopkins 2017 (https://arxiv.org/abs/1602.07703) */

  /*! The psi in eq (7) */
  float hll_riemann_solver_psi;

  /*! The epsilon in eq (14). This is a tolerance parameter. So, it must be 0
      <= epsilon <= 1. */
  float hll_riemann_solver_epsilon;

  /* CFL coefficient for integration on timesteps larger than the parabolic
     timestep */
  float C_CFL_chemistry;
};

#endif /* SWIFT_CHEMISTRY_PROPERTIES_GEAR_MF_DIFFUSION_H */
