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

/* Define the tolerable minimal negative mass for metals. Small negative masses
   can happen if the metal mass is close to 0. */
#define GEAR_NEGATIVE_METAL_MASS_FRACTION_TOLERANCE 0

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
 * @brief The Riemann solver type
 */
enum chemistry_riemann_solver {
  HLL,                       /* Regular HLL solver */
  HLL_parabolic_Hopkins2017, /* Hopkins (2017) HLL Riemann solver for
                                parabolic diffusion*/
  HLL_hyperbolic_Hopkins2017 /* Improved Hopkins (2017) HLL Riemann solver for
                                hyperbolic diffusion*/
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

  /* 1=Hopkins 2017, 2=HLL, 3=HLLC */
  int riemann_solver;
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
