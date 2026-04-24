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

#include "const.h"

/** @brief Size of the labels used for chemical species identifiers. */
/* Note: This is redundant with the one in
   src/feedback/GEAR/stellar_evolution_struct.h */
#define GEAR_LABELS_SIZE 10

/* --- Cell Slope Limiter (Hopkins 2015) --- */

/** @brief Minimum slope limiter factor (Barth-Jespersen-like).
 * Prevents the slope from being entirely zeroed out in smooth regions. */
#define GIZMO_SLOPE_LIMITER_BETA_MIN 1.0

/** @brief Maximum slope limiter factor.
 * Limits the interface reconstruction to prevent new local extrema. */
#define GIZMO_SLOPE_LIMITER_BETA_MAX 2.0

/* --- Positivity Preserving & Riemann Solvers --- */

/** @brief Tolerance for the cell-wide positivity limiter.
 * Allows a small fractional overshoot (20%) during diffusion to prevent
 * the limiter from being too aggressive and stalling the solution. */
#define GEAR_FVPM_DIFF_CELL_LIMITER_SHOOT_TOLERANGE 0.2

/** @brief Safety threshold for the Riemann solver wavespeed denominator.
 * Prevents division by zero or numerical "explosions" when left and
 * right wavespeeds are nearly identical. */
#define GEAR_FVPM_DIFF_WAVESPEED_ESTIMATE_DIFFERENCE_TOLERANCE 1e-8

/* --- GEAR FVPM Diffusion Flux Limiter Constants --- */

/** @brief Toggle for detailed logging of flux limiter events.
 * 0: Silent, 1+: Detailed per-interaction reports (high I/O overhead). */
#define GEAR_FVPM_DIFF_FLUX_LIMITER_VERBOSITY 0

/** @brief The "Noise Gate" threshold.
 * Fluxes smaller than 1e-15 relative to the source mass are treated as
 * numerical noise and zeroed out to prevent the "ratchet effect." */
#define GEAR_FVPM_DIFF_NOISE_GATE 1e-15

/** @brief Global safety factor for the rational flux limiter.
 * Ensures the attenuated flux remains well within the physical bounds. */
#define GEAR_FVPM_DIFF_LIMITER_SAFETY 0.5f

/** @brief Sink Capacity Constraint.
 * Restricts a single neighbor interaction from changing the sink particle's
 * metal mass by more than 25% to prevent neighbor-flooding oscillations. */
#define GEAR_FVPM_DIFF_LIMITER_SINK_STABILITY 0.25f

/** @brief Diffusion "Startup" Fraction.
 * Ensures pristine gas (Z=0) can be enriched. Allows a flux of up to 10%
 * of the source's mass to seed a metal-free neighbor. */
#define GEAR_FVPM_DIFF_LIMITER_STARTUP 0.1f

/* --- Error Handling & Debugging --- */

/** @brief Threshold for reporting persistent negative metal masses.
 * Number of consecutive timesteps a particle stays negative before
 * triggering a console warning. High values (200+) are for production. */
#define GEAR_FVPM_DIFF_NEGATIVITY_COUNTER_PRINT_LIMIT 200

/** @note DEBUG ONLY: Enforces one-sided updates for active pairs.
 * NOT MPI-compatible. Use only for serial/thread-testing of symmetry. */
/* #define GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE */

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
