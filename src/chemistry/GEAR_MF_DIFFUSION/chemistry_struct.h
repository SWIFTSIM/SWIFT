/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H

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
  constant_mode,   /* Constant */
  soundspeed_mode, /* Based on the gas soundspeed: tau = ||K||/(c_s^2 rho) */
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

/**
 * @brief Properties of the chemistry function for #part.
 */
struct chemistry_part_data {

  /*! Total mass of element in a particle. This is the primitive variable *
   * volume. */
  double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Metal mass flux computed with the Riemann solver */
  double diffused_metal_mass_fluxes[GEAR_CHEMISTRY_ELEMENT_COUNT];

#ifdef HYDRO_DOES_MASS_FLUX
  /* Note: This is only used by the MFV hydro scheme. */
  /*! Mass fluxes of the metals in a given element */
  double metal_mass_fluxes[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /*! Total metal mass diffused during the simulation for this particle */
  double diffused_metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

  /*! Condition number of matrix_E (eq C1) */
  float geometry_condition_number;

  /*! Particle chemistry time-step. */
  float flux_dt;

  /*! Isotropic diffusion coefficient. The matrix K is proportional to kappa.
   Note about units:
   - For the isotropic constant case, the units are : U_L^2/U_T
   - Smagorinsky/Gradient, units are : U_M/(U_L*U_T) */
  double kappa;

  /*! Density of the previous timestep. This is used to compute quantities in
     the density loop while hydro loops are updating rho. */
  float rho_prev;

  /* Gradients. */
  struct {
    /*! Gradient of the metals. It is used to compute the diffusion flux. */
    double Z[GEAR_CHEMISTRY_ELEMENT_COUNT][3];

    /*! Density gradient */
    float rho[3];

    /*! Fluid velocity gradients. */
    float v[3][3];

  } gradients;

  /* Cell-wise limiter to avoid creating new min or max */
  struct {
    /*! Extreme values of the fluid metal mass fraction  among the neighbours.
     */
    double Z[GEAR_CHEMISTRY_ELEMENT_COUNT][2];

    /*! Extreme values of the density among the neigbours. */
    float rho[2];

    /*! Extreme values of the fluid velocity among the neighbours. */
    float v[3][2];

    /*! Extreme values of the filtered velocity among the neighbours. */
    float v_tilde[3][2];

    /*! Maximal distance to all neighbouring faces. */
    float maxr;

  } limiter;

  /* Here are the filtered quantities, i.e. "smoothed" over the resolution
     scale h_bar = \gamma_k h */
  struct {

    /*! Filtered density (rho_bar) */
    float rho;

    /*! Filtered density (rho_bar) of the previous timstep. This is used to
       compute quantities in the density loop while we are updating rho_bar. */
    float rho_prev;

    /*! Filtered (rho * v). tilde(v) = filtered(rho*v) / filtered(rho) */
    float rho_v[3];

    /*! Gradient of tilde(v) */
    float grad_v_tilde[3][3];
  } filtered;

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  /* Relaxation time */
  double tau;

  /* Hyperbolic flux scheme variables */
  struct {
    /*! Diffusion flux at the last active timestep. Units: U_M/(U_L^2 U_T) */
    double F_diff[3];

    /*! Predicted diffusion flux */
    double F_diff_pred[3];

    /*! Time derivative of the diffusion flux */
    double dF_dt[3];
  } hyperbolic_flux[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Variables used for timestep calculation. */
  struct {
    /* Maximum signal velocity among all the neighbours of the particle. The
     * signal velocity encodes information about the relative fluid velocities
     * AND particle velocities of the neighbour and this particle, as well as
     * the sound speed of both particles. */
    float vmax;

  } timestepvars;
#endif
};

/**
 * @brief Properties of the chemistry function for #spart.
 */
struct chemistry_spart_data {

  /*! Fraction of the particle mass in a given element */
  double metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

/**
 * @brief Chemical abundances traced by the #bpart in the GEAR model.
 */
struct chemistry_bpart_data {};

/**
 * @brief Chemical abundances traced by the #sink in the GEAR model.
 */
struct chemistry_sink_data {

  /*! Total mass of element in a particle. */
  double metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];
};

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_MF_DIFFUSION_H */
