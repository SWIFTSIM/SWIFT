/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2024 Roduit Darwin (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_MFM_DIFFUSION_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_MFM_DIFFUSION_H

#define GEAR_LABELS_SIZE 10  // redumndant with the one defined in

/**
 * @brief The diffusion mode
 */
enum chemistry_diffusion_mode {
  isotropic_constant,
  isotropic_smagorinsky, /* Smagorinsky turbulent diffusion \propto |S| */
  anisotropic_gradient   /* Rennehan (2021) model \propto S */
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

  /*! Diffusion mode. 0: isotropic with constant coefficient, 1: Smagorinsky
      isotrpoic diffusion, 2: anistropic diffusion with the shear tensor. */
  enum chemistry_diffusion_mode diffusion_mode;

  /***************************************************************************/
  /* HLL Riemann solver parameters
   * See Hopkins 2017 (https://arxiv.org/abs/1602.07703) */

  /*! The psi in eq (7) */
  float hll_riemann_solver_psi;

  /*! Use "usual" HLL or Hopkins 2017 modified HLL */
  int use_hokpins2017_hll_riemann_solver;

  /*! The epsilon in eq (14). This is a tolerance parameter. So, it must be 0
      <= epsilon <= 1. */
  float hll_riemann_solver_epsilon;

  /***************************************************************************/
  /* Supertimestepping */

  /*! Do we want to use supertimestepping? */
  int use_supertimestepping;

  /* Number of substeps */
  int N_substeps;

  /* Nu parameter */
  float nu;

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

#ifdef HYDRO_DOES_MASS_FLUX
  /* Note: This is only used by the MFV hydro scheme. */
  /*! Mass fluxes of the metals in a given element */
  double metal_mass_fluxes[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

  double diffusion_flux[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* Gradients. */
  struct {
    /* Gradient of the metals. It is used to compute the diffusion flux.
     */
    double Z[GEAR_CHEMISTRY_ELEMENT_COUNT][3];

    /* Fluid velocity gradients. */
    float v[3][3];

  } gradients;

  struct {
    /* Extreme values of the fluid metal_density among the neighbours. */
    double metal_density[GEAR_CHEMISTRY_ELEMENT_COUNT][2];

    /* Extreme values of the fluid velocity among the neighbours. */
    float v[3][2];

    /* Extreme values of the filtered velocity among the neighbours. */
    float v_tilde[3][2];

    /* Maximal distance to all neighbouring faces. */
    float maxr;

  } limiter;

  /* Geometrical quantities used for MFM/V hydro. */
  struct {

    /* Volume of the particle. */
    float volume;

    /* Geometrical shear matrix used to calculate second order accurate
       gradients */
    float matrix_E[3][3];

    /* Correction factor for wcount. */
    float wcorr;

    /* Condition number of matrix_E (eq C1) */
    float condition_number;

  } geometry;

  /* Particle chemistry time-step. */
  float flux_dt;

  /* Isotropic diffusion coefficient. */
  float kappa;

  /* Density of the previous timestep */
  float rho_prev;

  /* Here are the filtered quantities, i.e. "smoothed" over the resolution scale
   */
  struct {
    float rho;

    /* rho * v*/
    float rho_v[3];

    /* v_tilde = filtered(rho_v) / filtered(rho) */
    float grad_v_tilde[3][3];
  } filtered;

  struct {
    int current_substep;
    float explicit_timestep;
  } timesteps;
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

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_MFM_DIFFUSION_H */
