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
 * @brief Global chemical abundance information.
 */
struct chemistry_global_data {

  /* Initial mass fraction */
  double initial_metallicities[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* Solar mass abundances read from the chemistry table */
  float solar_abundances[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /*! Name of the different elements */
  char elements_name[GEAR_CHEMISTRY_ELEMENT_COUNT * GEAR_LABELS_SIZE];

  /*! Diffusion normalisation constant: \kappa \propto C */
  float diffusion_coefficient;

  /*! Is the diffusion isotropic? */
  int use_isotropic_diffusion;

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

  /*! Smoothed fraction of the particle mass in a given element */
  double smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* This is U.  */
  /* struct { */
  /*   double metal_density; */
  /* } conserved[GEAR_CHEMISTRY_ELEMENT_COUNT]; */

  double diffusion_flux[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* Gradients. */
  struct {
    /* This is \nabla \otimes \vec{q}. It is used to compute the diffusion flux
     */
    double nabla_otimes_q[3];
  } gradients[GEAR_CHEMISTRY_ELEMENT_COUNT];

  struct {
    double metal_density[2];
  } limiter[GEAR_CHEMISTRY_ELEMENT_COUNT];

  /* we only need one maxr. Hence, put it outside the limiter struct */
  float limiter_maxr;

  /* Geometrical quantities used for MFM hydro. */
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
