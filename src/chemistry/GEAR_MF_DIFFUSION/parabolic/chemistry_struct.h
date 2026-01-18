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
#ifndef SWIFT_CHEMISTRY_STRUCT_GEAR_MF_PARABOLIC_DIFFUSION_H
#define SWIFT_CHEMISTRY_STRUCT_GEAR_MF_PARABOLIC_DIFFUSION_H

#include "../chemistry_properties.h"

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

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  /*! Total metal mass diffused during the simulation for this particle */
  double diffused_metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];
#endif

  /* Fluxes. */
  struct {
    /*! Metal mass flux computed with the Riemann solver */
    double metal_mass[GEAR_CHEMISTRY_ELEMENT_COUNT];

    /*! Particle chemistry time-step. */
    float dt;

  } flux;

  struct {
    /* Counts the successive number of times this particle was negative */
    unsigned int negativity_counter[GEAR_CHEMISTRY_ELEMENT_COUNT];
  } check;

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
    /*! Metal mass fraction gradient. It is used to compute the diffusion flux.
     */
    double Z[GEAR_CHEMISTRY_ELEMENT_COUNT][3];

    /*! Metal density gradient. */
    double rhoZ[GEAR_CHEMISTRY_ELEMENT_COUNT][3];

    /*! Fluid velocity gradients. */
    float v[3][3];

  } gradients;

  /* Cell-wise limiter to avoid creating new min or max */
  struct {
    /*! Extreme values of the fluid metal mass fraction among the neighbours.
     */
    double Z[GEAR_CHEMISTRY_ELEMENT_COUNT][2];

    /*! Extreme values of the fluid metal density among the neighbours. */
    double rhoZ[GEAR_CHEMISTRY_ELEMENT_COUNT][2];

    /*! Extreme values of the fluid velocity among the neighbours. */
    float v[3][2];

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

    /*! Gradient of tilde(v). Note: We do not need to slope limit this
        gradient because it is not used for conserved quantity
        reconstruction. This gradient is used to compute the traceless shear
        matrix. */
    float grad_v_tilde[3][3];
  } filtered;

  /*! Variables used for timestep calculation. */
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  struct {
    /* Mean separation distance from Lanson & Vila (2008) */
    float delxbar;

  } timestepvars;
#endif
};

#endif /* SWIFT_CHEMISTRY_STRUCT_GEAR_MF_PARABOLIC_DIFFUSION_H */
