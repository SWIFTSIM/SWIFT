/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"

/**
 * @brief Compute the particle parabolic timestep proportional to h^2.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static float chemistry_diffusion_timestep(
    const struct part *restrict p,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

  const float CFL_condition = chem_data->C_CFL_chemistry;

#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  /* Lanson & Vila (2008) particle size */
  /* Note that we computed 1/delxbar so we need to take the inverse here. And
     we need to convert to physical units. */
  const float psize = cosmo->a * 0.5f / p->chemistry_data.timestepvars.delxbar;
#else
  /* Gizmo's particle size definition */
  const float psize = cosmo->a * powf(p->geometry.volume / hydro_dimension_unit_sphere, hydro_dimension_inv);
#endif

  const float delta_x = psize;

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  /* CFL condition */
  const float dt_cfl =
      CFL_condition * delta_x / p->chemistry_data.timestepvars.vmax;
  return dt_cfl;
#else
  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute the diffusion matrix K */
  double K[3][3];
  chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0) {
    return FLT_MAX;
  }

  if (chem_data->diffusion_mode == isotropic_constant) {
    /* Isotropic constant diffusion has the simple expression: */
    return CFL_condition * delta_x * delta_x / norm_matrix_K;
  }
  /* From here, we are either isotropic smagorinksy or anisotropic */

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...), and q = (Z_1, Z_2,
     ...). Hence, the term norm(U)/norm(q) in eq (15) is abs(rho). */
  const float norm_U_over_norm_q = hydro_get_physical_density(p, cosmo);

  /* Some helpful variables */
  float norm_q = 0.0;
  float norm_nabla_q = 0.0;
  float expression = 0.0;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_q += chemistry_get_metal_mass_fraction(p, i) *
              chemistry_get_metal_mass_fraction(p, i);

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q = Grad Z */
      norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt. Add the missing a^{-1} to have physical gradients */
  norm_q = sqrtf(norm_q);  // Physical
  norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv;

  /* If ||q|| is 0 (metal density = 0), then there is no metal to diffuse.
     Hence, no timestep limit is needed. */
  if (norm_q == 0) {
    return FLT_MAX;
  }

  if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* We should use the same expression as the anisotropic diffusion below.
       However, it is overly conservative for isotropic diffusion and slows
       down the simulations by order of magnitudes...
       Instead, we use a similar expression to the isotropic constant
       diffusion. */
    expression = delta_x;
  } else {
    /* Compute the expression in the square bracket in eq (15). Notice that I
       rewrote it to avoid division by 0 when norm_nabla_q = 0.
       This expression ensures the correctness of parabolic diffusion. But at
       the cost of being not usable is cosmo simulations... The SN enrichment
       creates huge gradients that slow down everything. */
    expression = norm_q * delta_x / (norm_nabla_q * delta_x + norm_q);
  }

  const float dt = CFL_condition * expression * expression / norm_matrix_K *
                   norm_U_over_norm_q;
  return dt;
#endif
}
#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H  */
