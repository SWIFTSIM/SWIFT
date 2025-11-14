/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@ealumni.pfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GETTERS_H

#include "../chemistry_utils.h"
#include "chemistry_struct.h"
#include "const.h"
#include "cosmology.h"
#include "hydro.h"
#include "part.h"

/**
 * @brief Get the physical hyperbolic diffusion soundspeed.
 *
 * Note: The units are always U_L/U_T.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_physical_hyperbolic_soundspeed(
    const struct part *restrict p,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {
  if (chem_data->relaxation_time_mode == constant_mode) {
    double K[3][3];
    chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);
    const double norm_matrix_K = chemistry_get_matrix_norm(K);

    /* Here we simply use the formula c_hyp = sqrt(||K||/tau) */
    return sqrt(norm_matrix_K / p->chemistry_data.tau);
  } else {
    /* Note that 1/|S| ~ time --> we define this as our turbulent relaxation
       time. Also note that we do not regularize the shear tensor here.
       (Shall we?) */
    double S[3][3];
    chemistry_get_physical_shear_tensor(p, cosmo, S);

    /* TODO: Add the alpha parameter to the code */
    /* The formula is c_hyp = sqrt(||K||/(rho tau)). We simplify it by hand to
       reduce rounding errors: c_hyp = sqrt(C/alpha) * gamma_k * h * ||S|| */
    const double delta_x = kernel_gamma * p->h;
    const double C_diff = chem_data->diffusion_coefficient;
    const double alpha = chem_data->tau;
    const double c_hyp =
        sqrt(C_diff / alpha) * delta_x * chemistry_get_matrix_norm(S);
    return c_hyp;
  }
}

/**
 * @brief Get the physical hyperbolic diffusion relaxation time.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_compute_physical_tau(const struct part *restrict p,
                               const struct chemistry_global_data *chem_data,
                               const struct cosmology *cosmo) {
  if (chem_data->relaxation_time_mode == constant_mode) {
    /* Tau is constant and chosen in the parameter file. Hence return this
     * value. */
    return chem_data->tau;
  } else {
    /* Note that 1/|S| ~ time --> we define this as our turbulent relaxation
       time. Also note that we do not regularize the shear tensor here. */
    double S[3][3];
    chemistry_get_physical_shear_tensor(p, cosmo, S);
    const double S_norm_inv = 1.0 / chemistry_get_matrix_norm(S);

    return chem_data->tau * S_norm_inv;
  }
}

/**
 * @brief Get the gradients of diffusion flux for a given metal specie.
 *
 * Note: Gradients are comoving;
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dFx Gradient (of size 3) of the flux's x component.
 * @param dFy Gradient (of size 3) of the flux's y component.
 * @param dFz Gradient (of size 3) of the flux's z component.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_hyperbolic_flux_gradients(const struct part *restrict p,
                                        int metal, double dFx[3], double dFy[3],
                                        double dFz[3]) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  dFx[0] = chd->gradients.flux[metal][0][0];
  dFx[1] = chd->gradients.flux[metal][0][1];
  dFx[2] = chd->gradients.flux[metal][0][2];

  dFy[0] = chd->gradients.flux[metal][1][0];
  dFy[1] = chd->gradients.flux[metal][1][1];
  dFy[2] = chd->gradients.flux[metal][1][2];

  dFz[0] = chd->gradients.flux[metal][2][0];
  dFz[1] = chd->gradients.flux[metal][2][1];
  dFz[2] = chd->gradients.flux[metal][2][2];
}

/**
 * @brief Get the gradients of diffusion flux for a given metal specie.
 *
 * Note: Gradients are comoving;
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dFx Gradient (of size 3) of the flux's x component.
 * @param dFy Gradient (of size 3) of the flux's y component.
 * @param dFz Gradient (of size 3) of the flux's z component.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_physical_hyperbolic_flux_gradients(
    const struct part *restrict p, int metal, double dFx[3], double dFy[3],
    double dFz[3], const struct cosmology *cosmo) {
  chemistry_get_hyperbolic_flux_gradients(p, metal, dFx, dFy, dFz);

  /* grad_p = a^{-1} grad_c */
  for (int i = 0; i < 3; i++) {
    dFx[i] *= cosmo->a_inv;
    dFy[i] *= cosmo->a_inv;
    dFz[i] *= cosmo->a_inv;
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_GETTERS_H */
