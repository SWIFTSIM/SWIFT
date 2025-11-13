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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H

#include "../chemistry_gradients.h"
#include "../chemistry_riemann_checks.h"
#include "../chemistry_riemann_utils.h"
#include "../chemistry_struct.h"
#include "hydro.h"
#include "sign.h"

// TODO: Same
// Test again with the hyperbolic soundspeed. Maybe take the max of c_s and
// c_hyp.
/**
 * @brief HLL riemann solver.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param hyper_flux_L The flux of the hyperbolic conservation law of the left
 * (in physical units).
 * @param hyper_flux_R The flux of the hyperbolic conservation law of the right
 * (in physical units).
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param fluxes (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_solver_HLL(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const double hyper_flux_L[4][3],
    const double hyper_flux_R[4][3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double fluxes[4]) {

  /***************************************************************************/
  /* Estimate the eigenvalue of the Jacobian matrix dF/dU */
  /* Everything is in physical units here */

  /* PVRS wavespeed approximation */
  const double uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
  const double uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];
  const double rhoLinv = (WL[0] > 0.0f) ? 1.0f / WL[0] : 0.0f;
  const double rhoRinv = (WR[0] > 0.0f) ? 1.0f / WR[0] : 0.0f;
  const double aL = sqrtf(hydro_gamma * WL[4] * rhoLinv);
  const double aR = sqrtf(hydro_gamma * WR[4] * rhoRinv);

  /* Pressure estimate */
  const double rhobar = WL[0] + WR[0];
  const double abar = aL + aR;
  const double pPVRS =
      0.5f * ((WL[4] + WR[4]) - 0.25f * (uR - uL) * rhobar * abar);
  const double pstar = max(0.0f, pPVRS);

  /* Wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  double qL = 1.0f;
  if (pstar > WL[4] && WL[4] > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
			  (pstar / WL[4] - 1.0f));
  }
  double qR = 1.0f;
  if (pstar > WR[4] && WR[4] > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
			  (pstar / WR[4] - 1.0f));
  }
  const double lambda_minus_hydro = -aL * qL;
  const double lambda_plus_hydro = aR * qR;

  /* const double c_diff_L = */
      /* chemistry_get_physical_hyperbolic_soundspeed(pi, chem_data, cosmo); */
  /* const double c_diff_R = */
      /* chemistry_get_physical_hyperbolic_soundspeed(pj, chem_data, cosmo); */
  /* const double lambda_plus_diffusion = max3(0.0, c_diff_L, c_diff_R); */
  /* const double lambda_minus_diffusion = min3(0.0, c_diff_L, c_diff_R); */

  const double lambda_minus = lambda_minus_hydro;
  const double lambda_plus = lambda_plus_hydro;
  /* const double lambda_minus = min(lambda_minus_hydro, lambda_minus_diffusion); */
  /* const double lambda_plus = max(lambda_plus_hydro, lambda_plus_diffusion); */

  if (lambda_plus == 0.0 && lambda_minus == 0.0) {
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
  }

  /***************************************************************************/
  /* Now project the fluxes */
  /* No conversion to physical needed, everything is physical here */

  /* Project the fluxes to reduce to a 1D Problem with 1 quantity */
  /* const double Flux_L = F_diff_L[0] * n_unit[0] + F_diff_L[1] * n_unit[1] + */
			/* F_diff_L[2] * n_unit[2]; */
  /* const double Flux_R = F_diff_R[0] * n_unit[0] + F_diff_R[1] * n_unit[1] + */
			/* F_diff_R[2] * n_unit[2]; */

  /***************************************************************************/
  /* Now solve the Riemann problem */

  /* const double dU = UR - UL; */

  /* const double delta_lambda = lambda_plus - lambda_minus; */
  /* if (fabs(delta_lambda) < 1e-8) { */
  /*   *metal_flux = 0.0; */
  /*   return; */
  /* } */

  /* const double one_over_dl = 1.f / (lambda_plus - lambda_minus); */
  /* const double F_2 = */
  /*     (lambda_plus * Flux_L - lambda_minus * Flux_R) * one_over_dl; */
  /* double F_U = lambda_plus * lambda_minus * dU * one_over_dl; */

  /* if (lambda_minus > 0.0) { */
  /*   *metal_flux = Flux_L; */
  /* } else if (lambda_minus <= 0.0 && lambda_plus >= 0.0) { */
  /*   *metal_flux = F_2 + F_U; */
  /* } else if (lambda_plus < 0.0) { */
  /*   *metal_flux = Flux_R; */
  /* } */
}

// TODO: Update arguments to take U[4]...
/**
 * @brief Solve the Riemann problem for the diffusion equations and return the
 * flux at the interface.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi Left particle
 * @param pj Right particle
 * @param UL left diffusion state (metal density, in physical units)
 * @param UR right diffusion state (metal density, in physical units)
 * @param WL Left state hydrodynamics primitve variables (density, velocity[3],
 * pressure) (in physical units)
 * @param WR Right state hydrodynamics primitve variables (density,
 * velocity[3], pressure) (in physical units)
 * @param hyper_flux_L The flux of the hyperbolic conservation law of the left
 * (in physical units).
 * @param hyper_flux_R The flux of the hyperbolic conservation law of the right
 * (in physical units).
 * @param Anorm Norm of the face between the left and right particles (in
 * physical units)
 * @param n_unit The unit vector perpendicular to the "intercell" surface.
 * @param m Index of metal specie to update.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param fluxes (return) The resulting flux at the interface.
 */
__attribute__((always_inline)) INLINE static void
chemistry_riemann_solve_for_flux(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const double hyper_flux_L[4][3],
    const double hyper_flux_R[4][3], const float Anorm, const float n_unit[3],
    const int m, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, double fluxes[4]) {

  chemistry_riemann_check_input(WL, WR, UL, UR, n_unit);

  /* Handle pure vacuum */
  if ((!UL[0] && !UR[0]) || (!WL[0] && !WR[0])) {
    fluxes[0] = 0.0;
    fluxes[1] = 0.0;
    fluxes[2] = 0.0;
    fluxes[3] = 0.0;
    return;
  }

  /* No conversion to physical needed, everything is physical here */
  chemistry_riemann_solver_HLL(dx, pi, pj, UL, UR, WL, WR, hyper_flux_L, hyper_flux_R,
			       Anorm, n_unit, m, chem_data, cosmo, fluxes);

  chemistry_riemann_check_output(WL, WR, UL, UR, n_unit, fluxes);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_RIEMANN_HLL_H */
