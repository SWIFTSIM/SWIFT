/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H
#define SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H


/**
 * @brief Solve the Riemann problem for the RT equations and return the
 * flux at the interface.
 *
 * @param UL left state (radiation energy density, flux)
 * @param UR right state (radiation energy density, flux)
 * @param FLnorm the norm of the radiation flux of the left state
 * @param FRnorm the norm of the radiation flux of the right state
 * @param hyperFluxL the flux of the hyperbolic conservation law of the left
 * state
 * @param hyperFluxR the flux of the hyperbolic conservation law of the right
 * state
 * @param n_unit the unit vector perpendicular to the "intercell" surface.
 * @param flux_half (return) the resulting flux at the interface
 */
__attribute__((always_inline)) INLINE static void chemistry_riemann_solve_for_flux(
    const struct part* restrict p, const float UL, const float UR, const float* n_unit,
    const float Anorm, const float F_diff_i[3], const float F_diff_j[3], float metal_flux) {
										   
    /* const float UL[4], const float UR[4], const float FLnorm, */
    /* const float FRnorm, float hyperFluxL[4][3], float hyperFluxR[4][3], */
    /* const float n_unit[3], float flux_half[4]) { */

  /* /\* Compute lambda_+/- *\/ */
  /* const float lminus = min3(lambdaLmin, lambdaRmin, 0.f); */
  /* const float lplus = max3(lambdaLmax, lambdaRmax, 0.f); */

  /* /\* Sanity check: This should give the same results as GLF solver *\/ */
  /* if (lminus == 0.f && lplus == 0.f) { */
  /*   flux_half[0] = 0.f; */
  /*   flux_half[1] = 0.f; */
  /*   flux_half[2] = 0.f; */
  /*   flux_half[3] = 0.f; */
  /*   return; */
  /* } */

  /* Compute F_diff_ij^* */

  /* Compute F_diff_ij^* cdot A_ij */

  // Have a look in riemann_hllc.h 
}


#endif /* SWIFT_GEAR_CHEMISTRY_RIEMANN_HLL_H */
