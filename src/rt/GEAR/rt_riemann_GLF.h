/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_RIEMANN_GLF_H
#define SWIFT_GEAR_RT_RIEMANN_GLF_H

#include "rt_getters.h"
#include "rt_parameters.h"
#include "rt_unphysical.h"

/**
 * @file src/rt/GEAR/rt_riemann_GLF.h
 * @brief  The Global Lax-Friedrich (GLF) Riemann solver for the moments of the
 * radiative transfer equation following Rosdahl et al 2013
 * */

/**
 * @brief Solve the Riemann problem for the RT equations and return the
 * flux at the interface.
 *
 * @param UL left state
 * @param UR right state
 * @param flux_half the resulting flux at the interface
 * @param n_unit the unit vector perpendicular to the "intercell" surface.
 */
__attribute__((always_inline)) INLINE static void rt_riemann_solve_for_flux(
    const float UL[4], const float UR[4], float flux_half[4],
    const float n_unit[3]) {

  float hyperFluxL[4][3];
  rt_get_hyperbolic_flux(UL, hyperFluxL);
  float hyperFluxR[4][3];
  rt_get_hyperbolic_flux(UR, hyperFluxR);

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_check_unphysical_hyperbolic_flux(hyperFluxL);
  rt_check_unphysical_hyperbolic_flux(hyperFluxR);
#endif

  float fluxL[4];
  fluxL[0] = hyperFluxL[0][0] * n_unit[0] + hyperFluxL[0][1] * n_unit[1] +
             hyperFluxL[0][2] * n_unit[2];
  fluxL[1] = hyperFluxL[1][0] * n_unit[0] + hyperFluxL[1][1] * n_unit[1] +
             hyperFluxL[1][2] * n_unit[2];
  fluxL[2] = hyperFluxL[2][0] * n_unit[0] + hyperFluxL[2][1] * n_unit[1] +
             hyperFluxL[2][2] * n_unit[2];
  fluxL[3] = hyperFluxL[3][0] * n_unit[0] + hyperFluxL[3][1] * n_unit[1] +
             hyperFluxL[3][2] * n_unit[2];

  float fluxR[4];
  fluxR[0] = hyperFluxR[0][0] * n_unit[0] + hyperFluxR[0][1] * n_unit[1] +
             hyperFluxR[0][2] * n_unit[2];
  fluxR[1] = hyperFluxR[1][0] * n_unit[0] + hyperFluxR[1][1] * n_unit[1] +
             hyperFluxR[1][2] * n_unit[2];
  fluxR[2] = hyperFluxR[2][0] * n_unit[0] + hyperFluxR[2][1] * n_unit[1] +
             hyperFluxR[2][2] * n_unit[2];
  fluxR[3] = hyperFluxR[3][0] * n_unit[0] + hyperFluxR[3][1] * n_unit[1] +
             hyperFluxR[3][2] * n_unit[2];

  const float c_red = rt_params.reduced_speed_of_light;
  float cdU[4] = {c_red * (UR[0] - UL[0]), c_red * (UR[1] - UL[1]),
                  c_red * (UR[2] - UL[2]), c_red * (UR[3] - UL[3])};

  flux_half[0] = 0.5f * (fluxL[0] + fluxR[0] - cdU[0]);
  flux_half[1] = 0.5f * (fluxL[1] + fluxR[1] - cdU[1]);
  flux_half[2] = 0.5f * (fluxL[2] + fluxR[2] - cdU[2]);
  flux_half[3] = 0.5f * (fluxL[3] + fluxR[3] - cdU[3]);
}

#endif /* SWIFT_GEAR_RT_RIEMANN_GLF_H */
