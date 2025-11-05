/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_RIEMANN_HLL_H
#define SWIFT_GEAR_RT_RIEMANN_HLL_H

#define RT_RIEMANN_HLL_NPOINTS 100
#define RT_RIEMANN_HLL_DF (1.f / (float)(RT_RIEMANN_HLL_NPOINTS - 1))
#define RT_RIEMANN_HLL_ONE_OVER_DF ((float)RT_RIEMANN_HLL_NPOINTS - 1)
#define RT_RIEMANN_HLL_DTHETA (M_PI / ((float)(RT_RIEMANN_HLL_NPOINTS - 1)))
#define RT_RIEMANN_HLL_ONE_OVER_DTHETA \
  (((float)(RT_RIEMANN_HLL_NPOINTS - 1)) / M_PI)

#include "rt_getters.h"
#include "rt_parameters.h"
#include "rt_riemann_HLL_eigenvalues.h"
#include "rt_unphysical.h"

/**
 * @file src/rt/GEAR/rt_riemann_HLL.h
 * @brief  The Hartmann-Lax-van Leer Riemann solver for the moments of the
 * radiative transfer equation following following Gonzalez et al 2007
 * (ui.adsabs.harvard.edu/abs/2007A%26A...464..429G).
 * */

/**
 * @brief Interpolate the minimal and maximal eigenvalues from
 * the lookup table given the reduced flux f and angle w.r.t.
 * the surface theta.
 *
 * @param f reduced flux |F|/(c E)
 * @param theta angle between flux and surface
 * @param lambda_min minimal eigenvalue
 * @param lambda_max maximal eigenvalue
 */

__attribute__((always_inline)) INLINE static void
rt_riemann_interpolate_eigenvals(float f, float theta, float *lambda_min,
                                 float *lambda_max) {

  /* find lower table indices for f and theta */
  int f_ind = floor(RT_RIEMANN_HLL_ONE_OVER_DF * f);
  if (f_ind >= RT_RIEMANN_HLL_NPOINTS - 1) f_ind = RT_RIEMANN_HLL_NPOINTS - 2;
  int theta_ind = floor(RT_RIEMANN_HLL_ONE_OVER_DTHETA * theta);
  if (theta_ind >= RT_RIEMANN_HLL_NPOINTS - 1)
    theta_ind = RT_RIEMANN_HLL_NPOINTS - 2;

  /* Grab the data */
  const float Q11min = rt_riemann_HLL_eigenvals[f_ind][theta_ind][0];
  const float Q12min = rt_riemann_HLL_eigenvals[f_ind][theta_ind + 1][0];
  const float Q21min = rt_riemann_HLL_eigenvals[f_ind + 1][theta_ind][0];
  const float Q22min = rt_riemann_HLL_eigenvals[f_ind + 1][theta_ind + 1][0];
  const float Q11max = rt_riemann_HLL_eigenvals[f_ind][theta_ind][1];
  const float Q12max = rt_riemann_HLL_eigenvals[f_ind][theta_ind + 1][1];
  const float Q21max = rt_riemann_HLL_eigenvals[f_ind + 1][theta_ind][1];
  const float Q22max = rt_riemann_HLL_eigenvals[f_ind + 1][theta_ind + 1][1];

  /* (f - f1)/(f2 - f1)  =  (f - f_ind * Delta f)/Delta f */
  const float df1 = f * RT_RIEMANN_HLL_ONE_OVER_DF - (float)f_ind;
  /* (f2 - f)/(f2 - f1)  =  ((f_ind + 1) * Delta f - f)/Delta f */
  const float df2 = 1.f - df1;

  /* linear interpolation in f direction */
  const float fmid1_min = df2 * Q11min + df1 * Q21min;
  const float fmid2_min = df2 * Q12min + df1 * Q22min;
  const float fmid1_max = df2 * Q11max + df1 * Q21max;
  const float fmid2_max = df2 * Q12max + df1 * Q22max;

  /* Now second interpolation in theta direction */
  /* (theta - theta1)/(theta2 - theta1) */
  const float dtheta1 =
      theta * RT_RIEMANN_HLL_ONE_OVER_DTHETA - (float)theta_ind;
  /* (theta2 - theta)/(theta2 - theta1) */
  const float dtheta2 = 1.f - dtheta1;

  /* Make sure -1 < eigenvalue < 1 */
  float lmin = dtheta2 * fmid1_min + dtheta1 * fmid2_min;
  lmin = max(lmin, -1.f);
  lmin = min(lmin, 1.f);
  *lambda_min = lmin;
  float lmax = dtheta2 * fmid1_max + dtheta1 * fmid2_max;
  lmax = max(lmax, -1.f);
  lmax = min(lmax, 1.f);
  *lambda_max = lmax;
}

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
__attribute__((always_inline)) INLINE static void rt_riemann_solve_for_flux(
    const float UL[4], const float UR[4], const float FLnorm,
    const float FRnorm, float hyperFluxL[4][3], float hyperFluxR[4][3],
    const float n_unit[3], float flux_half[4]) {

  /* Compute reduced fluxes and angles between surface and flux.
   * These are based on physical fluxes, not hyperbolic fluxes. */
  const float c_red = rt_params.reduced_speed_of_light;
  const float c_red_inv = rt_params.reduced_speed_of_light_inverse;

  float fL = 0.f;
  float thetaL = 0.f;

  if (UL[0] > 0.f) {
    fL = FLnorm / UL[0] * c_red_inv;
    fL = min(fL, 1.f);
  }

  if (FLnorm > 0.f) {
    /* cos theta = F * n / (|F| |n|) */
    const float FLdotn =
        UL[1] * n_unit[0] + UL[2] * n_unit[1] + UL[3] * n_unit[2];
    float costhetaL = min(FLdotn / FLnorm, 1.f);
    costhetaL = max(costhetaL, -1.f);
    thetaL = acosf(costhetaL);
  }

  float fR = 0.f;
  float thetaR = 0.f;

  if (UR[0] > 0.f) {
    fR = FRnorm / UR[0] * c_red_inv;
    fR = min(fR, 1.f);
  }

  if (FRnorm > 0.f) {
    const float FRdotn =
        UR[1] * n_unit[0] + UR[2] * n_unit[1] + UR[3] * n_unit[2];
    float costhetaR = min(FRdotn / FRnorm, 1.f);
    costhetaR = max(costhetaR, -1.f);
    thetaR = acosf(costhetaR);
  }

  /* interpolate eigenvalues in lookup table */
  float lambdaLmin = 0;
  float lambdaLmax = 0;
  rt_riemann_interpolate_eigenvals(fL, thetaL, &lambdaLmin, &lambdaLmax);
  float lambdaRmin = 0;
  float lambdaRmax = 0;
  rt_riemann_interpolate_eigenvals(fR, thetaR, &lambdaRmin, &lambdaRmax);

  const float lminus = min3(lambdaLmin, lambdaRmin, 0.f);
  const float lplus = max3(lambdaLmax, lambdaRmax, 0.f);

  /* Sanity check: This should give the same results as GLF solver */
  /* const float lminus = -1.f; */
  /* const float lplus = +1.f; */

  if (lminus == 0.f && lplus == 0.f) {
    flux_half[0] = 0.f;
    flux_half[1] = 0.f;
    flux_half[2] = 0.f;
    flux_half[3] = 0.f;
    return;
  }

  /* Project the (hyperbolic) flux along the surface,
   * reduce the problem to 1D with 4 quantities*/
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

  const float one_over_dl = 1.f / (lplus - lminus);
  /* Remember that the eigenvalues are in units of c */
  const float lprod = lplus * lminus * c_red;

  flux_half[0] =
      (lplus * fluxL[0] - lminus * fluxR[0] + lprod * (UR[0] - UL[0])) *
      one_over_dl;
  flux_half[1] =
      (lplus * fluxL[1] - lminus * fluxR[1] + lprod * (UR[1] - UL[1])) *
      one_over_dl;
  flux_half[2] =
      (lplus * fluxL[2] - lminus * fluxR[2] + lprod * (UR[2] - UL[2])) *
      one_over_dl;
  flux_half[3] =
      (lplus * fluxL[3] - lminus * fluxR[3] + lprod * (UR[3] - UL[3])) *
      one_over_dl;
}

#endif /* SWIFT_GEAR_RT_RIEMANN_HLL_H */
