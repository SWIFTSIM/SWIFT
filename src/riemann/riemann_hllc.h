/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 *               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_RIEMANN_HLLC_H
#define SWIFT_RIEMANN_HLLC_H

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "error.h"
#include "minmax.h"
#include "riemann_checks.h"
#include "riemann_vacuum.h"

#ifndef EOS_IDEAL_GAS
#error \
    "The HLLC Riemann solver currently only supports and ideal gas equation of state. Either select this equation of state, or try using another Riemann solver!"
#endif

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    const float *WL, const float *WR, const float *n, const float *vij,
    float *totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(WL, WR, n, vij);
#endif

  /* Handle pure vacuum */
  if (!WL[0] && !WR[0]) {
    totflux[0] = 0.f;
    totflux[1] = 0.f;
    totflux[2] = 0.f;
    totflux[3] = 0.f;
    totflux[4] = 0.f;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
  const float uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
  const float aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  const float aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    riemann_solve_vacuum_flux(WL, WR, uL, uR, aL, aR, n, vij, totflux);
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = 0.5f * (WL[0] + WR[0]);
  const float abar = 0.5f * (aL + aR);
  const float pPVRS = 0.5f * (WL[4] + WR[4]) - 0.5f * (uR - uL) * rhobar * abar;
  const float pstar = max(0.f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.f;
  if (pstar > WL[4] && WL[4] > 0.f) {
    qL = sqrtf(1.f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                         (pstar / WL[4] - 1.f));
  }
  float qR = 1.f;
  if (pstar > WR[4] && WR[4] > 0.f) {
    qR = sqrtf(1.f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                         (pstar / WR[4] - 1.f));
  }
  const float SL = uL - aL * qL;
  const float SR = uR + aR * qR;
  const float Sstar =
      (WR[4] - WL[4] + WL[0] * uL * (SL - uL) - WR[0] * uR * (SR - uR)) /
      (WL[0] * (SL - uL) - WR[0] * (SR - uR));

  /* STEP 3: HLLC flux in a frame moving with the interface velocity */
  if (Sstar >= 0.f) {
    /* flux FL */
    totflux[0] = WL[0] * uL;
    /* these are the actual correct fluxes in the boosted lab frame
       (not rotated to interface frame) */
    totflux[1] = WL[0] * WL[1] * uL + WL[4] * n[0];
    totflux[2] = WL[0] * WL[2] * uL + WL[4] * n[1];
    totflux[3] = WL[0] * WL[3] * uL + WL[4] * n[2];
    const float v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
    const float eL = WL[4] * hydro_one_over_gamma_minus_one / WL[0] + 0.5f * v2;
    totflux[4] = WL[0] * eL * uL + WL[4] * uL;
    if (SL < 0.f) {

      float UstarL[5];

      /* add flux FstarL */
      UstarL[0] = 1.f;
      /* we need UstarL in the lab frame:
       * subtract the velocity in the interface frame from the lab frame
       * velocity and then add Sstar in interface frame */
      UstarL[1] = WL[1] + (Sstar - uL) * n[0];
      UstarL[2] = WL[2] + (Sstar - uL) * n[1];
      UstarL[3] = WL[3] + (Sstar - uL) * n[2];
      UstarL[4] = eL + (Sstar - uL) * (Sstar + WL[4] / (WL[0] * (SL - uL)));
      UstarL[0] *= WL[0] * (SL - uL) / (SL - Sstar);
      UstarL[1] *= WL[0] * (SL - uL) / (SL - Sstar);
      UstarL[2] *= WL[0] * (SL - uL) / (SL - Sstar);
      UstarL[3] *= WL[0] * (SL - uL) / (SL - Sstar);
      UstarL[4] *= WL[0] * (SL - uL) / (SL - Sstar);
      totflux[0] += SL * (UstarL[0] - WL[0]);
      totflux[1] += SL * (UstarL[1] - WL[0] * WL[1]);
      totflux[2] += SL * (UstarL[2] - WL[0] * WL[2]);
      totflux[3] += SL * (UstarL[3] - WL[0] * WL[3]);
      totflux[4] += SL * (UstarL[4] - WL[0] * eL);
    }
  } else {
    /* flux FR */
    totflux[0] = WR[0] * uR;
    totflux[1] = WR[0] * WR[1] * uR + WR[4] * n[0];
    totflux[2] = WR[0] * WR[2] * uR + WR[4] * n[1];
    totflux[3] = WR[0] * WR[3] * uR + WR[4] * n[2];
    const float v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
    const float eR = WR[4] * hydro_one_over_gamma_minus_one / WR[0] + 0.5f * v2;
    totflux[4] = WR[0] * eR * uR + WR[4] * uR;
    if (SR > 0.f) {

      float UstarR[5];

      /* add flux FstarR */
      UstarR[0] = 1.f;

      /* we need UstarR in the lab frame:
       * subtract the velocity in the interface frame from the lab frame
       * velocity and then add Sstar in interface frame */
      UstarR[1] = WR[1] + (Sstar - uR) * n[0];
      UstarR[2] = WR[2] + (Sstar - uR) * n[1];
      UstarR[3] = WR[3] + (Sstar - uR) * n[2];
      UstarR[4] = eR + (Sstar - uR) * (Sstar + WR[4] / (WR[0] * (SR - uR)));
      UstarR[0] *= WR[0] * (SR - uR) / (SR - Sstar);
      UstarR[1] *= WR[0] * (SR - uR) / (SR - Sstar);
      UstarR[2] *= WR[0] * (SR - uR) / (SR - Sstar);
      UstarR[3] *= WR[0] * (SR - uR) / (SR - Sstar);
      UstarR[4] *= WR[0] * (SR - uR) / (SR - Sstar);
      totflux[0] += SR * (UstarR[0] - WR[0]);
      totflux[1] += SR * (UstarR[1] - WR[0] * WR[1]);
      totflux[2] += SR * (UstarR[2] - WR[0] * WR[2]);
      totflux[3] += SR * (UstarR[3] - WR[0] * WR[3]);
      totflux[4] += SR * (UstarR[4] - WR[0] * eR);
    }
  }

  /* deboost to lab frame we add the flux contribution due to the
     movement of the interface the density flux is unchanged
     we add the extra velocity flux due to the absolute motion of the fluid
     similarly, we need to add the energy fluxes due to the absolute motion */
  const float v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];

  /* order is important: we first use the momentum fluxes to update the energy
     flux and then de-boost the momentum fluxes! */
  totflux[4] += vij[0] * totflux[1] + vij[1] * totflux[2] +
                vij[2] * totflux[3] + 0.5f * v2 * totflux[0];
  totflux[1] += vij[0] * totflux[0];
  totflux[2] += vij[1] * totflux[0];
  totflux[3] += vij[2] * totflux[0];

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(WL, WR, n, vij, totflux);
#endif
}

__attribute__((always_inline)) INLINE static void
riemann_solve_for_middle_state_flux(const float *WL, const float *WR,
                                    const float *n, const float *vij,
                                    float *totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(WL, WR, n, vij);
#endif

  /* Handle pure vacuum */
  if (!WL[0] && !WR[0]) {
    totflux[0] = 0.f;
    totflux[1] = 0.f;
    totflux[2] = 0.f;
    totflux[3] = 0.f;
    totflux[4] = 0.f;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
  const float uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
  const float aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  const float aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    totflux[0] = 0.f;
    totflux[1] = 0.f;
    totflux[2] = 0.f;
    totflux[3] = 0.f;
    totflux[4] = 0.f;
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = 0.5f * (WL[0] + WR[0]);
  const float abar = 0.5f * (aL + aR);
  const float pPVRS = 0.5f * (WL[4] + WR[4]) - 0.5f * (uR - uL) * rhobar * abar;
  const float pstar = max(0.f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.f;
  if (pstar > WL[4] && WL[4] > 0.f) {
    qL = sqrtf(1.f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                         (pstar / WL[4] - 1.f));
  }
  float qR = 1.f;
  if (pstar > WR[4] && WR[4] > 0.f) {
    qR = sqrtf(1.f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                         (pstar / WR[4] - 1.f));
  }
  const float SL = uL - aL * qL;
  const float SR = uR + aR * qR;
  const float Sstar =
      (WR[4] - WL[4] + WL[0] * uL * (SL - uL) - WR[0] * uR * (SR - uR)) /
      (WL[0] * (SL - uL) - WR[0] * (SR - uR));

  totflux[0] = 0.0f;
  totflux[1] = pstar * n[0];
  totflux[2] = pstar * n[1];
  totflux[3] = pstar * n[2];
  const float vface = vij[0] * n[0] + vij[1] * n[1] + vij[2] * n[2];
  totflux[4] = pstar * (Sstar + vface);

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_output(WL, WR, n, vij, totflux);
#endif
}

#endif /* SWIFT_RIEMANN_HLLC_H */
