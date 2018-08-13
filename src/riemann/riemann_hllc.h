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

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    const float *WL, const float *WR, const float *n, const float *vij,
    float *totflux) {

#ifdef SWIFT_DEBUG_CHECKS
  riemann_check_input(WL, WR, n, vij);
#endif

  /* Handle pure vacuum */
  if (!WL[0] && !WR[0]) {
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
  const float uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
  const float rhoLinv = 1.0f / WL[0];
  const float rhoRinv = 1.0f / WR[0];
  const float aL = sqrtf(hydro_gamma * WL[4] * rhoLinv);
  const float aR = sqrtf(hydro_gamma * WR[4] * rhoRinv);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    riemann_solve_vacuum_flux(WL, WR, uL, uR, aL, aR, n, vij, totflux);
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = WL[0] + WR[0];
  const float abar = aL + aR;
  const float pPVRS =
      0.5f * ((WL[4] + WR[4]) - 0.25f * (uR - uL) * rhobar * abar);
  const float pstar = max(0.0f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.0f;
  if (pstar > WL[4] && WL[4] > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WL[4] - 1.0f));
  }
  float qR = 1.0f;
  if (pstar > WR[4] && WR[4] > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WR[4] - 1.0f));
  }
  const float SLmuL = -aL * qL;
  const float SRmuR = aR * qR;
  const float Sstar =
      (WR[4] - WL[4] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
      (WL[0] * SLmuL - WR[0] * SRmuR);

  /* STEP 3: HLLC flux in a frame moving with the interface velocity */
  if (Sstar >= 0.0f) {
    const float rhoLuL = WL[0] * uL;
    const float v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
    const float eL =
        WL[4] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5f * v2;
    const float SL = SLmuL + uL;

    /* flux FL */
    totflux[0] = rhoLuL;
    /* these are the actual correct fluxes in the boosted lab frame
       (not rotated to interface frame) */
    totflux[1] = rhoLuL * WL[1] + WL[4] * n[0];
    totflux[2] = rhoLuL * WL[2] + WL[4] * n[1];
    totflux[3] = rhoLuL * WL[3] + WL[4] * n[2];
    totflux[4] = rhoLuL * eL + WL[4] * uL;

    if (SL < 0.0f) {

      const float starfac = SLmuL / (SL - Sstar) - 1.0f;
      const float rhoLSL = WL[0] * SL;
      const float SstarmuL = Sstar - uL;
      const float rhoLSLstarfac = rhoLSL * starfac;
      const float rhoLSLSstarmuL = rhoLSL * SstarmuL;

      totflux[0] += rhoLSLstarfac;
      totflux[1] += rhoLSLstarfac * WL[1] + rhoLSLSstarmuL * n[0];
      totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[1];
      totflux[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n[2];
      totflux[4] += rhoLSLstarfac * eL +
                    rhoLSLSstarmuL * (Sstar + WL[4] / (WL[0] * SLmuL));
    }
  } else {
    const float rhoRuR = WR[0] * uR;
    const float v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
    const float eR =
        WR[4] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5f * v2;
    const float SR = SRmuR + uR;

    /* flux FR */
    totflux[0] = rhoRuR;
    totflux[1] = rhoRuR * WR[1] + WR[4] * n[0];
    totflux[2] = rhoRuR * WR[2] + WR[4] * n[1];
    totflux[3] = rhoRuR * WR[3] + WR[4] * n[2];
    totflux[4] = rhoRuR * eR + WR[4] * uR;

    if (SR > 0.0f) {

      const float starfac = SRmuR / (SR - Sstar) - 1.0f;
      const float rhoRSR = WR[0] * SR;
      const float SstarmuR = Sstar - uR;
      const float rhoRSRstarfac = rhoRSR * starfac;
      const float rhoRSRSstarmuR = rhoRSR * SstarmuR;

      totflux[0] += rhoRSRstarfac;
      totflux[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n[0];
      totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[1];
      totflux[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n[2];
      totflux[4] += rhoRSRstarfac * eR +
                    rhoRSRSstarmuR * (Sstar + WR[4] / (WR[0] * SRmuR));
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
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  const float uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
  const float uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
  const float aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  const float aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    totflux[0] = 0.0f;
    totflux[1] = 0.0f;
    totflux[2] = 0.0f;
    totflux[3] = 0.0f;
    totflux[4] = 0.0f;
    return;
  }

  /* STEP 1: pressure estimate */
  const float rhobar = WL[0] + WR[0];
  const float abar = aL + aR;
  const float pPVRS =
      0.5f * ((WL[4] + WR[4]) - 0.25f * (uR - uL) * rhobar * abar);
  const float pstar = max(0.f, pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  float qL = 1.0f;
  if (pstar > WL[4] && WL[4] > 0.0f) {
    qL = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WL[4] - 1.0f));
  }
  float qR = 1.0f;
  if (pstar > WR[4] && WR[4] > 0.0f) {
    qR = sqrtf(1.0f + 0.5f * hydro_gamma_plus_one * hydro_one_over_gamma *
                          (pstar / WR[4] - 1.0f));
  }
  const float SLmuL = -aL * qL;
  const float SRmuR = aR * qR;
  const float Sstar =
      (WR[4] - WL[4] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
      (WL[0] * SLmuL - WR[0] * SRmuR);

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
