/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#include "riemann_vacuum.h"

#ifndef EOS_IDEAL_GAS
#error \
    "The HLLC Riemann solver currently only supports and ideal gas equation of state. Either select this equation of state, or try using another Riemann solver!"
#endif

__attribute__((always_inline)) INLINE static void riemann_solve_for_flux(
    float *WL, float *WR, float *n, float *vij, float *totflux) {

  float uL, uR, aL, aR;
  float rhobar, abar, pPVRS, pstar, qL, qR, SL, SR, Sstar;
  float v2, eL, eR;
  float UstarL[5], UstarR[5];

  /* Handle pure vacuum */
  if (!WL[0] && !WR[0]) {
    totflux[0] = 0.;
    totflux[1] = 0.;
    totflux[2] = 0.;
    totflux[3] = 0.;
    totflux[4] = 0.;
    return;
  }

  /* STEP 0: obtain velocity in interface frame */
  uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
  uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
  aL = sqrtf(hydro_gamma * WL[4] / WL[0]);
  aR = sqrtf(hydro_gamma * WR[4] / WR[0]);

  /* Handle vacuum: vacuum does not require iteration and is always exact */
  if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    riemann_solve_vacuum_flux(WL, WR, uL, uR, aL, aR, n, vij, totflux);
    return;
  }

  /* STEP 1: pressure estimate */
  rhobar = 0.5 * (WL[0] + WR[0]);
  abar = 0.5 * (aL + aR);
  pPVRS = 0.5 * (WL[4] + WR[4]) - 0.5 * (uR - uL) * rhobar * abar;
  pstar = max(0., pPVRS);

  /* STEP 2: wave speed estimates
     all these speeds are along the interface normal, since uL and uR are */
  qL = 1.;
  if (pstar > WL[4]) {
    qL = sqrtf(1. +
               0.5 * (hydro_gamma + 1.) / hydro_gamma * (pstar / WL[4] - 1.));
  }
  qR = 1.;
  if (pstar > WR[4]) {
    qR = sqrtf(1. +
               0.5 * (hydro_gamma + 1.) / hydro_gamma * (pstar / WR[4] - 1.));
  }
  SL = uL - aL * qL;
  SR = uR + aR * qR;
  Sstar = (WR[4] - WL[4] + WL[0] * uL * (SL - uL) - WR[0] * uR * (SR - uR)) /
          (WL[0] * (SL - uL) - WR[0] * (SR - uR));

  /* STEP 3: HLLC flux in a frame moving with the interface velocity */
  if (Sstar >= 0.) {
    /* flux FL */
    totflux[0] = WL[0] * uL;
    /* these are the actual correct fluxes in the boosted lab frame
       (not rotated to interface frame) */
    totflux[1] = WL[0] * WL[1] * uL + WL[4] * n[0];
    totflux[2] = WL[0] * WL[2] * uL + WL[4] * n[1];
    totflux[3] = WL[0] * WL[3] * uL + WL[4] * n[2];
    v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
    eL = WL[4] / hydro_gamma_minus_one / WL[0] + 0.5 * v2;
    totflux[4] = WL[0] * eL * uL + WL[4] * uL;
    if (SL < 0.) {
      /* add flux FstarL */
      UstarL[0] = 1.;
      /* we need UstarL in the lab frame:
         subtract the velocity in the interface frame from the lab frame
         velocity and then add Sstar in interface frame */
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
    v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
    eR = WR[4] / hydro_gamma_minus_one / WR[0] + 0.5 * v2;
    totflux[4] = WR[0] * eR * uR + WR[4] * uR;
    if (SR > 0.) {
      /* add flux FstarR */
      UstarR[0] = 1.;
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

  /* deboost to lab frame
     we add the flux contribution due to the movement of the interface
     the density flux is unchanged
     we add the extra velocity flux due to the absolute motion of the fluid
     similarly, we need to add the energy fluxes due to the absolute motion */
  v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];
  // order is important: we first use the momentum fluxes to update the energy
  // flux and then de-boost the momentum fluxes!
  totflux[4] += vij[0] * totflux[1] + vij[1] * totflux[2] +
                vij[2] * totflux[3] + 0.5 * v2 * totflux[0];
  totflux[1] += vij[0] * totflux[0];
  totflux[2] += vij[1] * totflux[0];
  totflux[3] += vij[2] * totflux[0];
}

#endif /* SWIFT_RIEMANN_HLLC_H */
