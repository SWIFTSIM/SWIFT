/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, float *xij_i, float *xij_j,
    float r) {

  float xij_i_norm;
  float phi_i, phi_j;
  float delta1, delta2;
  float phiminus, phiplus;
  float phimin, phimax;
  float phibar;
  /* free parameters, values from Hopkins */
  float psi1 = 0.5, psi2 = 0.25;
  float phi_mid0, phi_mid;
  int k;

  for (k = 0; k < 10; k++) {
    if (k < 5) {
      phi_i = Wi[k];
      phi_j = Wj[k];
      phi_mid0 = Wi[k] + dWi[k];
      xij_i_norm = sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] +
                         xij_i[2] * xij_i[2]);
    } else {
      phi_i = Wj[k - 5];
      phi_j = Wi[k - 5];
      phi_mid0 = Wj[k - 5] + dWj[k - 5];
      xij_i_norm = sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] +
                         xij_j[2] * xij_j[2]);
    }

    delta1 = psi1 * fabs(phi_i - phi_j);
    delta2 = psi2 * fabs(phi_i - phi_j);

    phimin = fmin(phi_i, phi_j);
    phimax = fmax(phi_i, phi_j);

    phibar = phi_i + xij_i_norm / r * (phi_j - phi_i);

    /* if sign(phimax+delta1) == sign(phimax) */
    if ((phimax + delta1) * phimax > 0.0f) {
      phiplus = phimax + delta1;
    } else {
      phiplus = phimax / (1.0f + delta1 / fabs(phimax));
    }

    /* if sign(phimin-delta1) == sign(phimin) */
    if ((phimin - delta1) * phimin > 0.0f) {
      phiminus = phimin - delta1;
    } else {
      phiminus = phimin / (1.0f + delta1 / fabs(phimin));
    }

    if (phi_i == phi_j) {
      phi_mid = phi_i;
    } else {
      if (phi_i < phi_j) {
        phi_mid = fmax(phiminus, fmin(phibar + delta2, phi_mid0));
      } else {
        phi_mid = fmin(phiplus, fmax(phibar - delta2, phi_mid0));
      }
    }

    if (k < 5) {
      dWi[k] = phi_mid - phi_i;
    } else {
      dWj[k - 5] = phi_mid - phi_i;
    }
  }
}
