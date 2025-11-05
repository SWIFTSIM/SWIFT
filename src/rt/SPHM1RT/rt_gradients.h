/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_RT_GRAD_SPHM1RT_H
#define SWIFT_RT_GRAD_SPHM1RT_H

/**
 * @file SPHM1RT/rt_gradients.h
 * @brief Compute the gradient or divergence according to SPH
 */

/**
 * @brief Compute the gradient according to SPH
 * Note that the differential is: 1/rho * grad(rho * uin)
 * @param uini quantity for particle i to be differentiated
 * @param uinj quantity for particle j to be differentiated
 * @param mi mass of particle i
 * @param mj mass of particle j
 * @param forcefi a correction factor for spatial variations of hi
 * @param forcefj a correction factor for spatial variations of hj
 * @param rhoi density of particle i
 * @param rhoj density of particle j
 * @param wi_dr derivative of kernel i
 * @param wj_dr derivative of kernel j
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param r Comoving distance between the two particles.
 * @param diffmode mode=0 for ``difference''; mode=1 for ``symmetric''; mode=2
 * for ``difference'' but ignore h correction
 * @param gradi gradient of uin for particle i
 * @param gradj gradient of uin for particle j
 */
__attribute__((always_inline)) INLINE static void radiation_gradient_SPH(
    const float uini, const float uinj, const float mi, const float mj,
    const float forcefi, const float forcefj, const float rhoi,
    const float rhoj, const float wi_dr, const float wj_dr, const float dx[3],
    const float r, int diffmode, float gradi[3], float gradj[3]) {
  const float rho_i_inv = 1.0f / rhoi;
  const float rho_j_inv = 1.0f / rhoj;
  const float r_inv = 1.0f / r;
  if (diffmode == 0) {
    const float fradpre_common = (rhoi * uini - rhoj * uinj);
    const float fradprei =
        mj * forcefi * rho_i_inv * rho_i_inv * fradpre_common * wi_dr;
    const float fradprej =
        mi * forcefj * rho_j_inv * rho_j_inv * fradpre_common * wj_dr;
    gradi[0] = -fradprei * dx[0] * r_inv;
    gradi[1] = -fradprei * dx[1] * r_inv;
    gradi[2] = -fradprei * dx[2] * r_inv;
    gradj[0] = -fradprej * dx[0] * r_inv;
    gradj[1] = -fradprej * dx[1] * r_inv;
    gradj[2] = -fradprej * dx[2] * r_inv;
  } else if (diffmode == 1) {
    const float fradpre_common =
        uinj * rho_j_inv * wj_dr * forcefj + uini * rho_i_inv * wi_dr * forcefi;
    const float fradprei = mj * fradpre_common;
    const float fradprej = mi * fradpre_common;
    gradi[0] = fradprei * dx[0] * r_inv;
    gradi[1] = fradprei * dx[1] * r_inv;
    gradi[2] = fradprei * dx[2] * r_inv;
    gradj[0] = -fradprej * dx[0] * r_inv;
    gradj[1] = -fradprej * dx[1] * r_inv;
    gradj[2] = -fradprej * dx[2] * r_inv;
  } else if (diffmode == 2) {
    const float fradpre_common = (rhoi * uini - rhoj * uinj);
    const float fradprei =
        mj * rho_i_inv * rho_i_inv * fradpre_common * (wi_dr + wj_dr) * 0.5f;
    const float fradprej =
        mi * rho_j_inv * rho_j_inv * fradpre_common * (wi_dr + wj_dr) * 0.5f;
    gradi[0] = -fradprei * dx[0] * r_inv;
    gradi[1] = -fradprei * dx[1] * r_inv;
    gradi[2] = -fradprei * dx[2] * r_inv;
    gradj[0] = -fradprej * dx[0] * r_inv;
    gradj[1] = -fradprej * dx[1] * r_inv;
    gradj[2] = -fradprej * dx[2] * r_inv;
  } else {
    error("diffmode should be 0 or 1 or 2");
  }
}

/**
 * @brief Compute the anistropic gradient according to SPH
 * Note that the differential is: 1/rho * grad(Ftensor * rho * uin)
 * @param uini quantity for particle i to be differentiated
 * @param uinj quantity for particle j to be differentiated
 * @param mi mass of particle i
 * @param mj mass of particle j
 * @param forcefi a correction factor for spatial variations of hi
 * @param forcefj a correction factor for spatial variations of hj
 * @param rhoi density of particle i
 * @param rhoj density of particle j
 * @param wi_dr derivative of kernel i
 * @param wj_dr derivative of kernel j
 * @param Fanisoi Ftensor i
 * @param Fanisoj Ftensor j
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param r Comoving distance between the two particles.
 * @param diffmode mode=0 for ``difference''; mode=1 for ``symmetric''; mode=2
 * for ``difference'' but ignore h correction; mode=3
 * for ``symmetric'' but ignore h correction;
 * @param gradi gradient uin for particle i
 * @param gradj gradient uin for particle j
 */
__attribute__((always_inline)) INLINE static void radiation_gradient_aniso_SPH(
    const float uini, const float uinj, const float mi, const float mj,
    const float forcefi, const float forcefj, const float rhoi,
    const float rhoj, const float wi_dr, const float wj_dr, float Fanisoi[3][3],
    float Fanisoj[3][3], const float dx[3], const float r, int diffmode,
    float gradi[3], float gradj[3]) {
  const float rho_i_inv = 1.0f / rhoi;
  const float rho_j_inv = 1.0f / rhoj;
  const float r_inv = 1.0f / r;
  if (diffmode == 0) {
    float tempi[3], tempj[3];
    tempi[0] =
        Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] + Fanisoi[0][2] * dx[2];
    tempi[0] *= rhoi * uini * r_inv;
    tempi[1] =
        Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] + Fanisoi[1][2] * dx[2];
    tempi[1] *= rhoi * uini * r_inv;
    tempi[2] =
        Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] + Fanisoi[2][2] * dx[2];
    tempi[2] *= rhoi * uini * r_inv;
    tempj[0] =
        Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] + Fanisoj[0][2] * dx[2];
    tempj[0] *= rhoj * uinj * r_inv;
    tempj[1] =
        Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] + Fanisoj[1][2] * dx[2];
    tempj[1] *= rhoj * uinj * r_inv;
    tempj[2] =
        Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] + Fanisoj[2][2] * dx[2];
    tempj[2] *= rhoj * uinj * r_inv;
    gradi[0] =
        -(tempi[0] - tempj[0]) * mj * rho_i_inv * rho_i_inv * forcefi * wi_dr;
    gradi[1] =
        -(tempi[1] - tempj[1]) * mj * rho_i_inv * rho_i_inv * forcefi * wi_dr;
    gradi[2] =
        -(tempi[2] - tempj[2]) * mj * rho_i_inv * rho_i_inv * forcefi * wi_dr;
    gradj[0] =
        -(tempi[0] - tempj[0]) * mi * rho_j_inv * rho_j_inv * forcefj * wj_dr;
    gradj[1] =
        -(tempi[1] - tempj[1]) * mi * rho_j_inv * rho_j_inv * forcefj * wj_dr;
    gradj[2] =
        -(tempi[2] - tempj[2]) * mi * rho_j_inv * rho_j_inv * forcefj * wj_dr;
  } else if (diffmode == 1) {
    float tempi[3], tempj[3];
    tempi[0] =
        Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] + Fanisoi[0][2] * dx[2];
    tempi[0] *= uini * rho_i_inv * forcefi * r_inv * wi_dr;
    tempi[1] =
        Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] + Fanisoi[1][2] * dx[2];
    tempi[1] *= uini * rho_i_inv * forcefi * r_inv * wi_dr;
    tempi[2] =
        Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] + Fanisoi[2][2] * dx[2];
    tempi[2] *= uini * rho_i_inv * forcefi * r_inv * wi_dr;
    tempj[0] =
        Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] + Fanisoj[0][2] * dx[2];
    tempj[0] *= uinj * rho_j_inv * forcefj * r_inv * wj_dr;
    tempj[1] =
        Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] + Fanisoj[1][2] * dx[2];
    tempj[1] *= uinj * rho_j_inv * forcefj * r_inv * wj_dr;
    tempj[2] =
        Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] + Fanisoj[2][2] * dx[2];
    tempj[2] *= uinj * rho_j_inv * forcefj * r_inv * wj_dr;
    gradi[0] = (tempi[0] + tempj[0]) * mj;
    gradi[1] = (tempi[1] + tempj[1]) * mj;
    gradi[2] = (tempi[2] + tempj[2]) * mj;
    gradj[0] = -(tempi[0] + tempj[0]) * mi;
    gradj[1] = -(tempi[1] + tempj[1]) * mi;
    gradj[2] = -(tempi[2] + tempj[2]) * mi;
  } else if (diffmode == 2) {
    float tempi[3], tempj[3];
    tempi[0] =
        Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] + Fanisoi[0][2] * dx[2];
    tempi[0] *= rhoi * uini * r_inv;
    tempi[1] =
        Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] + Fanisoi[1][2] * dx[2];
    tempi[1] *= rhoi * uini * r_inv;
    tempi[2] =
        Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] + Fanisoi[2][2] * dx[2];
    tempi[2] *= rhoi * uini * r_inv;
    tempj[0] =
        Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] + Fanisoj[0][2] * dx[2];
    tempj[0] *= rhoj * uinj * r_inv;
    tempj[1] =
        Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] + Fanisoj[1][2] * dx[2];
    tempj[1] *= rhoj * uinj * r_inv;
    tempj[2] =
        Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] + Fanisoj[2][2] * dx[2];
    tempj[2] *= rhoj * uinj * r_inv;
    gradi[0] = -(tempi[0] - tempj[0]) * mj * rho_i_inv * rho_i_inv *
               (wi_dr + wj_dr) * 0.5f;
    gradi[1] = -(tempi[1] - tempj[1]) * mj * rho_i_inv * rho_i_inv *
               (wi_dr + wj_dr) * 0.5f;
    gradi[2] = -(tempi[2] - tempj[2]) * mj * rho_i_inv * rho_i_inv *
               (wi_dr + wj_dr) * 0.5f;
    gradj[0] = -(tempi[0] - tempj[0]) * mi * rho_j_inv * rho_j_inv *
               (wi_dr + wj_dr) * 0.5f;
    gradj[1] = -(tempi[1] - tempj[1]) * mi * rho_j_inv * rho_j_inv *
               (wi_dr + wj_dr) * 0.5f;
    gradj[2] = -(tempi[2] - tempj[2]) * mi * rho_j_inv * rho_j_inv *
               (wi_dr + wj_dr) * 0.5f;
  } else if (diffmode == 3) {
    float tempi[3], tempj[3];
    tempi[0] =
        Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] + Fanisoi[0][2] * dx[2];
    tempi[0] *= uini * rho_i_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    tempi[1] =
        Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] + Fanisoi[1][2] * dx[2];
    tempi[1] *= uini * rho_i_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    tempi[2] =
        Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] + Fanisoi[2][2] * dx[2];
    tempi[2] *= uini * rho_i_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    tempj[0] =
        Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] + Fanisoj[0][2] * dx[2];
    tempj[0] *= uinj * rho_j_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    tempj[1] =
        Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] + Fanisoj[1][2] * dx[2];
    tempj[1] *= uinj * rho_j_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    tempj[2] =
        Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] + Fanisoj[2][2] * dx[2];
    tempj[2] *= uinj * rho_j_inv * r_inv * (wi_dr + wj_dr) * 0.5f;
    gradi[0] = (tempi[0] + tempj[0]) * mj;
    gradi[1] = (tempi[1] + tempj[1]) * mj;
    gradi[2] = (tempi[2] + tempj[2]) * mj;
    gradj[0] = -(tempi[0] + tempj[0]) * mi;
    gradj[1] = -(tempi[1] + tempj[1]) * mi;
    gradj[2] = -(tempi[2] + tempj[2]) * mi;
  } else {
    error("diffmode should be 0 or 1 or 2 or 3");
  }
}

/**
 * @brief Compute the divergence according to SPH
 * Note that the differential is: 1/rho * div(rho * fin)
 * @param fini quantity for particle i to be differentiated
 * @param finj quantity for particle j to be differentiated
 * @param mi mass of particle i
 * @param mj mass of particle j
 * @param forcefi a correction factor for spatial variations of hi
 * @param forcefj a correction factor for spatial variations of hj
 * @param rhoi density of particle i
 * @param rhoj density of particle j
 * @param wi_dr derivative of kernel i
 * @param wj_dr derivative of kernel j
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param r Comoving distance between the two particles.
 * @param diffmode mode=0 for ``difference''; mode=1 for ``symmetric''; mode=2
 * for ``difference'' but ignore h correction
 * @param divfi divergence of f for particle i
 * @param divfj divergence of f for particle j
 */
__attribute__((always_inline)) INLINE static void radiation_divergence_SPH(
    const float *fini, const float *finj, const float mi, const float mj,
    const float forcefi, const float forcefj, const float rhoi,
    const float rhoj, const float wi_dr, const float wj_dr, const float dx[3],
    const float r, int diffmode, float *divfi, float *divfj) {
  const float rho_i_inv = 1.0f / rhoi;
  const float rho_j_inv = 1.0f / rhoj;
  const float r_inv = 1.0f / r;
  if (diffmode == 0) {
    float drhof[3];
    const float faci = mj * forcefi * rho_i_inv * rho_i_inv * wi_dr;
    const float facj = mi * forcefj * rho_j_inv * rho_j_inv * wj_dr;
    /* Compute drhof dot r */
    drhof[0] = rhoi * fini[0] - rhoj * finj[0];
    drhof[1] = rhoi * fini[1] - rhoj * finj[1];
    drhof[2] = rhoi * fini[2] - rhoj * finj[2];
    const float drhofdr =
        drhof[0] * dx[0] + drhof[1] * dx[1] + drhof[2] * dx[2];
    *divfi = -faci * drhofdr * r_inv;
    *divfj = -facj * drhofdr * r_inv;
  } else if (diffmode == 1) {
    float uradpre_common = 0.0f;
    uradpre_common += fini[0] * rho_i_inv * forcefi * dx[0] * r_inv * wi_dr;
    uradpre_common += fini[1] * rho_i_inv * forcefi * dx[1] * r_inv * wi_dr;
    uradpre_common += fini[2] * rho_i_inv * forcefi * dx[2] * r_inv * wi_dr;
    uradpre_common += finj[0] * rho_j_inv * forcefj * dx[0] * r_inv * wj_dr;
    uradpre_common += finj[1] * rho_j_inv * forcefj * dx[1] * r_inv * wj_dr;
    uradpre_common += finj[2] * rho_j_inv * forcefj * dx[2] * r_inv * wj_dr;
    *divfi = mj * uradpre_common;
    *divfj = -mi * uradpre_common;
  } else if (diffmode == 2) {
    float drhof[3];
    const float faci = mj * rho_i_inv * rho_i_inv * (wi_dr + wj_dr) * 0.5f;
    const float facj = mi * rho_j_inv * rho_j_inv * (wi_dr + wj_dr) * 0.5f;
    /* Compute drhof dot r */
    drhof[0] = rhoi * fini[0] - rhoj * finj[0];
    drhof[1] = rhoi * fini[1] - rhoj * finj[1];
    drhof[2] = rhoi * fini[2] - rhoj * finj[2];
    const float drhofdr =
        drhof[0] * dx[0] + drhof[1] * dx[1] + drhof[2] * dx[2];
    *divfi = -faci * drhofdr * r_inv;
    *divfj = -facj * drhofdr * r_inv;
  } else {
    error("diffmode should be 0 or 1 or 2");
  }
}

/**
 * @brief Compute the anistropic divergence according to SPH
 * Note that the differential is: 1/rho * div(Ftensor * rho * fin)
 * @param fini quantity for particle i to be differentiated
 * @param finj quantity for particle j to be differentiated
 * @param mi mass of particle i
 * @param mj mass of particle j
 * @param forcefi a correction factor for spatial variations of hi
 * @param forcefj a correction factor for spatial variations of hj
 * @param rhoi density of particle i
 * @param rhoj density of particle j
 * @param wi_dr derivative of kernel i
 * @param wj_dr derivative of kernel j
 * @param Fanisoi Ftensor i
 * @param Fanisoj Ftensor j
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param r Comoving distance between the two particles.
 * @param diffmode mode=0 for ``difference''; mode=1 for ``symmetric''; mode=2
 * for ``difference'' but ignore h correction
 * @param divfi divergence of f for particle i
 * @param divfj divergence of f for particle j
 */
__attribute__((always_inline)) INLINE static void
radiation_divergence_aniso_SPH(const float *fini, const float *finj,
                               const float mi, const float mj,
                               const float forcefi, const float forcefj,
                               const float rhoi, const float rhoj,
                               const float wi_dr, const float wj_dr,
                               float Fanisoi[3][3], float Fanisoj[3][3],
                               const float dx[3], const float r, int diffmode,
                               float *divfi, float *divfj) {
  const float rho_i_inv = 1.0f / rhoi;
  const float rho_j_inv = 1.0f / rhoj;
  const float r_inv = 1.0f / r;
  if (diffmode == 0) {
    float tempi, tempj;
    tempi = (Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] +
             Fanisoi[0][2] * dx[2]) *
            fini[0];
    tempi += (Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] +
              Fanisoi[1][2] * dx[2]) *
             fini[1];
    tempi += (Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] +
              Fanisoi[2][2] * dx[2]) *
             fini[2];
    tempi *= rhoi * r_inv;
    tempj = (Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] +
             Fanisoj[0][2] * dx[2]) *
            finj[0];
    tempj += (Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] +
              Fanisoj[1][2] * dx[2]) *
             finj[1];
    tempj += (Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] +
              Fanisoj[2][2] * dx[2]) *
             finj[2];
    tempj *= rhoj * r_inv;
    *divfi = -(tempi - tempj) * mj * forcefi * rho_i_inv * rho_i_inv * wi_dr;
    *divfj = -(tempi - tempj) * mi * forcefj * rho_j_inv * rho_j_inv * wj_dr;
  } else if (diffmode == 1) {
    float tempi, tempj;
    tempi = (Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] +
             Fanisoi[0][2] * dx[2]) *
            fini[0];
    tempi += (Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] +
              Fanisoi[1][2] * dx[2]) *
             fini[1];
    tempi += (Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] +
              Fanisoi[2][2] * dx[2]) *
             fini[2];
    tempi *= rho_i_inv * forcefi * r_inv * wi_dr;
    tempj = (Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] +
             Fanisoj[0][2] * dx[2]) *
            finj[0];
    tempj += (Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] +
              Fanisoj[1][2] * dx[2]) *
             finj[1];
    tempj += (Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] +
              Fanisoj[2][2] * dx[2]) *
             finj[2];
    tempj *= rho_j_inv * forcefj * r_inv * wj_dr;
    *divfi = (tempi + tempj) * mj;
    *divfj = -(tempi + tempj) * mi;
  } else if (diffmode == 2) {
    float tempi, tempj;
    tempi = (Fanisoi[0][0] * dx[0] + Fanisoi[0][1] * dx[1] +
             Fanisoi[0][2] * dx[2]) *
            fini[0];
    tempi += (Fanisoi[1][0] * dx[0] + Fanisoi[1][1] * dx[1] +
              Fanisoi[1][2] * dx[2]) *
             fini[1];
    tempi += (Fanisoi[2][0] * dx[0] + Fanisoi[2][1] * dx[1] +
              Fanisoi[2][2] * dx[2]) *
             fini[2];
    tempi *= rhoi * r_inv;
    tempj = (Fanisoj[0][0] * dx[0] + Fanisoj[0][1] * dx[1] +
             Fanisoj[0][2] * dx[2]) *
            finj[0];
    tempj += (Fanisoj[1][0] * dx[0] + Fanisoj[1][1] * dx[1] +
              Fanisoj[1][2] * dx[2]) *
             finj[1];
    tempj += (Fanisoj[2][0] * dx[0] + Fanisoj[2][1] * dx[1] +
              Fanisoj[2][2] * dx[2]) *
             finj[2];
    tempj *= rhoj * r_inv;
    *divfi =
        -(tempi - tempj) * mj * rho_i_inv * rho_i_inv * (wi_dr + wj_dr) * 0.5f;
    *divfj =
        -(tempi - tempj) * mi * rho_j_inv * rho_j_inv * (wi_dr + wj_dr) * 0.5f;
  } else {
    error("diffmode should be 0 or 1 or 2");
  }
}

/**
 * @brief Compute the divergence according to SPH
 * Note that the differential is: 1/rho * div(rho * fin)
 * @param fini quantity for particle i to be differentiated
 * @param finj quantity for particle j to be differentiated
 * @param mi mass of particle i
 * @param mj mass of particle j
 * @param forcefi a correction factor for spatial variations of hi
 * @param forcefj a correction factor for spatial variations of hj
 * @param rhoi density of particle i
 * @param rhoj density of particle j
 * @param wi_dr derivative of kernel i
 * @param wj_dr derivative of kernel j
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param r Comoving distance between the two particles.
 * @param diffmode mode=0 for ``difference''
 * @param shearfi divergence of f for particle i
 * @param shearfj divergence of f for particle j
 */
__attribute__((always_inline)) INLINE static void radiation_gradient_vec_SPH(
    const float *fini, const float *finj, const float mi, const float mj,
    const float forcefi, const float forcefj, const float rhoi,
    const float rhoj, const float wi_dr, const float wj_dr, const float dx[3],
    const float r, int diffmode, float shearfi[3][3], float shearfj[3][3]) {
  const float rho_i_inv = 1.0f / rhoi;
  const float rho_j_inv = 1.0f / rhoj;
  const float r_inv = 1.0f / r;
  if (diffmode == 0) {
    float drhof[3];
    const float faci = mj * forcefi * rho_i_inv * rho_i_inv * wi_dr;
    const float facj = mi * forcefj * rho_j_inv * rho_j_inv * wj_dr;
    /* Compute drhof dot r */
    drhof[0] = rhoi * fini[0] - rhoj * finj[0];
    drhof[1] = rhoi * fini[1] - rhoj * finj[1];
    drhof[2] = rhoi * fini[2] - rhoj * finj[2];
    for (int k = 0; k < 3; k++) {
      shearfi[k][0] = -faci * r_inv * drhof[0] * dx[k];
      shearfi[k][1] = -faci * r_inv * drhof[1] * dx[k];
      shearfi[k][2] = -faci * r_inv * drhof[2] * dx[k];
      shearfj[k][0] = -facj * r_inv * drhof[0] * dx[k];
      shearfj[k][1] = -facj * r_inv * drhof[1] * dx[k];
      shearfj[k][2] = -facj * r_inv * drhof[2] * dx[k];
    }
  } else {
    error("diffmode should be 0");
  }
}

#endif /* SWIFT_RT_GRAD_SPHM1RT_H */
