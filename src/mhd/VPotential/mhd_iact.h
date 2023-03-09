/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_VECTOR_POTENTIAL_MHD_IACT_H
#define SWIFT_VECTOR_POTENTIAL_MHD_IACT_H

#include "periodic.h"

/**
 * @brief MHD-Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const float a, const float H) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  double dA[3];
  for (int i = 0; i < 3; ++i)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];

  const double dAdr = dA[0] * dx[0] + dA[1] * dx[1] + dA[2] * dx[2];
  pi->mhd_data.divA -= faci * dAdr;
  pj->mhd_data.divA -= facj * dAdr;
  /////
  // bi = dj ak - dk aj
  // bj = dk ai - di ak
  // bk = di aj - dj ai
  //
  for (int i = 0; i < 3; ++i) {
    pi->mhd_data.BPred[i] += faci * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);
    pj->mhd_data.BPred[i] += facj * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);
  }
  return;
}

/**
 * @brief MHD-Density interaction between two particles. (non-symmetric)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mhd_density(const float r2, const float dx[3],
                               const float hi, const float hj,
                               struct part *restrict pi,
                               const struct part *restrict pj,
                               const double mu_0, const float a,
                               const float H) {
  float wi, wi_dx;

  const float r = sqrtf(r2);

  /* Get the mass. */
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  double dA[3];
  for (int i = 0; i < 3; ++i)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];

  const double dAdr = dA[0] * dx[0] + dA[1] * dx[1] + dA[2] * dx[2];
  pi->mhd_data.divA -= faci * dAdr;
  for (int i = 0; i < 3; ++i)
    pi->mhd_data.BPred[i] += faci * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);

  return;
}

/**
 * @brief Calculate the MHD-gradient interaction between particle i and particle
 * j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const float a, const float H) {

  float wi, wi_dx, wj, wj_dx;
  const float r = sqrtf(r2);

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float wj_dr = hjd_inv * wj_dx;

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = pj->mass * wi_dr * r_inv / pi->rho;
  const float facj = pi->mass * wj_dr * r_inv / pj->rho;

  double dB[3];
  for (int i = 0; i < 3; ++i)
    dB[i] = pi->mhd_data.BPred[i] - pj->mhd_data.BPred[i];
  const double dBdr = dB[0] * dx[0] + dB[1] * dx[1] + dB[2] * dx[2];
  pi->mhd_data.divB -= faci * dBdr;
  pj->mhd_data.divB -= facj * dBdr;

  for (int i = 0; i < 3; i++) {
    pi->mhd_data.BSmooth[i] += pj->mass * wi * pj->mhd_data.BPred[i];
    pj->mhd_data.BSmooth[i] += pi->mass * wj * pi->mhd_data.BPred[i];
  }
  // pi->mhd_data.GauSmooth += pj->mass * wi * pj->mhd_data.Gau;
  // pj->mhd_data.GauSmooth += pi->mass * wj * pi->mhd_data.Gau;
  pi->mhd_data.Q0 += pj->mass * wi;
  pj->mhd_data.Q0 += pi->mass * wj;

  return;
}

/**
 * @brief Calculate the MHDgradient interaction between particle i and particle
 * j (non-symmetric)
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mhd_gradient(const float r2, const float dx[3],
                                const float hi, const float hj,
                                struct part *restrict pi,
                                const struct part *restrict pj,
                                const double mu_0, const float a,
                                const float H) {
  float wi, wi_dx;
  const float r = sqrtf(r2);

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = pj->mass * wi_dr * r_inv / pi->rho;

  double dB[3];
  for (int i = 0; i < 3; ++i)
    dB[i] = pi->mhd_data.BPred[i] - pj->mhd_data.BPred[i];
  const double dBdr = dB[0] * dx[0] + dB[1] * dx[1] + dB[2] * dx[2];
  pi->mhd_data.divB -= faci * dBdr;

  for (int i = 0; i < 3; i++)
    pi->mhd_data.BSmooth[i] += pj->mass * wi * pj->mhd_data.BPred[i];
  // pi->mhd_data.GauSmooth += pj->mass * wi * pj->mhd_data.Gau;
  pi->mhd_data.Q0 += pj->mass * wi;

  return;
}

/**
 * @brief MHD-Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const float a, const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
  const float rho_ij = rhoi + rhoj;

  const float a_fac =
      pow(a, 2.f * mhd_comoving_factor + 3.f * (hydro_gamma - 1.f));

  const float mag_faci = f_ij * wi_dr * r_inv / (rhoi * rhoi) / mu_0 * a_fac;
  const float mag_facj = f_ji * wj_dr * r_inv / (rhoj * rhoj) / mu_0 * a_fac;
  float Bi[3], Bj[3];
  float mm_i[3][3], mm_j[3][3];

  for (int i = 0; i < 3; i++) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }

  ///////////////////////////// FORCE MAXWELL TENSOR
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      mm_i[i][j] = Bi[i] * Bi[j];
      mm_j[i][j] = Bj[i] * Bj[j];
    }
  for (int j = 0; j < 3; j++) {
    mm_i[j][j] -= 0.5 * (Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2]);
    mm_j[j][j] -= 0.5 * (Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2]);
  }
  //////////////////////////// Apply to the Force and DIVB TERM SUBTRACTION
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      pi->a_hydro[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pj->a_hydro[i] -=
          mi * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->a_hydro[i] -= pi->mhd_data.Q0 * mj * Bi[i] *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
      pj->a_hydro[i] += pj->mhd_data.Q0 * mi * Bj[i] *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
    }
  /////////////////////////// VP evolution
  const float mag_VPIndi = wi_dr * r_inv / rhoi;
  const float mag_VPIndj = wj_dr * r_inv / rhoj;
  // Normal Gauge
  double dA[3];
  for (int i = 0; i < 3; i++)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];
  /*float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float SourceAi = dv[0] * pi->mhd_data.APred[0] +
                         dv[1] * pi->mhd_data.APred[1] +
                         dv[2] * pi->mhd_data.APred[2];
  const float SourceAj = dv[0] * pj->mhd_data.APred[0] +
                         dv[1] * pj->mhd_data.APred[1] +
                         dv[2] * pj->mhd_data.APred[2];
  */
  const float SourceAi =
      -(dA[0] * pi->v[0] + dA[1] * pi->v[1] + dA[2] * pi->v[2]);
  const float SourceAj =
      -(dA[0] * pj->v[0] + dA[1] * pj->v[1] + dA[2] * pj->v[2]);
  float SAi = SourceAi + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  float SAj = SourceAj + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);

  for (int i = 0; i < 3; i++) {
    pi->mhd_data.dAdt[i] += mj * mag_VPIndi * SAi * dx[i];
    pj->mhd_data.dAdt[i] += mi * mag_VPIndj * SAj * dx[i];
  }
  /// DISSSIPATION
  const float mag_Disi =
      (wi_dx + wj_dx) / 2.f * r_inv * rhoi / (rho_ij * rho_ij);
  const float mag_Disj =
      (wj_dx + wi_dx) / 2.f * r_inv * rhoj / (rho_ij * rho_ij);
  for (int i = 0; i < 3; i++) {
    pi->mhd_data.dAdt[i] += mj * 2.0 * pi->mhd_data.Deta * mag_Disi * dA[i];
    pj->mhd_data.dAdt[i] += mi * 2.0 * pj->mhd_data.Deta * mag_Disj * dA[i];
  }
  return;
}

/**
 * @brief MHD-Force interaction between two particles. non-symmetric version.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_mhd_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const double mu_0,
    const float a, const float H) {

  /* Cosmological factors entering the EoMs */
  // const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  // const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  // const float pressurei = pi->force.pressure;
  // const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
  const float rho_ij = rhoi + rhoj;

  const float a_fac =
      pow(a, 2.f * mhd_comoving_factor + 3.f * (hydro_gamma - 1.f));

  const float mag_faci = f_ij * wi_dr * r_inv / (rhoi * rhoi) / mu_0 * a_fac;
  const float mag_facj = f_ji * wj_dr * r_inv / (rhoj * rhoj) / mu_0 * a_fac;
  float Bi[3], Bj[3];
  float mm_i[3][3], mm_j[3][3];

  for (int i = 0; i < 3; i++) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }

  ///////////////////////////// FORCE MAXWELL TENSOR
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      mm_i[i][j] = Bi[i] * Bi[j];
      mm_j[i][j] = Bj[i] * Bj[j];
    }
  for (int j = 0; j < 3; j++) {
    mm_i[j][j] -= 0.5 * (Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2]);
    mm_j[j][j] -= 0.5 * (Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2]);
  }
  //////////////////////////// Apply to the Force and DIVB TERM SUBTRACTION
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      pi->a_hydro[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->a_hydro[i] -= pi->mhd_data.Q0 * mj * Bi[i] *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
    }
  /////////////////////////// VP INDUCTION
  const float mag_VPIndi = wi_dr * r_inv / rhoi;
  // Normal Gauge
  double dA[3];
  for (int i = 0; i < 3; i++)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];
  // float dv[3];
  // dv[0] = pi->v[0] - pj->v[0];
  // dv[1] = pi->v[1] - pj->v[1];
  // dv[2] = pi->v[2] - pj->v[2];
  // const float SourceAi = dv[0] * pi->mhd_data.APred[0] +
  //                        dv[1] * pi->mhd_data.APred[1] +
  //                        dv[2] * pi->mhd_data.APred[2];

  const float SourceAi =
      -(dA[0] * pi->v[0] + dA[1] * pi->v[1] + dA[2] * pi->v[2]);
  float SAi = SourceAi + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  for (int i = 0; i < 3; i++)
    pi->mhd_data.dAdt[i] += mj * mag_VPIndi * SAi * dx[i];
  /// DISSSIPATION
  const float mag_Disi =
      (wi_dx + wj_dx) / 2.f * r_inv * rhoi / (rho_ij * rho_ij);
  for (int i = 0; i < 3; i++)
    pi->mhd_data.dAdt[i] += mj * 2.0 * pi->mhd_data.Deta * mag_Disi * dA[i];
  return;
}
#endif /* SWIFT_VECTOR_POTENTIAL_MHD_H */
