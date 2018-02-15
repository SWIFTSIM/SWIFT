/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H

/**
 * @file PressureEnergy/hydro_iact.h
 * @brief Pressure-Energy implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (16), (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
 */

/**
 * @brief Density loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float dv[3], curlvr[3]; 

  /* Get r and r inverse. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Get the masses and energies */
  const float mj = pj->mass;
  const float int_energy_j = pj->u;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float hj_inv = 1.0f / hj;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  /* Note that wcount is \bar{n} */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx) * hj_inv;

  /* Compute the 'real' SPH density for other uses */
  const float wimj = wi * mj;
  pi->rho += wimj;

  /* Compute the contribution to the local pressure */
  /* The (gamma - 1) factor is added in hydro_end_density. */
  pi->pressure_bar += int_energy_j * wimj;
  pi->pressure_bar_dh -=
    mj * int_energy_j * hj_inv * (hydro_dimension * wi + ui * wi_dx);

  /* The following is lifted directly from the PressureEntropy code */
  const float fac = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  pi->density.div_v -= fac * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += fac * curlvr[0];
  pi->density.rot_v[1] += fac * curlvr[1];
  pi->density.rot_v[2] += fac * curlvr[2];
}

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* For now, just do the simple case */

  runner_iact_nonsym_density(r2, dx, hi, hj, pi, pj);
  runner_iact_nonsym_density(r2, dx, hj, hi, pj, pi);
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {
  
  float wi, wi_dx, wj, wj_dx;
  const float fac_mu = 1.f; /* Will change with cosmological integration */
  float dv[3];

  /* Masses */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get r and r inverse */
  const float r = sqrtf(r2);
  const float ri = 1.0f / r;

  /* Get velocity difference */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  const float dvdr = dv[0] * dx[0] +
                     dv[1] * dx[1] +
                     dv[2] * dx[2];

  /* Compute kernel function */
  const float hi_inv = 1.0f / hi;
  const float hj_inv = 1.0f / hj;
  const float hid_inv = pow_dimension_plus_one(hi_inv);
  const float hjd_inv = pow_dimension_plus_one(hj_inv);

  const float ui = r * hi_inv;
  const float uj = r * hj_inv;
  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  const float wi_dr = hid_inv * wi_dx;
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute gradient terms */

  const float f_ij = hydro_h_term(pi, pj, hi);
  const float f_ji = hydro_h_term(pj, pi, hj);

  /* Artificial viscosity */

  const float omega_ij = fminf(dvdr, 0.f);
  const float mu_ij = fac_mu * ri * omega_ij;

  const float v_sig = pi->force.soundspeed + pj->force.soundspeed - 3.f * mu_ij;

  const float rho_ij = 0.5f * (pi->rho + pj->rho);
  const float visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij *
                     (pi->force.balsara + pj->force.balsara) / rho_ij;
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr);

  /* Calculate the equation of motion H13:15 */

  const float sph_term = mj * pj->u * pi->u *
    hydro_gamma_minus_one * hydro_gamma_minus_one *
    ((f_ij/pi->pressure_bar) * wi_dr + (f_ji/pj->pressure_bar) * wj_dr);

  const float acc = (visc_term + sph_term) * ri;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= acc * dx[0];
  pi->a_hydro[1] -= acc * dx[1];
  pi->a_hydro[2] -= acc * dx[2];

  /* Get the time derivative for h */
  pi->force.h_dt -= mj * dvdr * ri / pj->rho * wi_dr;

  /* Update the signal velocity */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);

  /* Calcualte change in energy from viscosity */
  const float u_dt_visc = mj * visc_term * dvdr * ri;

  /* Calculate the change in internal energy from hydro H13:16 */
  const float u_dt_hydro =
    hydro_gamma_minus_one * hydro_gamma_minus_one *
    mj * pj->u * pi->u * (f_ij / pi->pressure_bar) * wi_dr * dvdr * ri;

  pi->u_dt += (u_dt_hydro + u_dt_visc);
}

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* For now, just do the simple case */

  runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj);
  runner_iact_nonsym_force(r2, dx, hj, hi, pj, pi);
}

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_IACT_H */
