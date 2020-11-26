/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_DIFFUSION_CHEMISTRY_IACT_H
#define SWIFT_GEAR_DIFFUSION_CHEMISTRY_IACT_H

/**
 * @file GEAR/chemistry_iact.h
 * @brief Smooth metal interaction + diffusion functions following the GEAR
 * version of smooth metalicity.
 */

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);
  const float mi = hydro_get_mass(pi);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_deval(uj, &wj, &wj_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  const float wj_dr = wj_dx * r_inv;
  const float mi_wj_dr = mi * wj_dr;

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chi->smoothed_metal_mass_fraction[i] += chj->metal_mass[i] * wi;
    chj->smoothed_metal_mass_fraction[i] += chi->metal_mass[i] * wj;
  }

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    chi->S[i][0] += (pj->v[0] - pi->v[0]) * dx[i] * mj_wi_dr;
    chi->S[i][1] += (pj->v[1] - pi->v[1]) * dx[i] * mj_wi_dr;
    chi->S[i][2] += (pj->v[2] - pi->v[2]) * dx[i] * mj_wi_dr;

    chj->S[i][0] -= (pj->v[0] - pi->v[0]) * dx[i] * mi_wj_dr;
    chj->S[i][1] -= (pj->v[1] - pi->v[1]) * dx[i] * mi_wj_dr;
    chj->S[i][2] -= (pj->v[2] - pi->v[2]) * dx[i] * mi_wj_dr;
  }
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  /* Compute contribution to the smooth metallicity */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chi->smoothed_metal_mass_fraction[i] += chj->metal_mass[i] * wi;
  }

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    chi->S[i][0] += (pj->v[0] - pi->v[0]) * dx[i] * mj_wi_dr;
    chi->S[i][1] += (pj->v[1] - pi->v[1]) * dx[i] * mj_wi_dr;
    chi->S[i][2] += (pj->v[2] - pi->v[2]) * dx[i] * mj_wi_dr;
  }
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, float time_base,
    integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* No need to diffuse if both particles are not diffusing. */
  if (chj->diff_coef > 0 && chi->diff_coef > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float mi = hydro_get_mass(pi);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);

    const float wi_dr = dwi_dx * r_inv;
    const float wj_dr = dwj_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;
    const float mi_dw_r = mi * wj_dr;

    /* Compute the diffusion coefficient <D> / <rho> in physical units. */
    const float coef = 2. * (chi->diff_coef + chj->diff_coef) / (rhoi + rhoj);

    const float coef_i = coef * mj_dw_r;
    const float coef_j = coef * mi_dw_r;

    /* Compute the time derivative */
    const float inv_mi = 1.f / hydro_get_mass(pi);
    const float inv_mj = 1.f / hydro_get_mass(pj);
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      const double mi_frac = chi->metal_mass[i] * inv_mi;
      const double mj_frac = chj->metal_mass[i] * inv_mj;
      const double dm = mi_frac - mj_frac;
      chi->metal_mass_dt[i] += coef_i * dm;
      chj->metal_mass_dt[i] -= coef_j * dm;
    }
  }
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, float time_base,
    integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diff_coef > 0 && chi->diff_coef > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, dwi_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);
    const float wi_dr = dwi_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;

    /* Compute the diffusion coefficient <D> / <rho> in physical units. */
    const float coef = 2. * (chi->diff_coef + chj->diff_coef) / (rhoi + rhoj);

    const float coef_i = coef * mj_dw_r;

    /* Compute the time derivative */
    const float inv_mi = 1.f / hydro_get_mass(pi);
    const float inv_mj = 1.f / hydro_get_mass(pj);
    for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
      const double mi_frac = chi->metal_mass[i] * inv_mi;
      const double mj_frac = chj->metal_mass[i] * inv_mj;
      const double dm = mi_frac - mj_frac;
      chi->metal_mass_dt[i] += coef_i * dm;
    }
  }
}

#endif /* SWIFT_GEAR_DIFFUSION_CHEMISTRY_IACT_H */
