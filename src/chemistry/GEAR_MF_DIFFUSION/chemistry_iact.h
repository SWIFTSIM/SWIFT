/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_IACT_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_IACT_H

#include "chemistry_flux.h"
#include "chemistry_getters.h"
#include "chemistry_gradients.h"

#define GIZMO_VOLUME_CORRECTION

/**
 * @file GEAR_MF_DIFFUSION/chemistry_iact.h
 * @brief Diffusion of metals with MFM.
 *
 * The description of the algorithms for diffusion are given in Hopkins 2017
 * (https://arxiv.org/abs/1602.07703) */

/**
 * @brief Do chemistry computation after the runner_iact_density (symmetric
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
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;
  const float r = sqrtf(r2);

  /*****************************************/
  /* Compute the filtered quantities */

  /* Compute the filtered rho = same as rho in SPH but with h_bar instead of h,
     where h_bar = compact support of the kernel */
  float hi_bar = hi * kernel_gamma;
  float hj_bar = hj * kernel_gamma;
  float wi_bar, wj_bar;
  kernel_eval(r / hi_bar, &wi_bar);
  kernel_eval(r / hj_bar, &wj_bar);

  /* j contributes to i and vice-versa */
  chi->filtered.rho += hydro_get_mass(pj) * wi_bar;
  chj->filtered.rho += hydro_get_mass(pi) * wj_bar;

  /* Some smoothing length multiples. */
  float h_bar_ij = 0.5 * (hi_bar + hj_bar);         /* arithmetic mean */
  const float h_inv_bar = 1.0f / h_bar_ij;          /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv_bar); /* 1/h^d */

  /* Take the previous value of bar{rho} since we are computing it now */
  float rho_i_bar = chi->filtered.rho_prev;
  float rho_j_bar = chj->filtered.rho_prev;
  float rho_bar_mean =
      2 * rho_i_bar * rho_j_bar / (rho_i_bar + rho_j_bar); /* harmonic mean */
  float w_filtered;
  kernel_eval(r * h_inv_bar, &w_filtered);

  /* Take the previous value of rho since it is being computed in the density
     loop */
  float rho_i = chi->rho_prev;
  float rho_j = chj->rho_prev;

  /* Avoid 0 division and NaN */
  if (rho_bar_mean != 0 && !isnan(rho_bar_mean)) {
    /* Now compute the filtered rho*v */
    chi->filtered.rho_v[0] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[0] - rho_i * pi->v[0]) *
                              w_filtered * h_inv_dim;
    chi->filtered.rho_v[1] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[1] - rho_i * pi->v[1]) *
                              w_filtered * h_inv_dim;
    chi->filtered.rho_v[2] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[2] - rho_i * pi->v[2]) *
                              w_filtered * h_inv_dim;

    /* Notice the - since the subtraction must be inverted */
    chj->filtered.rho_v[0] -= hydro_get_mass(pi) / rho_bar_mean *
                              (rho_j * pj->v[0] - rho_i * pi->v[0]) *
                              w_filtered * h_inv_dim;
    chj->filtered.rho_v[1] -= hydro_get_mass(pi) / rho_bar_mean *
                              (rho_j * pj->v[1] - rho_i * pi->v[1]) *
                              w_filtered * h_inv_dim;
    chj->filtered.rho_v[2] -= hydro_get_mass(pi) / rho_bar_mean *
                              (rho_j * pj->v[2] - rho_i * pi->v[2]) *
                              w_filtered * h_inv_dim;
  }
}

/**
 * @brief Do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * Compute MFM geometry variables if needed.
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
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;
  const float r = sqrtf(r2);

  /*****************************************/
  /* Compute the filtered quantities */

  /* Compute the filtered rho = same as rho in SPH but with h_bar instead of h,
     where h_bar = compact support of the kernel */
  float hi_bar = hi * kernel_gamma;
  float hj_bar = hj * kernel_gamma;
  float wi_bar;
  kernel_eval(r / hi_bar, &wi_bar);

  /* j contributes to i */
  chi->filtered.rho += hydro_get_mass(pj) * wi_bar;

  /* Some smoothing length multiples. */
  float h_bar_ij = 0.5 * (hi_bar + hj_bar);         /* arithmetic mean */
  const float h_inv_bar = 1.0f / h_bar_ij;          /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv_bar); /* 1/h^d */

  /* Take the previous value of \bar{rho} since we are computing it now */
  float rho_i_bar = chi->filtered.rho_prev;
  float rho_j_bar = chj->filtered.rho_prev;
  float rho_bar_mean =
      2 * rho_i_bar * rho_j_bar / (rho_i_bar + rho_j_bar); /* harmonic mean */
  float w_filtered;
  kernel_eval(r * h_inv_bar, &w_filtered);

  /* Take the previous value of rho since it is being computed in the density
     loop */
  float rho_i = chi->rho_prev;
  float rho_j = chj->rho_prev;

  /* Avoid 0 division and NaN */
  if (rho_bar_mean != 0 && !isnan(rho_bar_mean)) {
    /* Now compute the filtered rho*v */
    chi->filtered.rho_v[0] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[0] - rho_i * pi->v[0]) *
                              w_filtered * h_inv_dim;
    chi->filtered.rho_v[1] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[1] - rho_i * pi->v[1]) *
                              w_filtered * h_inv_dim;
    chi->filtered.rho_v[2] += hydro_get_mass(pj) / rho_bar_mean *
                              (rho_j * pj->v[2] - rho_i * pi->v[2]) *
                              w_filtered * h_inv_dim;
  }
}

/**
 * @brief Do metal diffusion computations in the gradient loop (symmetric
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
__attribute__((always_inline)) INLINE static void
runner_iact_gradient_diffusion(const float r2, const float dx[3],
                               const float hi, const float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, const float a,
                               const float H) {
  chemistry_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Do metal diffusion computations in the gradient loop (nonsymmetric
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
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_gradient_diffusion(const float r2, const float dx[3],
                                      const float hi, const float hj,
                                      struct part *restrict pi,
                                      struct part *restrict pj, const float a,
                                      const float H) {
  chemistry_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Common part of the flux calculations between particle i and j.
 *
 * Since the only difference between the symmetric and non-symmetric version
 * of the flux calculation  is in the update of the conserved variables at the
 * very end (which is not done for particle j if mode is 0), both
 * runner_iact_diffusion() and runner_iact_diffusion() call this method, with
 * an appropriate mode.
 *
 * This method calculates the surface area of the interface between particle i
 * and particle j, as well as the interface position and velocity. These are
 * then used to reconstruct and predict the primitive variables, which are then
 * fed to a Riemann solver that calculates a flux. This flux is used to update
 * the conserved variables of particle i or both particles.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The #cosmology.
 * @param mode 0 if non-symmetric interaction, 1 if symmetric
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_fluxes_common(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, int mode) {

  /* If the masses are null, then there is nothing to diffuse. */
  if (hydro_get_mass(pi) == 0.0 || hydro_get_mass(pj) == 0) {
    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Initialize local variables */
  float Bi[3][3];
  float Bj[3][3];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  const float Vi = pi->geometry.volume;
  const float Vj = pj->geometry.volume;

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute (square of) area */
  /* eqn. (7) */
  float Anorm2 = 0.0f;
  float A[3];
  if (fvpm_part_geometry_well_behaved(pi) &&
      fvpm_part_geometry_well_behaved(pj)) {
    /* in principle, we use Vi and Vj as weights for the left and right
     * contributions to the generalized surface vector.
     * However, if Vi and Vj are very different (because they have very
     * different smoothing lengths), then the expressions below are more
     * stable. */
    float Xi = Vi;
    float Xj = Vj;
#ifdef GIZMO_VOLUME_CORRECTION
    if (fabsf(Vi - Vj) / min(Vi, Vj) > 1.5f * hydro_dimension) {
      Xi = (Vi * hj + Vj * hi) / (hi + hj);
      Xj = Xi;
    }
#endif
    for (int k = 0; k < 3; k++) {
      /* we add a minus sign since dx is pi->x - pj->x */
      A[k] = -Xi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) *
                 wi * hi_inv_dim -
             Xj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) *
                 wj * hj_inv_dim;
      Anorm2 += A[k] * A[k];
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    const float hidp1 = pow_dimension_plus_one(hi_inv);
    const float hjdp1 = pow_dimension_plus_one(hj_inv);
    const float Anorm =
        -(hidp1 * Vi * Vi * wi_dx + hjdp1 * Vj * Vj * wj_dx) * r_inv;
    A[0] = -Anorm * dx[0];
    A[1] = -Anorm * dx[1];
    A[2] = -Anorm * dx[2];
    Anorm2 = Anorm * Anorm * r2;
  }

  /* if the interface has no area, nothing happens and we return */
  /* continuing results in dividing by zero and NaN's... */
  if (Anorm2 == 0.0f) {
    return;
  }

  /* Compute the area */
  const float Anorm_inv = 1.0f / sqrtf(Anorm2);
  const float Anorm = Anorm2 * Anorm_inv;

#ifdef SWIFT_DEBUG_CHECKS
  /* For stability reasons, we do require A and dx to have opposite
   * directions (basically meaning that the surface normal for the surface
   * always points from particle i to particle j, as it would in a real
   * moving-mesh code). If not, our scheme is no longer upwind and hence can
   * become unstable. */
  const float dA_dot_dx = A[0] * dx[0] + A[1] * dx[1] + A[2] * dx[2];
  /* In GIZMO, Phil Hopkins reverts to an SPH integration scheme if this
   * happens. We curently just ignore this case and display a message. */
  const float rdim = pow_dimension(r);
  if (dA_dot_dx > 1.e-6f * rdim) {
    message("Ill conditioned gradient matrix (%g %g %g %g %g)!", dA_dot_dx,
            Anorm, Vi, Vj, r);
  }
#endif

  /* Compute the normal vector of the interface */
  const float n_unit[3] = {A[0] * Anorm_inv, A[1] * Anorm_inv,
                           A[2] * Anorm_inv};

  /* Compute interface position (relative to pi, since we don't need
   * the actual position) eqn. (8) */
  const float xfac = -hi / (hi + hj);
  const float xij_i[3] = {xfac * dx[0], xfac * dx[1], xfac * dx[2]};

  /* Get the time step for the flux exchange. This is always the smallest time
   * step among the two particles. */
  const float mindt =
      (chj->flux_dt > 0.f) ? fminf(chi->flux_dt, chj->flux_dt) : chi->flux_dt;

  /*****************************************/
  /* Predict the velocity at the interface to compute fluxes */
  /* Get the hydro W_L and W_R */
  float vi[3] = {pi->v[0], pi->v[1], pi->v[2]};
  float vj[3] = {pj->v[0], pj->v[1], pj->v[2]};

  /* Compute interface velocity */
  const float vij[3] = {vi[0] + (vi[0] - vj[0]) * xfac,
                        vi[1] + (vi[1] - vj[1]) * xfac,
                        vi[2] + (vi[2] - vj[2]) * xfac};

  /* Get the primitive variable of Euler eq */
  float Wi[5] = {hydro_get_comoving_density(pi), vi[0], vi[1], vi[2],
                 hydro_get_comoving_pressure(pi)};
  float Wj[5] = {hydro_get_comoving_density(pj), vj[0], vj[1], vj[2],
                 hydro_get_comoving_pressure(pj)};

  chemistry_gradients_predict_hydro(pi, pj, dx, r, xij_i, Wi, Wj);

  /* Boost the primitive variables to the frame of reference of the interface */
  /* Note that velocities are indices 1-3 in W */
  /* Note: This is necessary to properly follow the fluid motion. */
  Wi[1] -= vij[0];
  Wi[2] -= vij[1];
  Wi[3] -= vij[2];
  Wj[1] -= vij[0];
  Wj[2] -= vij[1];
  Wj[3] -= vij[2];

  /* Convert to physical units */
  Wi[0] *= cosmo->a3_inv;
  Wi[1] /= cosmo->a;
  Wi[2] /= cosmo->a;
  Wi[3] /= cosmo->a;
  Wi[4] *= cosmo->a_factor_pressure;

  Wj[0] *= cosmo->a3_inv;
  Wj[1] /= cosmo->a;
  Wj[2] /= cosmo->a;
  Wj[3] /= cosmo->a;
  Wj[4] *= cosmo->a_factor_pressure;

  /* Helper variable */
  const float a2 = cosmo->a * cosmo->a;

  /*****************************************/
  /* Now solve the Riemann problem for each metal specie */
  for (int g = 0; g < GEAR_CHEMISTRY_ELEMENT_COUNT; g++) {

    /* Predict the diffusion state at the interface to compute fluxes */
    double Ui, Uj;
    chemistry_gradients_predict(pi, pj, g, dx, r, xij_i, &Ui, &Uj);

    /* Convert Ui and Uj to physical units */
    Ui *= cosmo->a3_inv;
    Uj *= cosmo->a3_inv;

    /* Solve the 1D Riemann problem at the interface A_ij _physical units_ */
    double totflux;
    chemistry_compute_flux(pi, pj, Ui, Uj, Wi, Wj, n_unit, a2 * Anorm, g,
                           &totflux, chem_data, cosmo);

    /* Limit the mass flux to 1/4 of the total mass. This avoids ending with
       more metal mass than the gas_mass. */
    /* if (fabs(totflux * mindt) > 0.0) { */
    /*   const double mi = hydro_get_mass(pi); */
    /*   const double mj = hydro_get_mass(pj); */
    /*   const double Zi = chemistry_get_metal_mass_fraction(pi, g); */
    /*   const double Zj = chemistry_get_metal_mass_fraction(pj, g); */
    /*   const double min_m_Z = min(mi, mj) * fabs(Zi - Zj); */
    /*   const double max_m_Z = max(mi * Zi, mj * Zj); */
    /*   double m_Z_lim = 0.25 * min(min_m_Z, max_m_Z); */

    /*   if (fabs(totflux * mindt) > m_Z_lim) { */
    /*     totflux *= m_Z_lim / fabs(totflux); */
    /*   } */
    /* } */

    /* When solving the Riemann problem, we assume pi is left state, and
     * pj is right state. The sign convention is that a positive total
     * flux is subtracted from the left state, and added to the right
     * state, based on how we chose the unit vector. By this convention,
     * the time integration results in conserved quantity += flux * dt */
    /* Unlike in SPH schemes, we do need to update inactive neighbours, so that
     * the fluxes are always exchanged symmetrically. Thanks to our sneaky use
     * of flux_dt, we can detect inactive neighbours through their negative time
     * step. */
    /* Update V*U. */
    chi->diffusion_flux[g] -= totflux * mindt;
    if (mode == 1 || (chj->flux_dt < 0.f)) {
      chj->diffusion_flux[g] += totflux * mindt;
    }
  }
}

/**
 * @brief Do metal diffusion computation in the <FORCE LOOP> (symmetric version)
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
 * @param chem_data The global properties of the chemistry scheme.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {

  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      1);
}

/**
 * @brief Do metal diffusion computation in the <FORCE LOOP>
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
 * @param chem_data The global properties of the chemistry scheme.
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {

  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      0);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_IACT_H */
