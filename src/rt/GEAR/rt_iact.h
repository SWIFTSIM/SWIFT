/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IACT_GEAR_H
#define SWIFT_RT_IACT_GEAR_H

#include "rt_debugging.h"
#include "rt_flux.h"
#include "rt_gradients.h"

/**
 * @file src/rt/GEAR/rt_iact.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * particle interactions.
 */

/**
 * @brief Preparation step for injection to gather necessary data.
 * This function gets called during the feedback force loop.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle (not updated).
 * @param cosmo The cosmological model.
 * @param rt_props Properties of the RT scheme.
 */

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_injection_prep(const float r2, const float *dx,
                                     const float hi, const float hj,
                                     struct spart *si, const struct part *pj,
                                     const struct cosmology *cosmo,
                                     const struct rt_props *rt_props) {

  /* If the star doesn't have any neighbours, we
   * have nothing to do here. */
  if (si->density.wcount == 0.f) return;

#ifdef SWIFT_RT_DEBUG_CHECKS
  si->rt_data.debug_iact_hydro_inject_prep += 1;
#endif

  /* Compute the weight of the neighbouring particle */
  const float hi_inv = 1.f / hi;
  const float r = sqrtf(r2);
  const float xi = r * hi_inv;
  float wi;
  kernel_eval(xi, &wi);
  const float hi_inv_dim = pow_dimension(hi_inv);
  /* psi(x_star - x_gas, h_star) */
  /* Note: skip the devision by si->density.wcount here. It'll cancel out by the
   * normalization anyway, and furthermore now that the injection prep is done
   * during the star density loop, si->density.wcount won't be computed at this
   * stage yet. */
  const float psi = wi * hi_inv_dim;

  /* Now add that weight to the appropriate octant */
  int octant_index = 0;

  if (dx[0] > 0.f) octant_index += 1;
  if (dx[1] > 0.f) octant_index += 2;
  if (dx[2] > 0.f) octant_index += 4;

  si->rt_data.octant_weights[octant_index] += psi;
}

/**
 * @brief Injection step interaction between star and hydro particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si Star particle.
 * @param pj Hydro particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param rt_props Properties of the RT scheme.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_inject(
    const float r2, float *dx, const float hi, const float hj,
    struct spart *restrict si, struct part *restrict pj, float a, float H,
    const struct rt_props *rt_props) {

  /* If the star doesn't have any neighbours, we
   * have nothing to do here. */
  if (si->density.wcount == 0.f) return;

#ifdef SWIFT_RT_DEBUG_CHECKS

  /* Do some checks and increase neighbour counts
   * before other potential early exits */
  if (si->rt_data.debug_iact_hydro_inject_prep == 0)
    error(
        "Injecting energy from star that wasn't called during injection prep");

  if (!si->rt_data.debug_emission_rate_set)
    error("Injecting energy from star without setting emission rate");

  si->rt_data.debug_iact_hydro_inject += 1;
  si->rt_data.debug_radiation_emitted_tot += 1ULL;

  pj->rt_data.debug_iact_stars_inject += 1;
  pj->rt_data.debug_radiation_absorbed_tot += 1ULL;

#endif

  /* Compute the weight of the neighbouring particle */
  const float hi_inv = 1.f / hi;
  const float r = sqrtf(r2);
  const float xi = r * hi_inv;
  float wi;
  kernel_eval(xi, &wi);
  const float hi_inv_dim = pow_dimension(hi_inv);
  /* psi(x_star - x_gas, h_star) */
  /* Skip the division by si->density.wcount to remain consistent */
  const float psi = wi * hi_inv_dim;

#if defined(HYDRO_DIMENSION_3D)
  const int maxind = 8;
#elif defined(HYDRO_DIMENSION_2D)
  const int maxind = 4;
#elif defined(HYDRO_DIMENSION_1D)
  const int maxind = 2;
#endif

  /* Get weight for particle, including isotropy correction */
  float nonempty_octants = 0.f;

  for (int i = 0; i < maxind; i++) {
    if (si->rt_data.octant_weights[i] > 0.f) nonempty_octants += 1.f;
  }

  int octant_index = 0;
  if (dx[0] > 0.f) octant_index += 1;
  if (dx[1] > 0.f) octant_index += 2;
  if (dx[2] > 0.f) octant_index += 4;

  const float octw = si->rt_data.octant_weights[octant_index];
  /* We might end up in this scenario due to roundoff errors */
  if (psi == 0.f || octw == 0.f) return;

  const float weight = psi / (nonempty_octants * octw);
  const float Vinv = 1.f / pj->geometry.volume;

  /* Nurse, the patient is ready now */
  for (int g = 0; g < RT_NGROUPS; g++) {
    /* Inject energy. */
    const float injected_energy_density =
        si->rt_data.emission_this_step[g] * weight * Vinv;
    pj->rt_data.radiation[g].energy_density += injected_energy_density;

    /* Don't inject flux. */
  }

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Take note how much energy we actually injected */
  for (int g = 0; g < RT_NGROUPS; g++) {
    const float injected_energy = si->rt_data.emission_this_step[g] * weight;
    if (isinf(injected_energy) || isnan(injected_energy))
      error(
          "Injecting abnormal energy spart %lld part %lld group %d | %.6e %.6e "
          "%.6e",
          si->id, pj->id, g, injected_energy, weight,

          si->rt_data.emission_this_step[g]);
    si->rt_data.debug_injected_energy[g] += injected_energy;
    si->rt_data.debug_injected_energy_tot[g] += injected_energy;
  }
  si->rt_data.debug_psi_sum += psi;
#endif
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param mode 0 if non-symmetric interaction, 1 if symmetric
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_flux_common(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H, int mode) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  const char *func_name = (mode == 1) ? "sym flux iact" : "nonsym flux iact";
  rt_debug_sequence_check(pi, 3, func_name);
  pi->rt_data.debug_calls_iact_transport_interaction += 1;

  if (mode == 1) {
    rt_debug_sequence_check(pj, 3, func_name);
    pj->rt_data.debug_calls_iact_transport_interaction += 1;
  }
#endif

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
  if (hydro_part_geometry_well_behaved(pi) &&
      hydro_part_geometry_well_behaved(pj)) {
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

  /* compute the normal vector of the interface */
  const float n_unit[3] = {A[0] * Anorm_inv, A[1] * Anorm_inv,
                           A[2] * Anorm_inv};

  /* Compute interface position (relative to pi, since we don't need
   * the actual position) eqn. (8) */
  const float xfac = -hi / (hi + hj);
  const float xij_i[3] = {xfac * dx[0], xfac * dx[1], xfac * dx[2]};

  struct rt_part_data *restrict rti = &pi->rt_data;
  struct rt_part_data *restrict rtj = &pj->rt_data;
  /* Get the time step for the flux exchange. This is always the smallest time
   * step among the two particles. */
  const float mindt =
      (rtj->flux_dt > 0.f) ? fminf(rti->flux_dt, rtj->flux_dt) : rti->flux_dt;

  for (int g = 0; g < RT_NGROUPS; g++) {

    /* radiation state to be used to compute the flux */
    float Ui[4], Uj[4];
    rt_gradients_predict(pi, pj, Ui, Uj, g, dx, r, xij_i);

    /* For first order method, skip the gradients */
    /* float Ui[4], Uj[4]; */
    /* rt_part_get_radiation_state_vector(pi, g, Ui); */
    /* rt_part_get_radiation_state_vector(pj, g, Uj); */
    /* No need to check for unphysical quantities, they
     * haven't been touched since
     * rt_injection_update_photon_densities */

    float totflux[4];

    rt_compute_flux(Ui, Uj, n_unit, Anorm, totflux);

    /* When solving the Riemann problem, we assume pi is left state, and
     * pj is right state. The sign convention is that a positive total
     * flux is subtracted from the left state, and added to the right
     * state, based on how we chose the unit vector. By this convention,
     * the time integration results in conserved quantity += flux * dt */
    /* Unlike in SPH schemes, we do need to update inactive neighbours, so that
     * the fluxes are always exchanged symmetrically. Thanks to our sneaky use
     * of flux_dt, we can detect inactive neighbours through their negative time
     * step. */
    rti->flux[g].energy -= totflux[0] * mindt;
    rti->flux[g].flux[0] -= totflux[1] * mindt;
    rti->flux[g].flux[1] -= totflux[2] * mindt;
    rti->flux[g].flux[2] -= totflux[3] * mindt;
    if (mode == 1 || (rtj->flux_dt < 0.f)) {
      rtj->flux[g].energy += totflux[0] * mindt;
      rtj->flux[g].flux[0] += totflux[1] * mindt;
      rtj->flux[g].flux[1] += totflux[2] * mindt;
      rtj->flux[g].flux[2] += totflux[3] * mindt;
    }
  }
}

/**
 * @brief Flux calculation between particle i and particle j
 *
 * This method calls runner_iact_rt_fluxes_common with mode 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_transport(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 1);
}

/**
 * @brief Flux calculation between particle i and particle j: non-symmetric
 * version
 *
 * This method calls runner_iact_rt_fluxes_common with mode 0.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_transport(float r2, const float *dx, float hi, float hj,
                                struct part *restrict pi,
                                struct part *restrict pj, float a, float H) {

  runner_iact_rt_flux_common(r2, dx, hi, hj, pi, pj, a, H, 0);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around rt_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_rt_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  rt_gradients_collect(r2, dx, hi, hj, pi, pj);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around rt_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_rt_gradient(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj, float a, float H) {

  rt_gradients_nonsym_collect(r2, dx, hi, hj, pi, pj);
}

#endif /* SWIFT_RT_IACT_GEAR_H */
