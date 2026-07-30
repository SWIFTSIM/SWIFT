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
#ifndef SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_IACT_H
#define SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_IACT_H

#include "chemistry_flux.h"
#include "chemistry_getters.h"
#include "chemistry_gradients.h"
#include "chemistry_properties.h"

#define GIZMO_VOLUME_CORRECTION

/* Passes of the two-pass FCT positivity limiter: the prep pass computes the
   prospective fluxes and accumulates the donor-side outflow sums, the apply
   pass computes the same fluxes, scales them by the donor's theta and applies
   them. */
#define CHEMISTRY_FCT_PASS_PREP 0
#define CHEMISTRY_FCT_PASS_APPLY 1

/* Status codes of chemistry_gear_fvpm_compute_pair_fluxes(), ordered by how
   far the computation got before an early exit. */
#define CHEMISTRY_PAIR_FLUX_NONE 0     /* no interaction at all */
#define CHEMISTRY_PAIR_FLUX_VMAX 1     /* vmax valid; interface has no area */
#define CHEMISTRY_PAIR_FLUX_GEOMETRY 2 /* vmax+delxbar valid; mindt == 0 */
#define CHEMISTRY_PAIR_FLUX_FLUXES 3   /* fluxes valid */

/**
 * @brief Computes the diffusion pair fluxes of all elements (side-effect
 * free).
 *
 * Defined out-of-line in chemistry.c and deliberately NOT inlined: the FCT
 * prep loop and the flux-applying force loop must obtain bit-identical fluxes
 * for the donor-side positivity cap to hold, which is only guaranteed if both
 * passes execute the exact same machine code.
 *
 * @return a CHEMISTRY_PAIR_FLUX_* status describing which outputs are valid.
 */
int chemistry_gear_fvpm_compute_pair_fluxes(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *restrict pi, const struct part *restrict pj,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo,
    double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4], float *mindt_out,
    float *vmax_out
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
    ,
    float *delxbar_i_out, float *delxbar_j_out
#endif
);

/**
 * @file GEAR_FVPM_DIFFUSION/chemistry_iact.h
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
  float rho_ij_bar = rho_i_bar + rho_j_bar;

  /* Harmonic mean */
  float rho_bar_mean = 2 * rho_i_bar * rho_j_bar / rho_ij_bar;
  float w_filtered;
  kernel_eval(r * h_inv_bar, &w_filtered);

  /* Take the previous value of rho since it is being computed in the density
     loop */
  float rho_i = chi->rho_prev;
  float rho_j = chj->rho_prev;

  /* Avoid 0 division and NaN */
  if (rho_bar_mean > 0.0 && rho_ij_bar > 0.0) {
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
  float rho_ij_bar = rho_i_bar + rho_j_bar;

  /* Harmonic mean */
  float rho_bar_mean = 2 * rho_i_bar * rho_j_bar / rho_ij_bar;
  float w_filtered;
  kernel_eval(r * h_inv_bar, &w_filtered);

  /* Take the previous value of rho since it is being computed in the density
     loop */
  float rho_i = chi->rho_prev;
  float rho_j = chj->rho_prev;

  /* Avoid 0 division and NaN */
  if (rho_bar_mean > 0.0 && rho_ij_bar > 0.0) {
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
 * @param interaction_mode 0 if non-symmetric interaction, 1 if symmetric.
 * @param fct_pass CHEMISTRY_FCT_PASS_PREP to only accumulate the donor-side
 * outflow sums, CHEMISTRY_FCT_PASS_APPLY to scale by the donor's theta and
 * apply the fluxes.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_fluxes_common(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, int interaction_mode, int fct_pass) {

  /* Compute the pair fluxes in the shared, non-inlined helper (see
     chemistry.c). The prep and apply passes live in different translation
     units, so a single out-of-line copy is what makes them provably execute
     identical code under fast-math/LTO, rather than two copies free to round
     differently and flip discontinuous limiter branches. */
  double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4];
  float mindt = 0.f, vmax = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  float delxbar_i = 0.f, delxbar_j = 0.f;
  const int status = chemistry_gear_fvpm_compute_pair_fluxes(
      r2, dx, hi, hj, pi, pj, chem_data, cosmo, totflux, &mindt, &vmax,
      &delxbar_i, &delxbar_j);
#else
  const int status = chemistry_gear_fvpm_compute_pair_fluxes(
      r2, dx, hi, hj, pi, pj, chem_data, cosmo, totflux, &mindt, &vmax);
#endif

  if (status == CHEMISTRY_PAIR_FLUX_NONE) return;

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  /* Store the signal velocity (only in the pass that applies fluxes, to
     leave the time-step estimate identical to the single-pass scheme) */
  if (fct_pass == CHEMISTRY_FCT_PASS_APPLY) {
    chi->timestepvars.vmax = max(chi->timestepvars.vmax, vmax);
    if (interaction_mode == 1) {
      chj->timestepvars.vmax = max(chj->timestepvars.vmax, vmax);
    }
  }
#else
  (void)vmax;
#endif

  if (status == CHEMISTRY_PAIR_FLUX_VMAX) return;

#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  /* Lanson & Vila (2008), denominator of equation (58) */
  if (fct_pass == CHEMISTRY_FCT_PASS_APPLY) {
    chi->timestepvars.delxbar += delxbar_i;
    if (interaction_mode == 1) {
      chj->timestepvars.delxbar += delxbar_j;
    }
  }
#endif

  if (status == CHEMISTRY_PAIR_FLUX_GEOMETRY) return;

  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {

    /* Two-pass FCT positivity limiter. The donor of a pair is the side losing
     * metal mass (positive flux is subtracted from the left state pi). The
     * prep pass accumulates each donor's prospective outflow, but only in the
     * evaluation where the donor's flux field would be updated below, to avoid
     * double counting in the nonsym double evaluation of a pair. The apply
     * pass scales the flux by the donor's theta (frozen since the
     * chemistry_fct_ghost task), so both evaluations of a pair apply the
     * identical scaled flux and the exchange stays conservative. */
    if (fct_pass == CHEMISTRY_FCT_PASS_PREP) {
      /* NOTE: the write to the (possibly inactive) neighbour chj below is
         race-safe only because self/pair hydro tasks lock BOTH cells
         (task_lock default hydro branch). Do not relax the prep task's
         locking to a single cell. */
      if (totflux[m][0] > 0.0) {
        chi->fct_sum_out[m] += totflux[m][0] * mindt;
      } else if (totflux[m][0] < 0.0 &&
                 (interaction_mode == 1 || (chj->flux.dt < 0.f))) {
        chj->fct_sum_out[m] += -totflux[m][0] * mindt;
      }
      /* Do not apply anything in the prep pass */
      continue;
    }

    if (totflux[m][0] > 0.0) {
      totflux[m][0] *= chi->fct_theta[m];
    } else if (totflux[m][0] < 0.0) {
      totflux[m][0] *= chj->fct_theta[m];
    }

    /* Update V*U ****************************************/
    /* When solving the Riemann problem, we assume pi is left state, and
     * pj is right state. The sign convention is that a positive total
     * flux is subtracted from the left state, and added to the right
     * state, based on how we chose the unit vector. By this convention,
     * the time integration results in conserved quantity += flux * dt */
    /* Unlike in SPH schemes, we do need to update inactive neighbours, so that
     * the fluxes are always exchanged symmetrically. Thanks to our sneaky use
     * of flux.dt, we can detect inactive neighbours through their negative time
     * step. */
    chemistry_part_update_fluxes_left(pi, m, totflux[m], mindt);
    if (interaction_mode == 1 || (chj->flux.dt < 0.f)) {
      chemistry_part_update_fluxes_right(pj, m, totflux[m], mindt);
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
                                      1, CHEMISTRY_FCT_PASS_APPLY);
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

#if !defined(GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE)
  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      0, CHEMISTRY_FCT_PASS_APPLY);
#else
  int local_mode = 0;
  const int pi_is_active = pi->chemistry_data.flux.dt > 0.f;
  const int pj_is_active = pj->chemistry_data.flux.dt > 0.f;

  if (pi_is_active && pj_is_active) {
    if (pi->id < pj->id) {
      /* pi has the duty to update both pi and pj symmetrically */
      local_mode = 1;
    } else {
      /* Nothing to do. The other particle will update symmetrically */
      return;
    }
  } else {
    /* pi or pj is active. No trick is needed here: only the active particle
       will update both particles. So we can use the default implementation. */
    local_mode = 0;
  }
  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      local_mode, CHEMISTRY_FCT_PASS_APPLY);
#endif
}

/**
 * @brief Prep pass of the FCT positivity limiter (symmetric version).
 *
 * Runs in the chemistry_fct_prep loop, between the extra ghost and the
 * chemistry_fct_ghost task. Computes the same fluxes as the flux exchange in
 * the force loop, but only accumulates the donor-side outflow sums; nothing
 * is applied. The signature matches runner_iact_diffusion() so the fct_prep
 * loop instantiation can alias one to the other.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_diffusion_fct_prep(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {

  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      1, CHEMISTRY_FCT_PASS_PREP);
}

/**
 * @brief Prep pass of the FCT positivity limiter (nonsymmetric version).
 *
 * Mirrors runner_iact_nonsym_diffusion() exactly (including the debug
 * one-sided update mode) so the prep pass sees the identical set of pair
 * evaluations as the apply pass in the force loop.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_diffusion_fct_prep(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {

#if !defined(GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE)
  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      0, CHEMISTRY_FCT_PASS_PREP);
#else
  int local_mode = 0;
  const int pi_is_active = pi->chemistry_data.flux.dt > 0.f;
  const int pj_is_active = pj->chemistry_data.flux.dt > 0.f;

  if (pi_is_active && pj_is_active) {
    if (pi->id < pj->id) {
      local_mode = 1;
    } else {
      return;
    }
  } else {
    local_mode = 0;
  }
  runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data, cosmo,
                                      local_mode, CHEMISTRY_FCT_PASS_PREP);
#endif
}

#endif /* SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_IACT_H */
