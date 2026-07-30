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

#include "atomic.h"
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
 * @param update_left Whether pi receives the flux.
 * @param update_right Whether pj receives the flux.
 * @param accum_ts_left Whether pi accumulates this pair's timestep variables
 * (vmax, delxbar). Not the same mask as update_left: the nonsym contract
 * gives the calling wrapper's active particle the timestep-variable credit
 * even on an evaluation where the OTHER (inactive) particle is the one
 * receiving the flux update (see runner_iact_nonsym_diffusion()).
 * @param accum_ts_right Whether pj accumulates this pair's timestep
 * variables. See accum_ts_left.
 * @param fct_pass CHEMISTRY_FCT_PASS_PREP to only accumulate the donor-side
 * outflow sums, CHEMISTRY_FCT_PASS_APPLY to scale by the donor's theta and
 * apply the fluxes.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_fluxes_common(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, int update_left, int update_right,
    int accum_ts_left, int accum_ts_right, int fct_pass) {

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
    if (accum_ts_left)
      chi->timestepvars.vmax = max(chi->timestepvars.vmax, vmax);
    if (accum_ts_right)
      chj->timestepvars.vmax = max(chj->timestepvars.vmax, vmax);
  }
#else
  (void)vmax;
#endif

  if (status == CHEMISTRY_PAIR_FLUX_VMAX) return;

#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  /* Lanson & Vila (2008), denominator of equation (58) */
  if (fct_pass == CHEMISTRY_FCT_PASS_APPLY) {
    if (accum_ts_left) chi->timestepvars.delxbar += delxbar_i;
    if (accum_ts_right) chj->timestepvars.delxbar += delxbar_j;
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
      if (totflux[m][0] > 0.0 && update_left) {
        chi->fct_sum_out[m] += totflux[m][0] * mindt;
      } else if (totflux[m][0] < 0.0 && update_right) {
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
    if (update_left)
      chemistry_part_update_fluxes_left(pi, m, totflux[m], mindt);
    if (update_right)
      chemistry_part_update_fluxes_right(pj, m, totflux[m], mindt);
  }
}

#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT
/* Step-1 instrumentation: how often the local mixed-band tie-break arm
   below actually decides to act versus skip. Defined once in chemistry.c
   (this header may be textually included into more than one translation
   unit); incremented atomically since the force loop runs multithreaded.
   This aggregate pair of counters CANNOT prove "no dropped pairs" on its
   own: an arm-1 (both_updatable_here) visit fully covers a pair without
   touching either counter, so acted > skipped is also consistent with a
   pair being dropped here but covered elsewhere by arm 1. The per-pair log
   below is the instrument that actually discriminates; these two counters
   are kept only as a cheap sanity readout. */
extern long long chemistry_fvpm_visit_acted;
extern long long chemistry_fvpm_visit_skipped;

/* Per-pair discriminating log (see chemistry.c for the full rationale and
   the post-processing that turns this into dropped/double_applied counts).
   NOT thread-safe: diagnostic use only, with --threads=1. */
enum chemistry_fvpm_visit_arm {
  CHEMISTRY_FVPM_VISIT_SYM = 0,  /* arm 1: both_updatable_here */
  CHEMISTRY_FVPM_VISIT_ACT = 1,  /* arm 3: mixed-band owner acts */
  CHEMISTRY_FVPM_VISIT_SKIP = 2, /* arm 3: mixed-band non-owner skips */
};

struct chemistry_fvpm_visit_record {
  const struct part *a; /* min(pi, pj) by address */
  const struct part *b; /* max(pi, pj) by address */
  enum chemistry_fvpm_visit_arm arm;
};

extern void chemistry_fvpm_visit_log_record(const struct part *pi,
                                            const struct part *pj,
                                            enum chemistry_fvpm_visit_arm arm);
extern void chemistry_fvpm_visit_log_analyze(void);
extern void chemistry_fvpm_visit_log_reset(void);
#endif

/**
 * @brief Strict lexicographic order on absolute (wrapped) particle
 * positions.
 *
 * Frame-independent (unlike the sid-shifted local coordinates used inside
 * the sorted pair loops) and rank-independent (proxy particles carry the
 * owner's absolute position verbatim), so two independent evaluations of
 * the same physical pair -- one per task that visits it -- agree on which
 * side "wins" without any shared mutable state. A tie (all three
 * coordinates bitwise equal) means the two particles are coincident; the
 * caller treats that as a signal to skip, since a coincident pair's FVPM
 * face has zero area anyway (Anorm2 == 0, see chemistry.c).
 */
__attribute__((always_inline)) INLINE static int chemistry_pos_less(
    const struct part *restrict a, const struct part *restrict b) {
  if (a->x[0] != b->x[0]) return a->x[0] < b->x[0];
  if (a->x[1] != b->x[1]) return a->x[1] < b->x[1];
  return a->x[2] < b->x[2];
}

/**
 * @brief Designate the owner of a locally-active mixed-depth-band pair: the
 * larger-h particle, with chemistry_pos_less() as an h-tie fallback.
 *
 * This is not an arbitrary tie-break: a pair interacts at all only if
 * r2 < gamma^2 * max(hi, hj)^2, so the partner always lies inside the
 * larger-h particle's own kernel, and SWIFT's hydro loop machinery must
 * enumerate every active particle's full own-kernel neighbour set for
 * density/force correctness. The larger-h side of an active mixed-band pair
 * is therefore GUARANTEED its own visit to this pair; designating it the
 * owner means the acting evaluation is always guaranteed to happen,
 * independent of whether the smaller-h side's visit does. Position-only
 * designation lacked this guarantee (the designated side could be the one
 * whose own kernel does not contain the pair) and measurably dropped a
 * small fraction of mixed-band pairs (Step-1 instrumentation, first
 * implementation of this session). h is per-particle, frozen for the step
 * after the density loop, communicated to proxies, and identical in the FCT
 * prep and apply passes -- the same stability properties chemistry_pos_less
 * relies on for x. Bitwise-equal h implies equal depth_h implies the same
 * band, which the both_updatable_here arm already covers -- so the
 * chemistry_pos_less() fallback here is expected to be dead in practice;
 * kept as a safety net.
 */
__attribute__((always_inline)) INLINE static int chemistry_flux_owner(
    const float hi, const float hj, const struct part *restrict pi,
    const struct part *restrict pj) {
  if (hi != hj) return hi > hj;
  return chemistry_pos_less(pi, pj);
}

/**
 * @brief Shared decision table for the chemistry diffusion flux exchange,
 * used by both the flux-applying force loop and (via fct_pass) the FCT prep
 * pass.
 *
 * The contract at every call site: @p pi is the particle this call is
 * entitled to update (the pre-existing nonsym convention); @p pj is its
 * neighbour. @p both_updatable_here is the call site's own doi&&doj-style
 * gate -- when true, both particles are guaranteed active, local and in the
 * same depth band, and this is the pair's ONE symmetric evaluation for the
 * step (established structural fact, see the plan). @p local_first and
 * @p local_second report whether the cell that @p pi (resp. @p pj) belongs
 * to at THIS call site is local to this rank -- named by argument position,
 * not by any fixed "cell i / cell j" identity, since several call sites
 * place the cell-j particle first.
 *
 * Decision table (order matters):
 *   1. both_updatable_here      -> symmetric update, today's sym path,
 *                                  no reordering.
 *   2. pj inactive (flux.dt<0)  -> nonsym + sentinel, today's path
 *                                  (pj still receives the flux so the
 *                                  exchange stays symmetric; the pair is
 *                                  visited exactly once this step).
 *   3. local_first && local_second (mixed depth band, both active, both
 *      local) -> h-primary ownership: chemistry_flux_owner(hi, hj, pi, pj)
 *      decides which of the two independent tasks visiting this pair
 *      performs the one symmetric update (the larger-h side, whose own
 *      visit to this pair is guaranteed by the hydro enumeration
 *      invariant -- see chemistry_flux_owner()); the loser skips entirely
 *      (an h-and-position tie -- coincident particles -- also skips, on
 *      both sides, consistently).
 *   4. otherwise (cross-rank: pj is a foreign/proxy copy) -> compute the
 *      flux in position-canonical input order (so both ranks evaluate the
 *      identical floating-point sequence) but apply the result only to the
 *      local particle; a tie skips (degenerate face) rather than picking
 *      an arbitrary side.
 *
 * A separate, purely diagnostic bypass
 * (GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE) reproduces the id-keyed
 * one-sided update used to establish the exact-zero conservation reference this
 * session's fix targets; it is NOT MPI-safe (see chemistry_properties.h) and
 * plays no part in the decision table above.
 *
 * @param fct_pass CHEMISTRY_FCT_PASS_PREP or CHEMISTRY_FCT_PASS_APPLY.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_flux_exchange_decide(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const int both_updatable_here, const int local_first,
    const int local_second, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo, const int fct_pass) {

  if (both_updatable_here) {
    /* Today's sym path: both particles already guaranteed active, local and
       in-band by the caller's own gate. No reordering -- see the plan for
       why reordering here would only add rounding noise for no benefit. */
#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT
    /* Logged only on the APPLY pass so each pair contributes one record per
       step (PREP visits the same pairs and would otherwise double them). */
    if (fct_pass == CHEMISTRY_FCT_PASS_APPLY)
      chemistry_fvpm_visit_log_record(pi, pj, CHEMISTRY_FVPM_VISIT_SYM);
#endif
    runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                        cosmo, /*update_left=*/1,
                                        /*update_right=*/1,
                                        /*accum_ts_left=*/1,
                                        /*accum_ts_right=*/1, fct_pass);
    return;
  }

#if defined(GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE)
  /* Reference-only bypass: id-keyed one-sided update, kept verbatim
     (including its id dependence) for the exact-zero conservation
     comparison run. Not part of the production decision table. */
  {
    const int pi_is_active = pi->chemistry_data.flux.dt > 0.f;
    const int pj_is_active = pj->chemistry_data.flux.dt > 0.f;
    if (pi_is_active && pj_is_active) {
      if (pi->id < pj->id) {
        runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                            cosmo, 1, 1, 1, 1, fct_pass);
      }
      /* else: the mirror call (from pj's own active dispatch) performs the
         full symmetric update; nothing to do here. */
    } else {
      const int pj_inactive = (pj->chemistry_data.flux.dt < 0.f);
      runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                          cosmo, 1, pj_inactive, 1, 0,
                                          fct_pass);
    }
  }
  return;
#else

  const int pj_inactive = (pj->chemistry_data.flux.dt < 0.f);
  if (pj_inactive) {
    /* Today's nonsym path: pj still receives the flux so the exchange stays
       symmetric, but this pair is visited exactly once this step (from
       pi's active side), so it gets none of the timestep-variable credit. */
    runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                        cosmo, /*update_left=*/1,
                                        /*update_right=*/1,
                                        /*accum_ts_left=*/1,
                                        /*accum_ts_right=*/0, fct_pass);
    return;
  }

  if (local_first && local_second) {
    /* Both particles active and local, but this call site alone does not
       cover the pair symmetrically (mixed depth band): exactly one of the
       two independent tasks that visit this pair must act -- the larger-h
       side (see chemistry_flux_owner()), whose own visit is guaranteed. */
    if (chemistry_flux_owner(hi, hj, pi, pj)) {
#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT
      atomic_add(&chemistry_fvpm_visit_acted, 1);
      if (fct_pass == CHEMISTRY_FCT_PASS_APPLY)
        chemistry_fvpm_visit_log_record(pi, pj, CHEMISTRY_FVPM_VISIT_ACT);
#endif
      runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                          cosmo, 1, 1, 1, 1, fct_pass);
    } else {
#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT
      atomic_add(&chemistry_fvpm_visit_skipped, 1);
      if (fct_pass == CHEMISTRY_FCT_PASS_APPLY)
        chemistry_fvpm_visit_log_record(pi, pj, CHEMISTRY_FVPM_VISIT_SKIP);
#endif
#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
      if (hi == hj && pi->x[0] == pj->x[0] && pi->x[1] == pj->x[1] &&
          pi->x[2] == pj->x[2])
        warning(
            "Coincident distinct particles (ids %lld/%lld) at pos "
            "[%g %g %g] in the local mixed-band chemistry flux exchange -- "
            "skipping degenerate face.",
            pi->id, pj->id, pi->x[0], pi->x[1], pi->x[2]);
#endif
      /* Loser: the mirror task's own evaluation of this pair wins the
         comparison the other way and acts instead. */
    }
    return;
  }

  /* Cross-rank: local_second is false here (the local_first && local_second
     case was handled above), and every call site's first argument is by
     construction the one this call is entitled to update -- so pi is local
     and pj is a foreign/proxy copy that must never be written. Compute the
     flux in position-canonical input order so both ranks evaluate the
     identical floating-point sequence, then apply the result only to
     whichever slot pi ends up in (never pj). */
#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  if (!local_first)
    error(
        "runner_iact_chemistry_flux_exchange's cross-rank arm reached with "
        "local_first false: the call site's first argument must always be "
        "the locally-owned particle.");
#endif
  const int pi_less_pj = chemistry_pos_less(pi, pj);
  const int pj_less_pi = chemistry_pos_less(pj, pi);
  if (!pi_less_pj && !pj_less_pi) {
#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
    warning(
        "Coincident distinct particles (ids %lld/%lld) at pos [%g %g %g] "
        "in the cross-rank chemistry flux exchange -- skipping degenerate "
        "face.",
        pi->id, pj->id, pi->x[0], pi->x[1], pi->x[2]);
#endif
    return;
  }
  if (pi_less_pj) {
    runner_iact_chemistry_fluxes_common(r2, dx, hi, hj, pi, pj, chem_data,
                                        cosmo, /*update_left=*/1,
                                        /*update_right=*/0,
                                        /*accum_ts_left=*/1,
                                        /*accum_ts_right=*/0, fct_pass);
  } else {
    const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
    runner_iact_chemistry_fluxes_common(r2, mdx, hj, hi, pj, pi, chem_data,
                                        cosmo, /*update_left=*/0,
                                        /*update_right=*/1,
                                        /*accum_ts_left=*/0,
                                        /*accum_ts_right=*/1, fct_pass);
  }
#endif /* !GEAR_FVPM_DIFF_DEBUG_FORCE_LOOP_ONESIDED_UPDATE */
}

/**
 * @brief Do metal diffusion computation in the <FORCE LOOP>: single
 * dispatcher covering the symmetric, nonsym, mixed-band and cross-rank
 * cases (see runner_iact_chemistry_flux_exchange_decide()).
 *
 * @param pi The particle this call site is entitled to update.
 * @param pj Its neighbour.
 * @param both_updatable_here Whether this call site's own gate already
 * guarantees both particles active, local and in the same depth band (the
 * pair's one symmetric evaluation for the step).
 * @param local_first Whether the cell @p pi belongs to is local to this
 * rank.
 * @param local_second Whether the cell @p pj belongs to is local to this
 * rank.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param t_current The current time (in integer).
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 * @param chem_data The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_flux_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const int both_updatable_here, const int local_first,
    const int local_second, const float a, const float H, const float time_base,
    const integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology, const struct chemistry_global_data *chem_data) {

  (void)a;
  (void)H;
  (void)time_base;
  (void)t_current;
  (void)with_cosmology;

  runner_iact_chemistry_flux_exchange_decide(
      r2, dx, hi, hj, pi, pj, both_updatable_here, local_first, local_second,
      chem_data, cosmo, CHEMISTRY_FCT_PASS_APPLY);
}

/**
 * @brief Prep pass of the FCT positivity limiter: mirrors
 * runner_iact_chemistry_flux_exchange() exactly (identical decision table,
 * identical position tie-break) so the prep pass sees the identical set of
 * pair evaluations, in the identical order, as the apply pass in the force
 * loop. The signature matches runner_iact_chemistry_flux_exchange() so the
 * fct_prep loop instantiation's object-macro redirect (runner_doiact_hydro.c)
 * can alias one to the other.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_flux_exchange_fct_prep(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const int both_updatable_here, const int local_first,
    const int local_second, const float a, const float H, const float time_base,
    const integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology, const struct chemistry_global_data *chem_data) {

  (void)a;
  (void)H;
  (void)time_base;
  (void)t_current;
  (void)with_cosmology;

  runner_iact_chemistry_flux_exchange_decide(
      r2, dx, hi, hj, pi, pj, both_updatable_here, local_first, local_second,
      chem_data, cosmo, CHEMISTRY_FCT_PASS_PREP);
}

/* The template's original diffusion hooks: kept (not removed) so private
   out-of-tree schemes sharing this template's call sites are unaffected,
   but empty here -- this scheme's logic now lives entirely in
   runner_iact_chemistry_flux_exchange() above. */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {}

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {}

__attribute__((always_inline)) INLINE static void
runner_iact_diffusion_fct_prep(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_diffusion_fct_prep(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data) {}

#endif /* SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_IACT_H */
