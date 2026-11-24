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

#include "chemistry_debug.h"
#include "chemistry_flux.h"
#include "chemistry_getters.h"
#include "chemistry_gradients.h"
#include "chemistry_properties.h"
#include "flux_exchange.h"

#define GIZMO_VOLUME_CORRECTION

/* Status codes of chemistry_compute_pair_fluxes(), ordered by how
   far the computation got before an early exit. */
#define CHEMISTRY_PAIR_FLUX_NONE 0     /* no interaction at all */
#define CHEMISTRY_PAIR_FLUX_VMAX 1     /* vmax valid; interface has no area */
#define CHEMISTRY_PAIR_FLUX_GEOMETRY 2 /* vmax+delxbar valid; flux_dt == 0 */
#define CHEMISTRY_PAIR_FLUX_FLUXES 3   /* fluxes valid */

/* Defined in chemistry.c; see there for the full doc. */
int chemistry_compute_pair_fluxes(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct part *restrict pi, const struct part *restrict pj,
    const double time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data,
    double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4], float *flux_dt_out,
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
 * @brief Add the outflow of a pair to the sum of its donor, in the prep pass
 * of the positivity limiter.
 *
 * The donor is the side that loses metal mass: a positive flux is taken from
 * pi. Only a side that this call updates is credited, so a pair that is
 * visited twice is not counted twice.
 *
 * @param chi The chemistry data of particle i.
 * @param chj The chemistry data of particle j.
 * @param totflux The fluxes of all elements across the interface.
 * @param flux_dt The time interval of the exchange.
 * @param update_left Whether this call updates pi.
 * @param update_right Whether this call updates pj.
 */
__attribute__((always_inline)) INLINE static void
chemistry_accumulate_donor_outflow(
    struct chemistry_part_data *restrict chi,
    struct chemistry_part_data *restrict chj,
    double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4], const float flux_dt,
    const int update_left, const int update_right) {

  /* The write to a possibly inactive neighbour is race-free only because the
     self and pair hydro tasks lock both cells. Do not relax the locking of
     the prep task to a single cell. */
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    if (totflux[m][0] > 0.0 && update_left) {
      chi->fct_sum_out[m] += totflux[m][0] * flux_dt;
    } else if (totflux[m][0] < 0.0 && update_right) {
      chj->fct_sum_out[m] += -totflux[m][0] * flux_dt;
    }
  }
}

/**
 * @brief Add the timestep variables of a pair (signal velocity, particle
 * size), in the apply pass.
 *
 * Only the apply pass adds them, so the time-step estimate is the same as
 * without the limiter. Their masks differ from the flux masks: when an
 * inactive neighbour receives the flux, only the active particle gets this
 * credit.
 *
 * @param chi The chemistry data of particle i.
 * @param chj The chemistry data of particle j.
 * @param status The CHEMISTRY_PAIR_FLUX_* status of the flux computation.
 * @param vmax The signal velocity of the pair.
 * @param delxbar_i The particle-size term of particle i.
 * @param delxbar_j The particle-size term of particle j.
 * @param accum_ts_left Whether pi gets the timestep variables.
 * @param accum_ts_right Whether pj gets the timestep variables.
 */
__attribute__((always_inline)) INLINE static void
chemistry_accumulate_timestep_vars(struct chemistry_part_data *restrict chi,
                                   struct chemistry_part_data *restrict chj,
                                   const int status, const float vmax,
                                   const float delxbar_i, const float delxbar_j,
                                   const int accum_ts_left,
                                   const int accum_ts_right) {

#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  if (accum_ts_left) chi->timestepvars.vmax = max(chi->timestepvars.vmax, vmax);
  if (accum_ts_right)
    chj->timestepvars.vmax = max(chj->timestepvars.vmax, vmax);
#endif

#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  /* An interface without area has no particle-size term. */
  if (status == CHEMISTRY_PAIR_FLUX_VMAX) return;

  /* Lanson & Vila (2008), denominator of equation (58) */
  if (accum_ts_left) chi->timestepvars.delxbar += delxbar_i;
  if (accum_ts_right) chj->timestepvars.delxbar += delxbar_j;
#endif
}

/**
 * @brief Apply the fluxes of a pair, in the apply pass.
 *
 * With the positivity limiter on, each flux is first scaled by the theta of
 * its donor, which is frozen since the chemistry flux ghost. Every visit of
 * a pair then applies the same scaled flux, so metal mass is conserved.
 *
 * @param pi Particle i, the left state.
 * @param pj Particle j, the right state.
 * @param totflux The fluxes of all elements across the interface, scaled in
 * place by theta.
 * @param flux_dt The time interval of the exchange.
 * @param update_left Whether pi receives the flux.
 * @param update_right Whether pj receives the flux.
 * @param chem_data The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_apply_pair_fluxes(
    struct part *restrict pi, struct part *restrict pj,
    double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4], const float flux_dt,
    const int update_left, const int update_right,
    const struct chemistry_global_data *chem_data) {

  const struct chemistry_part_data *chi = &pi->chemistry_data;
  const struct chemistry_part_data *chj = &pj->chemistry_data;

  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    if (chem_data->use_positivity_limiter) {
      if (totflux[m][0] > 0.0) {
        totflux[m][0] *= chi->fct_theta[m];
      } else if (totflux[m][0] < 0.0) {
        totflux[m][0] *= chj->fct_theta[m];
      }
    }

    /* A positive flux is subtracted from the left state and added to the
       right state. Unlike in SPH, an inactive neighbour is updated too, so
       that the exchange stays symmetric. */
    if (update_left)
      chemistry_part_update_fluxes_left(pi, m, totflux[m], flux_dt);
    if (update_right)
      chemistry_part_update_fluxes_right(pj, m, totflux[m], flux_dt);
  }
}

/**
 * @brief Is this particle inactive for the flux exchange of this step?
 *
 * This is the time-bin test of part_is_active_no_debug(), which the
 * interaction loops use to pick the branch, so it holds for foreign copies
 * and in the prep pass. active.h cannot be included here: it includes this
 * file through engine.h. A particle woken by the limiter in this step is
 * still inactive under its old bin.
 *
 * @param p The #part.
 * @param max_active_bin The highest active time-bin.
 */
__attribute__((always_inline)) INLINE static int
chemistry_part_is_flux_inactive(const struct part *p,
                                const timebin_t max_active_bin) {
  return p->time_bin > max_active_bin;
}

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
/**
 * @brief Add the integer interval of a pair exchange to the particles it
 * updates, and count it as a catch-up if it is longer than the current step
 * of every active member of the pair.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param t_current The current integer time.
 * @param update_left Whether pi receives the flux.
 * @param update_right Whether pj receives the flux.
 */
__attribute__((always_inline)) INLINE static void
chemistry_debug_count_pair_interval(struct part *restrict pi,
                                    struct part *restrict pj,
                                    const integertime_t t_current,
                                    const int update_left,
                                    const int update_right) {

  const integertime_t interval =
      t_current - chemistry_get_pair_ti_flux_start(pi, pj);
  const timebin_t max_active_bin = get_max_active_bin(t_current);
  integertime_t step_min = max_nr_timesteps;
  if (!chemistry_part_is_flux_inactive(pi, max_active_bin))
    step_min = min(step_min, get_integer_timestep(pi->time_bin));
  if (!chemistry_part_is_flux_inactive(pj, max_active_bin))
    step_min = min(step_min, get_integer_timestep(pj->time_bin));
  const int catch_up = interval > step_min;

  if (update_left) {
    pi->chemistry_data.flux_coverage_ticks += interval;
    pi->chemistry_data.flux_catch_up_count += catch_up;
  }
  if (update_right) {
    pj->chemistry_data.flux_coverage_ticks += interval;
    pj->chemistry_data.flux_catch_up_count += catch_up;
  }
}
#endif

/**
 * @brief Compute the metal diffusion flux of a pair once, and use it in the
 * current pass of the positivity limiter.
 *
 * The prep pass adds the donor outflow. The apply pass adds the timestep
 * variables and applies the fluxes. Both passes call the out-of-line
 * chemistry_compute_pair_fluxes() with the same inputs, so they see the same
 * fluxes.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param t_current The current time (in integer).
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 * @param chem_data The global properties of the chemistry scheme.
 * @param update_left Whether pi receives the flux.
 * @param update_right Whether pj receives the flux.
 * @param accum_ts_left Whether pi gets the timestep variables.
 * @param accum_ts_right Whether pj gets the timestep variables.
 * @param pass The pass of the positivity limiter.
 */
__attribute__((always_inline)) INLINE static void chemistry_exchange_pair_flux(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const double time_base,
    const integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology, const struct chemistry_global_data *chem_data,
    const int update_left, const int update_right, const int accum_ts_left,
    const int accum_ts_right, const enum chemistry_flux_pass pass) {

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  if (update_left) pi->chemistry_data.flux_exchange_count[pass]++;
  if (update_right) pj->chemistry_data.flux_exchange_count[pass]++;
  if (pass == chemistry_flux_pass_apply)
    chemistry_debug_count_pair_interval(pi, pj, t_current, update_left,
                                        update_right);
#endif

  /* The flux computation is out of line, so that the prep and the apply
     loops run the same compiled code. Two inlined copies in different loop
     bodies can round differently under fast-math. */
  double totflux[GEAR_CHEMISTRY_ELEMENT_COUNT][4];
  float flux_dt = 0.f, vmax = 0.f;
  float delxbar_i = 0.f, delxbar_j = 0.f;
#ifdef GIZMO_LANSON_VILA_PARTICLE_SIZE
  const int status = chemistry_compute_pair_fluxes(
      r2, dx, hi, hj, pi, pj, time_base, t_current, cosmo, with_cosmology,
      chem_data, totflux, &flux_dt, &vmax, &delxbar_i, &delxbar_j);
#else
  const int status = chemistry_compute_pair_fluxes(
      r2, dx, hi, hj, pi, pj, time_base, t_current, cosmo, with_cosmology,
      chem_data, totflux, &flux_dt, &vmax);
#endif

  if (status == CHEMISTRY_PAIR_FLUX_NONE) return;

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (pass == chemistry_flux_pass_prep) {
    if (status == CHEMISTRY_PAIR_FLUX_FLUXES)
      chemistry_accumulate_donor_outflow(chi, chj, totflux, flux_dt,
                                         update_left, update_right);
    return;
  }

  chemistry_accumulate_timestep_vars(chi, chj, status, vmax, delxbar_i,
                                     delxbar_j, accum_ts_left, accum_ts_right);

  if (status == CHEMISTRY_PAIR_FLUX_FLUXES)
    chemistry_apply_pair_fluxes(pi, pj, totflux, flux_dt, update_left,
                                update_right, chem_data);
}

/**
 * @brief Exchange the flux of a pair whose neighbour is a foreign copy.
 *
 * Both ranks compute the flux with the two particles in the same order,
 * sorted by position, so they get the same floating-point result. Only the
 * local particle pi is updated; the foreign copy pj is never written. A tie
 * in position skips the degenerate face.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi The local particle.
 * @param pj The foreign copy of its neighbour.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param t_current The current time (in integer).
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 * @param chem_data The global properties of the chemistry scheme.
 * @param pass The pass of the positivity limiter.
 */
__attribute__((always_inline)) INLINE static void
chemistry_exchange_across_ranks(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const double time_base,
    const integertime_t t_current, const struct cosmology *cosmo,
    const int with_cosmology, const struct chemistry_global_data *chem_data,
    const enum chemistry_flux_pass pass) {

  const int pi_less_pj = flux_exchange_precedes(pi->x, pj->x);
  const int pj_less_pi = flux_exchange_precedes(pj->x, pi->x);
  if (!pi_less_pj && !pj_less_pi) {
#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
    warning(
        "Coincident distinct particles (ids %lld/%lld) at pos [%g %g %g] "
        "in the cross-rank chemistry flux exchange: skipping degenerate "
        "face.",
        pi->id, pj->id, pi->x[0], pi->x[1], pi->x[2]);
#endif
    return;
  }
  if (pi_less_pj) {
    chemistry_exchange_pair_flux(r2, dx, hi, hj, pi, pj, time_base, t_current,
                                 cosmo, with_cosmology, chem_data,
                                 /*update_left=*/1, /*update_right=*/0,
                                 /*accum_ts_left=*/1, /*accum_ts_right=*/0,
                                 pass);
  } else {
    const float mdx[3] = {-dx[0], -dx[1], -dx[2]};
    chemistry_exchange_pair_flux(r2, mdx, hj, hi, pj, pi, time_base, t_current,
                                 cosmo, with_cosmology, chem_data,
                                 /*update_left=*/0, /*update_right=*/1,
                                 /*accum_ts_left=*/0, /*accum_ts_right=*/1,
                                 pass);
  }
}

/**
 * @brief Exchange the metal diffusion flux of a pair, in the prep pass (the
 * chemistry flux prep loop) or in the apply pass (the force loop).
 *
 * @p pi is the particle this call may update and @p pj its neighbour. The
 * cases, in this order:
 *   1. both_updatable_here: the gate of the call site covers the pair
 *      symmetrically. Both particles are updated.
 *   2. pj inactive: this is the only visit of the pair this step. Both
 *      particles are updated, and only pi gets the timestep variables. pi
 *      can be a foreign copy, whose update is then discarded.
 *   3. otherwise pj is an active foreign copy: only pi is updated, see
 *      chemistry_exchange_across_ranks(). Two local active particles always
 *      take case 1.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi The particle this call site may update.
 * @param pj Its neighbour.
 * @param both_updatable_here Whether the gate of this call site already
 * guarantees both particles active and local (the one symmetric visit of the
 * pair this step).
 * @param local_first Whether the cell @p pi belongs to is local to this
 * rank. Named by argument position, since several call sites pass the
 * particle of cell j first.
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
 * @param pass The pass of the positivity limiter.
 * @param max_active_bin The highest active time-bin.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_chemistry_flux_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj,
    const int both_updatable_here, const int local_first,
    const int local_second, const float a, const float H,
    const double time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data *chem_data,
    const enum chemistry_flux_pass pass, const timebin_t max_active_bin) {

#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
  if (both_updatable_here && !(local_first && local_second))
    error(
        "runner_iact_chemistry_flux_exchange called with both_updatable_here "
        "on a pair with a foreign particle (ids %lld/%lld).",
        pi->id, pj->id);
  if (max_active_bin != get_max_active_bin(t_current))
    error("max_active_bin=%d does not match t_current=%lld.", max_active_bin,
          t_current);
  if (chemistry_part_is_flux_inactive(pi, max_active_bin))
    error("The first particle of the flux exchange (id %lld) is inactive.",
          pi->id);
#endif

  if (both_updatable_here) {
    chemistry_exchange_pair_flux(r2, dx, hi, hj, pi, pj, time_base, t_current,
                                 cosmo, with_cosmology, chem_data,
                                 /*update_left=*/1, /*update_right=*/1,
                                 /*accum_ts_left=*/1, /*accum_ts_right=*/1,
                                 pass);
  } else if (chemistry_part_is_flux_inactive(pj, max_active_bin)) {
    chemistry_exchange_pair_flux(r2, dx, hi, hj, pi, pj, time_base, t_current,
                                 cosmo, with_cosmology, chem_data,
                                 /*update_left=*/1, /*update_right=*/1,
                                 /*accum_ts_left=*/1, /*accum_ts_right=*/0,
                                 pass);
  } else {
#ifdef SWIFT_CHEMISTRY_DEBUG_CHECKS
    if (local_first && local_second)
      error(
          "Unreachable: both local active particles (ids %lld/%lld) take the "
          "symmetric call.",
          pi->id, pj->id);
    if (!local_first)
      error(
          "runner_iact_chemistry_flux_exchange's cross-rank arm reached with "
          "local_first false: the call site's first argument must always be "
          "the locally-owned particle.");
#endif
    chemistry_exchange_across_ranks(r2, dx, hi, hj, pi, pj, time_base,
                                    t_current, cosmo, with_cosmology, chem_data,
                                    pass);
  }
}

/* The template's diffusion hooks are empty here: this scheme uses
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

#endif /* SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_IACT_H */
