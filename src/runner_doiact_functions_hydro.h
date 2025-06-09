/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_hydro.h"

/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR1_NAIVE(struct runner *r, const struct cell *restrict ci,
                   const struct cell *restrict cj, const int limit_min_h,
                   const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->dmin != cj->dmin) error("Cells of different size!");
#endif

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = PART_IS_ACTIVE(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = PART_IS_ACTIVE(pj, e);
      const char depth_j = pj->depth_h;
      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      const int doi = pi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = pj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

      /* Hit or miss? */
      if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
      if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions between a cell pair (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR2_NAIVE(struct runner *r, const struct cell *restrict ci,
                   const struct cell *restrict cj, const int limit_min_h,
                   const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = PART_IS_ACTIVE(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = PART_IS_ACTIVE(pj, e);
      const char depth_j = pj->depth_h;
      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      const int doi = pi_active && (depth_i >= min_depth) &&
                      (depth_i <= max_depth) && ((r2 < hig2) || (r2 < hjg2));
      const int doj = pj_active && (depth_j >= min_depth) &&
                      (depth_j <= max_depth) && ((r2 < hjg2) || (r2 < hig2));

      /* Hit or miss? */
      if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
      if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1_NAIVE(struct runner *r, const struct cell *c,
                   const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(c, e)) return;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  const int count = c->hydro.count;
  struct part *parts = c->hydro.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = PART_IS_ACTIVE(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = PART_IS_ACTIVE(pj, e);
      const char depth_j = pj->depth_h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = pi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = pj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        IACT(r2, dx, hi, hj, pi, pj, a, H);
        IACT_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                              t_current, cosmo, with_cosmology);
#endif
      } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions within a cell (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF2_NAIVE(struct runner *r, const struct cell *c,
                   const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(c, e)) return;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  const int count = c->hydro.count;
  struct part *parts = c->hydro.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = PART_IS_ACTIVE(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = PART_IS_ACTIVE(pj, e);
      const char depth_j = pj->depth_h;
      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = pi_active && (depth_i >= min_depth) &&
                      (depth_i <= max_depth) && ((r2 < hig2) || (r2 < hjg2));
      const int doj = pj_active && (depth_j >= min_depth) &&
                      (depth_j <= max_depth) && ((r2 < hig2) || (r2 < hjg2));

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        IACT(r2, dx, hi, hj, pi, pj, a, H);
        IACT_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                              t_current, cosmo, with_cosmology);
#endif
      } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET_NAIVE(struct runner *r, const struct cell *restrict ci,
                         struct part *restrict parts_i, const int *ind,
                         const int count, const struct cell *restrict cj,
                         const double shift[3]) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  /* Loop over the parts_i. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[ind[pid]];
    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, pj->h, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, pj->h, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, pj->h, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, pj->h, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(timer_dopair_subset_naive);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param flipped Flag to check whether the cells have been flipped or not.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET(struct runner *r, const struct cell *restrict ci,
                   struct part *restrict parts_i, const int *ind,
                   const int count, const struct cell *restrict cj,
                   const int sid, const int flipped, const double shift[3]) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  /* Pick-out the sorted lists. */
  const struct sort_entry *sort_j = cell_get_hydro_sorts(cj, sid);
  const float dxj = cj->hydro.dx_max_sort;

  /* Parts are on the left? */
  if (!flipped) {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the parts in cj. */
    } /* loop over the parts in ci. */
  }

  /* Parts are on the right. */
  else {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the parts in cj. */
    } /* loop over the parts in ci. */
  }

  TIMER_TOC(timer_dopair_subset);
}

/**
 * @brief Determine which version of DOPAIR_SUBSET needs to be called depending
 * on the
 * orientation of the cells or whether DOPAIR_SUBSET needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR_SUBSET_BRANCH(struct runner *r, const struct cell *restrict ci,
                          struct part *restrict parts_i, const int *ind,
                          const int count, struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->hydro.count == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Get the sorting index. */
  int sid = 0;
  for (int k = 0; k < 3; k++)
    sid = 3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0)   ? 0
                     : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2
                                                                : 1);

  /* Switch the cells around? */
  const int flipped = runner_flip[sid];
  sid = sortlistID[sid];

  /* Let's first lock the cell */
  lock_lock(&cj->hydro.extra_sort_lock);

  /* Is it sorted, if not we use the naive interactions. */
  const int is_sorted =
      (cj->hydro.sorted & (1 << sid)) &&
      (cj->hydro.dx_max_sort_old <= space_maxreldx * cj->dmin);

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  const int force_naive = 1;
#else
  const int force_naive = 0;
#endif

  if (force_naive || !is_sorted) {
    DOPAIR_SUBSET_NAIVE(r, ci, parts_i, ind, count, cj, shift);
  } else {
#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH)
    if (sort_is_face(sid))
      runner_dopair_subset_density_vec(r, ci, parts_i, ind, count, cj, sid,
                                       flipped, shift);
    else
      DOPAIR_SUBSET(r, ci, parts_i, ind, count, cj, sid, flipped, shift);
#else
    DOPAIR_SUBSET(r, ci, parts_i, ind, count, cj, sid, flipped, shift);
#endif
  }

  /* Now we can unlock */
  if (lock_unlock(&cj->hydro.extra_sort_lock) != 0)
    error("Impossible to unlock cell!");
}

/**
 * @brief Compute the interactions between a cell, but only for the
 *      given indices in c.
 *
 * @param r The #runner.
 * @param i The #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET(struct runner *r, const struct cell *c,
                   struct part *restrict parts, const int *ind,
                   const int count) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  const int count_cell = c->hydro.count;
  struct part *restrict parts_j = c->hydro.parts;
  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *pi = &parts[ind[pid]];
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_cell; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip oneself */
      if (pi == pj) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
        runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
        runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                     t_current, cosmo, with_cosmology);
#endif
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(timer_doself_subset);
}

/**
 * @brief Determine which version of DOSELF_SUBSET needs to be called depending
 * on the optimisation level.

 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET_BRANCH(struct runner *r, const struct cell *ci,
                          struct part *restrict parts, const int *ind,
                          const int count) {

#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH)
  runner_doself_subset_density_vec(r, ci, parts, ind, count);
#else
  DOSELF_SUBSET(r, ci, parts, ind, count);
#endif
}

/**
 * @brief Compute the interactions between a cell pair (non-symmetric).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1(struct runner *r, const struct cell *restrict ci,
             const struct cell *restrict cj, const int limit_min_h,
             const int limit_max_h, const int sid, const double shift[3]) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Get the cut_off shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
  /* Some constants used to checks that the parts are in the right frame */
  /* TODO MLADEN: coordinate 2. -> 2.02 with Matthieu */
  const float shift_threshold_x =
      2.02 * ci->width[0] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2.02 * ci->width[1] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2.02 * ci->width[2] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
#endif
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;

  /* Get some other useful values. */
  const double hi_max =
      min(h_max, ci->hydro.h_max_active) * kernel_gamma - rshift;
  const double hj_max = min(h_max, cj->hydro.h_max_active) * kernel_gamma;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  if (CELL_IS_ACTIVE(ci, e)) {

    /* Loop over the *active* parts in ci that are within range (on the axis)
       of any particle in cj. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      const char depth_i = pi->depth_h;
      const float hi = pi->h;

      /* Skip inactive particles */
      if (!PART_IS_ACTIVE(pi, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hi > ci->hydro.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_i < min_depth) continue;
      if (depth_i > max_depth) continue;

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the parts in cj. */
    } /* loop over the parts in ci. */
  } /* Cell ci is active */

  if (CELL_IS_ACTIVE(cj, e)) {

    /* Loop over the *active* parts in cj that are within range (on the axis)
       of any particle in ci. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *pj = &parts_j[sort_j[pjd].i];
      const char depth_j = pj->depth_h;
      const float hj = pj->h;

      /* Skip inactive particles */
      if (!PART_IS_ACTIVE(pj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hj > cj->hydro.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_j < min_depth) continue;
      if (depth_j > max_depth) continue;

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = pj->x[0] - cj->loc[0];
      const float pjy = pj->x[1] - cj->loc[1];
      const float pjz = pj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
          IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the parts in ci. */
    } /* loop over the parts in cj. */
  } /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

  /* Check that cells are drifted. */
  if (!CELL_ARE_PART_DRIFTED(ci, e) || !CELL_ARE_PART_DRIFTED(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID.
   * Note: this may swap the ci and cj pointers!! */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells (ci).");

  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells (cj).");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR1_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (!sort_is_corner(sid))
    runner_dopair1_density_vec(r, ci, cj, sid, shift);
  else
    DOPAIR1(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#else
  DOPAIR1(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#endif
}

/**
 * @brief Compute the interactions between a cell pair (symmetric)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 * @param sid The direction of the pair
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR2(struct runner *r, const struct cell *restrict ci,
             const struct cell *restrict cj, const int limit_min_h,
             const int limit_max_h, const int sid, const double shift[3]) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  /* Get the cut_off shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
  /* Some constants used to checks that the parts are in the right frame */
  /* TODO MLADEN: coordinate 2. -> 2.02 with Matthieu */
  const float shift_threshold_x =
      2.02 * ci->width[0] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2.02 * ci->width[1] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2.02 * ci->width[2] +
      2.02 * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;
#endif

  /* Get some other useful values. */
  const int local_i = ci->nodeID == e->nodeID;
  const int local_j = cj->nodeID == e->nodeID;
  const double hi_max = ci->hydro.h_max;
  const double hj_max = cj->hydro.h_max;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  /* Maximal displacement since last rebuild */
  const double dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Position on the axis of the particles closest to the interface */
  const double di_max = sort_i[count_i - 1].d;
  const double dj_min = sort_j[0].d;

  /* Shifts to apply to the particles to be in a good frame */
  const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                             cj->loc[2] + shift[2]};
  const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  int count_active_i = 0, count_active_j = 0;
  struct sort_entry *restrict sort_active_i = NULL;
  struct sort_entry *restrict sort_active_j = NULL;

  // MATTHIEU: temporary disable this optimization
  if (0 /*&& cell_is_all_active_hydro(ci, e)*/) {
    /* If everybody is active don't bother copying */
    sort_active_i = sort_i;
    count_active_i = count_i;
  } else if (CELL_IS_ACTIVE(ci, e)) {
    if (posix_memalign((void **)&sort_active_i, SWIFT_CACHE_ALIGNMENT,
                       sizeof(struct sort_entry) * count_i) != 0)
      error("Failed to allocate active sortlists.");

    /* Collect the active particles in ci */
    for (int k = 0; k < count_i; k++) {
      const struct part *p = &parts_i[sort_i[k].i];
      const char depth = p->depth_h;

      if (PART_IS_ACTIVE(p, e) && (depth >= min_depth) &&
          (depth <= max_depth)) {
        sort_active_i[count_active_i] = sort_i[k];
        count_active_i++;
      }
    }
  }

  // MATTHIEU: temporary disable this optimization
  if (0 /*&& cell_is_all_active_hydro(cj, e)*/) {
    /* If everybody is active don't bother copying */
    sort_active_j = sort_j;
    count_active_j = count_j;
  } else if (CELL_IS_ACTIVE(cj, e)) {
    if (posix_memalign((void **)&sort_active_j, SWIFT_CACHE_ALIGNMENT,
                       sizeof(struct sort_entry) * count_j) != 0)
      error("Failed to allocate active sortlists.");

    /* Collect the active particles in cj */
    for (int k = 0; k < count_j; k++) {
      const struct part *p = &parts_j[sort_j[k].i];
      const char depth = p->depth_h;

      if (PART_IS_ACTIVE(p, e) && (depth >= min_depth) &&
          (depth <= max_depth)) {
        sort_active_j[count_active_j] = sort_j[k];
        count_active_j++;
      }
    }
  }

  /* Loop over *all* the parts in ci starting from the centre until
     we are out of range of anything in cj (using the maximal hi). */
  for (int pid = count_i - 1;
       pid >= 0 &&
       sort_i[pid].d + hi_max * kernel_gamma + dx_max - rshift > dj_min;
       pid--) {

    /* Get a hold of the ith part in ci. */
    struct part *pi = &parts_i[sort_i[pid].i];
    const char depth_i = pi->depth_h;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const float hi = pi->h;

    /* Is there anything we need to interact with (for this specific hi) ? */
    const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
    if (di < dj_min) continue;

    /* Get some additional information about pi */
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix = pi->x[0] - shift_i[0];
    const float piy = pi->x[1] - shift_i[1];
    const float piz = pi->x[2] - shift_i[2];

    const int update_i = PART_IS_ACTIVE(pi, e) && local_i &&
                         (depth_i >= min_depth) && (depth_i <= max_depth);

    /* Do we need to only check active parts in cj
       (i.e. pi does not need updating) ? */
    if (!update_i) {

      /* Loop over the *active* parts in cj within range of pi */
      for (int pjd = 0; pjd < count_active_j && sort_active_j[pjd].d < di;
           pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_active_j[pjd].i];

        /* Skip inhibited particles.
         * Note we are looping over active particles but in the case where
         * the cell thinks all the particles are active (because of the
         * ti_end_max), particles may have nevertheless been inhibted by BH
         * swallowing in the mean time. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Get the position of pj in the right frame */
        const float pjx = pj->x[0] - shift_j[0];
        const float pjy = pj->x[1] - shift_j[1];
        const float pjz = pj->x[2] - shift_j[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");

        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
#endif

        /* Hit or miss?
           (note that we will do the other condition in the reverse loop) */
        if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
          IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the active parts in cj. */
    }

    else { /* pi is active, we may need to update pi and pj */

      /* Loop over *all* the parts in cj in range of pi. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];
        const char depth_j = pj->depth_h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Get the position of pj in the right frame */
        const float pjx = pj->x[0] - shift_j[0];
        const float pjy = pj->x[1] - shift_j[1];
        const float pjz = pj->x[2] - shift_j[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");

        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

#endif

        /* Hit or miss?
           (note that we will do the other condition in the reverse loop) */
        if (r2 < hig2) {

          const int doj = PART_IS_ACTIVE(pj, e) && local_j &&
                          (depth_j >= min_depth) && (depth_j <= max_depth);

          /* Does pj need to be updated too? */
          if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
            if (hi < h_min || hi >= h_max)
              error("Inappropriate h for this level!");
            if (hj < h_min || hj >= h_max)
              error("Inappropriate h for this level!");
#endif

            IACT(r2, dx, hi, hj, pi, pj, a, H);
            IACT_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_star_formation(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
            runner_iact_timebin(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                  t_current, cosmo, with_cosmology);
#endif
          } else {

#ifdef SWIFT_DEBUG_CHECKS
            if (hi < h_min || hi >= h_max)
              error("Inappropriate h for this level!");
#endif
            IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
            IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
            runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
            runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H,
                                         time_base, t_current, cosmo,
                                         with_cosmology);
#endif
          }
        }
      } /* loop over the parts in cj. */
    } /* Is pi active? */
  } /* Loop over all ci */

  /* Loop over *all* the parts in cj starting from the centre until
     we are out of range of anything in ci (using the maximal hj). */
  for (int pjd = 0;
       pjd < count_j &&
       sort_j[pjd].d - hj_max * kernel_gamma - dx_max < di_max - rshift;
       pjd++) {

    /* Get a hold of the jth part in cj. */
    struct part *pj = &parts_j[sort_j[pjd].i];
    const char depth_j = pj->depth_h;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pj, e)) continue;

    const float hj = pj->h;

    /* Is there anything we need to interact with (for this specific hj) ? */
    const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max;
    if (dj > di_max - rshift) continue;

    /* Get some additional information about pj */
    const float hjg2 = hj * hj * kernel_gamma2;
    const float pjx = pj->x[0] - shift_j[0];
    const float pjy = pj->x[1] - shift_j[1];
    const float pjz = pj->x[2] - shift_j[2];

    const int update_j = PART_IS_ACTIVE(pj, e) && local_j &&
                         (depth_j >= min_depth) && (depth_j <= max_depth);

    /* Do we need to only check active parts in ci
       (i.e. pj does not need updating) ? */
    if (!update_j) {

      /* Loop over the *active* parts in ci. */
      for (int pid = count_active_i - 1;
           pid >= 0 && sort_active_i[pid].d - rshift > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_active_i[pid].i];

        /* Skip inhibited particles.
         * Note we are looping over active particles but in the case where
         * the cell thinks all the particles are active (because of the
         * ti_end_max), particles may have nevertheless been inhibted by BH
         * swallowing in the mean time. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;

        /* Get the position of pi in the right frame */
        const float pix = pi->x[0] - shift_i[0];
        const float piy = pi->x[1] - shift_i[1];
        const float piz = pi->x[2] - shift_i[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
#endif

        /* Hit or miss?
           (note that we must avoid the r2 < hig2 cases we already processed) */
        if (r2 < hjg2 && r2 >= hig2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over the active parts in ci. */
    }

    else { /* pj is active, we may need to update pj and pi */

      /* Loop over *all* the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d - rshift > dj;
           pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];
        const char depth_i = pi->depth_h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;

        /* Get the position of pi in the right frame */
        const float pix = pi->x[0] - shift_i[0];
        const float piy = pi->x[1] - shift_i[1];
        const float piz = pi->x[2] - shift_i[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
#endif

        /* Hit or miss?
           (note that we must avoid the r2 < hig2 cases we already processed) */
        if (r2 < hjg2 && r2 >= hig2) {

          const int doi = PART_IS_ACTIVE(pi, e) && local_i &&
                          (depth_i >= min_depth) && (depth_i <= max_depth);

          /* Does pi need to be updated too? */
          if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
            if (hi < h_min || hi >= h_max)
              error("Inappropriate h for this level!");
            if (hj < h_min || hj >= h_max)
              error("Inappropriate h for this level!");
#endif

            IACT(r2, dx, hj, hi, pj, pi, a, H);
            IACT_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_star_formation(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
            runner_iact_timebin(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                  t_current, cosmo, with_cosmology);
#endif
          } else {

#ifdef SWIFT_DEBUG_CHECKS
            if (hj < h_min || hj >= h_max)
              error("Inappropriate h for this level!");
#endif

            IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
            IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
            runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
            runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H,
                                         time_base, t_current, cosmo,
                                         with_cosmology);
#endif
          }
        }
      } /* loop over the parts in ci. */
    } /* Is pj active? */
  } /* Loop over all cj */

  /* Clean-up if necessary */  // MATTHIEU: temporary disable this optimization
  if (CELL_IS_ACTIVE(ci, e))   // && !cell_is_all_active_hydro(ci, e))
    free(sort_active_i);
  if (CELL_IS_ACTIVE(cj, e))  // && !cell_is_all_active_hydro(cj, e))
    free(sort_active_j);

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Determine which version of DOPAIR2 needs to be called depending on the
 * orientation of the cells or whether DOPAIR2 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

  /* Check that cells are drifted. */
  if (!CELL_ARE_PART_DRIFTED(ci, e) || !CELL_ARE_PART_DRIFTED(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells (ci).");

  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells (cj).");

#ifdef SWIFT_USE_NAIVE_INTERACTIONS
  DOPAIR2_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  if (!sort_is_corner(sid))
    runner_dopair2_force_vec(r, ci, cj, sid, shift);
  else
    DOPAIR2(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#else
  DOPAIR2(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#endif
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1(struct runner *r, const struct cell *c, const int limit_min_h,
             const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  struct part *parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Set up a list of the particles for which we want to compute interactions */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++) {
    const struct part *p = &parts[k];
    const char depth = p->depth_h;
    if (PART_IS_ACTIVE(p, e) && (depth >= min_depth) && (depth <= max_depth)) {
      indt[countdt] = k;
      countdt += 1;
    }
  }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  /* Loop over *all* the particles (i.e. the ones to update and not to update).
   *
   * Note the additional condition to make the loop abort if all the active
   * particles have been processed. */
  for (int pid = 0; pid < count && firstdt < countdt; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];
    const char depth_i = pi->depth_h;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    /* Get the particle position and (square of) search radius. */
    const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle active and in the range of h we care about? */
    const int update_i = PART_IS_ACTIVE(pi, e) && (depth_i >= min_depth) &&
                         (depth_i <= max_depth);

    /* If false then it can only act as a neighbour of others */
    if (!update_i) {

      /* Loop over the particles we want to update. */
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. (by construction pi != pj) */
        struct part *restrict pj = &parts[indt[pjd]];

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        const float dx[3] = {(float)(pjx[0] - pix[0]), (float)(pjx[1] - pix[1]),
                             (float)(pjx[2] - pix[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
          IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over all the particles we want to update. */
    }

    /* Otherwise, interact with all candidates. */
    else {

      /* We caught a live one!
       * Move the start of the list of active ones by one slot as it will have
       * been fully processed after the following loop so no need to consider it
       * in the previous loop any more. */
      firstdt += 1;

      /* Loop over *all* the particles (i.e. the ones to update and not to
       * update) but starting from where we are in the overall list. */
      for (int pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle (by construction pi != pj). */
        struct part *restrict pj = &parts[pjd];
        const char depth_j = pj->depth_h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                       (float)(pix[2] - pjx[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Decide which of the two particles to update */

        /* We know pi is active and in the right range of h
         * -> Only check the distance to pj */
        const int doi = (r2 < hig2);

        /* We know nothing about pj
         * -> Check whether it is active
         * -> Check whether it is in the right range of h
         * -> Check the distance to pi */
        const int doj = (PART_IS_ACTIVE(pj, e)) && (depth_j >= min_depth) &&
                        (depth_j <= max_depth) && (r2 < hjg2);

        /* Hit or miss? */
        if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif
          /* Update both pi and pj */

          IACT(r2, dx, hi, hj, pi, pj, a, H);
          IACT_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                t_current, cosmo, with_cosmology);
#endif
        } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif
          /* Update only pi */

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update only pj */

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
          IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        } /* Hit or miss */
      } /* loop over all other particles. */
    } /* pi is active */
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Determine which version of DOSELF1 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1_BRANCH(struct runner *r, const struct cell *c,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    if (!limit_max_h && c->hydro.h_max_active * kernel_gamma > c->dmin)
      error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");
#endif

  /* Check that cells are drifted. */
  if (!CELL_ARE_PART_DRIFTED(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_NAIVE(r, c, limit_min_h, limit_max_h);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  runner_doself1_density_vec(r, c);
#else
  DOSELF1(r, c, limit_min_h, limit_max_h);
#endif
}

/**
 * @brief Compute the cell self-interaction (symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF2(struct runner *r, const struct cell *c, const int limit_min_h,
             const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  const double time_base = e->time_base;
  const integertime_t t_current = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
#endif

  TIMER_TIC;

  struct part *parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Set up a list of the particles for which we want to compute interactions */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++) {
    const struct part *p = &parts[k];
    const char depth = p->depth_h;
    if (PART_IS_ACTIVE(p, e) && (depth >= min_depth) && (depth <= max_depth)) {
      indt[countdt] = k;
      countdt += 1;
    }
  }

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;
  GET_MU0();

  /* Loop over *all* the particles (the ones to update and others!) in the cell.
   *
   * Note the additional condition to make the loop abort if all the active
   * particles have been processed. */
  for (int pid = 0; pid < count && firstdt < countdt; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];
    const char depth_i = pi->depth_h;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    /* Get the particle position and (square of) search radius. */
    const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle active and in the range of h we care about? */
    const int update_i = PART_IS_ACTIVE(pi, e) && (depth_i >= min_depth) &&
                         (depth_i <= max_depth);

    /* If false then it can only act as a neighbour of others */
    if (!update_i) {

      /* Loop over the active particles we want to update. */
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. (by construction pi != pj) */
        struct part *restrict pj = &parts[indt[pjd]];

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        const float dx[3] = {(float)(pjx[0] - pix[0]), (float)(pjx[1] - pix[1]),
                             (float)(pjx[2] - pix[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
          IACT_NONSYM_MHD(r2, dx, hj, hi, pj, pi, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_sink(r2, dx, hj, hi, pj, pi, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hj, hi, pj, pi, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        }
      } /* loop over all other particles. */
    }

    /* Otherwise, interact with all candidates. */
    else {

      /* We caught a live one!
       * Move the start of the list of active ones by one slot as it will have
       * been fully processed after the following loop so no need to consider it
       * in the previous loop any more. */
      firstdt += 1;

      /* Loop over *all* the particles (i.e. the ones to update and not to
       * update) but starting from where we are in the overall list. */
      for (int pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];
        const char depth_j = pj->depth_h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                       (float)(pix[2] - pjx[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Decide which of the two particles to update */

        /* We know pi is active and in the right range of h
         * -> Only check the distance to pj */
        const int doi = (r2 < hig2) || (r2 < hjg2);

        /* We know nothing about pj
         * -> Check whether it is active
         * -> Check whether it is in the right range of h
         * -> Check the distance to pi */
        const int doj = (PART_IS_ACTIVE(pj, e)) && (depth_j >= min_depth) &&
                        (depth_j <= max_depth) && ((r2 < hjg2) || (r2 < hig2));

        /* Hit or miss? */
        if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update both pi and pj */

          IACT(r2, dx, hi, hj, pi, pj, a, H);
          IACT_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                t_current, cosmo, with_cosmology);
#endif
        } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update only pi */

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
          IACT_NONSYM_MHD(r2, dx, hi, hj, pi, pj, mu_0, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_star_formation(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_sink(r2, dx, hi, hj, pi, pj, a, H);
#endif
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
          runner_iact_nonsym_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_rt_timebin(r2, dx, hi, hj, pi, pj, a, H);
          runner_iact_nonsym_diffusion(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                       t_current, cosmo, with_cosmology);
#endif
        } else if (doj) {

          /* Update only doj
           *
           * Note: This is impossible since if doj==True so does doi */

#ifdef SWIFT_DEBUG_CHECKS
          error("Impossible problem in the logic!!!");
#endif
        } /* Hit or miss */
      } /* loop over all other particles. */
    } /* pi is active */
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Determine which version of DOSELF2 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF2_BRANCH(struct runner *r, const struct cell *c,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    if (!limit_max_h && c->hydro.h_max_active * kernel_gamma > c->dmin)
      error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");

#endif

  /* Check that cells are drifted. */
  if (!CELL_ARE_PART_DRIFTED(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF2_NAIVE(r, c, limit_min_h, limit_max_h);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  runner_doself2_force_vec(r, c);
#else
  DOSELF2(r, c, limit_min_h, limit_max_h);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param recurse_below_h_max Are we currently recursing at a level where we
 * violated the h < cell size condition.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int recurse_below_h_max, const int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Should we even bother? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->hydro.count < space_recurse_size_pair_hydro ||
      !cj->split || cj->hydro.count < space_recurse_size_pair_hydro) {

    /* Do any of the cells need to be sorted first?
     * Since h_max might have changed, we may not have sorted at this level */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH(r, ci, cj, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Both ci and cj are split */

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max && (!cell_can_recurse_in_subpair_hydro_task(ci) ||
                                 !cell_can_recurse_in_subpair_hydro_task(cj))) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Do any of the cells need to be sorted first?
       * Since h_max might have changed, we may not have sorted at this level */
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }
      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }

      /* message("Multi-level PAIR! ci->count=%d cj->count=%d", ci->hydro.count,
       */
      /* 	      cj->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR1_BRANCH(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1(r, ci->progeny[pid], cj->progeny[pjd], recurse_below_h_max,
                    /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param recurse_below_h_max Are we currently recursing at a level where we
 * violated the h < cell size condition.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1(struct runner *r, struct cell *c, int recurse_below_h_max,
                 const int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  if (c->hydro.count == 0 || !CELL_IS_ACTIVE(c, r->e)) return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->hydro.count < space_recurse_size_self_hydro) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH(r, c, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself_hydro_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1(r, c->progeny[k], recurse_below_h_max, /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1(r, c->progeny[k], c->progeny[j], recurse_below_h_max,
                        /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}

/**
 * @brief Compute grouped sub-cell interactions for pairs (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj,
                 int recurse_below_h_max, const int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->hydro.count < space_recurse_size_pair_hydro ||
      !cj->split || cj->hydro.count < space_recurse_size_pair_hydro) {

    /* Do any of the cells need to be sorted first?
     * Since h_max might have changed, we may not have sorted at this level */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR2_BRANCH(r, ci, cj, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max &&
        (!cell_can_recurse_in_subpair2_hydro_task(ci) ||
         !cell_can_recurse_in_subpair2_hydro_task(cj))) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Do any of the cells need to be sorted first?
       * Since h_max might have changed, we may not have sorted at this level */
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }
      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }

      /* message("Multi-level PAIR! ci->count=%d cj->count=%d", ci->hydro.count,
       */
      /* 	      cj->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR2_BRANCH(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR2(r, ci->progeny[pid], cj->progeny[pjd], recurse_below_h_max,
                    /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF2(struct runner *r, struct cell *c, int recurse_below_h_max,
                 const int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  if (c->hydro.count == 0 || !CELL_IS_ACTIVE(c, r->e)) return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->hydro.count < space_recurse_size_self_hydro) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF2_BRANCH(r, c, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself2_hydro_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF2_BRANCH(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF2(r, c->progeny[k], recurse_below_h_max, /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR2(r, c->progeny[k], c->progeny[j], recurse_below_h_max,
                        /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}

/**
 * @brief Find which sub-cell of a cell contain the subset of particles given
 * by the list of indices.
 *
 * Will throw an error if the sub-cell can't be found.
 *
 * @param c The #cell
 * @param parts An array of #part.
 * @param ind Index of the #part's in the particle array to find in the subs.
 */
struct cell *FIND_SUB(const struct cell *const c,
                      const struct part *const parts, const int *const ind) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!c->split) error("Can't search for subs in a non-split cell");
#endif

  /* Find out in which sub-cell of ci the parts are.
   *
   * Note: We only need to check the first particle in the list */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL) {
      if (&parts[ind[0]] >= &c->progeny[k]->hydro.parts[0] &&
          &parts[ind[0]] <
              &c->progeny[k]->hydro.parts[c->progeny[k]->hydro.count]) {
        return c->progeny[k];
      }
    }
  }
  error("Invalid sub!");
  return NULL;
}

void DOSUB_PAIR_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                       const int *ind, const int count, struct cell *cj,
                       const int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_pair_hydro_task(ci) &&
      cell_can_recurse_in_pair_hydro_task(cj)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB(ci, parts, ind);

    /* Get the type of pair and flip ci/cj if needed. */
    double shift[3];
    const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
        DOSUB_PAIR_SUBSET(r, ci->progeny[pid], parts, ind, count,
                          cj->progeny[pjd],
                          /*gettimer=*/0);
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
        DOSUB_PAIR_SUBSET(r, cj->progeny[pjd], parts, ind, count,
                          ci->progeny[pid],
                          /*gettimer=*/0);
    }

  }
  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_hydro(ci, e)) {

    /* Do any of the cells need to be drifted first? */
    if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");

    DOPAIR_SUBSET_BRANCH(r, ci, parts, ind, count, cj);
  }

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}

void DOSUB_SELF_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                       const int *ind, const int count, const int gettimer) {

  const struct engine *e = r->e;

  /* Should we even bother? */
  if (ci->hydro.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse? */
  if (ci->split && cell_can_recurse_in_self_hydro_task(ci)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB(ci, parts, ind);

    /* Loop over all progeny. */
    DOSUB_SELF_SUBSET(r, sub, parts, ind, count, /*gettimer=*/0);
    for (int j = 0; j < 8; j++)
      if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
        DOSUB_PAIR_SUBSET(r, sub, parts, ind, count, ci->progeny[j],
                          /*gettimer=*/0);
  }

  /* Otherwise, compute self-interaction. */
  else
    DOSELF_SUBSET_BRANCH(r, ci, parts, ind, count);
}
