/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Gas-gas neighbour loop for sink particle formation.
 *
 * Implements pair and self interaction loops that accumulate quantities on
 * active gas particles using a fixed aperture radius @p r_cut, rather than
 * the per-particle smoothing-length radius used by the standard hydro density
 * loop.  The aperture radius is passed explicitly as an argument to every
 * function, making the loop generic and independent of any particular subgrid
 * sink model.
 *
 * Before including this file, define FUNCTION (e.g. "prep_sink_formation") and
 * FUNCTION_TASK_LOOP (e.g. TASK_LOOP_PREP_SINK_FORMATION).  The macro
 * expansion creates the following set of functions:
 *   runner_dopair1_naive_hydro_aperture_FUNCTION   (brute-force reference)
 *   runner_doself1_naive_hydro_aperture_FUNCTION   (brute-force reference)
 *   runner_dopair1_hydro_aperture_FUNCTION         (sorted, optimised)
 *   runner_dopair1_branch_hydro_aperture_FUNCTION  (dispatch)
 *   runner_doself1_hydro_aperture_FUNCTION         (optimised, active-list)
 *   runner_doself1_branch_hydro_aperture_FUNCTION  (dispatch)
 *   runner_dosub_pair1_hydro_aperture_FUNCTION     (recursive pair)
 *   runner_dosub_self1_hydro_aperture_FUNCTION     (recursive self)
 *
 * Interaction calls are dispatched via the IACT_NONSYM_HYDRO_APERTURE and
 * IACT_HYDRO_APERTURE macros defined in runner_doiact_hydro_aperture.h, which
 * expand to runner_iact_nonsym_hydro_aperture_FUNCTION and
 * runner_iact_hydro_aperture_FUNCTION respectively.
 *
 * Only non-symmetric (_1_) variants are provided.  The formation loop uses a
 * pure gather pattern — each active particle independently accumulates from
 * its neighbours — so symmetric (_2_) variants and subset-reiteration
 * functions (used by the standard hydro ghost loop) are not needed.
 *
 * Recursion termination:
 *   DOSUB_PAIR1/SELF1 recurse while r_cut < 0.5 * cell->dmin.  Once the
 *   aperture equals or exceeds half the cell's minimum dimension, further
 *   subdivision offers no pruning benefit (all sub-cell pairs would interact
 *   anyway), so the leaf pair/self function is called at that level instead
 *   of descending to the true leaf cells. */

#include "runner_doiact_hydro_aperture.h"

/* ============================================================
 * DOPAIR1_NAIVE_HYDRO_APERTURE
 * Brute-force non-symmetric pair interaction.
 * ============================================================ */

/**
 * @brief Compute the interactions between a cell pair (non-symmetric,
 *        brute force).
 *
 * Loops over all combinations of particles in @p ci and @p cj and calls the
 * non-symmetric interaction function for any pair within the fixed aperture
 * radius @p r_cut.  Used as a reference implementation and when
 * #SWIFT_USE_NAIVE_INTERACTIONS is defined.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR1_NAIVE_HYDRO_APERTURE(struct runner *r, const struct cell *restrict ci,
                               const struct cell *restrict cj,
                               const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->dmin != cj->dmin) error("Cells of different size!");
  if (r_cut > ci->dmin || r_cut > cj->dmin)
    error("Cell sizes (%e, %e) smaller than r_cut (%e)", ci->dmin, cj->dmin,
          r_cut);
#endif

  /* Get the relative distance between the cell pair, wrapping. */
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
    const float hi = pi->h;
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
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time. */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (pi_active && r2 < r_cut2) {
        IACT_NONSYM_HYDRO_APERTURE(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
      if (pj_active && r2 < r_cut2) {
        float mdx[3] = {-dx[0], -dx[1], -dx[2]};
        IACT_NONSYM_HYDRO_APERTURE(r2, mdx, hj, hi, pj, pi, a, H,
                                with_self_gravity, cosmo, e->sink_properties);
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR_HYDRO_APERTURE);
}

/* ============================================================
 * DOSELF1_NAIVE_HYDRO_APERTURE
 * Brute-force non-symmetric self interaction.
 * ============================================================ */

/**
 * @brief Compute cell self-interactions (non-symmetric, brute force).
 *
 * Loops over all distinct pairs within @p c and calls the non-symmetric
 * interaction function for any pair within the fixed aperture radius @p r_cut.
 * Used as a reference and when #SWIFT_USE_NAIVE_INTERACTIONS is defined.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF1_NAIVE_HYDRO_APERTURE(struct runner *r, const struct cell *c,
                               const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  struct part *restrict parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

#ifdef SWIFT_DEBUG_CHECKS
  if (r_cut > c->dmin)
    error("Cell size (%e) smaller than r_cut (%e)", c->dmin, r_cut);
#endif

  /* Loop over the parts in c. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in c. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = PART_IS_ACTIVE(pi, e);
    const float hi = pi->h;

    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = PART_IS_ACTIVE(pj, e);
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pi->x[0] - pj->x[0]), (float)(pi->x[1] - pj->x[1]),
                     (float)(pi->x[2] - pj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time. */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < r_cut2) {
        if (pi_active) {
          IACT_NONSYM_HYDRO_APERTURE(r2, dx, hi, hj, pi, pj, a, H,
                                  with_self_gravity, cosmo, e->sink_properties);
        }
        if (pj_active) {
          float mdx[3] = {-dx[0], -dx[1], -dx[2]};
          IACT_NONSYM_HYDRO_APERTURE(r2, mdx, hj, hi, pj, pi, a, H,
                                  with_self_gravity, cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOSELF_HYDRO_APERTURE);
}

/* ============================================================
 * DOPAIR1_HYDRO_APERTURE
 * Sorted non-symmetric pair interaction with fixed aperture r_cut.
 * ============================================================ */

/**
 * @brief Compute non-symmetric pair interactions using sorted particle lists.
 *
 * Uses two independent sweeps over the sorted lists of @p ci and @p cj:
 *
 *   - Active-ci pass: for each active pi in @p ci (descending along the pair
 *     axis), scan @p cj particles that fall within the fixed scan bound
 *     @c sort_i[pid].d + r_cut + dx_max - rshift.  The hit condition is
 *     @c r2 < r_cut^2.
 *
 *   - Active-cj pass: for each active pj in @p cj (ascending), scan @p ci
 *     particles within @c sort_j[pjd].d - r_cut - dx_max + rshift.
 *
 * The scan bounds replace @c hi*kernel_gamma with @p r_cut, giving a uniform
 * fixed-aperture search independent of smoothing length.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param sid The direction index of the pair.
 * @param shift The shift vector to apply to @p ci particles.
 */
void DOPAIR1_HYDRO_APERTURE(struct runner *r, const struct cell *restrict ci,
                         const struct cell *restrict cj, const float r_cut,
                         const int sid, const double shift[3]) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  /* Get the shift along the pair axis. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

  /* Fixed-aperture scan bounds (replace hi/hj * kernel_gamma in the hydro
   * loop).  hi_max accounts for the rshift projection so that the outer guard
   * in the ci loop correctly prunes particles whose sphere cannot reach cj. */
  const double hi_max = (double)r_cut - rshift;
  const double hj_max = (double)r_cut;

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Sorted axis extent of each cell. */
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;

  /* Maximum displacement of particles from their sorted positions. */
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

#ifdef SWIFT_DEBUG_CHECKS
  if (r_cut > ci->dmin || r_cut > cj->dmin)
    error("Cell sizes (%e, %e) smaller than r_cut (%e)", ci->dmin, cj->dmin,
          r_cut);
#endif

  /* ---- Update active particles in ci from all pj within r_cut. ---- */
  if (CELL_IS_ACTIVE(ci, e)) {

    /* Loop over the *active* parts in ci that can be within r_cut of any
     * part in cj (outer scan bound on the pair axis). */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];

      /* Skip inactive particles. */
      if (!PART_IS_ACTIVE(pi, e)) continue;

      const float hi = pi->h;

      /* Inner scan bound for pi: the largest d value in cj that could be
       * within r_cut of pi (accounting for particle drift). */
      const double di = sort_i[pid].d + r_cut + dx_max - rshift;
      if (di < dj_min) continue;

      /* Position of pi relative to cj origin (shifted for periodicity). */
      const float pix = (float)(pi->x[0] - (cj->loc[0] + shift[0]));
      const float piy = (float)(pi->x[1] - (cj->loc[1] + shift[1]));
      const float piz = (float)(pi->x[2] - (cj->loc[2] + shift[2]));

      /* Loop over the parts in cj that are within the inner scan bound. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Compute the pairwise distance. */
        const float pjx = (float)(pj->x[0] - cj->loc[0]);
        const float pjy = (float)(pj->x[1] - cj->loc[1]);
        const float pjz = (float)(pj->x[2] - cj->loc[2]);
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time. */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < r_cut2) {
          IACT_NONSYM_HYDRO_APERTURE(r2, dx, hi, hj, pi, pj, a, H,
                                  with_self_gravity, cosmo, e->sink_properties);
        }
      }
    }
  }

  /* ---- Update active particles in cj from all pi within r_cut. ---- */
  if (CELL_IS_ACTIVE(cj, e)) {

    /* Loop over the *active* parts in cj that can be within r_cut of any
     * part in ci (outer scan bound on the pair axis). */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *pj = &parts_j[sort_j[pjd].i];

      /* Skip inactive particles. */
      if (!PART_IS_ACTIVE(pj, e)) continue;

      const float hj = pj->h;

      /* Inner scan bound for pj: the smallest d value in ci that could
       * be within r_cut of pj. */
      const double dj = sort_j[pjd].d - r_cut - dx_max + rshift;
      if (dj - rshift > di_max) continue;

      /* Position of pj relative to cj origin. */
      const float pjx = (float)(pj->x[0] - cj->loc[0]);
      const float pjy = (float)(pj->x[1] - cj->loc[1]);
      const float pjz = (float)(pj->x[2] - cj->loc[2]);

      /* Loop over the parts in ci that are within the inner scan bound. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Get a hold of the ith part in ci. */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;

        /* Position of pi relative to cj origin. */
        const float pix = (float)(pi->x[0] - (cj->loc[0] + shift[0]));
        const float piy = (float)(pi->x[1] - (cj->loc[1] + shift[1]));
        const float piz = (float)(pi->x[2] - (cj->loc[2] + shift[2]));

        /* dx = pj - pi (pj is the receiver for the non-symmetric call). */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time. */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < r_cut2) {
          /* Note: dx points from pi to pj, so pj is the accumulator. */
          IACT_NONSYM_HYDRO_APERTURE(r2, dx, hj, hi, pj, pi, a, H,
                                  with_self_gravity, cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR_HYDRO_APERTURE);
}

/* ============================================================
 * DOPAIR1_BRANCH_HYDRO_APERTURE
 * Dispatch to sorted or naive DOPAIR1.
 * ============================================================ */

/**
 * @brief Dispatch to the appropriate DOPAIR1 variant for this cell pair.
 *
 * Computes the pair direction SID, checks that sorted lists are up to date,
 * and calls either the sorted DOPAIR1_HYDRO_APERTURE or the brute-force
 * DOPAIR1_NAIVE_HYDRO_APERTURE variant.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR1_BRANCH_HYDRO_APERTURE(struct runner *r, struct cell *ci,
                                struct cell *cj, const float r_cut) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure both cells have been drifted. */
  if (!CELL_ARE_PART_DRIFTED(ci, e) || !CELL_ARE_PART_DRIFTED(cj, e))
    error("Interacting undrifted cells.");
#endif

  /* Get the pair direction and apply the canonical cell ordering. */
#if !defined(SWIFT_USE_NAIVE_INTERACTIONS)
  double shift[3];
  const int sid = space_getsid_and_swap_cells(e->s, &ci, &cj, shift);
#endif

#ifndef SWIFT_USE_NAIVE_INTERACTIONS
  /* Sorted lists must be valid before calling the sorted variant. */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells (ci).");

  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells (cj).");
#endif

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR1_NAIVE_HYDRO_APERTURE(r, ci, cj, r_cut);
#else
  DOPAIR1_HYDRO_APERTURE(r, ci, cj, r_cut, sid, shift);
#endif
}

/* ============================================================
 * DOSELF1_HYDRO_APERTURE
 * Optimised non-symmetric self interaction with active-particle list.
 * ============================================================ */

/**
 * @brief Compute non-symmetric cell self-interactions (optimised).
 *
 * Builds a compact list @c indt[] of active particle indices to avoid
 * processing inactive-inactive pairs.  The main loop then uses two
 * regimes:
 *
 *   - If @c pi is passive: only iterate over active @c pj entries in
 *     @c indt[firstdt..countdt] that appear after @c pi in the array
 *     order, calling the non-symmetric interaction with @c pj as the
 *     accumulator.
 *
 *   - If @c pi is active: advance @c firstdt by 1 (pi has moved out of
 *     the pending-active window) and iterate all @c pj from @c pid+1 to
 *     @c count, calling the symmetric variant when both are active and
 *     the non-symmetric variant otherwise.
 *
 * Each unordered pair is visited exactly once, avoiding double-counting.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF1_HYDRO_APERTURE(struct runner *r, const struct cell *c,
                         const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  struct part *restrict parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Build a compact list of indices of active particles. */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate active-particle index array.");
  for (int k = 0; k < count; k++) {
    if (PART_IS_ACTIVE(&parts[k], e)) {
      indt[countdt] = k;
      countdt += 1;
    }
  }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

#ifdef SWIFT_DEBUG_CHECKS
  if (r_cut > c->dmin)
    error("Cell size (%e) smaller than r_cut (%e)", c->dmin, r_cut);
#endif

  /* Main loop: iterate until the active window is exhausted. */
  for (int pid = 0; pid < count && firstdt < countdt; pid++) {

    /* Get a hold of the ith part in c. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const float hi = pi->h;
    const int update_i = PART_IS_ACTIVE(pi, e);

    if (!update_i) {

      /* pi is passive.  Loop over active pj particles that come after pi in
       * array order (they are in indt[firstdt..countdt]).  We use pj as the
       * accumulator since pi cannot be updated. */
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        struct part *restrict pj = &parts[indt[pjd]];
        const float hj = pj->h;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time. */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* dx = pj - pi (pj is the receiver). */
        float dx[3] = {(float)(pj->x[0] - pi->x[0]),
                       (float)(pj->x[1] - pi->x[1]),
                       (float)(pj->x[2] - pi->x[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < r_cut2) {
          IACT_NONSYM_HYDRO_APERTURE(r2, dx, hj, hi, pj, pi, a, H,
                                  with_self_gravity, cosmo, e->sink_properties);
        }
      }

    } else {

      /* pi is active.  Advance the active window (pi has been removed from
       * the pending set) and process all pj from pid+1 to count. */
      firstdt += 1;

      for (int pjd = pid + 1; pjd < count; pjd++) {

        struct part *restrict pj = &parts[pjd];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time. */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* dx = pi - pj (pi is the receiver). */
        float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                       (float)(pi->x[1] - pj->x[1]),
                       (float)(pi->x[2] - pj->x[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < r_cut2) {
          if (PART_IS_ACTIVE(pj, e)) {
            /* Both active: symmetric call updates both pi and pj. */
            IACT_HYDRO_APERTURE(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                             cosmo, e->sink_properties);
          } else {
            /* Only pi is active. */
            IACT_NONSYM_HYDRO_APERTURE(r2, dx, hi, hj, pi, pj, a, H,
                                    with_self_gravity, cosmo,
                                    e->sink_properties);
          }
        }
      }
    }
  }

  free(indt);

  TIMER_TOC(TIMER_DOSELF_HYDRO_APERTURE);
}

/* ============================================================
 * DOSELF1_BRANCH_HYDRO_APERTURE
 * Dispatch to optimised or naive DOSELF1.
 * ============================================================ */

/**
 * @brief Dispatch to the appropriate DOSELF1 variant for this cell.
 *
 * Checks that the cell is active and drifted, then calls either the optimised
 * DOSELF1_HYDRO_APERTURE or the brute-force DOSELF1_NAIVE_HYDRO_APERTURE variant.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF1_BRANCH_HYDRO_APERTURE(struct runner *r, const struct cell *c,
                                const float r_cut) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!CELL_IS_ACTIVE(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure the cell has been drifted. */
  if (!CELL_ARE_PART_DRIFTED(c, e)) error("Interacting undrifted cell.");
#endif

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_NAIVE_HYDRO_APERTURE(r, c, r_cut);
#else
  DOSELF1_HYDRO_APERTURE(r, c, r_cut);
#endif
}

/* ============================================================
 * DOSUB_PAIR1_HYDRO_APERTURE
 * Recursive sub-cell pair interaction.
 * ============================================================ */

/**
 * @brief Recursively compute non-symmetric pair interactions for sub-cells.
 *
 * Recurses while @c r_cut < 0.5 * ci->dmin — i.e. while the fixed aperture
 * is smaller than a sub-cell.  Once the aperture equals or exceeds half the
 * cell's minimum dimension, all sub-cell progeny pairs would interact and
 * no pruning is possible, so DOPAIR1_BRANCH_HYDRO_APERTURE is called at this
 * level instead of descending further.
 *
 * This is the fixed-aperture analogue of the variable-h hydro loop's
 * @c recurse_below_h_max logic: the recursion stops when the search radius
 * covers a whole sub-cell.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_PAIR1_HYDRO_APERTURE(struct runner *r, struct cell *ci, struct cell *cj,
                             const float r_cut, const int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!CELL_IS_ACTIVE(ci, e) && !CELL_IS_ACTIVE(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the pair direction and apply the canonical cell ordering. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* We reached a leaf OR the aperture is larger than a sub-cell — no benefit
   * in recursing further since all sub-cell pairs would interact anyway.
   * For the fixed r_cut loop the threshold is r_cut >= 0.5 * ci->dmin, which
   * mirrors the cell_can_recurse_in_subpair_hydro_task condition that uses
   * h_max * kernel_gamma < 0.5 * dmin. */
  if (!ci->split || ci->hydro.count < space_recurse_size_pair_hydro ||
      !cj->split || cj->hydro.count < space_recurse_size_pair_hydro ||
      r_cut >= 0.5f * ci->dmin) {

    /* Ensure cells are sorted before calling the branch function. */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
      runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_requests_sort=*/0, /*clock=*/0);
    }
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
      runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_requests_sort=*/0, /*clock=*/0);
    }

    DOPAIR1_BRANCH_HYDRO_APERTURE(r, ci, cj, r_cut);

  } else {

    /* Both cells are split and the aperture fits within a sub-cell.
     * Recurse over the progeny pairs for the current pair direction. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_HYDRO_APERTURE(r, ci->progeny[pid], cj->progeny[pjd], r_cut,
                                /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR_HYDRO_APERTURE);
}

/* ============================================================
 * DOSUB_SELF1_HYDRO_APERTURE
 * Recursive sub-cell self interaction.
 * ============================================================ */

/**
 * @brief Recursively compute non-symmetric self interactions for sub-cells.
 *
 * Recurses while @c r_cut < 0.5 * c->dmin — i.e. while the fixed aperture
 * is smaller than a sub-cell.  At the leaf level (or once the aperture covers
 * the sub-cell scale), calls DOSELF1_BRANCH_HYDRO_APERTURE for the self part
 * and DOSUB_PAIR1_HYDRO_APERTURE for the cross-progeny part.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_SELF1_HYDRO_APERTURE(struct runner *r, struct cell *c,
                             const float r_cut, const int gettimer) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0 || !CELL_IS_ACTIVE(c, r->e)) return;

  /* We reached a leaf OR the aperture is larger than a sub-cell — call the
   * leaf self function at this level. */
  if (!c->split || c->hydro.count < space_recurse_size_self_hydro ||
      r_cut >= 0.5f * c->dmin) {

    DOSELF1_BRANCH_HYDRO_APERTURE(r, c, r_cut);

  } else {

    /* Cell is split and the aperture fits within a sub-cell.
     * Recurse: self interactions of each progeny + pair interactions between
     * every distinct pair of progeny. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        /* Self interactions within this progeny. */
        DOSUB_SELF1_HYDRO_APERTURE(r, c->progeny[k], r_cut, /*gettimer=*/0);

        /* Pair interactions with all subsequent progeny. */
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_HYDRO_APERTURE(r, c->progeny[k], c->progeny[j], r_cut,
                                    /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_HYDRO_APERTURE);
}
