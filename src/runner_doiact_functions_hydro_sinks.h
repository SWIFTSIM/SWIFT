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
 * Before including this file, define FUNCTION (e.g. "density") and
 * FUNCTION_TASK_LOOP (e.g. TASK_LOOP_SINK_FORMATION).  The macro expansion
 * creates the full set of functions:
 *   runner_dopair1_hydro_sinks_FUNCTION
 *   runner_doself1_hydro_sinks_FUNCTION
 *   runner_dosub_pair1_hydro_sinks_FUNCTION
 *   runner_dosub_self1_hydro_sinks_FUNCTION
 *   (and the symmetric _2_ variants that delegate to the _1_ variants)
 *
 * Interaction calls are dispatched via the IACT_NONSYM_HYDRO_SINKS and
 * IACT_HYDRO_SINKS macros defined in runner_doiact_hydro_sinks.h, which
 * expand to runner_iact_nonsym_hydro_sinks_FUNCTION and
 * runner_iact_hydro_sinks_FUNCTION respectively.
 *
 * Recursion termination:
 *   DOSUB_PAIR/SELF recurse while r_cut < 0.5 * cell->dmin.  Once the
 *   aperture equals or exceeds half the cell's minimum dimension, further
 *   subdivision offers no pruning benefit (all sub-cell pairs would interact
 *   anyway), so the leaf pair/self function is called at that level instead
 *   of descending to the true leaf cells. */

#include "runner_doiact_hydro_sinks.h"

/* ============================================================
 * DOPAIR1_NAIVE_HYDRO_SINKS
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
void DOPAIR1_NAIVE_HYDRO_SINKS(struct runner *r,
                                const struct cell *restrict ci,
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
        IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
      if (pj_active && r2 < r_cut2) {
        float mdx[3] = {-dx[0], -dx[1], -dx[2]};
        IACT_NONSYM_HYDRO_SINKS(r2, mdx, hj, hi, pj, pi, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR_HYDRO_SINKS);
}

/* ============================================================
 * DOPAIR2_NAIVE_HYDRO_SINKS
 * Brute-force symmetric pair interaction.
 * ============================================================ */

/**
 * @brief Compute the interactions between a cell pair (symmetric, brute force).
 *
 * Identical to DOPAIR1_NAIVE_HYDRO_SINKS but uses the symmetric interaction
 * call when both particles are active.  Used as a reference and when
 * #SWIFT_USE_NAIVE_INTERACTIONS is defined.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR2_NAIVE_HYDRO_SINKS(struct runner *r,
                                const struct cell *restrict ci,
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
      if (r2 < r_cut2) {
        if (pi_active && pj_active) {
          /* Both active: symmetric call updates both. */
          IACT_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                           cosmo, e->sink_properties);
        } else if (pi_active) {
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        } else if (pj_active) {
          float mdx[3] = {-dx[0], -dx[1], -dx[2]};
          IACT_NONSYM_HYDRO_SINKS(r2, mdx, hj, hi, pj, pi, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR_HYDRO_SINKS);
}

/* ============================================================
 * DOSELF1_NAIVE_HYDRO_SINKS
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
void DOSELF1_NAIVE_HYDRO_SINKS(struct runner *r, const struct cell *c,
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
      float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                     (float)(pi->x[1] - pj->x[1]),
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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
        if (pj_active) {
          float mdx[3] = {-dx[0], -dx[1], -dx[2]};
          IACT_NONSYM_HYDRO_SINKS(r2, mdx, hj, hi, pj, pi, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOSELF_HYDRO_SINKS);
}

/* ============================================================
 * DOSELF2_NAIVE_HYDRO_SINKS
 * Brute-force symmetric self interaction.
 * ============================================================ */

/**
 * @brief Compute cell self-interactions (symmetric, brute force).
 *
 * Identical to DOSELF1_NAIVE_HYDRO_SINKS but uses the symmetric interaction
 * call when both particles are active.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF2_NAIVE_HYDRO_SINKS(struct runner *r, const struct cell *c,
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
      float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                     (float)(pi->x[1] - pj->x[1]),
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
        if (pi_active && pj_active) {
          /* Both active: symmetric call updates both. */
          IACT_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                           cosmo, e->sink_properties);
        } else if (pi_active) {
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        } else if (pj_active) {
          float mdx[3] = {-dx[0], -dx[1], -dx[2]};
          IACT_NONSYM_HYDRO_SINKS(r2, mdx, hj, hi, pj, pi, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOSELF_HYDRO_SINKS);
}

/* ============================================================
 * DOPAIR_SUBSET_NOSORT_HYDRO_SINKS
 * Brute-force pair subset (no sort, computes shift internally).
 * ============================================================ */

/**
 * @brief Compute pair interactions for a subset of particles in @p ci (no
 *        sort).
 *
 * Interaction partners in @p cj are scanned without using the sorted particle
 * lists.  The periodic shift between cell centres is computed internally.
 *
 * @param r The #runner.
 * @param ci The first #cell (contains the active particle subset).
 * @param parts_i Particle array for @p ci.
 * @param ind Indices of the particles in the subset.
 * @param count Number of particles in the subset.
 * @param cj The second #cell (all particles interact as neighbours).
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR_SUBSET_NOSORT_HYDRO_SINKS(struct runner *r,
                                       struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       int *restrict ind, int count,
                                       struct cell *restrict cj,
                                       const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

  /* Get the periodic shift between the two cell centres. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the active particles in the subset. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part. */
    struct part *restrict pi = &parts_i[ind[pid]];
    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];
    const float hi = pi->h;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e))
      error("Inactive particle in subset function!");
#endif

    /* Loop over all the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(pix - pj->x[0]),
                           (float)(piy - pj->x[1]),
                           (float)(piz - pj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < r_cut2) {
        IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, pj->h, pi, pj, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
    }
  }

  TIMER_TOC(timer_dopair_subset);
}

/* ============================================================
 * DOPAIR_SUBSET_NAIVE_HYDRO_SINKS
 * Brute-force pair subset (with pre-computed shift).
 * ============================================================ */

/**
 * @brief Compute pair interactions for a subset of particles (brute force,
 *        pre-computed shift).
 *
 * As DOPAIR_SUBSET_NOSORT_HYDRO_SINKS but receives the periodic shift vector
 * from the caller rather than computing it internally.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param shift Periodic shift vector to apply to @p ci particles.
 */
void DOPAIR_SUBSET_NAIVE_HYDRO_SINKS(struct runner *r,
                                      const struct cell *restrict ci,
                                      struct part *restrict parts_i,
                                      const int *ind, const int count,
                                      const struct cell *restrict cj,
                                      const float r_cut,
                                      const double shift[3]) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

  /* Loop over the active particles in the subset. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part. */
    struct part *restrict pi = &parts_i[ind[pid]];
    const double pix = pi->x[0] - shift[0];
    const double piy = pi->x[1] - shift[1];
    const double piz = pi->x[2] - shift[2];
    const float hi = pi->h;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e))
      error("Inactive particle in subset function!");
#endif

    /* Loop over all the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(pix - pj->x[0]),
                           (float)(piy - pj->x[1]),
                           (float)(piz - pj->x[2])};
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
        IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, pj->h, pi, pj, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
    }
  }

  TIMER_TOC(timer_dopair_subset_naive);
}

/* ============================================================
 * DOPAIR_SUBSET_HYDRO_SINKS
 * Sorted pair subset using fixed aperture r_cut.
 * ============================================================ */

/**
 * @brief Compute pair interactions for a subset of particles using the sorted
 *        particle lists.
 *
 * Uses the axis-aligned sorted lists to bound the inner loop over @p cj
 * particles.  The scan bound replaces @c hi*kernel_gamma with @p r_cut
 * throughout, giving a fixed-aperture search.
 *
 * @param r The #runner.
 * @param ci The first #cell (contains the active particle subset).
 * @param parts_i Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param sid The direction index of the pair (identifies the sort axis).
 * @param flipped Whether the cells have been swapped from the canonical order.
 * @param shift Periodic shift vector to apply to @p ci particles.
 */
void DOPAIR_SUBSET_HYDRO_SINKS(struct runner *r, const struct cell *restrict ci,
                                struct part *restrict parts_i, const int *ind,
                                const int count,
                                const struct cell *restrict cj,
                                const float r_cut, const int sid,
                                const int flipped, const double shift[3]) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

  /* Sorted list along the pair direction for cj. */
  const struct sort_entry *sort_j = cell_get_hydro_sorts(cj, sid);
  const float dxj = cj->hydro.dx_max_sort;

  if (!flipped) {

    /* Loop over the active particles in the subset. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - shift[0];
      const double piy = pi->x[1] - shift[1];
      const double piz = pi->x[2] - shift[2];
      const float hi = pi->h;

      /* Scan bound: fixed r_cut replaces hi * kernel_gamma. */
      const double di = r_cut + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj within the scan bound. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float dx[3] = {(float)(pix - pj->x[0]),
                             (float)(piy - pj->x[1]),
                             (float)(piz - pj->x[2])};
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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, pj->h, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }

  } else {

    /* Flipped case: ci and cj are in reverse canonical order. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - shift[0];
      const double piy = pi->x[1] - shift[1];
      const double piz = pi->x[2] - shift[2];
      const float hi = pi->h;

      /* Scan bound for the flipped pair direction. */
      const double di = -r_cut - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj within the scan bound. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float dx[3] = {(float)(pix - pj->x[0]),
                             (float)(piy - pj->x[1]),
                             (float)(piz - pj->x[2])};
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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, pj->h, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(timer_dopair_subset);
}

/* ============================================================
 * DOPAIR_SUBSET_BRANCH_HYDRO_SINKS
 * Branch: choose sorted or naive pair subset depending on sort availability.
 * ============================================================ */

/**
 * @brief Choose the sorted or naive pair subset function depending on whether
 *        the sorted lists for @p cj are available and up to date.
 *
 * @param r The #runner.
 * @param ci The first #cell (contains the active particle subset).
 * @param parts_i Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR_SUBSET_BRANCH_HYDRO_SINKS(struct runner *r,
                                       const struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *ind, const int count,
                                       struct cell *restrict cj,
                                       const float r_cut) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->hydro.count == 0) return;

  /* Compute the periodic shift and the pair direction index. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Compute the SID without cell swapping. */
  int sid = 0;
  for (int k = 0; k < 3; k++) {
    const double d = cj->loc[k] - ci->loc[k] + shift[k];
    sid = 3 * sid + (d < 0 ? 0 : (d > 0 ? 2 : 1));
  }
  const int flipped = runner_flip[sid];
  sid = sortlistID[sid];

  /* Is the sorted list valid? */
  const int is_sorted =
      (cj->hydro.sorted & (1 << sid)) &&
      (cj->hydro.dx_max_sort_old <= space_maxreldx * cj->dmin);

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR_SUBSET_NAIVE_HYDRO_SINKS(r, ci, parts_i, ind, count, cj, r_cut,
                                   shift);
#else
  if (!is_sorted) {
    DOPAIR_SUBSET_NAIVE_HYDRO_SINKS(r, ci, parts_i, ind, count, cj, r_cut,
                                     shift);
  } else {
    DOPAIR_SUBSET_HYDRO_SINKS(r, ci, parts_i, ind, count, cj, r_cut, sid,
                               flipped, shift);
  }
#endif
}

/* ============================================================
 * DOSELF_SUBSET_HYDRO_SINKS
 * Self subset interaction.
 * ============================================================ */

/**
 * @brief Compute self interactions for a subset of particles in @p c.
 *
 * Each active particle in the subset interacts with all non-inhibited particles
 * in the full cell @p c (including itself, which is skipped by the pi == pj
 * guard).
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param parts Particle array for @p c.
 * @param ind Indices of the active particles in the subset.
 * @param count Number of particles in the subset.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF_SUBSET_HYDRO_SINKS(struct runner *r, const struct cell *c,
                                struct part *restrict parts, const int *ind,
                                const int count, const float r_cut) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float r_cut2 = r_cut * r_cut;

  TIMER_TIC;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);

  const int count_cell = c->hydro.count;
  struct part *restrict parts_all = c->hydro.parts;

  /* Loop over the active particles in the subset. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part. */
    struct part *restrict pi = &parts[ind[pid]];
    const float hi = pi->h;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over all the parts in c. */
    for (int pjd = 0; pjd < count_cell; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_all[pjd];

      /* Skip self-interaction and inhibited particles. */
      if (pi == pj) continue;
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pi->x[0] - pj->x[0]),
                     (float)(pi->x[1] - pj->x[1]),
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
        IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, pj->h, pi, pj, a, H, with_self_gravity,
                                cosmo, e->sink_properties);
      }
    }
  }

  TIMER_TOC(timer_doself_subset);
}

/* ============================================================
 * DOSELF_SUBSET_BRANCH_HYDRO_SINKS
 * Branch for self subset (dispatches to DOSELF_SUBSET_HYDRO_SINKS).
 * ============================================================ */

/**
 * @brief Dispatch to DOSELF_SUBSET_HYDRO_SINKS.
 *
 * @param r The #runner.
 * @param ci The #cell.
 * @param parts Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF_SUBSET_BRANCH_HYDRO_SINKS(struct runner *r,
                                       const struct cell *ci,
                                       struct part *restrict parts,
                                       const int *ind, const int count,
                                       const float r_cut) {
  DOSELF_SUBSET_HYDRO_SINKS(r, ci, parts, ind, count, r_cut);
}

/* ============================================================
 * DOPAIR1_HYDRO_SINKS
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
void DOPAIR1_HYDRO_SINKS(struct runner *r, const struct cell *restrict ci,
                          const struct cell *restrict cj,
                          const float r_cut, const int sid,
                          const double shift[3]) {

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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  /* ---- Update active particles in cj from all pi within r_cut. ---- */
  if (CELL_IS_ACTIVE(cj, e)) {

    /* Loop over the *active* parts in cj that can be within r_cut of any
     * part in ci (outer scan bound on the pair axis). */
    for (int pjd = 0;
         pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max; pjd++) {

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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hj, hi, pj, pi, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
        }
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR_HYDRO_SINKS);
}

/* ============================================================
 * DOPAIR1_BRANCH_HYDRO_SINKS
 * Dispatch to sorted or naive DOPAIR1.
 * ============================================================ */

/**
 * @brief Dispatch to the appropriate DOPAIR1 variant for this cell pair.
 *
 * Computes the pair direction SID, checks that sorted lists are up to date,
 * and calls either the sorted DOPAIR1_HYDRO_SINKS or the brute-force
 * DOPAIR1_NAIVE_HYDRO_SINKS variant.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR1_BRANCH_HYDRO_SINKS(struct runner *r, struct cell *ci,
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
  double shift[3];
  const int sid = space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

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
  DOPAIR1_NAIVE_HYDRO_SINKS(r, ci, cj, r_cut);
#else
  DOPAIR1_HYDRO_SINKS(r, ci, cj, r_cut, sid, shift);
#endif
}

/* ============================================================
 * DOPAIR2_HYDRO_SINKS / DOPAIR2_BRANCH_HYDRO_SINKS
 * Symmetric variants — delegate to the non-symmetric versions.
 * ============================================================ */

/**
 * @brief Compute symmetric pair interactions using sorted particle lists.
 *
 * Delegates to DOPAIR1_HYDRO_SINKS since the sink formation gather uses
 * non-symmetric semantics: each active particle independently accumulates
 * contributions from its neighbours, so a second symmetric sweep is
 * not required.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param sid The direction index of the pair.
 * @param shift The shift vector.
 */
void DOPAIR2_HYDRO_SINKS(struct runner *r, const struct cell *restrict ci,
                          const struct cell *restrict cj,
                          const float r_cut, const int sid,
                          const double shift[3]) {
  DOPAIR1_HYDRO_SINKS(r, ci, cj, r_cut, sid, shift);
}

/**
 * @brief Dispatch to the appropriate DOPAIR2 variant for this cell pair.
 *
 * Delegates to DOPAIR1_BRANCH_HYDRO_SINKS.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOPAIR2_BRANCH_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                 struct cell *cj, const float r_cut) {
  DOPAIR1_BRANCH_HYDRO_SINKS(r, ci, cj, r_cut);
}

/* ============================================================
 * DOSELF1_HYDRO_SINKS
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
void DOSELF1_HYDRO_SINKS(struct runner *r, const struct cell *c,
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
          IACT_NONSYM_HYDRO_SINKS(r2, dx, hj, hi, pj, pi, a, H, with_self_gravity,
                                  cosmo, e->sink_properties);
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
            IACT_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity, cosmo,
                             e->sink_properties);
          } else {
            /* Only pi is active. */
            IACT_NONSYM_HYDRO_SINKS(r2, dx, hi, hj, pi, pj, a, H, with_self_gravity, cosmo,
                                    e->sink_properties);
          }
        }
      }
    }
  }

  free(indt);

  TIMER_TOC(TIMER_DOSELF_HYDRO_SINKS);
}

/* ============================================================
 * DOSELF1_BRANCH_HYDRO_SINKS
 * Dispatch to sorted or naive DOSELF1.
 * ============================================================ */

/**
 * @brief Dispatch to the appropriate DOSELF1 variant for this cell.
 *
 * Checks that the cell is active and drifted, then calls either the optimised
 * DOSELF1_HYDRO_SINKS or the brute-force DOSELF1_NAIVE_HYDRO_SINKS variant.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF1_BRANCH_HYDRO_SINKS(struct runner *r, const struct cell *c,
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
  DOSELF1_NAIVE_HYDRO_SINKS(r, c, r_cut);
#else
  DOSELF1_HYDRO_SINKS(r, c, r_cut);
#endif
}

/* ============================================================
 * DOSELF2_HYDRO_SINKS / DOSELF2_BRANCH_HYDRO_SINKS
 * Symmetric variants — delegate to DOSELF1.
 * ============================================================ */

/**
 * @brief Compute symmetric cell self-interactions (optimised).
 *
 * Delegates to DOSELF1_HYDRO_SINKS since the sink formation gather uses
 * non-symmetric semantics.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF2_HYDRO_SINKS(struct runner *r, const struct cell *c,
                          const float r_cut) {
  DOSELF1_HYDRO_SINKS(r, c, r_cut);
}

/**
 * @brief Dispatch to the appropriate DOSELF2 variant for this cell.
 *
 * Delegates to DOSELF1_BRANCH_HYDRO_SINKS.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 */
void DOSELF2_BRANCH_HYDRO_SINKS(struct runner *r, const struct cell *c,
                                 const float r_cut) {
  DOSELF1_BRANCH_HYDRO_SINKS(r, c, r_cut);
}

/* ============================================================
 * FIND_SUB_HYDRO_SINKS
 * Find which sub-cell of c contains parts[ind[0]].
 * ============================================================ */

/**
 * @brief Find the sub-cell of @p c that contains the first particle in the
 *        subset @p ind.
 *
 * Walks the eight progeny of @p c and returns the one whose particle array
 * contains @c parts[ind[0]].
 *
 * @param c The parent #cell (must be split).
 * @param parts Particle array (same as @c c->hydro.parts base).
 * @param ind Indices into @p parts; @c ind[0] identifies the particle.
 * @return Pointer to the matching progeny #cell.
 */
struct cell *FIND_SUB_HYDRO_SINKS(const struct cell *const c,
                                   const struct part *const parts,
                                   const int *const ind) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!c->split) error("Can't search for sub-cells in a non-split cell.");
#endif

  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL) {
      const struct cell *const prog = c->progeny[k];
      if (&parts[ind[0]] >= &prog->hydro.parts[0] &&
          &parts[ind[0]] < &prog->hydro.parts[prog->hydro.count]) {
        return c->progeny[k];
      }
    }
  }
  error("Could not find a matching sub-cell — invalid particle index.");
  return NULL;
}

/* ============================================================
 * DOSUB_PAIR_SUBSET_HYDRO_SINKS
 * Recursive pair subset.
 * ============================================================ */

/**
 * @brief Recursively compute pair interactions for a subset of @p ci
 *        particles against all particles in @p cj.
 *
 * Recurses while @c r_cut < 0.5 * ci->dmin (the aperture fits within a
 * sub-cell), descending along the pair direction SID. At the leaf level,
 * calls DOPAIR_SUBSET_BRANCH_HYDRO_SINKS.
 *
 * @param r The #runner.
 * @param ci The first #cell (contains the active particle subset).
 * @param parts Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_PAIR_SUBSET_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                    struct part *parts, const int *ind,
                                    const int count, struct cell *cj,
                                    const float r_cut, const int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse if both cells are split and the aperture still fits within a
   * sub-cell (r_cut < 0.5 * dmin).  Beyond this threshold all sub-cell pairs
   * would interact anyway, so calling the leaf function at this level is both
   * correct and more efficient. */
  if (ci->split && cj->split && r_cut < 0.5f * ci->dmin &&
      r_cut < 0.5f * cj->dmin) {

    /* Find which progeny of ci contains the first particle in the subset. */
    struct cell *const sub = FIND_SUB_HYDRO_SINKS(ci, parts, ind);

    /* Determine the pair direction SID (apply canonical cell ordering). */
    double shift[3];
    struct cell *ci_tmp = ci, *cj_tmp = cj;
    const int sid = space_getsid_and_swap_cells(s, &ci_tmp, &cj_tmp, shift);

    /* Recurse over the progeny pairs along this direction. */
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci_tmp->progeny[pid] == sub && cj_tmp->progeny[pjd] != NULL)
        DOSUB_PAIR_SUBSET_HYDRO_SINKS(r, ci_tmp->progeny[pid], parts, ind,
                                       count, cj_tmp->progeny[pjd], r_cut,
                                       /*gettimer=*/0);
      if (cj_tmp->progeny[pjd] == sub && ci_tmp->progeny[pid] != NULL)
        DOSUB_PAIR_SUBSET_HYDRO_SINKS(r, cj_tmp->progeny[pjd], parts, ind,
                                       count, ci_tmp->progeny[pid], r_cut,
                                       /*gettimer=*/0);
    }

  } else {

    /* Reached a leaf or the aperture covers the sub-cell scale. */
#ifdef SWIFT_DEBUG_CHECKS
    if (!cell_are_part_drifted(cj, e)) error("Cell cj should be drifted!");
#endif
    DOPAIR_SUBSET_BRANCH_HYDRO_SINKS(r, ci, parts, ind, count, cj, r_cut);
  }

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}

/* ============================================================
 * DOSUB_SELF_SUBSET_HYDRO_SINKS
 * Recursive self subset.
 * ============================================================ */

/**
 * @brief Recursively compute self interactions for a subset of @p ci
 *        particles.
 *
 * Recurses while @c r_cut < 0.5 * ci->dmin, then at the leaf level calls
 * DOSELF_SUBSET_BRANCH_HYDRO_SINKS.
 *
 * @param r The #runner.
 * @param ci The #cell.
 * @param parts Particle array for @p ci.
 * @param ind Indices of the active particles.
 * @param count Number of particles in the subset.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_SELF_SUBSET_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                    struct part *parts, const int *ind,
                                    const int count, const float r_cut,
                                    const int gettimer) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse while the aperture fits within a sub-cell. */
  if (ci->split && r_cut < 0.5f * ci->dmin) {

    /* Find which progeny of ci contains the first particle in the subset. */
    struct cell *const sub = FIND_SUB_HYDRO_SINKS(ci, parts, ind);

    /* Self-interaction of the sub-cell containing the subset. */
    DOSUB_SELF_SUBSET_HYDRO_SINKS(r, sub, parts, ind, count, r_cut,
                                   /*gettimer=*/0);

    /* Pair interactions between the sub-cell and all other progeny. */
    for (int j = 0; j < 8; j++) {
      if (ci->progeny[j] != NULL && ci->progeny[j] != sub) {
        DOSUB_PAIR_SUBSET_HYDRO_SINKS(r, sub, parts, ind, count,
                                       ci->progeny[j], r_cut, /*gettimer=*/0);
      }
    }

  } else {

    /* Reached a leaf or the aperture covers the sub-cell scale. */
    DOSELF_SUBSET_BRANCH_HYDRO_SINKS(r, ci, parts, ind, count, r_cut);
  }

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}

/* ============================================================
 * DOSUB_PAIR1_HYDRO_SINKS
 * Recursive sub-cell pair interaction.
 * ============================================================ */

/**
 * @brief Recursively compute non-symmetric pair interactions for sub-cells.
 *
 * Recurses while @c r_cut < 0.5 * ci->dmin — i.e. while the fixed aperture
 * is smaller than a sub-cell.  Once the aperture equals or exceeds half the
 * cell's minimum dimension, all sub-cell progeny pairs would interact and
 * no pruning is possible, so DOPAIR1_BRANCH_HYDRO_SINKS is called at this
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
void DOSUB_PAIR1_HYDRO_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
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

    DOPAIR1_BRANCH_HYDRO_SINKS(r, ci, cj, r_cut);

  } else {

    /* Both cells are split and the aperture fits within a sub-cell.
     * Recurse over the progeny pairs for the current pair direction. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_HYDRO_SINKS(r, ci->progeny[pid], cj->progeny[pjd], r_cut,
                                 /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR_HYDRO_SINKS);
}

/* ============================================================
 * DOSUB_SELF1_HYDRO_SINKS
 * Recursive sub-cell self interaction.
 * ============================================================ */

/**
 * @brief Recursively compute non-symmetric self interactions for sub-cells.
 *
 * Recurses while @c r_cut < 0.5 * c->dmin — i.e. while the fixed aperture
 * is smaller than a sub-cell.  At the leaf level (or once the aperture covers
 * the sub-cell scale), calls DOSELF1_BRANCH_HYDRO_SINKS for the self part
 * and DOSUB_PAIR1_HYDRO_SINKS for the cross-progeny part.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_SELF1_HYDRO_SINKS(struct runner *r, struct cell *c,
                              const float r_cut, const int gettimer) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0 || !CELL_IS_ACTIVE(c, r->e)) return;

  /* We reached a leaf OR the aperture is larger than a sub-cell — call the
   * leaf self function at this level. */
  if (!c->split || c->hydro.count < space_recurse_size_self_hydro ||
      r_cut >= 0.5f * c->dmin) {

    DOSELF1_BRANCH_HYDRO_SINKS(r, c, r_cut);

  } else {

    /* Both cell is split and the aperture fits within a sub-cell.
     * Recurse: self interactions of each progeny + pair interactions between
     * every distinct pair of progeny. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        /* Self interactions within this progeny. */
        DOSUB_SELF1_HYDRO_SINKS(r, c->progeny[k], r_cut, /*gettimer=*/0);

        /* Pair interactions with all subsequent progeny. */
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_HYDRO_SINKS(r, c->progeny[k], c->progeny[j], r_cut,
                                     /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_HYDRO_SINKS);
}

/* ============================================================
 * DOSUB_PAIR2_HYDRO_SINKS / DOSUB_SELF2_HYDRO_SINKS
 * Symmetric sub-cell variants — delegate to non-symmetric.
 * ============================================================ */

/**
 * @brief Recursively compute symmetric pair interactions for sub-cells.
 *
 * Delegates to DOSUB_PAIR1_HYDRO_SINKS since sink formation uses
 * non-symmetric gather semantics.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_PAIR2_HYDRO_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                              const float r_cut, const int gettimer) {
  DOSUB_PAIR1_HYDRO_SINKS(r, ci, cj, r_cut, gettimer);
}

/**
 * @brief Recursively compute symmetric self interactions for sub-cells.
 *
 * Delegates to DOSUB_SELF1_HYDRO_SINKS since sink formation uses
 * non-symmetric gather semantics.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param r_cut The fixed aperture radius for interactions.
 * @param gettimer Whether to record a timer for this call.
 */
void DOSUB_SELF2_HYDRO_SINKS(struct runner *r, struct cell *c,
                              const float r_cut, const int gettimer) {
  DOSUB_SELF1_HYDRO_SINKS(r, c, r_cut, gettimer);
}
