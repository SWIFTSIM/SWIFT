/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Rob McGibbon (mcgibbon@strw.leidenuniv.nl)
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

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "runner_doiact_bh_stars.h"

/* Local headers. */
#include "active.h"
#include "black_holes_iact.h"
#include "cell.h"
#include "engine.h"
#include "runner.h"
#include "space_getsid.h"
#include "timers.h"

/* The loops in this file only exist for black hole models that carry the
 * star-density fields. For other models we only provide stubs that abort:
 * the corresponding tasks are never created in that configuration. */
#ifdef BLACK_HOLES_HAVE_STAR_DENSITY

/**
 * @brief Calculate the star density around the #bpart in a single cell.
 *
 * The black holes gather from the star particles (non-symmetric).
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_doself_bh_stars_density(struct runner *r, struct cell *c,
                                    int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->black_holes.count == 0 || c->stars.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  const int bcount = c->black_holes.count;
  const int scount = c->stars.count;
  struct bpart *restrict bparts = c->black_holes.parts;
  struct spart *restrict sparts = c->stars.parts;

  /* Loop over the bparts in c. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith bpart in c. */
    struct bpart *restrict bi = &bparts[bid];

    /* Skip inhibited particles */
    if (bpart_is_inhibited(bi, e)) continue;

    /* Skip inactive particles */
    if (!bpart_is_active(bi, e)) continue;

    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;
    const float bix[3] = {(float)(bi->x[0] - c->loc[0]),
                          (float)(bi->x[1] - c->loc[1]),
                          (float)(bi->x[2] - c->loc[2])};

    /* Loop over the sparts in c. */
    for (int sjd = 0; sjd < scount; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts[sjd];

      /* Early abort? */
      if (spart_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - c->loc[0]),
                            (float)(sj->x[1] - c->loc[1]),
                            (float)(sj->x[2] - c->loc[2])};
      const float dx[3] = {bix[0] - sjx[0], bix[1] - sjx[1], bix[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (bi->ti_drift != e->ti_current)
        error("Particle bi not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    } /* loop over the sparts in c. */
  } /* loop over the bparts in c. */

  if (timer) TIMER_TOC(timer_doself_bh_stars_density);
}

/**
 * @brief Compute the star density around the #bpart of one cell from the
 * star particles of another cell, without using the sort arrays.
 *
 * @param r The #runner.
 * @param ci The #cell holding the black holes.
 * @param cj The #cell holding the star particles.
 * @param shift The periodic wrapping to apply to the particles in ci.
 */
static void runner_dopair_bh_stars_density_naive_one_sided(
    struct runner *r, struct cell *restrict ci, struct cell *restrict cj,
    const double shift[3]) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->black_holes.count == 0 || cj->stars.count == 0) return;
  if (ci->nodeID != e->nodeID) return;
  if (!cell_is_active_black_holes(ci, e)) return;

  const int bcount_i = ci->black_holes.count;
  const int scount_j = cj->stars.count;
  struct bpart *restrict bparts_i = ci->black_holes.parts;
  struct spart *restrict sparts_j = cj->stars.parts;

  /* Loop over the bparts in ci. */
  for (int bid = 0; bid < bcount_i; bid++) {

    /* Get a hold of the ith bpart in ci. */
    struct bpart *restrict bi = &bparts_i[bid];

    /* Skip inhibited particles */
    if (bpart_is_inhibited(bi, e)) continue;

    /* Skip inactive particles */
    if (!bpart_is_active(bi, e)) continue;

    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;
    const double bix = bi->x[0] - shift[0];
    const double biy = bi->x[1] - shift[1];
    const double biz = bi->x[2] - shift[2];

    /* Loop over the sparts in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts_j[sjd];

      /* Skip inhibited particles. */
      if (spart_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(bix - sj->x[0]), (float)(biy - sj->x[1]),
                           (float)(biz - sj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (bi->ti_drift != e->ti_current)
        error("Particle bi not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    } /* loop over the sparts in cj. */
  } /* loop over the bparts in ci. */
}

/**
 * @brief Compute the star density around the #bpart of both cells of a pair.
 *
 * The black holes are not sorted (they are rare), so this simply runs the
 * naive one-sided loop once in each direction.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param timer 1 if the time is to be recorded.
 */
void runner_dopair_bh_stars_density(struct runner *r, struct cell *ci,
                                    struct cell *cj, int timer) {

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  const int do_ci = (ci->nodeID == e->nodeID) &&
                    (ci->black_holes.count != 0) && (cj->stars.count != 0) &&
                    cell_is_active_black_holes(ci, e);
  const int do_cj = (cj->nodeID == e->nodeID) &&
                    (cj->black_holes.count != 0) && (ci->stars.count != 0) &&
                    cell_is_active_black_holes(cj, e);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_bpart_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_bpart_drifted(cj, e) || !cell_are_spart_drifted(ci, e)))
    error("Interacting undrifted cells.");

  /* The black holes are not sorted, so we always use the naive loops. Note
   * that runner_dopair_bh_stars_density_naive_one_sided() applies the wrapping
   * shift to the cell holding the black holes, i.e. the ci argument. */
  runner_dopair_bh_stars_density_naive_one_sided(r, ci, cj, shift);
  const double shift_inv[3] = {-shift[0], -shift[1], -shift[2]};
  runner_dopair_bh_stars_density_naive_one_sided(r, cj, ci, shift_inv);

  if (timer) TIMER_TOC(timer_dopair_bh_stars_density);
}

/**
 * @brief Compute the star density around a subset of the #bpart of a cell,
 * using the star particles of the same cell.
 *
 * Used by the ghost task when the search radius of some black holes did not
 * converge. Always naive: this concerns very few particles.
 *
 * @param r The #runner.
 * @param ci The #cell.
 * @param bparts The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 */
void runner_doself_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts,
                                                  int *ind, const int bcount) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;

  const int scount_i = ci->stars.count;
  struct spart *restrict sparts_j = ci->stars.parts;

  /* Early abort? */
  if (scount_i == 0) return;

  /* Loop over the subset of bparts. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *bi = &bparts[ind[bid]];
    const float bix[3] = {(float)(bi->x[0] - ci->loc[0]),
                          (float)(bi->x[1] - ci->loc[1]),
                          (float)(bi->x[2] - ci->loc[2])};
    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bi, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the sparts in ci. */
    for (int sjd = 0; sjd < scount_i; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts_j[sjd];

      /* Skip inhibited particles */
      if (spart_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - ci->loc[0]),
                            (float)(sj->x[1] - ci->loc[1]),
                            (float)(sj->x[2] - ci->loc[2])};
      const float dx[3] = {bix[0] - sjx[0], bix[1] - sjx[1], bix[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    } /* loop over the sparts in ci. */
  } /* loop over the subset of bparts. */
}

/**
 * @brief Compute the star density around a subset of the #bpart of a cell,
 * using the star particles of another cell.
 *
 * Used by the ghost task when the search radius of some black holes did not
 * converge. Always naive: this concerns very few particles.
 *
 * @param r The #runner.
 * @param ci The #cell holding the black holes.
 * @param bparts_i The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The #cell holding the star particles.
 */
void runner_dopair_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts_i,
                                                  int *ind, const int bcount,
                                                  struct cell *cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;

  const int scount_j = cj->stars.count;
  struct spart *restrict sparts_j = cj->stars.parts;

  /* Early abort? */
  if (scount_j == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the subset of bparts. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *bi = &bparts_i[ind[bid]];

    const double bix = bi->x[0] - shift[0];
    const double biy = bi->x[1] - shift[1];
    const double biz = bi->x[2] - shift[2];
    const float hi = bi->h_star;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bi, e))
      error("Trying to correct the search radius of an inactive particle !");
#endif

    /* Loop over the sparts in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct spart *restrict sj = &sparts_j[sjd];

      /* Skip inhibited particles */
      if (spart_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(bix - sj->x[0]), (float)(biy - sj->x[1]),
                           (float)(biz - sj->x[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        runner_iact_nonsym_bh_stars_density(r2, dx, hi, bi, sj);
      }
    } /* loop over the sparts in cj. */
  } /* loop over the subset of bparts. */
}

/**
 * @brief Recursively compute the star density around a subset of the #bpart
 * of a cell.
 *
 * @param r The #runner.
 * @param ci The #cell holding the black holes of interest.
 * @param bparts The #bpart array of the (top) cell the subset belongs to.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The second #cell (NULL for a self interaction).
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_subset_bh_stars_density(struct runner *r, struct cell *ci,
                                          struct bpart *bparts, int *ind,
                                          const int bcount, struct cell *cj,
                                          int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (!cell_is_active_black_holes(ci, e) &&
      (cj == NULL || !cell_is_active_black_holes(cj, e)))
    return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&bparts[ind[0]] >= &ci->progeny[k]->black_holes.parts[0] &&
            &bparts[ind[0]] <
                &ci->progeny[k]
                     ->black_holes.parts[ci->progeny[k]->black_holes.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_bh_stars_task(ci)) {

      /* Loop over all progeny. */
      runner_dosub_subset_bh_stars_density(r, sub, bparts, ind, bcount, NULL,
                                           0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          runner_dosub_subset_bh_stars_density(r, sub, bparts, ind, bcount,
                                               ci->progeny[j], 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      runner_doself_subset_branch_bh_stars_density(r, ci, bparts, ind, bcount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_bh_stars_task(ci, cj) &&
        cell_can_recurse_in_pair_bh_stars_task(cj, ci)) {

      /* Get the type of pair and flip ci/cj if needed. */
      double shift[3] = {0.0, 0.0, 0.0};
      const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

      struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
          runner_dosub_subset_bh_stars_density(r, ci->progeny[pid], bparts, ind,
                                               bcount, cj->progeny[pjd], 0);
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
          runner_dosub_subset_bh_stars_density(r, cj->progeny[pjd], bparts, ind,
                                               bcount, ci->progeny[pid], 0);
      }
    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_black_holes(ci, e) && cj->stars.count > 0) {

      /* Do any of the cells need to be drifted first? */
      if (!cell_are_bpart_drifted(ci, e)) error("Cell should be drifted!");
      if (!cell_are_spart_drifted(cj, e)) error("Cell should be drifted!");

      runner_dopair_subset_branch_bh_stars_density(r, ci, bparts, ind, bcount,
                                                   cj);
    }

  } /* otherwise, pair interaction. */
}

/**
 * @brief Determine whether DOSELF needs to be called at all.
 *
 * @param r #runner
 * @param c #cell c
 */
void runner_doself_branch_bh_stars_density(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->black_holes.count == 0 || c->stars.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_black_holes(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->black_holes.h_star_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than the star search radius");

  runner_doself_bh_stars_density(r, c, 1);
}

/**
 * @brief Determine whether DOPAIR needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 */
void runner_dopair_branch_bh_stars_density(struct runner *r, struct cell *ci,
                                           struct cell *cj) {

  const struct engine *restrict e = r->e;

  const int do_ci =
      (ci->black_holes.count != 0 && cj->stars.count != 0 &&
       cell_is_active_black_holes(ci, e) && ci->nodeID == e->nodeID);
  const int do_cj =
      (cj->black_holes.count != 0 && ci->stars.count != 0 &&
       cell_is_active_black_holes(cj, e) && cj->nodeID == e->nodeID);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  runner_dopair_bh_stars_density(r, ci, cj, 1);
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_pair_bh_stars_density(struct runner *r, struct cell *ci,
                                        struct cell *cj, int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = ci->black_holes.count != 0 &&
                           cj->stars.count != 0 &&
                           cell_is_active_black_holes(ci, e);
  const int should_do_cj = cj->black_holes.count != 0 &&
                           ci->stars.count != 0 &&
                           cell_is_active_black_holes(cj, e);

  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_bh_stars_task(ci, cj) &&
      cell_can_recurse_in_pair_bh_stars_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        runner_dosub_pair_bh_stars_density(r, ci->progeny[pid],
                                           cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci = ci->black_holes.count != 0 && cj->stars.count != 0 &&
                      cell_is_active_black_holes(ci, e) &&
                      ci->nodeID == e->nodeID;
    const int do_cj = cj->black_holes.count != 0 && ci->stars.count != 0 &&
                      cell_is_active_black_holes(cj, e) &&
                      cj->nodeID == e->nodeID;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_bpart_drifted(ci, e))
        error("Interacting undrifted cells (bparts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_bpart_drifted(cj, e))
        error("Interacting undrifted cells (bparts).");
    }

    if (do_ci || do_cj) runner_dopair_branch_bh_stars_density(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(timer_dosub_pair_bh_stars_density);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_self_bh_stars_density(struct runner *r, struct cell *ci,
                                        int gettimer) {

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  const int should_do_ci = ci->black_holes.count != 0 &&
                           ci->stars.count != 0 &&
                           cell_is_active_black_holes(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_bh_stars_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        runner_dosub_self_bh_stars_density(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            runner_dosub_pair_bh_stars_density(r, ci->progeny[k],
                                               ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_bpart_drifted(ci, e)) error("Interacting undrifted cell.");

    if (!cell_are_spart_drifted(ci, e))
      error("Interacting undrifted cell (sparts).");

    runner_doself_branch_bh_stars_density(r, ci);
  }

  if (gettimer) TIMER_TOC(timer_dosub_self_bh_stars_density);
}

#else /* BLACK_HOLES_HAVE_STAR_DENSITY */

/* The selected black hole model does not support the star-density loops.
 * The tasks calling these functions are never created in that case. */

void runner_doself_bh_stars_density(struct runner *r, struct cell *c,
                                    int timer) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dopair_bh_stars_density(struct runner *r, struct cell *ci,
                                    struct cell *cj, int timer) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_doself_branch_bh_stars_density(struct runner *r, struct cell *c) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dopair_branch_bh_stars_density(struct runner *r, struct cell *ci,
                                           struct cell *cj) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_doself_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts,
                                                  int *ind, const int bcount) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dopair_subset_branch_bh_stars_density(struct runner *r,
                                                  struct cell *ci,
                                                  struct bpart *bparts_i,
                                                  int *ind, const int bcount,
                                                  struct cell *cj) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dosub_subset_bh_stars_density(struct runner *r, struct cell *ci,
                                          struct bpart *bparts, int *ind,
                                          const int bcount, struct cell *cj,
                                          int gettimer) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dosub_self_bh_stars_density(struct runner *r, struct cell *c,
                                        int gettimer) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}
void runner_dosub_pair_bh_stars_density(struct runner *r, struct cell *ci,
                                        struct cell *cj, int gettimer) {
  error("SWIFT was compiled with a black hole model without star-density support.");
}

#endif /* BLACK_HOLES_HAVE_STAR_DENSITY */
