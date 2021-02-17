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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "active.h"
#include "runner_doiact_sinks.h"

/**
 * @brief Calculate the number density of #part around the #sink
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_SINKS(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->sinks.count == 0) return;
  if (!cell_is_active_hydro(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->sinks.count;
  const int count = c->hydro.count;
  struct sink *restrict sinks = c->sinks.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the particles in ci. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts[pjd];
    const float hj = pj->h;

    /* Skip inactive particles */
    if (!part_is_active(pj, e)) continue;

    const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                          (float)(pj->x[1] - c->loc[1]),
                          (float)(pj->x[2] - c->loc[2])};

    /* Loop over the sinks in cj. */
    for (int sid = 0; sid < scount; sid++) {

      /* Get a hold of the ith sink in ci. */
      struct sink *restrict si = &sinks[sid];

      /* Early abort? */
      if (sink_is_inhibited(si, e)) continue;

      const float ri = si->r_cut;
      const float ri2 = ri * ri;

      const float six[3] = {(float)(si->x[0] - c->loc[0]),
                            (float)(si->x[1] - c->loc[1]),
                            (float)(si->x[2] - c->loc[2])};

      /* Compute the pairwise distance. */
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");

      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
#endif

      if (r2 < ri2) {
        IACT_SINK(r2, dx, ri, hj, si, pj, a, H);
      }
    } /* loop over the sinks in ci. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Calculate the number density of cj #part around the ci #sink
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_SINKS_NAIVE(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->hydro.count == 0 || ci->sinks.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->sinks.count;
  const int count_j = cj->hydro.count;
  struct sink *restrict sinks_i = ci->sinks.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts in cj. */
  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];
    const float hj = pj->h;

    /* Skip inactive particles */
    if (!part_is_active(pj, e)) continue;

    const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                          (float)(pj->x[1] - cj->loc[1]),
                          (float)(pj->x[2] - cj->loc[2])};

    /* Loop over the sinks in ci. */
    for (int sid = 0; sid < scount_i; sid++) {

      /* Get a hold of the ith sink in ci. */
      struct sink *restrict si = &sinks_i[sid];

      /* Skip inhibited particles. */
      if (sink_is_inhibited(si, e)) continue;

      /* Get the radius */
      const float ri = si->r_cut;
      const float ri2 = ri * ri;

      /* Compute the pairwise distance. */
      const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                            (float)(si->x[1] - (cj->loc[1] + shift[1])),
                            (float)(si->x[2] - (cj->loc[2] + shift[2]))};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");

      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
#endif

      if (r2 < ri2) {
        IACT_SINK(r2, dx, ri, hj, si, pj, a, H);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DO_SYM_PAIR1_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                        const int sid, const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  const int do_ci_sinks = (ci->nodeID == e->nodeID) && (ci->sinks.count != 0) &&
                          (cj->hydro.count != 0) && cell_is_active_hydro(ci, e);
  const int do_cj_sinks = (cj->nodeID == e->nodeID) && (cj->sinks.count != 0) &&
                          (ci->hydro.count != 0) && cell_is_active_hydro(cj, e);

  if (do_ci_sinks) {

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->sinks.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->sinks.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->sinks.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const int count_i = ci->sinks.count;
    const int count_j = cj->hydro.count;
    struct sink *restrict sinks_i = ci->sinks.parts;
    struct part *restrict parts_j = cj->hydro.parts;
    const double dj_min = sort_j[0].d;
    const float hydro_dx_max_rshift = cj->hydro.dx_max_sort - rshift;

    /* Loop over the sinks in ci. */
    for (int i = 0; i < count_i; i++) {

      /* Get a hold of the ith part in ci. */
      struct sink *restrict spi = &sinks_i[i];
      const float ri = spi->r_cut;

      /* Skip inhibited particles */
      if (sink_is_inhibited(spi, e)) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spi->x[0], spi->x[1], spi->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double di = dist + ri + hydro_dx_max_rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float ri2 = ri * ri;
      const float pix = spi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = spi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = spi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inactive particles. */
        if (!part_is_active(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
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

        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle spi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < ri2) {
          IACT_SINK(r2, dx, ri, hj, spi, pj, a, H);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* do_ci_sinks */

  if (do_cj_sinks) {
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->hydro.dx_max_part, cj->sinks.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->hydro.dx_max_part, cj->sinks.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->hydro.dx_max_part, cj->sinks.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const int count_i = ci->hydro.count;
    const int count_j = cj->sinks.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct sink *restrict sinks_j = cj->sinks.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float hydro_dx_max_rshift = ci->hydro.dx_max_sort - rshift;

    /* Loop over the sinks in cj. */
    for (int j = 0; j < count_j; j++) {

      /* Get a hold of the jth part in cj. */
      struct sink *spj = &sinks_j[j];
      const float rj = spj->r_cut;

      /* Skip inhibited particles */
      if (sink_is_inhibited(spj, e)) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spj->x[0], spj->x[1], spj->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double dj = dist - rj - hydro_dx_max_rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float rj2 = rj * rj;
      const float pjx = spj->x[0] - cj->loc[0];
      const float pjy = spj->x[1] - cj->loc[1];
      const float pjz = spj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inactive particles. */
        if (!part_is_active(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
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

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (spj->ti_drift != e->ti_current)
          error("Particle spj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < rj2) {

          IACT_SINK(r2, dx, rj, hi, spj, pi, a, H);
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */
}

void DOPAIR1_SINKS_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {

  const int do_ci_sinks = ci->nodeID == r->e->nodeID;
  const int do_cj_sinks = cj->nodeID == r->e->nodeID;
  if (do_ci_sinks && ci->sinks.count != 0 && cj->hydro.count != 0)
    DO_NONSYM_PAIR1_SINKS_NAIVE(r, ci, cj);
  if (do_cj_sinks && cj->sinks.count != 0 && ci->hydro.count != 0)
    DO_NONSYM_PAIR1_SINKS_NAIVE(r, cj, ci);
}

/**
 * @brief Determine which version of DOSELF1_SINKS needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_SINKS(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->sinks.r_cut_max_old > c->dmin)
    error("Cell smaller than the cut off radius");

  DOSELF1_SINKS(r, c, 1);
}

/**
 * @brief Determine which version of DOPAIR1_SINKS needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_SINKS needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_SINKS(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  const int ci_active = cell_is_active_hydro(ci, e);
  const int cj_active = cell_is_active_hydro(cj, e);
  const int do_ci_sinks = ci->nodeID == e->nodeID;
  const int do_cj_sinks = cj->nodeID == e->nodeID;
  const int do_ci = (ci->sinks.count != 0 && cj->hydro.count != 0 &&
                     ci_active && do_ci_sinks);
  const int do_cj = (cj->sinks.count != 0 && ci->hydro.count != 0 &&
                     cj_active && do_cj_sinks);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci && (!cell_are_sink_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_ci && (!(cj->hydro.sorted & (1 << sid)) ||
                cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

  if (do_cj && (!cell_are_part_drifted(ci, e) || !cell_are_sink_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_cj && (!(ci->hydro.sorted & (1 << sid)) ||
                ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_SINKS
  DOPAIR1_SINKS_NAIVE(r, ci, cj, 1);
#else
  DO_SYM_PAIR1_SINKS(r, ci, cj, sid, shift);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = ci->sinks.count != 0 && cj->hydro.count != 0 &&
                           cell_is_active_hydro(ci, e);
  const int should_do_cj = cj->sinks.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_hydro(cj, e);
  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_sinks_task(ci, cj) &&
      cell_can_recurse_in_pair_sinks_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_SINKS(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci_sinks = ci->nodeID == e->nodeID;
    const int do_cj_sinks = cj->nodeID == e->nodeID;
    const int do_ci = ci->sinks.count != 0 && cj->hydro.count != 0 &&
                      cell_is_active_hydro(ci, e) && do_ci_sinks;
    const int do_cj = cj->sinks.count != 0 && ci->hydro.count != 0 &&
                      cell_is_active_hydro(cj, e) && do_cj_sinks;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_sink_drifted(ci, e))
        error("Interacting undrifted cells (sinks).");

      if (!cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");

      /* Do any of the cells need to be sorted first? */
      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx)
        error("Interacting unsorted cell (parts). %i", cj->nodeID);
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (parts).");

      if (!cell_are_sink_drifted(cj, e))
        error("Interacting undrifted cells (sinks).");

      /* Do any of the cells need to be sorted first? */
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
        error("Interacting unsorted cell (parts).");
      }
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_SINKS(r, ci, cj);
  }
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SINKS(struct runner *r, struct cell *ci, int gettimer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (ci->hydro.count == 0 || ci->sinks.count == 0 ||
      !cell_is_active_hydro(ci, r->e))
    return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_SINKS(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_SINKS(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_sink_drifted(ci, r->e)) error("Interacting undrifted cell.");

    DOSELF1_BRANCH_SINKS(r, ci);
  }
}
