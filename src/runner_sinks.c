/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "sink.h"
#include "space_getsid.h"
#include "timers.h"

/**
 * @brief Calculate gas and sink interaction around #sinks
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_doself_sinks_swallow(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;
  if (!cell_is_active_sinks(c, e)) return;

  const int scount = c->sinks.count;
  const int count = c->hydro.count;
  struct sink *restrict sinks = c->sinks.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Do we actually have any gas neighbours? */
  if (c->hydro.count != 0) {

    /* Loop over the sinks in ci. */
    for (int sid = 0; sid < scount; sid++) {

      /* Get a hold of the ith sinks in ci. */
      struct sink *restrict si = &sinks[sid];

      /* Skip inactive particles */
      if (!sink_is_active(si, e)) continue;

      const float ri = si->r_cut;
      const float ri2 = ri * ri;
      const float six[3] = {(float)(si->x[0] - c->loc[0]),
                            (float)(si->x[1] - c->loc[1]),
                            (float)(si->x[2] - c->loc[2])};

      /* Loop over the parts (gas) in cj. */
      for (int pjd = 0; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];
        const float hj = pj->h;

        /* Early abort? */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                              (float)(pj->x[1] - c->loc[1]),
                              (float)(pj->x[2] - c->loc[2])};
        const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < ri2) {
          runner_iact_nonsym_sinks_gas_swallow(r2, dx, ri, hj, si, pj, with_cosmology, cosmo,
					       e->gravity_properties,e->sink_properties);
        }
      } /* loop over the parts in ci. */
    }   /* loop over the bparts in ci. */
  }     /* Do we have gas particles in the cell? */

  /* When doing sink swallowing, we need a quick loop also over the sink
   * neighbours */

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith sink in ci. */
    struct sink *restrict si = &sinks[sid];

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount; sjd++) {

      /* Skip self interaction */
      if (sid == sjd) continue;

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks[sjd];
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

      /* Early abort? */
      if (sink_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - c->loc[0]),
                            (float)(sj->x[1] - c->loc[1]),
                            (float)(sj->x[2] - c->loc[2])};
      const float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (si->ti_drift != e->ti_current)
        error("Particle bi not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle bj not drifted to current time");
#endif

      if (r2 < ri2 || r2 < rj2) {
        runner_iact_nonsym_sinks_sink_swallow(r2, dx, ri, rj, si, sj, with_cosmology,
					      cosmo, e->gravity_properties); 
      }
    } /* loop over the sinks in ci. */
  }   /* loop over the sinks in ci. */

  if (timer) TIMER_TOC(timer_doself_sink_swallow);
}

/**
 * @brief Calculate gas and sink interaction around #sinks
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void runner_do_nonsym_pair_sinks_naive_swallow(struct runner *r,
                                               struct cell *restrict ci,
                                               struct cell *restrict cj) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  /* Anything to do here? */
  if (ci->sinks.count == 0) return;
  if (!cell_is_active_sinks(ci, e)) return;

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

  /* Do we actually have any gas neighbours? */
  if (cj->hydro.count != 0) {

    /* Loop over the sinks in ci. */
    for (int sid = 0; sid < scount_i; sid++) {

      /* Get a hold of the ith bpart in ci. */
      struct sink *restrict si = &sinks_i[sid];

      /* Skip inactive particles */
      if (!sink_is_active(si, e)) continue;

      const float ri = si->r_cut;
      const float ri2 = ri * ri;
      const float six[3] = {(float)(si->x[0] - cj->loc[0]),
                            (float)(si->x[1] - cj->loc[1]),
                            (float)(si->x[2] - cj->loc[2])};

      /* Loop over the parts (gas) in cj. */
      for (int pjd = 0; pjd < count_j; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[pjd];
        const float hj = pj->h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                              (float)(pj->x[1] - cj->loc[1]),
                              (float)(pj->x[2] - cj->loc[2])};
        const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < ri2) {
          runner_iact_nonsym_sinks_gas_swallow(r2, dx, ri, hj, si, pj, with_cosmology, cosmo,
					       e->gravity_properties, e->sink_properties);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the sinks in ci. */
  }     /* Do we have gas particles in the cell? */

  /* When doing sink swallowing, we need a quick loop also over the sinks
   * neighbours */

  const int scount_j = cj->sinks.count;
  struct sink *restrict sinks_j = cj->sinks.parts;

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith bpart in ci. */
    struct sink *restrict si = &sinks_i[sid];

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks_j[sjd];
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

      /* Skip inhibited particles. */
      if (sink_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - cj->loc[0]),
                            (float)(sj->x[1] - cj->loc[1]),
                            (float)(sj->x[2] - cj->loc[2])};
      const float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < ri2 || r2 < rj2) {
        runner_iact_nonsym_sinks_sink_swallow(r2, dx, ri, rj, si, sj, with_cosmology,
					      cosmo, e->gravity_properties);
      }
    } /* loop over the sinks in cj. */
  }   /* loop over the sinks in ci. */
}

/**
 * @brief Calculate swallow for ci #sinks part around the cj #gas and sinks and
 *                              cj #sinks part around the ci #gas and sinks
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void runner_dopair_sinks_naive_swallow(struct runner *r,
                                       struct cell *restrict ci,
                                       struct cell *restrict cj, int timer) {

  TIMER_TIC;

  runner_do_nonsym_pair_sinks_naive_swallow(r, ci, cj);
  runner_do_nonsym_pair_sinks_naive_swallow(r, cj, ci);

  if (timer) TIMER_TOC(timer_dopair_sink_swallow);
}

/**
 * @brief Wrapper to runner_doself_sinks_swallow
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void runner_doself_branch_sinks_swallow(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_sinks(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->sinks.r_cut_max_old > c->dmin)
    error("Cell smaller than the cut off radius");

  runner_doself_sinks_swallow(r, c, 1);
}

/**
 * @brief Wrapper for runner_dopair_sinks_naive_swallow.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void runner_dopair_branch_sinks_swallow(struct runner *r, struct cell *ci,
                                        struct cell *cj) {

  const struct engine *restrict e = r->e;

  const int ci_active = cell_is_active_sinks(ci, e);
  const int cj_active = cell_is_active_sinks(cj, e);

  const int do_ci = (ci->sinks.count != 0 && cj->hydro.count != 0 && ci_active);
  const int do_cj = (cj->sinks.count != 0 && ci->hydro.count != 0 && cj_active);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci && (!cell_are_sink_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj && (!cell_are_part_drifted(ci, e) || !cell_are_sink_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* No sorted interactions here -> use the naive ones */
  runner_dopair_sinks_naive_swallow(r, ci, cj, 1);
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
void runner_dosub_pair_sinks_swallow(struct runner *r, struct cell *ci,
                                     struct cell *cj, int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother?
   * In the swallow case we care about sink-sink and sink-gas
   * interactions. */

  const int should_do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);
  const int should_do_cj = cj->sinks.count != 0 && cell_is_active_sinks(cj, e);

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
        runner_dosub_pair_sinks_swallow(r, ci->progeny[pid], cj->progeny[pjd],
                                        0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);
    const int do_cj = cj->sinks.count != 0 && cell_is_active_sinks(cj, e);

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_sink_drifted(ci, e))
        error("Interacting undrifted cells (sinks).");

      if (cj->hydro.count != 0 && !cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (parts).");

      if (!cell_are_sink_drifted(cj, e))
        error("Interacting undrifted cells (sinks).");
    }

    if (do_ci || do_cj) runner_dopair_branch_sinks_swallow(r, ci, cj);
  }

  if (timer) TIMER_TOC(timer_dosub_pair_sink_swallow);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void runner_dosub_self_sinks_swallow(struct runner *r, struct cell *ci,
                                     int timer) {

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother?
   * In the swallow case we care about sink-sink and sink-gas
   * interactions. */

  const int should_do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        runner_dosub_self_sinks_swallow(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            runner_dosub_pair_sinks_swallow(r, ci->progeny[k], ci->progeny[j],
                                            0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_sink_drifted(ci, e)) error("Interacting undrifted cell.");

    if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
      error("Interacting undrifted cells (parts).");

    runner_doself_branch_sinks_swallow(r, ci);
  }

  if (timer) TIMER_TOC(timer_dosub_self_sink_swallow);
}

/**
 * @brief Process all the gas particles in a cell that have been flagged for
 * swallowing by a sink.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of sink for the particle with the corresponding
 * ID. If found, the sink swallows the gas particle and the gas particle is
 * removed. If the cell is local, we may be looking for a foreign sink, in which
 * case, we do not update the sink (that will be done on its node) but just
 * remove the gas particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;

  struct sink *sinks = s->sinks;
  const size_t nr_sink = s->nr_sinks;
#ifdef WITH_MPI
  error("MPI is not implemented yet for sink particles.");
#endif

  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;

  timebin_t max_bin = e->max_active_bin;
  integertime_t ti_current = e->ti_current;
  integertime_t ti_beg_max = 0;
  int count = 0;

  /* Early abort?
   * (We only want cells for which we drifted the gas as these are
   * the only ones that could have gas particles that have been flagged
   * for swallowing) */
  if (c->hydro.count == 0 || c->hydro.ti_old_part != e->ti_current) {
    return;
  }

  /* Loop over the progeny ? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        runner_do_sinks_gas_swallow(r, cp, 0);

        /* Propagate the ti_beg_max from the leaves to the roots.
         * See bug fix below. */
        ti_beg_max = max(cp->hydro.ti_beg_max, ti_beg_max);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell
     * Note that the cell (and hence the parts) may be local or foreign. */
    const size_t nr_parts = c->hydro.count;
    for (size_t k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Ignore inhibited particles (they have already been removed!) */
      if (part_is_inhibited(p, e)) continue;

      /* Get the ID of the sink that will swallow this part */
      const long long swallow_id = sink_get_part_swallow_id(&p->sink_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (p->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the sink swallowing this particle */
        const long long sink_id = swallow_id;

        /* Have we found this particle's sink already? */
        int found = 0;

        /* Let's look for the hungry sink in the local list */
        for (size_t i = 0; i < nr_sink; ++i) {

          /* Get a handle on the bpart. */
          struct sink *sp = &sinks[i];

          if (sp->id == sink_id) {

            /* Lock the space as we are going to work directly on the spart list
             */
            lock_lock(&s->lock);

            /* Swallow the gas particle (i.e. update the sink properties) */
            sink_swallow_part(sp, p, xp, e->cosmology);

            /* Release the space as we are done updating the spart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

            /* If the gas particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              lock_lock(&e->s->lock);

              /* Re-check that the particle has not been removed
               * by another thread before we do the deed. */
              if (!part_is_inhibited(p, e)) {

                /* Finally, remove the gas particle from the system
                 * Recall that the gpart associated with it is also removed
                 * at the same time. */
                cell_remove_part(e, c, p, xp);
              }

              if (lock_unlock(&e->s->lock) != 0)
                error("Failed to unlock the space!");
            }

            /* In any case, prevent the particle from being re-swallowed */
            sink_mark_part_as_swallowed(&p->sink_data);

            found = 1;
            break;
          }

        } /* Loop over local sinks */

#ifdef WITH_MPI
        error("MPI is not implemented yet for sink particles.");
#endif

        /* If we have a local particle, we must have found the sink in one
         * of our list of sinks. */
        if (c->nodeID == e->nodeID && !found) {
          error("Gas particle %lld could not find sink %lld to be swallowed",
                p->id, swallow_id);
        }
      } /* Part was flagged for swallowing */

      /* Bug fix : Change the hydro.ti_beg_max when a sink eats the last gas
       * particle possessing the ti_beg_max of the cell. We set hydro.ti_beg_max
       * to the max ti_beg of the remaining gas particle. Why this fix ?
       * Otherwise, we fail the check from cell_check_timesteps. This bug is
       * rare because it needs that the swallowed gas is the last part with the
       * ti_beg_max of the cell.
       * The same is not done for ti_end_min since it may inactivate cells that
       * need to perform sinks tasks.
       */

      /* A part may habe been swallowed just before so continue if this is the
         case */
      if (part_is_inhibited(p, e)) continue;

      ++count;

      integertime_t ti_beg;

      if (p->time_bin <= max_bin) {
        ti_beg = get_integer_time_begin(ti_current + 1, p->time_bin);
      } else {
        ti_beg = get_integer_time_begin(ti_current + 1, p->time_bin);
      }

      ti_beg_max = max(ti_beg, ti_beg_max);
    } /* Loop over the parts */
  }   /* Cell is not split */

  /* Update ti_beg_max. See bug fix above. */
  if (ti_beg_max != c->hydro.ti_beg_max) {
    c->hydro.ti_beg_max = ti_beg_max;
  }
}

/**
 * @brief Processing of gas particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_self(struct runner *r, struct cell *c,
                                      int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_sinks(c, r->e) && !cell_is_active_hydro(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_sinks_gas_swallow(r, c, timer);
}

/**
 * @brief Processing of gas particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_pair(struct runner *r, struct cell *ci,
                                      struct cell *cj, int timer) {

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active sink */
  if (cell_is_active_sinks(cj, e)) runner_do_sinks_gas_swallow(r, ci, timer);
  if (cell_is_active_sinks(ci, e)) runner_do_sinks_gas_swallow(r, cj, timer);
}

/**
 * @brief Process all the sink particles in a cell that have been flagged for
 * swallowing by a sink.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of sinks for the particle with the corresponding
 * ID. If found, the sink swallows the sink particle and the sink particle is
 * removed. If the cell is local, we may be looking for a foreign sink, in which
 * case, we do not update the sink (that will be done on its node) but just
 * remove the sink particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  struct space *s = e->s;

  struct sink *sinks = s->sinks;
  const size_t nr_sink = s->nr_sinks;
#ifdef WITH_MPI
  error("MPI is not implemented yet for sink particles.");
#endif

  struct sink *cell_sinks = c->sinks.parts;

  /* Early abort?
   * (We only want cells for which we drifted the sink as these are
   * the only ones that could have sink particles that have been flagged
   * for swallowing) */
  if (c->sinks.count == 0 || c->sinks.ti_old_part != e->ti_current) {
    return;
  }

  /* Loop over the progeny ? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        runner_do_sinks_sink_swallow(r, cp, 0);
      }
    }
  } else {

    /* Loop over all the sinks particles in the cell
     * Note that the cell (and hence the sinks) may be local or foreign. */
    const size_t nr_cell_sinks = c->sinks.count;

    for (size_t k = 0; k < nr_cell_sinks; k++) {

      /* Get a handle on the part. */
      struct sink *const cell_sp = &cell_sinks[k];

      /* Ignore inhibited particles (they have already been removed!) */
      if (sink_is_inhibited(cell_sp, e)) continue;

      /* Get the ID of the sink that will swallow this sink */
      const long long swallow_id =
          sink_get_sink_swallow_id(&cell_sp->merger_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

#ifdef SWIFT_DEBUG_CHECKS
        if (cell_sp->ti_drift != e->ti_current)
          error("Trying to swallow an un-drifted particle.");
#endif

        /* ID of the sink swallowing this particle */
        const long long sink_id = swallow_id;

        /* Have we found this particle's sink already? */
        int found = 0;

        /* Let's look for the hungry sink in the local list */
        for (size_t i = 0; i < nr_sink; ++i) {

          /* Get a handle on the bpart. */
          struct sink *sp = &sinks[i];

          if (sp->id == sink_id) {

            /* Is the swallowing sink itself flagged for swallowing by
               another sink? */
            if (sink_get_sink_swallow_id(&sp->merger_data) != -1) {

              /* Pretend it was found and abort */
              sink_mark_sink_as_not_swallowed(&cell_sp->merger_data);
              found = 1;
              break;
            }

            /* Lock the space as we are going to work directly on the
             * space's bpart list */
            lock_lock(&s->lock);

            /* Swallow the sink particle (i.e. update the swallowing sink
             * properties with the properties of cell_sp) */
            sink_swallow_sink(sp, cell_sp, e->cosmology);

            /* Release the space as we are done updating the spart */
            if (lock_unlock(&s->lock) != 0)
              error("Failed to unlock the space.");

            // message("sink %lld swallowing sink particle %lld", sp->id,
            // cell_sp->id);

            /* If the sink particle is local, remove it */
            if (c->nodeID == e->nodeID) {

              /* Finally, remove the sink particle from the system
               * Recall that the gpart associated with it is also removed
               * at the same time. */
              cell_remove_sink(e, c, cell_sp);
            }

            /* In any case, prevent the particle from being re-swallowed */
            sink_mark_sink_as_merged(&cell_sp->merger_data);

            found = 1;
            break;
          }

        } /* Loop over local sinks */

#ifdef WITH_MPI
        error("MPI is not implemented yet for sink particles.");
#endif

        /* If we have a local particle, we must have found the sink in one
         * of our list of sinks. */
        if (c->nodeID == e->nodeID && !found) {
          error("sink particle %lld could not find sink %lld to be swallowed",
                cell_sp->id, swallow_id);
        }

      } /* Part was flagged for swallowing */
    }   /* Loop over the parts */
  }     /* Cell is not split */
}

/**
 * @brief Processing of sink particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_self(struct runner *r, struct cell *c,
                                       int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Running self task on foreign node");
  if (!cell_is_active_sinks(c, r->e) && !cell_is_active_hydro(c, r->e))
    error("Running self task on inactive cell");
#endif

  runner_do_sinks_sink_swallow(r, c, timer);
}

/**
 * @brief Processing of sink particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_pair(struct runner *r, struct cell *ci,
                                       struct cell *cj, int timer) {

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Running pair task on foreign node");
#endif

  /* Run the swallowing loop only in the cell that is the neighbour of the
   * active sink */
  if (cell_is_active_sinks(cj, e)) runner_do_sinks_sink_swallow(r, ci, timer);
  if (cell_is_active_sinks(ci, e)) runner_do_sinks_sink_swallow(r, cj, timer);
}
