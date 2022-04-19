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
#include "../config.h"

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "sink.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"
#include "space_getsid.h"



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
  //const integertime_t ti_current = e->ti_current;
  //const struct cosmology *cosmo = e->cosmology;
  //const int with_cosmology = e->policy & engine_policy_cosmology;
  //const int si_is_local = 1; /* SELF tasks are always local */

  /* Anything to do here? */
  if (c->sinks.count == 0) return;
  if (!cell_is_active_sinks(c, e)) return;

  const int scount = c->sinks.count;
  const int count = c->hydro.count;
  struct sink *restrict sinks = c->sinks.parts;
  struct part *restrict parts = c->hydro.parts;
  //struct xpart *restrict xparts = c->hydro.xparts;   used by iact

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
          runner_iact_nonsym_sinks_gas_swallow(r2, dx, ri, hj, si, pj);

          //if (si_is_local) {      /* SELF tasks are always local */
          //runner_iact_nonsym_bh_gas_repos(r2, dx, hi, pj->h, bi, pj, xpj, with_cosmology, cosmo,e->gravity_properties, e->black_holes_properties, e->entropy_floor, ti_current, e->time);
          //}
        }
      } /* loop over the parts in ci. */
    }   /* loop over the bparts in ci. */
  }     /* Do we have gas particles in the cell? */

  /* When doing BH swallowing, we need a quick loop also over the BH
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
        runner_iact_nonsym_sinks_sink_swallow(r2, dx, ri, rj, si, sj);

        //if (si_is_local) {
        //  runner_iact_nonsym_bh_bh_repos(r2, dx, hi, hj, bi, bj, cosmo, e->gravity_properties, e->black_holes_properties, ti_current);
        //}
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
void runner_do_nonsym_pair_sinks_naive_swallow(struct runner *r, struct cell *restrict ci,
                              struct cell *restrict cj) {



#ifdef SWIFT_DEBUG_CHECKS
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  //const integertime_t ti_current = e->ti_current;
  //const struct cosmology *cosmo = e->cosmology;
  //const int with_cosmology = e->policy & engine_policy_cosmology;
  //const int si_is_local = ci->nodeID == e->nodeID;

  /* Anything to do here? */
  if (ci->sinks.count == 0) return;
  if (!cell_is_active_sinks(ci, e)) return;

  const int scount_i = ci->sinks.count;
  const int count_j = cj->hydro.count;
  struct sink *restrict sinks_i = ci->sinks.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  //struct xpart *restrict xparts_j = cj->hydro.xparts;

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
          runner_iact_nonsym_sinks_gas_swallow(r2, dx, ri, hj, si, pj);
          //if (si_is_local) {
          //    runner_iact_nonsym_bh_gas_repos(r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo, e->gravity_properties, e->black_holes_properties, e->entropy_floor, ti_current, e->time);
          //}
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
        runner_iact_nonsym_sinks_sink_swallow(r2, dx, ri, rj, si, sj);

        //if (si_is_local) {
        //  runner_iact_nonsym_bh_bh_repos(r2, dx, hi, hj, bi, bj, cosmo,e->gravity_properties,e->black_holes_properties, ti_current);
        //}
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
void runner_dopair_sinks_naive_swallow(struct runner *r, struct cell *restrict ci,
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
void runner_dopair_branch_sinks_swallow(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  const int ci_active = cell_is_active_sinks(ci, e);
  const int cj_active = cell_is_active_sinks(cj, e);

  const int do_ci = (ci->sinks.count != 0 && cj->hydro.count != 0 && ci_active);
  const int do_cj = (cj->sinks.count != 0 && ci->hydro.count != 0 && cj_active);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;
  
  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_sink_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_part_drifted(ci, e) || !cell_are_sink_drifted(cj, e)))
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
void runner_dosub_pair_sinks_swallow(struct runner *r, struct cell *ci, struct cell *cj,
                    int timer) {

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
        runner_dosub_pair_sinks_swallow(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci = ci->black_holes.count != 0 && cell_is_active_black_holes(ci, e);
    const int do_cj = cj->black_holes.count != 0 && cell_is_active_black_holes(cj, e);

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
void runner_dosub_self_sinks_swallow(struct runner *r, struct cell *ci, int timer) {

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

    /* Should we even bother?
     * In the swallow case we care about sink-sink and sink-gas
     * interactions. */

  const int should_do_ci =
      ci->sinks.count != 0 && cell_is_active_sinks(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        runner_dosub_self_sinks_swallow(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            runner_dosub_pair_sinks_swallow(r, ci->progeny[k], ci->progeny[j], 0);
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
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the gas particle and the gas particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the gas particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow(struct runner *r, struct cell *c, int timer) {

}

/**
 * @brief Processing of gas particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_gas_swallow_self(struct runner *r, struct cell *c, int timer) {


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

 
}

/**
 * @brief Process all the BH particles in a cell that have been flagged for
 * swallowing by a black hole.
 *
 * This is done by recursing down to the leaf-level and skipping the sub-cells
 * that have not been drifted as they would not have any particles with
 * swallowing flag. We then loop over the particles with a flag and look into
 * the space-wide list of black holes for the particle with the corresponding
 * ID. If found, the BH swallows the BH particle and the BH particle is
 * removed. If the cell is local, we may be looking for a foreign BH, in which
 * case, we do not update the BH (that will be done on its node) but just remove
 * the BH particle.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow(struct runner *r, struct cell *c, int timer) {

}

/**
 * @brief Processing of bh particles to swallow - self task case.
 *
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_self(struct runner *r, struct cell *c, int timer) {


}

/**
 * @brief Processing of bh particles to swallow - pair task case.
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_sinks_sink_swallow_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer) {


}
