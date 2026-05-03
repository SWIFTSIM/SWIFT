/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "fof.h"
#include "gravity.h"
#include "hydro.h"
#include "runner_doiact_sinks.h"
#include "space.h"
#include "star_formation.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"

/**
 * @brief Recurse to find the appropriate level to ionize gas particles due to
 * stellar radiation around stars.
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer) {
#ifdef IONIZATION_FEEDBACK_LOOP
  struct engine *e = r->e;

  /* TODO: Or h_max_old? */
  /* TODO: Which cell h_hii_max shall I look at? The current code changes the
     value as we do deeper in the cell hierarchy. Ideally, we want to use the
     smallest hmax for a cell, process the stars and parts at this level for
     this cell, and go higher for cells/stars with higher h_hii */
  const float r_hii_max = c->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;
  const int can_recurse = c->split && (interaction_limit < 0.5f * c->dmin);

  if (c->stars.count == 0 || c->hydro.count == 0 || r_hii_max <= 0.f ||
      !cell_is_active_stars(c, e))
    return;

  message("[%lld], hydro super = %lld", c->cellID, c->hydro.super->cellID);

  TIMER_TIC;

  /*
   * Decision: Should we process this cell or go deeper?
   * We stop recursing if:
   * 1. The cell is a leaf (!c->split).
   * 2. The cell is small enough that we can't efficiently split the
   *    interaction radius further (similar to Mosaics logic).
   * 3. The R_HII_max covers a significant fraction of the cell.
   */
  /* Is the cell split and not smaller than the smoothing length? */
  if (can_recurse) {
    /* Keep recursing deeper into the super-cell hierarchy. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        runner_do_stars_hii_ionization_feedback(r, c->progeny[k], 0);
    }
  } else {
    /* We have reached the 'Working Level' */
    runner_do_stars_hii_ionization_feedback_branch(r, c);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Ionize gas particles due to stellar radiation around stars.
 *
 * @param r The #runner thread.
 * @param c The #cell.
 */
void runner_do_stars_hii_ionization_feedback_branch(struct runner *r,
                                                    struct cell *c) {

  struct engine *e = r->e;
  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  /* OR h_max_old? */
  const float r_hii_max = c->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;

  if (c->stars.count == 0 || c->hydro.count == 0 || r_hii_max <= 0.f ||
      !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > c->dmin)
    error("Cell smaller than HII interaction length");
#endif

  message("[c = %lld, super = %lld] interaction_limit = %e, cell_dmin = %e",
          c->cellID, c->hydro.super->cellID, interaction_limit, c->dmin);

  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[sid];

    /* Is this part within the timestep? */
    if (!spart_is_active(si, e) && !feedback_is_active(si, e) &&
        feedback_is_HII_ionization_active(si, e))
      return;

    /***************************************************/
    /* First loop over particles in the current cell */
    runner_do_stars_hii_ionization_feedback_self(r, c, si);

    /***************************************************/
    /* Now loop over particles in the neighboring cells */

    /* Climb up the cell hierarchy. */
    for (struct cell *finger = c; finger != NULL; finger = finger->parent) {

      /* These are defined at the super level... When we reach the progeny,
         these will be NULL... So we need to grab the finger first */
      for (struct link *l = finger->stars.density; l != NULL; l = l->next) {
        /* We have already handled the self case */
        if (l->t->type == task_type_self) continue;

        struct cell *c_in = l->t->cj;

        /* Now, find the correct level... */
        const int can_recurse_j =
            c_in->split && (interaction_limit < 0.5f * c_in->dmin);
        if (can_recurse_j) {
          /* Keep recursing deeper into the super-cell hierarchy. */
          for (int k = 0; k < 8; k++)
            if (c_in->progeny[k] != NULL)
              runner_do_stars_hii_ionization_feedback_pair(
                  r, c, c_in->progeny[k], si);
        } else {
          /* We have reached the 'Working Level' */
          runner_do_stars_hii_ionization_feedback_pair(r, c, c_in, si);
        } /* Recurse */
      } /* Neighbour search */
    } /* Climb up in the cell hierarchy */

    /***************************************************/
    /* Now sort the gas particles */

    /***************************************************/
    /* Now let's ionize the gas particles! */

  } /* Loop over sparts */
}

void runner_do_stars_hii_ionization_feedback_self(struct runner *r,
                                                  struct cell *c,
                                                  struct spart *si) {

  struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  const float hi = si->h_hii;
  const float hig2 = hi * hi * kernel_gamma2;
  const float six[3] = {(float)(si->x[0] - c->loc[0]),
                        (float)(si->x[1] - c->loc[1]),
                        (float)(si->x[2] - c->loc[2])};

  /* Loop over the parts in c. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = &xparts[pjd];

    /* const float hj = pj->h; */

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                          (float)(pj->x[1] - c->loc[1]),
                          (float)(pj->x[2] - c->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if (r2 < hig2) {
      /* Gather */

    }
  } /* Loop in current cell */
}

void runner_do_stars_hii_ionization_feedback_pair(struct runner *r,
                                                  struct cell *ci,
                                                  struct cell *cj,
                                                  struct spart *si) {

  struct engine *e = r->e;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  /* OR h_max_old? */
  const float r_hii_max = ci->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > cj->dmin)
    error("Cell smaller than HII interaction length");

  /* Call the function that will do the work */
  if (interaction_limit > cj->dmin) {
    warning(
        "[c = %lld, cj = %lld, star = %lld] cj size is too small. We need to "
        "go up in the cell hierarchy!",
        ci->cellID, cj->cellID, si->id);
  }

  message(
      "[c = %lld, cj = %lld, cj->super = %lld, star = %lld] Neighbour search, "
      "interaction_limit = %e, cell_dmin = %e",
      ci->cellID, cj->cellID, cj->hydro.super->cellID, si->id,
      interaction_limit, cj->dmin);
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  const float hi = si->h_hii;
  const float hig2 = hi * hi * kernel_gamma2;
  const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                        (float)(si->x[1] - (cj->loc[1] + shift[1])),
                        (float)(si->x[2] - (cj->loc[2] + shift[2]))};

  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];
    struct xpart *restrict xpj = &xparts_j[pjd];

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                          (float)(pj->x[1] - cj->loc[1]),
                          (float)(pj->x[2] - cj->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if (r2 < hig2) {
      /* Gather */
    }
  }
}
