/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_DOIACT_GRAV_H
#define SWIFT_RUNNER_DOIACT_GRAV_H

/* Includes. */
#include "cell.h"
#include "clocks.h"
#include "part.h"

/**
 * @brief Compute the sorted gravity interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void runner_dopair_grav_new(struct runner *r, struct cell *ci,
                            struct cell *cj) {

  struct engine *restrict e = r->e;
  int pid, pjd, k, sid;
  double rshift, shift[3] = {0.0, 0.0, 0.0}, nshift[3];
  struct entry *restrict sort_i, *restrict sort_j;
  struct gpart *restrict pi, *restrict pj, *restrict parts_i, *restrict parts_j;
  double pix[3];
  float dx[3], r2, h_max, di, dj;
  int count_i, count_j, cnj, cnj_new;
  const int ti_current = e->ti_current;
  struct multipole m;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct gpart *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the sort ID. */
  sid = space_getsid(e->s, &ci, &cj, shift);

  /* Make sure the cells are sorted. */
  runner_dogsort(r, ci, (1 << sid), 0);
  runner_dogsort(r, cj, (1 << sid), 0);

  /* Have the cells been sorted? */
  if (!(ci->gsorted & (1 << sid)) || !(cj->gsorted & (1 << sid)))
    error("Trying to interact unsorted cells.");

  /* Get the cutoff shift. */
  for (rshift = 0.0, k = 0; k < 3; k++)
    rshift += shift[k] * runner_shift[3 * sid + k];

  /* Pick-out the sorted lists. */
  sort_i = &ci->gsort[sid * (ci->count + 1)];
  sort_j = &cj->gsort[sid * (cj->count + 1)];

  /* Get some other useful values. */
  h_max =
      sqrtf(ci->h[0] * ci->h[0] + ci->h[1] * ci->h[1] + ci->h[2] * ci->h[2]) *
      const_theta_max;
  count_i = ci->gcount;
  count_j = cj->gcount;
  parts_i = ci->gparts;
  parts_j = cj->gparts;
  cnj = count_j;
  multipole_reset(&m);
  nshift[0] = -shift[0];
  nshift[1] = -shift[1];
  nshift[2] = -shift[2];

  /* Loop over the parts in ci. */
  for (pid = count_i - 1; pid >= 0; pid--) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[sort_i[pid].i];
    if (pi->ti_end > ti_current) continue;
    di = sort_i[pid].d + h_max - rshift;

    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < cnj && sort_j[pjd].d < di; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[sort_j[pjd].i];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#ifndef VECTORIZE

      // if ( pi->part->id == 3473472412525 || pj->part->id == 3473472412525 )
      //     message( "interacting particles pi=%lli and pj=%lli with r=%.3e in
      // cells %lli/%lli." , pi->part->id , pj->part->id , sqrtf(r2) , ((long
      // long int)ci) / sizeof(struct cell) , ((long long int)cj) /
      // sizeof(struct cell) );

      runner_iact_grav(r2, dx, pi, pj);

#else

      /* Add this interaction to the queue. */
      r2q[icount] = r2;
      dxq[3 * icount + 0] = dx[0];
      dxq[3 * icount + 1] = dx[1];
      dxq[3 * icount + 2] = dx[2];
      piq[icount] = pi;
      pjq[icount] = pj;
      icount += 1;

      /* Flush? */
      if (icount == VEC_SIZE) {
        runner_iact_vec_grav(r2q, dxq, piq, pjq);
        icount = 0;
      }

#endif

    } /* loop over the parts in cj. */

    /* Set the new limit. */
    cnj_new = pjd;

    /* Add trailing parts to the multipole. */
    for (pjd = cnj_new; pjd < cnj; pjd++) {

      /* Add the part to the multipole. */
      multipole_addpart(&m, &parts_j[sort_j[pjd].i]);

    } /* add trailing parts to the multipole. */

    /* Set the new cnj. */
    cnj = cnj_new;

    /* Interact the ith particle with the multipole. */
    multipole_iact_mp(&m, pi, nshift);

  } /* loop over the parts in ci. */

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      runner_iact_grav(r2q[k], &dxq[3 * k], piq[k], pjq[k]);
#endif

  /* Re-set the multipole. */
  multipole_reset(&m);

  /* Loop over the parts in cj and interact with the multipole in ci. */
  for (pid = count_i - 1, pjd = 0; pjd < count_j; pjd++) {

    /* Get the position of pj along the axis. */
    dj = sort_j[pjd].d - h_max + rshift;

    /* Add any left-over parts in cell_i to the multipole. */
    while (pid >= 0 && sort_i[pid].d < dj) {

      /* Add this particle to the multipole. */
      multipole_addpart(&m, &parts_i[sort_i[pid].i]);

      /* Decrease pid. */
      pid -= 1;
    }

    /* Interact pj with the multipole. */
    multipole_iact_mp(&m, &parts_j[sort_j[pjd].i], shift);

  } /* loop over the parts in cj and interact with the multipole. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the recursive upward sweep, i.e. construct the
 *        multipoles in a cell hierarchy.
 *
 * @param r The #runner.
 * @param c The top-level #cell.
 */

void runner_dograv_up(struct runner *r, struct cell *c) {

  /* Re-set this cell's multipole. */
  multipole_reset(&c->multipole);

  /* Split? */
  if (c->split) {

    /* Recurse. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_dograv_up(r, c->progeny[k]);

    /* Collect the multipoles from the progeny. */
    multipole_reset(&c->multipole);
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        multipole_merge(&c->multipole, &c->progeny[k]->multipole);

  }

  /* No, leaf node. */
  else

    /* Just collect the multipole. */
    multipole_init(&c->multipole, c->gparts, c->gcount);
}

/**
 * @brief Compute the recursive downward sweep, i.e. apply the multipole
 *        acceleration on all the particles.
 *
 * @param r The #runner.
 * @param c The top-level #cell.
 */

void runner_dograv_down(struct runner *r, struct cell *c) {

  struct multipole *m = &c->multipole;

  /* Split? */
  if (c->split) {

    /* Apply this cell's acceleration on the multipoles below. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct multipole *mp = &c->progeny[k]->multipole;
        mp->a[0] += m->a[0];
        mp->a[1] += m->a[1];
        mp->a[2] += m->a[2];
      }

    /* Recurse. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_dograv_down(r, c->progeny[k]);

  }

  /* No, leaf node. */
  else {

    /* Apply the multipole acceleration to all gparts. */
    for (int k = 0; k < c->gcount; k++) {
      struct gpart *p = &c->gparts[k];
      p->a_grav[0] += m->a[0];
      p->a_grav[1] += m->a[1];
      p->a_grav[2] += m->a[2];
    }
  }
}

/**
 * @brief Compute the multipole-multipole interaction between two cells.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void runner_dograv_mm(struct runner *r, struct cell *restrict ci,
                      struct cell *restrict cj) {

  struct engine *e = r->e;
  int k;
  double shift[3] = {0.0, 0.0, 0.0};
  float dx[3], theta;

  /* Compute the shift between the cells. */
  for (k = 0; k < 3; k++) {
    dx[k] = cj->loc[k] - ci->loc[k];
    if (r->e->s->periodic) {
      if (dx[k] < -e->s->dim[k] / 2)
        shift[k] = e->s->dim[k];
      else if (dx[k] > e->s->dim[k] / 2)
        shift[k] = -e->s->dim[k];
      dx[k] += shift[k];
    }
  }
  theta =
      sqrt((dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) /
           (ci->h[0] * ci->h[0] + ci->h[1] * ci->h[1] + ci->h[2] * ci->h[2]));

  /* Do an MM or an MP/PM? */
  if (theta > const_theta_max * 4) {

    /* Update the multipoles. */
    multipole_iact_mm(&ci->multipole, &cj->multipole, shift);

  } else {

    /* Interact the multipoles via their parts. */
    for (k = 0; k < ci->gcount; k++)
      multipole_iact_mp(&cj->multipole, &ci->gparts[k], shift);
    for (k = 0; k < cj->gcount; k++)
      multipole_iact_mp(&ci->multipole, &cj->gparts[k], shift);
  }
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void runner_dopair_grav(struct runner *r, struct cell *restrict ci,
                        struct cell *restrict cj) {

  struct engine *e = r->e;
  int pid, pjd, k, count_i = ci->gcount, count_j = cj->gcount;
  double shift[3] = {0.0, 0.0, 0.0};
  struct gpart *restrict parts_i = ci->gparts, *restrict parts_j = cj->gparts;
  struct gpart *restrict pi, *restrict pj;
  double pix[3];
  float dx[3], r2;
  const int ti_current = r->e->ti_current;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct gpart *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the relative distance between the pairs, wrapping. */
  if (e->s->periodic)
    for (k = 0; k < 3; k++) {
      if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
        shift[k] = e->s->dim[k];
      else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
        shift[k] = -e->s->dim[k];
    }

  /* Loop over the parts in ci. */
  for (pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[pid];
    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

/* Compute the interaction. */
#ifndef VECTORIZE

      // if ( pi->part->id == 3473472412525 || pj->part->id == 3473472412525 )
      //     message( "interacting particles pi=%lli and pj=%lli with r=%.3e in
      // cells %lli/%lli." , pi->part->id , pj->part->id , sqrtf(r2) , ((long
      // long int)ci) / sizeof(struct cell) , ((long long int)cj) /
      // sizeof(struct cell) );

      runner_iact_grav(r2, dx, pi, pj);

#else

      /* Add this interaction to the queue. */
      r2q[icount] = r2;
      dxq[3 * icount + 0] = dx[0];
      dxq[3 * icount + 1] = dx[1];
      dxq[3 * icount + 2] = dx[2];
      piq[icount] = pi;
      pjq[icount] = pj;
      icount += 1;

      /* Flush? */
      if (icount == VEC_SIZE) {
        runner_iact_vec_grav(r2q, dxq, piq, pjq);
        icount = 0;
      }

#endif

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      runner_iact_grav(r2q[k], &dxq[3 * k], piq[k], pjq[k]);
#endif

  TIMER_TOC(timer_dopair_grav);
}

/**
 * @brief Compute the interactions within a cell.
 *
 * @param r The #runner.
 * @param c The #cell.
 */

void runner_doself_grav(struct runner *r, struct cell *restrict c) {

  int pid, pjd, k, count = c->gcount;
  struct gpart *restrict parts = c->gparts;
  struct gpart *restrict pi, *restrict pj;
  double pix[3] = {0.0, 0.0, 0.0};
  float dx[3], r2;
  const int ti_current = r->e->ti_current;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct gpart *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (c->ti_end_min > ti_current) return;

  /* Loop over every part in c. */
  for (pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts[pid];
    for (k = 0; k < 3; k++) pix[k] = pi->x[k];

    /* Loop over every other part in c. */
    for (pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

/* Compute the interaction. */
#ifndef VECTORIZE

      // if ( pi->part->id == 3473472412525 || pj->part->id == 3473472412525 )
      //     message( "interacting particles pi=%lli and pj=%lli with r=%.3e." ,
      // pi->part->id , pj->part->id , sqrtf(r2) );

      runner_iact_grav(r2, dx, pi, pj);

#else

      /* Add this interaction to the queue. */
      r2q[icount] = r2;
      dxq[3 * icount + 0] = dx[0];
      dxq[3 * icount + 1] = dx[1];
      dxq[3 * icount + 2] = dx[2];
      piq[icount] = pi;
      pjq[icount] = pj;
      icount += 1;

      /* Flush? */
      if (icount == VEC_SIZE) {
        runner_iact_vec_grav(r2q, dxq, piq, pjq);
        icount = 0;
      }

#endif

    } /* loop over the remaining parts in c. */

  } /* loop over the parts in c. */

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      runner_iact_grav(r2q[k], &dxq[3 * k], piq[k], pjq[k]);
#endif

  TIMER_TOC(timer_doself_grav);
}

/**
 * @brief Compute a gravity sub-task.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Flag to record timer or not.
 */

void runner_dosub_grav(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer) {

  int j, k, periodic = r->e->s->periodic;
  struct space *s = r->e->s;

  TIMER_TIC

  /* Self-interaction? */
  if (cj == NULL) {

    /* If the cell is split, recurse. */
    if (ci->split) {

      /* Split this task into tasks on its progeny. */
      for (j = 0; j < 8; j++)
        if (ci->progeny[j] != NULL) {
          runner_dosub_grav(r, ci->progeny[j], NULL, 0);
          for (k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              runner_dosub_grav(r, ci->progeny[j], ci->progeny[k], 0);
        }

    }

    /* Otherwise, just make a pp task out of it. */
    else
      runner_doself_grav(r, ci);

  }

  /* Nope, pair. */
  else {

    /* Get the opening angle theta. */
    float dx[3], theta;
    for (k = 0; k < 3; k++) {
      dx[k] = fabs(ci->loc[k] - cj->loc[k]);
      if (periodic && dx[k] > 0.5 * s->dim[k]) dx[k] = -dx[k] + s->dim[k];
      if (dx[k] > 0.0f) dx[k] -= ci->h[k];
    }
    theta = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) /
            (ci->h[0] * ci->h[0] + ci->h[1] * ci->h[1] + ci->h[2] * ci->h[2]);

    /* Split the interaction? */
    if (theta < const_theta_max * const_theta_max) {

      /* Are both ci and cj split? */
      if (ci->split && cj->split) {

        /* Split this task into tasks on its progeny. */
        for (j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL) {
            for (k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL)
                runner_dosub_grav(r, ci->progeny[j], cj->progeny[k], 0);
          }

      }

      /* Otherwise, make a pp task out of it. */
      else
        runner_dopair_grav(r, ci, cj);

    }

    /* Otherwise, mm interaction is fine. */
    else
      runner_dograv_mm(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(timer_dosub_grav);
}
#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
