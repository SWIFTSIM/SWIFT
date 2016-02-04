/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "atomic.h"
#include "const.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "scheduler.h"
#include "space.h"
#include "task.h"
#include "timers.h"
#include "hydro.h"
#include "gravity.h"

/* Orientation of the cell pairs */
const float runner_shift[13 * 3] = {
    5.773502691896258e-01, 5.773502691896258e-01,  5.773502691896258e-01,
    7.071067811865475e-01, 7.071067811865475e-01,  0.0,
    5.773502691896258e-01, 5.773502691896258e-01,  -5.773502691896258e-01,
    7.071067811865475e-01, 0.0,                    7.071067811865475e-01,
    1.0,                   0.0,                    0.0,
    7.071067811865475e-01, 0.0,                    -7.071067811865475e-01,
    5.773502691896258e-01, -5.773502691896258e-01, 5.773502691896258e-01,
    7.071067811865475e-01, -7.071067811865475e-01, 0.0,
    5.773502691896258e-01, -5.773502691896258e-01, -5.773502691896258e-01,
    0.0,                   7.071067811865475e-01,  7.071067811865475e-01,
    0.0,                   1.0,                    0.0,
    0.0,                   7.071067811865475e-01,  -7.071067811865475e-01,
    0.0,                   0.0,                    1.0, };

/* Does the axis need flipping ? */
const char runner_flip[27] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Import the density loop functions. */
#define FUNCTION density
#include "runner_doiact.h"

/* Import the force loop functions. */
#undef FUNCTION
#define FUNCTION force
#include "runner_doiact.h"

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */

void runner_dosort_ascending(struct entry *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[10];
  int qpos, i, j, lo, hi, imin;
  struct entry temp;
  float pivot;

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    if (hi - lo < 15) {
      for (i = lo; i < hi; i++) {
        imin = i;
        for (j = i + 1; j <= hi; j++)
          if (sort[j].d < sort[imin].d) imin = j;
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      while (i <= j) {
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        if (i <= j) {
          if (i < j) {
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */

void runner_dosort(struct runner *r, struct cell *c, int flags, int clock) {

  struct entry *finger;
  struct entry *fingers[8];
  struct part *parts = c->parts;
  struct entry *sort;
  int j, k, count = c->count;
  int i, ind, off[8], inds[8], temp_i, missing;
  // float shift[3];
  float buff[8], px[3];

  TIMER_TIC

  /* Clean-up the flags, i.e. filter out what's already been sorted. */
  flags &= ~c->sorted;
  if (flags == 0) return;

  /* start by allocating the entry arrays. */
  if (c->sort == NULL || c->sortsize < count) {
    if (c->sort != NULL) free(c->sort);
    c->sortsize = count * 1.1;
    if ((c->sort = (struct entry *)malloc(sizeof(struct entry) *
                                          (c->sortsize + 1) * 13)) == NULL)
      error("Failed to allocate sort memory.");
  }
  sort = c->sort;

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    for (k = 0; k < 8; k++) {
      if (c->progeny[k] == NULL) continue;
      missing = flags & ~c->progeny[k]->sorted;
      if (missing) runner_dosort(r, c->progeny[k], missing, 0);
    }

    /* Loop over the 13 different sort arrays. */
    for (j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      for (off[0] = 0, k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      for (k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->count > 0) {
          fingers[k] = &c->progeny[k]->sort[j * (c->progeny[k]->count + 1)];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (i = 0; i < 7; i++)
        for (k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      finger = &sort[j * (count + 1)];
      for (ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      sort[j * (count + 1) + count].d = FLT_MAX;
      sort[j * (count + 1) + count].i = 0;

      /* Mark as sorted. */
      c->sorted |= (1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Fill the sort array. */
    for (k = 0; k < count; k++) {
      px[0] = parts[k].x[0];
      px[1] = parts[k].x[1];
      px[2] = parts[k].x[2];
      for (j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          sort[j * (count + 1) + k].i = k;
          sort[j * (count + 1) + k].d = px[0] * runner_shift[3 * j + 0] +
                                        px[1] * runner_shift[3 * j + 1] +
                                        px[2] * runner_shift[3 * j + 2];
        }
    }

    /* Add the sentinel and sort. */
    for (j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        sort[j * (count + 1) + count].d = FLT_MAX;
        sort[j * (count + 1) + count].i = 0;
        runner_dosort_ascending(&sort[j * (count + 1)], count);
        c->sorted |= (1 << j);
      }
  }

/* Verify the sorting. */
/* for ( j = 0 ; j < 13 ; j++ ) {
    if ( !( flags & (1 << j) ) )
        continue;
    finger = &sort[ j*(count + 1) ];
    for ( k = 1 ; k < count ; k++ ) {
        if ( finger[k].d < finger[k-1].d )
            error( "Sorting failed, ascending array." );
        if ( finger[k].i >= count )
            error( "Sorting failed, indices borked." );
        }
    } */

#ifdef TIMER_VERBOSE
  message(
      "runner %02i: %i parts at depth %i (flags = %i%i%i%i%i%i%i%i%i%i%i%i%i) "
      "took %.3f ms.",
      r->id, count, c->depth, (flags & 0x1000) >> 12, (flags & 0x800) >> 11,
      (flags & 0x400) >> 10, (flags & 0x200) >> 9, (flags & 0x100) >> 8,
      (flags & 0x80) >> 7, (flags & 0x40) >> 6, (flags & 0x20) >> 5,
      (flags & 0x10) >> 4, (flags & 0x8) >> 3, (flags & 0x4) >> 2,
      (flags & 0x2) >> 1, (flags & 0x1) >> 0,
      ((double)TIMER_TOC(timer_dosort)) / CPU_TPS * 1000);
  fflush(stdout);
#else
  if (clock) TIMER_TOC(timer_dosort);
#endif
}

void runner_dogsort(struct runner *r, struct cell *c, int flags, int clock) {

  struct entry *finger;
  struct entry *fingers[8];
  struct gpart *gparts = c->gparts;
  struct entry *gsort;
  int j, k, count = c->gcount;
  int i, ind, off[8], inds[8], temp_i, missing;
  // float shift[3];
  float buff[8], px[3];

  TIMER_TIC

  /* Clean-up the flags, i.e. filter out what's already been sorted. */
  flags &= ~c->gsorted;
  if (flags == 0) return;

  /* start by allocating the entry arrays. */
  if (c->gsort == NULL || c->gsortsize < count) {
    if (c->gsort != NULL) free(c->gsort);
    c->gsortsize = count * 1.1;
    if ((c->gsort = (struct entry *)malloc(sizeof(struct entry) *
                                           (c->gsortsize + 1) * 13)) == NULL)
      error("Failed to allocate sort memory.");
  }
  gsort = c->gsort;

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    for (k = 0; k < 8; k++) {
      if (c->progeny[k] == NULL) continue;
      missing = flags & ~c->progeny[k]->gsorted;
      if (missing) runner_dogsort(r, c->progeny[k], missing, 0);
    }

    /* Loop over the 13 different sort arrays. */
    for (j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      for (off[0] = 0, k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->gcount;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      for (k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->gcount > 0) {
          fingers[k] = &c->progeny[k]->gsort[j * (c->progeny[k]->gcount + 1)];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (i = 0; i < 7; i++)
        for (k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      finger = &gsort[j * (count + 1)];
      for (ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      gsort[j * (count + 1) + count].d = FLT_MAX;
      gsort[j * (count + 1) + count].i = 0;

      /* Mark as sorted. */
      c->gsorted |= (1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Fill the sort array. */
    for (k = 0; k < count; k++) {
      px[0] = gparts[k].x[0];
      px[1] = gparts[k].x[1];
      px[2] = gparts[k].x[2];
      for (j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          gsort[j * (count + 1) + k].i = k;
          gsort[j * (count + 1) + k].d = px[0] * runner_shift[3 * j + 0] +
                                         px[1] * runner_shift[3 * j + 1] +
                                         px[2] * runner_shift[3 * j + 2];
        }
    }

    /* Add the sentinel and sort. */
    for (j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        gsort[j * (count + 1) + count].d = FLT_MAX;
        gsort[j * (count + 1) + count].i = 0;
        runner_dosort_ascending(&gsort[j * (count + 1)], count);
        c->gsorted |= (1 << j);
      }
  }

/* Verify the sorting. */
/* for ( j = 0 ; j < 13 ; j++ ) {
    if ( !( flags & (1 << j) ) )
        continue;
    finger = &c->gsort[ j*(count + 1) ];
    for ( k = 1 ; k < count ; k++ ) {
        if ( finger[k].d < finger[k-1].d )
            error( "Sorting failed, ascending array." );
        if ( finger[k].i < 0 || finger[k].i >= count )
            error( "Sorting failed, indices borked." );
        }
    } */

#ifdef TIMER_VERBOSE
  message(
      "runner %02i: %i parts at depth %i (flags = %i%i%i%i%i%i%i%i%i%i%i%i%i) "
      "took %.3f ms.",
      r->id, count, c->depth, (flags & 0x1000) >> 12, (flags & 0x800) >> 11,
      (flags & 0x400) >> 10, (flags & 0x200) >> 9, (flags & 0x100) >> 8,
      (flags & 0x80) >> 7, (flags & 0x40) >> 6, (flags & 0x20) >> 5,
      (flags & 0x10) >> 4, (flags & 0x8) >> 3, (flags & 0x4) >> 2,
      (flags & 0x2) >> 1, (flags & 0x1) >> 0,
      ((double)TIMER_TOC(timer_dosort)) / CPU_TPS * 1000);
  fflush(stdout);
#else
  if (clock) TIMER_TOC(timer_dosort);
#endif
}

/**
 * @brief Initialize the particles before the density calculation
 *
 * @param r The runner thread.
 * @param c The cell.
 */

void runner_doinit(struct runner *r, struct cell *c, int timer) {

  struct part *p, *parts = c->parts;
  const int count = c->count;
  const float t_end = r->e->time;

  TIMER_TIC;
  
  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_doinit(r, c->progeny[k], 0);
    return;
  }
  else {
  
    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {
      
      /* Get a direct pointer on the part. */
      p = &parts[i];

      if(p->id == 1000) message("init 1000!");
      if(p->id == 515050) message("init 515050!");
      
      if (p->t_end <= t_end) {
	
	/* Get ready for a density calculation */
	hydro_init_part(p);
      }
    }
  }

  if (timer) {
#ifdef TIMER_VERBOSE
    message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
            c->depth, ((double)TIMER_TOC(timer_init)) / CPU_TPS * 1000);
    fflush(stdout);
#else
    TIMER_TOC(timer_init);
#endif
  }

  
}

/**
 * @brief Intermediate task between density and force
 *
 * @param r The runner thread.
 * @param c The cell.
 */

void runner_doghost(struct runner *r, struct cell *c) {

  struct part *p, *parts = c->parts;
  struct xpart *xp, *xparts = c->xparts;
  struct cell *finger;
  int redo, count = c->count;
  int *pid;
  float h_corr;
  float t_end = r->e->time;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_doghost(r, c->progeny[k]);
    return;
  }

  /* Init the IDs that have to be updated. */
  if ((pid = (int *)alloca(sizeof(int) * count)) == NULL)
    error("Call to alloca failed.");
  for (int k = 0; k < count; k++) pid[k] = k;

  /* While there are particles that need to be updated... */
  for (int num_reruns = 0; count > 0 && num_reruns < const_smoothing_max_iter;
       num_reruns++) {

    /* Reset the redo-count. */
    redo = 0;

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      p = &parts[pid[i]];
      xp = &xparts[pid[i]];

      if(p->id == 1000) message("ghost 1000");
      if(p->id == 515050) message("ghost 515050");
      
      /* Is this part within the timestep? */
      if (p->t_end <= t_end) {

	/* Finish the density calculation */
	hydro_end_density(p, t_end);

        /* If no derivative, double the smoothing length. */
        if (p->density.wcount_dh == 0.0f) h_corr = p->h;

        /* Otherwise, compute the smoothing length update (Newton step). */
        else {
          h_corr = (kernel_nwneigh - p->density.wcount) / p->density.wcount_dh;

          /* Truncate to the range [ -p->h/2 , p->h ]. */
          h_corr = fminf(h_corr, p->h);
          h_corr = fmaxf(h_corr, -p->h * 0.5f);
        }

        /* Did we get the right number density? */
        if (p->density.wcount > kernel_nwneigh + const_delta_nwneigh ||
           p->density. wcount < kernel_nwneigh - const_delta_nwneigh) {

          /* Ok, correct then */
          p->h += h_corr;

          /* Flag for another round of fun */
          pid[redo] = pid[i];
          redo += 1;

	  /* Re-initialise everything */
	  hydro_init_part(p);

	  /* Off we go ! */
          continue;
        }

        /* We now have a particle whose smoothing length has converged */
	//if(p->id == 1000)
	//  printParticle(parts, 1000, count);
	
        /* As of here, particle force variables will be set. Do _NOT_
           try to read any particle density variables! */

	/* Compute variables required for the force loop */
	hydro_prepare_force(p, xp);
	
	/* Prepare the particle for the force loop over neighbours */
	hydro_reset_acceleration(p);

      }
    }

    /* We now need to treat the particles whose smoothing length had not
     * converged again */

    /* Re-set the counter for the next loop (potentially). */
    count = redo;
    if (count > 0) {

      /* Climb up the cell hierarchy. */
      for (finger = c; finger != NULL; finger = finger->parent) {

        /* Run through this cell's density interactions. */
        for (struct link *l = finger->density; l != NULL; l = l->next) {

          /* Self-interaction? */
          if (l->t->type == task_type_self)
            runner_doself_subset_density(r, finger, parts, pid, count);

          /* Otherwise, pair interaction? */
          else if (l->t->type == task_type_pair) {

            /* Left or right? */
            if (l->t->ci == finger)
              runner_dopair_subset_density(r, finger, parts, pid, count,
                                           l->t->cj);
            else
              runner_dopair_subset_density(r, finger, parts, pid, count,
                                           l->t->ci);

          }

          /* Otherwise, sub interaction? */
          else if (l->t->type == task_type_sub) {

            /* Left or right? */
            if (l->t->ci == finger)
              runner_dosub_subset_density(r, finger, parts, pid, count,
                                          l->t->cj, -1, 1);
            else
              runner_dosub_subset_density(r, finger, parts, pid, count,
                                          l->t->ci, -1, 1);
          }
        }
      }
    }
  }

  if (count)
    message("Smoothing length failed to converge on %i particles.", count);

#ifdef TIMER_VERBOSE
  message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
          c->depth, ((double)TIMER_TOC(timer_doghost)) / CPU_TPS * 1000);
  fflush(stdout);
#else
  TIMER_TOC(timer_doghost);
#endif
}

/**
 * @brief Drift particles forward in time
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_dodrift(struct runner *r, struct cell *c, int timer) {

  const int nr_parts = c->count;
  const float dt = r->e->time - r->e->timeOld;
  struct part *restrict p, *restrict parts = c->parts;
  struct xpart *restrict xp, *restrict xparts = c->xparts;
  float w;
  float dx_max = 0.f, h_max = 0.f;
  
  TIMER_TIC

  /* No children? */
  if (!c->split) {

    /* Loop over all the particles in the cell */
    for (int k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      p = &parts[k];
      xp = &xparts[k];

      /* Get local copies of particle data. */
      const float h = p->h;
      const float ih = 1.0f / h;

      /* Drift... */
      p->x[0] += xp->v_full[0] * dt;
      p->x[1] += xp->v_full[1] * dt;
      p->x[2] += xp->v_full[2] * dt;

      /* Predict velocities */
      p->v[0] += p->a[0] * dt;
      p->v[1] += p->a[1] * dt;
      p->v[2] += p->a[2] * dt;

      /* Predict smoothing length */
      w = p->force.h_dt * ih * dt;
      if (fabsf(w) < 0.01f) /* 1st order expansion of exp(w) */
	p->h *=
	  1.0f +
	  w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
      else
	p->h *= expf(w);
      
      /* Predict density */
      w = -3.0f * p->force.h_dt * ih * dt;
      if (fabsf(w) < 0.1f)
	p->rho *=
	  1.0f +
	  w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
      else
	p->rho *= expf(w);
      
      /* Predict the values of the extra fields */
      hydro_predict_extra(p, xp, dt);

      /* Compute motion since last cell construction */
      const float dx = sqrtf((p->x[0] - xp->x_old[0]) * (p->x[0] - xp->x_old[0]) +
			     (p->x[1] - xp->x_old[1]) * (p->x[1] - xp->x_old[1]) +
			     (p->x[2] - xp->x_old[2]) * (p->x[2] - xp->x_old[2]));
      dx_max = fmaxf(dx_max, dx);

      /* Maximal smoothing length */
      h_max = fmaxf(p->h, h_max);
    }
  }

  /* Otherwise, aggregate data from children. */
  else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_dodrift(r, cp, 0);

	dx_max = fmaxf(dx_max, cp->dx_max);
	h_max = fmaxf(h_max, cp->h_max);
      }
  }

  /* Store the values */
  c->h_max = h_max;
  c->dx_max = dx_max;
  
  if (timer) {
#ifdef TIMER_VERBOSE
    message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
            c->depth, ((double)TIMER_TOC(timer_drift)) / CPU_TPS * 1000);
    fflush(stdout);
#else
    TIMER_TOC(timer_drift);
#endif
  }
}

/**
 * @brief Combined second and first kick for fixed dt.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer The timer
 */

void runner_dokick(struct runner *r, struct cell *c, int timer) {

  const float dt_max_timeline = r->e->timeEnd - r->e->timeBegin;
  const float global_dt_min = r->e->dt_min, global_dt_max = r->e->dt_max;
  const float t_current = r->e->time;
  const int count = c->count;
  const int is_fixdt = (r->e->policy & engine_policy_fixdt) == engine_policy_fixdt;

  float new_dt;
  float t_start, t_end, dt;
  float dt_timeline;

  int updated = 0;
  float t_end_min = FLT_MAX, t_end_max = 0.f;
  double ekin = 0.0, epot = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f}, ang[3] = {0.0f, 0.0f, 0.0f};
  float m, x[3], v_full[3];
  struct part *restrict p, *restrict parts = c->parts;
  struct xpart *restrict xp, *restrict xparts = c->xparts;

  TIMER_TIC

  /* No children? */
  if (!c->split) {

    /* Loop over the particles and kick the active ones. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      p = &parts[k];
      xp = &xparts[k];

      m = p->mass;
      x[0] = p->x[0];
      x[1] = p->x[1];
      x[2] = p->x[2];

      /* If particle needs to be kicked */
      if ( is_fixdt || p->t_end <= t_current ) {

        /* First, finish the force loop */
	hydro_end_force(p);
	  
	if( is_fixdt ) {

	  /* Now we have a time step, proceed with the kick */
	  new_dt = global_dt_max;
	  
	} else {
	
	  /* Compute the next timestep */
	  const float new_dt_hydro = hydro_compute_timestep(p, xp);
	  const float new_dt_grav = gravity_compute_timestep(p, xp);
	  
	  new_dt = fminf(new_dt_hydro, new_dt_grav);

	  /* Recover the current timestep */
	  const float current_dt = p->t_end - p->t_begin;
	  
	  /* Limit timestep increase */
	  if (current_dt > 0.0f) new_dt = fminf(new_dt, 2.0f * current_dt);
	
	  /* Limit timestep within the allowed range */
	  new_dt = fminf(new_dt, global_dt_max);
	  new_dt = fmaxf(new_dt, global_dt_min);
	
	  /* Put this timestep on the time line */
	  dt_timeline = dt_max_timeline;
	  while (new_dt < dt_timeline) dt_timeline /= 2.;
	  
	  /* Now we have a time step, proceed with the kick */
	  new_dt = dt_timeline;
	}

        /* Compute the time step for this kick */
        t_start = 0.5f * (p->t_begin + p->t_end);
        t_end = p->t_end + 0.5f * new_dt;
        dt = t_end - t_start;

        /* Move particle forward in time */
        p->t_begin = p->t_end;
        p->t_end = p->t_begin + new_dt;

        /* Kick particles in momentum space */
        xp->v_full[0] += p->a[0] * dt;
        xp->v_full[1] += p->a[1] * dt;
        xp->v_full[2] += p->a[2] * dt;

        p->v[0] = xp->v_full[0] - 0.5f * new_dt * p->a[0];
        p->v[1] = xp->v_full[1] - 0.5f * new_dt * p->a[1];
        p->v[2] = xp->v_full[2] - 0.5f * new_dt * p->a[2];
      }

      /* Now collect quantities for statistics */

      v_full[0] = xp->v_full[0];
      v_full[1] = xp->v_full[1];
      v_full[2] = xp->v_full[2];

      /* Collect momentum */
      mom[0] += m * v_full[0];
      mom[1] += m * v_full[1];
      mom[2] += m * v_full[2];

      /* Collect angular momentum */
      ang[0] += m * (x[1] * v_full[2] - x[2] * v_full[1]);
      ang[1] += m * (x[2] * v_full[0] - x[0] * v_full[2]);
      ang[2] += m * (x[0] * v_full[1] - x[1] * v_full[0]);

      /* Collect total energy. */
      ekin += 0.5 * m * (v_full[0] * v_full[0] + v_full[1] * v_full[1] +
                         v_full[2] * v_full[2]);
      epot += m * xp->u_hdt;

      /* Minimal time for next end of time-step */
      t_end_min = fminf(p->t_end, t_end_min);
      t_end_max = fmaxf(p->t_end, t_end_max);

      /* Number of updated particles */
      updated++;
    }

  }

  /* Otherwise, aggregate data from children. */
  else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_dokick(r, cp, 0);

	updated += cp->updated;
        ekin += cp->ekin;
        epot += cp->epot;
        mom[0] += cp->mom[0];
        mom[1] += cp->mom[1];
        mom[2] += cp->mom[2];
        ang[0] += cp->ang[0];
        ang[1] += cp->ang[1];
        ang[2] += cp->ang[2];
        t_end_min = fminf(cp->t_end_min, t_end_min);
        t_end_max = fmaxf(cp->t_end_max, t_end_max);
      }
  }

  /* Store the values. */
  c->updated = updated;
  c->ekin = ekin;
  c->epot = epot;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];
  c->t_end_min = t_end_min;
  c->t_end_max = t_end_max;

  if (timer) {
#ifdef TIMER_VERBOSE
    message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
            c->depth, ((double)TIMER_TOC(timer_kick)) / CPU_TPS * 1000);
    fflush(stdout);
#else
    TIMER_TOC(timer_kick);
#endif
  }
}

/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */

void *runner_main(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  struct task *t = NULL;
  struct cell *ci, *cj;
  struct part *parts;
  int k, nr_parts;

  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    engine_barrier(e, r->id);

    /* Re-set the pointer to the previous task, as there is none. */
    struct task* prev = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, prev);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL) break;
      }

      /* Get the cells. */
      ci = t->ci;
      cj = t->cj;
      t->rid = r->cpuid;

      /* Different types of tasks... */
      switch (t->type) {
        case task_type_self:
          if (t->subtype == task_subtype_density)
            runner_doself1_density(r, ci);
          else if (t->subtype == task_subtype_force)
	    runner_doself2_force(r, ci);
          else
            error("Unknown task subtype.");
          break;
        case task_type_pair:
          if (t->subtype == task_subtype_density)
            runner_dopair1_density(r, ci, cj);
          else if (t->subtype == task_subtype_force)
	    runner_dopair2_force(r, ci, cj);
          else
            error("Unknown task subtype.");
          break;
        case task_type_sort:
          runner_dosort(r, ci, t->flags, 1);
          break;
        case task_type_sub:
          if (t->subtype == task_subtype_density)
            runner_dosub1_density(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_force)
            runner_dosub2_force(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_grav)
            runner_dosub_grav(r, ci, cj, 1);
          else
            error("Unknown task subtype.");
          break;
        case task_type_init:
          runner_doinit(r, ci, 1);
          break;
        case task_type_ghost:
          runner_doghost(r, ci);
          break;
        case task_type_drift:
          runner_dodrift(r, ci, 1);
          break;
        case task_type_kick:
          runner_dokick(r, ci, 1);
          break;
        case task_type_send:
          break;
        case task_type_recv:
          parts = ci->parts;
          nr_parts = ci->count;
          ci->t_end_min = ci->t_end_max = FLT_MAX;
          for (k = 0; k < nr_parts; k++) parts[k].t_end = FLT_MAX;
          break;
        case task_type_grav_pp:
          if (t->cj == NULL)
            runner_doself_grav(r, t->ci);
          else
            runner_dopair_grav(r, t->ci, t->cj);
          break;
        case task_type_grav_mm:
          runner_dograv_mm(r, t->ci, t->cj);
          break;
        case task_type_grav_up:
          runner_dograv_up(r, t->ci);
          break;
        case task_type_grav_down:
          runner_dograv_down(r, t->ci);
          break;
        case task_type_psort:
          space_do_parts_sort();
          break;
        case task_type_split_cell:
          space_split(e->s, t->ci);
          break;
        case task_type_rewait:
	  scheduler_do_rewait((struct task *)t->ci, (struct task *)t->cj, t->flags, t->rank);
          break;
        default:
          error("Unknown task type.");
      }

      /* We're done with this task, see if we get a next one. */
      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
