/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "const.h"
#include "engine.h"
#include "error.h"
#include "scheduler.h"
#include "space.h"
#include "task.h"
#include "timers.h"

/* Include the right variant of the SPH interactions */
#ifdef LEGACY_GADGET2_SPH
#include "runner_iact_legacy.h"
#else
#include "runner_iact.h"
#endif
#include "runner_iact_grav.h"

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/* The counters. */
int runner_counter[runner_counter_count];

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
 * @brief Intermediate task between density and force
 *
 * @param r The runner thread.
 * @param c The cell.
 */

void runner_doghost(struct runner *r, struct cell *c) {

  struct part *p, *parts = c->parts;
  struct cell *finger;
  int i, k, redo, count = c->count;
  int *pid;
  float h, ih, ih2, ih4, h_corr, rho, wcount, rho_dh, wcount_dh, u, fc;
  float normDiv_v, normCurl_v;
#ifndef LEGACY_GADGET2_SPH
  float alpha_dot, tau, S;
#endif
  float dt_step = r->e->dt_step;
  TIMER_TIC

  /* Recurse? */
  if (c->split) {
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_doghost(r, c->progeny[k]);
    return;
  }

  /* Init the IDs that have to be updated. */
  if ((pid = (int *)alloca(sizeof(int) * count)) == NULL)
    error("Call to alloca failed.");
  for (k = 0; k < count; k++) pid[k] = k;

  /* While there are particles that need to be updated... */
  while (count > 0) {

    /* Reset the redo-count. */
    redo = 0;

    /* Loop over the parts in this cell. */
    __builtin_prefetch(&parts[pid[0]], 0, 1);
    __builtin_prefetch(&parts[pid[0]].rho_dh, 0, 1);
    __builtin_prefetch(&parts[pid[1]], 0, 1);
    __builtin_prefetch(&parts[pid[1]].rho_dh, 0, 1);
    __builtin_prefetch(&parts[pid[2]], 0, 1);
    __builtin_prefetch(&parts[pid[2]].rho_dh, 0, 1);
    for (i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      __builtin_prefetch(&parts[pid[i + 3]], 0, 1);
      __builtin_prefetch(&parts[pid[i + 3]].rho_dh, 0, 1);
      p = &parts[pid[i]];

      /* Is this part within the timestep? */
      if (p->dt <= dt_step) {

        /* Some smoothing length multiples. */
        h = p->h;
        ih = 1.0f / h;
        ih2 = ih * ih;
        ih4 = ih2 * ih2;

        /* Final operation on the density. */
        p->rho = rho = ih * ih2 * (p->rho + p->mass * kernel_root);
        p->rho_dh = rho_dh = (p->rho_dh - 3.0f * p->mass * kernel_root) * ih4;
        wcount = (p->density.wcount + kernel_root) *
                 (4.0f / 3.0 * M_PI * kernel_gamma3);
        wcount_dh =
            p->density.wcount_dh * ih * (4.0f / 3.0 * M_PI * kernel_gamma3);

        /* If no derivative, double the smoothing length. */
        if (wcount_dh == 0.0f) h_corr = p->h;

        /* Otherwise, compute the smoothing length update (Newton step). */
        else {
          h_corr = (kernel_nwneigh - wcount) / wcount_dh;

          /* Truncate to the range [ -p->h/2 , p->h ]. */
          h_corr = fminf(h_corr, h);
          h_corr = fmaxf(h_corr, -h / 2.f);
        }

        /* Apply the correction to p->h and to the compact part. */
        p->h += h_corr;

        /* Did we get the right number density? */
        if (wcount > kernel_nwneigh + const_delta_nwneigh ||
            wcount < kernel_nwneigh - const_delta_nwneigh) {
          // message( "particle %lli (h=%e,depth=%i) has bad wcount=%.3f." ,
          // p->id , p->h , c->depth , wcount ); fflush(stdout);
          // p->h += ( p->density.wcount + kernel_root - kernel_nwneigh ) /
          // p->density.wcount_dh;
          pid[redo] = pid[i];
          redo += 1;
          p->density.wcount = 0.0;
          p->density.wcount_dh = 0.0;
          p->rho = 0.0;
          p->rho_dh = 0.0;
          p->density.div_v = 0.0;
          for (k = 0; k < 3; k++) p->density.curl_v[k] = 0.0;
          continue;
        }

        /* Pre-compute some stuff for the balsara switch. */
        normDiv_v = fabs(p->density.div_v / rho * ih4);
        normCurl_v = sqrtf(p->density.curl_v[0] * p->density.curl_v[0] +
                           p->density.curl_v[1] * p->density.curl_v[1] +
                           p->density.curl_v[2] * p->density.curl_v[2]) /
                     rho * ih4;

        /* As of here, particle force variables will be set. Do _NOT_
           try to read any particle density variables! */

        /* Compute this particle's sound speed. */
        u = p->u;
        p->force.c = fc =
            sqrtf(const_hydro_gamma * (const_hydro_gamma - 1.0f) * u);

        /* Compute the P/Omega/rho2. */
        p->force.POrho2 =
            u * (const_hydro_gamma - 1.0f) / (rho + h * rho_dh / 3.0f);

        /* Balsara switch */
        p->force.balsara =
            normDiv_v / (normDiv_v + normCurl_v + 0.0001f * fc * ih);

#ifndef LEGACY_GADGET2_SPH
        /* Viscosity parameter decay time */
        tau = h / (2.f * const_viscosity_length * p->force.c);

        /* Viscosity source term */
        S = fmaxf(-normDiv_v, 0.f);

        /* Compute the particle's viscosity parameter time derivative */
        alpha_dot = (const_viscosity_alpha_min - p->alpha) / tau +
                    (const_viscosity_alpha_max - p->alpha) * S;

        /* Update particle's viscosity paramter */
        p->alpha += alpha_dot * p->dt;
#endif

        /* Reset the acceleration. */
        for (k = 0; k < 3; k++) p->a[k] = 0.0f;

        /* Reset the time derivatives. */
        p->force.u_dt = 0.0f;
        p->force.h_dt = 0.0f;
        p->force.v_sig = 0.0f;
      }
    }

    /* Re-set the counter for the next loop (potentially). */
    count = redo;
    if (count > 0) {

      // error( "Bad smoothing length, fixing this isn't implemented yet." );

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

#ifdef TIMER_VERBOSE
  message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
          c->depth, ((double)TIMER_TOC(timer_doghost)) / CPU_TPS * 1000);
  fflush(stdout);
#else
  TIMER_TOC(timer_doghost);
#endif
}

/**
 * @brief Compute the second kick of the given cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 */

void runner_dokick2(struct runner *r, struct cell *c) {

  int j, k, count = 0, nr_parts = c->count;
  float dt_min = FLT_MAX, dt_max = 0.0f;
  double ekin = 0.0, epot = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f}, ang[3] = {0.0f, 0.0f, 0.0f};
  float x[3], v_hdt[3], u_hdt, h, pdt, m;
  float dt_step = r->e->dt_step, dt = r->e->dt, hdt, idt;
  float dt_cfl, dt_h_change, dt_u_change, dt_new;
  float h_dt, u_dt;
  struct part *restrict p, *restrict parts = c->parts;
  struct xpart *restrict xp, *restrict xparts = c->xparts;

  TIMER_TIC

  /* Init idt to avoid compiler stupidity. */
  idt = (dt > 0) ? 1.0f / dt : 0.0f;
  hdt = dt / 2;

  /* Loop over the particles and kick them. */
  __builtin_prefetch(&parts[0], 0, 1);
  __builtin_prefetch(&parts[0].rho_dh, 0, 1);
  __builtin_prefetch(&xparts[0], 0, 1);
  __builtin_prefetch(&parts[1], 0, 1);
  __builtin_prefetch(&parts[1].rho_dh, 0, 1);
  __builtin_prefetch(&xparts[1], 0, 1);
  __builtin_prefetch(&parts[2], 0, 1);
  __builtin_prefetch(&parts[2].rho_dh, 0, 1);
  __builtin_prefetch(&xparts[2], 0, 1);
  for (k = 0; k < nr_parts; k++) {

    /* Get a handle on the part. */
    __builtin_prefetch(&parts[k + 3], 0, 1);
    __builtin_prefetch(&parts[k + 3].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[k + 3], 0, 1);
    p = &parts[k];
    xp = &xparts[k];

    /* Get local copies of particle data. */
    pdt = p->dt;
    m = p->mass;
    x[0] = p->x[0];
    x[1] = p->x[1];
    x[2] = p->x[2];
    v_hdt[0] = xp->v_hdt[0];
    v_hdt[1] = xp->v_hdt[1];
    v_hdt[2] = xp->v_hdt[2];
    u_hdt = xp->u_hdt;

    /* Update the particle's data (if active). */
    if (pdt <= dt_step) {

      /* Increase the number of particles updated. */
      count += 1;

      /* Scale the derivatives as they're freshly computed. */
      h = p->h;
      h_dt = p->force.h_dt *= h * 0.333333333f;
      xp->omega = 1.0f + h * p->rho_dh / p->rho * 0.3333333333f;

      /* Compute the new time step. */
      u_dt = p->force.u_dt;
      dt_cfl = const_cfl * h / p->force.v_sig;
      dt_h_change =
          (h_dt != 0.0f) ? fabsf(const_ln_max_h_change * h / h_dt) : FLT_MAX;
      dt_u_change =
          (u_dt != 0.0f) ? fabsf(const_max_u_change * p->u / u_dt) : FLT_MAX;
      dt_new = fminf(dt_cfl, fminf(dt_h_change, dt_u_change));
      if (pdt == 0.0f)
        p->dt = pdt = dt_new;
      else
        p->dt = pdt = fminf(dt_new, 2.0f * pdt);

      /* Update positions and energies at the full step. */
      p->v[0] = v_hdt[0] + hdt * p->a[0];
      p->v[1] = v_hdt[1] + hdt * p->a[1];
      p->v[2] = v_hdt[2] + hdt * p->a[2];
      p->u = u_hdt + hdt * u_dt;

      /* Set the new particle-specific time step. */
      if (dt > 0.0f) {
        float dt_curr = dt;
        j = (int)(pdt * idt);
        while (j > 1) {
          dt_curr *= 2.0f;
          j >>= 1;
        }
        xp->dt_curr = dt_curr;
      }
    }

    /* Get the smallest/largest dt. */
    dt_min = fminf(dt_min, pdt);
    dt_max = fmaxf(dt_max, pdt);

    /* Collect total energy. */
    ekin += 0.5 * m *
            (v_hdt[0] * v_hdt[0] + v_hdt[1] * v_hdt[1] + v_hdt[2] * v_hdt[2]);
    epot += m * u_hdt;

    /* Collect momentum */
    mom[0] += m * v_hdt[0];
    mom[1] += m * v_hdt[1];
    mom[2] += m * v_hdt[2];

    /* Collect angular momentum */
    ang[0] += m * (x[1] * v_hdt[2] - x[2] * v_hdt[1]);
    ang[1] += m * (x[2] * v_hdt[0] - x[0] * v_hdt[2]);
    ang[2] += m * (x[0] * v_hdt[1] - x[1] * v_hdt[0]);

    /* Collect entropic function */
    // lent += u * pow( p->rho, 1.f-const_gamma );
  }

#ifdef TIMER_VERBOSE
  message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
          c->depth, ((double)TIMER_TOC(timer_kick2)) / CPU_TPS * 1000);
  fflush(stdout);
#else
  TIMER_TOC(timer_kick2);
#endif

  /* Store the computed values in the cell. */
  c->dt_min = dt_min;
  c->dt_max = dt_max;
  c->updated = count;
  c->ekin = ekin;
  c->epot = epot;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];
}

/**
 * @brief Mapping function to set dt_min and dt_max, do the first
 * kick.
 */

void runner_dokick1(struct runner *r, struct cell *c) {

  int j, k;
  struct engine *e = r->e;
  float pdt, dt_step = e->dt_step, dt = e->dt, hdt = dt / 2;
  float dt_min, dt_max, h_max, dx, dx_max;
  float a[3], v[3], u, u_dt, h, h_dt, w, rho;
  double x[3], x_old[3];
  struct part *restrict p, *restrict parts = c->parts;
  struct xpart *restrict xp, *restrict xparts = c->xparts;

  /* No children? */
  if (!c->split) {

    /* Init the min/max counters. */
    dt_min = FLT_MAX;
    dt_max = 0.0f;
    h_max = 0.0f;
    dx_max = 0.0f;

    /* Loop over parts. */
    __builtin_prefetch(&parts[0], 0, 1);
    __builtin_prefetch(&parts[0].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[0], 0, 1);
    __builtin_prefetch(&parts[1], 0, 1);
    __builtin_prefetch(&parts[1].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[1], 0, 1);
    __builtin_prefetch(&parts[2], 0, 1);
    __builtin_prefetch(&parts[2].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[2], 0, 1);
    for (k = 0; k < c->count; k++) {

      /* Get a handle on the kth particle. */
      __builtin_prefetch(&parts[k + 3], 0, 1);
      __builtin_prefetch(&parts[k + 3].rho_dh, 0, 1);
      __builtin_prefetch(&xparts[k + 3], 0, 1);
      p = &parts[k];
      xp = &xparts[k];

      /* Load the data locally. */
      a[0] = p->a[0];
      a[1] = p->a[1];
      a[2] = p->a[2];
      v[0] = p->v[0];
      v[1] = p->v[1];
      v[2] = p->v[2];
      x[0] = p->x[0];
      x[1] = p->x[1];
      x[2] = p->x[2];
      x_old[0] = xp->x_old[0];
      x_old[1] = xp->x_old[1];
      x_old[2] = xp->x_old[2];
      h = p->h;
      u = p->u;
      h_dt = p->force.h_dt;
      u_dt = p->force.u_dt;
      pdt = p->dt;

      /* Store the min/max dt. */
      dt_min = fminf(dt_min, pdt);
      dt_max = fmaxf(dt_max, pdt);

      /* Update the half-step velocities from the current velocities. */
      xp->v_hdt[0] = v[0] + hdt * a[0];
      xp->v_hdt[1] = v[1] + hdt * a[1];
      xp->v_hdt[2] = v[2] + hdt * a[2];
      xp->u_hdt = u + hdt * u_dt;

      /* Move the particles with the velocities at the half-step. */
      p->x[0] = x[0] += dt * xp->v_hdt[0];
      p->x[1] = x[1] += dt * xp->v_hdt[1];
      p->x[2] = x[2] += dt * xp->v_hdt[2];
      dx = sqrtf((x[0] - x_old[0]) * (x[0] - x_old[0]) +
                 (x[1] - x_old[1]) * (x[1] - x_old[1]) +
                 (x[2] - x_old[2]) * (x[2] - x_old[2]));
      dx_max = fmaxf(dx_max, dx);

      /* Update positions and energies at the half-step. */
      p->v[0] = v[0] + dt * a[0];
      p->v[1] = v[1] + dt * a[1];
      p->v[2] = v[2] + dt * a[2];
      w = u_dt / u * dt;
      if (fabsf(w) < 0.01f)
        p->u = u *=
            1.0f +
            w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
      else
        p->u = u *= expf(w);
      w = h_dt / h * dt;
      if (fabsf(w) < 0.01f)
        p->h = h *=
            1.0f +
            w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
      else
        p->h = h *= expf(w);
      h_max = fmaxf(h_max, h);

      /* Integrate other values if this particle will not be updated. */
      /* Init fields for density calculation. */
      if (pdt > dt_step) {
        float w = -3.0f * h_dt / h * dt;
        if (fabsf(w) < 0.1f)
          rho = p->rho *=
              1.0f +
              w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w)));
        else
          rho = p->rho *= expf(w);
        p->force.POrho2 = u * (const_hydro_gamma - 1.0f) / (rho * xp->omega);
      } else {
        p->density.wcount = 0.0f;
        p->density.wcount_dh = 0.0f;
        p->rho = 0.0f;
        p->rho_dh = 0.0f;
        p->density.div_v = 0.0f;
        for (j = 0; j < 3; ++j) p->density.curl_v[j] = 0.0f;
      }
    }

  }

  /* Otherwise, agregate data from children. */
  else {

    /* Init with the first non-null child. */
    dt_min = FLT_MAX;
    dt_max = 0.0f;
    h_max = 0.0f;
    dx_max = 0.0f;

    /* Loop over the progeny. */
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        if (c->count < space_subsize) runner_dokick1(r, c->progeny[k]);
        dt_min = fminf(dt_min, c->progeny[k]->dt_min);
        dt_max = fmaxf(dt_max, c->progeny[k]->dt_max);
        h_max = fmaxf(h_max, c->progeny[k]->h_max);
        dx_max = fmaxf(dx_max, c->progeny[k]->dx_max);
      }
  }

  /* Store the values. */
  c->dt_min = dt_min;
  c->dt_max = dt_max;
  c->h_max = h_max;
  c->dx_max = dx_max;
}

/**
 * @brief Combined second and first kick for fixed dt.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer The timer
 */

void runner_dokick(struct runner *r, struct cell *c, int timer) {

  int k, count = 0, nr_parts = c->count, updated;
  float dt_min = FLT_MAX, dt_max = 0.0f;
  float h_max, dx, dx_max;
  double ekin = 0.0, epot = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f}, ang[3] = {0.0f, 0.0f, 0.0f};
  float x[3], x_old[3], v_hdt[3], a[3], u, u_hdt, h, pdt, m, w;
  float dt = r->e->dt, hdt = 0.5f * dt;
  float dt_cfl, dt_h_change, dt_u_change, dt_new;
  float h_dt, u_dt;
  struct part *restrict p, *restrict parts = c->parts;
  struct xpart *restrict xp, *restrict xparts = c->xparts;

  TIMER_TIC

  /* No children? */
  if (!c->split) {

    /* Init the min/max counters. */
    dt_min = FLT_MAX;
    dt_max = 0.0f;
    h_max = 0.0f;
    dx_max = 0.0f;

    /* Loop over the particles and kick them. */
    __builtin_prefetch(&parts[0], 0, 1);
    __builtin_prefetch(&parts[0].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[0], 0, 1);
    __builtin_prefetch(&parts[1], 0, 1);
    __builtin_prefetch(&parts[1].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[1], 0, 1);
    __builtin_prefetch(&parts[2], 0, 1);
    __builtin_prefetch(&parts[2].rho_dh, 0, 1);
    __builtin_prefetch(&xparts[2], 0, 1);
    for (k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      __builtin_prefetch(&parts[k + 3], 0, 1);
      __builtin_prefetch(&parts[k + 3].rho_dh, 0, 1);
      __builtin_prefetch(&xparts[k + 3], 0, 1);
      p = &parts[k];
      xp = &xparts[k];

      /* Get local copies of particle data. */
      pdt = p->dt;
      u_dt = p->force.u_dt;
      h = p->h;
      m = p->mass;
      x[0] = p->x[0];
      x[1] = p->x[1];
      x[2] = p->x[2];
      a[0] = p->a[0];
      a[1] = p->a[1];
      a[2] = p->a[2];
      x_old[0] = xp->x_old[0];
      x_old[1] = xp->x_old[1];
      x_old[2] = xp->x_old[2];
      v_hdt[0] = xp->v_hdt[0];
      v_hdt[1] = xp->v_hdt[1];
      v_hdt[2] = xp->v_hdt[2];
      u_hdt = xp->u_hdt;

      /* Scale the derivatives if they're freshly computed. */
      h_dt = p->force.h_dt *= h * 0.333333333f;
      count += 1;
      xp->omega = 1.0f + h * p->rho_dh / p->rho * 0.3333333333f;

      /* Update the particle's time step. */
      dt_cfl = const_cfl * h / p->force.v_sig;
      dt_h_change =
          (h_dt != 0.0f) ? fabsf(const_ln_max_h_change * h / h_dt) : FLT_MAX;
      dt_u_change =
          (u_dt != 0.0f) ? fabsf(const_max_u_change * p->u / u_dt) : FLT_MAX;
      dt_new = fminf(dt_cfl, fminf(dt_h_change, dt_u_change));
      if (pdt == 0.0f)
        p->dt = pdt = dt_new;
      else
        p->dt = pdt = fminf(dt_new, 2.0f * pdt);

      /* Get the smallest/largest dt. */
      dt_min = fminf(dt_min, pdt);
      dt_max = fmaxf(dt_max, pdt);

      /* Step and store the velocity and internal energy. */
      xp->v_hdt[0] = (v_hdt[0] += dt * a[0]);
      xp->v_hdt[1] = (v_hdt[1] += dt * a[1]);
      xp->v_hdt[2] = (v_hdt[2] += dt * a[2]);
      xp->u_hdt = (u_hdt += dt * u_dt);

      /* Move the particles with the velocitie at the half-step. */
      p->x[0] = x[0] += dt * v_hdt[0];
      p->x[1] = x[1] += dt * v_hdt[1];
      p->x[2] = x[2] += dt * v_hdt[2];
      dx = sqrtf((x[0] - x_old[0]) * (x[0] - x_old[0]) +
                 (x[1] - x_old[1]) * (x[1] - x_old[1]) +
                 (x[2] - x_old[2]) * (x[2] - x_old[2]));
      dx_max = fmaxf(dx_max, dx);

      /* Update positions and energies at the next full step. */
      p->v[0] = v_hdt[0] + hdt * a[0];
      p->v[1] = v_hdt[1] + hdt * a[1];
      p->v[2] = v_hdt[2] + hdt * a[2];
      w = u_dt / u_hdt * hdt;
      if (fabsf(w) < 0.01f)
        p->u = u =
            u_hdt *
            (1.0f +
             w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w))));
      else
        p->u = u = u_hdt * expf(w);
      w = h_dt / h * dt;
      if (fabsf(w) < 0.01f)
        p->h = h *=
            (1.0f +
             w * (1.0f + w * (0.5f + w * (1.0f / 6.0f + 1.0f / 24.0f * w))));
      else
        p->h = h *= expf(w);
      h_max = fmaxf(h_max, h);

      /* Collect momentum */
      mom[0] += m * v_hdt[0];
      mom[1] += m * v_hdt[1];
      mom[2] += m * v_hdt[2];

      /* Collect angular momentum */
      ang[0] += m * (x[1] * v_hdt[2] - x[2] * v_hdt[1]);
      ang[1] += m * (x[2] * v_hdt[0] - x[0] * v_hdt[2]);
      ang[2] += m * (x[0] * v_hdt[1] - x[1] * v_hdt[0]);

      /* Collect total energy. */
      ekin += 0.5 * m *
              (v_hdt[0] * v_hdt[0] + v_hdt[1] * v_hdt[1] + v_hdt[2] * v_hdt[2]);
      epot += m * u_hdt;

      /* Init fields for density calculation. */
      p->density.wcount = 0.0f;
      p->density.wcount_dh = 0.0f;
      p->rho = 0.0f;
      p->rho_dh = 0.0f;
      p->density.div_v = 0.0f;
      p->density.curl_v[0] = 0.0f;
      p->density.curl_v[1] = 0.0f;
      p->density.curl_v[2] = 0.0f;
    }

  }

  /* Otherwise, agregate data from children. */
  else {

    /* Init with the first non-null child. */
    dt_min = FLT_MAX;
    dt_max = 0.0f;
    h_max = 0.0f;
    dx_max = 0.0f;
    updated = 0;
    ekin = 0.0;
    epot = 0.0;
    mom[0] = 0.0f;
    mom[1] = 0.0f;
    mom[2] = 0.0f;
    ang[0] = 0.0f;
    ang[1] = 0.0f;
    ang[2] = 0.0f;

    /* Loop over the progeny. */
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_dokick(r, cp, 0);
        dt_min = fminf(dt_min, cp->dt_min);
        dt_max = fmaxf(dt_max, cp->dt_max);
        h_max = fmaxf(h_max, cp->h_max);
        dx_max = fmaxf(dx_max, cp->dx_max);
        updated += cp->count;
        ekin += cp->ekin;
        epot += cp->epot;
        mom[0] += cp->mom[0];
        mom[1] += cp->mom[1];
        mom[2] += cp->mom[2];
        ang[0] += cp->ang[0];
        ang[1] += cp->ang[1];
        ang[2] += cp->ang[2];
      }
  }

  /* Store the values. */
  c->dt_min = dt_min;
  c->dt_max = dt_max;
  c->h_max = h_max;
  c->dx_max = dx_max;
  c->updated = count;
  c->ekin = ekin;
  c->epot = epot;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];

  if (timer) {
#ifdef TIMER_VERBOSE
    message("runner %02i: %i parts at depth %i took %.3f ms.", r->id, c->count,
            c->depth, ((double)TIMER_TOC(timer_kick2)) / CPU_TPS * 1000);
    fflush(stdout);
#else
    TIMER_TOC(timer_kick2);
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
  struct cell *ci, *cj, *super;
  struct part *parts;
  int k, nr_parts;

  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    engine_barrier(e, r->id);

    /* Re-set the pointer to the previous super cell. */
    super = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, super);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL) break;
      }

      /* Get the cells. */
      ci = t->ci;
      cj = t->cj;
      t->rid = r->cpuid;

      /* Set super to the first cell that I own. */
      if (ci != NULL && ci->super != NULL && ci->super->owner == r->qid)
        super = ci->super;
      else if (cj != NULL && cj->super != NULL && cj->super->owner == r->qid)
        super = cj->super;
      else
        super = NULL;

      /* Prefetch? */
      if (runner_prefetch && t->type != task_type_kick1 &&
          t->type != task_type_kick2 && t->type != task_type_ghost) {
        for (int k = 0; k < ci->count; k++)
          __builtin_prefetch(&ci->parts[k], 1, 3);
        if (cj != NULL)
          for (int k = 0; k < cj->count; k++)
            __builtin_prefetch(&cj->parts[k], 1, 3);
      }

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
        case task_type_ghost:
          runner_doghost(r, ci);
          break;
        case task_type_kick1:
          runner_dokick1(r, ci);
          break;
        case task_type_kick2:
          if (e->policy & engine_policy_fixdt)
            runner_dokick(r, ci, 1);
          else
            runner_dokick2(r, ci);
          break;
        case task_type_send:
          break;
        case task_type_recv:
          parts = ci->parts;
          nr_parts = ci->count;
          for (k = 0; k < nr_parts; k++) parts[k].dt = FLT_MAX;
          ci->dt_min = ci->dt_max = FLT_MAX;
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
        default:
          error("Unknown task type.");
      }

      /* We're done with this task, see if we get a next one. */
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
