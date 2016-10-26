/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "atomic.h"
#include "const.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "lock.h"
#include "minmax.h"
#include "runner.h"
#include "threadpool.h"
#include "tools.h"

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize = space_subsize_default;
int space_maxsize = space_maxsize_default;
int space_maxcount = space_maxcount_default;

/* Map shift vector to sortlist. */
const int sortlistID[27] = {
    /* ( -1 , -1 , -1 ) */ 0,
    /* ( -1 , -1 ,  0 ) */ 1,
    /* ( -1 , -1 ,  1 ) */ 2,
    /* ( -1 ,  0 , -1 ) */ 3,
    /* ( -1 ,  0 ,  0 ) */ 4,
    /* ( -1 ,  0 ,  1 ) */ 5,
    /* ( -1 ,  1 , -1 ) */ 6,
    /* ( -1 ,  1 ,  0 ) */ 7,
    /* ( -1 ,  1 ,  1 ) */ 8,
    /* (  0 , -1 , -1 ) */ 9,
    /* (  0 , -1 ,  0 ) */ 10,
    /* (  0 , -1 ,  1 ) */ 11,
    /* (  0 ,  0 , -1 ) */ 12,
    /* (  0 ,  0 ,  0 ) */ 0,
    /* (  0 ,  0 ,  1 ) */ 12,
    /* (  0 ,  1 , -1 ) */ 11,
    /* (  0 ,  1 ,  0 ) */ 10,
    /* (  0 ,  1 ,  1 ) */ 9,
    /* (  1 , -1 , -1 ) */ 8,
    /* (  1 , -1 ,  0 ) */ 7,
    /* (  1 , -1 ,  1 ) */ 6,
    /* (  1 ,  0 , -1 ) */ 5,
    /* (  1 ,  0 ,  0 ) */ 4,
    /* (  1 ,  0 ,  1 ) */ 3,
    /* (  1 ,  1 , -1 ) */ 2,
    /* (  1 ,  1 ,  0 ) */ 1,
    /* (  1 ,  1 ,  1 ) */ 0};

/**
 * @brief Interval stack necessary for parallel particle sorting.
 */
struct qstack {
  volatile ptrdiff_t i, j;
  volatile int min, max;
  volatile int ready;
};

/**
 * @brief Parallel particle-sorting stack
 */
struct parallel_sort {
  struct part *parts;
  struct gpart *gparts;
  struct xpart *xparts;
  int *ind;
  struct qstack *stack;
  unsigned int stack_size;
  volatile unsigned int first, last, waiting;
};

/**
 * @brief Information required to compute the particle cell indices.
 */
struct index_data {
  struct space *s;
  struct cell *cells;
  int *ind;
};

/**
 * @brief Get the shift-id of the given pair of cells, swapping them
 *      if need be.
 *
 * @param s The space
 * @param ci Pointer to first #cell.
 * @param cj Pointer second #cell.
 * @param shift Vector from ci to cj.
 *
 * @return The shift ID and set shift, may or may not swap ci and cj.
 */
int space_getsid(struct space *s, struct cell **ci, struct cell **cj,
                 double *shift) {

  /* Get the relative distance between the pairs, wrapping. */
  const int periodic = s->periodic;
  double dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = (*cj)->loc[k] - (*ci)->loc[k];
    if (periodic && dx[k] < -s->dim[k] / 2)
      shift[k] = s->dim[k];
    else if (periodic && dx[k] > s->dim[k] / 2)
      shift[k] = -s->dim[k];
    else
      shift[k] = 0.0;
    dx[k] += shift[k];
  }

  /* Get the sorting index. */
  int sid = 0;
  for (int k = 0; k < 3; k++)
    sid = 3 * sid + ((dx[k] < 0.0) ? 0 : ((dx[k] > 0.0) ? 2 : 1));

  /* Switch the cells around? */
  if (runner_flip[sid]) {
    struct cell *temp = *ci;
    *ci = *cj;
    *cj = temp;
    for (int k = 0; k < 3; k++) shift[k] = -shift[k];
  }
  sid = sortlistID[sid];

  /* Return the sort ID. */
  return sid;
}

/**
 * @brief Recursively dismantle a cell tree.
 *
 * @param s The #space.
 * @param c The #cell to recycle.
 */
void space_rebuild_recycle(struct space *s, struct cell *c) {

  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        space_rebuild_recycle(s, c->progeny[k]);
        space_recycle(s, c->progeny[k]);
        c->progeny[k] = NULL;
      }
}

/**
 * @brief Re-build the top-level cell grid.
 *
 * @param s The #space.
 * @param cell_max Maximum cell edge length.
 * @param verbose Print messages to stdout or not.
 */
void space_regrid(struct space *s, double cell_max, int verbose) {

  const size_t nr_parts = s->nr_parts;
  const ticks tic = getticks();
  const int ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the cells and get the current h_max. */
  // tic = getticks();
  float h_max = s->cell_min / kernel_gamma / space_stretch;
  if (nr_parts > 0) {
    if (s->cells_top != NULL) {
      for (int k = 0; k < s->nr_cells; k++) {
        if (s->cells_top[k].nodeID == engine_rank &&
            s->cells_top[k].h_max > h_max) {
          h_max = s->cells_top[k].h_max;
        }
      }
    } else {
      for (size_t k = 0; k < nr_parts; k++) {
        if (s->parts[k].h > h_max) h_max = s->parts[k].h;
      }
    }
  }

/* If we are running in parallel, make sure everybody agrees on
   how large the largest cell should be. */
#ifdef WITH_MPI
  {
    float buff;
    if (MPI_Allreduce(&h_max, &buff, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    h_max = buff;
  }
#endif
  if (verbose) message("h_max is %.3e (cell_max=%.3e).", h_max, cell_max);

  /* Get the new putative cell dimensions. */
  const int cdim[3] = {
      floor(s->dim[0] / fmax(h_max * kernel_gamma * space_stretch, cell_max)),
      floor(s->dim[1] / fmax(h_max * kernel_gamma * space_stretch, cell_max)),
      floor(s->dim[2] / fmax(h_max * kernel_gamma * space_stretch, cell_max))};

  /* Check if we have enough cells for periodicity. */
  if (s->periodic && (cdim[0] < 3 || cdim[1] < 3 || cdim[2] < 3))
    error(
        "Must have at least 3 cells in each spatial dimension when periodicity "
        "is switched on.\nThis error is often caused by any of the "
        "followings:\n"
        " - too few particles to generate a sensible grid,\n"
        " - the initial value of 'SPH:max_smoothing_length' is too large,\n"
        " - the (minimal) time-step is too large leading to particles with "
        "predicted smoothing lengths too large for the box size,\n"
        " - particle with velocities so large that they move by more than two "
        "box sizes per time-step.\n");

  /* Check if we have enough cells for gravity. */
  if (s->gravity && (cdim[0] < 8 || cdim[1] < 8 || cdim[2] < 8))
    error(
        "Must have at least 8 cells in each spatial dimension when gravity "
        "is switched on.");

/* In MPI-Land, changing the top-level cell size requires that the
 * global partition is recomputed and the particles redistributed.
 * Be prepared to do that. */
#ifdef WITH_MPI
  double oldwidth[3];
  double oldcdim[3];
  int *oldnodeIDs = NULL;
  if (cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] || cdim[2] < s->cdim[2]) {

    /* Capture state of current space. */
    oldcdim[0] = s->cdim[0];
    oldcdim[1] = s->cdim[1];
    oldcdim[2] = s->cdim[2];
    oldwidth[0] = s->width[0];
    oldwidth[1] = s->width[1];
    oldwidth[2] = s->width[2];

    if ((oldnodeIDs = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate temporary nodeIDs.");

    int cid = 0;
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {
          cid = cell_getid(oldcdim, i, j, k);
          oldnodeIDs[cid] = s->cells_top[cid].nodeID;
        }
      }
    }
  }

#endif

  /* Do we need to re-build the upper-level cells? */
  // tic = getticks();
  if (s->cells_top == NULL || cdim[0] < s->cdim[0] || cdim[1] < s->cdim[1] ||
      cdim[2] < s->cdim[2]) {

    /* Free the old cells, if they were allocated. */
    if (s->cells_top != NULL) {
      for (int k = 0; k < s->nr_cells; k++) {
        space_rebuild_recycle(s, &s->cells_top[k]);
        if (s->cells_top[k].sort != NULL) free(s->cells_top[k].sort);
      }
      free(s->cells_top);
      s->maxdepth = 0;
    }

    /* Set the new cell dimensions only if smaller. */
    for (int k = 0; k < 3; k++) {
      s->cdim[k] = cdim[k];
      s->width[k] = s->dim[k] / cdim[k];
      s->iwidth[k] = 1.0 / s->width[k];
    }
    const float dmin = min(s->width[0], min(s->width[1], s->width[2]));

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells = cdim[0] * cdim[1] * cdim[2];
    if (posix_memalign((void *)&s->cells_top, cell_align,
                       s->nr_cells * sizeof(struct cell)) != 0)
      error("Failed to allocate cells.");
    bzero(s->cells_top, s->nr_cells * sizeof(struct cell));
    for (int k = 0; k < s->nr_cells; k++)
      if (lock_init(&s->cells_top[k].lock) != 0)
        error("Failed to init spinlock.");

    /* Set the cell location and sizes. */
    for (int i = 0; i < cdim[0]; i++)
      for (int j = 0; j < cdim[1]; j++)
        for (int k = 0; k < cdim[2]; k++) {
          struct cell *restrict c = &s->cells_top[cell_getid(cdim, i, j, k)];
          c->loc[0] = i * s->width[0];
          c->loc[1] = j * s->width[1];
          c->loc[2] = k * s->width[2];
          c->width[0] = s->width[0];
          c->width[1] = s->width[1];
          c->width[2] = s->width[2];
          c->dmin = dmin;
          c->depth = 0;
          c->count = 0;
          c->gcount = 0;
          c->super = c;
          c->ti_old = ti_current;
          lock_init(&c->lock);
        }

    /* Be verbose about the change. */
    if (verbose)
      message("set cell dimensions to [ %i %i %i ].", cdim[0], cdim[1],
              cdim[2]);
    fflush(stdout);

#ifdef WITH_MPI
    if (oldnodeIDs != NULL) {
      /* We have changed the top-level cell dimension, so need to redistribute
       * cells around the nodes. We repartition using the old space node
       * positions as a grid to resample. */
      if (s->e->nodeID == 0)
        message(
            "basic cell dimensions have increased - recalculating the "
            "global partition.");

      if (!partition_space_to_space(oldwidth, oldcdim, oldnodeIDs, s)) {

        /* Failed, try another technique that requires no settings. */
        message("Failed to get a new partition, trying less optimal method");
        struct partition initial_partition;
#ifdef HAVE_METIS
        initial_partition.type = INITPART_METIS_NOWEIGHT;
#else
        initial_partition.type = INITPART_VECTORIZE;
#endif
        partition_initial_partition(&initial_partition, s->e->nodeID,
                                    s->e->nr_nodes, s);
      }

      /* Re-distribute the particles to their new nodes. */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);

      /* Finished with these. */
      free(oldnodeIDs);
    }
#endif

    // message( "rebuilding upper-level cells took %.3f %s." ,
    // clocks_from_ticks(double)(getticks() - tic), clocks_getunit());

  } /* re-build upper-level cells? */

  else { /* Otherwise, just clean up the cells. */

    /* Free the old cells, if they were allocated. */
    for (int k = 0; k < s->nr_cells; k++) {
      space_rebuild_recycle(s, &s->cells_top[k]);
      s->cells_top[k].sorts = NULL;
      s->cells_top[k].nr_tasks = 0;
      s->cells_top[k].density = NULL;
      s->cells_top[k].gradient = NULL;
      s->cells_top[k].force = NULL;
      s->cells_top[k].grav = NULL;
      s->cells_top[k].dx_max = 0.0f;
      s->cells_top[k].sorted = 0;
      s->cells_top[k].count = 0;
      s->cells_top[k].gcount = 0;
      s->cells_top[k].init = NULL;
      s->cells_top[k].extra_ghost = NULL;
      s->cells_top[k].ghost = NULL;
      s->cells_top[k].kick = NULL;
      s->cells_top[k].cooling = NULL;
      s->cells_top[k].sourceterms = NULL;
      s->cells_top[k].super = &s->cells_top[k];
#if WITH_MPI
      s->cells_top[k].recv_xv = NULL;
      s->cells_top[k].recv_rho = NULL;
      s->cells_top[k].recv_gradient = NULL;
      s->cells_top[k].recv_ti = NULL;

      s->cells_top[k].send_xv = NULL;
      s->cells_top[k].send_rho = NULL;
      s->cells_top[k].send_gradient = NULL;
      s->cells_top[k].send_ti = NULL;
#endif
    }
    s->maxdepth = 0;
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Re-build the cells as well as the tasks.
 *
 * @param s The #space in which to update the cells.
 * @param cell_max Maximal cell size.
 * @param verbose Print messages to stdout or not
 *
 */
void space_rebuild(struct space *s, double cell_max, int verbose) {

  const ticks tic = getticks();

  /* Be verbose about this. */
  // message("re)building space..."); fflush(stdout);

  /* Re-grid if necessary, or just re-set the cell data. */
  space_regrid(s, cell_max, verbose);

  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  struct cell *restrict cells_top = s->cells_top;
  const int ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the particles and get their cell index. */
  const size_t ind_size = s->size_parts;
  int *ind;
  if ((ind = (int *)malloc(sizeof(int) * ind_size)) == NULL)
    error("Failed to allocate temporary particle indices.");
  if (ind_size > 0) space_parts_get_cell_index(s, ind, cells_top, verbose);
  for (size_t i = 0; i < s->nr_parts; ++i) cells_top[ind[i]].count++;

  /* Run through the gravity particles and get their cell index. */
  const size_t gind_size = s->size_gparts;
  int *gind;
  if ((gind = (int *)malloc(sizeof(int) * gind_size)) == NULL)
    error("Failed to allocate temporary g-particle indices.");
  if (gind_size > 0) space_gparts_get_cell_index(s, gind, cells_top, verbose);
  for (size_t i = 0; i < s->nr_gparts; ++i) cells_top[gind[i]].gcount++;

#ifdef WITH_MPI

  /* Move non-local parts to the end of the list. */
  const int local_nodeID = s->e->nodeID;
  for (size_t k = 0; k < nr_parts;) {
    if (cells_top[ind[k]].nodeID != local_nodeID) {
      cells_top[ind[k]].count -= 1;
      nr_parts -= 1;
      const struct part tp = s->parts[k];
      s->parts[k] = s->parts[nr_parts];
      s->parts[nr_parts] = tp;
      if (s->parts[k].gpart != NULL) {
        s->parts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->parts[nr_parts].gpart != NULL) {
        s->parts[nr_parts].gpart->id_or_neg_offset = -nr_parts;
      }
      const struct xpart txp = s->xparts[k];
      s->xparts[k] = s->xparts[nr_parts];
      s->xparts[nr_parts] = txp;
      const int t = ind[k];
      ind[k] = ind[nr_parts];
      ind[nr_parts] = t;
    } else {
      /* Increment when not exchanging otherwise we need to retest "k".*/
      k++;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all parts are in the correct places. */
  for (size_t k = 0; k < nr_parts; k++) {
    if (cells_top[ind[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local parts to send list");
    }
  }
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    if (cells_top[ind[k]].nodeID == local_nodeID) {
      error("Failed to remove local parts from send list");
    }
  }
#endif

  /* Move non-local gparts to the end of the list. */
  for (size_t k = 0; k < nr_gparts;) {
    if (cells_top[gind[k]].nodeID != local_nodeID) {
      cells_top[gind[k]].gcount -= 1;
      nr_gparts -= 1;
      const struct gpart tp = s->gparts[k];
      s->gparts[k] = s->gparts[nr_gparts];
      s->gparts[nr_gparts] = tp;
      if (s->gparts[k].id_or_neg_offset <= 0) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }
      if (s->gparts[nr_gparts].id_or_neg_offset <= 0) {
        s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      }
      const int t = gind[k];
      gind[k] = gind[nr_gparts];
      gind[nr_gparts] = t;
    } else {
      /* Increment when not exchanging otherwise we need to retest "k".*/
      k++;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all gparts are in the correct place (untested). */
  for (size_t k = 0; k < nr_gparts; k++) {
    if (cells_top[gind[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local gparts to send list");
    }
  }
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    if (cells_top[gind[k]].nodeID == local_nodeID) {
      error("Failed to remove local gparts from send list");
    }
  }
#endif

  /* Exchange the strays, note that this potentially re-allocates
     the parts arrays. */
  size_t nr_parts_exchanged = s->nr_parts - nr_parts;
  size_t nr_gparts_exchanged = s->nr_gparts - nr_gparts;
  engine_exchange_strays(s->e, nr_parts, &ind[nr_parts], &nr_parts_exchanged,
                         nr_gparts, &gind[nr_gparts], &nr_gparts_exchanged);

  /* Set the new particle counts. */
  s->nr_parts = nr_parts + nr_parts_exchanged;
  s->nr_gparts = nr_gparts + nr_gparts_exchanged;

  /* Re-allocate the index array if needed.. */
  if (s->nr_parts > ind_size) {
    int *ind_new;
    if ((ind_new = (int *)malloc(sizeof(int) * s->nr_parts)) == NULL)
      error("Failed to allocate temporary particle indices.");
    memcpy(ind_new, ind, sizeof(int) * nr_parts);
    free(ind);
    ind = ind_new;
  }

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih[3] = {s->iwidth[0], s->iwidth[1], s->iwidth[2]};

  /* Assign each particle to its cell. */
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    const struct part *const p = &s->parts[k];
    ind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cells_top[ind[k]].count += 1;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[ind[k]].nodeID != local_nodeID)
      error("Received part that does not belong to me (nodeID=%i).",
            cells_top[ind[k]].nodeID);
#endif
  }
  nr_parts = s->nr_parts;

#endif /* WITH_MPI */

  /* Sort the parts according to their cells. */
  space_parts_sort(s, ind, nr_parts, 0, s->nr_cells - 1, verbose);

  /* Re-link the gparts. */
  if (nr_parts > 0 && nr_gparts > 0) part_relink_gparts(s->parts, nr_parts, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify space_sort_struct. */
  for (size_t k = 1; k < nr_parts; k++) {
    if (ind[k - 1] > ind[k]) {
      error("Sort failed!");
    } else if (ind[k] != cell_getid(s->cdim, s->parts[k].x[0] * s->iwidth[0],
                                    s->parts[k].x[1] * s->iwidth[1],
                                    s->parts[k].x[2] * s->iwidth[2])) {
      error("Incorrect indices!");
    }
  }
#endif

  /* We no longer need the indices as of here. */
  free(ind);

#ifdef WITH_MPI

  /* Re-allocate the index array if needed.. */
  if (s->nr_gparts > gind_size) {
    int *gind_new;
    if ((gind_new = (int *)malloc(sizeof(int) * s->nr_gparts)) == NULL)
      error("Failed to allocate temporary g-particle indices.");
    memcpy(gind_new, gind, sizeof(int) * nr_gparts);
    free(gind);
    gind = gind_new;
  }

  /* Assign each particle to its cell. */
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    const struct gpart *const p = &s->gparts[k];
    gind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cells_top[gind[k]].gcount += 1;

#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[ind[k]].nodeID != s->e->nodeID)
      error("Received part that does not belong to me (nodeID=%i).",
            cells_top[ind[k]].nodeID);
#endif
  }
  nr_gparts = s->nr_gparts;

#endif

  /* Sort the gparts according to their cells. */
  space_gparts_sort(s, gind, nr_gparts, 0, s->nr_cells - 1, verbose);

  /* Re-link the parts. */
  if (nr_parts > 0 && nr_gparts > 0)
    part_relink_parts(s->gparts, nr_gparts, s->parts);

  /* We no longer need the indices as of here. */
  free(gind);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  for (size_t k = 0; k < nr_gparts; ++k) {

    if (s->gparts[k].id_or_neg_offset < 0) {

      const struct part *part = &s->parts[-s->gparts[k].id_or_neg_offset];

      if (part->gpart != &s->gparts[k]) error("Linking problem !");

      if (s->gparts[k].x[0] != part->x[0] || s->gparts[k].x[1] != part->x[1] ||
          s->gparts[k].x[2] != part->x[2])
        error("Linked particles are not at the same position !");
    }
  }
  for (size_t k = 0; k < nr_parts; ++k) {

    if (s->parts[k].gpart != NULL &&
        s->parts[k].gpart->id_or_neg_offset != -(ptrdiff_t)k) {
      error("Linking problem !");
    }
  }
#endif

  /* Hook the cells up to the parts. */
  // tic = getticks();
  struct part *finger = s->parts;
  struct xpart *xfinger = s->xparts;
  struct gpart *gfinger = s->gparts;
  for (int k = 0; k < s->nr_cells; k++) {
    struct cell *restrict c = &cells_top[k];
    c->ti_old = ti_current;
    c->parts = finger;
    c->xparts = xfinger;
    c->gparts = gfinger;
    finger = &finger[c->count];
    xfinger = &xfinger[c->count];
    gfinger = &gfinger[c->gcount];
  }
  // message( "hooking up cells took %.3f %s." ,
  // clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* At this point, we have the upper-level cells, old or new. Now make
     sure that the parts in each cell are ok. */
  space_split(s, cells_top, s->nr_cells, verbose);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Split particles between cells of a hierarchy
 *
 * This is done in parallel using threads in the #threadpool.
 *
 * @param s The #space.
 * @param cells The cell hierarchy.
 * @param nr_cells The number of cells.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, struct cell *cells, int nr_cells,
                 int verbose) {

  const ticks tic = getticks();

  threadpool_map(&s->e->threadpool, space_split_mapper, cells, nr_cells,
                 sizeof(struct cell), 1, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Runs through the top-level cells and checks whether tasks associated
 * with them can be split. If not, try to sanitize the cells.
 *
 * @param s The #space to act upon.
 */
void space_sanitize(struct space *s) {

  for (int k = 0; k < s->nr_cells; k++) {

    struct cell *c = &s->cells_top[k];
    const double min_width = c->dmin;

    /* Do we have a problem ? */
    if (c->h_max * kernel_gamma * space_stretch > min_width * 0.5 &&
        c->count > space_maxcount) {

      /* Ok, clean-up the mess */
      cell_sanitize(c);
    }
  }
}

/**
 * @brief #threadpool mapper function to compute the particle cell indices.
 *
 * @param map_data Pointer towards the particles.
 * @param nr_parts The number of particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_parts_get_cell_index_mapper(void *map_data, int nr_parts,
                                       void *extra_data) {

  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(parts - s->parts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    struct part *restrict p = &parts[k];

    const double old_pos_x = p->x[0];
    const double old_pos_y = p->x[1];
    const double old_pos_z = p->x[2];

    /* Put it back into the simulation volume */
    const double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    const double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    const double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);
    ind[k] = index;

    /* Update the position */
    p->x[0] = pos_x;
    p->x[1] = pos_y;
    p->x[2] = pos_z;
  }
}

/**
 * @brief #threadpool mapper function to compute the g-particle cell indices.
 *
 * @param map_data Pointer towards the g-particles.
 * @param nr_gparts The number of g-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_gparts_get_cell_index_mapper(void *map_data, int nr_gparts,
                                        void *extra_data) {

  /* Unpack the data */
  struct gpart *restrict gparts = (struct gpart *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(gparts - s->gparts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  for (int k = 0; k < nr_gparts; k++) {

    /* Get the particle */
    struct gpart *restrict gp = &gparts[k];

    const double old_pos_x = gp->x[0];
    const double old_pos_y = gp->x[1];
    const double old_pos_z = gp->x[2];

    /* Put it back into the simulation volume */
    const double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    const double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    const double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);
    ind[k] = index;

    /* Update the position */
    gp->x[0] = pos_x;
    gp->x[1] = pos_y;
    gp->x[2] = pos_z;
  }
}

/**
 * @brief Computes the cell index of all the particles and update the cell
 * count.
 *
 * @param s The #space.
 * @param ind The array of indices to fill.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_parts_get_cell_index(struct space *s, int *ind, struct cell *cells,
                                int verbose) {

  const ticks tic = getticks();

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = ind;

  threadpool_map(&s->e->threadpool, space_parts_get_cell_index_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), 1000, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the g-particles and update the cell
 * gcount.
 *
 * @param s The #space.
 * @param gind The array of indices to fill.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_gparts_get_cell_index(struct space *s, int *gind, struct cell *cells,
                                 int verbose) {

  const ticks tic = getticks();

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = gind;

  threadpool_map(&s->e->threadpool, space_gparts_get_cell_index_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart), 1000, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Sort the particles and condensed particles according to the given
 * indices.
 *
 * @param s The #space.
 * @param ind The indices with respect to which the parts are sorted.
 * @param N The number of parts
 * @param min Lowest index.
 * @param max highest index.
 * @param verbose Are we talkative ?
 */
void space_parts_sort(struct space *s, int *ind, size_t N, int min, int max,
                      int verbose) {

  const ticks tic = getticks();

  /* Populate a parallel_sort structure with the input data */
  struct parallel_sort sort_struct;
  sort_struct.parts = s->parts;
  sort_struct.xparts = s->xparts;
  sort_struct.ind = ind;
  sort_struct.stack_size = 2 * (max - min + 1) + 10 + s->e->nr_threads;
  if ((sort_struct.stack =
           malloc(sizeof(struct qstack) * sort_struct.stack_size)) == NULL)
    error("Failed to allocate sorting stack.");
  for (unsigned int i = 0; i < sort_struct.stack_size; i++)
    sort_struct.stack[i].ready = 0;

  /* Add the first interval. */
  sort_struct.stack[0].i = 0;
  sort_struct.stack[0].j = N - 1;
  sort_struct.stack[0].min = min;
  sort_struct.stack[0].max = max;
  sort_struct.stack[0].ready = 1;
  sort_struct.first = 0;
  sort_struct.last = 1;
  sort_struct.waiting = 1;

  /* Launch the sorting tasks with a stride of zero such that the same
     map data is passed to each thread. */
  threadpool_map(&s->e->threadpool, space_parts_sort_mapper, &sort_struct,
                 s->e->threadpool.num_threads, 0, 1, NULL);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify space_sort_struct. */
  for (size_t i = 1; i < N; i++)
    if (ind[i - 1] > ind[i])
      error("Sorting failed (ind[%zu]=%i,ind[%zu]=%i), min=%i, max=%i.", i - 1,
            ind[i - 1], i, ind[i], min, max);
  message("Sorting succeeded.");
#endif

  /* Clean up. */
  free(sort_struct.stack);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_parts_sort_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  /* Unpack the mapping data. */
  struct parallel_sort *sort_struct = (struct parallel_sort *)map_data;

  /* Pointers to the sorting data. */
  int *ind = sort_struct->ind;
  struct part *parts = sort_struct->parts;
  struct xpart *xparts = sort_struct->xparts;

  /* Main loop. */
  while (sort_struct->waiting) {

    /* Grab an interval off the queue. */
    int qid = atomic_inc(&sort_struct->first) % sort_struct->stack_size;

    /* Wait for the entry to be ready, or for the sorting do be done. */
    while (!sort_struct->stack[qid].ready)
      if (!sort_struct->waiting) return;

    /* Get the stack entry. */
    ptrdiff_t i = sort_struct->stack[qid].i;
    ptrdiff_t j = sort_struct->stack[qid].j;
    int min = sort_struct->stack[qid].min;
    int max = sort_struct->stack[qid].max;
    sort_struct->stack[qid].ready = 0;

    /* Loop over sub-intervals. */
    while (1) {

      /* Bring beer. */
      const int pivot = (min + max) / 2;
      /* message("Working on interval [%i,%i] with min=%i, max=%i, pivot=%i.",
              i, j, min, max, pivot); */

      /* One pass of QuickSort's partitioning. */
      ptrdiff_t ii = i;
      ptrdiff_t jj = j;
      while (ii < jj) {
        while (ii <= j && ind[ii] <= pivot) ii++;
        while (jj >= i && ind[jj] > pivot) jj--;
        if (ii < jj) {
          size_t temp_i = ind[ii];
          ind[ii] = ind[jj];
          ind[jj] = temp_i;
          struct part temp_p = parts[ii];
          parts[ii] = parts[jj];
          parts[jj] = temp_p;
          struct xpart temp_xp = xparts[ii];
          xparts[ii] = xparts[jj];
          xparts[jj] = temp_xp;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify space_sort_struct. */
      for (int k = i; k <= jj; k++)
        if (ind[k] > pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.",
                  k, ind[k], pivot, i, j);
          error("Partition failed (<=pivot).");
        }
      for (int k = jj + 1; k <= j; k++)
        if (ind[k] <= pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.",
                  k, ind[k], pivot, i, j);
          error("Partition failed (>pivot).");
        }
#endif

      /* Split-off largest interval. */
      if (jj - i > j - jj + 1) {

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          qid = atomic_inc(&sort_struct->last) % sort_struct->stack_size;
          while (sort_struct->stack[qid].ready)
            ;
          sort_struct->stack[qid].i = i;
          sort_struct->stack[qid].j = jj;
          sort_struct->stack[qid].min = min;
          sort_struct->stack[qid].max = pivot;
          if (atomic_inc(&sort_struct->waiting) >= sort_struct->stack_size)
            error("Qstack overflow.");
          sort_struct->stack[qid].ready = 1;
        }

        /* Recurse on the right? */
        if (jj + 1 < j && pivot + 1 < max) {
          i = jj + 1;
          min = pivot + 1;
        } else
          break;

      } else {

        /* Recurse on the right? */
        if (pivot + 1 < max) {
          qid = atomic_inc(&sort_struct->last) % sort_struct->stack_size;
          while (sort_struct->stack[qid].ready)
            ;
          sort_struct->stack[qid].i = jj + 1;
          sort_struct->stack[qid].j = j;
          sort_struct->stack[qid].min = pivot + 1;
          sort_struct->stack[qid].max = max;
          if (atomic_inc(&sort_struct->waiting) >= sort_struct->stack_size)
            error("Qstack overflow.");
          sort_struct->stack[qid].ready = 1;
        }

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          j = jj;
          max = pivot;
        } else
          break;
      }

    } /* loop over sub-intervals. */

    atomic_dec(&sort_struct->waiting);

  } /* main loop. */
}

/**
 * @brief Sort the g-particles according to the given indices.
 *
 * @param s The #space.
 * @param ind The indices with respect to which the gparts are sorted.
 * @param N The number of gparts
 * @param min Lowest index.
 * @param max highest index.
 * @param verbose Are we talkative ?
 */
void space_gparts_sort(struct space *s, int *ind, size_t N, int min, int max,
                       int verbose) {

  const ticks tic = getticks();

  /*Populate a global parallel_sort structure with the input data */
  struct parallel_sort sort_struct;
  sort_struct.gparts = s->gparts;
  sort_struct.ind = ind;
  sort_struct.stack_size = 2 * (max - min + 1) + 10 + s->e->nr_threads;
  if ((sort_struct.stack =
           malloc(sizeof(struct qstack) * sort_struct.stack_size)) == NULL)
    error("Failed to allocate sorting stack.");
  for (unsigned int i = 0; i < sort_struct.stack_size; i++)
    sort_struct.stack[i].ready = 0;

  /* Add the first interval. */
  sort_struct.stack[0].i = 0;
  sort_struct.stack[0].j = N - 1;
  sort_struct.stack[0].min = min;
  sort_struct.stack[0].max = max;
  sort_struct.stack[0].ready = 1;
  sort_struct.first = 0;
  sort_struct.last = 1;
  sort_struct.waiting = 1;

  /* Launch the sorting tasks with a stride of zero such that the same
     map data is passed to each thread. */
  threadpool_map(&s->e->threadpool, space_gparts_sort_mapper, &sort_struct,
                 s->e->threadpool.num_threads, 0, 1, NULL);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify space_sort_struct. */
  for (size_t i = 1; i < N; i++)
    if (ind[i - 1] > ind[i])
      error("Sorting failed (ind[%zu]=%i,ind[%zu]=%i), min=%i, max=%i.", i - 1,
            ind[i - 1], i, ind[i], min, max);
  message("Sorting succeeded.");
#endif

  /* Clean up. */
  free(sort_struct.stack);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_gparts_sort_mapper(void *map_data, int num_elements,
                              void *extra_data) {

  /* Unpack the mapping data. */
  struct parallel_sort *sort_struct = (struct parallel_sort *)map_data;

  /* Pointers to the sorting data. */
  int *ind = sort_struct->ind;
  struct gpart *gparts = sort_struct->gparts;

  /* Main loop. */
  while (sort_struct->waiting) {

    /* Grab an interval off the queue. */
    int qid = atomic_inc(&sort_struct->first) % sort_struct->stack_size;

    /* Wait for the entry to be ready, or for the sorting do be done. */
    while (!sort_struct->stack[qid].ready)
      if (!sort_struct->waiting) return;

    /* Get the stack entry. */
    ptrdiff_t i = sort_struct->stack[qid].i;
    ptrdiff_t j = sort_struct->stack[qid].j;
    int min = sort_struct->stack[qid].min;
    int max = sort_struct->stack[qid].max;
    sort_struct->stack[qid].ready = 0;

    /* Loop over sub-intervals. */
    while (1) {

      /* Bring beer. */
      const int pivot = (min + max) / 2;
      /* message("Working on interval [%i,%i] with min=%i, max=%i, pivot=%i.",
              i, j, min, max, pivot); */

      /* One pass of QuickSort's partitioning. */
      ptrdiff_t ii = i;
      ptrdiff_t jj = j;
      while (ii < jj) {
        while (ii <= j && ind[ii] <= pivot) ii++;
        while (jj >= i && ind[jj] > pivot) jj--;
        if (ii < jj) {
          size_t temp_i = ind[ii];
          ind[ii] = ind[jj];
          ind[jj] = temp_i;
          struct gpart temp_p = gparts[ii];
          gparts[ii] = gparts[jj];
          gparts[jj] = temp_p;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify space_sort_struct. */
      for (int k = i; k <= jj; k++)
        if (ind[k] > pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.",
                  k, ind[k], pivot, i, j);
          error("Partition failed (<=pivot).");
        }
      for (int k = jj + 1; k <= j; k++)
        if (ind[k] <= pivot) {
          message("sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.",
                  k, ind[k], pivot, i, j);
          error("Partition failed (>pivot).");
        }
#endif

      /* Split-off largest interval. */
      if (jj - i > j - jj + 1) {

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          qid = atomic_inc(&sort_struct->last) % sort_struct->stack_size;
          while (sort_struct->stack[qid].ready)
            ;
          sort_struct->stack[qid].i = i;
          sort_struct->stack[qid].j = jj;
          sort_struct->stack[qid].min = min;
          sort_struct->stack[qid].max = pivot;
          if (atomic_inc(&sort_struct->waiting) >= sort_struct->stack_size)
            error("Qstack overflow.");
          sort_struct->stack[qid].ready = 1;
        }

        /* Recurse on the right? */
        if (jj + 1 < j && pivot + 1 < max) {
          i = jj + 1;
          min = pivot + 1;
        } else
          break;

      } else {

        /* Recurse on the right? */
        if (pivot + 1 < max) {
          qid = atomic_inc(&sort_struct->last) % sort_struct->stack_size;
          while (sort_struct->stack[qid].ready)
            ;
          sort_struct->stack[qid].i = jj + 1;
          sort_struct->stack[qid].j = j;
          sort_struct->stack[qid].min = pivot + 1;
          sort_struct->stack[qid].max = max;
          if (atomic_inc(&sort_struct->waiting) >= sort_struct->stack_size)
            error("Qstack overflow.");
          sort_struct->stack[qid].ready = 1;
        }

        /* Recurse on the left? */
        if (jj > i && pivot > min) {
          j = jj;
          max = pivot;
        } else
          break;
      }

    } /* loop over sub-intervals. */

    atomic_dec(&sort_struct->waiting);

  } /* main loop. */
}

/**
 * @brief Mapping function to free the sorted indices buffers.
 */
void space_map_clearsort(struct cell *c, void *data) {

  if (c->sort != NULL) {
    free(c->sort);
    c->sort = NULL;
  }
}

/**
 * @brief Map a function to all particles in a cell recursively.
 *
 * @param c The #cell we are working in.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */
static void rec_map_parts(struct cell *c,
                          void (*fun)(struct part *p, struct cell *c,
                                      void *data),
                          void *data) {
  /* No progeny? */
  if (!c->split)
    for (int k = 0; k < c->count; k++) fun(&c->parts[k], c, data);

  /* Otherwise, recurse. */
  else
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) rec_map_parts(c->progeny[k], fun, data);
}

/**
 * @brief Map a function to all particles in a space.
 *
 * @param s The #space we are working in.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */
void space_map_parts(struct space *s,
                     void (*fun)(struct part *p, struct cell *c, void *data),
                     void *data) {

  /* Call the recursive function on all higher-level cells. */
  for (int cid = 0; cid < s->nr_cells; cid++)
    rec_map_parts(&s->cells_top[cid], fun, data);
}

/**
 * @brief Map a function to all particles in a cell recursively.
 *
 * @param c The #cell we are working in.
 * @param fun Function pointer to apply on the cells.
 */
static void rec_map_parts_xparts(struct cell *c,
                                 void (*fun)(struct part *p, struct xpart *xp,
                                             struct cell *c)) {

  /* No progeny? */
  if (!c->split)
    for (int k = 0; k < c->count; k++) fun(&c->parts[k], &c->xparts[k], c);

  /* Otherwise, recurse. */
  else
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) rec_map_parts_xparts(c->progeny[k], fun);
}

/**
 * @brief Map a function to all particles (#part and #xpart) in a space.
 *
 * @param s The #space we are working in.
 * @param fun Function pointer to apply on the particles in the cells.
 */
void space_map_parts_xparts(struct space *s,
                            void (*fun)(struct part *p, struct xpart *xp,
                                        struct cell *c)) {

  /* Call the recursive function on all higher-level cells. */
  for (int cid = 0; cid < s->nr_cells; cid++)
    rec_map_parts_xparts(&s->cells_top[cid], fun);
}

/**
 * @brief Map a function to all particles in a cell recursively.
 *
 * @param c The #cell we are working in.
 * @param full Map to all cells, including cells with sub-cells.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */
static void rec_map_cells_post(struct cell *c, int full,
                               void (*fun)(struct cell *c, void *data),
                               void *data) {
  /* Recurse. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        rec_map_cells_post(c->progeny[k], full, fun, data);

  /* No progeny? */
  if (full || !c->split) fun(c, data);
}

/**
 * @brief Map a function to all particles in a aspace.
 *
 * @param s The #space we are working in.
 * @param full Map to all cells, including cells with sub-cells.
 * @param fun Function pointer to apply on the cells.
 * @param data Data passed to the function fun.
 */
void space_map_cells_post(struct space *s, int full,
                          void (*fun)(struct cell *c, void *data), void *data) {

  /* Call the recursive function on all higher-level cells. */
  for (int cid = 0; cid < s->nr_cells; cid++)
    rec_map_cells_post(&s->cells_top[cid], full, fun, data);
}

static void rec_map_cells_pre(struct cell *c, int full,
                              void (*fun)(struct cell *c, void *data),
                              void *data) {

  /* No progeny? */
  if (full || !c->split) fun(c, data);

  /* Recurse. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        rec_map_cells_pre(c->progeny[k], full, fun, data);
}

/**
 * @brief Calls function fun on the cells in the space s
 *
 * @param s The #space
 * @param full If true calls the function on all cells and not just on leaves
 * @param fun The function to call.
 * @param data Additional data passed to fun() when called
 */
void space_map_cells_pre(struct space *s, int full,
                         void (*fun)(struct cell *c, void *data), void *data) {

  /* Call the recursive function on all higher-level cells. */
  for (int cid = 0; cid < s->nr_cells; cid++)
    rec_map_cells_pre(&s->cells_top[cid], full, fun, data);
}

/**
 * @brief #threadpool mapper function to split cells if they contain
 *        too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void space_split_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *restrict cells_top = (struct cell *)map_data;
  struct engine *e = s->e;

  for (int ind = 0; ind < num_cells; ind++) {

    struct cell *c = &cells_top[ind];

    const int count = c->count;
    const int gcount = c->gcount;
    const int depth = c->depth;
    int maxdepth = 0;
    float h_max = 0.0f;
    int ti_end_min = max_nr_timesteps, ti_end_max = 0;
    struct cell *temp;
    struct part *parts = c->parts;
    struct gpart *gparts = c->gparts;
    struct xpart *xparts = c->xparts;

    /* Check the depth. */
    while (depth > (maxdepth = s->maxdepth)) {
      atomic_cas(&s->maxdepth, maxdepth, depth);
    }

    /* If the depth is too large, we have a problem and should stop. */
    if (maxdepth > space_cell_maxdepth) {
      error("Exceeded maximum depth (%d) when splitting cells, aborting",
            space_cell_maxdepth);
    }

    /* Split or let it be? */
    if (count > space_splitsize || gcount > space_splitsize) {

      /* No longer just a leaf. */
      c->split = 1;

      /* Create the cell's progeny. */
      for (int k = 0; k < 8; k++) {
        temp = space_getcell(s);
        temp->count = 0;
        temp->gcount = 0;
        temp->ti_old = e->ti_current;
        temp->loc[0] = c->loc[0];
        temp->loc[1] = c->loc[1];
        temp->loc[2] = c->loc[2];
        temp->width[0] = c->width[0] / 2;
        temp->width[1] = c->width[1] / 2;
        temp->width[2] = c->width[2] / 2;
        temp->dmin = c->dmin / 2;
        if (k & 4) temp->loc[0] += temp->width[0];
        if (k & 2) temp->loc[1] += temp->width[1];
        if (k & 1) temp->loc[2] += temp->width[2];
        temp->depth = c->depth + 1;
        temp->split = 0;
        temp->h_max = 0.0;
        temp->dx_max = 0.f;
        temp->nodeID = c->nodeID;
        temp->parent = c;
        temp->super = NULL;
        c->progeny[k] = temp;
      }

      /* Split the cell data. */
      cell_split(c, c->parts - s->parts);

      /* Remove any progeny with zero parts. */
      for (int k = 0; k < 8; k++)
        if (c->progeny[k]->count == 0 && c->progeny[k]->gcount == 0) {
          space_recycle(s, c->progeny[k]);
          c->progeny[k] = NULL;
        } else {
          space_split_mapper(c->progeny[k], 1, s);
          h_max = max(h_max, c->progeny[k]->h_max);
          ti_end_min = min(ti_end_min, c->progeny[k]->ti_end_min);
          ti_end_max = max(ti_end_max, c->progeny[k]->ti_end_max);
          if (c->progeny[k]->maxdepth > maxdepth)
            maxdepth = c->progeny[k]->maxdepth;
        }

    }

    /* Otherwise, collect the data for this cell. */
    else {

      /* Clear the progeny. */
      bzero(c->progeny, sizeof(struct cell *) * 8);
      c->split = 0;
      maxdepth = c->depth;

      /* Get dt_min/dt_max. */
      for (int k = 0; k < count; k++) {
        struct part *p = &parts[k];
        struct xpart *xp = &xparts[k];
        const float h = p->h;
        const int ti_end = p->ti_end;
        xp->x_diff[0] = 0.f;
        xp->x_diff[1] = 0.f;
        xp->x_diff[2] = 0.f;
        if (h > h_max) h_max = h;
        if (ti_end < ti_end_min) ti_end_min = ti_end;
        if (ti_end > ti_end_max) ti_end_max = ti_end;
      }
      for (int k = 0; k < gcount; k++) {
        struct gpart *gp = &gparts[k];
        const int ti_end = gp->ti_end;
        gp->x_diff[0] = 0.f;
        gp->x_diff[1] = 0.f;
        gp->x_diff[2] = 0.f;
        if (ti_end < ti_end_min) ti_end_min = ti_end;
        if (ti_end > ti_end_max) ti_end_max = ti_end;
      }
    }

    /* Set the values for this cell. */
    c->h_max = h_max;
    c->ti_end_min = ti_end_min;
    c->ti_end_max = ti_end_max;
    c->maxdepth = maxdepth;

    /* Set ownership according to the start of the parts array. */
    if (s->nr_parts > 0)
      c->owner =
          ((c->parts - s->parts) % s->nr_parts) * s->nr_queues / s->nr_parts;
    else if (s->nr_gparts > 0)
      c->owner = ((c->gparts - s->gparts) % s->nr_gparts) * s->nr_queues /
                 s->nr_gparts;
    else
      c->owner = 0; /* Ok, there is really nothing on this rank... */
  }
}

/**
 * @brief Return a used cell to the buffer od unused sub-cells.
 *
 * @param s The #space.
 * @param c The #cell.
 */
void space_recycle(struct space *s, struct cell *c) {

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Clear the cell. */
  if (lock_destroy(&c->lock) != 0) error("Failed to destroy spinlock.");

  /* Clear this cell's sort arrays. */
  if (c->sort != NULL) free(c->sort);

  /* Clear the cell data. */
  bzero(c, sizeof(struct cell));

  /* Hook this cell into the buffer. */
  c->next = s->cells_sub;
  s->cells_sub = c;
  s->tot_cells -= 1;

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Get a new empty (sub-)#cell.
 *
 * If there are cells in the buffer, use the one at the end of the linked list.
 * If we have no cells, allocate a new chunk of memory and pick one from there.
 *
 * @param s The #space.
 */
struct cell *space_getcell(struct space *s) {

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Is the buffer empty? */
  if (s->cells_sub == NULL) {
    if (posix_memalign((void *)&s->cells_sub, cell_align,
                       space_cellallocchunk * sizeof(struct cell)) != 0)
      error("Failed to allocate more cells.");

    /* Zero everything for good measure */
    bzero(s->cells_sub, space_cellallocchunk * sizeof(struct cell));

    /* Constructed a linked list */
    for (int k = 0; k < space_cellallocchunk - 1; k++)
      s->cells_sub[k].next = &s->cells_sub[k + 1];
    s->cells_sub[space_cellallocchunk - 1].next = NULL;
  }

  /* Pick off the next cell. */
  struct cell *c = s->cells_sub;
  s->cells_sub = c->next;
  s->tot_cells += 1;

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);

  /* Init some things in the cell we just got. */
  bzero(c, sizeof(struct cell));
  c->nodeID = -1;
  if (lock_init(&c->lock) != 0 || lock_init(&c->glock) != 0)
    error("Failed to initialize cell spinlocks.");

  return c;
}

/**
 * @brief Initialises all the particles by setting them into a valid state
 *
 * Calls hydro_first_init_part() on all the particles
 */
void space_init_parts(struct space *s) {

  const size_t nr_parts = s->nr_parts;
  struct part *restrict p = s->parts;
  struct xpart *restrict xp = s->xparts;

  for (size_t i = 0; i < nr_parts; ++i) {

#ifdef HYDRO_DIMENSION_2D
    p[i].x[2] = 0.f;
    p[i].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    p[i].x[1] = p[i].x[2] = 0.f;
    p[i].v[1] = p[i].v[2] = 0.f;
#endif

    hydro_first_init_part(&p[i], &xp[i]);
  }
}

/**
 * @brief Initialises all the extra particle data
 *
 * Calls cooling_init_xpart() on all the particles
 */
void space_init_xparts(struct space *s) {

  const size_t nr_parts = s->nr_parts;
  struct part *restrict p = s->parts;
  struct xpart *restrict xp = s->xparts;

  for (size_t i = 0; i < nr_parts; ++i) {

    cooling_init_part(&p[i], &xp[i]);
  }
}

/**
 * @brief Initialises all the g-particles by setting them into a valid state
 *
 * Calls gravity_first_init_gpart() on all the particles
 */
void space_init_gparts(struct space *s) {

  const size_t nr_gparts = s->nr_gparts;
  struct gpart *restrict gp = s->gparts;

  for (size_t i = 0; i < nr_gparts; ++i) {

#ifdef HYDRO_DIMENSION_2D
    gp[i].x[2] = 0.f;
    gp[i].v_full[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    gp[i].x[1] = gp[i].x[2] = 0.f;
    gp[i].v_full[1] = gp[i].v_full[2] = 0.f;
#endif

    gravity_first_init_gpart(&gp[i]);
  }
}

/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param s The #space to initialize.
 * @param params The parsed parameter file.
 * @param dim Spatial dimensions of the domain.
 * @param parts Array of Gas particles.
 * @param gparts Array of Gravity particles.
 * @param Npart The number of Gas particles in the space.
 * @param Ngpart The number of Gravity particles in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param gravity flag whether we are doing gravity or not.
 * @param verbose Print messages to stdout or not.
 * @param dry_run If 1, just initialise stuff, don't do anything with the parts.
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_splitsize
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */
void space_init(struct space *s, const struct swift_params *params,
                double dim[3], struct part *parts, struct gpart *gparts,
                size_t Npart, size_t Ngpart, int periodic, int gravity,
                int verbose, int dry_run) {

  /* Clean-up everything */
  bzero(s, sizeof(struct space));

  /* Store everything in the space. */
  s->dim[0] = dim[0];
  s->dim[1] = dim[1];
  s->dim[2] = dim[2];
  s->periodic = periodic;
  s->gravity = gravity;
  s->nr_parts = Npart;
  s->size_parts = Npart;
  s->parts = parts;
  s->nr_gparts = Ngpart;
  s->size_gparts = Ngpart;
  s->gparts = gparts;
  s->cell_min = parser_get_param_double(params, "SPH:max_smoothing_length");
  s->nr_queues = 1; /* Temporary value until engine construction */

  /* Get the constants for the scheduler */
  space_maxsize = parser_get_opt_param_int(params, "Scheduler:cell_max_size",
                                           space_maxsize_default);
  space_subsize = parser_get_opt_param_int(params, "Scheduler:cell_sub_size",
                                           space_subsize_default);
  space_splitsize = parser_get_opt_param_int(
      params, "Scheduler:cell_split_size", space_splitsize_default);
  space_maxcount = parser_get_opt_param_int(params, "Scheduler:cell_max_count",
                                            space_maxcount_default);
  if (verbose)
    message("max_size set to %d, sub_size set to %d, split_size set to %d",
            space_maxsize, space_subsize, space_splitsize);

  /* Check that we have enough cells */
  if (s->cell_min * 3 > dim[0] || s->cell_min * 3 > dim[1] ||
      s->cell_min * 3 > dim[2])
    error(
        "Maximal smoothing length (%e) too large. Needs to be "
        "smaller than 1/3 the simulation box size [%e %e %e]",
        s->cell_min, dim[0], dim[1], dim[2]);

  /* Apply h scaling */
  const double scaling =
      parser_get_opt_param_double(params, "InitialConditions:h_scaling", 1.0);
  if (scaling != 1.0 && !dry_run) {
    message("Re-scaling smoothing lengths by a factor %e", scaling);
    for (size_t k = 0; k < Npart; k++) parts[k].h *= scaling;
  }

  /* Apply shift */
  double shift[3] = {0.0, 0.0, 0.0};
  shift[0] =
      parser_get_opt_param_double(params, "InitialConditions:shift_x", 0.0);
  shift[1] =
      parser_get_opt_param_double(params, "InitialConditions:shift_y", 0.0);
  shift[2] =
      parser_get_opt_param_double(params, "InitialConditions:shift_z", 0.0);
  if ((shift[0] != 0. || shift[1] != 0. || shift[2] != 0.) && !dry_run) {
    message("Shifting particles by [%e %e %e]", shift[0], shift[1], shift[2]);
    for (size_t k = 0; k < Npart; k++) {
      parts[k].x[0] += shift[0];
      parts[k].x[1] += shift[1];
      parts[k].x[2] += shift[2];
    }
    for (size_t k = 0; k < Ngpart; k++) {
      gparts[k].x[0] += shift[0];
      gparts[k].x[1] += shift[1];
      gparts[k].x[2] += shift[2];
    }
  }

  if (!dry_run) {

    /* Check that all the part positions are reasonable, wrap if periodic. */
    if (periodic) {
      for (size_t k = 0; k < Npart; k++)
        for (int j = 0; j < 3; j++) {
          while (parts[k].x[j] < 0) parts[k].x[j] += dim[j];
          while (parts[k].x[j] >= dim[j]) parts[k].x[j] -= dim[j];
        }
    } else {
      for (size_t k = 0; k < Npart; k++)
        for (int j = 0; j < 3; j++)
          if (parts[k].x[j] < 0 || parts[k].x[j] >= dim[j])
            error("Not all particles are within the specified domain.");
    }

    /* Same for the gparts */
    if (periodic) {
      for (size_t k = 0; k < Ngpart; k++)
        for (int j = 0; j < 3; j++) {
          while (gparts[k].x[j] < 0) gparts[k].x[j] += dim[j];
          while (gparts[k].x[j] >= dim[j]) gparts[k].x[j] -= dim[j];
        }
    } else {
      for (size_t k = 0; k < Ngpart; k++)
        for (int j = 0; j < 3; j++)
          if (gparts[k].x[j] < 0 || gparts[k].x[j] >= dim[j])
            error("Not all g-particles are within the specified domain.");
    }
  }

  /* Allocate the extra parts array. */
  if (Npart > 0) {
    if (posix_memalign((void *)&s->xparts, xpart_align,
                       Npart * sizeof(struct xpart)) != 0)
      error("Failed to allocate xparts.");
    bzero(s->xparts, Npart * sizeof(struct xpart));
  }

  /* Set the particles in a state where they are ready for a run */
  space_init_parts(s);
  space_init_xparts(s);
  space_init_gparts(s);

  /* Init the space lock. */
  if (lock_init(&s->lock) != 0) error("Failed to create space spin-lock.");

  /* Build the cells and the tasks. */
  if (!dry_run) space_regrid(s, s->cell_min, verbose);
}

/**
 * @brief Cleans-up all the cell links in the space
 *
 * Expensive funtion. Should only be used for debugging purposes.
 */
void space_link_cleanup(struct space *s) {

  /* Recursively apply the cell link cleaning routine */
  space_map_cells_pre(s, 1, cell_clean_links, NULL);
}

/**
 * @brief Frees up the memory allocated for this #space
 */
void space_clean(struct space *s) {

  for (int i = 0; i < s->nr_cells; ++i) cell_clean(&s->cells_top[i]);
  free(s->cells_top);
  free(s->parts);
  free(s->xparts);
  free(s->gparts);
}
