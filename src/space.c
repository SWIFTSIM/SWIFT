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
#include "memswap.h"
#include "minmax.h"
#include "multipole.h"
#include "runner.h"
#include "sort_part.h"
#include "stars.h"
#include "threadpool.h"
#include "tools.h"

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize_pair = space_subsize_pair_default;
int space_subsize_self = space_subsize_self_default;
int space_subsize_self_grav = space_subsize_self_grav_default;
int space_maxsize = space_maxsize_default;

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
  struct spart *sparts;
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
 * @param cell_rec_begin Pointer to the start of the list of cells to recycle.
 * @param cell_rec_end Pointer to the end of the list of cells to recycle.
 * @param multipole_rec_begin Pointer to the start of the list of multipoles to
 * recycle.
 * @param multipole_rec_end Pointer to the end of the list of multipoles to
 * recycle.
 */
void space_rebuild_recycle_rec(struct space *s, struct cell *c,
                               struct cell **cell_rec_begin,
                               struct cell **cell_rec_end,
                               struct gravity_tensors **multipole_rec_begin,
                               struct gravity_tensors **multipole_rec_end) {
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        space_rebuild_recycle_rec(s, c->progeny[k], cell_rec_begin,
                                  cell_rec_end, multipole_rec_begin,
                                  multipole_rec_end);

        c->progeny[k]->next = *cell_rec_begin;
        *cell_rec_begin = c->progeny[k];

        if (s->gravity) {
          c->progeny[k]->multipole->next = *multipole_rec_begin;
          *multipole_rec_begin = c->progeny[k]->multipole;
        }

        if (*cell_rec_end == NULL) *cell_rec_end = *cell_rec_begin;
        if (s->gravity && *multipole_rec_end == NULL)
          *multipole_rec_end = *multipole_rec_begin;

        c->progeny[k]->multipole = NULL;
        c->progeny[k] = NULL;
      }
}

void space_rebuild_recycle_mapper(void *map_data, int num_elements,
                                  void *extra_data) {

  struct space *s = (struct space *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int k = 0; k < num_elements; k++) {
    struct cell *c = &cells[k];
    struct cell *cell_rec_begin = NULL, *cell_rec_end = NULL;
    struct gravity_tensors *multipole_rec_begin = NULL,
                           *multipole_rec_end = NULL;
    space_rebuild_recycle_rec(s, c, &cell_rec_begin, &cell_rec_end,
                              &multipole_rec_begin, &multipole_rec_end);
    if (cell_rec_begin != NULL)
      space_recycle_list(s, cell_rec_begin, cell_rec_end, multipole_rec_begin,
                         multipole_rec_end);
    c->sorts = NULL;
    c->nr_tasks = 0;
    c->density = NULL;
    c->gradient = NULL;
    c->force = NULL;
    c->grav = NULL;
    c->dx_max_part = 0.0f;
    c->dx_max_gpart = 0.0f;
    c->dx_max_sort = 0.0f;
    c->sorted = 0;
    c->count = 0;
    c->gcount = 0;
    c->scount = 0;
    c->init_grav = NULL;
    c->extra_ghost = NULL;
    c->ghost_in = NULL;
    c->ghost_out = NULL;
    c->ghost = NULL;
    c->kick1 = NULL;
    c->kick2 = NULL;
    c->timestep = NULL;
    c->drift_part = NULL;
    c->drift_gpart = NULL;
    c->cooling = NULL;
    c->sourceterms = NULL;
    c->grav_ghost[0] = NULL;
    c->grav_ghost[1] = NULL;
    c->grav_long_range = NULL;
    c->grav_down = NULL;
    c->super = c;
    c->parts = NULL;
    c->xparts = NULL;
    c->gparts = NULL;
    c->sparts = NULL;
    for (int i = 0; i < 13; i++)
      if (c->sort[i] != NULL) {
        free(c->sort[i]);
        c->sort[i] = NULL;
      }
#if WITH_MPI
    c->recv_xv = NULL;
    c->recv_rho = NULL;
    c->recv_gradient = NULL;
    c->recv_ti = NULL;

    c->send_xv = NULL;
    c->send_rho = NULL;
    c->send_gradient = NULL;
    c->send_ti = NULL;
#endif
  }
}

/**
 * @brief Free up any allocated cells.
 */
void space_free_cells(struct space *s) {
  threadpool_map(&s->e->threadpool, space_rebuild_recycle_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), 0, s);
  s->maxdepth = 0;
}

/**
 * @brief Re-build the top-level cell grid.
 *
 * @param s The #space.
 * @param verbose Print messages to stdout or not.
 */
void space_regrid(struct space *s, int verbose) {

  const size_t nr_parts = s->nr_parts;
  const ticks tic = getticks();
  const integertime_t ti_old = (s->e != NULL) ? s->e->ti_old : 0;

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
  if (verbose) message("h_max is %.3e (cell_min=%.3e).", h_max, s->cell_min);

  /* Get the new putative cell dimensions. */
  const int cdim[3] = {
      floor(s->dim[0] /
            fmax(h_max * kernel_gamma * space_stretch, s->cell_min)),
      floor(s->dim[1] /
            fmax(h_max * kernel_gamma * space_stretch, s->cell_min)),
      floor(s->dim[2] /
            fmax(h_max * kernel_gamma * space_stretch, s->cell_min))};

  /* Check if we have enough cells for periodicity. */
  if (s->periodic && (cdim[0] < 3 || cdim[1] < 3 || cdim[2] < 3))
    error(
        "Must have at least 3 cells in each spatial dimension when periodicity "
        "is switched on.\nThis error is often caused by any of the "
        "followings:\n"
        " - too few particles to generate a sensible grid,\n"
        " - the initial value of 'Scheduler:max_top_level_cells' is too "
        "small,\n"
        " - the (minimal) time-step is too large leading to particles with "
        "predicted smoothing lengths too large for the box size,\n"
        " - particles with velocities so large that they move by more than two "
        "box sizes per time-step.\n");

  /* Check if we have enough cells for periodic gravity. */
  if (s->gravity && s->periodic && (cdim[0] < 8 || cdim[1] < 8 || cdim[2] < 8))
    error(
        "Must have at least 8 cells in each spatial dimension when periodic "
        "gravity is switched on.\nThis error is often caused by any of the "
        "followings:\n"
        " - too few particles to generate a sensible grid,\n"
        " - the initial value of 'Scheduler:max_top_level_cells' is too "
        "small,\n"
        " - the (minimal) time-step is too large leading to particles with "
        "predicted smoothing lengths too large for the box size,\n"
        " - particles with velocities so large that they move by more than two "
        "box sizes per time-step.\n");

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

/* Be verbose about this. */
#ifdef SWIFT_DEBUG_CHECKS
    message("(re)griding space cdim=(%d %d %d)", cdim[0], cdim[1], cdim[2]);
    fflush(stdout);
#endif

    /* Free the old cells, if they were allocated. */
    if (s->cells_top != NULL) {
      space_free_cells(s);
      free(s->cells_top);
      free(s->multipoles_top);
    }

    /* Also free the task arrays, these will be regenerated and we can use the
     * memory while copying the particle arrays. */
    if (s->e != NULL) scheduler_free_tasks(&s->e->sched);

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
      error("Failed to allocate top-level cells.");
    bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

    /* Allocate the multipoles for the top-level cells. */
    if (s->gravity) {
      if (posix_memalign((void *)&s->multipoles_top, multipole_align,
                         s->nr_cells * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate top-level multipoles.");
      bzero(s->multipoles_top, s->nr_cells * sizeof(struct gravity_tensors));
    }

    /* Set the cells' locks */
    for (int k = 0; k < s->nr_cells; k++) {
      if (lock_init(&s->cells_top[k].lock) != 0)
        error("Failed to init spinlock for hydro.");
      if (lock_init(&s->cells_top[k].glock) != 0)
        error("Failed to init spinlock for gravity.");
      if (lock_init(&s->cells_top[k].mlock) != 0)
        error("Failed to init spinlock for multipoles.");
      if (lock_init(&s->cells_top[k].slock) != 0)
        error("Failed to init spinlock for stars.");
    }

    /* Set the cell location and sizes. */
    for (int i = 0; i < cdim[0]; i++)
      for (int j = 0; j < cdim[1]; j++)
        for (int k = 0; k < cdim[2]; k++) {
          const size_t cid = cell_getid(cdim, i, j, k);
          struct cell *restrict c = &s->cells_top[cid];
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
          c->scount = 0;
          c->super = c;
          c->ti_old_part = ti_old;
          c->ti_old_gpart = ti_old;
          c->ti_old_multipole = ti_old;
          if (s->gravity) c->multipole = &s->multipoles_top[cid];
        }

    /* Be verbose about the change. */
    if (verbose)
      message("set cell dimensions to [ %i %i %i ].", cdim[0], cdim[1],
              cdim[2]);

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
#endif /* WITH_MPI */

    // message( "rebuilding upper-level cells took %.3f %s." ,
    // clocks_from_ticks(double)(getticks() - tic), clocks_getunit());

  }      /* re-build upper-level cells? */
  else { /* Otherwise, just clean up the cells. */

    /* Free the old cells, if they were allocated. */
    space_free_cells(s);
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Re-build the cells as well as the tasks.
 *
 * @param s The #space in which to update the cells.
 * @param verbose Print messages to stdout or not
 *
 */
void space_rebuild(struct space *s, int verbose) {

  const ticks tic = getticks();

/* Be verbose about this. */
#ifdef SWIFT_DEBUG_CHECKS
  if (s->e->nodeID == 0 || verbose) message("(re)building space");
  fflush(stdout);
#endif

  /* Re-grid if necessary, or just re-set the cell data. */
  space_regrid(s, verbose);

  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;
  struct cell *restrict cells_top = s->cells_top;
  const integertime_t ti_old = (s->e != NULL) ? s->e->ti_old : 0;

  /* Run through the particles and get their cell index. Allocates
     an index that is larger than the number of particles to avoid
     re-allocating after shuffling. */
  const size_t ind_size = s->size_parts + 100;
  int *ind;
  if ((ind = (int *)malloc(sizeof(int) * ind_size)) == NULL)
    error("Failed to allocate temporary particle indices.");
  if (s->size_parts > 0) space_parts_get_cell_index(s, ind, cells_top, verbose);

  /* Run through the gravity particles and get their cell index. */
  const size_t gind_size = s->size_gparts + 100;
  int *gind;
  if ((gind = (int *)malloc(sizeof(int) * gind_size)) == NULL)
    error("Failed to allocate temporary g-particle indices.");
  if (s->size_gparts > 0)
    space_gparts_get_cell_index(s, gind, cells_top, verbose);

  /* Run through the star particles and get their cell index. */
  const size_t sind_size = s->size_sparts + 100;
  int *sind;
  if ((sind = (int *)malloc(sizeof(int) * sind_size)) == NULL)
    error("Failed to allocate temporary s-particle indices.");
  if (s->size_sparts > 0)
    space_sparts_get_cell_index(s, sind, cells_top, verbose);

#ifdef WITH_MPI
  const int local_nodeID = s->e->nodeID;

  /* Move non-local parts to the end of the list. */
  for (size_t k = 0; k < nr_parts;) {
    if (cells_top[ind[k]].nodeID != local_nodeID) {
      nr_parts -= 1;
      /* Swap the particle */
      const struct part tp = s->parts[k];
      s->parts[k] = s->parts[nr_parts];
      s->parts[nr_parts] = tp;
      /* Swap the link with the gpart */
      if (s->parts[k].gpart != NULL) {
        s->parts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->parts[nr_parts].gpart != NULL) {
        s->parts[nr_parts].gpart->id_or_neg_offset = -nr_parts;
      }
      /* Swap the xpart */
      const struct xpart txp = s->xparts[k];
      s->xparts[k] = s->xparts[nr_parts];
      s->xparts[nr_parts] = txp;
      /* Swap the index */
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

  /* Move non-local sparts to the end of the list. */
  for (size_t k = 0; k < nr_sparts;) {
    if (cells_top[sind[k]].nodeID != local_nodeID) {
      nr_sparts -= 1;
      /* Swap the particle */
      const struct spart tp = s->sparts[k];
      s->sparts[k] = s->sparts[nr_sparts];
      s->sparts[nr_sparts] = tp;
      /* Swap the link with the gpart */
      if (s->sparts[k].gpart != NULL) {
        s->sparts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->sparts[nr_sparts].gpart != NULL) {
        s->sparts[nr_sparts].gpart->id_or_neg_offset = -nr_sparts;
      }
      /* Swap the index */
      const int t = sind[k];
      sind[k] = sind[nr_sparts];
      sind[nr_sparts] = t;
    } else {
      /* Increment when not exchanging otherwise we need to retest "k".*/
      k++;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all sparts are in the correct place (untested). */
  for (size_t k = 0; k < nr_sparts; k++) {
    if (cells_top[sind[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local sparts to send list");
    }
  }
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    if (cells_top[sind[k]].nodeID == local_nodeID) {
      error("Failed to remove local sparts from send list");
    }
  }
#endif

  /* Move non-local gparts to the end of the list. */
  for (size_t k = 0; k < nr_gparts;) {
    if (cells_top[gind[k]].nodeID != local_nodeID) {
      nr_gparts -= 1;
      /* Swap the particle */
      const struct gpart tp = s->gparts[k];
      s->gparts[k] = s->gparts[nr_gparts];
      s->gparts[nr_gparts] = tp;
      /* Swap the link with part/spart */
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_star) {
        s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }
      if (s->gparts[nr_gparts].type == swift_type_gas) {
        s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      } else if (s->gparts[nr_gparts].type == swift_type_star) {
        s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      }
      /* Swap the index */
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
  size_t nr_sparts_exchanged = s->nr_sparts - nr_sparts;
  engine_exchange_strays(s->e, nr_parts, &ind[nr_parts], &nr_parts_exchanged,
                         nr_gparts, &gind[nr_gparts], &nr_gparts_exchanged,
                         nr_sparts, &sind[nr_sparts], &nr_sparts_exchanged);

  /* Set the new particle counts. */
  s->nr_parts = nr_parts + nr_parts_exchanged;
  s->nr_gparts = nr_gparts + nr_gparts_exchanged;
  s->nr_sparts = nr_sparts + nr_sparts_exchanged;

  /* Re-allocate the index array for the parts if needed.. */
  if (s->nr_parts + 1 > ind_size) {
    int *ind_new;
    if ((ind_new = (int *)malloc(sizeof(int) * (s->nr_parts + 1))) == NULL)
      error("Failed to allocate temporary particle indices.");
    memcpy(ind_new, ind, sizeof(int) * nr_parts);
    free(ind);
    ind = ind_new;
  }

  /* Re-allocate the index array for the sparts if needed.. */
  if (s->nr_sparts + 1 > sind_size) {
    int *sind_new;
    if ((sind_new = (int *)malloc(sizeof(int) * (s->nr_sparts + 1))) == NULL)
      error("Failed to allocate temporary s-particle indices.");
    memcpy(sind_new, sind, sizeof(int) * nr_sparts);
    free(sind);
    sind = sind_new;
  }

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih[3] = {s->iwidth[0], s->iwidth[1], s->iwidth[2]};

  /* Assign each received part to its cell. */
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    const struct part *const p = &s->parts[k];
    ind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[ind[k]].nodeID != local_nodeID)
      error("Received part that does not belong to me (nodeID=%i).",
            cells_top[ind[k]].nodeID);
#endif
  }
  nr_parts = s->nr_parts;

  /* Assign each received spart to its cell. */
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    const struct spart *const sp = &s->sparts[k];
    sind[k] =
        cell_getid(cdim, sp->x[0] * ih[0], sp->x[1] * ih[1], sp->x[2] * ih[2]);
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[sind[k]].nodeID != local_nodeID)
      error("Received s-part that does not belong to me (nodeID=%i).",
            cells_top[sind[k]].nodeID);
#endif
  }
  nr_sparts = s->nr_sparts;

#endif /* WITH_MPI */

  /* Sort the parts according to their cells. */
  if (nr_parts > 0)
    space_parts_sort(s, ind, nr_parts, 0, s->nr_cells - 1, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < nr_parts; k++) {
    const struct part *p = &s->parts[k];

    /* New cell index */
    const int new_ind =
        cell_getid(s->cdim, p->x[0] * s->iwidth[0], p->x[1] * s->iwidth[1],
                   p->x[2] * s->iwidth[2]);

    /* New cell of this part */
    const struct cell *c = &s->cells_top[new_ind];

    if (ind[k] != new_ind)
      error("part's new cell index not matching sorted index.");

    if (p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] < c->loc[1] || p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] < c->loc[2] || p->x[2] > c->loc[2] + c->width[2])
      error("part not sorted into the right top-level cell!");
  }
#endif

  /* Sort the sparts according to their cells. */
  if (nr_sparts > 0)
    space_sparts_sort(s, sind, nr_sparts, 0, s->nr_cells - 1, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < nr_sparts; k++) {
    const struct spart *sp = &s->sparts[k];

    /* New cell index */
    const int new_sind =
        cell_getid(s->cdim, sp->x[0] * s->iwidth[0], sp->x[1] * s->iwidth[1],
                   sp->x[2] * s->iwidth[2]);

    /* New cell of this spart */
    const struct cell *c = &s->cells_top[new_sind];

    if (sind[k] != new_sind)
      error("spart's new cell index not matching sorted index.");

    if (sp->x[0] < c->loc[0] || sp->x[0] > c->loc[0] + c->width[0] ||
        sp->x[1] < c->loc[1] || sp->x[1] > c->loc[1] + c->width[1] ||
        sp->x[2] < c->loc[2] || sp->x[2] > c->loc[2] + c->width[2])
      error("spart not sorted into the right top-level cell!");
  }
#endif

  /* Re-link the gparts to their (s-)particles. */
  if (nr_parts > 0 && nr_gparts > 0)
    part_relink_gparts_to_parts(s->parts, nr_parts, 0);
  if (nr_sparts > 0 && nr_gparts > 0)
    part_relink_gparts_to_sparts(s->sparts, nr_sparts, 0);

  /* Extract the cell counts from the sorted indices. */
  size_t last_index = 0;
  ind[nr_parts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_parts; k++) {
    if (ind[k] < ind[k + 1]) {
      cells_top[ind[k]].count = k - last_index + 1;
      last_index = k + 1;
    }
  }

  /* Extract the cell counts from the sorted indices. */
  size_t last_sindex = 0;
  sind[nr_sparts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_sparts; k++) {
    if (sind[k] < sind[k + 1]) {
      cells_top[sind[k]].scount = k - last_sindex + 1;
      last_sindex = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  free(ind);
  free(sind);

#ifdef WITH_MPI

  /* Re-allocate the index array for the gparts if needed.. */
  if (s->nr_gparts + 1 > gind_size) {
    int *gind_new;
    if ((gind_new = (int *)malloc(sizeof(int) * (s->nr_gparts + 1))) == NULL)
      error("Failed to allocate temporary g-particle indices.");
    memcpy(gind_new, gind, sizeof(int) * nr_gparts);
    free(gind);
    gind = gind_new;
  }

  /* Assign each received gpart to its cell. */
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    const struct gpart *const p = &s->gparts[k];
    gind[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);

#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[gind[k]].nodeID != s->e->nodeID)
      error("Received g-part that does not belong to me (nodeID=%i).",
            cells_top[gind[k]].nodeID);
#endif
  }
  nr_gparts = s->nr_gparts;

#endif /* WITH_MPI */

  /* Sort the gparts according to their cells. */
  if (nr_gparts > 0)
    space_gparts_sort(s, gind, nr_gparts, 0, s->nr_cells - 1, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < nr_gparts; k++) {
    const struct gpart *gp = &s->gparts[k];

    /* New cell index */
    const int new_gind =
        cell_getid(s->cdim, gp->x[0] * s->iwidth[0], gp->x[1] * s->iwidth[1],
                   gp->x[2] * s->iwidth[2]);

    /* New cell of this gpart */
    const struct cell *c = &s->cells_top[new_gind];

    if (gind[k] != new_gind)
      error("gpart's new cell index not matching sorted index.");

    if (gp->x[0] < c->loc[0] || gp->x[0] > c->loc[0] + c->width[0] ||
        gp->x[1] < c->loc[1] || gp->x[1] > c->loc[1] + c->width[1] ||
        gp->x[2] < c->loc[2] || gp->x[2] > c->loc[2] + c->width[2])
      error("gpart not sorted into the right top-level cell!");
  }
#endif

  /* Re-link the parts. */
  if (nr_parts > 0 && nr_gparts > 0)
    part_relink_parts_to_gparts(s->gparts, nr_gparts, s->parts);

  /* Re-link the sparts. */
  if (nr_sparts > 0 && nr_gparts > 0)
    part_relink_sparts_to_gparts(s->gparts, nr_gparts, s->sparts);

  /* Extract the cell counts from the sorted indices. */
  size_t last_gindex = 0;
  gind[nr_gparts] = s->nr_cells;
  for (size_t k = 0; k < nr_gparts; k++) {
    if (gind[k] < gind[k + 1]) {
      cells_top[gind[k]].gcount = k - last_gindex + 1;
      last_gindex = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  free(gind);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  if ((nr_gparts > 0 && nr_parts > 0) || (nr_gparts > 0 && nr_sparts > 0))
    part_verify_links(s->parts, s->gparts, s->sparts, nr_parts, nr_gparts,
                      nr_sparts, verbose);
#endif

  /* Hook the cells up to the parts. */
  // tic = getticks();
  struct part *finger = s->parts;
  struct xpart *xfinger = s->xparts;
  struct gpart *gfinger = s->gparts;
  struct spart *sfinger = s->sparts;
  for (int k = 0; k < s->nr_cells; k++) {
    struct cell *restrict c = &cells_top[k];
    c->ti_old_part = ti_old;
    c->ti_old_gpart = ti_old;
    c->ti_old_multipole = ti_old;
    if (c->nodeID == engine_rank) {
      c->parts = finger;
      c->xparts = xfinger;
      c->gparts = gfinger;
      c->sparts = sfinger;
      finger = &finger[c->count];
      xfinger = &xfinger[c->count];
      gfinger = &gfinger[c->gcount];
      sfinger = &sfinger[c->scount];
    }
  }
  // message( "hooking up cells took %.3f %s." ,
  // clocks_from_ticks(getticks() - tic), clocks_getunit());

  /* At this point, we have the upper-level cells, old or new. Now make
     sure that the parts in each cell are ok. */
  space_split(s, cells_top, s->nr_cells, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the multipole construction went OK */
  if (s->gravity)
    for (int k = 0; k < s->nr_cells; k++)
      cell_check_multipole(&s->cells_top[k], NULL);
#endif

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
                 sizeof(struct cell), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief #threadpool mapper function to sanitize the cells
 *
 * @param map_data Pointers towards the top-level cells.
 * @param num_cells The number of top-level cells.
 * @param extra_data Unused parameters.
 */
void space_sanitize_mapper(void *map_data, int num_cells, void *extra_data) {
  /* Unpack the inputs. */
  struct cell *cells_top = (struct cell *)map_data;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[ind];
    cell_sanitize(c, 0);
  }
}

/**
 * @brief Runs through the top-level cells and sanitize their h values
 *
 * @param s The #space to act upon.
 */
void space_sanitize(struct space *s) {

  if (s->e->nodeID == 0) message("Cleaning up unreasonable values of h");

  threadpool_map(&s->e->threadpool, space_sanitize_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), 0, NULL);
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

#ifdef SWIFT_DEBUG_CHECKS
    if (pos_x > dim_x || pos_y > dim_y || pos_z > pos_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

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
 * @brief #threadpool mapper function to compute the s-particle cell indices.
 *
 * @param map_data Pointer towards the s-particles.
 * @param nr_sparts The number of s-particles to treat.
 * @param extra_data Pointers to the space and index list
 */
void space_sparts_get_cell_index_mapper(void *map_data, int nr_sparts,
                                        void *extra_data) {

  /* Unpack the data */
  struct spart *restrict sparts = (struct spart *)map_data;
  struct index_data *data = (struct index_data *)extra_data;
  struct space *s = data->s;
  int *const ind = data->ind + (ptrdiff_t)(sparts - s->sparts);

  /* Get some constants */
  const double dim_x = s->dim[0];
  const double dim_y = s->dim[1];
  const double dim_z = s->dim[2];
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih_x = s->iwidth[0];
  const double ih_y = s->iwidth[1];
  const double ih_z = s->iwidth[2];

  for (int k = 0; k < nr_sparts; k++) {

    /* Get the particle */
    struct spart *restrict sp = &sparts[k];

    const double old_pos_x = sp->x[0];
    const double old_pos_y = sp->x[1];
    const double old_pos_z = sp->x[2];

    /* Put it back into the simulation volume */
    const double pos_x = box_wrap(old_pos_x, 0.0, dim_x);
    const double pos_y = box_wrap(old_pos_y, 0.0, dim_y);
    const double pos_z = box_wrap(old_pos_z, 0.0, dim_z);

    /* Get its cell index */
    const int index =
        cell_getid(cdim, pos_x * ih_x, pos_y * ih_y, pos_z * ih_z);
    ind[k] = index;

    /* Update the position */
    sp->x[0] = pos_x;
    sp->x[1] = pos_y;
    sp->x[2] = pos_z;
  }
}

/**
 * @brief Computes the cell index of all the particles.
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
                 s->nr_parts, sizeof(struct part), 0, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the g-particles.
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
                 s->gparts, s->nr_gparts, sizeof(struct gpart), 0, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the s-particles.
 *
 * @param s The #space.
 * @param sind The array of indices to fill.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_sparts_get_cell_index(struct space *s, int *sind, struct cell *cells,
                                 int verbose) {

  const ticks tic = getticks();

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = sind;

  threadpool_map(&s->e->threadpool, space_sparts_get_cell_index_mapper,
                 s->sparts, s->nr_sparts, sizeof(struct spart), 0, &data);

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
  if (s->e->nodeID == 0 || verbose) message("Sorting succeeded.");
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
          memswap(&ind[ii], &ind[jj], sizeof(int));
          memswap(&parts[ii], &parts[jj], sizeof(struct part));
          memswap(&xparts[ii], &xparts[jj], sizeof(struct xpart));
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify space_sort_struct. */
      if (i != j) {
        for (int k = i; k <= jj; k++) {
          if (ind[k] > pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.", k,
                ind[k], pivot, i, j);
            error("Partition failed (<=pivot).");
          }
        }
        for (int k = jj + 1; k <= j; k++) {
          if (ind[k] <= pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.", k,
                ind[k], pivot, i, j);
            error("Partition failed (>pivot).");
          }
        }
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
 * @brief Sort the s-particles according to the given indices.
 *
 * @param s The #space.
 * @param ind The indices with respect to which the #spart are sorted.
 * @param N The number of parts
 * @param min Lowest index.
 * @param max highest index.
 * @param verbose Are we talkative ?
 */
void space_sparts_sort(struct space *s, int *ind, size_t N, int min, int max,
                       int verbose) {

  const ticks tic = getticks();

  /* Populate a parallel_sort structure with the input data */
  struct parallel_sort sort_struct;
  sort_struct.sparts = s->sparts;
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
  threadpool_map(&s->e->threadpool, space_sparts_sort_mapper, &sort_struct,
                 s->e->threadpool.num_threads, 0, 1, NULL);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify space_sort_struct. */
  for (size_t i = 1; i < N; i++)
    if (ind[i - 1] > ind[i])
      error("Sorting failed (ind[%zu]=%i,ind[%zu]=%i), min=%i, max=%i.", i - 1,
            ind[i - 1], i, ind[i], min, max);
  if (s->e->nodeID == 0 || verbose) message("Sorting succeeded.");
#endif

  /* Clean up. */
  free(sort_struct.stack);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_sparts_sort_mapper(void *map_data, int num_elements,
                              void *extra_data) {

  /* Unpack the mapping data. */
  struct parallel_sort *sort_struct = (struct parallel_sort *)map_data;

  /* Pointers to the sorting data. */
  int *ind = sort_struct->ind;
  struct spart *sparts = sort_struct->sparts;

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
          memswap(&ind[ii], &ind[jj], sizeof(int));
          memswap(&sparts[ii], &sparts[jj], sizeof(struct spart));
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify space_sort_struct. */
      if (i != j) {
        for (int k = i; k <= jj; k++) {
          if (ind[k] > pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li "
                "min=%i max=%i.",
                k, ind[k], pivot, i, j, min, max);
            error("Partition failed (<=pivot).");
          }
        }
        for (int k = jj + 1; k <= j; k++) {
          if (ind[k] <= pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.", k,
                ind[k], pivot, i, j);
            error("Partition failed (>pivot).");
          }
        }
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
  if (s->e->nodeID == 0 || verbose) message("Sorting succeeded.");
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
          memswap(&ind[ii], &ind[jj], sizeof(int));
          memswap(&gparts[ii], &gparts[jj], sizeof(struct gpart));
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Verify space_sort_struct. */
      if (i != j) {
        for (int k = i; k <= jj; k++) {
          if (ind[k] > pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.", k,
                ind[k], pivot, i, j);
            error("Partition failed (<=pivot).");
          }
        }
        for (int k = jj + 1; k <= j; k++) {
          if (ind[k] <= pivot) {
            message(
                "sorting failed at k=%i, ind[k]=%i, pivot=%i, i=%li, j=%li.", k,
                ind[k], pivot, i, j);
            error("Partition failed (>pivot).");
          }
        }
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

  for (int i = 0; i < 13; i++)
    if (c->sort[i] != NULL) {
      free(c->sort[i]);
      c->sort[i] = NULL;
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
 * @brief Recursively split a cell.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param buff A buffer for particle sorting, should be of size at least
 *        c->count or @c NULL.
 * @param sbuff A buffer for particle sorting, should be of size at least
 *        c->scount or @c NULL.
 * @param gbuff A buffer for particle sorting, should be of size at least
 *        c->gcount or @c NULL.
 */
void space_split_recursive(struct space *s, struct cell *c,
                           struct cell_buff *buff, struct cell_buff *sbuff,
                           struct cell_buff *gbuff) {

  const int count = c->count;
  const int gcount = c->gcount;
  const int scount = c->scount;
  const int depth = c->depth;
  int maxdepth = 0;
  float h_max = 0.0f;
  integertime_t ti_end_min = max_nr_timesteps, ti_end_max = 0, ti_beg_max = 0;
  struct part *parts = c->parts;
  struct gpart *gparts = c->gparts;
  struct spart *sparts = c->sparts;
  struct xpart *xparts = c->xparts;
  struct engine *e = s->e;

  /* If the buff is NULL, allocate it, and remember to free it. */
  const int allocate_buffer = (buff == NULL && gbuff == NULL && sbuff == NULL);
  if (allocate_buffer) {
    if (count > 0) {
      if (posix_memalign((void *)&buff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * count) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < count; k++) {
        buff[k].x[0] = parts[k].x[0];
        buff[k].x[1] = parts[k].x[1];
        buff[k].x[2] = parts[k].x[2];
      }
    }
    if (gcount > 0) {
      if (posix_memalign((void *)&gbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * gcount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < gcount; k++) {
        gbuff[k].x[0] = gparts[k].x[0];
        gbuff[k].x[1] = gparts[k].x[1];
        gbuff[k].x[2] = gparts[k].x[2];
      }
    }
    if (scount > 0) {
      if (posix_memalign((void *)&sbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * scount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < scount; k++) {
        sbuff[k].x[0] = sparts[k].x[0];
        sbuff[k].x[1] = sparts[k].x[1];
        sbuff[k].x[2] = sparts[k].x[2];
      }
    }
  }

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
  if (count > space_splitsize || gcount > space_splitsize ||
      scount > space_splitsize) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->count = 0;
      cp->gcount = 0;
      cp->scount = 0;
      cp->ti_old_part = c->ti_old_part;
      cp->ti_old_gpart = c->ti_old_gpart;
      cp->ti_old_multipole = c->ti_old_multipole;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->width[0] = c->width[0] / 2;
      cp->width[1] = c->width[1] / 2;
      cp->width[2] = c->width[2] / 2;
      cp->dmin = c->dmin / 2;
      if (k & 4) cp->loc[0] += cp->width[0];
      if (k & 2) cp->loc[1] += cp->width[1];
      if (k & 1) cp->loc[2] += cp->width[2];
      cp->depth = c->depth + 1;
      cp->split = 0;
      cp->h_max = 0.f;
      cp->dx_max_part = 0.f;
      cp->dx_max_gpart = 0.f;
      cp->dx_max_sort = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->super = NULL;
    }

    /* Split the cell data. */
    cell_split(c, c->parts - s->parts, c->sparts - s->sparts, buff, sbuff,
               gbuff);

    /* Remove any progeny with zero parts. */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k]->count == 0 && c->progeny[k]->gcount == 0 &&
          c->progeny[k]->scount == 0) {
        space_recycle(s, c->progeny[k]);
        c->progeny[k] = NULL;
      } else {
        space_split_recursive(s, c->progeny[k], progeny_buff, progeny_sbuff,
                              progeny_gbuff);
        progeny_buff += c->progeny[k]->count;
        progeny_gbuff += c->progeny[k]->gcount;
        progeny_sbuff += c->progeny[k]->scount;
        h_max = max(h_max, c->progeny[k]->h_max);
        ti_end_min = min(ti_end_min, c->progeny[k]->ti_end_min);
        ti_end_max = max(ti_end_max, c->progeny[k]->ti_end_max);
        ti_beg_max = max(ti_beg_max, c->progeny[k]->ti_beg_max);
        if (c->progeny[k]->maxdepth > maxdepth)
          maxdepth = c->progeny[k]->maxdepth;
      }
    }

    /* Deal with multipole */
    if (s->gravity) {

      /* Reset everything */
      gravity_reset(c->multipole);

      /* Compute CoM of all progenies */
      double CoM[3] = {0., 0., 0.};
      double mass = 0.;

      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct gravity_tensors *m = c->progeny[k]->multipole;
          CoM[0] += m->CoM[0] * m->m_pole.M_000;
          CoM[1] += m->CoM[1] * m->m_pole.M_000;
          CoM[2] += m->CoM[2] * m->m_pole.M_000;
          mass += m->m_pole.M_000;
        }
      }
      c->multipole->CoM[0] = CoM[0] / mass;
      c->multipole->CoM[1] = CoM[1] / mass;
      c->multipole->CoM[2] = CoM[2] / mass;

      /* Now shift progeny multipoles and add them up */
      struct multipole temp;
      double r_max = 0.;
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct cell *cp = c->progeny[k];
          const struct multipole *m = &cp->multipole->m_pole;

          /* Contribution to multipole */
          gravity_M2M(&temp, m, c->multipole->CoM, cp->multipole->CoM);
          gravity_multipole_add(&c->multipole->m_pole, &temp);

          /* Upper limit of max CoM<->gpart distance */
          const double dx = c->multipole->CoM[0] - cp->multipole->CoM[0];
          const double dy = c->multipole->CoM[1] - cp->multipole->CoM[1];
          const double dz = c->multipole->CoM[2] - cp->multipole->CoM[2];
          const double r2 = dx * dx + dy * dy + dz * dz;
          r_max = max(r_max, cp->multipole->r_max + sqrt(r2));
        }
      }
      /* Alternative upper limit of max CoM<->gpart distance */
      const double dx = c->multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
                            ? c->multipole->CoM[0] - c->loc[0]
                            : c->loc[0] + c->width[0] - c->multipole->CoM[0];
      const double dy = c->multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
                            ? c->multipole->CoM[1] - c->loc[1]
                            : c->loc[1] + c->width[1] - c->multipole->CoM[1];
      const double dz = c->multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
                            ? c->multipole->CoM[2] - c->loc[2]
                            : c->loc[2] + c->width[2] - c->multipole->CoM[2];

      /* Take minimum of both limits */
      c->multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));
      c->multipole->r_max_rebuild = c->multipole->r_max;
      c->multipole->CoM_rebuild[0] = c->multipole->CoM[0];
      c->multipole->CoM_rebuild[1] = c->multipole->CoM[1];
      c->multipole->CoM_rebuild[2] = c->multipole->CoM[2];
    } /* Deal with gravity */
  }

  /* Otherwise, collect the data for this cell. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
    maxdepth = c->depth;

    timebin_t time_bin_min = num_time_bins, time_bin_max = 0;

    /* parts: Get dt_min/dt_max and h_max. */
    for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (parts[k].time_bin == time_bin_inhibited)
        error("Inhibited particle present in space_split()");
#endif
      time_bin_min = min(time_bin_min, parts[k].time_bin);
      time_bin_max = max(time_bin_max, parts[k].time_bin);
      h_max = max(h_max, parts[k].h);
    }
    /* parts: Reset x_diff */
    for (int k = 0; k < count; k++) {
      xparts[k].x_diff[0] = 0.f;
      xparts[k].x_diff[1] = 0.f;
      xparts[k].x_diff[2] = 0.f;
    }
    /* gparts: Get dt_min/dt_max, reset x_diff. */
    for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (gparts[k].time_bin == time_bin_inhibited)
        error("Inhibited s-particle present in space_split()");
#endif
      time_bin_min = min(time_bin_min, gparts[k].time_bin);
      time_bin_max = max(time_bin_max, gparts[k].time_bin);

      gparts[k].x_diff[0] = 0.f;
      gparts[k].x_diff[1] = 0.f;
      gparts[k].x_diff[2] = 0.f;
    }
    /* sparts: Get dt_min/dt_max */
    for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[k].time_bin == time_bin_inhibited)
        error("Inhibited g-particle present in space_split()");
#endif
      time_bin_min = min(time_bin_min, sparts[k].time_bin);
      time_bin_max = max(time_bin_max, sparts[k].time_bin);
    }

    /* Convert into integer times */
    ti_end_min = get_integer_time_end(e->ti_current, time_bin_min);
    ti_end_max = get_integer_time_end(e->ti_current, time_bin_max);
    ti_beg_max = get_integer_time_begin(e->ti_current + 1, time_bin_max);

    /* Construct the multipole and the centre of mass*/
    if (s->gravity) {
      if (gcount > 0) {
        gravity_P2M(c->multipole, c->gparts, c->gcount);
        const double dx = c->multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
                              ? c->multipole->CoM[0] - c->loc[0]
                              : c->loc[0] + c->width[0] - c->multipole->CoM[0];
        const double dy = c->multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
                              ? c->multipole->CoM[1] - c->loc[1]
                              : c->loc[1] + c->width[1] - c->multipole->CoM[1];
        const double dz = c->multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
                              ? c->multipole->CoM[2] - c->loc[2]
                              : c->loc[2] + c->width[2] - c->multipole->CoM[2];
        c->multipole->r_max = sqrt(dx * dx + dy * dy + dz * dz);
      } else {
        gravity_multipole_init(&c->multipole->m_pole);
        c->multipole->CoM[0] = c->loc[0] + c->width[0] / 2.;
        c->multipole->CoM[1] = c->loc[1] + c->width[1] / 2.;
        c->multipole->CoM[2] = c->loc[2] + c->width[2] / 2.;
        c->multipole->r_max = 0.;
      }
      c->multipole->r_max_rebuild = c->multipole->r_max;
      c->multipole->CoM_rebuild[0] = c->multipole->CoM[0];
      c->multipole->CoM_rebuild[1] = c->multipole->CoM[1];
      c->multipole->CoM_rebuild[2] = c->multipole->CoM[2];
    }
  }

  /* Set the values for this cell. */
  c->h_max = h_max;
  c->ti_end_min = ti_end_min;
  c->ti_end_max = ti_end_max;
  c->ti_beg_max = ti_beg_max;
  c->maxdepth = maxdepth;

  /* Set ownership according to the start of the parts array. */
  if (s->nr_parts > 0)
    c->owner =
        ((c->parts - s->parts) % s->nr_parts) * s->nr_queues / s->nr_parts;
  else if (s->nr_sparts > 0)
    c->owner =
        ((c->sparts - s->sparts) % s->nr_sparts) * s->nr_queues / s->nr_sparts;
  else if (s->nr_gparts > 0)
    c->owner =
        ((c->gparts - s->gparts) % s->nr_gparts) * s->nr_queues / s->nr_gparts;
  else
    c->owner = 0; /* Ok, there is really nothing on this rank... */

  /* Clean up. */
  if (allocate_buffer) {
    if (buff != NULL) free(buff);
    if (gbuff != NULL) free(gbuff);
    if (sbuff != NULL) free(sbuff);
  }
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

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[ind];
    space_split_recursive(s, c, NULL, NULL, NULL);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    if (!checkCellhdxmax(&cells_top[ind], &depth))
      message("    at cell depth %d", depth);
  }
#endif
}

/**
 * @brief Return a used cell to the buffer of unused sub-cells.
 *
 * @param s The #space.
 * @param c The #cell.
 */
void space_recycle(struct space *s, struct cell *c) {

  /* Clear the cell. */
  if (lock_destroy(&c->lock) != 0 || lock_destroy(&c->glock) != 0 ||
      lock_destroy(&c->mlock) != 0 || lock_destroy(&c->slock) != 0)
    error("Failed to destroy spinlocks.");

  /* Clear this cell's sort arrays. */
  for (int i = 0; i < 13; i++)
    if (c->sort[i] != NULL) free(c->sort[i]);

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Hook the multipole back in the buffer */
  if (s->gravity) {
    c->multipole->next = s->multipoles_sub;
    s->multipoles_sub = c->multipole;
  }

  /* Hook this cell into the buffer. */
  c->next = s->cells_sub;
  s->cells_sub = c;
  s->tot_cells -= 1;

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Return a list of used cells to the buffer of unused sub-cells.
 *
 * @param s The #space.
 * @param cell_list_begin Pointer to the first #cell in the linked list of
 *        cells joined by their @c next pointers.
 * @param cell_list_end Pointer to the last #cell in the linked list of
 *        cells joined by their @c next pointers. It is assumed that this
 *        cell's @c next pointer is @c NULL.
 * @param multipole_list_begin Pointer to the first #multipole in the linked
 * list of
 *        multipoles joined by their @c next pointers.
 * @param multipole_list_end Pointer to the last #multipole in the linked list
 * of
 *        multipoles joined by their @c next pointers. It is assumed that this
 *        multipole's @c next pointer is @c NULL.
 */
void space_recycle_list(struct space *s, struct cell *cell_list_begin,
                        struct cell *cell_list_end,
                        struct gravity_tensors *multipole_list_begin,
                        struct gravity_tensors *multipole_list_end) {

  int count = 0;

  /* Clean up the list of cells. */
  for (struct cell *c = cell_list_begin; c != NULL; c = c->next) {
    /* Clear the cell. */
    if (lock_destroy(&c->lock) != 0 || lock_destroy(&c->glock) != 0 ||
        lock_destroy(&c->mlock) != 0 || lock_destroy(&c->slock) != 0)
      error("Failed to destroy spinlocks.");

    /* Clear this cell's sort arrays. */
    for (int i = 0; i < 13; i++)
      if (c->sort[i] != NULL) free(c->sort[i]);

    /* Count this cell. */
    count += 1;
  }

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Hook the cells into the buffer. */
  cell_list_end->next = s->cells_sub;
  s->cells_sub = cell_list_begin;
  s->tot_cells -= count;

  /* Hook the multipoles into the buffer. */
  if (s->gravity) {
    multipole_list_end->next = s->multipoles_sub;
    s->multipoles_sub = multipole_list_begin;
  }

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
 * @param nr_cells Number of #cell to pick up.
 * @param cells Array of @c nr_cells #cell pointers in which to store the
 *        new cells.
 */
void space_getcells(struct space *s, int nr_cells, struct cell **cells) {

  /* Lock the space. */
  lock_lock(&s->lock);

  /* For each requested cell... */
  for (int j = 0; j < nr_cells; j++) {

    /* Is the cell buffer empty? */
    if (s->cells_sub == NULL) {
      if (posix_memalign((void *)&s->cells_sub, cell_align,
                         space_cellallocchunk * sizeof(struct cell)) != 0)
        error("Failed to allocate more cells.");

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->cells_sub[k].next = &s->cells_sub[k + 1];
      s->cells_sub[space_cellallocchunk - 1].next = NULL;
    }

    /* Is the multipole buffer empty? */
    if (s->gravity && s->multipoles_sub == NULL) {
      if (posix_memalign(
              (void *)&s->multipoles_sub, multipole_align,
              space_cellallocchunk * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate more multipoles.");

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->multipoles_sub[k].next = &s->multipoles_sub[k + 1];
      s->multipoles_sub[space_cellallocchunk - 1].next = NULL;
    }

    /* Pick off the next cell. */
    cells[j] = s->cells_sub;
    s->cells_sub = cells[j]->next;
    s->tot_cells += 1;

    /* Hook the multipole */
    if (s->gravity) {
      cells[j]->multipole = s->multipoles_sub;
      s->multipoles_sub = cells[j]->multipole->next;
    }
  }

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);

  /* Init some things in the cell we just got. */
  for (int j = 0; j < nr_cells; j++) {
    struct gravity_tensors *temp = cells[j]->multipole;
    bzero(cells[j], sizeof(struct cell));
    cells[j]->multipole = temp;
    cells[j]->nodeID = -1;
    if (lock_init(&cells[j]->lock) != 0 || lock_init(&cells[j]->glock) != 0 ||
        lock_init(&cells[j]->mlock) != 0 || lock_init(&cells[j]->slock) != 0)
      error("Failed to initialize cell spinlocks.");
  }
}

void space_synchronize_particle_positions_mapper(void *map_data, int nr_gparts,
                                                 void *extra_data) {
  /* Unpack the data */
  struct gpart *restrict gparts = (struct gpart *)map_data;
  struct space *s = (struct space *)extra_data;

  for (int k = 0; k < nr_gparts; k++) {

    /* Get the particle */
    const struct gpart *restrict gp = &gparts[k];

    if (gp->type == swift_type_dark_matter)
      continue;

    else if (gp->type == swift_type_gas) {

      /* Get it's gassy friend */
      struct part *p = &s->parts[-gp->id_or_neg_offset];
      struct xpart *xp = &s->xparts[-gp->id_or_neg_offset];

      /* Synchronize positions and velocities */
      p->x[0] = gp->x[0];
      p->x[1] = gp->x[1];
      p->x[2] = gp->x[2];

      xp->v_full[0] = gp->v_full[0];
      xp->v_full[1] = gp->v_full[1];
      xp->v_full[2] = gp->v_full[2];
    }

    else if (gp->type == swift_type_star) {

      /* Get it's stellar friend */
      struct spart *sp = &s->sparts[-gp->id_or_neg_offset];

      /* Synchronize positions */
      sp->x[0] = gp->x[0];
      sp->x[1] = gp->x[1];
      sp->x[2] = gp->x[2];
    }
  }
}

void space_synchronize_particle_positions(struct space *s) {

  if ((s->nr_gparts > 0 && s->nr_parts > 0) ||
      (s->nr_gparts > 0 && s->nr_sparts > 0))
    threadpool_map(&s->e->threadpool,
                   space_synchronize_particle_positions_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart), 0, (void *)s);
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

#ifdef SWIFT_DEBUG_CHECKS
    p->ti_drift = 0;
    p->ti_kick = 0;
#endif
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

#ifdef SWIFT_DEBUG_CHECKS
    gp->ti_drift = 0;
    gp->ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the s-particles by setting them into a valid state
 *
 * Calls star_first_init_spart() on all the particles
 */
void space_init_sparts(struct space *s) {

  const size_t nr_sparts = s->nr_sparts;
  struct spart *restrict sp = s->sparts;

  for (size_t i = 0; i < nr_sparts; ++i) {

#ifdef HYDRO_DIMENSION_2D
    sp[i].x[2] = 0.f;
    sp[i].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    sp[i].x[1] = sp[i].x[2] = 0.f;
    sp[i].v[1] = sp[i].v[2] = 0.f;
#endif

    star_first_init_spart(&sp[i]);

#ifdef SWIFT_DEBUG_CHECKS
    sp->ti_drift = 0;
    sp->ti_kick = 0;
#endif
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
 * @param sparts Array of star particles.
 * @param Npart The number of Gas particles in the space.
 * @param Ngpart The number of Gravity particles in the space.
 * @param Nspart The number of star particles in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param replicate How many replications along each direction do we want ?
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
                struct spart *sparts, size_t Npart, size_t Ngpart,
                size_t Nspart, int periodic, int replicate, int gravity,
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
  s->nr_sparts = Nspart;
  s->size_sparts = Nspart;
  s->sparts = sparts;
  s->nr_queues = 1; /* Temporary value until engine construction */

  /* Are we replicating the space ? */
  if (replicate < 1)
    error("Value of 'InitialConditions:replicate' (%d) is too small",
          replicate);
  if (replicate > 1) {
    space_replicate(s, replicate, verbose);
    parts = s->parts;
    gparts = s->gparts;
    sparts = s->sparts;
    Npart = s->nr_parts;
    Ngpart = s->nr_gparts;
    Nspart = s->nr_sparts;
  }

  /* Decide on the minimal top-level cell size */
  const double dmax = max(max(s->dim[0], s->dim[1]), s->dim[2]);
  int maxtcells =
      parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                               space_max_top_level_cells_default);
  s->cell_min = 0.99 * dmax / maxtcells;

  /* Check that it is big enough. */
  const double dmin = min(min(s->dim[0], s->dim[1]), s->dim[2]);
  int needtcells = 3 * dmax / dmin;
  if (maxtcells < needtcells)
    error(
        "Scheduler:max_top_level_cells is too small %d, needs to be at "
        "least %d",
        maxtcells, needtcells);

  /* Get the constants for the scheduler */
  space_maxsize = parser_get_opt_param_int(params, "Scheduler:cell_max_size",
                                           space_maxsize_default);
  space_subsize_pair = parser_get_opt_param_int(
      params, "Scheduler:cell_sub_size_pair", space_subsize_pair_default);
  space_subsize_self = parser_get_opt_param_int(
      params, "Scheduler:cell_sub_size_self", space_subsize_self_default);
  space_subsize_self_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_grav",
                               space_subsize_self_grav_default);
  space_splitsize = parser_get_opt_param_int(
      params, "Scheduler:cell_split_size", space_splitsize_default);

  if (verbose)
    message(
        "max_size set to %d, sub_size_pair set to %d, sub_size_self set to %d, "
        "split_size set to %d",
        space_maxsize, space_subsize_pair, space_subsize_self, space_splitsize);

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
    for (size_t k = 0; k < Nspart; k++) {
      sparts[k].x[0] += shift[0];
      sparts[k].x[1] += shift[1];
      sparts[k].x[2] += shift[2];
    }
  }

  if (!dry_run) {

    /* Check that all the part positions are reasonable, wrap if periodic. */
    if (periodic) {
      for (size_t k = 0; k < Npart; k++)
        for (int j = 0; j < 3; j++) {
          while (parts[k].x[j] < 0) parts[k].x[j] += s->dim[j];
          while (parts[k].x[j] >= s->dim[j]) parts[k].x[j] -= s->dim[j];
        }
    } else {
      for (size_t k = 0; k < Npart; k++)
        for (int j = 0; j < 3; j++)
          if (parts[k].x[j] < 0 || parts[k].x[j] >= s->dim[j])
            error("Not all particles are within the specified domain.");
    }

    /* Same for the gparts */
    if (periodic) {
      for (size_t k = 0; k < Ngpart; k++)
        for (int j = 0; j < 3; j++) {
          while (gparts[k].x[j] < 0) gparts[k].x[j] += s->dim[j];
          while (gparts[k].x[j] >= s->dim[j]) gparts[k].x[j] -= s->dim[j];
        }
    } else {
      for (size_t k = 0; k < Ngpart; k++)
        for (int j = 0; j < 3; j++)
          if (gparts[k].x[j] < 0 || gparts[k].x[j] >= s->dim[j])
            error("Not all g-particles are within the specified domain.");
    }

    /* Same for the sparts */
    if (periodic) {
      for (size_t k = 0; k < Nspart; k++)
        for (int j = 0; j < 3; j++) {
          while (sparts[k].x[j] < 0) sparts[k].x[j] += s->dim[j];
          while (sparts[k].x[j] >= s->dim[j]) sparts[k].x[j] -= s->dim[j];
        }
    } else {
      for (size_t k = 0; k < Nspart; k++)
        for (int j = 0; j < 3; j++)
          if (sparts[k].x[j] < 0 || sparts[k].x[j] >= s->dim[j])
            error("Not all s-particles are within the specified domain.");
    }
  }

  /* Allocate the extra parts array for the gas particles. */
  if (Npart > 0) {
    if (posix_memalign((void *)&s->xparts, xpart_align,
                       Npart * sizeof(struct xpart)) != 0)
      error("Failed to allocate xparts.");
    bzero(s->xparts, Npart * sizeof(struct xpart));
  }

  hydro_space_init(&s->hs, s);

  /* Set the particles in a state where they are ready for a run */
  space_init_parts(s);
  space_init_xparts(s);
  space_init_gparts(s);
  space_init_sparts(s);

  /* Init the space lock. */
  if (lock_init(&s->lock) != 0) error("Failed to create space spin-lock.");

  /* Build the cells recursively. */
  if (!dry_run) space_regrid(s, verbose);
}

/**
 * @brief Replicate the content of a space along each axis.
 *
 * Should only be called during initialisation.
 *
 * @param s The #space to replicate.
 * @param replicate The number of copies along each axis.
 * @param verbose Are we talkative ?
 */
void space_replicate(struct space *s, int replicate, int verbose) {

  if (replicate < 1) error("Invalid replicate value: %d", replicate);

  message("Replicating space %d times along each axis.", replicate);

  const int factor = replicate * replicate * replicate;

  /* Store the current values */
  const size_t nr_parts = s->nr_parts;
  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_sparts = s->nr_sparts;
  const size_t nr_dm = nr_gparts - nr_parts - nr_sparts;

  s->size_parts = s->nr_parts = nr_parts * factor;
  s->size_gparts = s->nr_gparts = nr_gparts * factor;
  s->size_sparts = s->nr_sparts = nr_sparts * factor;

  /* Allocate space for new particles */
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;

  if (posix_memalign((void *)&parts, part_align,
                     s->nr_parts * sizeof(struct part)) != 0)
    error("Failed to allocate new part array.");

  if (posix_memalign((void *)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate new gpart array.");

  if (posix_memalign((void *)&sparts, spart_align,
                     s->nr_sparts * sizeof(struct spart)) != 0)
    error("Failed to allocate new spart array.");

  /* Replicate everything */
  for (int i = 0; i < replicate; ++i) {
    for (int j = 0; j < replicate; ++j) {
      for (int k = 0; k < replicate; ++k) {
        const size_t offset = i * replicate * replicate + j * replicate + k;

        /* First copy the data */
        memcpy(parts + offset * nr_parts, s->parts,
               nr_parts * sizeof(struct part));
        memcpy(sparts + offset * nr_sparts, s->sparts,
               nr_sparts * sizeof(struct spart));
        memcpy(gparts + offset * nr_gparts, s->gparts,
               nr_gparts * sizeof(struct gpart));

        /* Shift the positions */
        const double shift[3] = {i * s->dim[0], j * s->dim[1], k * s->dim[2]};

        for (size_t n = offset * nr_parts; n < (offset + 1) * nr_parts; ++n) {
          parts[n].x[0] += shift[0];
          parts[n].x[1] += shift[1];
          parts[n].x[2] += shift[2];
        }
        for (size_t n = offset * nr_gparts; n < (offset + 1) * nr_gparts; ++n) {
          gparts[n].x[0] += shift[0];
          gparts[n].x[1] += shift[1];
          gparts[n].x[2] += shift[2];
        }
        for (size_t n = offset * nr_sparts; n < (offset + 1) * nr_sparts; ++n) {
          sparts[n].x[0] += shift[0];
          sparts[n].x[1] += shift[1];
          sparts[n].x[2] += shift[2];
        }

        /* Set the correct links (recall gpart are sorted by type at start-up):
           first DM (unassociated gpart), then gas, then stars */
        if (nr_parts > 0 && nr_gparts > 0) {
          const size_t offset_part = offset * nr_parts;
          const size_t offset_gpart = offset * nr_gparts + nr_dm;

          for (size_t n = 0; n < nr_parts; ++n) {
            parts[offset_part + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_part + n);
          }
        }
        if (nr_sparts > 0 && nr_gparts > 0) {
          const size_t offset_spart = offset * nr_sparts;
          const size_t offset_gpart = offset * nr_gparts + nr_dm + nr_parts;

          for (size_t n = 0; n < nr_sparts; ++n) {
            sparts[offset_spart + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_spart + n);
          }
        }
      }
    }
  }

  /* Replace the content of the space */
  free(s->parts);
  free(s->gparts);
  free(s->sparts);
  s->parts = parts;
  s->gparts = gparts;
  s->sparts = sparts;

  /* Finally, update the domain size */
  s->dim[0] *= replicate;
  s->dim[1] *= replicate;
  s->dim[2] *= replicate;

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that everything is correct */
  part_verify_links(s->parts, s->gparts, s->sparts, s->nr_parts, s->nr_gparts,
                    s->nr_sparts, verbose);
#endif
}

/**
 * @brief Cleans-up all the cell links in the space
 *
 * Expensive funtion. Should only be used for debugging purposes.
 *
 * @param s The #space to clean.
 */
void space_link_cleanup(struct space *s) {

  /* Recursively apply the cell link cleaning routine */
  space_map_cells_pre(s, 1, cell_clean_links, NULL);
}

/**
 * @brief Checks that all cells have been drifted to a given point in time
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to check.
 * @param ti_drift The (integer) time.
 * @param multipole Are we also checking the multipoles ?
 */
void space_check_drift_point(struct space *s, integertime_t ti_drift,
                             int multipole) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Recursively check all cells */
  space_map_cells_pre(s, 1, cell_check_part_drift_point, &ti_drift);
  space_map_cells_pre(s, 1, cell_check_gpart_drift_point, &ti_drift);
  if (multipole)
    space_map_cells_pre(s, 1, cell_check_multipole_drift_point, &ti_drift);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

void space_check_top_multipoles_drift_point(struct space *s,
                                            integertime_t ti_drift) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_check_multipole_drift_point(&s->cells_top[i], &ti_drift);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that all particles and local cells have a non-zero time-step.
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to check.
 */
void space_check_timesteps(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_check_timesteps(&s->cells_top[i]);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Resets all the individual cell task counters to 0.
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to reset.
 */
void space_reset_task_counters(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_reset_task_counters(&s->cells_top[i]);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Frees up the memory allocated for this #space
 */
void space_clean(struct space *s) {

  for (int i = 0; i < s->nr_cells; ++i) cell_clean(&s->cells_top[i]);
  free(s->cells_top);
  free(s->multipoles_top);
  free(s->parts);
  free(s->xparts);
  free(s->gparts);
  free(s->sparts);
}
