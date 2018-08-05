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
#include "chemistry.h"
#include "const.h"
#include "cooling.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "lock.h"
#include "memswap.h"
#include "minmax.h"
#include "multipole.h"
#include "restart.h"
#include "sort_part.h"
#include "stars.h"
#include "threadpool.h"
#include "tools.h"

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize_pair_hydro = space_subsize_pair_hydro_default;
int space_subsize_self_hydro = space_subsize_self_hydro_default;
int space_subsize_pair_grav = space_subsize_pair_grav_default;
int space_subsize_self_grav = space_subsize_self_grav_default;
int space_subdepth_grav = space_subdepth_grav_default;
int space_maxsize = space_maxsize_default;
#ifdef SWIFT_DEBUG_CHECKS
int last_cell_id;
#endif

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
  int *cell_counts;
};

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
    c->dx_max_sort = 0.0f;
    c->sorted = 0;
    c->count = 0;
    c->gcount = 0;
    c->scount = 0;
    c->init_grav = NULL;
    c->init_grav_out = NULL;
    c->extra_ghost = NULL;
    c->ghost_in = NULL;
    c->ghost_out = NULL;
    c->ghost = NULL;
    c->kick1 = NULL;
    c->kick2 = NULL;
    c->timestep = NULL;
    c->end_force = NULL;
    c->drift_part = NULL;
    c->drift_gpart = NULL;
    c->cooling = NULL;
    c->sourceterms = NULL;
    c->grav_long_range = NULL;
    c->grav_down_in = NULL;
    c->grav_down = NULL;
    c->grav_mesh = NULL;
    c->super = c;
    c->super_hydro = c;
    c->super_gravity = c;
    c->parts = NULL;
    c->xparts = NULL;
    c->gparts = NULL;
    c->sparts = NULL;
    c->do_sub_sort = 0;
    c->do_grav_sub_drift = 0;
    c->do_sub_drift = 0;
    if (s->gravity) bzero(c->multipole, sizeof(struct gravity_tensors));
    for (int i = 0; i < 13; i++)
      if (c->sort[i] != NULL) {
        free(c->sort[i]);
        c->sort[i] = NULL;
      }
#if WITH_MPI
    c->recv_xv = NULL;
    c->recv_rho = NULL;
    c->recv_gradient = NULL;
    c->recv_grav = NULL;
    c->recv_ti = NULL;

    c->send_xv = NULL;
    c->send_rho = NULL;
    c->send_gradient = NULL;
    c->send_grav = NULL;
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
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

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
      (int)floor(s->dim[0] /
                 fmax(h_max * kernel_gamma * space_stretch, s->cell_min)),
      (int)floor(s->dim[1] /
                 fmax(h_max * kernel_gamma * space_stretch, s->cell_min)),
      (int)floor(s->dim[2] /
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

  /* Are we about to allocate new top level cells without a regrid?
   * Can happen when restarting the application. */
  int no_regrid = (s->cells_top == NULL && oldnodeIDs == NULL);
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
      free(s->local_cells_top);
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
    const float dmin = min3(s->width[0], s->width[1], s->width[2]);

    /* Allocate the highest level of cells. */
    s->tot_cells = s->nr_cells = cdim[0] * cdim[1] * cdim[2];
    if (posix_memalign((void **)&s->cells_top, cell_align,
                       s->nr_cells * sizeof(struct cell)) != 0)
      error("Failed to allocate top-level cells.");
    bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

    /* Allocate the multipoles for the top-level cells. */
    if (s->gravity) {
      if (posix_memalign((void **)&s->multipoles_top, multipole_align,
                         s->nr_cells * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate top-level multipoles.");
      bzero(s->multipoles_top, s->nr_cells * sizeof(struct gravity_tensors));
    }

    /* Allocate the indices of local cells */
    if (posix_memalign((void **)&s->local_cells_top, SWIFT_STRUCT_ALIGNMENT,
                       s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells.");
    bzero(s->local_cells_top, s->nr_cells * sizeof(int));

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
          c->super_hydro = c;
          c->super_gravity = c;
          c->ti_old_part = ti_current;
          c->ti_old_gpart = ti_current;
          c->ti_old_multipole = ti_current;
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

    } else if (no_regrid && s->e != NULL) {
      /* If we have created the top-levels cells and not done an initial
       * partition (can happen when restarting), then the top-level cells
       * are not assigned to a node, we must do that and then associate the
       * particles with the cells. Note requires that
       * partition_store_celllist() was called once before, or just before
       * dumping the restart files.*/
      partition_restore_celllist(s, s->e->reparttype);

      /* Now re-distribute the particles, should just add to cells? */
      engine_redistribute(s->e);

      /* Make the proxies. */
      engine_makeproxies(s->e);
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
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the particles and get their cell index. Allocates
     an index that is larger than the number of particles to avoid
     re-allocating after shuffling. */
  const size_t ind_size = s->size_parts + 100;
  int *ind = (int *)malloc(sizeof(int) * ind_size);
  if (ind == NULL) error("Failed to allocate temporary particle indices.");
  int *cell_part_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_part_counts == NULL)
    error("Failed to allocate cell part count buffer.");
  if (s->size_parts > 0)
    space_parts_get_cell_index(s, ind, cell_part_counts, cells_top, verbose);

  /* Run through the gravity particles and get their cell index. */
  const size_t gind_size = s->size_gparts + 100;
  int *gind = (int *)malloc(sizeof(int) * gind_size);
  if (gind == NULL) error("Failed to allocate temporary g-particle indices.");
  int *cell_gpart_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_gpart_counts == NULL)
    error("Failed to allocate cell gpart count buffer.");
  if (s->size_gparts > 0)
    space_gparts_get_cell_index(s, gind, cell_gpart_counts, cells_top, verbose);

  /* Run through the star particles and get their cell index. */
  const size_t sind_size = s->size_sparts + 100;
  int *sind = (int *)malloc(sizeof(int) * sind_size);
  if (sind == NULL) error("Failed to allocate temporary s-particle indices.");
  int *cell_spart_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_spart_counts == NULL)
    error("Failed to allocate cell gpart count buffer.");
  if (s->size_sparts > 0)
    space_sparts_get_cell_index(s, sind, cell_spart_counts, cells_top, verbose);

#ifdef WITH_MPI
  const int local_nodeID = s->e->nodeID;

  /* Move non-local parts to the end of the list. */
  for (size_t k = 0; k < nr_parts;) {
    if (cells_top[ind[k]].nodeID != local_nodeID) {
      nr_parts -= 1;
      /* Swap the particle */
      memswap(&s->parts[k], &s->parts[nr_parts], sizeof(struct part));
      /* Swap the link with the gpart */
      if (s->parts[k].gpart != NULL) {
        s->parts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->parts[nr_parts].gpart != NULL) {
        s->parts[nr_parts].gpart->id_or_neg_offset = -nr_parts;
      }
      /* Swap the xpart */
      memswap(&s->xparts[k], &s->xparts[nr_parts], sizeof(struct xpart));
      /* Swap the index */
      memswap(&ind[k], &ind[nr_parts], sizeof(int));
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
      memswap(&s->sparts[k], &s->sparts[nr_sparts], sizeof(struct spart));
      /* Swap the link with the gpart */
      if (s->sparts[k].gpart != NULL) {
        s->sparts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->sparts[nr_sparts].gpart != NULL) {
        s->sparts[nr_sparts].gpart->id_or_neg_offset = -nr_sparts;
      }
      /* Swap the index */
      memswap(&sind[k], &sind[nr_sparts], sizeof(int));
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
      memswap(&s->gparts[k], &s->gparts[nr_gparts], sizeof(struct gpart));
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
      memswap(&gind[k], &gind[nr_gparts], sizeof(int));
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

  /* Clear non-local cell counts. */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID != local_nodeID) {
      cell_part_counts[k] = 0;
      cell_spart_counts[k] = 0;
      cell_gpart_counts[k] = 0;
    }
  }

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
    cell_part_counts[ind[k]]++;
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
    cell_spart_counts[sind[k]]++;
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
    space_parts_sort(s->parts, s->xparts, ind, cell_part_counts, s->nr_cells,
                     0);

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
    space_sparts_sort(s->sparts, sind, cell_spart_counts, s->nr_cells, 0);

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
  free(cell_part_counts);
  free(sind);
  free(cell_spart_counts);

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
    cell_gpart_counts[gind[k]]++;
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
    space_gparts_sort(s->gparts, s->parts, s->sparts, gind, cell_gpart_counts,
                      s->nr_cells);

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
  free(cell_gpart_counts);

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
    c->ti_old_part = ti_current;
    c->ti_old_gpart = ti_current;
    c->ti_old_multipole = ti_current;
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

  /* Clean up any stray sort indices in the cell buffer. */
  space_free_buff_sort_indices(s);

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

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;

  /* Loop over the parts. */
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

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    ind[k] = index;
    cell_counts[index]++;

    /* Compute minimal mass */
    min_mass = min(min_mass, hydro_get_mass(p));

    /* Compute sum of velocity norm */
    sum_vel_norm += p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2];

    /* Update the position */
    p->x[0] = pos_x;
    p->x[1] = pos_y;
    p->x[2] = pos_z;
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_part_mass, min_mass);
  atomic_add_f(&s->sum_part_vel_norm, sum_vel_norm);
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

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;

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

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    ind[k] = index;
    cell_counts[index]++;

    /* Compute minimal mass */
    if (gp->type == swift_type_dark_matter) {
      min_mass = min(min_mass, gp->mass);
      sum_vel_norm += gp->v_full[0] * gp->v_full[0] +
                      gp->v_full[1] * gp->v_full[1] +
                      gp->v_full[2] * gp->v_full[2];
    }

    /* Update the position */
    gp->x[0] = pos_x;
    gp->x[1] = pos_y;
    gp->x[2] = pos_z;
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_gpart_mass, min_mass);
  atomic_add_f(&s->sum_gpart_vel_norm, sum_vel_norm);
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

  /* Init the local count buffer. */
  int *cell_counts = (int *)calloc(sizeof(int), s->nr_cells);
  if (cell_counts == NULL)
    error("Failed to allocate temporary cell count buffer.");

  /* Init the local collectors */
  float min_mass = FLT_MAX;
  float sum_vel_norm = 0.f;

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

#ifdef SWIFT_DEBUG_CHECKS
    if (index < 0 || index >= cdim[0] * cdim[1] * cdim[2])
      error("Invalid index=%d cdim=[%d %d %d] p->x=[%e %e %e]", index, cdim[0],
            cdim[1], cdim[2], pos_x, pos_y, pos_z);

    if (pos_x >= dim_x || pos_y >= dim_y || pos_z >= dim_z || pos_x < 0. ||
        pos_y < 0. || pos_z < 0.)
      error("Particle outside of simulation box. p->x=[%e %e %e]", pos_x, pos_y,
            pos_z);
#endif

    ind[k] = index;
    cell_counts[index]++;

    /* Compute minimal mass */
    min_mass = min(min_mass, sp->mass);

    /* Compute sum of velocity norm */
    sum_vel_norm +=
        sp->v[0] * sp->v[0] + sp->v[1] * sp->v[1] + sp->v[2] * sp->v[2];

    /* Update the position */
    sp->x[0] = pos_x;
    sp->x[1] = pos_y;
    sp->x[2] = pos_z;
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write back the minimal part mass and velocity sum */
  atomic_min_f(&s->min_spart_mass, min_mass);
  atomic_add_f(&s->sum_spart_vel_norm, sum_vel_norm);
}

/**
 * @brief Computes the cell index of all the particles.
 *
 * Also computes the minimal mass of all #part.
 *
 * @param s The #space.
 * @param ind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_parts_get_cell_index(struct space *s, int *ind, int *cell_counts,
                                struct cell *cells, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_part_mass = FLT_MAX;
  s->sum_part_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = ind;
  data.cell_counts = cell_counts;

  threadpool_map(&s->e->threadpool, space_parts_get_cell_index_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), 0, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the g-particles.
 *
 * Also computes the minimal mass of all dark-matter #gpart.
 *
 * @param s The #space.
 * @param gind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_gparts_get_cell_index(struct space *s, int *gind, int *cell_counts,
                                 struct cell *cells, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_gpart_mass = FLT_MAX;
  s->sum_gpart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = gind;
  data.cell_counts = cell_counts;

  threadpool_map(&s->e->threadpool, space_gparts_get_cell_index_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart), 0, &data);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Computes the cell index of all the s-particles.
 *
 * Also computes the minimal mass of all #spart.
 *
 * @param s The #space.
 * @param sind The array of indices to fill.
 * @param cell_counts The cell counters to update.
 * @param cells The array of #cell to update.
 * @param verbose Are we talkative ?
 */
void space_sparts_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                 struct cell *cells, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_spart_mass = FLT_MAX;
  s->sum_spart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.cells = cells;
  data.ind = sind;
  data.cell_counts = cell_counts;

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
 * @param parts The array of #part to sort.
 * @param xparts The corresponding #xpart array to sort as well.
 * @param ind The indices with respect to which the parts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of count).
 * @param parts_offset Offset of the #part array from the global #part array.
 */
void space_parts_sort(struct part *parts, struct xpart *xparts, int *ind,
                      int *counts, int num_bins, ptrdiff_t parts_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (posix_memalign((void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct part temp_part = parts[k];
      struct xpart temp_xpart = xparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&parts[j], &temp_part, sizeof(struct part));
        memswap(&xparts[j], &temp_xpart, sizeof(struct xpart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (parts[j].gpart)
          parts[j].gpart->id_or_neg_offset = -(j + parts_offset);
      }
      parts[k] = temp_part;
      xparts[k] = temp_xpart;
      ind[k] = target_cid;
      if (parts[k].gpart)
        parts[k].gpart->id_or_neg_offset = -(k + parts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  free(offsets);
}

/**
 * @brief Sort the s-particles according to the given indices.
 *
 * @param sparts The array of #spart to sort.
 * @param ind The indices with respect to which the #spart are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param sparts_offset Offset of the #spart array from the global #spart.
 * array.
 */
void space_sparts_sort(struct spart *sparts, int *ind, int *counts,
                       int num_bins, ptrdiff_t sparts_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (posix_memalign((void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct spart temp_spart = sparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&sparts[j], &temp_spart, sizeof(struct spart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (sparts[j].gpart)
          sparts[j].gpart->id_or_neg_offset = -(j + sparts_offset);
      }
      sparts[k] = temp_spart;
      ind[k] = target_cid;
      if (sparts[k].gpart)
        sparts[k].gpart->id_or_neg_offset = -(k + sparts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  free(offsets);
}

/**
 * @brief Sort the g-particles according to the given indices.
 *
 * @param gparts The array of #gpart to sort.
 * @param parts Global #part array for re-linking.
 * @param sparts Global #spart array for re-linking.
 * @param ind The indices with respect to which the gparts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 */
void space_gparts_sort(struct gpart *gparts, struct part *parts,
                       struct spart *sparts, int *ind, int *counts,
                       int num_bins) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (posix_memalign((void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct gpart temp_gpart = gparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&gparts[j], &temp_gpart, sizeof(struct gpart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (gparts[j].type == swift_type_gas) {
          parts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        } else if (gparts[j].type == swift_type_star) {
          sparts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        }
      }
      gparts[k] = temp_gpart;
      ind[k] = target_cid;
      if (gparts[k].type == swift_type_gas) {
        parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_star) {
        sparts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  free(offsets);
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
  const int with_gravity = s->gravity;
  const int depth = c->depth;
  int maxdepth = 0;
  float h_max = 0.0f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  struct part *parts = c->parts;
  struct gpart *gparts = c->gparts;
  struct spart *sparts = c->sparts;
  struct xpart *xparts = c->xparts;
  struct engine *e = s->e;
  const integertime_t ti_current = e->ti_current;

  /* If the buff is NULL, allocate it, and remember to free it. */
  const int allocate_buffer = (buff == NULL && gbuff == NULL && sbuff == NULL);
  if (allocate_buffer) {
    if (count > 0) {
      if (posix_memalign((void **)&buff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * count) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < count; k++) {
        buff[k].x[0] = parts[k].x[0];
        buff[k].x[1] = parts[k].x[1];
        buff[k].x[2] = parts[k].x[2];
      }
    }
    if (gcount > 0) {
      if (posix_memalign((void **)&gbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * gcount) != 0)
        error("Failed to allocate temporary indices.");
      for (int k = 0; k < gcount; k++) {
        gbuff[k].x[0] = gparts[k].x[0];
        gbuff[k].x[1] = gparts[k].x[1];
        gbuff[k].x[2] = gparts[k].x[2];
      }
    }
    if (scount > 0) {
      if (posix_memalign((void **)&sbuff, SWIFT_STRUCT_ALIGNMENT,
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
  if ((with_gravity && gcount > space_splitsize) ||
      (!with_gravity &&
       (count > space_splitsize || scount > space_splitsize))) {

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
      cp->dx_max_sort = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->super = NULL;
      cp->super_hydro = NULL;
      cp->super_gravity = NULL;
      cp->do_sub_sort = 0;
      cp->do_grav_sub_drift = 0;
      cp->do_sub_drift = 0;
#ifdef SWIFT_DEBUG_CHECKS
      cp->cellID = last_cell_id++;
#endif
    }

    /* Split the cell's partcle data. */
    cell_split(c, c->parts - s->parts, c->sparts - s->sparts, buff, sbuff,
               gbuff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff;

    for (int k = 0; k < 8; k++) {

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->count == 0 && cp->gcount == 0 && cp->scount == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {

        /* Recurse */
        space_split_recursive(s, cp, progeny_buff, progeny_sbuff,
                              progeny_gbuff);

        /* Update the pointers in the buffers */
        progeny_buff += cp->count;
        progeny_gbuff += cp->gcount;
        progeny_sbuff += cp->scount;

        /* Update the cell-wide properties */
        h_max = max(h_max, cp->h_max);
        ti_hydro_end_min = min(ti_hydro_end_min, cp->ti_hydro_end_min);
        ti_hydro_end_max = max(ti_hydro_end_max, cp->ti_hydro_end_max);
        ti_hydro_beg_max = max(ti_hydro_beg_max, cp->ti_hydro_beg_max);
        ti_gravity_end_min = min(ti_gravity_end_min, cp->ti_gravity_end_min);
        ti_gravity_end_max = max(ti_gravity_end_max, cp->ti_gravity_end_max);
        ti_gravity_beg_max = max(ti_gravity_beg_max, cp->ti_gravity_beg_max);

        /* Increase the depth */
        if (cp->maxdepth > maxdepth) maxdepth = cp->maxdepth;
      }
    }

    /* Deal with the multipole */
    if (s->gravity) {

      /* Reset everything */
      gravity_reset(c->multipole);

      /* Compute CoM and bulk velocity from all progenies */
      double CoM[3] = {0., 0., 0.};
      double vel[3] = {0., 0., 0.};
      float max_delta_vel[3] = {0.f, 0.f, 0.f};
      float min_delta_vel[3] = {0.f, 0.f, 0.f};
      double mass = 0.;

      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct gravity_tensors *m = c->progeny[k]->multipole;

          mass += m->m_pole.M_000;

          CoM[0] += m->CoM[0] * m->m_pole.M_000;
          CoM[1] += m->CoM[1] * m->m_pole.M_000;
          CoM[2] += m->CoM[2] * m->m_pole.M_000;

          vel[0] += m->m_pole.vel[0] * m->m_pole.M_000;
          vel[1] += m->m_pole.vel[1] * m->m_pole.M_000;
          vel[2] += m->m_pole.vel[2] * m->m_pole.M_000;

          max_delta_vel[0] = max(m->m_pole.max_delta_vel[0], max_delta_vel[0]);
          max_delta_vel[1] = max(m->m_pole.max_delta_vel[1], max_delta_vel[1]);
          max_delta_vel[2] = max(m->m_pole.max_delta_vel[2], max_delta_vel[2]);

          min_delta_vel[0] = min(m->m_pole.min_delta_vel[0], min_delta_vel[0]);
          min_delta_vel[1] = min(m->m_pole.min_delta_vel[1], min_delta_vel[1]);
          min_delta_vel[2] = min(m->m_pole.min_delta_vel[2], min_delta_vel[2]);
        }
      }

      /* Final operation on the CoM and bulk velocity */
      const double inv_mass = 1. / mass;
      c->multipole->CoM[0] = CoM[0] * inv_mass;
      c->multipole->CoM[1] = CoM[1] * inv_mass;
      c->multipole->CoM[2] = CoM[2] * inv_mass;
      c->multipole->m_pole.vel[0] = vel[0] * inv_mass;
      c->multipole->m_pole.vel[1] = vel[1] * inv_mass;
      c->multipole->m_pole.vel[2] = vel[2] * inv_mass;

      /* Min max velocity along each axis */
      c->multipole->m_pole.max_delta_vel[0] = max_delta_vel[0];
      c->multipole->m_pole.max_delta_vel[1] = max_delta_vel[1];
      c->multipole->m_pole.max_delta_vel[2] = max_delta_vel[2];
      c->multipole->m_pole.min_delta_vel[0] = min_delta_vel[0];
      c->multipole->m_pole.min_delta_vel[1] = min_delta_vel[1];
      c->multipole->m_pole.min_delta_vel[2] = min_delta_vel[2];

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

      /* Store the value at rebuild time */
      c->multipole->r_max_rebuild = c->multipole->r_max;
      c->multipole->CoM_rebuild[0] = c->multipole->CoM[0];
      c->multipole->CoM_rebuild[1] = c->multipole->CoM[1];
      c->multipole->CoM_rebuild[2] = c->multipole->CoM[2];

      /* We know the first-order multipole (dipole) is 0. */
      c->multipole->m_pole.M_100 = 0.f;
      c->multipole->m_pole.M_010 = 0.f;
      c->multipole->m_pole.M_001 = 0.f;

    } /* Deal with gravity */
  }   /* Split or let it be? */

  /* Otherwise, collect the data from the particles this cell. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
    maxdepth = c->depth;

    timebin_t hydro_time_bin_min = num_time_bins, hydro_time_bin_max = 0;
    timebin_t gravity_time_bin_min = num_time_bins, gravity_time_bin_max = 0;

    /* parts: Get dt_min/dt_max and h_max. */
    for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (parts[k].time_bin == time_bin_inhibited)
        error("Inhibited particle present in space_split()");
#endif
      hydro_time_bin_min = min(hydro_time_bin_min, parts[k].time_bin);
      hydro_time_bin_max = max(hydro_time_bin_max, parts[k].time_bin);
      h_max = max(h_max, parts[k].h);
    }

    /* xparts: Reset x_diff */
    for (int k = 0; k < count; k++) {
      xparts[k].x_diff[0] = 0.f;
      xparts[k].x_diff[1] = 0.f;
      xparts[k].x_diff[2] = 0.f;
    }

    /* gparts: Get dt_min/dt_max. */
    for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (gparts[k].time_bin == time_bin_inhibited)
        error("Inhibited g-particle present in space_split()");
#endif
      gravity_time_bin_min = min(gravity_time_bin_min, gparts[k].time_bin);
      gravity_time_bin_max = max(gravity_time_bin_max, gparts[k].time_bin);
    }

    /* sparts: Get dt_min/dt_max */
    for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[k].time_bin == time_bin_inhibited)
        error("Inhibited s-particle present in space_split()");
#endif
      gravity_time_bin_min = min(gravity_time_bin_min, sparts[k].time_bin);
      gravity_time_bin_max = max(gravity_time_bin_max, sparts[k].time_bin);
    }

    /* Convert into integer times */
    ti_hydro_end_min = get_integer_time_end(ti_current, hydro_time_bin_min);
    ti_hydro_end_max = get_integer_time_end(ti_current, hydro_time_bin_max);
    ti_hydro_beg_max =
        get_integer_time_begin(ti_current + 1, hydro_time_bin_max);
    ti_gravity_end_min = get_integer_time_end(ti_current, gravity_time_bin_min);
    ti_gravity_end_max = get_integer_time_end(ti_current, gravity_time_bin_max);
    ti_gravity_beg_max =
        get_integer_time_begin(ti_current + 1, gravity_time_bin_max);

    /* Construct the multipole and the centre of mass*/
    if (s->gravity) {
      if (gcount > 0) {

        gravity_P2M(c->multipole, c->gparts, c->gcount);

      } else {

        /* No gparts in that leaf cell */

        /* Set the values to something sensible */
        gravity_multipole_init(&c->multipole->m_pole);
        if (c->nodeID == engine_rank) {
          c->multipole->CoM[0] = c->loc[0] + c->width[0] / 2.;
          c->multipole->CoM[1] = c->loc[1] + c->width[1] / 2.;
          c->multipole->CoM[2] = c->loc[2] + c->width[2] / 2.;
          c->multipole->r_max = 0.;
        }
      }

      /* Store the value at rebuild time */
      c->multipole->r_max_rebuild = c->multipole->r_max;
      c->multipole->CoM_rebuild[0] = c->multipole->CoM[0];
      c->multipole->CoM_rebuild[1] = c->multipole->CoM[1];
      c->multipole->CoM_rebuild[2] = c->multipole->CoM[2];
    }
  }

  /* Set the values for this cell. */
  c->h_max = h_max;
  c->ti_hydro_end_min = ti_hydro_end_min;
  c->ti_hydro_end_max = ti_hydro_end_max;
  c->ti_hydro_beg_max = ti_hydro_beg_max;
  c->ti_gravity_end_min = ti_gravity_end_min;
  c->ti_gravity_end_max = ti_gravity_end_max;
  c->ti_gravity_beg_max = ti_gravity_beg_max;
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
      if (posix_memalign((void **)&s->cells_sub, cell_align,
                         space_cellallocchunk * sizeof(struct cell)) != 0)
        error("Failed to allocate more cells.");

      /* Clear the newly-allocated cells. */
      bzero(s->cells_sub, sizeof(struct cell) * space_cellallocchunk);

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->cells_sub[k].next = &s->cells_sub[k + 1];
      s->cells_sub[space_cellallocchunk - 1].next = NULL;
    }

    /* Is the multipole buffer empty? */
    if (s->gravity && s->multipoles_sub == NULL) {
      if (posix_memalign(
              (void **)&s->multipoles_sub, multipole_align,
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
    for (int k = 0; k < 13; k++)
      if (cells[j]->sort[k] != NULL) free(cells[j]->sort[k]);
    struct gravity_tensors *temp = cells[j]->multipole;
    bzero(cells[j], sizeof(struct cell));
    cells[j]->multipole = temp;
    cells[j]->nodeID = -1;
    if (lock_init(&cells[j]->lock) != 0 || lock_init(&cells[j]->glock) != 0 ||
        lock_init(&cells[j]->mlock) != 0 || lock_init(&cells[j]->slock) != 0)
      error("Failed to initialize cell spinlocks.");
  }
}

/**
 * @brief Free sort arrays in any cells in the cell buffer.
 *
 * @param s The #space.
 */
void space_free_buff_sort_indices(struct space *s) {
  for (struct cell *finger = s->cells_sub; finger != NULL;
       finger = finger->next) {
    for (int k = 0; k < 13; k++)
      if (finger->sort[k] != NULL) {
        free(finger->sort[k]);
        finger->sort[k] = NULL;
      }
  }
}

/**
 * @brief Construct the list of top-level cells that have any tasks in
 * their hierarchy.
 *
 * This assumes the list has been pre-allocated at a regrid.
 *
 * @param s The #space.
 */
void space_list_cells_with_tasks(struct space *s) {

  /* Let's rebuild the list of local top-level cells */
  s->nr_local_cells = 0;
  for (int i = 0; i < s->nr_cells; ++i)
    if (cell_has_tasks(&s->cells_top[i])) {
      s->local_cells_top[s->nr_local_cells] = i;
      s->nr_local_cells++;
    }
  if (s->e->verbose)
    message("Have %d local cells (total=%d)", s->nr_local_cells, s->nr_cells);
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

void space_first_init_parts_mapper(void *restrict map_data, int count,
                                   void *restrict extra_data) {

  struct part *restrict p = (struct part *)map_data;
  const struct space *restrict s = (struct space *)extra_data;
  const struct engine *e = s->e;

  const ptrdiff_t delta = p - s->parts;
  struct xpart *restrict xp = s->xparts + delta;

  /* Extract some constants */
  const struct cosmology *cosmo = s->e->cosmology;
  const struct phys_const *phys_const = s->e->physical_constants;
  const struct unit_system *us = s->e->internal_units;
  const float a_factor_vel = cosmo->a;

  const struct hydro_props *hydro_props = s->e->hydro_properties;
  const float u_init = hydro_props->initial_internal_energy;
  const float u_min = hydro_props->minimal_internal_energy;

  const struct chemistry_global_data *chemistry = e->chemistry;
  const struct cooling_function_data *cool_func = e->cooling_func;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {
    p[k].v[0] *= a_factor_vel;
    p[k].v[1] *= a_factor_vel;
    p[k].v[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    p[k].x[2] = 0.f;
    p[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    p[k].x[1] = p[k].x[2] = 0.f;
    p[k].v[1] = p[k].v[2] = 0.f;
#endif
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    hydro_first_init_part(&p[k], &xp[k]);

    /* Overwrite the internal energy? */
    if (u_init > 0.f) hydro_set_init_internal_energy(&p[k], u_init);
    if (u_min > 0.f) hydro_set_init_internal_energy(&p[k], u_min);

    /* Also initialise the chemistry */
    chemistry_first_init_part(phys_const, us, cosmo, chemistry, &p[k], &xp[k]);

    /* And the cooling */
    cooling_first_init_part(phys_const, us, cosmo, cool_func, &p[k], &xp[k]);

#ifdef SWIFT_DEBUG_CHECKS
    /* Check part->gpart->part linkeage. */
    if (p[k].gpart && p[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> part link");

    /* Initialise the time-integration check variables */
    p[k].ti_drift = 0;
    p[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the particles by setting them into a valid state
 *
 * Calls hydro_first_init_part() on all the particles
 * Calls chemistry_first_init_part() on all the particles
 * Calls cooling_first_init_part() on all the particles
 */
void space_first_init_parts(struct space *s, int verbose) {

  const ticks tic = getticks();
  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_parts_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_gparts_mapper(void *restrict map_data, int count,
                                    void *restrict extra_data) {

  struct gpart *restrict gp = (struct gpart *)map_data;
  const struct space *restrict s = (struct space *)extra_data;

  const struct cosmology *cosmo = s->e->cosmology;
  const float a_factor_vel = cosmo->a;
  const struct gravity_props *grav_props = s->e->gravity_properties;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {
    gp[k].v_full[0] *= a_factor_vel;
    gp[k].v_full[1] *= a_factor_vel;
    gp[k].v_full[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    gp[k].x[2] = 0.f;
    gp[k].v_full[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    gp[k].x[1] = gp[k].x[2] = 0.f;
    gp[k].v_full[1] = gp[k].v_full[2] = 0.f;
#endif
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    gravity_first_init_gpart(&gp[k], grav_props);

#ifdef SWIFT_DEBUG_CHECKS
    /* Initialise the time-integration check variables */
    gp[k].ti_drift = 0;
    gp[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the g-particles by setting them into a valid state
 *
 * Calls gravity_first_init_gpart() on all the particles
 */
void space_first_init_gparts(struct space *s, int verbose) {

  const ticks tic = getticks();
  if (s->nr_gparts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_gparts_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_sparts_mapper(void *restrict map_data, int count,
                                    void *restrict extra_data) {

  struct spart *restrict sp = (struct spart *)map_data;
  const struct space *restrict s = (struct space *)extra_data;

#ifdef SWIFT_DEBUG_CHECKS
  const ptrdiff_t delta = sp - s->sparts;
#endif

  const struct cosmology *cosmo = s->e->cosmology;
  const float a_factor_vel = cosmo->a;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {

    sp[k].v[0] *= a_factor_vel;
    sp[k].v[1] *= a_factor_vel;
    sp[k].v[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    sp[k].x[2] = 0.f;
    sp[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    sp[k].x[1] = sp[k].x[2] = 0.f;
    sp[k].v[1] = sp[k].v[2] = 0.f;
#endif
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    star_first_init_spart(&sp[k]);

#ifdef SWIFT_DEBUG_CHECKS
    if (sp[k].gpart && sp[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> spart link");

    /* Initialise the time-integration check variables */
    sp[k].ti_drift = 0;
    sp[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the s-particles by setting them into a valid state
 *
 * Calls star_first_init_spart() on all the particles
 */
void space_first_init_sparts(struct space *s, int verbose) {
  const ticks tic = getticks();
  if (s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_sparts_mapper, s->sparts,
                   s->nr_sparts, sizeof(struct spart), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_parts_mapper(void *restrict map_data, int count,
                             void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct hydro_space *restrict hs = (struct hydro_space *)extra_data;
  for (int k = 0; k < count; k++) hydro_init_part(&parts[k], hs);
}

/**
 * @brief Calls the #part initialisation function on all particles in the space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_parts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_init_parts_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), 0, &s->hs);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_gparts_mapper(void *restrict map_data, int count,
                              void *restrict extra_data) {

  struct gpart *gparts = (struct gpart *)map_data;
  for (int k = 0; k < count; k++) gravity_init_gpart(&gparts[k]);
}

/**
 * @brief Calls the #gpart initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_gparts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_gparts > 0)
    threadpool_map(&s->e->threadpool, space_init_gparts_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart), 0, NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_convert_quantities_mapper(void *restrict map_data, int count,
                                     void *restrict extra_data) {
  struct space *s = (struct space *)extra_data;
  const struct cosmology *cosmo = s->e->cosmology;
  struct part *restrict parts = (struct part *)map_data;
  const ptrdiff_t index = parts - s->parts;
  struct xpart *restrict xparts = s->xparts + index;

  for (int k = 0; k < count; k++)
    hydro_convert_quantities(&parts[k], &xparts[k], cosmo);
}

/**
 * @brief Calls the #part quantities conversion function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_convert_quantities(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_convert_quantities_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param s The #space to initialize.
 * @param params The parsed parameter file.
 * @param cosmo The current cosmological model.
 * @param dim Spatial dimensions of the domain.
 * @param parts Array of Gas particles.
 * @param gparts Array of Gravity particles.
 * @param sparts Array of star particles.
 * @param Npart The number of Gas particles in the space.
 * @param Ngpart The number of Gravity particles in the space.
 * @param Nspart The number of star particles in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param replicate How many replications along each direction do we want?
 * @param generate_gas_in_ics Are we generating gas particles from the gparts?
 * @param self_gravity flag whether we are doing gravity or not?
 * @param verbose Print messages to stdout or not.
 * @param dry_run If 1, just initialise stuff, don't do anything with the parts.
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_splitsize
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */
void space_init(struct space *s, struct swift_params *params,
                const struct cosmology *cosmo, double dim[3],
                struct part *parts, struct gpart *gparts, struct spart *sparts,
                size_t Npart, size_t Ngpart, size_t Nspart, int periodic,
                int replicate, int generate_gas_in_ics, int self_gravity,
                int verbose, int dry_run) {

  /* Clean-up everything */
  bzero(s, sizeof(struct space));

  /* Store everything in the space. */
  s->dim[0] = dim[0];
  s->dim[1] = dim[1];
  s->dim[2] = dim[2];
  s->periodic = periodic;
  s->gravity = self_gravity;
  s->nr_parts = Npart;
  s->size_parts = Npart;
  s->parts = parts;
  s->nr_gparts = Ngpart;
  s->size_gparts = Ngpart;
  s->gparts = gparts;
  s->nr_sparts = Nspart;
  s->size_sparts = Nspart;
  s->sparts = sparts;
  s->min_part_mass = FLT_MAX;
  s->min_gpart_mass = FLT_MAX;
  s->min_spart_mass = FLT_MAX;
  s->sum_part_vel_norm = 0.f;
  s->sum_gpart_vel_norm = 0.f;
  s->sum_spart_vel_norm = 0.f;
  s->nr_queues = 1; /* Temporary value until engine construction */

  /* Are we generating gas from the DM-only ICs? */
  if (generate_gas_in_ics) {
    space_generate_gas(s, cosmo, verbose);
    parts = s->parts;
    gparts = s->gparts;
    Npart = s->nr_parts;
    Ngpart = s->nr_gparts;

#ifdef SWIFT_DEBUG_CHECKS
    part_verify_links(parts, gparts, sparts, Npart, Ngpart, Nspart, 1);
#endif
  }

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

#ifdef SWIFT_DEBUG_CHECKS
    part_verify_links(parts, gparts, sparts, Npart, Ngpart, Nspart, 1);
#endif
  }

  /* Decide on the minimal top-level cell size */
  const double dmax = max3(s->dim[0], s->dim[1], s->dim[2]);
  int maxtcells =
      parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                               space_max_top_level_cells_default);
  s->cell_min = 0.99 * dmax / maxtcells;

  /* Check that it is big enough. */
  const double dmin = min3(s->dim[0], s->dim[1], s->dim[2]);
  int needtcells = 3 * dmax / dmin;
  if (maxtcells < needtcells)
    error(
        "Scheduler:max_top_level_cells is too small %d, needs to be at "
        "least %d",
        maxtcells, needtcells);

  /* Get the constants for the scheduler */
  space_maxsize = parser_get_opt_param_int(params, "Scheduler:cell_max_size",
                                           space_maxsize_default);
  space_subsize_pair_hydro =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_pair_hydro",
                               space_subsize_pair_hydro_default);
  space_subsize_self_hydro =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_hydro",
                               space_subsize_self_hydro_default);
  space_subsize_pair_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_pair_grav",
                               space_subsize_pair_grav_default);
  space_subsize_self_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_grav",
                               space_subsize_self_grav_default);
  space_splitsize = parser_get_opt_param_int(
      params, "Scheduler:cell_split_size", space_splitsize_default);
  space_subdepth_grav = parser_get_opt_param_int(
      params, "Scheduler:cell_subdepth_grav", space_subdepth_grav_default);

  if (verbose) {
    message("max_size set to %d split_size set to %d", space_maxsize,
            space_splitsize);
    message("subdepth_grav set to %d", space_subdepth_grav);
    message("sub_size_pair_hydro set to %d, sub_size_self_hydro set to %d",
            space_subsize_pair_hydro, space_subsize_self_hydro);
    message("sub_size_pair_grav set to %d, sub_size_self_grav set to %d",
            space_subsize_pair_grav, space_subsize_self_grav);
  }

  /* Apply h scaling */
  const double scaling = parser_get_opt_param_double(
      params, "InitialConditions:smoothing_length_scaling", 1.0);
  if (scaling != 1.0 && !dry_run) {
    message("Re-scaling smoothing lengths by a factor %e", scaling);
    for (size_t k = 0; k < Npart; k++) parts[k].h *= scaling;
  }

  /* Apply shift */
  double shift[3] = {0.0, 0.0, 0.0};
  parser_get_opt_param_double_array(params, "InitialConditions:shift", 3,
                                    shift);
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
    if (posix_memalign((void **)&s->xparts, xpart_align,
                       Npart * sizeof(struct xpart)) != 0)
      error("Failed to allocate xparts.");
    bzero(s->xparts, Npart * sizeof(struct xpart));
  }

  hydro_space_init(&s->hs, s);

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

  if (verbose)
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

  if (posix_memalign((void **)&parts, part_align,
                     s->nr_parts * sizeof(struct part)) != 0)
    error("Failed to allocate new part array.");

  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate new gpart array.");

  if (posix_memalign((void **)&sparts, spart_align,
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

void space_generate_gas(struct space *s, const struct cosmology *cosmo,
                        int verbose) {

  if (verbose) message("Generating gas particles from gparts");

  /* Store the current values */
  const size_t nr_parts = s->nr_parts;
  const size_t nr_gparts = s->nr_gparts;

  if (nr_parts != 0)
    error("Generating gas particles from DM but gas already exists!");

  s->size_parts = s->nr_parts = nr_gparts;
  s->size_gparts = s->nr_gparts = 2 * nr_gparts;

  /* Allocate space for new particles */
  struct part *parts = NULL;
  struct gpart *gparts = NULL;

  if (posix_memalign((void **)&parts, part_align,
                     s->nr_parts * sizeof(struct part)) != 0)
    error("Failed to allocate new part array.");

  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate new gpart array.");

  /* Start by copying the gparts */
  memcpy(gparts, s->gparts, nr_gparts * sizeof(struct gpart));
  memcpy(gparts + nr_gparts, s->gparts, nr_gparts * sizeof(struct gpart));

  /* And zero the parts */
  bzero(parts, s->nr_parts * sizeof(struct part));

  /* Compute some constants */
  const double mass_ratio = cosmo->Omega_b / cosmo->Omega_m;
  const double bg_density = cosmo->Omega_m * cosmo->critical_density;
  const double bg_density_inv = 1. / bg_density;

  /* Update the particle properties */
  for (size_t i = 0; i < nr_gparts; ++i) {

    struct part *p = &parts[i];
    struct gpart *gp_gas = &gparts[i];
    struct gpart *gp_dm = &gparts[nr_gparts + i];

    /* Set the IDs */
    p->id = gp_gas->id_or_neg_offset * 2 + 1;
    gp_dm->id_or_neg_offset *= 2;

    if (gp_dm->id_or_neg_offset <= 0) error("DM particle ID overflowd");

    if (p->id <= 0) error("gas particle ID overflowd");

    /* Set the links correctly */
    p->gpart = gp_gas;
    gp_gas->id_or_neg_offset = -i;
    gp_gas->type = swift_type_gas;

    /* Compute positions shift */
    const double d = cbrt(gp_dm->mass * bg_density_inv);
    const double shift_dm = d * mass_ratio;
    const double shift_gas = d * (1. - mass_ratio);

    /* Set the masses */
    gp_dm->mass *= (1. - mass_ratio);
    gp_gas->mass *= mass_ratio;
    hydro_set_mass(p, gp_gas->mass);

    /* Set the new positions */
    gp_dm->x[0] += shift_dm;
    gp_dm->x[1] += shift_dm;
    gp_dm->x[2] += shift_dm;
    gp_gas->x[0] += shift_gas;
    gp_gas->x[1] += shift_gas;
    gp_gas->x[2] += shift_gas;
    p->x[0] = gp_gas->x[0];
    p->x[1] = gp_gas->x[1];
    p->x[2] = gp_gas->x[2];

    /* Also copy the velocities */
    p->v[0] = gp_gas->v_full[0];
    p->v[1] = gp_gas->v_full[1];
    p->v[2] = gp_gas->v_full[2];

    /* Set the smoothing length to the mean inter-particle separation */
    p->h = 30. * d;
  }

  /* Replace the content of the space */
  free(s->gparts);
  s->parts = parts;
  s->gparts = gparts;
}

/**
 * @brief Verify that the matter content matches the cosmology model.
 *
 * @param s The #space.
 * @param cosmo The current cosmology model.
 * @param rank The MPI rank of this #space.
 */
void space_check_cosmology(struct space *s, const struct cosmology *cosmo,
                           int rank) {

  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;

  /* Sum up the mass in this space */
  double mass = 0.;
  for (size_t i = 0; i < nr_gparts; ++i) {
    mass += gparts[i].mass;
  }

/* Reduce the total mass */
#ifdef WITH_MPI
  double total_mass;
  MPI_Reduce(&mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  double total_mass = mass;
#endif

  if (rank == 0) {

    const double volume = s->dim[0] * s->dim[1] * s->dim[2];

    /* Current Hubble constant */
    const double H = cosmo->H;

    /* z=0 Hubble parameter */
    const double H0 = cosmo->H0;

    /* Critical density at z=0 */
    const double rho_crit0 = cosmo->critical_density * H0 * H0 / (H * H);

    /* Compute the mass density */
    const double Omega_m = (total_mass / volume) / rho_crit0;

    if (fabs(Omega_m - cosmo->Omega_m) > 1e-3)
      error(
          "The matter content of the simulation does not match the cosmology "
          "in the parameter file cosmo.Omega_m=%e Omega_m=%e",
          cosmo->Omega_m, Omega_m);
  }
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
  free(s->local_cells_top);
  free(s->parts);
  free(s->xparts);
  free(s->gparts);
  free(s->sparts);
}

/**
 * @brief Write the space struct and its contents to the given FILE as a
 * stream of bytes.
 *
 * @param s the space
 * @param stream the file stream
 */
void space_struct_dump(struct space *s, FILE *stream) {

  restart_write_blocks(s, sizeof(struct space), 1, stream, "space",
                       "space struct");

  /* More things to write. */
  if (s->nr_parts > 0) {
    restart_write_blocks(s->parts, s->nr_parts, sizeof(struct part), stream,
                         "parts", "parts");
    restart_write_blocks(s->xparts, s->nr_parts, sizeof(struct xpart), stream,
                         "xparts", "xparts");
  }
  if (s->nr_gparts > 0)
    restart_write_blocks(s->gparts, s->nr_gparts, sizeof(struct gpart), stream,
                         "gparts", "gparts");

  if (s->nr_sparts > 0)
    restart_write_blocks(s->sparts, s->nr_sparts, sizeof(struct spart), stream,
                         "sparts", "sparts");
}

/**
 * @brief Re-create a space struct and its contents from the given FILE
 *        stream.
 *
 * @param s the space
 * @param stream the file stream
 */
void space_struct_restore(struct space *s, FILE *stream) {

  restart_read_blocks(s, sizeof(struct space), 1, stream, NULL, "space struct");

  /* Things that should be reconstructed in a rebuild. */
  s->cells_top = NULL;
  s->cells_sub = NULL;
  s->multipoles_top = NULL;
  s->multipoles_sub = NULL;
  s->local_cells_top = NULL;
  s->grav_top_level = NULL;
#ifdef WITH_MPI
  s->parts_foreign = NULL;
  s->size_parts_foreign = 0;
  s->gparts_foreign = NULL;
  s->size_gparts_foreign = 0;
  s->sparts_foreign = NULL;
  s->size_sparts_foreign = 0;
#endif

  /* More things to read. */
  s->parts = NULL;
  s->xparts = NULL;
  if (s->nr_parts > 0) {

    /* Need the memory for these. */
    if (posix_memalign((void **)&s->parts, part_align,
                       s->size_parts * sizeof(struct part)) != 0)
      error("Failed to allocate restore part array.");
    if (posix_memalign((void **)&s->xparts, xpart_align,
                       s->size_parts * sizeof(struct xpart)) != 0)
      error("Failed to allocate restore xpart array.");

    restart_read_blocks(s->parts, s->nr_parts, sizeof(struct part), stream,
                        NULL, "parts");
    restart_read_blocks(s->xparts, s->nr_parts, sizeof(struct xpart), stream,
                        NULL, "xparts");
  }
  s->gparts = NULL;
  if (s->nr_gparts > 0) {
    if (posix_memalign((void **)&s->gparts, gpart_align,
                       s->size_gparts * sizeof(struct gpart)) != 0)
      error("Failed to allocate restore gpart array.");

    restart_read_blocks(s->gparts, s->nr_gparts, sizeof(struct gpart), stream,
                        NULL, "gparts");
  }

  s->sparts = NULL;
  if (s->nr_sparts > 0) {
    if (posix_memalign((void **)&s->sparts, spart_align,
                       s->size_sparts * sizeof(struct spart)) != 0)
      error("Failed to allocate restore spart array.");

    restart_read_blocks(s->sparts, s->nr_sparts, sizeof(struct spart), stream,
                        NULL, "sparts");
  }

  /* Need to reconnect the gravity parts to their hydro and star particles. */
  /* Re-link the parts. */
  if (s->nr_parts > 0 && s->nr_gparts > 0)
    part_relink_parts_to_gparts(s->gparts, s->nr_gparts, s->parts);

  /* Re-link the sparts. */
  if (s->nr_sparts > 0 && s->nr_gparts > 0)
    part_relink_sparts_to_gparts(s->gparts, s->nr_gparts, s->sparts);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that everything is correct */
  part_verify_links(s->parts, s->gparts, s->sparts, s->nr_parts, s->nr_gparts,
                    s->nr_sparts, 1);
#endif
}
