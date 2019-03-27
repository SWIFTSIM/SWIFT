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
#include "tracers.h"

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize_pair_hydro = space_subsize_pair_hydro_default;
int space_subsize_self_hydro = space_subsize_self_hydro_default;
int space_subsize_pair_grav = space_subsize_pair_grav_default;
int space_subsize_self_grav = space_subsize_self_grav_default;
int space_subsize_pair_stars = space_subsize_pair_stars_default;
int space_subsize_self_stars = space_subsize_self_stars_default;
int space_subdepth_diff_grav = space_subdepth_diff_grav_default;
int space_maxsize = space_maxsize_default;

/*! Number of extra #part we allocate memory for per top-level cell */
int space_extra_parts = space_extra_parts_default;

/*! Number of extra #spart we allocate memory for per top-level cell */
int space_extra_sparts = space_extra_sparts_default;

/*! Number of extra #gpart we allocate memory for per top-level cell */
int space_extra_gparts = space_extra_gparts_default;

/*! Expected maximal number of strays received at a rebuild */
int space_expected_max_nr_strays = space_expected_max_nr_strays_default;
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
  int *ind;
  int *cell_counts;
  size_t count_inhibited_part;
  size_t count_inhibited_gpart;
  size_t count_inhibited_spart;
  size_t count_extra_part;
  size_t count_extra_gpart;
  size_t count_extra_spart;
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

        if (s->with_self_gravity) {
          c->progeny[k]->grav.multipole->next = *multipole_rec_begin;
          *multipole_rec_begin = c->progeny[k]->grav.multipole;
        }

        if (*cell_rec_end == NULL) *cell_rec_end = *cell_rec_begin;
        if (s->with_self_gravity && *multipole_rec_end == NULL)
          *multipole_rec_end = *multipole_rec_begin;

        c->progeny[k]->grav.multipole = NULL;
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
    c->hydro.sorts = NULL;
    c->stars.sorts = NULL;
    c->nr_tasks = 0;
    c->grav.nr_mm_tasks = 0;
    c->hydro.density = NULL;
    c->hydro.gradient = NULL;
    c->hydro.force = NULL;
    c->hydro.limiter = NULL;
    c->grav.grav = NULL;
    c->grav.mm = NULL;
    c->hydro.dx_max_part = 0.0f;
    c->hydro.dx_max_sort = 0.0f;
    c->stars.dx_max_part = 0.f;
    c->stars.dx_max_sort = 0.f;
    c->hydro.sorted = 0;
    c->stars.sorted = 0;
    c->hydro.count = 0;
    c->hydro.count_total = 0;
    c->hydro.updated = 0;
    c->hydro.inhibited = 0;
    c->grav.count = 0;
    c->grav.count_total = 0;
    c->grav.updated = 0;
    c->grav.inhibited = 0;
    c->stars.count = 0;
    c->stars.count_total = 0;
    c->stars.updated = 0;
    c->stars.inhibited = 0;
    c->grav.init = NULL;
    c->grav.init_out = NULL;
    c->hydro.extra_ghost = NULL;
    c->hydro.ghost_in = NULL;
    c->hydro.ghost_out = NULL;
    c->hydro.ghost = NULL;
    c->stars.ghost_in = NULL;
    c->stars.ghost_out = NULL;
    c->stars.ghost = NULL;
    c->stars.density = NULL;
    c->stars.feedback = NULL;
    c->kick1 = NULL;
    c->kick2 = NULL;
    c->timestep = NULL;
    c->timestep_limiter = NULL;
    c->end_force = NULL;
    c->hydro.drift = NULL;
    c->grav.drift = NULL;
    c->grav.drift_out = NULL;
    c->hydro.cooling = NULL;
    c->grav.long_range = NULL;
    c->grav.down_in = NULL;
    c->grav.down = NULL;
    c->grav.mesh = NULL;
    c->super = c;
    c->hydro.super = c;
    c->grav.super = c;
    c->hydro.parts = NULL;
    c->hydro.xparts = NULL;
    c->grav.parts = NULL;
    c->stars.parts = NULL;
    c->hydro.do_sub_sort = 0;
    c->stars.do_sub_sort = 0;
    c->grav.do_sub_drift = 0;
    c->hydro.do_sub_drift = 0;
    c->hydro.do_sub_limiter = 0;
    c->hydro.do_limiter = 0;
    c->hydro.ti_end_min = -1;
    c->hydro.ti_end_max = -1;
    c->grav.ti_end_min = -1;
    c->grav.ti_end_max = -1;
    c->stars.ti_end_min = -1;
#ifdef SWIFT_DEBUG_CHECKS
    c->cellID = 0;
#endif
    if (s->with_self_gravity)
      bzero(c->grav.multipole, sizeof(struct gravity_tensors));
    for (int i = 0; i < 13; i++) {
      if (c->hydro.sort[i] != NULL) {
        free(c->hydro.sort[i]);
        c->hydro.sort[i] = NULL;
      }
      if (c->stars.sort[i] != NULL) {
        free(c->stars.sort[i]);
        c->stars.sort[i] = NULL;
      }
    }
#if WITH_MPI
    c->mpi.tag = -1;

    c->mpi.hydro.recv_xv = NULL;
    c->mpi.hydro.recv_rho = NULL;
    c->mpi.hydro.recv_gradient = NULL;
    c->mpi.grav.recv = NULL;
    c->mpi.recv_ti = NULL;
    c->mpi.limiter.recv = NULL;

    c->mpi.hydro.send_xv = NULL;
    c->mpi.hydro.send_rho = NULL;
    c->mpi.hydro.send_gradient = NULL;
    c->mpi.grav.send = NULL;
    c->mpi.send_ti = NULL;
    c->mpi.limiter.send = NULL;
#endif
  }
}

/**
 * @brief Free up any allocated cells.
 */
void space_free_cells(struct space *s) {

  ticks tic = getticks();

  threadpool_map(&s->e->threadpool, space_rebuild_recycle_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), 0, s);
  s->maxdepth = 0;

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Re-build the top-level cell grid.
 *
 * @param s The #space.
 * @param verbose Print messages to stdout or not.
 */
void space_regrid(struct space *s, int verbose) {

  const size_t nr_parts = s->nr_parts;
  const size_t nr_sparts = s->nr_sparts;
  const ticks tic = getticks();
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;

  /* Run through the cells and get the current h_max. */
  // tic = getticks();
  float h_max = s->cell_min / kernel_gamma / space_stretch;
  if (nr_parts > 0) {

    /* Can we use the list of local non-empty top-level cells? */
    if (s->local_cells_with_particles_top != NULL) {
      for (int k = 0; k < s->nr_local_cells_with_particles; ++k) {
        const struct cell *c =
            &s->cells_top[s->local_cells_with_particles_top[k]];
        if (c->hydro.h_max > h_max) {
          h_max = c->hydro.h_max;
        }
        if (c->stars.h_max > h_max) {
          h_max = c->stars.h_max;
        }
      }

      /* Can we instead use all the top-level cells? */
    } else if (s->cells_top != NULL) {
      for (int k = 0; k < s->nr_cells; k++) {
        const struct cell *c = &s->cells_top[k];
        if (c->nodeID == engine_rank && c->hydro.h_max > h_max) {
          h_max = c->hydro.h_max;
        }
        if (c->nodeID == engine_rank && c->stars.h_max > h_max) {
          h_max = c->stars.h_max;
        }
      }

      /* Last option: run through the particles */
    } else {
      for (size_t k = 0; k < nr_parts; k++) {
        if (s->parts[k].h > h_max) h_max = s->parts[k].h;
      }
      for (size_t k = 0; k < nr_sparts; k++) {
        if (s->sparts[k].h > h_max) h_max = s->sparts[k].h;
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
  const int no_regrid = (s->cells_top == NULL && oldnodeIDs == NULL);
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
      free(s->local_cells_with_tasks_top);
      free(s->local_cells_top);
      free(s->cells_with_particles_top);
      free(s->local_cells_with_particles_top);
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
    if (s->with_self_gravity) {
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

    /* Allocate the indices of local cells with tasks */
    if (posix_memalign((void **)&s->local_cells_with_tasks_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level cells with tasks.");
    bzero(s->local_cells_with_tasks_top, s->nr_cells * sizeof(int));

    /* Allocate and initialise array of cell indices. */
    if (posix_memalign((void **)&s->cell_index, 32,
                       s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate list of cells for FOF search.");

    /* Set cell index into list of top-level cells. */
    for (int i = 0; i < s->nr_cells; i++) s->cell_index[i] = i;

    /* Allocate the indices of cells with particles */
    if (posix_memalign((void **)&s->cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of top-level cells with particles.");
    bzero(s->cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Allocate the indices of local cells with particles */
    if (posix_memalign((void **)&s->local_cells_with_particles_top,
                       SWIFT_STRUCT_ALIGNMENT, s->nr_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level cells with "
          "particles.");
    bzero(s->local_cells_with_particles_top, s->nr_cells * sizeof(int));

    /* Set the cells' locks */
    for (int k = 0; k < s->nr_cells; k++) {
      if (lock_init(&s->cells_top[k].hydro.lock) != 0)
        error("Failed to init spinlock for hydro.");
      if (lock_init(&s->cells_top[k].grav.plock) != 0)
        error("Failed to init spinlock for gravity.");
      if (lock_init(&s->cells_top[k].grav.mlock) != 0)
        error("Failed to init spinlock for multipoles.");
      if (lock_init(&s->cells_top[k].stars.lock) != 0)
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
          c->split = 0;
          c->hydro.count = 0;
          c->grav.count = 0;
          c->stars.count = 0;
          c->super = c;
          c->hydro.super = c;
          c->grav.super = c;
          c->hydro.ti_old_part = ti_current;
          c->grav.ti_old_part = ti_current;
          c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
          c->mpi.tag = -1;
          c->mpi.hydro.recv_xv = NULL;
          c->mpi.hydro.recv_rho = NULL;
          c->mpi.hydro.recv_gradient = NULL;
          c->mpi.hydro.send_xv = NULL;
          c->mpi.hydro.send_rho = NULL;
          c->mpi.hydro.send_gradient = NULL;
          c->mpi.grav.recv = NULL;
          c->mpi.grav.send = NULL;
#endif  // WITH_MPI
          if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
#ifdef SWIFT_DEBUG_CHECKS
          c->cellID = -last_cell_id;
          last_cell_id++;
#endif
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
#if defined(HAVE_PARMETIS) || defined(HAVE_METIS)
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
 * @brief Allocate memory for the extra particles used for on-the-fly creation.
 *
 * This rarely actually allocates memory. Most of the time, we convert
 * pre-allocated memory inot extra particles.
 *
 * This function also sets the extra particles' location to their top-level
 * cells. They can then be sorted into their correct memory position later on.
 *
 * @param s The current #space.
 * @param verbose Are we talkative?
 */
void space_allocate_extras(struct space *s, int verbose) {

  const int local_nodeID = s->e->nodeID;

  /* Anything to do here? (Abort if we don't want extras)*/
  if (space_extra_parts == 0 && space_extra_gparts == 0 &&
      space_extra_sparts == 0)
    return;

  /* The top-level cells */
  const struct cell *cells = s->cells_top;
  const double half_cell_width[3] = {0.5 * cells[0].width[0],
                                     0.5 * cells[0].width[1],
                                     0.5 * cells[0].width[2]};

  /* The current number of particles (including spare ones) */
  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;

  /* The current number of actual particles */
  size_t nr_actual_parts = nr_parts - s->nr_extra_parts;
  size_t nr_actual_gparts = nr_gparts - s->nr_extra_gparts;
  size_t nr_actual_sparts = nr_sparts - s->nr_extra_sparts;

  /* The number of particles we allocated memory for (MPI overhead) */
  size_t size_parts = s->size_parts;
  size_t size_gparts = s->size_gparts;
  size_t size_sparts = s->size_sparts;

  int local_cells = 0;
  for (int i = 0; i < s->nr_cells; ++i)
    if (s->cells_top[i].nodeID == local_nodeID) local_cells++;

  /* Number of extra particles we want for each type */
  const size_t expected_num_extra_parts = local_cells * space_extra_parts;
  const size_t expected_num_extra_gparts = local_cells * space_extra_gparts;
  const size_t expected_num_extra_sparts = local_cells * space_extra_sparts;

  if (verbose) {
    message("Currently have %zd/%zd/%zd real particles.", nr_actual_parts,
            nr_actual_gparts, nr_actual_sparts);
    message("Currently have %zd/%zd/%zd spaces for extra particles.",
            s->nr_extra_parts, s->nr_extra_gparts, s->nr_extra_sparts);
    message("Requesting space for future %zd/%zd/%zd part/gpart/sparts.",
            expected_num_extra_parts, expected_num_extra_gparts,
            expected_num_extra_sparts);
  }

  if (expected_num_extra_parts < s->nr_extra_parts)
    error("Reduction in top-level cells number not handled.");
  if (expected_num_extra_gparts < s->nr_extra_gparts)
    error("Reduction in top-level cells number not handled.");
  if (expected_num_extra_sparts < s->nr_extra_sparts)
    error("Reduction in top-level cells number not handled.");

  /* Do we have enough space for the extra gparts (i.e. we haven't used up any)
   * ? */
  if (nr_gparts + expected_num_extra_gparts > size_gparts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_gparts + expected_num_extra_gparts > size_gparts) {

      size_gparts = (nr_actual_gparts + expected_num_extra_gparts) *
                    engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating gparts array from %zd to %zd", s->size_gparts,
                size_gparts);

      /* Create more space for parts */
      struct gpart *gparts_new = NULL;
      if (posix_memalign((void **)&gparts_new, gpart_align,
                         sizeof(struct gpart) * size_gparts) != 0)
        error("Failed to allocate new gpart data");
      const ptrdiff_t delta = gparts_new - s->gparts;
      memcpy(gparts_new, s->gparts, sizeof(struct gpart) * s->size_gparts);
      free(s->gparts);
      s->gparts = gparts_new;

      /* Update the counter */
      s->size_gparts = size_gparts;

      /* We now need to reset all the part and spart pointers */
      for (size_t i = 0; i < nr_parts; ++i) {
        if (s->parts[i].time_bin != time_bin_not_created)
          s->parts[i].gpart += delta;
      }
      for (size_t i = 0; i < nr_sparts; ++i) {
        if (s->sparts[i].time_bin != time_bin_not_created)
          s->sparts[i].gpart += delta;
      }
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_gparts; i < nr_actual_gparts + expected_num_extra_gparts;
         ++i) {
      bzero(&s->gparts[i], sizeof(struct gpart));
      s->gparts[i].time_bin = time_bin_not_created;
      s->gparts[i].type = swift_type_dark_matter;
      s->gparts[i].id_or_neg_offset = -1;
    }

      /* Put the spare particles in their correct cell */
#ifdef WITH_MPI
    error("Need to do this correctly over MPI for only the local cells.");
#endif
    int count_in_cell = 0, current_cell = 0;
    size_t count_extra_gparts = 0;
    for (size_t i = 0; i < nr_actual_gparts + expected_num_extra_gparts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->gparts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->gparts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->gparts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->gparts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
        count_extra_gparts++;
      }

      /* Once we have reached the number of extra gpart per cell, we move to the
       * next */
      if (count_in_cell == space_extra_gparts) {
        ++current_cell;
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_gparts != expected_num_extra_gparts)
      error("Constructed the wrong number of extra gparts (%zd vs. %zd)",
            count_extra_gparts, expected_num_extra_gparts);
#endif

    /* Update the counters */
    s->nr_gparts = nr_actual_gparts + expected_num_extra_gparts;
    s->nr_extra_gparts = expected_num_extra_gparts;
  }

  /* Do we have enough space for the extra parts (i.e. we haven't used up any) ?
   */
  if (expected_num_extra_parts > s->nr_extra_parts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_parts + expected_num_extra_parts > size_parts) {

      size_parts = (nr_actual_parts + expected_num_extra_parts) *
                   engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating parts array from %zd to %zd", s->size_parts,
                size_parts);

      /* Create more space for parts */
      struct part *parts_new = NULL;
      if (posix_memalign((void **)&parts_new, part_align,
                         sizeof(struct part) * size_parts) != 0)
        error("Failed to allocate new part data");
      memcpy(parts_new, s->parts, sizeof(struct part) * s->size_parts);
      free(s->parts);
      s->parts = parts_new;

      /* Same for xparts */
      struct xpart *xparts_new = NULL;
      if (posix_memalign((void **)&xparts_new, xpart_align,
                         sizeof(struct xpart) * size_parts) != 0)
        error("Failed to allocate new xpart data");
      memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->size_parts);
      free(s->xparts);
      s->xparts = xparts_new;

      /* Update the counter */
      s->size_parts = size_parts;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_parts; i < nr_actual_parts + expected_num_extra_parts;
         ++i) {
      bzero(&s->parts[i], sizeof(struct part));
      bzero(&s->xparts[i], sizeof(struct xpart));
      s->parts[i].time_bin = time_bin_not_created;
      s->parts[i].id = -1;
    }

      /* Put the spare particles in their correct cell */
#ifdef WITH_MPI
    error("Need to do this correctly over MPI for only the local cells.");
#endif
    int count_in_cell = 0, current_cell = 0;
    size_t count_extra_parts = 0;
    for (size_t i = 0; i < nr_actual_parts + expected_num_extra_parts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->parts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->parts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->parts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->parts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
        count_extra_parts++;
      }

      /* Once we have reached the number of extra part per cell, we move to the
       * next */
      if (count_in_cell == space_extra_parts) {
        ++current_cell;
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_parts != expected_num_extra_parts)
      error("Constructed the wrong number of extra parts (%zd vs. %zd)",
            count_extra_parts, expected_num_extra_parts);
#endif

    /* Update the counters */
    s->nr_parts = nr_actual_parts + expected_num_extra_parts;
    s->nr_extra_parts = expected_num_extra_parts;
  }

  /* Do we have enough space for the extra sparts (i.e. we haven't used up any)
   * ? */
  if (nr_actual_sparts + expected_num_extra_sparts > nr_sparts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_sparts + expected_num_extra_sparts > size_sparts) {

      size_sparts = (nr_actual_sparts + expected_num_extra_sparts) *
                    engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating sparts array from %zd to %zd", s->size_sparts,
                size_sparts);

      /* Create more space for parts */
      struct spart *sparts_new = NULL;
      if (posix_memalign((void **)&sparts_new, spart_align,
                         sizeof(struct spart) * size_sparts) != 0)
        error("Failed to allocate new spart data");
      memcpy(sparts_new, s->sparts, sizeof(struct spart) * s->size_sparts);
      free(s->sparts);
      s->sparts = sparts_new;

      /* Update the counter */
      s->size_sparts = size_sparts;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_sparts; i < nr_actual_sparts + expected_num_extra_sparts;
         ++i) {
      bzero(&s->sparts[i], sizeof(struct spart));
      s->sparts[i].time_bin = time_bin_not_created;
      s->sparts[i].id = -42;
    }

      /* Put the spare particles in their correct cell */
#ifdef WITH_MPI
    error("Need to do this correctly over MPI for only the local cells.");
#endif
    int count_in_cell = 0, current_cell = 0;
    size_t count_extra_sparts = 0;
    for (size_t i = 0; i < nr_actual_sparts + expected_num_extra_sparts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->sparts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->sparts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->sparts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->sparts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
        count_extra_sparts++;
      }

      /* Once we have reached the number of extra spart per cell, we move to the
       * next */
      if (count_in_cell == space_extra_sparts) {
        ++current_cell;
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_sparts != expected_num_extra_sparts)
      error("Constructed the wrong number of extra sparts (%zd vs. %zd)",
            count_extra_sparts, expected_num_extra_sparts);
#endif

    /* Update the counters */
    s->nr_sparts = nr_actual_sparts + expected_num_extra_sparts;
    s->nr_extra_sparts = expected_num_extra_sparts;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  if ((nr_gparts > 0 && nr_parts > 0) || (nr_gparts > 0 && nr_sparts > 0))
    part_verify_links(s->parts, s->gparts, s->sparts, nr_parts, nr_gparts,
                      nr_sparts, verbose);
#endif
}

/**
 * @brief Re-build the cells as well as the tasks.
 *
 * @param s The #space in which to update the cells.
 * @param repartitioned Did we just repartition?
 * @param verbose Print messages to stdout or not
 */
void space_rebuild(struct space *s, int repartitioned, int verbose) {

  const ticks tic = getticks();

/* Be verbose about this. */
#ifdef SWIFT_DEBUG_CHECKS
  if (s->e->nodeID == 0 || verbose) message("(re)building space");
  fflush(stdout);
#endif

  /* Re-grid if necessary, or just re-set the cell data. */
  space_regrid(s, verbose);

  /* Allocate extra space for particles that will be created */
  if (s->with_star_formation) space_allocate_extras(s, verbose);

  struct cell *cells_top = s->cells_top;
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;
  const int local_nodeID = s->e->nodeID;

  /* The current number of particles */
  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;

  /* The number of particles we allocated memory for */
  size_t size_parts = s->size_parts;
  size_t size_gparts = s->size_gparts;
  size_t size_sparts = s->size_sparts;

  /* Counter for the number of inhibited particles found on the node */
  size_t count_inhibited_parts = 0;
  size_t count_inhibited_gparts = 0;
  size_t count_inhibited_sparts = 0;

  /* Counter for the number of extra particles found on the node */
  size_t count_extra_parts = 0;
  size_t count_extra_gparts = 0;
  size_t count_extra_sparts = 0;

  /* Number of particles we expect to have after strays exchange */
  const size_t h_index_size = size_parts + space_expected_max_nr_strays;
  const size_t g_index_size = size_gparts + space_expected_max_nr_strays;
  const size_t s_index_size = size_sparts + space_expected_max_nr_strays;

  /* Allocate arrays to store the indices of the cells where particles
     belong. We allocate extra space to allow for particles we may
     receive from other nodes */
  int *h_index = (int *)malloc(sizeof(int) * h_index_size);
  int *g_index = (int *)malloc(sizeof(int) * g_index_size);
  int *s_index = (int *)malloc(sizeof(int) * s_index_size);
  if (h_index == NULL || g_index == NULL || s_index == NULL)
    error("Failed to allocate temporary particle indices.");

  /* Allocate counters of particles that will land in each cell */
  int *cell_part_counts = (int *)malloc(sizeof(int) * s->nr_cells);
  int *cell_gpart_counts = (int *)malloc(sizeof(int) * s->nr_cells);
  int *cell_spart_counts = (int *)malloc(sizeof(int) * s->nr_cells);
  if (cell_part_counts == NULL || cell_gpart_counts == NULL ||
      cell_spart_counts == NULL)
    error("Failed to allocate cell particle count buffer.");

  /* Initialise the counters, including buffer space for future particles */
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_part_counts[i] = 0;
    cell_gpart_counts[i] = 0;
    cell_spart_counts[i] = 0;
  }

  /* Run through the particles and get their cell index. */
  if (nr_parts > 0)
    space_parts_get_cell_index(s, h_index, cell_part_counts,
                               &count_inhibited_parts, &count_extra_parts,
                               verbose);
  if (nr_gparts > 0)
    space_gparts_get_cell_index(s, g_index, cell_gpart_counts,
                                &count_inhibited_gparts, &count_extra_gparts,
                                verbose);
  if (nr_sparts > 0)
    space_sparts_get_cell_index(s, s_index, cell_spart_counts,
                                &count_inhibited_sparts, &count_extra_sparts,
                                verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Some safety checks */
  if (repartitioned && count_inhibited_parts)
    error("We just repartitioned but still found inhibited parts.");
  if (repartitioned && count_inhibited_sparts)
    error("We just repartitioned but still found inhibited sparts.");
  if (repartitioned && count_inhibited_gparts)
    error("We just repartitioned but still found inhibited gparts.");

  if (count_extra_parts != s->nr_extra_parts)
    error(
        "Number of extra parts in the part array not matching the space "
        "counter.");
  if (count_extra_gparts != s->nr_extra_gparts)
    error(
        "Number of extra gparts in the gpart array not matching the space "
        "counter.");
  if (count_extra_sparts != s->nr_extra_sparts)
    error(
        "Number of extra sparts in the spart array not matching the space "
        "counter.");
#endif

  /* Move non-local parts and inhibited parts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_parts > 0)) {
    for (size_t k = 0; k < nr_parts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (h_index[k] == -1 || cells_top[h_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
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
        memswap(&h_index[k], &h_index[nr_parts], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all parts are in the correct places. */
  size_t check_count_inhibited_part = 0;
  for (size_t k = 0; k < nr_parts; k++) {
    if (h_index[k] == -1 || cells_top[h_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local parts to send list");
    }
  }
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    if (h_index[k] != -1 && cells_top[h_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local parts from send list");
    }
    if (h_index[k] == -1) ++check_count_inhibited_part;
  }
  if (check_count_inhibited_part != count_inhibited_parts)
    error("Counts of inhibited particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local sparts and inhibited sparts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_sparts > 0)) {
    for (size_t k = 0; k < nr_sparts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (s_index[k] == -1 || cells_top[s_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
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
        memswap(&s_index[k], &s_index[nr_sparts], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all sparts are in the correct place. */
  size_t check_count_inhibited_spart = 0;
  for (size_t k = 0; k < nr_sparts; k++) {
    if (s_index[k] == -1 || cells_top[s_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local sparts to send list");
    }
  }
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    if (s_index[k] != -1 && cells_top[s_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local sparts from send list");
    }
    if (s_index[k] == -1) ++check_count_inhibited_spart;
  }
  if (check_count_inhibited_spart != count_inhibited_sparts)
    error("Counts of inhibited s-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local gparts and inhibited parts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_gparts > 0)) {
    for (size_t k = 0; k < nr_gparts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (g_index[k] == -1 || cells_top[g_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_gparts -= 1;

        /* Swap the particle */
        memswap(&s->gparts[k], &s->gparts[nr_gparts], sizeof(struct gpart));

        /* Swap the link with part/spart */
        if (s->gparts[k].type == swift_type_gas) {
          s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        } else if (s->gparts[k].type == swift_type_stars) {
          s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        }
        if (s->gparts[nr_gparts].type == swift_type_gas) {
          s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        } else if (s->gparts[nr_gparts].type == swift_type_stars) {
          s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        }

        /* Swap the index */
        memswap(&g_index[k], &g_index[nr_gparts], sizeof(int));
      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all gparts are in the correct place. */
  size_t check_count_inhibited_gpart = 0;
  for (size_t k = 0; k < nr_gparts; k++) {
    if (g_index[k] == -1 || cells_top[g_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local gparts to send list");
    }
  }
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    if (g_index[k] != -1 && cells_top[g_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local gparts from send list");
    }
    if (g_index[k] == -1) ++check_count_inhibited_gpart;
  }
  if (check_count_inhibited_gpart != count_inhibited_gparts)
    error("Counts of inhibited g-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

#ifdef WITH_MPI

  /* Exchange the strays, note that this potentially re-allocates
     the parts arrays. This can be skipped if we just repartitioned space
     as there should be no strays in that case */
  if (!repartitioned) {

    size_t nr_parts_exchanged = s->nr_parts - nr_parts;
    size_t nr_gparts_exchanged = s->nr_gparts - nr_gparts;
    size_t nr_sparts_exchanged = s->nr_sparts - nr_sparts;
    engine_exchange_strays(s->e, nr_parts, &h_index[nr_parts],
                           &nr_parts_exchanged, nr_gparts, &g_index[nr_gparts],
                           &nr_gparts_exchanged, nr_sparts, &s_index[nr_sparts],
                           &nr_sparts_exchanged);

    /* Set the new particle counts. */
    s->nr_parts = nr_parts + nr_parts_exchanged;
    s->nr_gparts = nr_gparts + nr_gparts_exchanged;
    s->nr_sparts = nr_sparts + nr_sparts_exchanged;

  } else {
#ifdef SWIFT_DEBUG_CHECKS
    if (s->nr_parts != nr_parts)
      error("Number of parts changing after repartition");
    if (s->nr_sparts != nr_sparts)
      error("Number of sparts changing after repartition");
    if (s->nr_gparts != nr_gparts)
      error("Number of gparts changing after repartition");
#endif
  }

  /* Clear non-local cell counts. */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID != local_nodeID) {
      cell_part_counts[k] = 0;
      cell_spart_counts[k] = 0;
      cell_gpart_counts[k] = 0;
    }
  }

  /* Re-allocate the index array for the parts if needed.. */
  if (s->nr_parts + 1 > h_index_size) {
    int *ind_new;
    if ((ind_new = (int *)malloc(sizeof(int) * (s->nr_parts + 1))) == NULL)
      error("Failed to allocate temporary particle indices.");
    memcpy(ind_new, h_index, sizeof(int) * nr_parts);
    free(h_index);
    h_index = ind_new;
  }

  /* Re-allocate the index array for the sparts if needed.. */
  if (s->nr_sparts + 1 > s_index_size) {
    int *sind_new;
    if ((sind_new = (int *)malloc(sizeof(int) * (s->nr_sparts + 1))) == NULL)
      error("Failed to allocate temporary s-particle indices.");
    memcpy(sind_new, s_index, sizeof(int) * nr_sparts);
    free(s_index);
    s_index = sind_new;
  }

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih[3] = {s->iwidth[0], s->iwidth[1], s->iwidth[2]};

  /* Assign each received part to its cell. */
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    const struct part *const p = &s->parts[k];
    h_index[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cell_part_counts[h_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[h_index[k]].nodeID != local_nodeID)
      error("Received part that does not belong to me (nodeID=%i).",
            cells_top[h_index[k]].nodeID);
#endif
  }
  nr_parts = s->nr_parts;

  /* Assign each received spart to its cell. */
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    const struct spart *const sp = &s->sparts[k];
    s_index[k] =
        cell_getid(cdim, sp->x[0] * ih[0], sp->x[1] * ih[1], sp->x[2] * ih[2]);
    cell_spart_counts[s_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[s_index[k]].nodeID != local_nodeID)
      error("Received s-part that does not belong to me (nodeID=%i).",
            cells_top[s_index[k]].nodeID);
#endif
  }
  nr_sparts = s->nr_sparts;

#else /* WITH_MPI */

  /* Update the part and spart counters */
  s->nr_parts = nr_parts;
  s->nr_sparts = nr_sparts;

#endif /* WITH_MPI */

  /* Sort the parts according to their cells. */
  if (nr_parts > 0)
    space_parts_sort(s->parts, s->xparts, h_index, cell_part_counts,
                     s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < nr_parts; k++) {
    const struct part *p = &s->parts[k];

    if (p->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_ind =
        cell_getid(s->cdim, p->x[0] * s->iwidth[0], p->x[1] * s->iwidth[1],
                   p->x[2] * s->iwidth[2]);

    /* New cell of this part */
    const struct cell *c = &s->cells_top[new_ind];

    if (h_index[k] != new_ind)
      error("part's new cell index not matching sorted index.");

    if (p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] < c->loc[1] || p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] < c->loc[2] || p->x[2] > c->loc[2] + c->width[2])
      error("part not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Sort the sparts according to their cells. */
  if (nr_sparts > 0)
    space_sparts_sort(s->sparts, s_index, cell_spart_counts, s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < nr_sparts; k++) {
    const struct spart *sp = &s->sparts[k];

    if (sp->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_sind =
        cell_getid(s->cdim, sp->x[0] * s->iwidth[0], sp->x[1] * s->iwidth[1],
                   sp->x[2] * s->iwidth[2]);

    /* New cell of this spart */
    const struct cell *c = &s->cells_top[new_sind];

    if (s_index[k] != new_sind)
      error("spart's new cell index not matching sorted index.");

    if (sp->x[0] < c->loc[0] || sp->x[0] > c->loc[0] + c->width[0] ||
        sp->x[1] < c->loc[1] || sp->x[1] > c->loc[1] + c->width[1] ||
        sp->x[2] < c->loc[2] || sp->x[2] > c->loc[2] + c->width[2])
      error("spart not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_index = 0;
  h_index[nr_parts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_parts; k++) {
    if (h_index[k] < h_index[k + 1]) {
      cells_top[h_index[k]].hydro.count =
          k - last_index + 1 - space_extra_parts;
      last_index = k + 1;
    }
  }

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_sindex = 0;
  s_index[nr_sparts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_sparts; k++) {
    if (s_index[k] < s_index[k + 1]) {
      cells_top[s_index[k]].stars.count =
          k - last_sindex + 1 - space_extra_sparts;
      last_sindex = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  free(h_index);
  free(cell_part_counts);
  free(s_index);
  free(cell_spart_counts);

#ifdef WITH_MPI

  /* Re-allocate the index array for the gparts if needed.. */
  if (s->nr_gparts + 1 > g_index_size) {
    int *gind_new;
    if ((gind_new = (int *)malloc(sizeof(int) * (s->nr_gparts + 1))) == NULL)
      error("Failed to allocate temporary g-particle indices.");
    memcpy(gind_new, g_index, sizeof(int) * nr_gparts);
    free(g_index);
    g_index = gind_new;
  }

  /* Assign each received gpart to its cell. */
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    const struct gpart *const p = &s->gparts[k];
    g_index[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cell_gpart_counts[g_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[g_index[k]].nodeID != s->e->nodeID)
      error("Received g-part that does not belong to me (nodeID=%i).",
            cells_top[g_index[k]].nodeID);
#endif
  }
  nr_gparts = s->nr_gparts;

#else /* WITH_MPI */

  /* Update the gpart counter */
  s->nr_gparts = nr_gparts;

#endif /* WITH_MPI */

  /* Mark that there are no inhibited particles left */
  s->nr_inhibited_parts = 0;
  s->nr_inhibited_gparts = 0;
  s->nr_inhibited_sparts = 0;

  /* Sort the gparts according to their cells. */
  if (nr_gparts > 0)
    space_gparts_sort(s->gparts, s->parts, s->sparts, g_index,
                      cell_gpart_counts, s->nr_cells);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < nr_gparts; k++) {
    const struct gpart *gp = &s->gparts[k];

    if (gp->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_gind =
        cell_getid(s->cdim, gp->x[0] * s->iwidth[0], gp->x[1] * s->iwidth[1],
                   gp->x[2] * s->iwidth[2]);

    /* New cell of this gpart */
    const struct cell *c = &s->cells_top[new_gind];

    if (g_index[k] != new_gind)
      error("gpart's new cell index not matching sorted index.");

    if (gp->x[0] < c->loc[0] || gp->x[0] > c->loc[0] + c->width[0] ||
        gp->x[1] < c->loc[1] || gp->x[1] > c->loc[1] + c->width[1] ||
        gp->x[2] < c->loc[2] || gp->x[2] > c->loc[2] + c->width[2])
      error("gpart not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_gindex = 0;
  g_index[nr_gparts] = s->nr_cells;
  for (size_t k = 0; k < nr_gparts; k++) {
    if (g_index[k] < g_index[k + 1]) {
      cells_top[g_index[k]].grav.count =
          k - last_gindex + 1 - space_extra_gparts;
      last_gindex = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  free(g_index);
  free(cell_gpart_counts);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  if ((nr_gparts > 0 && nr_parts > 0) || (nr_gparts > 0 && nr_sparts > 0))
    part_verify_links(s->parts, s->gparts, s->sparts, nr_parts, nr_gparts,
                      nr_sparts, verbose);
#endif

  /* Hook the cells up to the parts. Make list of local and non-empty cells */
  ticks tic2 = getticks();
  struct part *finger = s->parts;
  struct xpart *xfinger = s->xparts;
  struct gpart *gfinger = s->gparts;
  struct spart *sfinger = s->sparts;
  s->nr_cells_with_particles = 0;
  s->nr_local_cells_with_particles = 0;
  s->nr_local_cells = 0;
  for (int k = 0; k < s->nr_cells; k++) {
    struct cell *restrict c = &cells_top[k];
    c->hydro.ti_old_part = ti_current;
    c->grav.ti_old_part = ti_current;
    c->grav.ti_old_multipole = ti_current;

#ifdef SWIFT_DEBUG_CHECKS
    c->cellID = -last_cell_id;
    last_cell_id++;
#endif

    const int is_local = (c->nodeID == engine_rank);
    const int has_particles =
        (c->hydro.count > 0) || (c->grav.count > 0) || (c->stars.count > 0);

    if (is_local) {
      c->hydro.parts = finger;
      c->hydro.xparts = xfinger;
      c->grav.parts = gfinger;
      c->stars.parts = sfinger;

      c->hydro.count_total = c->hydro.count + space_extra_parts;
      c->grav.count_total = c->grav.count + space_extra_gparts;
      c->stars.count_total = c->stars.count + space_extra_sparts;

      finger = &finger[c->hydro.count_total];
      xfinger = &xfinger[c->hydro.count_total];
      gfinger = &gfinger[c->grav.count_total];
      sfinger = &sfinger[c->stars.count_total];

      /* Add this cell to the list of local cells */
      s->local_cells_top[s->nr_local_cells] = k;
      s->nr_local_cells++;
    }

    if (is_local && has_particles) {

      /* Add this cell to the list of non-empty cells */
      s->local_cells_with_particles_top[s->nr_local_cells_with_particles] = k;
      s->nr_local_cells_with_particles++;
    }
  }
  if (verbose) {
    message("Have %d local top-level cells with particles (total=%d)",
            s->nr_local_cells_with_particles, s->nr_cells);
    message("Have %d local top-level cells (total=%d)", s->nr_local_cells,
            s->nr_cells);
    message("hooking up cells took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());
  }

  /* Re-order the extra particles such that they are at the end of their cell's
     memory pool. */
  if (s->with_star_formation) space_reorder_extras(s, verbose);

  /* At this point, we have the upper-level cells. Now recursively split each
     cell to get the full AMR grid. */
  space_split(s, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the multipole construction went OK */
  if (s->with_self_gravity)
    for (int k = 0; k < s->nr_cells; k++)
      cell_check_multipole(&s->cells_top[k]);
#endif

  /* Clean up any stray sort indices in the cell buffer. */
  space_free_buff_sort_indices(s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * This is done in parallel using threads in the #threadpool.
 * Only do this for the local non-empty top-level cells.
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, int verbose) {

  const ticks tic = getticks();

  threadpool_map(&s->e->threadpool, space_split_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int), 0, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_reorder_extra_parts_mapper(void *map_data, int num_cells,
                                      void *extra_data) {

  struct cell *cells_top = (struct cell *)map_data;
  struct space *s = (struct space *)extra_data;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[ind];
    cell_reorder_extra_parts(c, c->hydro.parts - s->parts);
  }
}

void space_reorder_extra_gparts_mapper(void *map_data, int num_cells,
                                       void *extra_data) {

  struct cell *cells_top = (struct cell *)map_data;
  struct space *s = (struct space *)extra_data;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[ind];
    cell_reorder_extra_gparts(c, s->parts, s->sparts);
  }
}

void space_reorder_extra_sparts_mapper(void *map_data, int num_cells,
                                       void *extra_data) {

  struct cell *cells_top = (struct cell *)map_data;
  struct space *s = (struct space *)extra_data;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[ind];
    cell_reorder_extra_sparts(c, c->stars.parts - s->sparts);
  }
}

/**
 * @brief Re-orders the particles in each cell such that the extra particles
 * for on-the-fly creation are located at the end of their respective cells.
 *
 * This assumes that all the particles (real and extra) have already been sorted
 * in their correct top-level cell.
 *
 * @param s The #space to act upon.
 * @param verbose Are we talkative?
 */
void space_reorder_extras(struct space *s, int verbose) {

#ifdef WITH_MPI
  if (space_extra_parts || space_extra_gparts || space_extra_sparts)
    error("Need an MPI-proof version of this.");
#endif

  /* Re-order the gas particles */
  if (space_extra_parts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_parts_mapper,
                   s->cells_top, s->nr_cells, sizeof(struct cell), 0, s);

  /* Re-order the gravity particles */
  if (space_extra_gparts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_gparts_mapper,
                   s->cells_top, s->nr_cells, sizeof(struct cell), 0, s);

  /* Re-order the star particles */
  if (space_extra_sparts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_sparts_mapper,
                   s->cells_top, s->nr_cells, sizeof(struct cell), 0, s);
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
  size_t count_inhibited_part = 0;
  size_t count_extra_part = 0;

  /* Loop over the parts. */
  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    struct part *restrict p = &parts[k];

    const double old_pos_x = p->x[0];
    const double old_pos_y = p->x[1];
    const double old_pos_z = p->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

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

    if (p->time_bin == time_bin_inhibited) {
      /* Is this particle to be removed? */
      ind[k] = -1;
      ++count_inhibited_part;
    } else if (p->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_part;
    } else {
      /* Normal case: list its top-level cell index */
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
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra parts */
  if (count_inhibited_part)
    atomic_add(&data->count_inhibited_part, count_inhibited_part);
  if (count_extra_part) atomic_add(&data->count_extra_part, count_extra_part);

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
  size_t count_inhibited_gpart = 0;
  size_t count_extra_gpart = 0;

  for (int k = 0; k < nr_gparts; k++) {

    /* Get the particle */
    struct gpart *restrict gp = &gparts[k];

    const double old_pos_x = gp->x[0];
    const double old_pos_y = gp->x[1];
    const double old_pos_z = gp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

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

    if (gp->time_bin == time_bin_inhibited) {
      /* Is this particle to be removed? */
      ind[k] = -1;
      ++count_inhibited_gpart;
    } else if (gp->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_gpart;
    } else {
      /* List its top-level cell index */
      ind[k] = index;
      cell_counts[index]++;

      if (gp->type == swift_type_dark_matter) {

        /* Compute minimal mass */
        min_mass = min(min_mass, gp->mass);

        /* Compute sum of velocity norm */
        sum_vel_norm += gp->v_full[0] * gp->v_full[0] +
                        gp->v_full[1] * gp->v_full[1] +
                        gp->v_full[2] * gp->v_full[2];
      }

      /* Update the position */
      gp->x[0] = pos_x;
      gp->x[1] = pos_y;
      gp->x[2] = pos_z;
    }
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra gparts */
  if (count_inhibited_gpart)
    atomic_add(&data->count_inhibited_gpart, count_inhibited_gpart);
  if (count_extra_gpart)
    atomic_add(&data->count_extra_gpart, count_extra_gpart);

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
  size_t count_inhibited_spart = 0;
  size_t count_extra_spart = 0;

  for (int k = 0; k < nr_sparts; k++) {

    /* Get the particle */
    struct spart *restrict sp = &sparts[k];

    const double old_pos_x = sp->x[0];
    const double old_pos_y = sp->x[1];
    const double old_pos_z = sp->x[2];

#ifdef SWIFT_DEBUG_CHECKS
    if (!s->periodic) {
      if (old_pos_x < 0. || old_pos_x > dim_x)
        error("Particle outside of volume along X.");
      if (old_pos_y < 0. || old_pos_y > dim_y)
        error("Particle outside of volume along Y.");
      if (old_pos_z < 0. || old_pos_z > dim_z)
        error("Particle outside of volume along Z.");
    }
#endif

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

    /* Is this particle to be removed? */
    if (sp->time_bin == time_bin_inhibited) {
      ind[k] = -1;
      ++count_inhibited_spart;
    } else if (sp->time_bin == time_bin_not_created) {
      /* Is this a place-holder for on-the-fly creation? */
      ind[k] = index;
      cell_counts[index]++;
      ++count_extra_spart;
    } else {
      /* List its top-level cell index */
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
  }

  /* Write the counts back to the global array. */
  for (int k = 0; k < s->nr_cells; k++)
    if (cell_counts[k]) atomic_add(&data->cell_counts[k], cell_counts[k]);
  free(cell_counts);

  /* Write the count of inhibited and extra sparts */
  if (count_inhibited_spart)
    atomic_add(&data->count_inhibited_spart, count_inhibited_spart);
  if (count_extra_spart)
    atomic_add(&data->count_extra_spart, count_extra_spart);

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
 * @param count_inhibited_parts (return) The number of #part to remove.
 * @param count_extra_parts (return) The number of #part for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_parts_get_cell_index(struct space *s, int *ind, int *cell_counts,
                                size_t *count_inhibited_parts,
                                size_t *count_extra_parts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_part_mass = FLT_MAX;
  s->sum_part_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = ind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;

  threadpool_map(&s->e->threadpool, space_parts_get_cell_index_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), 0, &data);

  *count_inhibited_parts = data.count_inhibited_part;
  *count_extra_parts = data.count_extra_part;

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
 * @param count_inhibited_gparts (return) The number of #gpart to remove.
 * @param count_extra_gparts (return) The number of #gpart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_gparts_get_cell_index(struct space *s, int *gind, int *cell_counts,
                                 size_t *count_inhibited_gparts,
                                 size_t *count_extra_gparts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_gpart_mass = FLT_MAX;
  s->sum_gpart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = gind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;

  threadpool_map(&s->e->threadpool, space_gparts_get_cell_index_mapper,
                 s->gparts, s->nr_gparts, sizeof(struct gpart), 0, &data);

  *count_inhibited_gparts = data.count_inhibited_gpart;
  *count_extra_gparts = data.count_extra_gpart;

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
 * @param count_inhibited_sparts (return) The number of #spart to remove.
 * @param count_extra_sparts (return) The number of #spart for on-the-fly
 * creation.
 * @param verbose Are we talkative ?
 */
void space_sparts_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                 size_t *count_inhibited_sparts,
                                 size_t *count_extra_sparts, int verbose) {

  const ticks tic = getticks();

  /* Re-set the counters */
  s->min_spart_mass = FLT_MAX;
  s->sum_spart_vel_norm = 0.f;

  /* Pack the extra information */
  struct index_data data;
  data.s = s;
  data.ind = sind;
  data.cell_counts = cell_counts;
  data.count_inhibited_part = 0;
  data.count_inhibited_gpart = 0;
  data.count_inhibited_spart = 0;
  data.count_extra_part = 0;
  data.count_extra_gpart = 0;
  data.count_extra_spart = 0;

  threadpool_map(&s->e->threadpool, space_sparts_get_cell_index_mapper,
                 s->sparts, s->nr_sparts, sizeof(struct spart), 0, &data);

  *count_inhibited_sparts = data.count_inhibited_spart;
  *count_extra_sparts = data.count_extra_spart;

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
void space_parts_sort(struct part *parts, struct xpart *xparts,
                      int *restrict ind, int *restrict counts, int num_bins,
                      ptrdiff_t parts_offset) {
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
void space_sparts_sort(struct spart *sparts, int *restrict ind,
                       int *restrict counts, int num_bins,
                       ptrdiff_t sparts_offset) {
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
                       struct spart *sparts, int *restrict ind,
                       int *restrict counts, int num_bins) {
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
        } else if (gparts[j].type == swift_type_stars) {
          sparts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        }
      }
      gparts[k] = temp_gpart;
      ind[k] = target_cid;
      if (gparts[k].type == swift_type_gas) {
        parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_stars) {
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

  for (int i = 0; i < 13; i++) {
    if (c->hydro.sort[i] != NULL) {
      free(c->hydro.sort[i]);
      c->hydro.sort[i] = NULL;
    }
    if (c->stars.sort[i] != NULL) {
      free(c->stars.sort[i]);
      c->stars.sort[i] = NULL;
    }
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
    for (int k = 0; k < c->hydro.count; k++) fun(&c->hydro.parts[k], c, data);

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
    for (int k = 0; k < c->hydro.count; k++)
      fun(&c->hydro.parts[k], &c->hydro.xparts[k], c);

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
 *        c->hydro.count or @c NULL.
 * @param sbuff A buffer for particle sorting, should be of size at least
 *        c->stars.count or @c NULL.
 * @param gbuff A buffer for particle sorting, should be of size at least
 *        c->grav.count or @c NULL.
 */
void space_split_recursive(struct space *s, struct cell *c,
                           struct cell_buff *buff, struct cell_buff *sbuff,
                           struct cell_buff *gbuff) {

  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int with_self_gravity = s->with_self_gravity;
  const int depth = c->depth;
  int maxdepth = 0;
  float h_max = 0.0f;
  float stars_h_max = 0.f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct xpart *xparts = c->hydro.xparts;
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
#ifdef SWIFT_DEBUG_CHECKS
        if (parts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (parts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
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
#ifdef SWIFT_DEBUG_CHECKS
        if (gparts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (gparts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
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
#ifdef SWIFT_DEBUG_CHECKS
        if (sparts[k].time_bin == time_bin_inhibited)
          error("Inhibited particle present in space_split()");
        if (sparts[k].time_bin == time_bin_not_created)
          error("Extra particle present in space_split()");
#endif
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
  if ((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize))) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->hydro.count = 0;
      cp->grav.count = 0;
      cp->stars.count = 0;
      cp->hydro.count_total = 0;
      cp->grav.count_total = 0;
      cp->stars.count_total = 0;
      cp->hydro.ti_old_part = c->hydro.ti_old_part;
      cp->grav.ti_old_part = c->grav.ti_old_part;
      cp->grav.ti_old_multipole = c->grav.ti_old_multipole;
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
      cp->hydro.h_max = 0.f;
      cp->hydro.dx_max_part = 0.f;
      cp->hydro.dx_max_sort = 0.f;
      cp->stars.h_max = 0.f;
      cp->stars.dx_max_part = 0.f;
      cp->stars.dx_max_sort = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->super = NULL;
      cp->hydro.super = NULL;
      cp->grav.super = NULL;
      cp->hydro.do_sub_sort = 0;
      cp->stars.do_sub_sort = 0;
      cp->grav.do_sub_drift = 0;
      cp->hydro.do_sub_drift = 0;
      cp->hydro.do_sub_limiter = 0;
      cp->hydro.do_limiter = 0;
#ifdef WITH_MPI
      cp->mpi.tag = -1;
#endif  // WITH_MPI
#ifdef SWIFT_DEBUG_CHECKS
      cp->cellID = last_cell_id++;
#endif
    }

    /* Split the cell's partcle data. */
    cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts, buff,
               sbuff, gbuff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff;

    for (int k = 0; k < 8; k++) {

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->hydro.count == 0 && cp->grav.count == 0 && cp->stars.count == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {

        /* Recurse */
        space_split_recursive(s, cp, progeny_buff, progeny_sbuff,
                              progeny_gbuff);

        /* Update the pointers in the buffers */
        progeny_buff += cp->hydro.count;
        progeny_gbuff += cp->grav.count;
        progeny_sbuff += cp->stars.count;

        /* Update the cell-wide properties */
        h_max = max(h_max, cp->hydro.h_max);
        stars_h_max = max(stars_h_max, cp->stars.h_max);
        ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
        ti_hydro_end_max = max(ti_hydro_end_max, cp->hydro.ti_end_max);
        ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
        ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
        ti_gravity_end_max = max(ti_gravity_end_max, cp->grav.ti_end_max);
        ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);
        ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);

        /* Increase the depth */
        if (cp->maxdepth > maxdepth) maxdepth = cp->maxdepth;
      }
    }

    /* Deal with the multipole */
    if (s->with_self_gravity) {

      /* Reset everything */
      gravity_reset(c->grav.multipole);

      /* Compute CoM and bulk velocity from all progenies */
      double CoM[3] = {0., 0., 0.};
      double vel[3] = {0., 0., 0.};
      float max_delta_vel[3] = {0.f, 0.f, 0.f};
      float min_delta_vel[3] = {0.f, 0.f, 0.f};
      double mass = 0.;

      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct gravity_tensors *m = c->progeny[k]->grav.multipole;

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
      c->grav.multipole->CoM[0] = CoM[0] * inv_mass;
      c->grav.multipole->CoM[1] = CoM[1] * inv_mass;
      c->grav.multipole->CoM[2] = CoM[2] * inv_mass;
      c->grav.multipole->m_pole.vel[0] = vel[0] * inv_mass;
      c->grav.multipole->m_pole.vel[1] = vel[1] * inv_mass;
      c->grav.multipole->m_pole.vel[2] = vel[2] * inv_mass;

      /* Min max velocity along each axis */
      c->grav.multipole->m_pole.max_delta_vel[0] = max_delta_vel[0];
      c->grav.multipole->m_pole.max_delta_vel[1] = max_delta_vel[1];
      c->grav.multipole->m_pole.max_delta_vel[2] = max_delta_vel[2];
      c->grav.multipole->m_pole.min_delta_vel[0] = min_delta_vel[0];
      c->grav.multipole->m_pole.min_delta_vel[1] = min_delta_vel[1];
      c->grav.multipole->m_pole.min_delta_vel[2] = min_delta_vel[2];

      /* Now shift progeny multipoles and add them up */
      struct multipole temp;
      double r_max = 0.;
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct cell *cp = c->progeny[k];
          const struct multipole *m = &cp->grav.multipole->m_pole;

          /* Contribution to multipole */
          gravity_M2M(&temp, m, c->grav.multipole->CoM,
                      cp->grav.multipole->CoM);
          gravity_multipole_add(&c->grav.multipole->m_pole, &temp);

          /* Upper limit of max CoM<->gpart distance */
          const double dx =
              c->grav.multipole->CoM[0] - cp->grav.multipole->CoM[0];
          const double dy =
              c->grav.multipole->CoM[1] - cp->grav.multipole->CoM[1];
          const double dz =
              c->grav.multipole->CoM[2] - cp->grav.multipole->CoM[2];
          const double r2 = dx * dx + dy * dy + dz * dz;
          r_max = max(r_max, cp->grav.multipole->r_max + sqrt(r2));
        }
      }

      /* Alternative upper limit of max CoM<->gpart distance */
      const double dx =
          c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
              ? c->grav.multipole->CoM[0] - c->loc[0]
              : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
      const double dy =
          c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
              ? c->grav.multipole->CoM[1] - c->loc[1]
              : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
      const double dz =
          c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
              ? c->grav.multipole->CoM[2] - c->loc[2]
              : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];

      /* Take minimum of both limits */
      c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));

      /* Store the value at rebuild time */
      c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
      c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
      c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
      c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];

      /* We know the first-order multipole (dipole) is 0. */
      c->grav.multipole->m_pole.M_100 = 0.f;
      c->grav.multipole->m_pole.M_010 = 0.f;
      c->grav.multipole->m_pole.M_001 = 0.f;

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
    timebin_t stars_time_bin_min = num_time_bins;

    /* parts: Get dt_min/dt_max and h_max. */
    for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (parts[k].time_bin == time_bin_not_created)
        error("Extra particle present in space_split()");
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
      if (gparts[k].time_bin == time_bin_not_created)
        error("Extra g-particle present in space_split()");
      if (gparts[k].time_bin == time_bin_inhibited)
        error("Inhibited g-particle present in space_split()");
#endif
      gravity_time_bin_min = min(gravity_time_bin_min, gparts[k].time_bin);
      gravity_time_bin_max = max(gravity_time_bin_max, gparts[k].time_bin);
    }

    /* sparts: Get dt_min/dt_max */
    for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[k].time_bin == time_bin_not_created)
        error("Extra s-particle present in space_split()");
      if (sparts[k].time_bin == time_bin_inhibited)
        error("Inhibited s-particle present in space_split()");
#endif
      gravity_time_bin_min = min(gravity_time_bin_min, sparts[k].time_bin);
      gravity_time_bin_max = max(gravity_time_bin_max, sparts[k].time_bin);
      stars_time_bin_min = min(stars_time_bin_min, sparts[k].time_bin);

      stars_h_max = max(stars_h_max, sparts[k].h);

      /* Reset x_diff */
      sparts[k].x_diff[0] = 0.f;
      sparts[k].x_diff[1] = 0.f;
      sparts[k].x_diff[2] = 0.f;
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
    ti_stars_end_min = get_integer_time_end(ti_current, stars_time_bin_min);

    /* Construct the multipole and the centre of mass*/
    if (s->with_self_gravity) {
      if (gcount > 0) {

        gravity_P2M(c->grav.multipole, c->grav.parts, c->grav.count);

      } else {

        /* No gparts in that leaf cell */

        /* Set the values to something sensible */
        gravity_multipole_init(&c->grav.multipole->m_pole);
        if (c->nodeID == engine_rank) {
          c->grav.multipole->CoM[0] = c->loc[0] + c->width[0] / 2.;
          c->grav.multipole->CoM[1] = c->loc[1] + c->width[1] / 2.;
          c->grav.multipole->CoM[2] = c->loc[2] + c->width[2] / 2.;
          c->grav.multipole->r_max = 0.;
        }
      }

      /* Store the value at rebuild time */
      c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
      c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
      c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
      c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];
    }
  }

  /* Set the values for this cell. */
  c->hydro.h_max = h_max;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_end_max = ti_gravity_end_max;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.h_max = stars_h_max;
  c->maxdepth = maxdepth;

  /* Set ownership according to the start of the parts array. */
  if (s->nr_parts > 0)
    c->owner = ((c->hydro.parts - s->parts) % s->nr_parts) * s->nr_queues /
               s->nr_parts;
  else if (s->nr_sparts > 0)
    c->owner = ((c->stars.parts - s->sparts) % s->nr_sparts) * s->nr_queues /
               s->nr_sparts;
  else if (s->nr_gparts > 0)
    c->owner = ((c->grav.parts - s->gparts) % s->nr_gparts) * s->nr_queues /
               s->nr_gparts;
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
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    space_split_recursive(s, c, NULL, NULL, NULL);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
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
  if (lock_destroy(&c->lock) != 0 || lock_destroy(&c->grav.plock) != 0 ||
      lock_destroy(&c->mlock) != 0 || lock_destroy(&c->stars.lock) != 0)
    error("Failed to destroy spinlocks.");

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Hook the multipole back in the buffer */
  if (s->with_self_gravity) {
    c->grav.multipole->next = s->multipoles_sub;
    s->multipoles_sub = c->grav.multipole;
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
    if (lock_destroy(&c->lock) != 0 || lock_destroy(&c->grav.plock) != 0 ||
        lock_destroy(&c->mlock) != 0 || lock_destroy(&c->stars.lock) != 0)
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
  if (s->with_self_gravity) {
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
    if (s->with_self_gravity && s->multipoles_sub == NULL) {
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
    if (s->with_self_gravity) {
      cells[j]->grav.multipole = s->multipoles_sub;
      s->multipoles_sub = cells[j]->grav.multipole->next;
    }
  }

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);

  /* Init some things in the cell we just got. */
  for (int j = 0; j < nr_cells; j++) {
    for (int k = 0; k < 13; k++) {
      if (cells[j]->hydro.sort[k] != NULL) free(cells[j]->hydro.sort[k]);
      if (cells[j]->stars.sort[k] != NULL) free(cells[j]->stars.sort[k]);
    }
    struct gravity_tensors *temp = cells[j]->grav.multipole;
    bzero(cells[j], sizeof(struct cell));
    cells[j]->grav.multipole = temp;
    cells[j]->nodeID = -1;
    if (lock_init(&cells[j]->hydro.lock) != 0 ||
        lock_init(&cells[j]->grav.plock) != 0 ||
        lock_init(&cells[j]->grav.mlock) != 0 ||
        lock_init(&cells[j]->stars.lock) != 0)
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
    for (int k = 0; k < 13; k++) {
      if (finger->hydro.sort[k] != NULL) {
        free(finger->hydro.sort[k]);
        finger->hydro.sort[k] = NULL;
      }
      if (finger->stars.sort[k] != NULL) {
        free(finger->stars.sort[k]);
        finger->stars.sort[k] = NULL;
      }
    }
  }
}

/**
 * @brief Construct the list of top-level cells that have any tasks in
 * their hierarchy on this MPI rank. Also construct the list of top-level
 * cells on any rank that have > 0 particles (of any kind).
 *
 * This assumes the list has been pre-allocated at a regrid.
 *
 * @param s The #space.
 */
void space_list_useful_top_level_cells(struct space *s) {

  const ticks tic = getticks();

  s->nr_local_cells_with_tasks = 0;
  s->nr_cells_with_particles = 0;

  for (int i = 0; i < s->nr_cells; ++i) {
    struct cell *c = &s->cells_top[i];

    if (cell_has_tasks(c)) {
      s->local_cells_with_tasks_top[s->nr_local_cells_with_tasks] = i;
      s->nr_local_cells_with_tasks++;
    }

    const int has_particles =
        (c->hydro.count > 0) || (c->grav.count > 0) || (c->stars.count > 0) ||
        (c->grav.multipole != NULL && c->grav.multipole->m_pole.M_000 > 0.f);

    if (has_particles) {
      s->cells_with_particles_top[s->nr_cells_with_particles] = i;
      s->nr_cells_with_particles++;
    }
  }
  if (s->e->verbose) {
    message("Have %d local top-level cells with tasks (total=%d)",
            s->nr_local_cells_with_tasks, s->nr_cells);
    message("Have %d top-level cells with particles (total=%d)",
            s->nr_cells_with_particles, s->nr_cells);
  }

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
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

    else if (gp->type == swift_type_stars) {

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

  const ticks tic = getticks();

  if ((s->nr_gparts > 0 && s->nr_parts > 0) ||
      (s->nr_gparts > 0 && s->nr_sparts > 0))
    threadpool_map(&s->e->threadpool,
                   space_synchronize_particle_positions_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart), 0, (void *)s);

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
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
  const float hydro_h_min_ratio = e->hydro_properties->h_min_ratio;

  const struct gravity_props *grav_props = s->e->gravity_properties;
  const int with_gravity = e->policy & engine_policy_self_gravity;

  const struct chemistry_global_data *chemistry = e->chemistry;
  const struct cooling_function_data *cool_func = e->cooling_func;

  /* Check that the smoothing lengths are non-zero */
  for (int k = 0; k < count; k++) {
    if (p[k].h <= 0.)
      error("Invalid value of smoothing length for part %lld h=%e", p[k].id,
            p[k].h);

    if (with_gravity) {
      const struct gpart *gp = p[k].gpart;
      const float softening = gravity_get_softening(gp, grav_props);
      p->h = max(p->h, softening * hydro_h_min_ratio);
    }
  }

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
#ifdef WITH_LOGGER
    logger_part_data_init(&xp[k].logger_data);
#endif

    /* Overwrite the internal energy? */
    if (u_init > 0.f) hydro_set_init_internal_energy(&p[k], u_init);

    /* Also initialise the chemistry */
    chemistry_first_init_part(phys_const, us, cosmo, chemistry, &p[k], &xp[k]);

    /* And the cooling */
    cooling_first_init_part(phys_const, us, cosmo, cool_func, &p[k], &xp[k]);

    /* And the tracers */
    tracers_first_init_xpart(&p[k], &xp[k], us, phys_const, cosmo, hydro_props,
                             cool_func);

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
  const struct engine *e = s->e;

#ifdef SWIFT_DEBUG_CHECKS
  const ptrdiff_t delta = sp - s->sparts;
#endif

  const int with_feedback = (e->policy & engine_policy_feedback);

  const struct cosmology *cosmo = e->cosmology;
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

  /* Check that the smoothing lengths are non-zero */
  for (int k = 0; k < count; k++) {
    if (with_feedback && sp[k].h <= 0.)
      error("Invalid value of smoothing length for spart %lld h=%e", sp[k].id,
            sp[k].h);
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    stars_first_init_spart(&sp[k]);

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
 * Calls stars_first_init_spart() on all the particles
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

void space_init_sparts_mapper(void *restrict map_data, int scount,
                              void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  for (int k = 0; k < scount; k++) stars_init_spart(&sparts[k]);
}

/**
 * @brief Calls the #spart initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_sparts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, space_init_sparts_mapper, s->sparts,
                   s->nr_sparts, sizeof(struct spart), 0, NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_convert_quantities_mapper(void *restrict map_data, int count,
                                     void *restrict extra_data) {
  struct space *s = (struct space *)extra_data;
  const struct cosmology *cosmo = s->e->cosmology;
  const struct hydro_props *hydro_props = s->e->hydro_properties;
  struct part *restrict parts = (struct part *)map_data;
  const ptrdiff_t index = parts - s->parts;
  struct xpart *restrict xparts = s->xparts + index;

  /* Loop over all the particles ignoring the extra buffer ones for on-the-fly
   * creation */
  for (int k = 0; k < count; k++)
    if (parts[k].time_bin <= num_time_bins)
      hydro_convert_quantities(&parts[k], &xparts[k], cosmo, hydro_props);
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
 * @param sparts Array of stars particles.
 * @param Npart The number of Gas particles in the space.
 * @param Ngpart The number of Gravity particles in the space.
 * @param Nspart The number of stars particles in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param replicate How many replications along each direction do we want?
 * @param generate_gas_in_ics Are we generating gas particles from the gparts?
 * @param hydro flag whether we are doing hydro or not?
 * @param self_gravity flag whether we are doing gravity or not?
 * @param star_formation flag whether we are doing star formation or not?
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
                int replicate, int generate_gas_in_ics, int hydro,
                int self_gravity, int star_formation, int verbose,
                int dry_run) {

  /* Clean-up everything */
  bzero(s, sizeof(struct space));

  /* Store everything in the space. */
  s->dim[0] = dim[0];
  s->dim[1] = dim[1];
  s->dim[2] = dim[2];
  s->periodic = periodic;
  s->with_self_gravity = self_gravity;
  s->with_hydro = hydro;
  s->with_star_formation = star_formation;
  s->nr_parts = Npart;
  s->nr_gparts = Ngpart;
  s->nr_sparts = Nspart;
  s->size_parts = Npart;
  s->size_gparts = Ngpart;
  s->size_sparts = Nspart;
  s->nr_inhibited_parts = 0;
  s->nr_inhibited_gparts = 0;
  s->nr_inhibited_sparts = 0;
  s->nr_extra_parts = 0;
  s->nr_extra_gparts = 0;
  s->nr_extra_sparts = 0;
  s->parts = parts;
  s->gparts = gparts;
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
    space_generate_gas(s, cosmo, periodic, dim, verbose);
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
  space_subsize_pair_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_pair_stars",
                               space_subsize_pair_stars_default);
  space_subsize_self_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_stars",
                               space_subsize_self_stars_default);
  space_splitsize = parser_get_opt_param_int(
      params, "Scheduler:cell_split_size", space_splitsize_default);
  space_subdepth_diff_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_subdepth_diff_grav",
                               space_subdepth_diff_grav_default);
  space_extra_parts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_parts", space_extra_parts_default);
  space_extra_sparts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_sparts", space_extra_sparts_default);
  space_extra_gparts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_gparts", space_extra_gparts_default);

  if (verbose) {
    message("max_size set to %d split_size set to %d", space_maxsize,
            space_splitsize);
    message("subdepth_grav set to %d", space_subdepth_diff_grav);
    message("sub_size_pair_hydro set to %d, sub_size_self_hydro set to %d",
            space_subsize_pair_hydro, space_subsize_self_hydro);
    message("sub_size_pair_grav set to %d, sub_size_self_grav set to %d",
            space_subsize_pair_grav, space_subsize_self_grav);
    message("sub_size_pair_stars set to %d, sub_size_self_stars set to %d",
            space_subsize_pair_stars, space_subsize_self_stars);
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

#ifdef SWIFT_DEBUG_CHECKS
  last_cell_id = 1;
#endif

  /* Do we want any spare particles for on the fly cration? */
  if (!star_formation) space_extra_sparts = 0;

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

/**
 * @brief Duplicate all the dark matter particles to create the same number
 * of gas particles with mass ratios given by the cosmology.
 *
 * Note that this function alters the dark matter particle masses and positions.
 * Velocities are unchanged. We also leave the thermodynamic properties of the
 * gas un-initialised as they will be given a value from the parameter file at a
 * later stage.
 *
 * @param s The #space to create the particles in.
 * @param cosmo The current #cosmology model.
 * @param periodic Are we using periodic boundary conditions?
 * @param dim The size of the box (for periodic wrapping).
 * @param verbose Are we talkative?
 */
void space_generate_gas(struct space *s, const struct cosmology *cosmo,
                        int periodic, const double dim[3], int verbose) {

  /* Check that this is a sensible ting to do */
  if (!s->with_hydro)
    error(
        "Cannot generate gas from ICs if we are running without "
        "hydrodynamics. Need to run with -s and the corresponding "
        "hydrodynamics parameters in the YAML file.");

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
  const double bg_density = cosmo->Omega_m * cosmo->critical_density_0;
  const double bg_density_inv = 1. / bg_density;

  /* Update the particle properties */
  for (size_t i = 0; i < nr_gparts; ++i) {

    struct part *p = &parts[i];
    struct gpart *gp_gas = &gparts[i];
    struct gpart *gp_dm = &gparts[nr_gparts + i];

    /* Set the IDs */
    p->id = gp_gas->id_or_neg_offset * 2 + 1;
    gp_dm->id_or_neg_offset *= 2;

    if (gp_dm->id_or_neg_offset < 0)
      error("DM particle ID overflowd (DM id=%lld gas id=%lld)",
            gp_dm->id_or_neg_offset, p->id);

    if (p->id < 0) error("gas particle ID overflowd (id=%lld)", p->id);

    /* Set the links correctly */
    p->gpart = gp_gas;
    gp_gas->id_or_neg_offset = -i;
    gp_gas->type = swift_type_gas;

    /* Compute positions shift */
    const double d = cbrt(gp_dm->mass * bg_density_inv);
    const double shift_dm = 0.5 * d * mass_ratio;
    const double shift_gas = 0.5 * d * (1. - mass_ratio);

    /* Set the masses */
    gp_dm->mass *= (1. - mass_ratio);
    gp_gas->mass *= mass_ratio;
    hydro_set_mass(p, gp_gas->mass);

    /* Set the new positions */
    gp_dm->x[0] += shift_dm;
    gp_dm->x[1] += shift_dm;
    gp_dm->x[2] += shift_dm;
    gp_gas->x[0] -= shift_gas;
    gp_gas->x[1] -= shift_gas;
    gp_gas->x[2] -= shift_gas;

    /* Make sure the positions are identical between linked particles */
    p->x[0] = gp_gas->x[0];
    p->x[1] = gp_gas->x[1];
    p->x[2] = gp_gas->x[2];

    /* Box-wrap the whole thing to be safe */
    if (periodic) {
      gp_dm->x[0] = box_wrap(gp_dm->x[0], 0., dim[0]);
      gp_dm->x[1] = box_wrap(gp_dm->x[1], 0., dim[1]);
      gp_dm->x[2] = box_wrap(gp_dm->x[2], 0., dim[2]);
      gp_gas->x[0] = box_wrap(gp_gas->x[0], 0., dim[0]);
      gp_gas->x[1] = box_wrap(gp_gas->x[1], 0., dim[1]);
      gp_gas->x[2] = box_wrap(gp_gas->x[2], 0., dim[2]);
      p->x[0] = box_wrap(p->x[0], 0., dim[0]);
      p->x[1] = box_wrap(p->x[1], 0., dim[1]);
      p->x[2] = box_wrap(p->x[2], 0., dim[2]);
    }

    /* Also copy the velocities */
    p->v[0] = gp_gas->v_full[0];
    p->v[1] = gp_gas->v_full[1];
    p->v[2] = gp_gas->v_full[2];

    /* Set the smoothing length to the mean inter-particle separation */
    p->h = d;

    /* Note that the thermodynamic properties (u, S, ...) will be set later */
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
 * @brief #threadpool mapper function for the limiter debugging check
 */
void space_check_limiter_mapper(void *map_data, int nr_parts,
                                void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;

  /* Verify that all limited particles have been treated */
  for (int k = 0; k < nr_parts; k++) {

    if (parts[k].time_bin == time_bin_inhibited) continue;

    if (parts[k].wakeup == time_bin_awake)
      error("Particle still woken up! id=%lld", parts[k].id);

    if (parts[k].gpart != NULL)
      if (parts[k].time_bin != parts[k].gpart->time_bin)
        error("Gpart not on the same time-bin as part");
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that all particles have their wakeup flag in a correct state.
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to check.
 */
void space_check_limiter(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS

  threadpool_map(&s->e->threadpool, space_check_limiter_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), 1000, NULL);
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
  free(s->local_cells_with_tasks_top);
  free(s->cell_index);
  free(s->cells_with_particles_top);
  free(s->local_cells_with_particles_top);
  free(s->parts);
  free(s->xparts);
  free(s->gparts);
  free(s->sparts);
  free(s->fof_data.group_index);
  free(s->fof_data.group_size);
  free(s->fof_data.group_mass);
  free(s->fof_data.group_CoM);
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
  s->local_cells_with_tasks_top = NULL;
  s->cells_with_particles_top = NULL;
  s->local_cells_with_particles_top = NULL;
  s->nr_local_cells_with_tasks = 0;
  s->nr_cells_with_particles = 0;
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

  /* Need to reconnect the gravity parts to their hydro and stars particles. */
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
