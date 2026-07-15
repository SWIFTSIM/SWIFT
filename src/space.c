/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fftw3.h>
#include <time.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "black_holes.h"
#include "const.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "kernel_hydro.h"
#include "lock.h"
#include "mesh_gravity.h"
#include "mhd.h"
#include "minmax.h"
#include "proxy.h"
#include "restart.h"
#include "row_major_id.h"
#include "rt.h"
#include "sort_part.h"
#include "space_unique_id.h"
#include "star_formation.h"
#include "stars.h"
#include "threadpool.h"
#include "tools.h"
#include "tracers.h"

/* Split size. */
int space_splitsize = space_splitsize_default;
int space_subsize_pair_hydro = space_subsize_pair_hydro_default;
int space_subsize_self_hydro = space_subsize_self_hydro_default;
int space_subsize_pair_stars = space_subsize_pair_stars_default;
int space_subsize_self_stars = space_subsize_self_stars_default;
int space_subsize_pair_grav = space_subsize_pair_grav_default;
int space_subsize_self_grav = space_subsize_self_grav_default;
int space_subdepth_diff_grav = space_subdepth_diff_grav_default;
int space_maxsize = space_maxsize_default;
int space_grid_split_threshold = space_grid_split_threshold_default;

/* Recursion sizes */
int space_recurse_size_self_hydro = space_recurse_size_self_hydro_default;
int space_recurse_size_pair_hydro = space_recurse_size_pair_hydro_default;
int space_recurse_size_self_stars = space_recurse_size_self_stars_default;
int space_recurse_size_pair_stars = space_recurse_size_pair_stars_default;
int space_recurse_size_self_black_holes =
    space_recurse_size_self_black_holes_default;
int space_recurse_size_pair_black_holes =
    space_recurse_size_pair_black_holes_default;
int space_recurse_size_self_sinks = space_recurse_size_self_sinks_default;
int space_recurse_size_pair_sinks = space_recurse_size_pair_sinks_default;

/*! Number of extra #part we allocate memory for per top-level cell */
int space_extra_parts = space_extra_parts_default;

/*! Number of extra #spart we allocate memory for per top-level cell */
int space_extra_sparts = space_extra_sparts_default;

/*! Number of extra #bpart we allocate memory for per top-level cell */
int space_extra_bparts = space_extra_bparts_default;

/*! Number of extra #gpart we allocate memory for per top-level cell */
int space_extra_gparts = space_extra_gparts_default;

/*! Number of extra #sink we allocate memory for per top-level cell */
int space_extra_sinks = space_extra_sinks_default;

/*! Maximum number of particles per ghost */
int engine_max_parts_per_ghost = engine_max_parts_per_ghost_default;
int engine_max_sparts_per_ghost = engine_max_sparts_per_ghost_default;
int engine_max_parts_per_cooling = engine_max_parts_per_cooling_default;

/*! Allocation margins */
double engine_redistribute_alloc_margin =
    engine_redistribute_alloc_margin_default;
double engine_foreign_alloc_margin = engine_foreign_alloc_margin_default;

/*! Maximal depth at which the stars resort task can be pushed */
int engine_star_resort_task_depth = engine_star_resort_task_depth_default;

/*! Expected maximal number of strays received at a rebuild */
int space_expected_max_nr_strays = space_expected_max_nr_strays_default;

/*! Counter for cell IDs (when debugging + max vals for unique IDs exceeded) */
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
unsigned long long last_cell_id;
unsigned long long last_leaf_cell_id;
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
 * @brief Free any memory in use for foreign particles.
 *
 * @param s The #space.
 * @param clear_cell_pointers Are we also setting all the foreign cell particle
 * pointers to NULL?
 */
void space_free_foreign_parts(struct space *s, const int clear_cell_pointers) {

#ifdef WITH_MPI
  if (s->parts_foreign != NULL) {
    swift_free("parts_foreign", s->parts_foreign);
    s->size_parts_foreign = 0;
    s->parts_foreign = NULL;
  }
  if (s->gparts_foreign != NULL) {
    swift_free("gparts_foreign", s->gparts_foreign);
    s->size_gparts_foreign = 0;
    s->gparts_foreign = NULL;
  }
  if (s->gparts_fof_foreign != NULL) {
    swift_free("gparts_fof_foreign", s->gparts_fof_foreign);
    s->size_gparts_foreign = 0;
    s->gparts_fof_foreign = NULL;
  }
  if (s->sparts_foreign != NULL) {
    swift_free("sparts_foreign", s->sparts_foreign);
    s->size_sparts_foreign = 0;
    s->sparts_foreign = NULL;
  }
  if (s->bparts_foreign != NULL) {
    swift_free("bparts_foreign", s->bparts_foreign);
    s->size_bparts_foreign = 0;
    s->bparts_foreign = NULL;
  }
  if (s->sinks_foreign != NULL) {
    swift_free("sinks_foreign", s->sinks_foreign);
    s->size_sinks_foreign = 0;
    s->sinks_foreign = NULL;
  }
  if (clear_cell_pointers) {
    for (int k = 0; k < s->e->nr_proxies; k++) {
      for (int j = 0; j < s->e->proxies[k].nr_cells_in; j++) {
        cell_unlink_foreign_particles(s->e->proxies[k].cells_in[j]);
      }
    }
  }
#endif
}

void space_reorder_extra_parts_mapper(void *map_data, int num_cells,
                                      void *extra_data) {
  int *local_cells = (int *)map_data;
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells[ind]];
    cell_reorder_extra_parts(c, c->hydro.parts - s->parts);
  }
}

void space_reorder_extra_gparts_mapper(void *map_data, int num_cells,
                                       void *extra_data) {

  int *local_cells = (int *)map_data;
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells[ind]];
    cell_reorder_extra_gparts(c, s->parts, s->sparts, s->sinks, s->bparts);
  }
}

void space_reorder_extra_sparts_mapper(void *map_data, int num_cells,
                                       void *extra_data) {

  int *local_cells = (int *)map_data;
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells[ind]];
    cell_reorder_extra_sparts(c, c->stars.parts - s->sparts);
  }
}

void space_reorder_extra_sinks_mapper(void *map_data, int num_cells,
                                      void *extra_data) {

  int *local_cells = (int *)map_data;
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells[ind]];
    cell_reorder_extra_sinks(c, c->sinks.parts - s->sinks);
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

  /* Re-order the gas particles */
  if (space_extra_parts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_parts_mapper,
                   s->local_cells_top, s->nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, s);

  /* Re-order the gravity particles */
  if (space_extra_gparts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_gparts_mapper,
                   s->local_cells_top, s->nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, s);

  /* Re-order the star particles */
  if (space_extra_sparts)
    threadpool_map(&s->e->threadpool, space_reorder_extra_sparts_mapper,
                   s->local_cells_top, s->nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, s);

  /* Re-order the black hole particles */
  if (space_extra_bparts)
    error("Missing implementation of BH extra reordering");

  /* Re-order the sink particles */
  if (space_extra_sinks)
    threadpool_map(&s->e->threadpool, space_reorder_extra_sinks_mapper,
                   s->local_cells_top, s->nr_local_cells, sizeof(int),
                   threadpool_auto_chunk_size, s);
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
                 s->nr_cells, sizeof(struct cell), threadpool_auto_chunk_size,
                 /*extra_data=*/NULL);
}

/**
 * @brief Mapping function to free the sorted indices buffers.
 */
void space_map_clearsort(struct cell *c, void *data) {

  cell_free_hydro_sorts(c);
  cell_free_stars_sorts(c);
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
 * @brief Get a new empty (sub-)#cell.
 *
 * If there are cells in the buffer, use the one at the end of the linked list.
 * If we have no cells, allocate a new chunk of memory and pick one from there.
 *
 * @param s The #space.
 * @param nr_cells Number of #cell to pick up.
 * @param cells Array of @c nr_cells #cell pointers in which to store the
 *        new cells.
 * @param tpid ID of threadpool threadpool associated with cells_sub.
 */
void space_getcells(struct space *s, int nr_cells, struct cell **cells,
                    const short int tpid, int index, int ghost) {

  /* For each requested cell... */
  for (int j = 0; j < nr_cells; j++) {

    /* Is the cell buffer empty? */
    if (s->cells_sub[tpid] == NULL) {
      if (swift_memalign("cells_sub", (void **)&s->cells_sub[tpid], cell_align,
                         space_cellallocchunk * sizeof(struct cell)) != 0)
        error("Failed to allocate more cells.");

      /* This allocation is never correctly freed, that is ok. */
      swift_ignore_leak(s->cells_sub[tpid]);

      /* Clear the newly-allocated cells. */
      bzero(s->cells_sub[tpid], sizeof(struct cell) * space_cellallocchunk);

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->cells_sub[tpid][k].next = &s->cells_sub[tpid][k + 1];
      s->cells_sub[tpid][space_cellallocchunk - 1].next = NULL;
    }

    /* Is the multipole buffer empty? */
    if (s->with_self_gravity && s->multipoles_sub[tpid] == NULL) {
      if (swift_memalign(
              "multipoles_sub", (void **)&s->multipoles_sub[tpid],
              multipole_align,
              space_cellallocchunk * sizeof(struct gravity_tensors)) != 0)
        error("Failed to allocate more multipoles.");

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->multipoles_sub[tpid][k].next = &s->multipoles_sub[tpid][k + 1];
      s->multipoles_sub[tpid][space_cellallocchunk - 1].next = NULL;
    }

    /* Pick off the next cell. */
    cells[j+index] = s->cells_sub[tpid];
    s->cells_sub[tpid] = cells[j+index]->next;

    /* Hook the multipole */
    if (s->with_self_gravity) {
      cells[j+index]->grav.multipole = s->multipoles_sub[tpid];
      s->multipoles_sub[tpid] = cells[j+index]->grav.multipole->next;
    }
  }

  /* Unlock the space. */
  if (!ghost) atomic_add(&s->tot_cells, nr_cells);

  /* Init some things in the cell we just got. */
  for (int j = 0; j < nr_cells; j++) {
    cell_free_hydro_sorts(cells[j+index]);
    cell_free_stars_sorts(cells[j+index]);

    struct gravity_tensors *temp = cells[j+index]->grav.multipole;
    bzero(cells[j+index], sizeof(struct cell));
    cells[j+index]->grav.multipole = temp;
    cells[j+index]->nodeID = -1;
    cells[j+index]->tpid = tpid;
    if (lock_init(&cells[j+index]->hydro.lock) != 0 ||
        lock_init(&cells[j+index]->hydro.extra_sort_lock) != 0 ||
        lock_init(&cells[j+index]->grav.plock) != 0 ||
        lock_init(&cells[j+index]->grav.mlock) != 0 ||
        lock_init(&cells[j+index]->stars.lock) != 0 ||
        lock_init(&cells[j+index]->sinks.lock) != 0 ||
        lock_init(&cells[j+index]->sinks.sink_formation_lock) != 0 ||
        lock_init(&cells[j+index]->black_holes.lock) != 0 ||
        lock_init(&cells[j+index]->stars.star_formation_lock) != 0 ||
        lock_init(&cells[j+index]->grav.star_formation_lock) != 0)
      error("Failed to initialize cell spinlocks.");
  }
}

/**
 * @brief Free sort arrays in any cells in the cell buffer.
 *
 * @param s The #space.
 */
void space_free_buff_sort_indices(struct space *s) {
  for (short int tpid = 0; tpid < s->e->nr_pool_threads; ++tpid) {
    for (struct cell *finger = s->cells_sub[tpid]; finger != NULL;
         finger = finger->next) {
      cell_free_hydro_sorts(finger);
      cell_free_stars_sorts(finger);
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
        (c->black_holes.count > 0) || (c->sinks.count > 0) ||
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

void space_synchronize_part_positions_mapper(void *map_data, int nr_parts,
                                             void *extra_data) {
  /* Unpack the data */
  const struct part *parts = (struct part *)map_data;
  struct space *s = (struct space *)extra_data;
  const ptrdiff_t offset = parts - s->parts;
  const struct xpart *xparts = s->xparts + offset;

  for (int k = 0; k < nr_parts; k++) {

    /* Get the particle */
    const struct part *p = &parts[k];
    const struct xpart *xp = &xparts[k];

    /* Skip unimportant particles */
    if (p->time_bin == time_bin_not_created ||
        p->time_bin == time_bin_inhibited)
      continue;

    /* Get its gravity friend */
    struct gpart *gp = p->gpart;

#ifdef SWIFT_DEBUG_CHECKS
    if (gp == NULL) error("Unlinked particle!");
#endif

    /* Synchronize positions, velocities and masses */
    gp->x[0] = p->x[0];
    gp->x[1] = p->x[1];
    gp->x[2] = p->x[2];

    gp->v_full[0] = xp->v_full[0];
    gp->v_full[1] = xp->v_full[1];
    gp->v_full[2] = xp->v_full[2];

    gp->mass = hydro_get_mass(p);
  }
}

void space_synchronize_spart_positions_mapper(void *map_data, int nr_sparts,
                                              void *extra_data) {
  /* Unpack the data */
  const struct spart *sparts = (struct spart *)map_data;

  for (int k = 0; k < nr_sparts; k++) {

    /* Get the particle */
    const struct spart *sp = &sparts[k];

    /* Skip unimportant particles */
    if (sp->time_bin == time_bin_not_created ||
        sp->time_bin == time_bin_inhibited)
      continue;

    /* Get its gravity friend */
    struct gpart *gp = sp->gpart;

#ifdef SWIFT_DEBUG_CHECKS
    if (gp == NULL) error("Unlinked particle!");
#endif

    /* Synchronize positions, velocities and masses */
    gp->x[0] = sp->x[0];
    gp->x[1] = sp->x[1];
    gp->x[2] = sp->x[2];

    gp->v_full[0] = sp->v[0];
    gp->v_full[1] = sp->v[1];
    gp->v_full[2] = sp->v[2];

    gp->mass = sp->mass;
  }
}

void space_synchronize_bpart_positions_mapper(void *map_data, int nr_bparts,
                                              void *extra_data) {
  /* Unpack the data */
  const struct bpart *bparts = (struct bpart *)map_data;

  for (int k = 0; k < nr_bparts; k++) {

    /* Get the particle */
    const struct bpart *bp = &bparts[k];

    /* Skip unimportant particles */
    if (bp->time_bin == time_bin_not_created ||
        bp->time_bin == time_bin_inhibited)
      continue;

    /* Get its gravity friend */
    struct gpart *gp = bp->gpart;

#ifdef SWIFT_DEBUG_CHECKS
    if (gp == NULL) error("Unlinked particle!");
#endif

    /* Synchronize positions, velocities and masses */
    gp->x[0] = bp->x[0];
    gp->x[1] = bp->x[1];
    gp->x[2] = bp->x[2];

    gp->v_full[0] = bp->v[0];
    gp->v_full[1] = bp->v[1];
    gp->v_full[2] = bp->v[2];

    gp->mass = bp->mass;
  }
}

void space_synchronize_sink_positions_mapper(void *map_data, int nr_sinks,
                                             void *extra_data) {
  /* Unpack the data */
  const struct sink *sinks = (struct sink *)map_data;

  for (int k = 0; k < nr_sinks; k++) {

    /* Get the particle */
    const struct sink *sink = &sinks[k];

    /* Skip unimportant particles */
    if (sink->time_bin == time_bin_not_created ||
        sink->time_bin == time_bin_inhibited)
      continue;

    /* Get its gravity friend */
    struct gpart *gp = sink->gpart;

#ifdef SWIFT_DEBUG_CHECKS
    if (gp == NULL) error("Unlinked particle!");
#endif

    /* Synchronize positions, velocities and masses */
    gp->x[0] = sink->x[0];
    gp->x[1] = sink->x[1];
    gp->x[2] = sink->x[2];

    gp->v_full[0] = sink->v[0];
    gp->v_full[1] = sink->v[1];
    gp->v_full[2] = sink->v[2];

    gp->mass = sink->mass;
  }
}

/**
 * @brief Make sure the baryon particles are at the same position and
 * have the same velocity and mass as their #gpart friends.
 *
 * We copy the baryon particle properties to the #gpart type-by-type.
 *
 * @param s The #space.
 */
void space_synchronize_particle_positions(struct space *s) {

  const ticks tic = getticks();

  if (s->nr_gparts > 0 && s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_synchronize_part_positions_mapper,
                   s->parts, s->nr_parts, sizeof(struct part),
                   threadpool_auto_chunk_size, (void *)s);

  if (s->nr_gparts > 0 && s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, space_synchronize_spart_positions_mapper,
                   s->sparts, s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);

  if (s->nr_gparts > 0 && s->nr_bparts > 0)
    threadpool_map(&s->e->threadpool, space_synchronize_bpart_positions_mapper,
                   s->bparts, s->nr_bparts, sizeof(struct bpart),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);

  if (s->nr_gparts > 0 && s->nr_sinks > 0)
    threadpool_map(&s->e->threadpool, space_synchronize_sink_positions_mapper,
                   s->sinks, s->nr_sinks, sizeof(struct sink),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_convert_quantities_mapper(void *restrict map_data, int count,
                                     void *restrict extra_data) {
  struct space *s = (struct space *)extra_data;
  const struct cosmology *cosmo = s->e->cosmology;
  const struct hydro_props *hydro_props = s->e->hydro_properties;
  const struct pressure_floor_props *floor = s->e->pressure_floor_props;
  struct part *restrict parts = (struct part *)map_data;
  const ptrdiff_t index = parts - s->parts;
  struct xpart *restrict xparts = s->xparts + index;

  /* Loop over all the particles ignoring the extra buffer ones for on-the-fly
   * creation */
  for (int k = 0; k < count; k++) {
    if (parts[k].time_bin <= num_time_bins) {
      hydro_convert_quantities(&parts[k], &xparts[k], cosmo, hydro_props,
                               floor);
      mhd_convert_quantities(&parts[k], &xparts[k], cosmo, hydro_props);
    }
  }
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
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_convert_rt_quantities_mapper(void *restrict map_data, int count,
                                        void *restrict extra_data) {
  struct space *s = (struct space *)extra_data;
  const struct engine *restrict e = s->e;
  const int with_rt = (e->policy & engine_policy_rt);
  if (!with_rt) return;

  const struct rt_props *restrict rt_props = e->rt_props;
  const struct hydro_props *restrict hydro_props = e->hydro_properties;
  const struct phys_const *restrict phys_const = e->physical_constants;
  const struct unit_system *restrict iu = e->internal_units;
  const struct cosmology *restrict cosmo = e->cosmology;

  struct part *restrict parts = (struct part *)map_data;

  /* Loop over all the particles ignoring the extra buffer ones for on-the-fly
   * creation */
  for (int k = 0; k < count; k++) {
    if (parts[k].time_bin <= num_time_bins)
      rt_convert_quantities(&parts[k], rt_props, hydro_props, phys_const, iu,
                            cosmo);
  }
}

/**
 * @brief Calls the #part RT quantities conversion function on all particles in
 * the space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_convert_rt_quantities(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_convert_rt_quantities_mapper,
                   s->parts, s->nr_parts, sizeof(struct part),
                   threadpool_auto_chunk_size, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_post_init_parts_mapper(void *restrict map_data, int count,
                                  void *restrict extra_data) {
  struct space *s = (struct space *)extra_data;
  const struct engine *restrict e = s->e;

  const struct hydro_props *restrict hydro_props = e->hydro_properties;
  const struct phys_const *restrict phys_const = e->physical_constants;
  const struct unit_system *us = s->e->internal_units;
  const struct cosmology *restrict cosmo = e->cosmology;
  const struct cooling_function_data *cool_func = e->cooling_func;

  struct part *restrict p = (struct part *)map_data;
  const ptrdiff_t delta = p - s->parts;
  struct xpart *restrict xp = s->xparts + delta;

  /* Loop over all the particles ignoring the extra buffer ones for on-the-fly
   * creation
   * Here we can initialize the cooling properties of the (x-)particles
   * using quantities (like the density) defined only after the neighbour loop.
   *
   * */

  for (int k = 0; k < count; k++) {
    cooling_post_init_part(phys_const, us, hydro_props, cosmo, cool_func, &p[k],
                           &xp[k]);
  }
}

/**
 * @brief Calls the #part post-initialisation function on all particles in the
 * space.
 * Here we can initialize the cooling properties of the (x-)particles
 * using quantities (like the density) defined only after the initial neighbour
 * loop.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_post_init_parts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_post_init_parts_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_collect_sum_part_mass(void *restrict map_data, int count,
                                 void *restrict extra_data) {

  struct space *s = (struct space *)extra_data;
  const struct part *parts = (const struct part *)map_data;

  /* Local collection */
  double sum = 0.;
  for (int i = 0; i < count; ++i) sum += hydro_get_mass(&parts[i]);

  /* Store back */
  atomic_add_d(&s->initial_mean_mass_particles[0], sum);
  atomic_add(&s->initial_count_particles[0], count);
}

void space_collect_sum_gpart_mass(void *restrict map_data, int count,
                                  void *restrict extra_data) {

  struct space *s = (struct space *)extra_data;
  const struct gpart *gparts = (const struct gpart *)map_data;

  /* Local collection */
  double sum_DM = 0., sum_DM_background = 0., sum_nu = 0.;
  long long count_DM = 0, count_DM_background = 0, count_nu = 0;
  for (int i = 0; i < count; ++i) {
    /* Skip inexistant particles */
    if (gpart_is_inhibited(&gparts[i], s->e) ||
        gparts[i].time_bin == time_bin_not_created)
      continue;

    if (gparts[i].type == swift_type_dark_matter) {
      sum_DM += gparts[i].mass;
      count_DM++;
    }
    if (gparts[i].type == swift_type_dark_matter_background) {
      sum_DM_background += gparts[i].mass;
      count_DM_background++;
    }
    if (gparts[i].type == swift_type_neutrino) {
      sum_nu += gparts[i].mass;
      count_nu++;
    }
  }

  /* Store back */
  atomic_add_d(&s->initial_mean_mass_particles[1], sum_DM);
  atomic_add(&s->initial_count_particles[1], count_DM);
  atomic_add_d(&s->initial_mean_mass_particles[2], sum_DM_background);
  atomic_add(&s->initial_count_particles[2], count_DM_background);
  atomic_add_d(&s->initial_mean_mass_particles[6], sum_nu);
  atomic_add(&s->initial_count_particles[6], count_nu);
}

void space_collect_sum_sink_mass(void *restrict map_data, int count,
                                 void *restrict extra_data) {

  struct space *s = (struct space *)extra_data;
  const struct sink *sinks = (const struct sink *)map_data;

  /* Local collection */
  double sum = 0.;
  for (int i = 0; i < count; ++i) sum += sinks[i].mass;

  /* Store back */
  atomic_add_d(&s->initial_mean_mass_particles[3], sum);
  atomic_add(&s->initial_count_particles[3], count);
}

void space_collect_sum_spart_mass(void *restrict map_data, int count,
                                  void *restrict extra_data) {

  struct space *s = (struct space *)extra_data;
  const struct spart *sparts = (const struct spart *)map_data;

  /* Local collection */
  double sum = 0.;
  for (int i = 0; i < count; ++i) sum += sparts[i].mass;

  /* Store back */
  atomic_add_d(&s->initial_mean_mass_particles[4], sum);
  atomic_add(&s->initial_count_particles[4], count);
}

void space_collect_sum_bpart_mass(void *restrict map_data, int count,
                                  void *restrict extra_data) {

  struct space *s = (struct space *)extra_data;
  const struct bpart *bparts = (const struct bpart *)map_data;

  /* Local collection */
  double sum = 0.;
  for (int i = 0; i < count; ++i) sum += bparts[i].mass;

  /* Store back */
  atomic_add_d(&s->initial_mean_mass_particles[5], sum);
  atomic_add(&s->initial_count_particles[5], count);
}

/**
 * @brief Collect the mean mass of each particle type in the #space.
 */
void space_collect_mean_masses(struct space *s, int verbose) {

  /* Init counters */
  for (int i = 0; i < swift_type_count; ++i)
    s->initial_mean_mass_particles[i] = 0.;
  for (int i = 0; i < swift_type_count; ++i) s->initial_count_particles[i] = 0;

  /* Collect each particle type */
  threadpool_map(&s->e->threadpool, space_collect_sum_part_mass, s->parts,
                 s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                 s);
  threadpool_map(&s->e->threadpool, space_collect_sum_gpart_mass, s->gparts,
                 s->nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 s);
  threadpool_map(&s->e->threadpool, space_collect_sum_spart_mass, s->sparts,
                 s->nr_sparts, sizeof(struct spart), threadpool_auto_chunk_size,
                 s);
  threadpool_map(&s->e->threadpool, space_collect_sum_sink_mass, s->sinks,
                 s->nr_sinks, sizeof(struct sink), threadpool_auto_chunk_size,
                 s);
  threadpool_map(&s->e->threadpool, space_collect_sum_bpart_mass, s->bparts,
                 s->nr_bparts, sizeof(struct bpart), threadpool_auto_chunk_size,
                 s);

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, s->initial_mean_mass_particles, swift_type_count,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, s->initial_count_particles, swift_type_count,
                MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /* Get means
   *
   * Note: the Intel compiler vectorizes this loop and creates FPEs from
   * the masked bit of the vector... Silly ICC... */
  /* TK comment: the following also has problems with gnu_7.3.0 and optimization
   */
#if defined(__ICC)
#pragma novector
#endif
  for (int i = 0; i < swift_type_count; ++i)
    if (s->initial_count_particles[i] > 0)
      s->initial_mean_mass_particles[i] /=
          (double)s->initial_count_particles[i];
}

/**
 * @brief Split the space into cells given the array of particles.
 *
 * @param s The #space to initialize.
 * @param params The parsed parameter file.
 * @param cosmo The current cosmological model.
 * @param dim Spatial dimensions of the domain.
 * @param hydro_properties The properties of the hydro scheme.
 * @param parts Array of Gas particles.
 * @param gparts Array of Gravity particles.
 * @param sinks Array of sink particles.
 * @param sparts Array of stars particles.
 * @param bparts Array of black hole particles.
 * @param Npart The number of Gas particles in the space.
 * @param Ngpart The number of Gravity particles in the space.
 * @param Nsink The number of sink particles in the space.
 * @param Nspart The number of stars particles in the space.
 * @param Nbpart The number of black hole particles in the space.
 * @param periodic flag whether the domain is periodic or not.
 * @param replicate How many replications along each direction do we want?
 * @param remap_ids Are we remapping the IDs from 1 to N?
 * @param generate_gas_in_ics Are we generating gas particles from the gparts?
 * @param hydro flag whether we are doing hydro or not?
 * @param self_gravity flag whether we are doing gravity or not?
 * @param star_formation flag whether we are doing star formation or not?
 * @param with_sink flag whether we are doing sink particles or not?
 * @param DM_background Are we running with some DM background particles?
 * @param verbose Print messages to stdout or not.
 * @param dry_run If 1, just initialise stuff, don't do anything with the parts.
 * @param nr_nodes The number of MPI rank.
 *
 * Makes a grid of edge length > r_max and fills the particles
 * into the respective cells. Cells containing more than #space_splitsize
 * parts with a cutoff below half the cell width are then split
 * recursively.
 */
void space_init(struct space *s, struct swift_params *params,
                const struct cosmology *cosmo, double dim[3],
                const struct hydro_props *hydro_properties, struct part *parts,
                struct gpart *gparts, struct sink *sinks, struct spart *sparts,
                struct bpart *bparts, size_t Npart, size_t Ngpart, size_t Nsink,
                size_t Nspart, size_t Nbpart, size_t Nnupart, int periodic,
                int replicate, int remap_ids, int generate_gas_in_ics,
                int hydro, int self_gravity, int star_formation, int with_sink,
                int with_DM, int with_DM_background, int neutrinos, int verbose,
                int dry_run, int nr_nodes) {

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
  s->with_sink = with_sink;
  s->with_DM = with_DM;
  s->with_DM_background = with_DM_background;
  s->with_neutrinos = neutrinos;
  s->nr_parts = Npart;
  s->nr_gparts = Ngpart;
  s->nr_sparts = Nspart;
  s->nr_bparts = Nbpart;
  s->nr_sinks = Nsink;
  s->nr_nuparts = Nnupart;
  s->size_parts = Npart;
  s->size_gparts = Ngpart;
  s->size_sparts = Nspart;
  s->size_bparts = Nbpart;
  s->size_sinks = Nsink;
  s->nr_inhibited_parts = 0;
  s->nr_inhibited_gparts = 0;
  s->nr_inhibited_sparts = 0;
  s->nr_inhibited_bparts = 0;
  s->nr_inhibited_sinks = 0;
  s->nr_extra_parts = 0;
  s->nr_extra_gparts = 0;
  s->nr_extra_sparts = 0;
  s->nr_extra_bparts = 0;
  s->nr_extra_sinks = 0;
  s->parts = parts;
  s->gparts = gparts;
  s->sparts = sparts;
  s->bparts = bparts;
  s->sinks = sinks;
  s->min_part_mass = FLT_MAX;
  s->min_gpart_mass = FLT_MAX;
  s->min_sink_mass = FLT_MAX;
  s->min_spart_mass = FLT_MAX;
  s->min_bpart_mass = FLT_MAX;
  s->sum_part_vel_norm = 0.f;
  s->sum_gpart_vel_norm = 0.f;
  s->sum_sink_vel_norm = 0.f;
  s->sum_spart_vel_norm = 0.f;
  s->sum_bpart_vel_norm = 0.f;
  s->nr_queues = 1; /* Temporary value until engine construction */

  /* do a quick check that the box size has valid values */
#if defined HYDRO_DIMENSION_1D
  if (dim[0] <= 0.) error("Invalid box size: [%f]", dim[0]);
#elif defined HYDRO_DIMENSION_2D
  if (dim[0] <= 0. || dim[1] <= 0.)
    error("Invalid box size: [%f, %f]", dim[0], dim[1]);
#else
  if (dim[0] <= 0. || dim[1] <= 0. || dim[2] <= 0.)
    error("Invalid box size: [%f, %f, %f]", dim[0], dim[1], dim[2]);
#endif

  /* Initiate some basic randomness */
  srand(42);

  /* Are we remapping the IDs to the range [1, NumPart]? */
  if (remap_ids) {
    space_remap_ids(s, nr_nodes, verbose);
  }

  /* Are we generating gas from the DM-only ICs? */
  if (generate_gas_in_ics) {
    space_generate_gas(s, cosmo, hydro_properties, periodic, with_DM_background,
                       neutrinos, dim, verbose);
    parts = s->parts;
    gparts = s->gparts;
    Npart = s->nr_parts;
    Ngpart = s->nr_gparts;

#ifdef SWIFT_DEBUG_CHECKS
    if (!dry_run)
      part_verify_links(parts, gparts, sinks, sparts, bparts, Npart, Ngpart,
                        Nsink, Nspart, Nbpart, 1);
#endif
  }

  /* Are we replicating the space ? */
  if (replicate < 1)
    error("Value of 'InitialConditions:replicate' (%d) is too small",
          replicate);
  if (replicate > 1) {
    if (with_DM_background)
      error("Can't replicate the space if background DM particles are in use.");

    space_replicate(s, replicate, verbose);
    parts = s->parts;
    gparts = s->gparts;
    sparts = s->sparts;
    bparts = s->bparts;
    sinks = s->sinks;
    Npart = s->nr_parts;
    Ngpart = s->nr_gparts;
    Nspart = s->nr_sparts;
    Nbpart = s->nr_bparts;
    Nsink = s->nr_sinks;

#ifdef SWIFT_DEBUG_CHECKS
    part_verify_links(parts, gparts, sinks, sparts, bparts, Npart, Ngpart,
                      Nsink, Nspart, Nbpart, 1);
#endif
  }

  /* Decide on the minimal top-level cell size */
  const double dmax = max3(s->dim[0], s->dim[1], s->dim[2]);
  int maxtcells =
      parser_get_opt_param_int(params, "Scheduler:max_top_level_cells",
                               space_max_top_level_cells_default);
  /* We're at risk of rounding errors if tol ~ 1 - 1/maxtcells. */
  /* Make sure it's (much) closer to 1.0 than this. */
  /* But also ensure it's not so small that we round in the other direction */
  const float tol = max(1.0 - 1.0 / (maxtcells * maxtcells), 0.99);
  s->cell_min = tol * dmax / maxtcells;

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
  space_subsize_pair_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_pair_stars",
                               space_subsize_pair_stars_default);
  space_subsize_self_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_stars",
                               space_subsize_self_stars_default);
  space_subsize_pair_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_pair_grav",
                               space_subsize_pair_grav_default);
  space_subsize_self_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_sub_size_self_grav",
                               space_subsize_self_grav_default);
  space_splitsize = parser_get_opt_param_int(
      params, "Scheduler:cell_split_size", space_splitsize_default);
  space_grid_split_threshold = parser_get_opt_param_int(
      params, "Scheduler:grid_split_threshold", space_grid_split_threshold);
  space_subdepth_diff_grav =
      parser_get_opt_param_int(params, "Scheduler:cell_subdepth_diff_grav",
                               space_subdepth_diff_grav_default);

  space_recurse_size_self_hydro =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_self_hydro",
                               space_recurse_size_self_hydro_default);
  space_recurse_size_pair_hydro =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_pair_hydro",
                               space_recurse_size_pair_hydro_default);
  space_recurse_size_self_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_self_stars",
                               space_recurse_size_self_stars_default);
  space_recurse_size_pair_stars =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_pair_stars",
                               space_recurse_size_pair_stars_default);
  space_recurse_size_self_black_holes = parser_get_opt_param_int(
      params, "Scheduler:cell_recurse_size_self_black_holes",
      space_recurse_size_self_black_holes_default);
  space_recurse_size_pair_black_holes = parser_get_opt_param_int(
      params, "Scheduler:cell_recurse_size_pair_black_holes",
      space_recurse_size_pair_black_holes_default);
  space_recurse_size_self_sinks =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_self_sinks",
                               space_recurse_size_self_sinks_default);
  space_recurse_size_pair_sinks =
      parser_get_opt_param_int(params, "Scheduler:cell_recurse_size_pair_sinks",
                               space_recurse_size_pair_sinks_default);

  space_extra_parts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_parts", space_extra_parts_default);
  space_extra_sparts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_sparts", space_extra_sparts_default);
  space_extra_gparts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_gparts", space_extra_gparts_default);
  space_extra_bparts = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_bparts", space_extra_bparts_default);
  space_extra_sinks = parser_get_opt_param_int(
      params, "Scheduler:cell_extra_sinks", space_extra_sinks_default);

  engine_max_parts_per_ghost =
      parser_get_opt_param_int(params, "Scheduler:engine_max_parts_per_ghost",
                               engine_max_parts_per_ghost_default);
  engine_max_sparts_per_ghost =
      parser_get_opt_param_int(params, "Scheduler:engine_max_sparts_per_ghost",
                               engine_max_sparts_per_ghost_default);

  engine_max_parts_per_cooling =
      parser_get_opt_param_int(params, "Scheduler:engine_max_parts_per_cooling",
                               engine_max_parts_per_cooling_default);

  engine_redistribute_alloc_margin = parser_get_opt_param_double(
      params, "Scheduler:engine_redist_alloc_margin",
      engine_redistribute_alloc_margin_default);

  engine_foreign_alloc_margin = parser_get_opt_param_double(
      params, "Scheduler:engine_foreign_alloc_margin",
      engine_foreign_alloc_margin_default);

  if (verbose) {
    message("max_size set to %d split_size set to %d", space_maxsize,
            space_splitsize);
    message("subdepth_grav set to %d", space_subdepth_diff_grav);
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

  /* Read in imposed star smoothing length */
  s->initial_spart_h = parser_get_opt_param_float(
      params, "InitialConditions:stars_smoothing_length", -1.f);
  if (s->initial_spart_h != -1.f) {
    message("Imposing a star smoothing length of %e", s->initial_spart_h);
  }
  /* Read in imposed black hole smoothing length */
  s->initial_bpart_h = parser_get_opt_param_float(
      params, "InitialConditions:black_holes_smoothing_length", -1.f);
  if (s->initial_bpart_h != -1.f) {
    message("Imposing a BH smoothing length of %e", s->initial_bpart_h);
  }

  /* Apply shift */
  double shift[3] = {0.0, 0.0, 0.0};
  parser_get_opt_param_double_array(params, "InitialConditions:shift", 3,
                                    shift);
  memcpy(s->initial_shift, shift, 3 * sizeof(double));
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
    for (size_t k = 0; k < Nbpart; k++) {
      bparts[k].x[0] += shift[0];
      bparts[k].x[1] += shift[1];
      bparts[k].x[2] += shift[2];
    }
    for (size_t k = 0; k < Nsink; k++) {
      sinks[k].x[0] += shift[0];
      sinks[k].x[1] += shift[1];
      sinks[k].x[2] += shift[2];
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

    /* Same for the bparts */
    if (periodic) {
      for (size_t k = 0; k < Nbpart; k++)
        for (int j = 0; j < 3; j++) {
          while (bparts[k].x[j] < 0) bparts[k].x[j] += s->dim[j];
          while (bparts[k].x[j] >= s->dim[j]) bparts[k].x[j] -= s->dim[j];
        }
    } else {
      for (size_t k = 0; k < Nbpart; k++)
        for (int j = 0; j < 3; j++)
          if (bparts[k].x[j] < 0 || bparts[k].x[j] >= s->dim[j])
            error("Not all b-particles are within the specified domain.");
    }

    /* Same for the sinks */
    if (periodic) {
      for (size_t k = 0; k < Nsink; k++)
        for (int j = 0; j < 3; j++) {
          while (sinks[k].x[j] < 0) sinks[k].x[j] += s->dim[j];
          while (sinks[k].x[j] >= s->dim[j]) sinks[k].x[j] -= s->dim[j];
        }
    } else {
      for (size_t k = 0; k < Nsink; k++)
        for (int j = 0; j < 3; j++)
          if (sinks[k].x[j] < 0 || sinks[k].x[j] >= s->dim[j])
            error("Not all sink-particles are within the specified domain.");
    }
  }

  /* Allocate the extra parts array for the gas particles. */
  if (Npart > 0) {
    if (swift_memalign("xparts", (void **)&s->xparts, xpart_align,
                       Npart * sizeof(struct xpart)) != 0)
      error("Failed to allocate xparts.");
    bzero(s->xparts, Npart * sizeof(struct xpart));
  }

  hydro_space_init(&s->hs, s);

  /* Init the space lock. */
  if (lock_init(&s->lock) != 0) error("Failed to create space spin-lock.");

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  last_cell_id = 1ULL;
  last_leaf_cell_id = 1ULL;
#endif

  /* Do we want any spare particles for on the fly creation?
     This condition should be the same than in engine_config.c */
  if (!(star_formation || with_sink) ||
      !swift_star_formation_model_creates_stars) {
    space_extra_sparts = 0;
    space_extra_gparts = 0;
    space_extra_sinks = 0;
  }

  const int create_sparts =
      (star_formation && swift_star_formation_model_creates_stars) || with_sink;
  if (create_sparts && space_extra_sparts == 0) {
    error(
        "Running with star formation but without spare star particles. "
        "Increase 'Scheduler:cell_extra_sparts'.");
  }

  if (with_sink && space_extra_gparts == 0) {
    error(
        "Running with star formation from sink but without spare g-particles. "
        "Increase 'Scheduler:cell_extra_gparts'.");
  }
  if (with_sink && space_extra_sinks == 0) {
    error(
        "Running with star formation from sink but without spare "
        "sink-particles. "
        "Increase 'Scheduler:cell_extra_sinks'.");
  }

  /* Build the cells recursively. */
  //if (!dry_run) space_regrid(s, verbose);
  space_regrid(s, verbose);

  /* Compute the max id for the generation of unique id. */
  if (create_sparts) {
    space_init_unique_id(s, nr_nodes);
  }
  
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
  const size_t nr_bparts = s->nr_bparts;
  const size_t nr_sinks = s->nr_sinks;
  const size_t nr_nuparts = s->nr_nuparts;
  const size_t nr_dm = nr_gparts - nr_parts - nr_sparts - nr_bparts;

  s->size_parts = s->nr_parts = nr_parts * factor;
  s->size_gparts = s->nr_gparts = nr_gparts * factor;
  s->size_sparts = s->nr_sparts = nr_sparts * factor;
  s->size_bparts = s->nr_bparts = nr_bparts * factor;
  s->size_sinks = s->nr_sinks = nr_sinks * factor;
  s->nr_nuparts = nr_nuparts * factor;

  /* Allocate space for new particles */
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;
  struct sink *sinks = NULL;

  if (swift_memalign("parts", (void **)&parts, part_align,
                     s->nr_parts * sizeof(struct part)) != 0)
    error("Failed to allocate new part array.");

  if (swift_memalign("gparts", (void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate new gpart array.");

  if (swift_memalign("sparts", (void **)&sparts, spart_align,
                     s->nr_sparts * sizeof(struct spart)) != 0)
    error("Failed to allocate new spart array.");

  if (swift_memalign("sinks", (void **)&sinks, sink_align,
                     s->nr_sinks * sizeof(struct sink)) != 0)
    error("Failed to allocate new sink array.");

  if (swift_memalign("bparts", (void **)&bparts, bpart_align,
                     s->nr_bparts * sizeof(struct bpart)) != 0)
    error("Failed to allocate new bpart array.");

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
        memcpy(bparts + offset * nr_bparts, s->bparts,
               nr_bparts * sizeof(struct bpart));
        memcpy(gparts + offset * nr_gparts, s->gparts,
               nr_gparts * sizeof(struct gpart));
        memcpy(sinks + offset * nr_sinks, s->sinks,
               nr_sinks * sizeof(struct sink));

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
        for (size_t n = offset * nr_bparts; n < (offset + 1) * nr_bparts; ++n) {
          bparts[n].x[0] += shift[0];
          bparts[n].x[1] += shift[1];
          bparts[n].x[2] += shift[2];
        }
        for (size_t n = offset * nr_sinks; n < (offset + 1) * nr_sinks; ++n) {
          sinks[n].x[0] += shift[0];
          sinks[n].x[1] += shift[1];
          sinks[n].x[2] += shift[2];
        }

        /* Set the correct links (recall gpart are sorted by type at start-up):
           first DM (unassociated gpart), then gas, then sinks, then stars */
        if (nr_parts > 0 && nr_gparts > 0) {
          const size_t offset_part = offset * nr_parts;
          const size_t offset_gpart = offset * nr_gparts + nr_dm;

          for (size_t n = 0; n < nr_parts; ++n) {
            parts[offset_part + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_part + n);
          }
        }
        if (nr_sinks > 0 && nr_gparts > 0) {
          const size_t offset_sink = offset * nr_sinks;
          const size_t offset_gpart = offset * nr_gparts + nr_dm + nr_parts;

          for (size_t n = 0; n < nr_sinks; ++n) {
            sinks[offset_sink + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_sink + n);
          }
        }
        if (nr_sparts > 0 && nr_gparts > 0) {
          const size_t offset_spart = offset * nr_sparts;
          const size_t offset_gpart =
              offset * nr_gparts + nr_dm + nr_parts + nr_sinks;

          for (size_t n = 0; n < nr_sparts; ++n) {
            sparts[offset_spart + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_spart + n);
          }
        }
        if (nr_bparts > 0 && nr_gparts > 0) {
          const size_t offset_bpart = offset * nr_bparts;
          const size_t offset_gpart =
              offset * nr_gparts + nr_dm + nr_parts + nr_sinks + nr_sparts;

          for (size_t n = 0; n < nr_bparts; ++n) {
            bparts[offset_bpart + n].gpart = &gparts[offset_gpart + n];
            gparts[offset_gpart + n].id_or_neg_offset = -(offset_bpart + n);
          }
        }
      }
    }
  }

  /* Replace the content of the space */
  swift_free("parts", s->parts);
  swift_free("gparts", s->gparts);
  swift_free("sparts", s->sparts);
  swift_free("bparts", s->bparts);
  swift_free("sinks", s->sinks);
  s->parts = parts;
  s->gparts = gparts;
  s->sparts = sparts;
  s->bparts = bparts;
  s->sinks = sinks;

  /* Finally, update the domain size */
  s->dim[0] *= replicate;
  s->dim[1] *= replicate;
  s->dim[2] *= replicate;

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that everything is correct */
  part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                    s->nr_parts, s->nr_gparts, s->nr_sinks, s->nr_sparts,
                    s->nr_bparts, verbose);
#endif
}

/**
 * @brief Remaps the IDs of the particles to the range [1, N]
 *
 * The IDs are unique accross all MPI ranks and are generated
 * in ther order DM, gas, sinks, stars, BHs.
 *
 * @param s The current #space object.
 * @param nr_nodes The number of MPI ranks used in the run.
 * @param verbose Are we talkative?
 */
void space_remap_ids(struct space *s, int nr_nodes, int verbose) {

  if (verbose) message("Remapping all the IDs");

  long long local_nr_dm_background = 0;
  long long local_nr_nuparts = 0;
  for (size_t i = 0; i < s->nr_gparts; ++i) {
    if (s->gparts[i].type == swift_type_neutrino)
      local_nr_nuparts++;
    else if (s->gparts[i].type == swift_type_dark_matter_background)
      local_nr_dm_background++;
  }

  /* Get the current local number of particles */
  long long local_nr_parts = s->nr_parts;
  long long local_nr_sinks = s->nr_sinks;
  long long local_nr_gparts = s->nr_gparts;
  long long local_nr_sparts = s->nr_sparts;
  long long local_nr_bparts = s->nr_bparts;
  long long local_nr_baryons =
      local_nr_parts + local_nr_sinks + local_nr_sparts + local_nr_bparts;
  long long local_nr_dm = local_nr_gparts > 0
                              ? local_nr_gparts - local_nr_baryons -
                                    local_nr_nuparts - local_nr_dm_background
                              : 0;

  /* Get the global offsets */
  long long offset_parts = 0;
  long long offset_sinks = 0;
  long long offset_sparts = 0;
  long long offset_bparts = 0;
  long long offset_dm = 0;
  long long offset_dm_background = 0;
  long long offset_nuparts = 0;
#ifdef WITH_MPI
  MPI_Exscan(&local_nr_parts, &offset_parts, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_sinks, &offset_sinks, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_sparts, &offset_sparts, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_bparts, &offset_bparts, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_dm, &offset_dm, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_dm_background, &offset_dm_background, 1,
             MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan(&local_nr_nuparts, &offset_nuparts, 1, MPI_LONG_LONG_INT, MPI_SUM,
             MPI_COMM_WORLD);
#endif

  /* Total number of particles of each kind */
  long long total_dm = offset_dm + local_nr_dm;
  long long total_parts = offset_parts + local_nr_parts;
  long long total_sinks = offset_sinks + local_nr_sinks;
  long long total_sparts = offset_sparts + local_nr_sparts;
  long long total_bparts = offset_bparts + local_nr_bparts;
  long long total_nuparts = offset_nuparts + local_nr_nuparts;
  // long long total_dm_backgroud = offset_dm_background +
  // local_nr_dm_background;

#ifdef WITH_MPI
  /* The last rank now has the correct total, let's broadcast this back */
  MPI_Bcast(&total_dm, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  MPI_Bcast(&total_parts, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  MPI_Bcast(&total_sinks, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  MPI_Bcast(&total_sparts, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  MPI_Bcast(&total_bparts, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  MPI_Bcast(&total_nuparts, 1, MPI_LONG_LONG_INT, nr_nodes - 1, MPI_COMM_WORLD);
  // MPI_Bcast(&total_dm_background, 1, MPI_LONG_LONG_INT, nr_nodes - 1,
  // MPI_COMM_WORLD);
#endif

  /* Let's order the particles
   * IDs will be DM then gas then sinks than stars then BHs then nus then
   * DM background. Note that we leave a large gap (10x the number of particles)
   * in-between the regular particles and the background ones. This allow for
   * particle splitting to keep a compact set of ids. */
  offset_dm += 1;
  offset_parts += 1 + total_dm;
  offset_sinks += 1 + total_dm + total_parts;
  offset_sparts += 1 + total_dm + total_parts + total_sinks;
  offset_bparts += 1 + total_dm + total_parts + total_sinks + total_sparts;
  offset_nuparts +=
      1 + total_dm + total_parts + total_sinks + total_sparts + total_bparts;
  offset_dm_background +=
      1 + 10 * (total_dm + total_parts + total_sinks + total_sparts +
                total_bparts + total_nuparts);

  /* We can now remap the IDs in the range [offset offset + local_nr] */
  for (long long i = 0; i < local_nr_parts; ++i) {
    s->parts[i].id = offset_parts + i;
  }
  for (long long i = 0; i < local_nr_sinks; ++i) {
    s->sinks[i].id = offset_sinks + i;
  }
  for (long long i = 0; i < local_nr_sparts; ++i) {
    s->sparts[i].id = offset_sparts + i;
  }
  for (long long i = 0; i < local_nr_bparts; ++i) {
    s->bparts[i].id = offset_bparts + i;
  }
  long long count_dm = 0;
  long long count_dm_background = 0;
  long long count_nu = 0;
  for (size_t i = 0; i < s->nr_gparts; ++i) {
    if (s->gparts[i].type == swift_type_dark_matter) {
      s->gparts[i].id_or_neg_offset = offset_dm + count_dm;
      count_dm++;
    } else if (s->gparts[i].type == swift_type_neutrino) {
      s->gparts[i].id_or_neg_offset = offset_nuparts + count_nu;
      count_nu++;
    } else if (s->gparts[i].type == swift_type_dark_matter_background) {
      s->gparts[i].id_or_neg_offset =
          offset_dm_background + count_dm_background;
      count_dm_background++;
    }
  }
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
 * Background DM particles are not duplicated.
 *
 * @param s The #space to create the particles in.
 * @param cosmo The current #cosmology model.
 * @param hydro_properties The properties of the hydro scheme.
 * @param periodic Are we using periodic boundary conditions?
 * @param with_background Are we using background DM particles?
 * @param dim The size of the box (for periodic wrapping).
 * @param verbose Are we talkative?
 */
void space_generate_gas(struct space *s, const struct cosmology *cosmo,
                        const struct hydro_props *hydro_properties,
                        const int periodic, const int with_background,
                        const int with_neutrinos, const double dim[3],
                        const int verbose) {

  /* Check that this is a sensible ting to do */
  if (!s->with_hydro)
    error(
        "Cannot generate gas from ICs if we are running without "
        "hydrodynamics. Need to run with -s and the corresponding "
        "hydrodynamics parameters in the YAML file.");

  if (hydro_properties->initial_internal_energy == 0.)
    error(
        "Cannot generate gas from ICs if the initial temperature is set to 0. "
        "Need to set 'SPH:initial_temperature' to a sensible value.");

  if (cosmo->Omega_b == 0.)
    error("Cannot generate gas from ICs if Omega_b is set to 0.");

  if (verbose) message("Generating gas particles from gparts");

  /* Store the current values */
  const size_t current_nr_parts = s->nr_parts;
  const size_t current_nr_gparts = s->nr_gparts;

  /* Basic checks for unwanted modes */
  if (current_nr_parts != 0)
    error("Generating gas particles from DM but gas already exist!");

  if (s->nr_sparts != 0)
    error("Generating gas particles from DM but stars already exist!");

  if (s->nr_bparts != 0)
    error("Generating gas particles from DM but BHs already exist!");

  if (s->nr_sinks != 0)
    error("Generating gas particles from DM but sinks already exist!");

  /* Pull out information about particle splitting */
  const int particle_splitting = hydro_properties->particle_splitting;
  const float splitting_mass_threshold =
      hydro_properties->particle_splitting_mass_threshold;

  /* Start by counting the number of background, neutrino & zoom DM particles */
  size_t nr_background_gparts = 0;
  size_t nr_neutrino_gparts = 0;
  if (with_background) {
    for (size_t i = 0; i < current_nr_gparts; ++i)
      if (s->gparts[i].type == swift_type_dark_matter_background)
        ++nr_background_gparts;
  }
  if (with_neutrinos) {
    for (size_t i = 0; i < current_nr_gparts; ++i)
      if (s->gparts[i].type == swift_type_neutrino) ++nr_neutrino_gparts;
  }
  const size_t nr_zoom_gparts =
      current_nr_gparts - nr_background_gparts - nr_neutrino_gparts;

  if (nr_zoom_gparts == 0)
    error("Can't generate gas from ICs if there are no high res. particles");

  /* New particle counts after replication */
  s->size_parts = s->nr_parts = nr_zoom_gparts;
  s->size_gparts = s->nr_gparts =
      2 * nr_zoom_gparts + nr_background_gparts + nr_neutrino_gparts;

  /* Allocate space for new particles */
  struct part *parts = NULL;
  struct gpart *gparts = NULL;

  if (swift_memalign("parts", (void **)&parts, part_align,
                     s->nr_parts * sizeof(struct part)) != 0)
    error("Failed to allocate new part array.");

  if (swift_memalign("gparts", (void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate new gpart array.");

  /* And zero the parts */
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));
  bzero(parts, s->nr_parts * sizeof(struct part));

  /* Compute some constants */
  const double Omega_m = cosmo->Omega_cdm + cosmo->Omega_b;
  const double mass_ratio = cosmo->Omega_b / Omega_m;
  const double bg_density = Omega_m * cosmo->critical_density_0;
  const double bg_density_inv = 1. / bg_density;

  // message("%zd", current_nr_gparts);

  /* Update the particle properties */
  size_t j = 0;
  for (size_t i = 0; i < current_nr_gparts; ++i) {

    /* For the neutrino DM particles, just copy the data */
    if (s->gparts[i].type == swift_type_neutrino) {

      memcpy(&gparts[i], &s->gparts[i], sizeof(struct gpart));

      /* For the background DM particles, copy the data and give a better ID */
    } else if (s->gparts[i].type == swift_type_dark_matter_background) {

      memcpy(&gparts[i], &s->gparts[i], sizeof(struct gpart));

      /* Multiply the ID by two to match the convention of even IDs for DM. */
      gparts[i].id_or_neg_offset *= 2;

    } else {

      /* For the zoom DM particles, there is a lot of work to do */

      struct part *p = &parts[j];
      struct gpart *gp_gas = &gparts[current_nr_gparts + j];
      struct gpart *gp_dm = &gparts[i];

      /* Start by copying over the gpart */
      memcpy(gp_gas, &s->gparts[i], sizeof(struct gpart));
      memcpy(gp_dm, &s->gparts[i], sizeof(struct gpart));

      /* Set the IDs */
      p->id = gp_gas->id_or_neg_offset * 2 + 1;
      gp_dm->id_or_neg_offset *= 2;

      if (gp_dm->id_or_neg_offset < 0)
        error("DM particle ID overflowd (DM id=%lld gas id=%lld)",
              gp_dm->id_or_neg_offset, p->id);

      if (p->id < 0) error("gas particle ID overflowd (id=%lld)", p->id);

      /* Set the links correctly */
      p->gpart = gp_gas;
      gp_gas->id_or_neg_offset = -j;
      gp_gas->type = swift_type_gas;

      /* Compute positions shift */
      const double d = cbrt(gp_dm->mass * bg_density_inv);
      const double shift_dm = 0.5 * d * mass_ratio;
      const double shift_gas = 0.5 * d * (1. - mass_ratio);

      /* Set the masses */
      gp_dm->mass *= (1. - mass_ratio);
      gp_gas->mass *= mass_ratio;
      hydro_set_mass(p, gp_gas->mass);

      /* Verify that we are not generating a gas particle larger than the
         threashold for particle splitting */
      if (particle_splitting && gp_gas->mass > splitting_mass_threshold)
        error("Generating a gas particle above the threshold for splitting");

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

      /* Move on to the next free gas slot */
      ++j;
    }
  }

  /* Replace the content of the space */
  swift_free("gparts", s->gparts);
  s->parts = parts;
  s->gparts = gparts;
}

/**
 * @brief Verify that the matter content matches the cosmology model.
 *
 * @param s The #space.
 * @param cosmo The current cosmology model.
 * @param with_hydro Are we running with hydro switched on?
 * @param rank The MPI rank of this #space.
 * @param check_neutrinos Should neutrino masses be checked?
 */
void space_check_cosmology(struct space *s, const struct cosmology *cosmo,
                           const int with_hydro, const int rank,
                           const int check_neutrinos) {

  struct gpart *gparts = s->gparts;
  const size_t nr_gparts = s->nr_gparts;

  /* Sum up the mass in this space */
  int has_background_particles = 0;
  double mass_cdm = 0.;
  double mass_b = 0.;
  double mass_nu = 0.;
  for (size_t i = 0; i < nr_gparts; ++i) {

    /* Skip extra particles */
    if (gparts[i].time_bin == time_bin_not_created) continue;

    switch (gparts[i].type) {
      case swift_type_dark_matter:
      case swift_type_dark_matter_background:
        mass_cdm += gparts[i].mass;
        break;
      case swift_type_neutrino:
        mass_nu += gparts[i].mass;
        break;
      case swift_type_gas:
      case swift_type_stars:
      case swift_type_black_hole:
      case swift_type_sink:
        mass_b += gparts[i].mass;
        break;
      default:
        error("Invalid particle type");
    }

    if (gparts[i].type == swift_type_dark_matter_background)
      has_background_particles = 1;
  }

/* Reduce the total mass */
#ifdef WITH_MPI
  double total_mass_cdm;
  double total_mass_b;
  double total_mass_nu;
  MPI_Reduce(&mass_cdm, &total_mass_cdm, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&mass_b, &total_mass_b, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass_nu, &total_mass_nu, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
#else
  double total_mass_cdm = mass_cdm;
  double total_mass_b = mass_b;
  double total_mass_nu = mass_nu;
#endif

  if (rank == 0) {

    const double volume = s->dim[0] * s->dim[1] * s->dim[2];

    /* Current Hubble constant */
    const double H = cosmo->H;

    /* z=0 Hubble parameter */
    const double H0 = cosmo->H0;

    /* Critical density at z=0 */
    const double rho_crit0 = cosmo->critical_density * H0 * H0 / (H * H);

    /* Compute the mass densities */
    const double Omega_particles_cdm = (total_mass_cdm / volume) / rho_crit0;
    const double Omega_particles_b = (total_mass_b / volume) / rho_crit0;
    const double Omega_particles_nu = (total_mass_nu / volume) / rho_crit0;

    const double Omega_particles_m = Omega_particles_cdm + Omega_particles_b;

    /* Expected matter density */
    const double Omega_m = cosmo->Omega_cdm + cosmo->Omega_b;

    if (with_hydro && !has_background_particles &&
        fabs(Omega_particles_cdm - cosmo->Omega_cdm) > 1e-3)
      error(
          "The cold dark matter content of the simulation does not match the "
          "cosmology in the parameter file: cosmo.Omega_cdm = %e particles "
          "Omega_cdm = %e",
          cosmo->Omega_cdm, Omega_particles_cdm);

    if (with_hydro && !has_background_particles &&
        fabs(Omega_particles_b - cosmo->Omega_b) > 1e-3)
      error(
          "The baryon content of the simulation does not match the cosmology "
          "in the parameter file: cosmo.Omega_b = %e particles Omega_b = %e",
          cosmo->Omega_b, Omega_particles_b);

    if (check_neutrinos && fabs(Omega_particles_nu - cosmo->Omega_nu_0) > 1e-3)
      error(
          "The massive neutrino content of the simulation does not match the "
          "cosmology in the parameter file: cosmo.Omega_nu = %e particles "
          "Omega_nu = %e",
          cosmo->Omega_nu_0, Omega_particles_nu);

    if (fabs(Omega_particles_m - Omega_m) > 1e-3)
      error(
          "The total matter content of the simulation does not match the "
          "cosmology in the parameter file: cosmo.Omega_m = %e particles "
          "Omega_m = %e \n cosmo: Omega_b=%e Omega_cdm=%e \n "
          "particles: Omega_b=%e Omega_cdm=%e",
          Omega_m, Omega_particles_m, cosmo->Omega_b, cosmo->Omega_cdm,
          Omega_particles_b, Omega_particles_cdm);
  }
}

/**
 * @brief Compute the max id of any #part in this space.
 *
 * This function is inefficient. Don't call often.
 * Background particles are ignored.
 *
 * @param s The #space.
 */
long long space_get_max_parts_id(struct space *s) {

  long long max_id = -1;
  for (size_t i = 0; i < s->nr_parts; ++i) max_id = max(max_id, s->parts[i].id);
  for (size_t i = 0; i < s->nr_sinks; ++i) max_id = max(max_id, s->sinks[i].id);
  for (size_t i = 0; i < s->nr_sparts; ++i)
    max_id = max(max_id, s->sparts[i].id);
  for (size_t i = 0; i < s->nr_bparts; ++i)
    max_id = max(max_id, s->bparts[i].id);

  /* Note: We Explicitly do *NOT* consider background particles */
  for (size_t i = 0; i < s->nr_gparts; ++i)
    if (s->gparts[i].type == swift_type_dark_matter ||
        s->gparts[i].type == swift_type_neutrino)
      max_id = max(max_id, s->gparts[i].id_or_neg_offset);
  return max_id;
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
  space_map_cells_pre(s, 1, cell_check_spart_drift_point, &ti_drift);
  space_map_cells_pre(s, 1, cell_check_bpart_drift_point, &ti_drift);
  space_map_cells_pre(s, 1, cell_check_sink_drift_point, &ti_drift);
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
void space_check_timesteps(const struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < s->nr_cells; ++i) {
    if (s->cells_top[i].nodeID == engine_rank) {
      cell_check_timesteps(&s->cells_top[i], s->e->ti_current,
                           s->e->max_active_bin);
    }
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
  const struct space *s = (struct space *)extra_data;
  const int with_timestep_limiter =
      (s->e->policy & engine_policy_timestep_limiter);
  const int with_timestep_sync = (s->e->policy & engine_policy_timestep_sync);

  /* Verify that all limited particles have been treated */
  for (int k = 0; k < nr_parts; k++) {

    if (parts[k].time_bin == time_bin_inhibited) continue;

    if (parts[k].time_bin < 0) error("Particle has negative time-bin!");

    if (with_timestep_limiter &&
        parts[k].limiter_data.wakeup != time_bin_not_awake)
      error("Particle still woken up! id=%lld wakeup=%d", parts[k].id,
            parts[k].limiter_data.wakeup);

    if (with_timestep_sync && parts[k].limiter_data.to_be_synchronized != 0)
      error("Synchronized particle not treated! id=%lld synchronized=%d",
            parts[k].id, parts[k].limiter_data.to_be_synchronized);

    if (parts[k].gpart != NULL) {
      if (parts[k].time_bin != parts[k].gpart->time_bin) {
        error("Gpart not on the same time-bin as part %i %i", parts[k].time_bin,
              parts[k].gpart->time_bin);
      }
    }
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
                 s->nr_parts, sizeof(struct part), 1000, s);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief #threadpool mapper function for the swallow debugging check
 */
void space_check_part_swallow_mapper(void *map_data, int nr_parts,
                                     void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;

  /* Verify that all particles have been swallowed or are untouched */
  for (int k = 0; k < nr_parts; k++) {

    if (parts[k].time_bin == time_bin_inhibited) continue;

    const long long swallow_id =
        black_holes_get_part_swallow_id(&parts[k].black_holes_data);

    if (swallow_id != -1)
      error("Particle has not been swallowed! id=%lld", parts[k].id);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief #threadpool mapper function for the swallow debugging check
 */
void space_check_bpart_swallow_mapper(void *map_data, int nr_bparts,
                                      void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Unpack the data */
  struct bpart *restrict bparts = (struct bpart *)map_data;

  /* Verify that all particles have been swallowed or are untouched */
  for (int k = 0; k < nr_bparts; k++) {

    if (bparts[k].time_bin == time_bin_inhibited) continue;

    const long long swallow_id =
        black_holes_get_bpart_swallow_id(&bparts[k].merger_data);

    if (swallow_id != -1)
      error("BH particle has not been swallowed! id=%lld", bparts[k].id);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief #threadpool mapper function for the swallow debugging check
 */
void space_check_part_sink_swallow_mapper(void *map_data, int nr_parts,
                                          void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;

  /* Verify that all particles have been swallowed or are untouched */
  for (int k = 0; k < nr_parts; k++) {

    if (parts[k].time_bin == time_bin_inhibited) continue;

    const long long swallow_id = sink_get_part_swallow_id(&parts[k].sink_data);

    if (swallow_id != -1)
      error("Particle has not been swallowed! id=%lld", parts[k].id);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief #threadpool mapper function for the swallow debugging check
 */
void space_check_sink_sink_swallow_mapper(void *map_data, int nr_sinks,
                                          void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS
  /* Unpack the data */
  struct sink *restrict sinks = (struct sink *)map_data;

  /* Verify that all particles have been swallowed or are untouched */
  for (int k = 0; k < nr_sinks; k++) {

    if (sinks[k].time_bin == time_bin_inhibited) continue;

    const long long swallow_id =
        sink_get_sink_swallow_id(&sinks[k].merger_data);

    if (swallow_id != -1)
      error("Sink particle has not been swallowed! id=%lld", sinks[k].id);
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that all particles have their swallow flag in a "no swallow"
 * state.
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to check.
 */
void space_check_swallow(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS

  threadpool_map(&s->e->threadpool, space_check_part_swallow_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                 /*extra_data=*/NULL);

  threadpool_map(&s->e->threadpool, space_check_bpart_swallow_mapper, s->bparts,
                 s->nr_bparts, sizeof(struct bpart), threadpool_auto_chunk_size,
                 /*extra_data=*/NULL);

  threadpool_map(&s->e->threadpool, space_check_part_sink_swallow_mapper,
                 s->parts, s->nr_parts, sizeof(struct part),
                 threadpool_auto_chunk_size,
                 /*extra_data=*/NULL);

  threadpool_map(&s->e->threadpool, space_check_sink_sink_swallow_mapper,
                 s->sinks, s->nr_sinks, sizeof(struct sink),
                 threadpool_auto_chunk_size,
                 /*extra_data=*/NULL);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

void space_check_sort_flags_mapper(void *map_data, int nr_cells,
                                   void *extra_data) {

#ifdef SWIFT_DEBUG_CHECKS

  const struct space *s = (struct space *)extra_data;
  int *local_cells_top = map_data;

  for (int ind = 0; ind < nr_cells; ++ind) {
    const struct cell *c = &s->cells_top[local_cells_top[ind]];

    cell_check_sort_flags(c);
  }

#endif
}

/**
 * @brief Checks that all cells have cleared their sort flags.
 *
 * Should only be used for debugging purposes.
 *
 * @param s The #space to check.
 */
void space_check_sort_flags(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS

  threadpool_map(&s->e->threadpool, space_check_sort_flags_mapper,
                 s->local_cells_with_tasks_top, s->nr_local_cells_with_tasks,
                 sizeof(int), 1, s);
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
 * @brief Call the post-snapshot tracer on all the particles.
 *
 * @param s The #space.
 */
void space_after_snap_tracer(struct space *s, int verbose) {
  for (size_t i = 0; i < s->nr_parts; ++i) {
    tracers_after_snapshot_part(&s->parts[i], &s->xparts[i]);
  }
  for (size_t i = 0; i < s->nr_sparts; ++i) {
    tracers_after_snapshot_spart(&s->sparts[i]);
  }
  for (size_t i = 0; i < s->nr_bparts; ++i) {
    tracers_after_snapshot_bpart(&s->bparts[i]);
  }
}

/**
 * @brief Marks a top-level cell as having been updated by one of the
 * time-step updating tasks.
 */
void space_mark_cell_as_updated(struct space *s, const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c != c->top) error("Function can only be called at the top level!");
#endif

  /* Get the offest into the top-level cell array */
  const size_t delta = c - s->cells_top;

  /* Mark as updated */
  atomic_inc(&s->cells_top_updated[delta]);
}

/**
 * @brief Frees up the memory allocated for this #space
 */
void space_clean(struct space *s) {

  for (int i = 0; i < s->nr_cells; ++i) cell_clean(&s->cells_top[i]);
  swift_free("cells_top", s->cells_top);
  swift_free("cells_top_updated", s->cells_top_updated);
  swift_free("multipoles_top", s->multipoles_top);
  swift_free("local_cells_top", s->local_cells_top);
  swift_free("local_cells_with_tasks_top", s->local_cells_with_tasks_top);
  swift_free("cells_with_particles_top", s->cells_with_particles_top);
  swift_free("local_cells_with_particles_top",
             s->local_cells_with_particles_top);
  swift_free("parts", s->parts);
  swift_free("xparts", s->xparts);
  swift_free("gparts", s->gparts);
  swift_free("sparts", s->sparts);
  swift_free("bparts", s->bparts);
  swift_free("sinks", s->sinks);
#ifdef WITH_MPI
  swift_free("parts_foreign", s->parts_foreign);
  swift_free("sparts_foreign", s->sparts_foreign);
  swift_free("gparts_foreign", s->gparts_foreign);
  swift_free("bparts_foreign", s->bparts_foreign);
  swift_free("sinks_foreign", s->sinks_foreign);
#endif
  free(s->cells_sub);
  free(s->multipoles_sub);

  if (lock_destroy(&s->unique_id.lock) != 0)
    error("Failed to destroy spinlocks.");
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

  /* Now all our globals. */
  restart_write_blocks(&space_splitsize, sizeof(int), 1, stream,
                       "space_splitsize", "space_splitsize");
  restart_write_blocks(&space_maxsize, sizeof(int), 1, stream, "space_maxsize",
                       "space_maxsize");
  restart_write_blocks(&space_grid_split_threshold, sizeof(int), 1, stream,
                       "space_grid_split_threshold",
                       "space_grid_split_threshold");
  restart_write_blocks(&space_subsize_pair_hydro, sizeof(int), 1, stream,
                       "space_subsize_pair_hydro", "space_subsize_pair_hydro");
  restart_write_blocks(&space_subsize_self_hydro, sizeof(int), 1, stream,
                       "space_subsize_self_hydro", "space_subsize_self_hydro");
  restart_write_blocks(&space_subsize_pair_stars, sizeof(int), 1, stream,
                       "space_subsize_pair_stars", "space_subsize_pair_stars");
  restart_write_blocks(&space_subsize_self_stars, sizeof(int), 1, stream,
                       "space_subsize_self_stars", "space_subsize_self_stars");
  restart_write_blocks(&space_subsize_pair_grav, sizeof(int), 1, stream,
                       "space_subsize_pair_grav", "space_subsize_pair_grav");
  restart_write_blocks(&space_subsize_self_grav, sizeof(int), 1, stream,
                       "space_subsize_self_grav", "space_subsize_self_grav");
  restart_write_blocks(&space_subdepth_diff_grav, sizeof(int), 1, stream,
                       "space_subdepth_diff_grav", "space_subdepth_diff_grav");
  restart_write_blocks(&space_extra_parts, sizeof(int), 1, stream,
                       "space_extra_parts", "space_extra_parts");
  restart_write_blocks(&space_extra_gparts, sizeof(int), 1, stream,
                       "space_extra_gparts", "space_extra_gparts");
  restart_write_blocks(&space_extra_sinks, sizeof(int), 1, stream,
                       "space_extra_sinks", "space_extra_sinks");
  restart_write_blocks(&space_extra_sparts, sizeof(int), 1, stream,
                       "space_extra_sparts", "space_extra_sparts");
  restart_write_blocks(&space_extra_bparts, sizeof(int), 1, stream,
                       "space_extra_bparts", "space_extra_bparts");
  restart_write_blocks(&space_recurse_size_self_hydro, sizeof(int), 1, stream,
                       "space_recurse_size_self_hydro",
                       "space_recurse_size_self_hydro");
  restart_write_blocks(&space_recurse_size_pair_hydro, sizeof(int), 1, stream,
                       "space_recurse_size_pair_hydro",
                       "space_recurse_size_pair_hydro");
  restart_write_blocks(&space_recurse_size_self_stars, sizeof(int), 1, stream,
                       "space_recurse_size_self_stars",
                       "space_recurse_size_self_stars");
  restart_write_blocks(&space_recurse_size_pair_stars, sizeof(int), 1, stream,
                       "space_recurse_size_pair_stars",
                       "space_recurse_size_pair_stars");
  restart_write_blocks(&space_recurse_size_self_black_holes, sizeof(int), 1,
                       stream, "space_recurse_size_self_black_holes",
                       "space_recurse_size_self_black_holes");
  restart_write_blocks(&space_recurse_size_pair_black_holes, sizeof(int), 1,
                       stream, "space_recurse_size_pair_black_holes",
                       "space_recurse_size_pair_black_holes");
  restart_write_blocks(&space_recurse_size_self_sinks, sizeof(int), 1, stream,
                       "space_recurse_size_self_sinks",
                       "space_recurse_size_self_sinks");
  restart_write_blocks(&space_recurse_size_pair_sinks, sizeof(int), 1, stream,
                       "space_recurse_size_pair_sinks",
                       "space_recurse_size_pair_sinks");
  restart_write_blocks(&space_expected_max_nr_strays, sizeof(int), 1, stream,
                       "space_expected_max_nr_strays",
                       "space_expected_max_nr_strays");
  restart_write_blocks(&engine_max_parts_per_ghost, sizeof(int), 1, stream,
                       "engine_max_parts_per_ghost",
                       "engine_max_parts_per_ghost");
  restart_write_blocks(&engine_max_sparts_per_ghost, sizeof(int), 1, stream,
                       "engine_max_sparts_per_ghost",
                       "engine_max_sparts_per_ghost");
  restart_write_blocks(&engine_max_parts_per_cooling, sizeof(int), 1, stream,
                       "engine_max_parts_per_cooling",
                       "engine_max_parts_per_cooling");
  restart_write_blocks(&engine_star_resort_task_depth, sizeof(int), 1, stream,
                       "engine_star_resort_task_depth",
                       "engine_star_resort_task_depth");
  restart_write_blocks(&engine_redistribute_alloc_margin, sizeof(double), 1,
                       stream, "engine_redistribute_alloc_margin",
                       "engine_redistribute_alloc_margin");
  restart_write_blocks(&engine_foreign_alloc_margin, sizeof(double), 1, stream,
                       "engine_foreign_alloc_margin",
                       "engine_foreign_alloc_margin");

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

  if (s->nr_sinks > 0)
    restart_write_blocks(s->sinks, s->nr_sinks, sizeof(struct sink), stream,
                         "sinks", "sinks");

  if (s->nr_sparts > 0)
    restart_write_blocks(s->sparts, s->nr_sparts, sizeof(struct spart), stream,
                         "sparts", "sparts");
  if (s->nr_bparts > 0)
    restart_write_blocks(s->bparts, s->nr_bparts, sizeof(struct bpart), stream,
                         "bparts", "bparts");
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

  /* Now all our globals. */
  restart_read_blocks(&space_splitsize, sizeof(int), 1, stream, NULL,
                      "space_splitsize");
  restart_read_blocks(&space_maxsize, sizeof(int), 1, stream, NULL,
                      "space_maxsize");
  restart_read_blocks(&space_grid_split_threshold, sizeof(int), 1, stream, NULL,
                      "space_grid_split_threshold");
  restart_read_blocks(&space_subsize_pair_hydro, sizeof(int), 1, stream, NULL,
                      "space_subsize_pair_hydro");
  restart_read_blocks(&space_subsize_self_hydro, sizeof(int), 1, stream, NULL,
                      "space_subsize_self_hydro");
  restart_read_blocks(&space_subsize_pair_stars, sizeof(int), 1, stream, NULL,
                      "space_subsize_pair_stars");
  restart_read_blocks(&space_subsize_self_stars, sizeof(int), 1, stream, NULL,
                      "space_subsize_self_stars");
  restart_read_blocks(&space_subsize_pair_grav, sizeof(int), 1, stream, NULL,
                      "space_subsize_pair_grav");
  restart_read_blocks(&space_subsize_self_grav, sizeof(int), 1, stream, NULL,
                      "space_subsize_self_grav");
  restart_read_blocks(&space_subdepth_diff_grav, sizeof(int), 1, stream, NULL,
                      "space_subdepth_diff_grav");
  restart_read_blocks(&space_extra_parts, sizeof(int), 1, stream, NULL,
                      "space_extra_parts");
  restart_read_blocks(&space_extra_gparts, sizeof(int), 1, stream, NULL,
                      "space_extra_gparts");
  restart_read_blocks(&space_extra_sinks, sizeof(int), 1, stream, NULL,
                      "space_extra_sinks");
  restart_read_blocks(&space_extra_sparts, sizeof(int), 1, stream, NULL,
                      "space_extra_sparts");
  restart_read_blocks(&space_extra_bparts, sizeof(int), 1, stream, NULL,
                      "space_extra_bparts");
  restart_read_blocks(&space_recurse_size_self_hydro, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_self_hydro");
  restart_read_blocks(&space_recurse_size_pair_hydro, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_pair_hydro");
  restart_read_blocks(&space_recurse_size_self_stars, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_self_stars");
  restart_read_blocks(&space_recurse_size_pair_stars, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_pair_stars");
  restart_read_blocks(&space_recurse_size_self_black_holes, sizeof(int), 1,
                      stream, NULL, "space_recurse_size_self_black_holes");
  restart_read_blocks(&space_recurse_size_pair_black_holes, sizeof(int), 1,
                      stream, NULL, "space_recurse_size_pair_black_holes");
  restart_read_blocks(&space_recurse_size_self_sinks, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_self_sinks");
  restart_read_blocks(&space_recurse_size_pair_sinks, sizeof(int), 1, stream,
                      NULL, "space_recurse_size_pair_sinks");
  restart_read_blocks(&space_expected_max_nr_strays, sizeof(int), 1, stream,
                      NULL, "space_expected_max_nr_strays");
  restart_read_blocks(&engine_max_parts_per_ghost, sizeof(int), 1, stream, NULL,
                      "engine_max_parts_per_ghost");
  restart_read_blocks(&engine_max_sparts_per_ghost, sizeof(int), 1, stream,
                      NULL, "engine_max_sparts_per_ghost");
  restart_read_blocks(&engine_max_parts_per_cooling, sizeof(int), 1, stream,
                      NULL, "engine_max_parts_per_cooling");
  restart_read_blocks(&engine_star_resort_task_depth, sizeof(int), 1, stream,
                      NULL, "engine_star_resort_task_depth");
  restart_read_blocks(&engine_redistribute_alloc_margin, sizeof(double), 1,
                      stream, NULL, "engine_redistribute_alloc_margin");
  restart_read_blocks(&engine_foreign_alloc_margin, sizeof(double), 1, stream,
                      NULL, "engine_foreign_alloc_margin");

  /* Things that should be reconstructed in a rebuild. */
  s->cells_top = NULL;
  s->cells_top_updated = NULL;
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
  s->bparts_foreign = NULL;
  s->size_bparts_foreign = 0;
  s->sinks_foreign = NULL;
  s->size_sinks_foreign = 0;
#endif

  /* More things to read. */
  s->parts = NULL;
  s->xparts = NULL;
  if (s->nr_parts > 0) {

    /* Need the memory for these. */
    if (swift_memalign("parts", (void **)&s->parts, part_align,
                       s->size_parts * sizeof(struct part)) != 0)
      error("Failed to allocate restore part array.");

    if (swift_memalign("xparts", (void **)&s->xparts, xpart_align,
                       s->size_parts * sizeof(struct xpart)) != 0)
      error("Failed to allocate restore xpart array.");

    restart_read_blocks(s->parts, s->nr_parts, sizeof(struct part), stream,
                        NULL, "parts");
    restart_read_blocks(s->xparts, s->nr_parts, sizeof(struct xpart), stream,
                        NULL, "xparts");
  }
  s->gparts = NULL;
  if (s->nr_gparts > 0) {
    if (swift_memalign("gparts", (void **)&s->gparts, gpart_align,
                       s->size_gparts * sizeof(struct gpart)) != 0)
      error("Failed to allocate restore gpart array.");

    restart_read_blocks(s->gparts, s->nr_gparts, sizeof(struct gpart), stream,
                        NULL, "gparts");
  }

  s->sinks = NULL;
  if (s->nr_sinks > 0) {
    if (swift_memalign("sinks", (void **)&s->sinks, sink_align,
                       s->size_sinks * sizeof(struct sink)) != 0)
      error("Failed to allocate restore sink array.");

    restart_read_blocks(s->sinks, s->nr_sinks, sizeof(struct sink), stream,
                        NULL, "sinks");
  }

  s->sparts = NULL;
  if (s->nr_sparts > 0) {
    if (swift_memalign("sparts", (void **)&s->sparts, spart_align,
                       s->size_sparts * sizeof(struct spart)) != 0)
      error("Failed to allocate restore spart array.");

    restart_read_blocks(s->sparts, s->nr_sparts, sizeof(struct spart), stream,
                        NULL, "sparts");
  }
  s->bparts = NULL;
  if (s->nr_bparts > 0) {
    if (swift_memalign("bparts", (void **)&s->bparts, bpart_align,
                       s->size_bparts * sizeof(struct bpart)) != 0)
      error("Failed to allocate restore bpart array.");

    restart_read_blocks(s->bparts, s->nr_bparts, sizeof(struct bpart), stream,
                        NULL, "bparts");
  }

  /* Need to reconnect the gravity parts to their hydro, star and BH particles.
   * Note that we can't use the threadpool here as we have not restored it yet.
   */

  /* Re-link the parts. */
  if (s->nr_parts > 0 && s->nr_gparts > 0)
    part_relink_parts_to_gparts(s->gparts, s->nr_gparts, s->parts);

  /* Re-link the sinks. */
  if (s->nr_sinks > 0 && s->nr_gparts > 0)
    part_relink_sinks_to_gparts(s->gparts, s->nr_gparts, s->sinks);

  /* Re-link the sparts. */
  if (s->nr_sparts > 0 && s->nr_gparts > 0)
    part_relink_sparts_to_gparts(s->gparts, s->nr_gparts, s->sparts);

  /* Re-link the bparts. */
  if (s->nr_bparts > 0 && s->nr_gparts > 0)
    part_relink_bparts_to_gparts(s->gparts, s->nr_gparts, s->bparts);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that everything is correct */
  part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                    s->nr_parts, s->nr_gparts, s->nr_sinks, s->nr_sparts,
                    s->nr_bparts, 1);
#endif
}

#define root_cell_id 0
/**
 * @brief write a single cell in a csv file.
 *
 * @param s The #space.
 * @param f The file to use (already open).
 * @param c The current #cell.
 */
void space_write_cell(const struct space *s, FILE *f, const struct cell *c) {
#ifdef SWIFT_CELL_GRAPH

  if (c == NULL) return;

  /* Get parent ID */
  long long parent = root_cell_id;
  if (c->parent != NULL) parent = c->parent->cellID;

  /* Get super ID */
  char superID[100] = "";
  if (c->super != NULL) sprintf(superID, "%lld", c->super->cellID);

  /* Get hydro super ID */
  char hydro_superID[100] = "";
  if (c->hydro.super != NULL)
    sprintf(hydro_superID, "%lld", c->hydro.super->cellID);

  /* Write line for current cell */
  fprintf(f, "%lld,%lld,%i,", c->cellID, parent, c->nodeID);
  fprintf(f, "%i,%i,%i,%i,%s,%s,%g,%g,%g,%g,%g,%g, ", c->hydro.count,
          c->stars.count, c->grav.count, c->sinks.count, superID, hydro_superID,
          c->loc[0], c->loc[1], c->loc[2], c->width[0], c->width[1],
          c->width[2]);
  fprintf(f, "%g, %g, %i, %i\n", c->hydro.h_max, c->stars.h_max, c->depth,
          c->maxdepth);

  /* Write children */
  for (int i = 0; i < 8; i++) {
    space_write_cell(s, f, c->progeny[i]);
  }
#endif
}

/**
 * @brief Write a csv file containing the cell hierarchy
 *
 * @param s The #space.
 * @param j The file number.
 */
void space_write_cell_hierarchy(const struct space *s, int j) {

#ifdef SWIFT_CELL_GRAPH

  /* Open file */
  char filename[200];
  sprintf(filename, "cell_hierarchy_%04i_%04i.csv", j, engine_rank);
  FILE *f = fopen(filename, "w");
  if (f == NULL) error("Error opening task level file.");

  const int root_id = root_cell_id;
  /* Write header */
  if (engine_rank == 0) {
    fprintf(f, "name,parent,mpi_rank,");
    fprintf(f,
            "hydro_count,stars_count,gpart_count,sinks_count,super,hydro_super,"
            "loc1,loc2,loc3,width1,width2,width3,");
    fprintf(f, "hydro_h_max,stars_h_max,depth,maxdepth\n");

    /* Write root data */
    fprintf(f, "%i, ,-1,", root_id);
    fprintf(f, "%li,%li,%li, , , , , , , , ,", s->nr_parts, s->nr_sparts,
            s->nr_gparts);
    fprintf(f, " , , ,\n");
  }

  /* Write all the top level cells (and their children) */
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->nodeID == engine_rank) space_write_cell(s, f, c);
  }

  /* Cleanup */
  fclose(f);
#endif
}

/**
 * @brief Check the unskip flags for the cell and its progenies.
 *
 * @param c The current #cell.
 */
void space_recurse_check_unskip_flag(const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS

  /* Check the current cell. */
  if (cell_get_flag(c, cell_flag_unskip_self_grav_processed)) {
    error(
        "A cell is still containing a self unskip flag for the gravity "
        "depth=%d node=%d cellID=%lld",
        c->depth, c->nodeID, c->cellID);
  }
  if (cell_get_flag(c, cell_flag_unskip_pair_grav_processed)) {
    error(
        "A cell is still containing a pair unskip flag for the gravity "
        "depth=%d node=%d cellID=%lld",
        c->depth, c->nodeID, c->cellID);
  }

  /* Recurse */
  for (int i = 0; i < 8; i++) {
    if (c->progeny[i] != NULL) space_recurse_check_unskip_flag(c->progeny[i]);
  }
#endif
}

/**
 * @brief Loop over all the cells and ensure that the unskip
 * flag have been cleared.
 *
 * @param s The #space
 */
void space_check_unskip_flags(const struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS

  for (int i = 0; i < s->nr_cells; i++) {
    const struct cell *c = &s->cells_top[i];
    space_recurse_check_unskip_flag(c);
  }
#else
  error("This function should not be called without the debugging checks.");
#endif
}

struct cic_mapper_data {
  const struct cell* cells;
  double* rho;
  double* potential;
  int N;
  int use_local_patches;
  double fac;
  double dim[3];
  float const_G;
  struct neutrino_model* nu_model;
};

/**
 * @brief Compute the forces and potential on a series of uniform grids.
 *
 * The FMG algorithm is used: on every level the density is calculated, 
 * and at every level the multigrid method (i.e. with V-cycles) is applied
 * to get the potential. The potential guess on the next finer grid is set 
 * by prolongating the potential from a grid coarser.
 *
 * @param e The #engine.
 * @param N_min 1D size of the smallest grid on which to solve for the potential.
 * @param N_max 1D size of the largest grid on which to solve for the potential.
 */
void space_apply_FMG(const struct engine *e, const int N_min, const int N_max) {
  const struct space *s = e->s;
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  if (dim[0] != dim[1] || dim[0] != dim[2]) {
    message("Noncubic domain for density calculation. Did not fill density array");
    return;
  }
  if (N_max%2 != 0){
    error("Uneven domain for density calculation. Did not apply FMG");
  }

  /* Get array of the different grid sizes */
  int n = N_min;
  int level = 0;
  int grid_sizes[32];
  while (n<=N_max) {
    grid_sizes[level++] = n;
    n *= 2;
  }

  const int N_levels = level;
  double *rho[N_levels];
  double *pot[N_levels];
  double mean_density[N_levels];

  struct cic_mapper_data data;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];

  /* Initialise arrays for density and potential calculation. */
  for (int i=0; i<N_levels; i++) {
    rho[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    pot[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    
    /* Assign densities to all arrays */ 
    data.N = grid_sizes[i];
    data.rho = rho[i];
    data.fac = grid_sizes[i]/box_size;
    gpart_to_mesh_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);
    mean_density[i] = get_mean_density(data.rho, grid_sizes[i], 0);
    mean_density[i] *= 4.*M_PI * data.fac*data.fac*data.fac;
  }
  
  /* Set initial guess on the coarsest grid and solve directly */
  int cdim[3] = {N_min, N_min, N_min};
  set_initial_guess(pot[0], cdim, /*MG=*/0);
  apply_GS(rho[0], pot[0], cdim, mean_density[0], box_size);

  /* Solve on all finer grids by prolongating from the previous grid and performing the multigrid method */
  for (int i = 1; i<N_levels; i++) {
    message("Going to the level with grid size %d \n", grid_sizes[i]);
    prolongate_solution(pot[i-1], pot[i], grid_sizes[i-1], grid_sizes[i]); 
    
    int N = grid_sizes[i];
    cdim[0] = N;
    cdim[1] = N;
    cdim[2] = N;
    apply_multigrid(rho[i], pot[i], cdim, mean_density[i], box_size, N_min, N, 15);
  }

  int get_cell_acc = 0;
  if (get_cell_acc) {
    double const_G = s->e->physical_constants->const_newton_G;
    const double dx = s->dim[0]/N_max;

    int cdim_max[3] = {N_max, N_max, N_max};

    double *acc[3];
    for (int d = 0; d < 3; d++) {
      acc[d] = malloc((size_t)N_max * N_max * N_max * sizeof(double));
    }

    for (int i = 0; i < N_max; i++) {
      int ip = (i + 1) % N_max;
      int im = (i - 1 + N_max) % N_max;
      for (int j = 0; j < N_max; j++) {
        int jp = (j + 1) % N_max;
        int jm = (j - 1 + N_max) % N_max;
        for (int k = 0; k < N_max; k++) {

          int kp = (k + 1) % N_max;
          int km = (k - 1 + N_max) % N_max;

          int cid    = cell_getid(cdim_max, i,  j,  k);

          int cid_xp = cell_getid(cdim_max, ip, j,  k);
          int cid_xm = cell_getid(cdim_max, im, j,  k);

          int cid_yp = cell_getid(cdim_max, i,  jp, k);
          int cid_ym = cell_getid(cdim_max, i,  jm, k);

          int cid_zp = cell_getid(cdim_max, i,  j,  kp);
          int cid_zm = cell_getid(cdim_max, i,  j,  km);

          acc[0][cid] =
              -const_G *
              (pot[N_levels-1][cid_xp] - pot[N_levels-1][cid_xm])
              / (2.0 * dx);

          acc[1][cid] =
              -const_G *
              (pot[N_levels-1][cid_yp] - pot[N_levels-1][cid_ym])
              / (2.0 * dx);

          acc[2][cid] =
              -const_G *
              (pot[N_levels-1][cid_zp] - pot[N_levels-1][cid_zm])
              / (2.0 * dx);
        }
      }
    }
    for (int d = 0; d < 3; d++) {
      free(acc[d]);
    }
  }

  /* Now calculate the accelerations from the potential */
  data.rho = NULL;
  data.potential = pot[N_levels -1];
  data.fac = N_max / box_size;
  data.dim[0] = s->dim[0];
  data.dim[1] = s->dim[1];
  data.dim[2] = s->dim[2];
  data.const_G = s->e->physical_constants->const_newton_G;
  mesh_to_gpart_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);

  /* Free up memory */
  for (int i=0; i<N_levels; i++) {
    free(rho[i]);
    free(pot[i]);
  }
}

/**
 * @brief Set an initial guess on a grid by prolongation from a grid coarser.
 *
 * The chosen prolongation operator is trilinear interpolation. 
 *
 * @param pot_coarse Array containing the potential on the previous grid.
 * @param pot_fine Array to set the initial guess for the fine grid in. 
 * @param N 1D size of the coarse grid.
 * @param N_double 1D size of the fine grid.
 */
void prolongate_solution(const double *pot_coarse, double *pot_fine, const int N, const int N_double) {
  int cdimh[3] = {N_double, N_double, N_double};
  int cdim[3] = {N,N,N};

  for (int i = 0; i<N_double; i++) {
    for (int j = 0; j<N_double; j++) {
      for (int k = 0; k<N_double; k++) {
        int iH = i/2;
        int jH = j/2;
        int kH = k/2;
        double dx = (i%2) * 0.5;
        double dy = (j%2) * 0.5;
        double dz = (k%2) * 0.5;

        int iHp = (iH +1)%N;
        int jHp = (jH +1)%N;
        int kHp = (kH +1)%N;

        pot_fine[cell_getid(cdimh,i,j,k)] = ((1.-dx)*(1.-dy)*(1.-dz)*pot_coarse[cell_getid(cdim,iH,jH,kH)] + dx*(1.-dy)*(1.-dz)*pot_coarse[cell_getid(cdim,iHp, jH, kH)]
                                  + (1.-dx)*dy*(1.-dz)*pot_coarse[cell_getid(cdim,iH,jHp,kH)] + (1.-dx)*(1.-dy)*dz*pot_coarse[cell_getid(cdim,iH,jH,kHp)]
                                  + dx*dy*(1.-dz)*pot_coarse[cell_getid(cdim,iHp,jHp,kH)] + dx*(1.-dy)*dz*pot_coarse[cell_getid(cdim,iHp,jH, kHp)]
                                  + (1.-dx)*dy*dz*pot_coarse[cell_getid(cdim,iH,jHp,kHp)]+ dx*dy*dz*pot_coarse[cell_getid(cdim,iHp,jHp,kHp)]);
      }
    }
  }
}

/**
 * @brief Compute the forces and potential on a uniform grid.
 *
 * Interpolates the gparts onto a mesh, performs Gauss-Seidel relaxation
 * to compute the potential on the grid, and interpolates the potential 
 * and acceleration back to the particles. We use CIC for the interpolation.
 *
 * This function calls the appropriate relaxation method (multigrid acceleration
 * or directly solving) based on the parameter multigrid.
 *
 * @param e The #engine.
 * @param N 1D size of the grid on which to solve for the potential.
 * @param multigrid Are we doing multigrid acceleration?
 */
void space_get_density(struct engine *e, const int N, int multigrid) {
  int cdim[3] = {N,N,N};
  const struct space *s = e->s;
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double fac = ((double) N)/box_size;

  /* Check if the space and mesh satisfy the conditions for Gauss-Seidel */
  if (dim[0] != dim[1] || dim[0] != dim[2]) {
    error("Noncubic domain for density calculation. Did not fill density array");
  }
  if (N%2 != 0){
    error("Uneven domain for density calculation. Did not apply Gauss-Seidel");
  }

  /* Allocate memory for density and potential arrays (on the finest level) */
  double *density_array = NULL;
  density_array = (double*)calloc(N * N * N, sizeof(double));
  if (density_array == NULL)
    error("Error allocating memory for the density mesh.");
  memuse_log_allocation("mesh.density", density_array, 1,
                        sizeof(double) * N * N * N);
  double *potential_array = NULL;
  potential_array = (double*)calloc(N * N * N, sizeof(double));
  if (potential_array == NULL){
    error("Error allocating memory for the potential mesh.");
  }
  memuse_log_allocation("mesh.potential", potential_array, 1, sizeof(double)*N*N*N);
       
  /* Prepare for and initialise CIC to get the density at the mesh */
  struct cic_mapper_data data; 
  message("The number of cells is %d", N);
  data.N = N;
  data.rho = density_array;
  data.fac = fac;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];
  gpart_to_mesh_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);
 
  /* Have to do something about this */
  //get_pm_potential(&data, N, box_size, &e->threadpool, cdim);  //PM method to get the potential is implemented here
 
  /* Rescale density to get overdensity and rescale mean density for consistency with FFT */
  double mean_density = get_mean_density(density_array, N, 0);
  mean_density *= 4.* M_PI * fac*fac*fac;
  double sum = 0.0;
  double rms = 0.0;
  for (int i = 0; i < N*N*N; i++) {
    rms += fabs(density_array[i]-1);
    sum += density_array[i]-1;
  }
  if (sum < -0.1 || sum > 0.1) message("Warning! Mean of the overdensity is nonzero");

  set_initial_guess(potential_array, cdim, /*MG=*/0);

  if (multigrid) {
    const int N_min = 16; //Size of the smallest grid to be used in multigrid acceleration 
    apply_multigrid(density_array, potential_array, cdim, mean_density, box_size, N_min, N, 15);
  } else {
    apply_GS(density_array, potential_array, cdim, mean_density, box_size);
  }

  free(density_array);
  free(potential_array);
}

/**
 * @brief Compute potential on the grid using multigrid acceleration.
 *
 * Solves for the potential on the finest grid, passes the residual to 
 * the coarser grids and recursively solves for it on every level
 * using V-cycles.
 *
 * @param rho Array containing the density on the finest grid.
 * @param pot Array containing (approximate) potential values on the finest grid.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param box_size Side length of the box.
 * @param N_min 1D size of the coarsest grid.
 * @param N_max 1D size of the finest grid.
 * @param V_max Maximum number of V-cycles that may be performed.
 */
void apply_multigrid(const double *rho, double *pot, int cdim[3], const double mean_density, const double box_size, const int N_min, const int N_max, const int V_max) {
  message("Applying the multigrid method...");
  //int N_print = 256;

  /* Allocate the memory for the residual array on the finest level */
  double *residual_array = NULL;
  residual_array = (double*)malloc(sizeof(double) * N_max * N_max * N_max);
  if (residual_array == NULL){
    error("Error allocating memory for the residual array.");
  }
  memuse_log_allocation("residual.array", residual_array, 1, sizeof(double)*N_max*N_max*N_max);

  double delta = box_size/N_max; //Width of a grid cell of the finest level
  double residual; 
  residual = get_residual(pot, rho, cdim, mean_density, delta);
  message("The first residual is %lf", residual);
  const double tolerance = 10e-10; //Choose reasonable value here
  int counter = 0;
  const int fine_steps = 10; //Choose reasonable value here
  double sum = 0;

  int V_cycles=0;

  while (residual > tolerance) {
    message("Performing V-cycle %d", V_cycles);
    /* Pre-smoothing */
    for (int i=0; i<fine_steps; i++) {
      perform_red_black_sweep(pot, rho, cdim, mean_density, delta); 
    }
    residual = get_residual(pot, rho, cdim, mean_density, delta);
    get_residual_array(pot, rho, cdim, mean_density, residual_array, delta);

    /* Transfer residual array to get coarse-grid correction */
    solve_coarser_problem_recursive(pot, residual_array, cdim, delta, N_min, N_max);

    /* Post-smoothing if needed */
    residual = get_residual(pot, rho, cdim, mean_density, delta);
    message("Back on the finest grid the residual is %E", residual);
    if (residual > tolerance) {
      for (int i=0; i<fine_steps; i++) {
        perform_red_black_sweep(pot, rho, cdim, mean_density, delta); 
      }
    }
    residual = get_residual(pot, rho, cdim, mean_density, delta);  
    V_cycles +=1;
    message("After %d V-cycle(s) the residual is %E", V_cycles, residual);
  }

  message("Performed %d V-cycle(s) in total", V_cycles);
  
  /* Post-smoothing until convergence. Should not be necessary! */
  while (residual > tolerance) {
    perform_red_black_sweep(pot, rho, cdim, mean_density, delta);
    residual = get_residual(pot, rho, cdim, mean_density, delta);
    counter +=1;
    if (counter % 100 ==0){
      sum = 0;
      for (int i = 0; i < N_max*N_max*N_max; i++) {
          sum += pot[i];
      }
    }
  }
  residual = get_residual(pot, rho, cdim, mean_density, delta);
  message("Needed to do %d step(s) in post-smoothing after which the residual was %E", counter, residual);
  free(residual_array);

  /* Set mean of potential to zero */
  sum = 0.;
  for (int i = 0; i < N_max*N_max*N_max; i++) {
    sum += pot[i];
  }
  for (int i = 0; i<N_max*N_max*N_max; i++) {
    pot[i] -= sum/(N_max*N_max*N_max); 
  }

  /* Verify that the potential mean is now zero */
  double potential_sum = 0.;
  for (int i =0; i<N_max*N_max*N_max; i++) {
    potential_sum += pot[i];
  }
  if (potential_sum < -0.1 || potential_sum > 0.1) error("The mean of the potential is %lf", potential_sum/(N_max*N_max*N_max));
}

/**
 * @brief Directly solve for the potential with Gauss-Seidel relaxation.
 * 
 * Only one grid is used and no V-cycles are done.
 * 
 * @param rho Array containing the density on the grid.
 * @param pot Array containing (approximate) potential values on the grid.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param box_size Side length of the box.
 */
void apply_GS(const double *rho, double *pot, int cdim[3], double mean_density, double box_size) {
  message("Applying Gauss-Seidel...");
  int N = cdim[0];
  double delta = box_size/N;
  double residual = get_residual(pot, rho, cdim, mean_density, delta);
  double tolerance = 10e-12; //Choose reasonable value here
  int counter = 0;
  double sum = 0.;

  while (residual >= tolerance) {
    perform_red_black_sweep(pot, rho, cdim, mean_density, delta);
    residual = get_residual(pot, rho, cdim, mean_density, delta);
    counter +=1;
    if (counter %50 == 0) message("Did %d steps and the residual is %E", counter, residual);
  }
  message("We needed %d steps to converge", counter);

  /* Set mean of potential to zero */
  sum = 0.;
  for (int i = 0; i < N*N*N; i++) {
      sum += pot[i];
  }
  for (int i = 0; i<N*N*N; i++) {
    pot[i] -= sum/(N*N*N);
  }

  double potential_sum = 0;
  for (int i =0; i<N*N*N; i++)
    potential_sum += fabs(pot[i]);  
}

/**
 * @brief Executes the V-cycles. 
 *
 * At the coarsest level, the residual equation is solved exactly. On all intermediate
 * levels pre- and post-smoothing is done to solve the residual equation before the 
 * residual of this problem is restricted to the next coarser grid. 
 *
 * @param pot Array containing the approximate solution on one grid finer.
 * @param residual Array containing the density on the grid.
 * @param cdim 3D size of one grid finer.
 * @param delta Cell width of one grid finer.
 * @param N_stop Size of coarsest grid in the V-cycle.
 * @param N_start Size of finest grid in the V-cycle.
 */
void solve_coarser_problem_recursive(double *pot, const double *residual, int cdim[3], double delta, const int N_stop, const int N_start) {
  int N = cdim[0]; //Grid size of the current level we are on
  delta = delta*2.0; //Cells are twice as big on the coarser grid
  N = N/2; 

  /* Below arrays are the equivalent of the density, potential and residual arrays on the finest grid */
  /* Allocate the memory for the arrays of the problem on this level */
  double *coarser_residual = NULL;
  coarser_residual = (double*)malloc(sizeof(double) * N * N * N);
  if (coarser_residual == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.grid", coarser_residual, 1, sizeof(double)*N*N*N);

  double *coarser_solution = NULL;
  coarser_solution = (double*)calloc(N * N * N, sizeof(double)); //Set initial guess to zero as we are essentially solving for an error.
  if (coarser_solution == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.solution", coarser_solution, 1, sizeof(double)*N*N*N);

  double *new_residual_array = NULL; 
  new_residual_array = (double*)malloc(sizeof(double) * N * N * N);
  if (new_residual_array == NULL){
    error("Error allocating memory for the new coarser residual.");
  }
  memuse_log_allocation("coarser.newresidual", new_residual_array, 1, sizeof(double)*N*N*N);

  restrict_residual(coarser_residual, residual, cdim); //Transfer residual to the coarser grid
  int cdimH[] = {N, N, N}; 
  double multiplier = 1.0; 
  double tolerance = 10e-8; //Choose a reasonable value here
  int counter = 0;

  /* Multiply the residual array on H by -1 in order to solve the right equation */
  for (int i =0; i<N*N*N; i++) 
    coarser_residual[i] = -1.0*coarser_residual[i];
  double coarser_residual_abs = get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
  message("The first residual on the grid with size %d is %lf", N, coarser_residual_abs);

  /* Solve the equation exactly if we are on the coarsest grid */
  if (N==N_stop) {
    while (coarser_residual_abs >= tolerance) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
      coarser_residual_abs = get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
      counter +=1;
    }
    coarser_residual_abs = get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
    message("The total number of steps on the coarsest grid is %d and the residual is %lf", counter, coarser_residual_abs);
    prolongate_residual(coarser_solution, pot, cdimH);
  }

  /* Do some smoothing and proceed to coarser grids */
  else { 
    counter = 0;
    int coarse_steps = 50;
    for (int i=0; i<coarse_steps; i++) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
      coarser_residual_abs =  get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
      counter +=1;
      //if (coarser_residual_abs < tolerance) break;
    }
    get_residual_array(coarser_solution, coarser_residual, cdimH, multiplier, new_residual_array, delta);

    coarser_residual_abs = get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
    message("The residual after pre-smoothing is %lf", coarser_residual_abs);
    solve_coarser_problem_recursive(coarser_solution, new_residual_array, cdimH, delta, N_stop, N_start);
    /* Post-smoothing */
    coarser_residual_abs = get_residual(coarser_solution, coarser_residual, cdimH, multiplier, delta);
    message("Back on the grid with size %d the residual is %lf", N, coarser_residual_abs);
    for (int i=0; i<coarse_steps; i++) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
    }
    prolongate_residual(coarser_solution, pot, cdimH);
  }

  /* The coarser-grid correction has now been added to the finer-grid solution for the potential, so discard used arrays. */
  free(coarser_residual);
  free(coarser_solution);
  free(new_residual_array);
}

/**
 * @brief Add the converged solution of the coarser grid to the solution on the finer grid.
 * 
 * The chosen prolongation operator is trilinear interpolation. 
 *
 * @param coarser_solution Solution on the coarser grid.
 * @param pot Approximate solution on the finer grid. 
 * @param cdim 3D size of the coarser grid.
 */
void prolongate_residual(const double *coarser_solution, double *pot, int cdim[3]) {
  int N = cdim[0];
  int N_double = cdim[0]*2;
  int cdimh[] = {N_double, N_double, N_double};

  for (int i = 0; i<N_double; i++) {
    for (int j = 0; j<N_double; j++) {
      for (int k = 0; k<N_double; k++) {
        int iH = i/2;
        int jH = j/2;
        int kH = k/2;
        double dx = (i%2) * 0.5;
        double dy = (j%2) * 0.5;
        double dz = (k%2) * 0.5;

        int iHp = (iH +1)%N;
        int jHp = (jH +1)%N;
        int kHp = (kH +1)%N;

        pot[cell_getid(cdimh,i,j,k)] += ((1.-dx)*(1.-dy)*(1.-dz)*coarser_solution[cell_getid(cdim,iH,jH,kH)] + dx*(1.-dy)*(1.-dz)*coarser_solution[cell_getid(cdim,iHp, jH, kH)]
                                  + (1.-dx)*dy*(1.-dz)*coarser_solution[cell_getid(cdim,iH,jHp,kH)] + (1.-dx)*(1.-dy)*dz*coarser_solution[cell_getid(cdim,iH,jH,kHp)]
                                  + dx*dy*(1.-dz)*coarser_solution[cell_getid(cdim,iHp,jHp,kH)] + dx*(1.-dy)*dz*coarser_solution[cell_getid(cdim,iHp,jH, kHp)]
                                  + (1.-dx)*dy*dz*coarser_solution[cell_getid(cdim,iH,jHp,kHp)]+ dx*dy*dz*coarser_solution[cell_getid(cdim,iHp,jHp,kHp)]);
      }
    }
  }

  /* Prepare for the next layer */
  cdim[0] *= 2.;
  cdim[1] *= 2.;
  cdim[2] *= 2.;
}

/**
 * @brief Restrict the residual to a coarser grid by taking the average of the children.
 *
 * @param H_array Array to restrict the residual to.
 * @param residual Residual being passed to a coarser grid.
 * @param cdim 3D size of the fine grid.
 */
void restrict_residual(double *H_array, const double *residual, int cdim[3]) {
  int cdimH[] = {cdim[0]/2, cdim[1]/2, cdim[2]/2};
  for (int i=0; i<cdimH[0]; i++) {
    for (int j=0; j<cdimH[0]; j++) {
      for (int k=0; k<cdimH[0]; k++) {
        H_array[cell_getid(cdimH,i,j,k)] = (residual[cell_getid(cdim,2*i,2*j,2*k)] + residual[cell_getid(cdim,2*i+1,2*j,2*k)]
                                            +residual[cell_getid(cdim,2*i,2*j+1,2*k)] + residual[cell_getid(cdim,2*i,2*j,2*k+1)]
                                            +residual[cell_getid(cdim,2*i+1,2*j+1,2*k)] + residual[cell_getid(cdim,2*i+1,2*j,2*k+1)]
                                            +residual[cell_getid(cdim,2*i,2*j+1,2*k+1)] + residual[cell_getid(cdim,2*i+1,2*j+1,2*k+1)])/8.0;
      }
    }
  }
}

/**
 * @brief Get the residual at every grid point of the current approximation to the potential.
 *
 * Calculates the difference between the Laplacian of the potential and the right 
 * hand side of the Poisson equation.
 *
 * @param pot Array containing the potential on the grid.
 * @param rho Array containing the density on the grid.
 * @param cdim 3D size of the grid.
 * @param multiplier Multiplicative term on the right hand side of the Poisson equation.
 * @param residual Array containing the residual at every grid point.
 * @param delta Width of the grid cells.
 */
void get_residual_array(const double *pot, const double *rho, int cdim[3], double multiplier, double *residual, double delta) {
  const int N = cdim[0];
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0; 
  for (int k=0; k<N; k++){
    int k_plus = (k+1) % N;
    int k_min = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      int j_plus = (j+1) % N;
      int j_min = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        int i_plus = (i+1) % N;
        int i_min = (i-1>=0) ? (i-1) % N : (i-1) % N + N;
        residual[cell_getid(cdim,i,j,k)] = ((pot[cell_getid(cdim,i_min,j,k)]+ pot[cell_getid(cdim,i_plus,j,k)]+
                                                pot[cell_getid(cdim,i,j_min,k)]+pot[cell_getid(cdim,i,j_plus,k)]+
                                                pot[cell_getid(cdim,i,j,k_min)]+pot[cell_getid(cdim,i,j,k_plus)]
                                                - 6.0*pot[cell_getid(cdim,i,j,k)])/(delta*delta) - multiplier*(rho[cell_getid(cdim,i,j,k)]-rhs_correction));
      }
    }
  }
}

/**
 * @brief Get the residual of the current approximation to the potential on the grid
 *
 * Calculates the difference between the Laplacian of the potential and the right 
 * hand side of the Poisson equation.
 *
 * @param pot Array containing the potential on the grid.
 * @param rho Array containing the density on the grid.
 * @param cdim 3D size of the grid.
 * @param multiplier Multiplicative term on the right hand side of the Poisson equation.
 * @param delta Width of the grid cells.
 */
double get_residual(const double *pot, const double *rho, int cdim[3], double multiplier, double delta){
  double residual = 0.;
  int N = cdim[0];
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  for (int k=0; k<N; k++){
    int k_plus = (k+1) % N;
    int k_min = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      int j_plus = (j+1) % N;
      int j_min = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        int i_plus = (i+1) % N;
        int i_min = (i-1>=0) ? (i-1) % N : (i-1) % N + N;
        double res = ((pot[cell_getid(cdim,i_min,j,k)]+ pot[cell_getid(cdim,i_plus,j,k)]+
                                                pot[cell_getid(cdim,i,j_min,k)]+pot[cell_getid(cdim,i,j_plus,k)]+
                                                pot[cell_getid(cdim,i,j,k_min)]+pot[cell_getid(cdim,i,j,k_plus)]
                                                - 6.0*pot[cell_getid(cdim,i,j,k)])/(delta*delta) - multiplier*(rho[cell_getid(cdim,i,j,k)]-rhs_correction));
      residual += res*res;
      }
    }
  }
  return sqrt(residual/(N*N*N));
}

/**
 * @brief Get a new approximation to the potential on the grid.
 *
 * Performs a red-black sweep of the grid: first one half of the cells is
 * updated and next the other half.
 *
 * @param pot Array containing the potential on the grid.
 * @param rho Array containing the density on the grid.
 * @param cdim 3D size of the grid.
 * @param multiplier Multiplicative term on the right hand side of the Poisson equation.
 * @param delta Width of the grid cells.
 */
void perform_red_black_sweep(double *pot, const double *rho, int cdim[3], double multiplier, double delta) {
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  for (int col=0; col<2; col++){
    for (int k=0; k<cdim[2]; k++){
      int k_plus = (k+1) % cdim[2];
      int k_min = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
      for (int j=0; j<cdim[1]; j++){
        int j_plus = (j+1) % cdim[1];
        int j_min = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
        for (int i=0; i<cdim[0]; i++) {
          int i_plus = (i+1) % cdim[0];
          int i_min = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
          if ((i+j+k)%2 != col)
            continue;
          pot[cell_getid(cdim,i,j,k)] = ((pot[cell_getid(cdim,i_min,j,k)]+ pot[cell_getid(cdim,i_plus,j,k)]+ 
                                              pot[cell_getid(cdim,i,j_min,k)]+pot[cell_getid(cdim,i,j_plus,k)] + 
                                              pot[cell_getid(cdim,i,j,k_min)]+pot[cell_getid(cdim,i,j,k_plus)] - 
                                              delta*delta*multiplier*(rho[cell_getid(cdim,i,j,k)]-rhs_correction))/6.0);
          
        }
      }
    }
  }
}

/**
 * @brief Set initial guess for the potential
 *
 * @param pot Array of potential values on the grid.
 * @param cdim 3D size of the grid (equal in all dimensions).
 * @param MG Are we doing an MG calculation?
 */
void set_initial_guess(double *pot, const int cdim[3], int MG) {
  //Perhaps change this for the non-f(R) multigrid schemes still
  const int N = cdim[0];
  srand(time(NULL));
  for (int i = 0; i<N*N*N; i++){
    pot[i] = 0.; //Corresponds to fR = bar(fR) everywhere.
  }
}

void get_pm_potential(struct cic_mapper_data* data, const int N, const double box_size, struct threadpool* tp, int cdim[3]) {
  const int N_half = N / 2;
  double* rho_copy = malloc(N*N*N *sizeof(double));
  //rho_copy = data->rho;
  memcpy(rho_copy, data->rho, N*N*N* sizeof(double));
  if (rho_copy == NULL) error("Error allocating memory for density mesh");

  /* Allocates some memory for the mesh in Fourier space */
  fftw_complex* restrict frho =
      (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N_half + 1));
  if (frho == NULL)
    error("Error allocating memory for transform of density mesh");
  memuse_log_allocation("fftw_frho", frho, 1,
                        sizeof(fftw_complex) * N * N * (N_half + 1));

  /* Prepare the FFT library */
  fftw_plan forward_plan = fftw_plan_dft_r2c_3d(
      N, N, N, rho_copy, frho, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
  fftw_plan inverse_plan = fftw_plan_dft_c2r_3d(
      N, N, N, frho, rho_copy, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  /* Fourier transform to go to magic-land */
  fftw_execute(forward_plan);
  double r_s = 0.0;

  const int deconvolve = 0;
  const int discrete_symbol = 1;

  message("Going to apply Green function");
  /* Now de-convolve the CIC kernel and apply the Green function */
  mesh_apply_Green_function(tp, frho, /*slice_offset=*/0, /*slice_width=*/N,
                            /* mesh_size=*/N, r_s, box_size, deconvolve, discrete_symbol);        
  fftw_execute(inverse_plan);

  double conversion_factor;
  if (discrete_symbol)
    conversion_factor = (M_PI*M_PI/(N*N));
  else
    conversion_factor = 1.;
  for (int i =0; i<N*N*N; i++)
    rho_copy[i] = rho_copy[i] * conversion_factor;
  
  double sum=0.0;
  for (int i=0; i<N*N*N; i++) {
    sum += fabs(rho_copy[i]);
  }
  message("The mean absolute value of the potential is %lf", sum/(N*N*N));

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  free(rho_copy);
}

/**
 * @brief Compute the mean density on the mesh.
 *
 * If this function is called to solve the f_R field equation, 
 * the mean density is subtracted from the density array such that this
 * represents the overdensity. Else, the density is normalised such that it has
 * 1 as mean.
 *
 * @param rho Array containing the density in each cell of the grid.
 * @param N Size (1D) of the grid.
 * @param MG Are we doing an MG calculation?
 */
double get_mean_density(double *rho, const int N, int MG) {
  double sum = 0.;
  double mean_density = 0.;

  /* Find the mean density */
  for (int i = 0; i < N*N*N; i++) {
    sum += rho[i];
  }
  mean_density = sum/(N*N*N);

  /* Normalise to mean 1 or 0 */
  for (int i = 0; i < N*N*N; i++) {
    rho[i] = rho[i] / mean_density;
    if (MG) rho[i] -= 1.;
  }

  return mean_density;
}

/**
 * @brief Compute the field u = ln(f_R/mean(f_R)) on a series of uniform grids. 
 *
 * The FAS algorithm is used in combination with a multigrid scheme. At every level 
 * a solution for u is found using V-cycles. The initial guess for u on the next 
 * finer grid is set by prolongating u from a grid coarser.
 *
 * @param s The #space.
 * @param rho The density or overdensity on the grid.
 * @param u The field values u on the grid to solve for. 
 * @param MG Contains variables relevant for the MG calculation. 
 * @param N_min 1D size of the smallest grid on which to solve for the potential.
 * @param N_max 1D size of the largest grid on which to solve for the potential.
 * @param test Are we testing? (1 = uniform density, 2 = single particle, 3 = sine wave, 4 = two particles)
 */
void space_get_fR_contribution(const struct space *s, struct threadpool *tp, double *rho, double *u, struct MG_variables *MG, int N_min, const int N_max, const int test) {
  message("Doing test case %d", test);
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double fac = ((double) N_max)/box_size;
  int V_max = 15; //Maximum number of V-cycles

  /* Check if the space and mesh satisfy the conditions for Gauss-Seidel */
  if (dim[0] != dim[1] || dim[0] != dim[2]) {
    error("Noncubic domain for density calculation. Did not fill density array");
  }
  if (N_max%2 != 0){
    error("Uneven domain for density calculation. Did not apply Gauss-Seidel");
  }

  /* Get array of the different grid sizes */
  int n = N_min;
  int level = 0;
  int grid_sizes[32];
  while (n<=N_max) {
    grid_sizes[level++] = n;
    n *= 2;
  }

  const int N_levels = level;
  double *rho_levels[N_levels];
  double *u_levels[N_levels];
  double mean_density[N_levels];

  /* Store information of the largest grid */
  rho_levels[N_levels-1] = rho;
  u_levels[N_levels-1] = u;

  /* Decide if we need to calculate the mean density */
  if (!MG->overdensity) {
    mean_density[N_levels-1] = get_mean_density(rho_levels[N_levels-1], grid_sizes[N_levels-1], 1);
    if (test != 2 && test != 4) mean_density[N_levels-1] *= fac*fac*fac;
    message("The mean density is %lf", mean_density[N_levels-1]);
  }
  else {
    mean_density[N_levels-1] = 0.;
  }

  struct cic_mapper_data data;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];

  /* Initialise arrays for density and potential calculation. Except on the largest grid; there we already have them */
  for (int i=0; i<N_levels-1; i++) {
    rho_levels[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    u_levels[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    
    int N = grid_sizes[i];
    int cdim[3] = {N,N,N};
    double fac_level = grid_sizes[i]/box_size;
    double mean_density_set = 20.;

    /* Decide on density assignment based on the test. If we are not testing, just do CIC interpolation with the parts */
    switch (test) {
      case 1:
        /* Change the density field to homogeneous */
        for (int j=0; j<N*N*N; j++) {
          rho_levels[i][j] = mean_density[N_levels-1] /(fac*fac*fac);
        }
        break;
      case 2:
        /* Change the density field to represent a single particle at the centre of the box */
        for (int i2=0; i2<N; i2++) {
          for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
              if (i2==N/2 && j==N/2 && k==N/2) continue;
              rho_levels[i][cell_getid(cdim, i2, j, k)] = mean_density_set * (1. - 1e-4);
            }
          }
        }
        rho_levels[i][cell_getid(cdim, N/2, N/2, N/2)] = mean_density_set * (1. + 1e-4*(N*N*N-1.));
        break;
      case 3:
        /* Change the density field to represent a 1D sinusoid in the box */
        double fR_mod = MG->fR_bar;
        double delta = s->dim[0]/N;
        for (int i2=0; i2<N; i2++) {
          for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
              double x_dist = ((double) i2) * delta;
              rho_levels[i][cell_getid(cdim, i2, j, k)] = peak_overdensity(MG, x_dist, fR_mod, s->dim[0]);
            }
          }
        }
        break;
      case 4:
        /* Change the density field to represent two particles aligned along the x-axis of the box */
        for (int i2=0; i2<N; i2++) {
          for (int j=0; j<N; j++) {
            for (int k=0; k<N; k++) {
              if (i2==N/2 + 8*(i+1) && j==N/2 && k==N/2) continue;
              if (i2==N/2 - 8*(i+1) && j==N/2 && k==N/2) continue;
              rho_levels[i][cell_getid(cdim, i2, j, k)] = mean_density_set * (1. - 1e-4);
            }
          }
        }
        rho_levels[i][cell_getid(cdim, N/2 + 8*(i+1), N/2, N/2)] = mean_density_set * (1. + 1e-4*(N*N*N-2.));
        rho_levels[i][cell_getid(cdim, N/2 - 8*(i+1), N/2, N/2)] = mean_density_set * (1. + 1e-4*(N*N*N-2.));
        break;
      default:
        /* CIC assignment based on the particles */
        data.N = N;
        data.rho = rho_levels[i];
        data.fac = fac_level;
        gpart_to_mesh_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);        
    }

    /* Decide if we need to calculate the mean density */
    if (!MG->overdensity) {
      mean_density[i] = get_mean_density(rho_levels[i], grid_sizes[i], 1);
      if (test != 2 && test != 4) mean_density[i] *= fac_level*fac_level*fac_level;
    }
    else {
      for (int j=0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++) {
      }
      mean_density[i] = 0.;
    }
  }
  
  /* Set initial guess on the coarsest grid and solve directly */
  int cdim[3] = {N_min, N_min, N_min};
  set_initial_guess(u_levels[0], cdim, /*MG=*/1);
  apply_NGS(tp, rho_levels[0], u_levels[0], MG, cdim, mean_density[0], box_size);
  message("Finished applying Newton-Gauss-Seidel");

  /* Solve on all finer grids by prolongating from the previous grid and performing the multigrid method */
  for (int i = 1; i<N_levels; i++) {
    message("Going to the level with grid size %d \n", grid_sizes[i]);
    double mean_density_copy[i+1]; //Copy of mean_density containing the densities of the grid at level i to level 0
    for (int j=0; j<i+1; j++) {
      mean_density_copy[j] = mean_density[i-j];
    }
    prolongate_solution(u_levels[i-1], u_levels[i], grid_sizes[i-1], grid_sizes[i]); 
    
    int N = grid_sizes[i];
    cdim[0] = N;
    cdim[1] = N;
    cdim[2] = N;
    apply_multigrid_fR(tp, rho_levels[i], u_levels[i], MG, cdim, mean_density_copy, box_size, N_min, N, V_max);
  }
}

/**
 * @brief Directly solve for the potential with Newton-Gauss-Seidel relaxation.
 * 
 * Only one grid is used and no V-cycles are done.
 * 
 * @param rho Array containing the density on the grid.
 * @param u The field values u on the grid to solve for. 
 * @param MG Contains variables relevant for the MG calculation. 
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param box_size Side length of the box.
 */
void apply_NGS(struct threadpool *tp, const double *rho, double *u, struct MG_variables *MG, int cdim[3], double mean_density, double box_size) {
  message("Applying Newton-Gauss-Seidel...");
  int N = cdim[0];
  double delta = box_size/N;
  double residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density, delta);
  sleep(5);
  message("Before smoothing the residual is %E", residual);

  /* Get the mean */
  double sum = 0.;
  for (int i=0; i<N*N*N; i++) {
    sum += MG->fR_bar * exp(u[i]);
  }
  message("Before smoothing the mean is %E", sum/(N*N*N));

  double tolerance = 10e-9; //Choose reasonable value here
  int counter = 0;

  while (residual >= tolerance) {
    perform_red_black_sweep_fR(u, rho, MG, cdim, mean_density, delta);
    residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density, delta);
    if (counter%50 == 0) message("The residual is %E", residual);
    counter +=1;
  }
  message("We needed %d steps to converge", counter);
  sleep(5);
}

/**
 * @brief Get a new approximation to the field u on the grid.
 *
 * Performs a red-black sweep of the grid: first one half of the cells is
 * updated and next the other half.
 *
 * @param u Array containing the approximation to u on the grid.
 * @param rho Array containing the density on the grid.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param delta Width of the grid cells.
 */
void perform_red_black_sweep_fR(double *u, const double *rho, struct MG_variables *MG, int cdim[3], double mean_density, double delta) {
  int nbs[6];

  for (int col=0; col<2; col++){
    for (int k=0; k<cdim[2]; k++){
      nbs[4] = (k+1) % cdim[2];
      nbs[5] = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
      for (int j=0; j<cdim[1]; j++){
        nbs[2] = (j+1) % cdim[1];
        nbs[3] = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
        for (int i=0; i<cdim[0]; i++) {
          nbs[0] = (i+1) % cdim[0];
          nbs[1] = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
          if ((i+j+k)%2 != col)
            continue;
          
          /* Do we use the equation with density or with overdensity? */
          double density_term;
          if (MG->overdensity) density_term = 8. * M_PI * MG->G * rho[cell_getid(cdim, i,j,k)]/MG->a;
          else density_term = 8. * M_PI * MG->G * mean_density * rho[cell_getid(cdim, i,j,k)]/MG->a;
          
          double Laplacian_exp = get_Laplacian(MG, u, cdim, nbs, i, j, k);
          double field_term = MG->a*MG->a*MG->R* (1. - exp(-(1./2.)*u[cell_getid(cdim, i,j,k)]));
          double derivative_term = get_derivative(u, MG, cdim, nbs, i, j, k, delta);

          u[cell_getid(cdim,i,j,k)] -= (Laplacian_exp/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term + density_term))/derivative_term;
        }
      }
    }
  }
}

/**
 * @brief Get the derivative of the Laplacian operator with respect to u_{i,j,k}
 * 
 * @param u The approximation to the field u on the grid. 
 * @param MG Contains variables relevant for the MG calculation. 
 * @param cdim 3D size of the grid.
 * @param nbs Encodes for the six neighbours x+1, x-1, y+1, etc. of the cell (i,j,k).
 * @param i x-index of the current grid cell.
 * @param j y-index of the current grid cell.
 * @param k z-index of the current grid cell.
 * @param delta Width of a grid cell.
 */
double get_derivative(double *u, struct MG_variables *MG, int cdim[3], int nbs[6], int i, int j, int k, double delta) {
  double exp_term = (exp(u[cell_getid(cdim, nbs[0], j, k)]) + exp(u[cell_getid(cdim, nbs[1], j, k)]) + exp(u[cell_getid(cdim, i, nbs[2], k)])
                      + exp(u[cell_getid(cdim, i, nbs[3], k)]) + exp(u[cell_getid(cdim, i, j, nbs[4])]) + exp(u[cell_getid(cdim, i, j, nbs[5])])
                      + 6. * exp(u[cell_getid(cdim, i, j, k)]));
  double field_term = exp(u[cell_getid(cdim, i, j, k)]) * (u[cell_getid(cdim, nbs[0], j, k)] + u[cell_getid(cdim, nbs[1], j, k)] 
                      + u[cell_getid(cdim, i, nbs[2], k)] + u[cell_getid(cdim, i, nbs[3], k)] + u[cell_getid(cdim, i, j, nbs[4])] 
                      + u[cell_getid(cdim, i, j, nbs[5])] - 6. * u[cell_getid(cdim, i, j, k)]);
  double model_term = (1./2.) * MG->R*MG->a*MG->a * (1./(3.*MG->c*MG->c*MG->fR_bar)) * exp((-1./2.) * u[cell_getid(cdim, i, j, k)]);

  return (1./(2.*delta*delta)) * (field_term - exp_term) + model_term;
}

struct residual_mapper_data {
  const double *rho;
  const double *u;
  int N;
  struct MG_variables *MG;
  double *residual;
  double delta;
  double mean_density;
};

struct mesh_plane {
  int k;
};

void residual_mapper(void* map_data, int num, void* extra) {
  /* Unpack the shared information */
  struct residual_mapper_data *data = (struct residual_mapper_data*)extra;
  const double *rho = data->rho;
  const double *u = data->u;
  const int N = data->N;
  struct MG_variables *MG = data->MG;
  double *residual = data->residual;
  const double delta = data->delta;
  const double mean_density = data->mean_density;

  int cdim[3] = {N, N, N};

  /* Pointer to the chunk to be processed */
  const struct mesh_plane *planes = (const struct mesh_plane*)map_data;

  /* Loop over the elements assigned to this thread */
  for (int p=0; p<num; p++) {
    const int k = planes[p].k; //Which plane are we working on?
    message("Doing k-plane %d", k);
    int nbs[6];
    nbs[4] = (k+1) % N;
    nbs[5] = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      nbs[2] = (j+1) % N;
      nbs[3] = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        nbs[0] = (i+1) % N;
        nbs[1] = (i-1>=0) ? (i-1) % N : (i-1) % N + N;

        /* Do we use the equation with density or with overdensity? */
        double density_term;
        if (MG->overdensity) density_term = 8. * M_PI * MG->G * rho[cell_getid(cdim, i,j,k)]/MG->a;
        else density_term = 8. * M_PI * MG->G * mean_density * rho[cell_getid(cdim, i,j,k)]/MG->a;

        double Laplacian_exp = get_Laplacian(MG, u, cdim, nbs, i, j, k);
        double field_term = MG->a*MG->a*MG->R* (1. - exp(-(1./2.) * u[cell_getid(cdim, i,j,k)]));
        double res = Laplacian_exp/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term + density_term); 

        atomic_add_d(residual, (res*res)/(N*N*N));
      }
    }
  }
  
}

/**
 * @brief Get the residual of the current approximation to the field u on the grid
 *
 * Calculates the difference between the elliptic operator applied to u and the right 
 * hand side (density term) of the fR equation.
 *
 * @param u Array containing the approximation to u on the grid.
 * @param rho Array containing the density on the grid.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param delta Width of the grid cells.
 * @param verbose Are we talkative?
 */
double get_residual_fR(struct threadpool *tp, const double *u, const double *rho, struct MG_variables *MG, int cdim[3], double mean_density, double delta) {
  double residual = 0.;
  int N = cdim[0];

  struct residual_mapper_data data;
  data.rho = rho;
  data.u = u;
  data.N = N;
  data.MG = MG;
  data.residual = &residual;
  data.delta = delta;
  data.mean_density = mean_density;

  struct mesh_plane *planes = malloc(N*sizeof(struct mesh_plane));
  for (int k=0; k<N; k++) {
    planes[k].k = k;
  }

  threadpool_map(tp, residual_mapper, planes, N, sizeof(struct mesh_plane), threadpool_auto_chunk_size, (void*)&data);
  free(planes);
  return sqrt(residual);
}

/**
 * @brief Get the residual of the current approximation to the field u on the grid
 * at every grid point.
 *
 * Calculates the difference between the elliptic operator applied to u and the right 
 * hand side (density term) of the fR equation.
 *
 * @param u Array containing the approximation to u on the grid.
 * @param rho Array containing the density on the grid.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param delta Width of the grid cells.
 * @param verbose Are we talkative?
 */
void get_residual_array_fR(const double *u, const double *rho, struct MG_variables *MG, int cdim[3], double mean_density, double *residual_array, double delta) {
  int N = cdim[0];
  int nbs[6];

  for (int k=0; k<N; k++){
    nbs[4] = (k+1) % N;
    nbs[5] = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      nbs[2] = (j+1) % N;
      nbs[3] = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        nbs[0] = (i+1) % N;
        nbs[1] = (i-1>=0) ? (i-1) % N : (i-1) % N + N;

        /* Do we use the equation with density or with overdensity? */
        double density_term;
        if (MG->overdensity) density_term = 8. * M_PI * MG->G * rho[cell_getid(cdim, i,j,k)]/MG->a;
        else density_term = 8. * M_PI * MG->G * mean_density * rho[cell_getid(cdim, i,j,k)]/MG->a;

        double Laplacian_exp = get_Laplacian(MG, u, cdim, nbs, i, j, k);
        double field_term = MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*u[cell_getid(cdim, i,j,k)]));
        double res = Laplacian_exp/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term + density_term); 

        residual_array[cell_getid(cdim, i, j, k)] = res;
      }
    }
  }
}

/**
 * @brief Get the Laplacian operator of e^u following Oyaizu (2008)
 *
 * @param MG Contains variables relevant for the MG calculation.
 * @param u Array containing the approximation to u on the grid.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on the grid.
 * @param delta Width of the grid cells.
 * @param verbose Are we talkative?
 */
double get_Laplacian(struct MG_variables *MG, const double *u, int cdim[3], int nbs[6], int i, int j, int k) {
  double half_exp[6];
  half_exp[0] = (1./2.) * (exp(u[cell_getid(cdim, nbs[0], j, k)]) + exp(u[cell_getid(cdim, i, j, k)]));
  half_exp[1] = (1./2.) * (exp(u[cell_getid(cdim, i, j, k)]) + exp(u[cell_getid(cdim, nbs[1], j, k)]));

  half_exp[2] = (1./2.) * (exp(u[cell_getid(cdim, i, nbs[2], k)]) + exp(u[cell_getid(cdim, i, j, k)]));
  half_exp[3] = (1./2.) * (exp(u[cell_getid(cdim, i, j, k)]) + exp(u[cell_getid(cdim, i, nbs[3], k)]));

  half_exp[4] = (1./2.) * (exp(u[cell_getid(cdim, i, j, nbs[4])]) + exp(u[cell_getid(cdim, i, j, k)]));
  half_exp[5] = (1./2.) * (exp(u[cell_getid(cdim, i, j, k)]) + exp(u[cell_getid(cdim, i, j, nbs[5])]));

  double i_comp = half_exp[0] * u[cell_getid(cdim, nbs[0], j, k)] - u[cell_getid(cdim, i, j, k)] * (half_exp[0]+half_exp[1]) + half_exp[1] * u[cell_getid(cdim, nbs[1], j, k)];
  double j_comp = half_exp[2] * u[cell_getid(cdim, i, nbs[2], k)] - u[cell_getid(cdim, i, j, k)] * (half_exp[2]+half_exp[3]) + half_exp[3] * u[cell_getid(cdim, i, nbs[3], k)];
  double k_comp = half_exp[4] * u[cell_getid(cdim, i, j, nbs[4])] - u[cell_getid(cdim, i, j, k)] * (half_exp[4]+half_exp[5]) + half_exp[5] * u[cell_getid(cdim, i, j, nbs[5])];

  return i_comp + j_comp + k_comp;
}

/**
 * @brief Compute potential on the grid using multigrid acceleration.
 *
 * Solves for the field u on the finest grid, passes the solution
 * and residual to the coarser grids and recursively solves for it 
 * on every level using V-cycles and the FAS algorithm. When a solution
 * is found on a coarse grid, we prolongate u_H - R(u_h).
 *
 * @param rho Array containing the density on the finest grid.
 * @param u Array containing (approximate) u values on the finest grid.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param mean_density Mean density on all grids.
 * @param box_size Side length of the box.
 * @param N_min 1D size of the coarsest grid.
 * @param N_max 1D size of the finest grid.
 * @param V_max Maximum number of V-cycles that may be performed.
 */
void apply_multigrid_fR(struct threadpool *tp, const double *rho, double *u, struct MG_variables *MG, int cdim[3], const double *mean_density, const double box_size, const int N_min, const int N_max, const int V_max) {
  message("Applying the multigrid method for the grid with size %d...", N_max);

  /* Allocate the memory for the residual array on the finest level */
  double *residual_array = NULL;
  residual_array = (double*)malloc(sizeof(double) * N_max * N_max * N_max);
  if (residual_array == NULL){
    error("Error allocating memory for the residual array.");
  }
  memuse_log_allocation("residual.array", residual_array, 1, sizeof(double)*N_max*N_max*N_max);

  double field_sum = 0;
  for (int i =0; i<N_max*N_max*N_max; i++) {
    field_sum += fabs(u[i]);  
  }
  message("The msq of the field is %lf", field_sum/(N_max*N_max*N_max));

  double delta = box_size/N_max; //Width of a grid cell of the finest level
  double residual; 
  residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
  message("The first residual is %E", residual);
  double tolerance = 10e-9; //Choose reasonable value here
  int counter = 0;
  int fine_steps = 10; //Choose reasonable value here

  int V_cycles=0;
  int depth = 0;

  while (residual > tolerance && V_cycles < V_max) { 
    message("Performing V-cycle %d", V_cycles);
    /* Pre-smoothing */
    for (int i=0; i<fine_steps; i++) {
      perform_red_black_sweep_fR(u, rho, MG, cdim, mean_density[0], delta);
      residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
      if (residual<tolerance) break;
    }
    residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
    get_residual_array_fR(u, rho, MG, cdim, mean_density[0], residual_array, delta);
    if (residual<tolerance) break;
    
    /* Transfer residual array to get coarse-grid correction */
    message("After pre-smoothing the residual is %E. Going to recurse with V-cycles.", residual);
    FAS_recursive(u, residual_array, MG, cdim, delta, N_min, &depth);

    /* Post-smoothing if needed */
    residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
    message("Back on the finest grid the residual is %E", residual);
    if (residual > tolerance) {
      for (int i=0; i<fine_steps; i++) {
        perform_red_black_sweep_fR(u, rho, MG, cdim, mean_density[0], delta);
        residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
      }
    }
    residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
    V_cycles +=1;
    message("After %d V-cycle(s) the residual is %E", V_cycles, residual);
  }

  message("Performed %d V-cycle(s) in total", V_cycles);
  
  /* Post-smoothing until convergence. Should not be necessary! */
  while (residual > tolerance) {
    perform_red_black_sweep_fR(u, rho, MG, cdim, mean_density[0], delta);
    residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
    counter +=1;
    message("The counter is %d and the residual %E", counter, residual);
  }
  residual = get_residual_fR(tp, u, rho, MG, cdim, mean_density[0], delta);
  message("Needed to do %d step(s) in post-smoothing after which the residual was %E", counter, residual);
  free(residual_array);
}

/**
 * @brief Compute recursively the solutions of the coarser-grid equations.
 *
 * Solves for the solution u_H of the coarser-grid equation L_H(u_H) = L_H(R(u_h)) + R(f_h - L_h(u_h)).
 * On the coarsest grid this equation is solved exactly. Otherwise, the NGS smoothing is applied for a number
 * of steps, after which the residual of this equation is passed to the next coarser grid to 
 * recursively find a better approximation. When a solution is found on a coarse grid, we prolongate u_H - R(u_h).
 *
 * @param u Array containing (approximate) u values on the finest grid.
 * @param residual Array containing the residual at every fine grid point.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param delta Width of a fine grid cell.
 * @param N_stop 1D size of the coarsest grid.
 * @param depth How many levels are we into the V-cycle?
 */
void FAS_recursive(double *u, const double *residual, struct MG_variables *MG, int cdim[3], double delta, const int N_stop, int *depth) {
  *depth += 1;
  int N = cdim[0]; //Grid size of the current level we are on
  delta = delta*2.0; //Cells are twice as big on the coarser grid
  N = N/2; 

  /* Array for storing R(L_h(u_h) - f_h), the restriction of the residual on the finer grid */
  double *restricted_residual = NULL;
  restricted_residual = (double*)malloc(sizeof(double) * N * N * N);
  if (restricted_residual == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.grid", restricted_residual, 1, sizeof(double)*N*N*N);

  /* Array for storing R(u_h), the restriction of the solution on the finer grid */
  double *restricted_solution = NULL;
  restricted_solution = (double*)malloc(sizeof(double) * N * N * N);
  if (restricted_solution == NULL) {
    error("Error allocating memory for the coarser-grid equation.");
  }
  memuse_log_allocation("coarser.equation", restricted_solution, 1, sizeof(double)*N*N*N);

  /* Array for storing u_H, the solution of the equation on the coarser grid */
  double *coarser_solution = NULL;
  coarser_solution = (double*)malloc(N * N * N*sizeof(double)); 
  if (coarser_solution == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.solution", coarser_solution, 1, sizeof(double)*N*N*N);

  /* Array for storing L_H - f_H, the residual of the (approximate) solution of the equation on the coarser grid */
  double *coarser_residual = NULL; 
  coarser_residual = (double*)malloc(sizeof(double) * N * N * N);
  if (coarser_residual == NULL){
    error("Error allocating memory for the new coarser residual.");
  }
  memuse_log_allocation("coarser.newresidual", coarser_residual, 1, sizeof(double)*N*N*N);

  /* Restrict residual and solution of the finer grid*/
  restrict_residual(restricted_residual, residual, cdim);
  restrict_residual(restricted_solution, u, cdim); 

  /* Set initial guess on the coarser grid to be R(u_h) */
  for (int i=0; i<N*N*N; i++) {
    coarser_solution[i] = restricted_solution[i];
  }

  int cdimH[] = {N, N, N}; 
  double tolerance = 10e-9; //Choose a reasonable value here
  int counter = 0;
  double coarser_residual_abs = get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
  message("The first residual on the grid with size %d is %E", N, coarser_residual_abs);

  /* Solve the equation exactly if we are on the coarsest grid */
  if (N==N_stop) {
    while (coarser_residual_abs >= tolerance) {
      perform_red_black_sweep_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta); 
      coarser_residual_abs = get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
      counter +=1;
      if (counter%50 == 0) message("Did %d steps and the residual is %E", counter, coarser_residual_abs);
    }
    coarser_residual_abs = get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
    message("The total number of steps on the coarsest grid is %d and the residual is %E", counter, coarser_residual_abs);

    /* Prepare array for prolongation */
    for (int i=0; i<N*N*N; i++) {
      coarser_solution[i] -= restricted_solution[i];
    }
    prolongate_residual(coarser_solution, u, cdimH);
  }

  /* Do some smoothing and proceed to coarser grids */
  else { 
    counter = 0;
    int coarse_steps = 25;
    for (int i=0; i<coarse_steps; i++) {
      perform_red_black_sweep_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH,delta); 
      coarser_residual_abs =  get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
      counter +=1;
      //if (coarser_residual_abs < tolerance) break;
    }
    get_residual_array_coarser(coarser_solution, restricted_residual, restricted_solution, coarser_residual, MG, cdimH, delta);

    coarser_residual_abs = get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
    message("The residual after pre-smoothing is %E", coarser_residual_abs);

    FAS_recursive(coarser_solution, coarser_residual, MG, cdimH, delta, N_stop, depth);

    /* Post-smoothing */
    coarser_residual_abs = get_residual_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta);
    message("Back on the grid with size %d the residual is %E", N, coarser_residual_abs);
    for (int i=0; i<coarse_steps; i++) {
      perform_red_black_sweep_coarser(coarser_solution, restricted_residual, restricted_solution, MG, cdimH, delta); 
    }
    /* Prepare array for prolongation */
    for (int i=0; i<N*N*N; i++) {
      coarser_solution[i] -= restricted_solution[i];
    }
    prolongate_residual(coarser_solution, u, cdimH);
  }

  /* The coarser-grid correction has now been added to the finer-grid solution for the potential, so discard used arrays. */
  free(restricted_residual);
  free(coarser_solution);
  free(coarser_residual);
  free(restricted_solution);

  *depth -= 1;

}

/**
 * @brief Compute the residual of the coarse-grid equation at every grid point.
 * 
 * We are considering the equation L_H(u_H) = L_H(R(u_h)) + R(f_h - L_h(u_h)).
 *
 * @param coarser_solution Array containing the approximate solution at the current level.
 * @param restricted_residual Array containing the fine grid residual restricted to the current level.
 * @param restricted_solution Array containing the fine grid approximate solution restricted to the current level.
 * @param coarser_residual Array containing the residual of the coarse-grid equation at every grid point.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param delta Width of a fine grid cell.
 */
void get_residual_array_coarser(const double *coarser_solution, const double *restricted_residual, const double *restricted_solution, double *coarser_residual, struct MG_variables *MG, int cdim[3], double delta) {
  int N = cdim[0];
  int nbs[6];

  for (int k=0; k<N; k++){
    nbs[4] = (k+1) % N;
    nbs[5] = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      nbs[2] = (j+1) % N;
      nbs[3] = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        nbs[0] = (i+1) % N;
        nbs[1] = (i-1>=0) ? (i-1) % N : (i-1) % N + N;

        double Laplacian_exp[2] = {get_Laplacian(MG, coarser_solution, cdim, nbs, i, j, k), get_Laplacian(MG, restricted_solution, cdim, nbs, i, j, k)};
        double field_term[2] = {MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*coarser_solution[cell_getid(cdim, i,j,k)])), MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*restricted_solution[cell_getid(cdim, i,j,k)]))};
        double res = (Laplacian_exp[0] - Laplacian_exp[1])/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term[0] - field_term[1]) + restricted_residual[cell_getid(cdim, i,j,k)];

        coarser_residual[cell_getid(cdim, i, j, k)] = res;
      }
    }
  }
}

/**
 * @brief Compute a new approximation to the coarse-grid equation..
 * 
 * We are considering the equation L_H(u_H) = L_H(R(u_h)) + R(f_h - L_h(u_h)) and use Newton-Gauss-Seidel relaxation
 * in combination with a red-black sweep.
 *
 * @param coarser_solution Array containing the approximate solution at the current level.
 * @param restricted_residual Array containing the fine grid residual restricted to the current level.
 * @param restricted_solution Array containing the fine grid approximate solution restricted to the current level.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param delta Width of a fine grid cell.
 */
void perform_red_black_sweep_coarser(double *coarser_solution, const double *restricted_residual, const double *restricted_solution, struct MG_variables *MG, int cdim[3], double delta) {
  int nbs[6];

  for (int col=0; col<2; col++){
    for (int k=0; k<cdim[2]; k++){
      nbs[4] = (k+1) % cdim[2];
      nbs[5] = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
      for (int j=0; j<cdim[1]; j++){
        nbs[2] = (j+1) % cdim[1];
        nbs[3] = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
        for (int i=0; i<cdim[0]; i++) {
          nbs[0] = (i+1) % cdim[0];
          nbs[1] = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
          if ((i+j+k)%2 != col)
            continue;
          
          double Laplacian_exp[2] = {get_Laplacian(MG, coarser_solution, cdim, nbs, i, j, k), get_Laplacian(MG, restricted_solution, cdim, nbs, i, j, k)};
          double field_term[2] = {MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*coarser_solution[cell_getid(cdim, i,j,k)])), MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*restricted_solution[cell_getid(cdim, i,j,k)]))};
          double derivative_term = get_derivative(coarser_solution, MG, cdim, nbs, i, j, k, delta);

          coarser_solution[cell_getid(cdim,i,j,k)] -= ((Laplacian_exp[0] - Laplacian_exp[1])/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term[0] - field_term[1]) + (restricted_residual[cell_getid(cdim, i, j, k)]))/derivative_term;
        }
      }
    }
  }
}

/**
 * @brief Compute the residual of the coarse-grid equation.
 * 
 * We are considering the equation L_H(u_H) = L_H(R(u_h)) + R(f_h - L_h(u_h)).
 *
 * @param coarser_solution Array containing the approximate solution at the current level.
 * @param restricted_residual Array containing the fine grid residual restricted to the current level.
 * @param restricted_solution Array containing the fine grid approximate solution restricted to the current level.
 * @param MG Contains variables relevant for the MG calculation.
 * @param cdim 3D size of the grid.
 * @param delta Width of a fine grid cell.
 */
double get_residual_coarser(const double *coarser_solution, const double *restricted_residual, const double *restricted_solution, struct MG_variables *MG, int cdim[3], double delta) {
  double residual = 0.;
  int N = cdim[0];
  int nbs[6];

  for (int k=0; k<N; k++){
    nbs[4] = (k+1) % N;
    nbs[5] = (k-1>=0) ? (k-1) % N : (k-1) % N + N;
    for (int j=0; j<N; j++){
      nbs[2] = (j+1) % N;
      nbs[3] = (j-1>=0) ? (j-1) % N : (j-1) % N + N;
      for (int i=0; i<N; i++) {
        nbs[0] = (i+1) % N;
        nbs[1] = (i-1>=0) ? (i-1) % N : (i-1) % N + N;

        double Laplacian_exp[2] = {get_Laplacian(MG, coarser_solution, cdim, nbs, i, j, k), get_Laplacian(MG, restricted_solution, cdim, nbs, i, j, k)};
        double field_term[2] = {MG->R*MG->a*MG->a* (1. - exp(-(1./2.)*coarser_solution[cell_getid(cdim, i,j,k)])), MG->R*MG->a*MG->a*(1. - exp(-(1./2.)*restricted_solution[cell_getid(cdim, i,j,k)]))};
        double res = (Laplacian_exp[0] - Laplacian_exp[1])/(delta*delta) + (1./(3.*MG->c*MG->c*MG->fR_bar)) * (field_term[0] - field_term[1]) + restricted_residual[cell_getid(cdim, i,j,k)]; 

        residual += res*res;
      }
    }
  }
  return sqrt(residual/(N*N*N));
}

double peak_overdensity(struct MG_variables *MG, double delta_x, double fR_mean, double box_size) {
  double period = (2.*M_PI)/box_size;
  double term1 = 3. * fR_mean * sin(period*delta_x) * period*period *(MG->c*MG->c);
  double term2 = MG->a*MG->a*MG->R * (sqrt((1./(2.-sin(period*delta_x)))) - 1.);
  return ((MG->a)*(term2 - term1))/(8.*M_PI*MG->G);
}