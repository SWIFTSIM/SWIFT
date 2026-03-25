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

//struct AMR_levels {
  //struct cell **cells;
  //double mean_density; 
  //int cell_count; 
  //int ghost_count;
  //int depth; 
  //int cdim;
//};

void free_gparts_in_cells(struct cell *c, int *level) {
  // Free gparts array if it exists
  if (c->grav.parts != NULL) {
    c->grav.parts = NULL;
    c->grav.count = 0;
  }
  // Recursively free gparts in daughter cells
  if (c->split) {
    //message("Cell split");
    for (int i = 0; i < 8; ++i) {
      if (c->progeny[i] == NULL) {
        message("Child %d does not exist", i);
        continue;
      } 
      //message("Going to remove parts from child %d", i);
      *level += 1;
      free_gparts_in_cells(c->progeny[i], level);
    }
  }
  *level -= 1;
}

void init_test_single_particle(struct engine *e, int desired_depth) {
  message("Reached function");
  /* Clean current array of particles */
  struct space *s = e->s;
  int nr_cells = s->nr_cells;
  free(s->gparts);
  for (int i=0; i<nr_cells; i++) {
    int level = 0;
    free_gparts_in_cells(&(s->cells_top[i]), &level);
  }
  s->gparts = NULL;
  s->nr_gparts = 0;  

  /* Set the single particle to the desired state */
  message("Getting particle");
  s->nr_gparts = 1;
  s->gparts = calloc(1, sizeof(struct gpart));
  struct gpart *part = &s->gparts[0];

  //Put it in the middle of the box
  message("Initialising particle");
  part->mass = 50.;
  part->x[0] = s->dim[0]/2;
  part->x[1] = s->dim[0]/2;
  part->x[2] = s->dim[0]/2;

  /* Let the cells rebuild */
  message("Going to rebuild");
  engine_rebuild(e, 1, 0);

  /* Find the cell with the particle and split it desired_depth times */
  for (int i=0; i<nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    c->split = 0;
    if (c->split) error("Cell is split when it should not be");
    c->maxdepth = 0;
    int nr_gparts = c->grav.count;
    if (nr_gparts > 0) { 
      message("The gpart was found in cell %d with location (%lf, %lf, %lf)", i, c->loc[0], c->loc[1], c->loc[2]);
      int curr_depth = 0;
      get_progeny(s, c, desired_depth, &curr_depth);
      c->maxdepth = desired_depth;
    }
  }

  /* Split neighbouring cells desired_depth times */
  int i_loc = 16;
  int j_loc = 16;
  int k_loc = 16;
  int offset[3] = {-1, 0, 1};
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
        struct cell *c = &s->cells_top[cell_getid(s->cdim, i_loc-offset[i], j_loc-offset[j], k_loc-offset[k])];
        message("Going to do the cell with loc (%lf, %lf, %lf)", c->loc[0], c->loc[1], c->loc[2]);
        if (i==1 && j==1 && k==1) continue;
        int curr_depth = 0;
        get_progeny(s, c, desired_depth, &curr_depth);
        c->maxdepth = desired_depth;
      }
    }
  }


  //sleep(15);
}

void get_progeny(struct space *s, struct cell *c, int desired_depth, int *curr_depth) {
  if (*curr_depth < desired_depth) {
    space_getcells(s, 8, c->progeny, 0, 0, 0);
    c->split = 1;
    for (int k=0; k<8; k++) {
      struct cell *child = c->progeny[k];
      child->loc[0] = c->loc[0];
      child->loc[1] = c->loc[1];
      child->loc[2] = c->loc[2];

      double width = c->width[0]/2;
      child->width[0] = width;
      child->width[1] = width;
      child->width[2] = width;

      if (k & 4) {
        child->loc[0] += child->width[0];
        c->progeny[k]->neighbours[1] = c->progeny[k-4];
        c->progeny[k-4]->neighbours[0] = c->progeny[k];
      }
      if (k & 2) {
        child->loc[1] += child->width[1];
        c->progeny[k]->neighbours[3] = c->progeny[k-2];
        c->progeny[k-2]->neighbours[2] = c->progeny[k];
      }
      if (k & 1) {
        child->loc[2] += child->width[2];
        c->progeny[k]->neighbours[5] = c->progeny[k-1];
        c->progeny[k-1]->neighbours[4] = c->progeny[k];
      }
      child->parent = c;
    }
    /* Put the particles in the cell */
    if (c->hydro.count != 0 || c->grav.count != 0 || c->stars.count != 0 ||
    c->black_holes.count != 0 || c->sinks.count != 0) {
      divide_particles(c, s);
    }
    for (int k=0; k<8; k++) {
      *curr_depth += 1;
      get_progeny(s, c->progeny[k], desired_depth, curr_depth);
    }
  }
  *curr_depth -= 1;
}

void space_get_AMR_density(struct space *s, struct engine *e, int level_check, int desired_depth) {
  const double box_size = s->dim[0];
  struct cell *cells_top = s->cells_top;
  //int maxdepth = s->maxdepth;

  int split;

  int cell_number=0;
  /* Find the largest tree */
  int max_depth =0;
  for (int i =0; i<s->nr_cells; i++) {
    if (cells_top[i].maxdepth>max_depth) {
      max_depth = cells_top[i].maxdepth;
      cell_number = i;
    }
  }

  message("The largest depth is %d and is found in cell number %d", max_depth, cell_number);
  
  if (level_check > max_depth) {
    message("Nog enough levels found to perform check. Set the max level to the max depth.");
    level_check = max_depth;
  }

  int nr_gparts = s->nr_gparts;
  message("There are supposed to be %d particles", nr_gparts);

    /* Initialise the patch structure */
  struct AMR_levels levels[max_depth+1];
  memset(levels, 0, sizeof levels);

  int cdim = s->cdim[0];

  /* Set the top level */
  levels[0].cells = malloc(s->nr_cells * sizeof(struct cell *));
  for (int i = 0; i < s->nr_cells; i++) {
    levels[0].cells[i] = &cells_top[i];
  }
  levels[0].cell_count = s->nr_cells;
  levels[0].depth = 0;
  levels[0].ghost_count = 0;
  levels[0].cdim = cdim;

  split=0;
  for (int j=0; j<levels[0].cell_count; j++) {
    if (levels[0].cells[j]->split) {
      split = 1;
      break;
    }
  }

  /* Set the other levels */
  int level_i=0;
  while(split && level_i <= max_depth) {
    cdim *= 2;
    levels[level_i+1].cells = NULL;
    levels[level_i+1].cell_count = 0;
    levels[level_i+1].depth = level_i+1;
    levels[level_i+1].cdim = cdim;
    extract_AMR_patches(s, cells_top, level_i+1, &(levels[level_i+1].cells), &(levels[level_i+1].cell_count));
    level_i++;

    split = 0;
    for (int j=0; j<levels[level_i].cell_count; j++) {
      if (levels[level_i].cells[j]->split) {
        split = 1;
        break;
      }
    }
    if (split == 1 && level_i > max_depth) error("Found a split cell at the max depth");
  }

  /* Re-adjust the maximum depth and the array */
  max_depth = level_i;
  
  /* Set the minimum depth */ 
  int min_depth = 0;
  for (int i=1; i<max_depth+1; i++) {
    int N_curr = levels[i].cdim;
    if (levels[i].cell_count < N_curr * N_curr * N_curr) break;
    min_depth++;
  }

  for (int i=0; i<min_depth;i++) {
    /* Sort the uniform grids in row major order based on the location */
    qsort(levels[i+1].cells, levels[i+1].cell_count, sizeof(struct cell *), sort_cell);
  }

  int old_cell_count[max_depth+1];

  /* Link the uniform levels */
  for (int i=0; i<min_depth+1;i++) {
    link_uniform_level(&levels[i]);
    old_cell_count[i] = levels[i].cell_count;
  }

  for (int i=min_depth+1; i<max_depth+1; i++) {
    message("Linking nonuniform level %d", i);
    link_nonuniform_level(s, &levels[i], 0, levels[i].cell_count);
    old_cell_count[i] = levels[i].cell_count;
  }

  int cells_old = 1;
  message("The total number of cells before is %d", s->tot_cells);
  int cells_new = 0;
  int passes = 0;
  int new_cell_count[max_depth+1];

  /* Do the smoothing */
  message("Smoothing tree");
  while (cells_old != cells_new) {
    cells_old = s->tot_cells;
    build_refinement_map(s, min_depth, max_depth, levels);
    modify_tree(s, min_depth, max_depth, levels, new_cell_count);
    cells_new = s->tot_cells;
    passes +=1;
    message("The total number of cells after %d passes is %d", passes, s->tot_cells);

    for (int i = min_depth; i <= max_depth; i++) {
      for (int j = 0; j < levels[i].cell_count; j++) {
        levels[i].cells[j]->refine = 0;
        //levels[i].cells[j]->refine2 = 0;
      }
    }

    for (int i=min_depth+1; i<max_depth+1; i++) {
      link_nonuniform_level(s, &levels[i], old_cell_count[i], new_cell_count[i]);
      old_cell_count[i] += new_cell_count[i];
    }
    /* Did any of the nonuniform levels become uniform? */
    message("The largest uniform grid was at %d", min_depth);
    int temp_depth = min_depth;
    for (int i=min_depth+1; i<max_depth; i++) {
      int N = levels[i].cdim;
      if (levels[i].cell_count == N*N*N) temp_depth +=1;
      else break;
    }
    for (int i=min_depth; i<temp_depth;i++) {
      /* Sort the newly uniform grids in row major order based on the location */
      qsort(levels[i+1].cells, levels[i+1].cell_count, sizeof(struct cell *), sort_cell);
    }
    min_depth = temp_depth;
  }

  /* Set ghost values all to zero */
  for (int i=min_depth; i<max_depth+1; i++) {
    for (int j=0; j<levels[i].cell_count; j++) {
      levels[i].cells[j]->ghost=0;
    }
  }

  //for (int i=0; i<levels[3].cell_count; i++) {
    //struct cell *c = levels[3].cells[i];
    //int x = (int) (c->loc[0]/c->width[0]);
    //int y = (int) (c->loc[1]/c->width[1]);
    //int z = (int) (c->loc[2]/c->width[2]);
    //if (x == 1 && y == 10 && z == 18) {
      //if (c->split) {
        //for (int j=0; j<8; j++) {
          //message("The child %d has address %p", j, c->progeny[j]);
        //}
      //}
      //else {
        //message("Cell not split!");
      //}
    //}
  //}
  //sleep(10);

  //for (int i=min_depth; i<max_depth;i++) {
    /* Sort the uniform grids in row major order based on the location */
    //message("Sorting level %d", i);
    //qsort(levels[i+1].cells, levels[i+1].cell_count, sizeof(struct cell *), sort_cell);
  //}

  /* Loop over particles and perform tree searches to get the density everywhere*/
  assign_densities(cells_top, &e->threadpool, s->nr_cells, box_size, level_check);

  //min_depth = 0;
  int max_gridsize = perform_uniform_calculation(s, min_depth, max_depth, levels);
  message("max_gridsize = %d", max_gridsize);
  //max_gridsize = 8;

  //int cdim2[3] = {32, 32, 32};
  //struct cell *c = levels[min_depth].cells[cell_getid(cdim2, 16, 16, 16)];
  //message("The potential of the cell at (%lf, %lf, %lf) is %lf", c->loc[0], c->loc[1], c->loc[2], c->CIC_potential);
  
  //max_gridsize = 8;

  /* Check that the densities are the same for the existing cells */
  //double fac = box_size/16;
  //int cdim2[3] = {16,16,16};
  
  //message("Exporting AMR density data");
  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/AMR_16_density.txt", "w");
  //for (int k=0; k<levels[1].cell_count; k++) {
    //for (int j=0; j<16; j++) {
      //for (int i=0; i<16; i++) {
        //if (i>=2 || j>=2 || k>=2) fprintf(fptr, "%.15g %.15g %.15g %.15g \n", levels[1].cells[cell_getid(cdim2, i,j,k)]->CIC_density, (double) i*fac, (double) j*fac, (double) k*fac);
        //fprintf(fptr, "%.15g %.15g %.15g %.15g \n", levels[1].cells[k]->CIC_density, levels[1].cells[k]->loc[0], levels[1].cells[k]->loc[1], levels[1].cells[k]->loc[2]);
      //}
    //}
  //}
  //fclose(fptr); 
  //sleep(10);

  //struct AMR_levels levels[max_depth+1];
  

  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/uniform_check/FMG_16_atcells.txt", "w");
  //for (int j=0; j<levels[1].cell_count; j++) {
    //fprintf(fptr, "%lf \n", levels[1].cells[j]->CIC_potential);
  //}
  //fclose(fptr);
  //message("Exported the FMG solution at cells");

  //sleep(10);

  perform_nonuniform_calculation(s, min_depth, max_depth, levels, max_gridsize, box_size);

  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/AMR_routine_tests/nonuniform_routine_uniform_grid/AMR_16m_50s.txt", "w");
  //for (int j=0; j<16*16*16; j++) {
    //fprintf(fptr, "%lf \n", levels[1].cells[j]->CIC_potential);
  //}
  //fclose(fptr);
  //message("Exported the AMR solution");

  //perform_nonuniform_calculation(s, cells_top, min_depth, max_depth, max_gridsize, box_size);

  //test_density_assignment(cells_top, s->nr_cells, box_size, nr_gparts, s);
  //message("Exporting AMR potential particle data");
  //FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/missing_1stoct_16.txt", "w");
  //for (int i=0; i<levels[1].cell_count; i++) {
    //fprintf(file_swift_density, "%.15g %.15g %.15g %.15g \n", levels[1].cells[i]->CIC_potential, levels[1].cells[i]->loc[0], levels[1].cells[i]->loc[1], levels[1].cells[i]->loc[2]);
  //}
  //fclose(file_swift_density);
  //message("Exported the cell potential data.");
  //sleep(5);

  get_cell_accelerations(s, min_depth, max_depth, levels);

  potential_to_fake_gparts(s, min_depth, max_depth, levels, desired_depth);
  
  //potential_to_gparts(s, min_depth, max_depth, levels);

  /* Free memory of the pointer array to the cells */
  for (int i=0; i<max_depth+1; i++) {
    free(levels[i].cells);
  }
}

/* Three-point stencil to get the accelerations in the cells */
void get_cell_accelerations(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]) {
  for (int i=min_depth; i<max_depth+1; i++) {
    for (int j=0; j<levels[i].cell_count; j++) {
      struct cell *c = levels[i].cells[j];
      if (c->ghost) {
        message("Accidentally accessing a ghost cell! Check cell numberings");
        continue;
      }
      for (int k=0; k<3; k++) {
        c->CIC_acc[k] = 0.;
      }
      c->CIC_acc[0] = c->neighbours[1]->CIC_potential - c->neighbours[0]->CIC_potential;
      c->CIC_acc[1] = c->neighbours[3]->CIC_potential - c->neighbours[2]->CIC_potential;
      c->CIC_acc[2] = c->neighbours[5]->CIC_potential - c->neighbours[4]->CIC_potential;
    }
  }
}

void potential_to_fake_gparts(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int desired_depth) {
  /* Add function to create 'fake' particles and add them to the cells */
  //desired_depth = 2;
  const int N_parts_old = s->nr_gparts;
  int N_parts_new = 30*9;
  s->nr_gparts += N_parts_new;
  //s->gparts = realloc(s->gparts, s->nr_gparts * sizeof(*s->gparts));

  /* Create array of particle locations */
  double part_loc[N_parts_new][3];
  double r_parts = 5.;
  for (int i=0; i<9; i++) {
    r_parts += 0.5 * (int) i;
    generate_particles(s, 30, part_loc, r_parts,30*i);
  }
  
  struct gpart_ref p_ref[N_parts_old + N_parts_new];

  message("Going to add the old particles");
  int counter_added= 0;
  for (int i=0; i<levels[desired_depth].cell_count; i++) {
    struct cell *c = levels[desired_depth].cells[i];
    if (c->grav.count > 0) {
      message("The grav count of cell %d is %d", i, c->grav.count);
      for (int j=0; j<c->grav.count; j++) {
        p_ref[counter_added].cell = c;
        p_ref[counter_added].index = j;
        p_ref[counter_added].level = 0;
        counter_added += 1;
      }
    }
  }
  if (counter_added>N_parts_old) error("Somehow added too many particles");

  /* Pass particle locations to the gparts */
  counter_added = 0;
  int depth = 0;
  for (int i=N_parts_old; i<N_parts_old + N_parts_new; i++) {
    for (int j=0; j<s->nr_cells; j++) {
      struct cell *c = &s->cells_top[j];
      //message("Accessing cell %d",j);
      double cx0 = c->loc[0];
      //message("Accessed cell %d", j);
      double cx1 = c->loc[0] + c->width[0];
      double cy0 = c->loc[1];
      double cy1 = c->loc[1] + c->width[1];
      double cz0 = c->loc[2];
      double cz1 = c->loc[2] + c->width[2];
      if (part_loc[i-N_parts_old][0] >= cx0 && part_loc[i-N_parts_old][0] < cx1 && part_loc[i-N_parts_old][1] >= cy0 && part_loc[i-N_parts_old][1] <cy1 && part_loc[i-N_parts_old][2] >= cz0 && part_loc[i-N_parts_old][2] < cz1) {
        if (i== 126) message("Found overlap for %d with cell at (%lf, %lf, %lf)", i, cx0, cy0, cz0);
        //if (j==16912) message("Placing particle at (%lf,%lf, %lf) in this cell", part_loc[i-N_parts_old][0], part_loc[i-N_parts_old][1], part_loc[i-N_parts_old][2]);
        const int nr_old = c->grav.count;
        c->grav.count += 1;
        if (c->split && depth<desired_depth) {
          if (i== 126) message("Going to recurse for %d", j);
          particle_to_cells_recursive(part_loc[i-N_parts_old], c->progeny, 8, &counter_added, N_parts_old, N_parts_new, p_ref, &depth, desired_depth);
        }
        else {
          if (i== 126) message("Not going to recurse for %d", j);
        //message("Adding particle at (%lf, %lf, %lf) to cell with location (%lf, %lf, %lf) and width (%lf, %lf, %lf)", part->x[0], part->x[1], part->x[2], c->loc[0], c->loc[1], c->loc[2], c->width[0], c->width[1], c->width[2]);
          if (nr_old == 0) c->grav.parts = malloc(sizeof(struct gpart));
          else {
            //if (j == 16912) message("We found %d particles in this cell", c->grav.count);
            c->grav.parts = realloc(c->grav.parts, (nr_old+1) * sizeof(struct gpart));
          }
          p_ref[N_parts_old + counter_added].cell = c;
          p_ref[N_parts_old + counter_added].index = c->grav.count - 1;
          p_ref[N_parts_old + counter_added].level = 0;
          c->grav.parts[nr_old].x[0] = part_loc[i-N_parts_old][0];
          c->grav.parts[nr_old].x[1] = part_loc[i-N_parts_old][1];
          c->grav.parts[nr_old].x[2] = part_loc[i-N_parts_old][2];
          //message("Added particle to cell. Particle %d in cell has location (%lf, %lf, %lf)", c->grav.count + 1, c->grav.parts[nr_old].x[0], c->grav.parts[nr_old].x[1], c->grav.parts[nr_old].x[2]);
          counter_added +=1;
        }        
        break;
      }
    }
  }
  message("Finished adding particles to cells. The added counter is %d", counter_added);

  //for (int i=0; i<N_parts_new + N_parts_old; i++) {
    //message("particle %d in the pointer array has loc (%lf, %lf, %lf)", i, gparts[i]->x[0], gparts[i]->x[1], gparts[i]->x[2]);
  //}

  /* Check if cells have more than one particle */
  /*for (int j=0; j<s->nr_cells; j++) {
    struct cell *c = levels[0].cells[j];
    if (c->grav.count > 1) {
      message("Cell %d has %d particles", j, c->grav.count);
    }
  }*/

  /* Pass potential to the particles (real and fake ones) */
  potential_to_gparts(s, min_depth, max_depth, levels);

  message("Writing accelerations to file");
  FILE *accelerations;
  accelerations = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/single_particle_test/trilinear/potential_5_fullnormalised_morecells_32.txt", "w");
  for (int i=0; i<N_parts_old + N_parts_new; i++) {
    struct gpart gpart = p_ref[i].cell->grav.parts[p_ref[i].index];
    //message("Going to write %d", i);
    //fprintf(accelerations, "%.15g %.15g %.15g %.15g %.15g %.15g \n", gparts[i]->a_grav_mesh[0], gparts[i]->a_grav_mesh[1], gparts[i]->a_grav_mesh[2], gparts[i]->x[0], gparts[i]->x[1], gparts[i]->x[2]);
    fprintf(accelerations, "%.15g %.15g %.15g %.15g \n", gpart.potential_mesh, gpart.x[0], gpart.x[1], gpart.x[2]);
  }
  fclose(accelerations);
  message("Done writing accelerations to file");
}

void particle_to_cells_recursive(double part_loc[3], struct cell **cells, int nr_cells, int *counter_added, int N_parts_old, int N_parts_new, struct gpart_ref p_ref[N_parts_old + N_parts_new], int *depth, int desired_depth) {
  //message("The number of cells is %d", nr_cells);
  *depth += 1;
  for (int i=0; i<nr_cells; i++) {
    struct cell *c = cells[i];
    //message("Accessing cell %d",i);
    double cx0 = c->loc[0];
    //message("Accessed cell %d", i);
    double cx1 = c->loc[0] + c->width[0];
    double cy0 = c->loc[1];
    double cy1 = c->loc[1] + c->width[1];
    double cz0 = c->loc[2];
    double cz1 = c->loc[2] + c->width[2];
    if (part_loc[0] >= cx0 && part_loc[0] < cx1 && part_loc[1] >= cy0 && part_loc[1] <cy1 && part_loc[2] >= cz0 && part_loc[2] < cz1) {
      if (*counter_added == 125) message("Found overlap for %d at (%lf, %lf, %lf) with the cell at (%lf, %lf, %lf)", *counter_added, part_loc[0], part_loc[1], part_loc[2], cx0, cy0, cz0);
      const int nr_old = c->grav.count;
      c->grav.count += 1;
      //message("Adding particle at (%lf, %lf, %lf) to cell with location (%lf, %lf, %lf) and width (%lf, %lf, %lf)", part->x[0], part->x[1], part->x[2], c->loc[0], c->loc[1], c->loc[2], c->width[0], c->width[1], c->width[2]);
      if (c->split && *depth<desired_depth) {
        particle_to_cells_recursive(part_loc, c->progeny, 8, counter_added, N_parts_old, N_parts_new, p_ref, depth, desired_depth);
      }
      else {
        if (nr_old == 0) c->grav.parts = malloc(sizeof(struct gpart));
        else c->grav.parts = realloc(c->grav.parts, (nr_old+1) * sizeof(struct gpart));
        c->grav.parts[nr_old].x[0] = part_loc[0];
        c->grav.parts[nr_old].x[1] = part_loc[1];
        c->grav.parts[nr_old].x[2] = part_loc[2];
        if (*counter_added == 125) message("Added particle to cell. Particle %d in cell has location (%lf, %lf, %lf)", c->grav.count + 1, c->grav.parts[nr_old].x[0], c->grav.parts[nr_old].x[1], c->grav.parts[nr_old].x[2]);
        p_ref[N_parts_old + *counter_added].cell = c;
        p_ref[N_parts_old + *counter_added].index = c->grav.count - 1;
        p_ref[N_parts_old + *counter_added].level = *depth;
        *counter_added +=1;
        //message("Added particle to cell. Particle %d in cell has location (%lf, %lf, %lf)", c->grav.count + 1, c->grav.parts[nr_old].x[0], c->grav.parts[nr_old].x[1], c->grav.parts[nr_old].x[2]);
      }
      break;
    }
  }
  *depth -= 1;
}

void generate_particles(struct space *s, int N_parts_new, double positions[N_parts_new][3], double r_parts, int start) {
  const double golden_angle = M_PI * (3.0 - sqrt(5.0));
  double centre_loc = s->dim[0]/2.;

  for (int i=0; i<N_parts_new; i++) {
    double cos_theta = 1.0 - 2.0 * ((double)i + 0.5) / (double)N_parts_new;
    double sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta * cos_theta));

    double phi = golden_angle * (double)i;

    double ux = cos(phi) * sin_theta;
    double uy = sin(phi) * sin_theta;
    double uz = cos_theta;

    positions[i+start][0] = centre_loc + r_parts * ux;
    positions[i+start][1] = centre_loc + r_parts * uy;
    positions[i+start][2] = centre_loc + r_parts * uz;
  }
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
/*
struct cell_basic {
  double loc[3];
  double width;
  double CIC_density; 
  double CIC_potential;
  double update_mask;
};*/

void link_nonuniform_level(struct space *s, struct AMR_levels *level, int start_index, int link_nr) {
  double fac = s->dim[0]/level->cdim;
  int cdim[3] = {level->cdim, level->cdim, level->cdim};

  for (int i=start_index;i<start_index + link_nr;i++) {
    /* Box wrap the multipole's position */
    const double search_loc_x = box_wrap(level->cells[i]->loc[0], 0., s->dim[0]);
    const double search_loc_y = box_wrap(level->cells[i]->loc[1], 0., s->dim[1]);
    const double search_loc_z = box_wrap(level->cells[i]->loc[2], 0., s->dim[2]);

    int cell_i = (int) ((search_loc_x+0.0001)/fac);
    int cell_j = (int) ((search_loc_y+0.0001)/fac);
    int cell_k = (int) ((search_loc_z+0.0001)/fac);

    int i_plus = (cell_i + 1) % cdim[0];
    int i_min = (cell_i-1>=0) ? (cell_i-1) % cdim[0] : (cell_i-1) % cdim[0] + cdim[0];
    int j_plus = (cell_j + 1) % cdim[0];
    int j_min = (cell_j-1>=0) ? (cell_j-1) % cdim[0] : (cell_j-1) % cdim[0] + cdim[0];
    int k_plus = (cell_k + 1) % cdim[0];
    int k_min = (cell_k-1>=0) ? (cell_k-1) % cdim[0] : (cell_k-1) % cdim[0] + cdim[0];

    int search_id[6] = {cell_getid(cdim, i_plus, cell_j, cell_k), cell_getid(cdim, i_min, cell_j, cell_k),
                        cell_getid(cdim, cell_i, j_plus, cell_k), cell_getid(cdim, cell_i, j_min, cell_k),
                        cell_getid(cdim, cell_i, cell_j, k_plus), cell_getid(cdim, cell_i, cell_j, k_min)};

    struct cell *parent = level->cells[i]->parent;

    int neighbour_counter = 0;
    int no_neighbour_counter = 0;

    for (int j=0; j<6; j++) {
      if (level->cells[i]->neighbours[j] == NULL) {
        if (cell_i==0 && cell_j==0 && cell_k==2) message("Searching for neighbour %d for cell %d", j, i);
        int pre_smoothing = 1;
        int neighbour = find_neighbour(level, parent, level->cells[i], search_id[j], fac, j, pre_smoothing, 0);
        if (cell_i==0 && cell_j==0 && cell_k==2 && neighbour) message("Found neighbour %d at (%lf, %lf, %lf)", j, level->cells[i]->neighbours[j]->loc[0], level->cells[i]->neighbours[j]->loc[1], level->cells[i]->neighbours[j]->loc[2]);
        no_neighbour_counter += 1;
        //message("For cell %d, neighbour %d was not linked", i, j);
      }
      else {
        neighbour_counter +=1;
        //message("For cell %d, neighbour %d already linked", i, j);
      }
    }
  }
}

void potential_to_gparts(struct space *s, int min_depth, int max_depth, struct AMR_levels *levels) {
  //message("Calling this function");
  /* Loop over the levels. At each level loop over the cells and then over the particles in each cell */
  //for (int i=max_depth-2; i<max_depth+1; i++) {
  for (int i=0; i<max_depth+1; i++) {
    struct AMR_levels *level = &(levels[max_depth-i]); //Start at deepest level
    //message("Doing level %d", max_depth-i);
    int nr_cells = level->cell_count; //Cell count should be excluding ghost cells
    message("The cell count is %d", nr_cells);
    for (int j=0; j<nr_cells; j++) {
      struct cell *c = level->cells[j];
      if (i==0 && !c->split && c->grav.count>0) message("Calling cell %d", j);
      if (c->split) {
      //if (c->split && level->depth<2) {
        if (i==0) message("Cell %d is split", j);
        continue;
      } ///Means we have already done CIC with the particles on another level
      int nr_gparts = c->grav.count;
      //message("The particle count is %d", nr_gparts);
      //if (nr_gparts>0) message("Going to do cell %d with loc (%lf, %lf, %lf) \n", j, c->loc[0], c->loc[1], c->loc[2]);
      for (int k=0; k<nr_gparts; k++) {
        struct gpart *gp = &(c->grav.parts[k]);
        message("Going to get the potential for the particle at (%lf, %lf, %lf)", gp->x[0], gp->x[1], gp->x[2]);
        get_AMR_potential(s, max_depth, max_depth-i, levels, gp, j);
        //get_AMR_potential(s, max_depth, max_depth-i, levels, gp, j);
        message("The assigned potential is %lf", c->grav.parts[k].potential_mesh);
      }
    }
  }
}

void get_AMR_potential(struct space *s, int max_depth, int current_depth, struct AMR_levels *levels, struct gpart *gp, int cell_nr) {
  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  double boxsize = dim[0];
  int stop_depth = current_depth;
  //message("Going to get the potential for the particle at (%lf, %lf, %lf)", gp->x[0], gp->x[1], gp->x[2]);
  /* Box wrap the particle's position */
  const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
  const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
  const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

  /* Some local accumulators */
  //double p = 0.;
  //double a[3] = {0.};

  /* CIC for the potential: Get overlap with the cell it is in and 7 others */
  struct cell *home_cell = levels[current_depth].cells[cell_nr];
  /* Check if the particle cloud lies entirely within this level */
  for (int i=0; i<current_depth+1; i++) {
    if (home_cell->neighbours[0] == NULL || home_cell->neighbours[2] == NULL || home_cell->neighbours[4] == NULL) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    else if (home_cell->neighbours[0]->ghost || home_cell->neighbours[2]->ghost || home_cell->neighbours[4]->ghost) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    else if (home_cell->neighbours[0]->neighbours[2] == NULL || home_cell->neighbours[0]->neighbours[4] == NULL || home_cell->neighbours[2]->neighbours[4] == NULL) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    else if (home_cell->neighbours[0]->neighbours[2]->ghost || home_cell->neighbours[0]->neighbours[4]->ghost || home_cell->neighbours[2]->neighbours[4]->ghost) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    else if (home_cell->neighbours[0]->neighbours[2]->neighbours[4] == NULL) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    else if (home_cell->neighbours[0]->neighbours[2]->neighbours[4]->ghost) {
      home_cell = levels[current_depth].cells[cell_nr]->parent;
      stop_depth -= 1;
      if (stop_depth < 0) error("Cells that should exist not found!");
      continue;
    }
    break;
  }
  message("Found the stop depth %d", stop_depth);
  /*if (stop_depth == 3) {
    message("Particle at (%lf, %lf, %lf) has home cell (%lf, %lf, %lf)", pos_x, pos_y, pos_z, home_cell->loc[0], home_cell->loc[1], home_cell->loc[2]);
    message("The neighbour 1 is at (%lf, %lf, %lf)", home_cell->neighbours[0]->loc[0], home_cell->neighbours[0]->loc[1], home_cell->neighbours[0]->loc[2]);
    message("The neighbour 2 is at (%lf, %lf, %lf)", home_cell->neighbours[2]->loc[0], home_cell->neighbours[2]->loc[1], home_cell->neighbours[2]->loc[2]);
    message("The neighbour 3 is at (%lf, %lf, %lf)", home_cell->neighbours[4]->loc[0], home_cell->neighbours[4]->loc[1], home_cell->neighbours[4]->loc[2]);
    message("The neighbour 4 is at (%lf, %lf, %lf)", home_cell->neighbours[2]->neighbours[4]->loc[0], home_cell->neighbours[2]->neighbours[4]->loc[1], home_cell->neighbours[2]->neighbours[4]->loc[2]);
    message("The neighbour 5 is at (%lf, %lf, %lf)", home_cell->neighbours[0]->neighbours[4]->loc[0], home_cell->neighbours[0]->neighbours[4]->loc[1], home_cell->neighbours[0]->neighbours[4]->loc[2]);
    message("The neighbour 6 is at (%lf, %lf, %lf)", home_cell->neighbours[0]->neighbours[2]->loc[0], home_cell->neighbours[0]->neighbours[2]->loc[1], home_cell->neighbours[0]->neighbours[2]->loc[2]);
    message("The neighbour 7 is at (%lf, %lf, %lf)", home_cell->neighbours[0]->neighbours[2]->neighbours[4]->loc[0], home_cell->neighbours[0]->neighbours[2]->neighbours[4]->loc[1], home_cell->neighbours[0]->neighbours[2]->neighbours[4]->loc[2]);
  }*/

  double x[3] = {pos_x, pos_y, pos_z};
  double width[3] = {home_cell->width[0], home_cell->width[1], home_cell->width[2]};
  //CIC_get_AMR(s, gp, x, width, boxsize);

  //double overlap = 0.;

  gp->potential_mesh = 0.;
  gp->a_grav_mesh[0] = 0.f;
  gp->a_grav_mesh[1] = 0.f;
  gp->a_grav_mesh[2] = 0.f;

  //overlap = get_overlap(pbox, home_cell->loc, home_cell->width);
  //if (overlap <= 0.) error("Overlap with home cell not valid");
  //message("Passing (%lf, %lf, %lf)", x[0], x[1], x[2]);
  int calc_acc = 0;
  double pot_part = 0.;
  CIC_get_AMR(s, gp, x, width, boxsize, stop_depth, calc_acc, &pot_part);
  calc_acc = 1;
  double acc_part[3] = {0.,0.,0.};
  CIC_get_AMR(s, gp, x, width, boxsize, stop_depth, calc_acc, acc_part);
  //message("Adding %lf to potential with value %lf", acc, gp->potential_mesh);
  gp->potential_mesh += pot_part;
  gp->a_grav_mesh[0] = (1./2.) * acc_part[0];
  gp->a_grav_mesh[1] = (1./2.) * acc_part[1];
  gp->a_grav_mesh[2] = (1./2.) * acc_part[2];

  /* 3-point stencil along each axis for the accelerations */
  /*x[0] = fmod(x[0] + width[0], boxsize);
  gp->a_grav_mesh[0] -= (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);
  x[0] = (x[0] - width[0] > 0) ? x[0] - width[0] : x[0] - width[0] + boxsize;
  gp->a_grav_mesh[0] += (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);

  x[1] = fmod(x[1] + width[1], boxsize);
  gp->a_grav_mesh[1] -= (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);
  x[1] = (x[1] - width[1] > 0) ? x[1] - width[1] : x[1] - width[1] + boxsize;
  gp->a_grav_mesh[1] += (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);

  x[2] = fmod(x[2] + width[2], boxsize);
  gp->a_grav_mesh[2] -= (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);
  x[2] = (x[2] - width[2] > 0) ? x[2] - width[2] : x[2] - width[2] + boxsize;
  gp->a_grav_mesh[2] += (1./2.) * CIC_get_AMR(s, gp, x, width, boxsize);*/

  //message("We added (%lf, %lf, %lf)", gp->a_grav_mesh[0], gp->a_grav_mesh[1], gp->a_grav_mesh[2]);

  /* Don't forget G */
  double const_G = s->e->physical_constants->const_newton_G;
  //message("Before, we had %lf", gp->potential_mesh);
  gp->a_grav_mesh[0] *= const_G;
  gp->a_grav_mesh[1] *= const_G;
  gp->a_grav_mesh[2] *= const_G;
  gp->potential_mesh *= const_G;
  //message("After, we have %lf", gp->potential_mesh);
}

void CIC_get_AMR(struct space *s, struct gpart *gp, double x[3], double width[3], double boxsize, int stop_depth, int calc_acc, double *acc) {
  double pbox[6] = {x[0], x[0] + width[0], x[1], x[1] + width[1], x[2], x[2] + width[2]};
  double pbox_shift[6];

  double x_shift[2] = {-boxsize, 0.};
  double y_shift[2] = {-boxsize, 0.};
  double z_shift[2] = {-boxsize, 0.};

  for (int i=0; i<2;i++) {
    for (int j=0; j<2;j++) {
      for (int k=0; k<2;k++) {
        pbox_shift[0] = pbox[0]+x_shift[i];
        pbox_shift[1] = pbox[1]+x_shift[i];
        pbox_shift[2] = pbox[2]+y_shift[j];
        pbox_shift[3] = pbox[3]+y_shift[j];
        pbox_shift[4] = pbox[4]+z_shift[k];
        pbox_shift[5] = pbox[5]+z_shift[k];
        if (pbox_shift[1]<0. || pbox_shift[3] <0. || pbox_shift[5]<0.) continue;
        //if (calc_acc) search_tree(gp, pbox_shift, s->cells_top, width[0], s->nr_cells, acc, stop_depth, 1, 0);
        //else search_tree(gp, pbox_shift, s->cells_top, width[0], s->nr_cells, &acc, stop_depth, 1, 0);
        search_tree(gp, pbox_shift, s->cells_top, width[0], s->nr_cells, acc, stop_depth, 1, 0, calc_acc);
        //message("We have found the value %lf", acc);
      }
    }
  }
  //return acc;
}

int sort_cell(const void *a, const void *b) {
    const struct cell *A = *(const struct cell * const *)a;
    const struct cell *B = *(const struct cell * const *)b;

    int ixA = (int)(A->loc[0] / A->width[0]);
    int ixB = (int)(B->loc[0] / B->width[0]);
    if (ixA != ixB) return ixA - ixB;

    int iyA = (int)(A->loc[1] / A->width[1]);
    int iyB = (int)(B->loc[1] / B->width[1]);
    if (iyA != iyB) return iyA - iyB;

    int izA = (int)(A->loc[2] / A->width[2]);
    int izB = (int)(B->loc[2] / B->width[2]);
    return izA - izB;
}

/* For now only for one depth. Add more depths later. */
void perform_nonuniform_calculation(struct space *s, int min_depth, int max_depth, struct AMR_levels *levels, int max_gridsize, double box_size) {
  //int smooth_tree = 1;

  for (int i=0; i<max_depth+1; i++) {
    get_patch_density(s, &levels[i]);
  }

  /* Link the uniform levels */
  //for (int i=0; i<min_depth+1;i++) {
    //link_uniform_level(&levels[i]);
  //}
  for (int i=0; i<min_depth; i++) {
    get_patch_potential(s, &levels[i+1], &levels[i]);
  }

  /* Maybe write some code to get the potential for the uniform levels at this point and not earlier */
  //for (int i=0;i<1;i++) {
  for (int i=min_depth;i<max_depth;i++) { //PLACEHOLDER
    get_patch_potential(s, &levels[i+1], &levels[i]);
    if (i==max_depth-1) {
      message("We found the patch potential %lf", levels[max_depth].mean_potential);
      //sleep(10);
    }
    max_gridsize *= 2;
    int nr_cells = levels[i+1].cell_count;
    message("The number of cells is %d", nr_cells);
    for (int j=0;j<nr_cells;j++) {
      (*(levels[i+1].cells[j])).mask_value = 1.;
    }
    printf("\n");
    message("Setting patch guess for depth %d. min_depth = %d", i+1, min_depth);
    //for (int j=0; j<levels[i].cell_count; j++) {
      //if (levels[i].cells[j]->split) message("The potential of coarse cell %d is %lf", j, levels[i].cells[j]->CIC_potential);
    //}
    set_patch_guess(s, &levels[i], &levels[i+1], nr_cells, max_gridsize, box_size, min_depth);

    double msq = 0.;
    for (int j=0; j<levels[i+1].cell_count + levels[i+1].ghost_count; j++) {
      msq += fabs(levels[i+1].cells[j]->CIC_potential);
    }
    msq = msq/(levels[i+1].cell_count + levels[i+1].ghost_count);
    message("The initial msq potential of level %d after prolongating from level %d is %lf", i+1, i, msq);
    
    /* Extrapolate to all levels. */
    for (int j=0; j<i+1; j++) {
      //message("Extrapolating to level %d", i-j);
      extrapolate_mask_values(&levels[i-j]); //I think we can just do this using the tree. Before summing we must first overwrite the values that were stored before, as they were relevant for the calculation on the previous level
    }

    //for (int j=0; j<levels[0].cell_count; j++) {
      //if (levels[0].cells[j]->mask_value != 1) {
        //message("Cell %d had mask value %lf", j, levels[0].cells[j]->mask_value);
      //}
    //}
    /*
    message("The ghost count is %d", levels[i+1].ghost_count);
    double fac = s->dim[0]/levels[i+1].cdim;
    for (int j=0; j<levels[i+1].ghost_count; j++) {
      struct cell *ghost = levels[i+1].cells[levels[i+1].cell_count+j];
      int cell_i = (ghost->loc[0]+0.001)/fac;
      int cell_j = (ghost->loc[1]+0.001)/fac;
      int cell_k = (ghost->loc[2]+0.001)/fac;

      if (cell_i == 0 && cell_j == 0 && cell_k == 0) ghost->CIC_potential = -2498.85353095973;
      if (cell_i == 0 && cell_j == 0 && cell_k == 1) ghost->CIC_potential = 435.977450810689;
      if (cell_i == 0 && cell_j == 1 && cell_k == 0) ghost->CIC_potential = -4523.97096905245;
      if (cell_i == 0 && cell_j == 1 && cell_k == 1) ghost->CIC_potential = -1569.15188828691;
      if (cell_i == 1 && cell_j == 0 && cell_k == 0) ghost->CIC_potential = -1297.60186602418;
      if (cell_i == 1 && cell_j == 0 && cell_k == 1) ghost->CIC_potential = 1506.67695673741;
      if (cell_i == 1 && cell_j == 1 && cell_k == 0) ghost->CIC_potential = -1961.62197957554;
      if (cell_i == 1 && cell_j == 1 && cell_k == 1) ghost->CIC_potential = 449.393180675248;

    }*/

    //nr_cells = levels[1].cell_count;
    ///message("Exporting AMR density and potential data");
      //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/FFT_16_initguess.txt", "w");
      //for (int j = 0; j < nr_cells; j++) {
        //fprintf(file_swift_density2, "%lf %.15g %.15g %.15g \n", levels[1].cells[j]->CIC_potential, levels[1].cells[j]->loc[0], levels[1].cells[j]->loc[1], levels[1].cells[j]->loc[2]);
      //}
      //fprintf(file_swift_density2, "%lf\n", 70.);
      //fclose(file_swift_density2);
      //message("Exported the AMR uniform potential data.");
      //sleep(15);

    perform_multigrid_acceleration(s, min_depth, max_depth, levels, levels[i+1].depth);
    msq = 0.;
    for (int j=0; j<levels[i+1].cell_count + levels[i+1].ghost_count; j++) {
      msq += fabs(levels[i+1].cells[j]->CIC_potential);
    }
    msq = msq/(levels[i+1].cell_count + levels[i+1].ghost_count);
    message("The final msq potential of level %d is %lf", i+1, msq);

    //for (int j=0; j<levels[i+1].cell_count; j++) {
      //struct cell *c = levels[i+1].cells[j];
      //message("The potential of cell %d at (%lf, %lf, %lf) is %lf", j, c->loc[0], c->loc[1], c->loc[2], c->CIC_potential);
    //}
  }

  //int nr_cells = levels[1].cell_count;
  //message("Exporting AMR density and potential data");
  //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/AMR_16_nooct_prolongated_trilinear_converged.txt", "w");
  //for (int j = 0; j < nr_cells; j++) {
    //fprintf(file_swift_density2, "%lf %.15g %.15g %.15g \n", levels[1].cells[j]->CIC_potential, levels[1].cells[j]->loc[0], levels[1].cells[j]->loc[1], levels[1].cells[j]->loc[2]);
  //}
  //fclose(file_swift_density2);
  //message("Exported the final AMR uniform potential data.");
  //sleep(15);

  //int nr_cells = levels[2].cell_count;
  //double fac = box_size/levels[1].cdim;
  //message("Exporting AMR density and potential data");
  //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/uniform_check/AMR_32_uniform.txt", "w");
  //for (int j = 0; j < nr_cells; j++) {
    //if (!levels[2].cells[j]->ghost) {
      //double search_loc_x = levels[2].cells[j]->loc[0];
      //double search_loc_y = levels[2].cells[j]->loc[1];
      //double search_loc_z = levels[2].cells[j]->loc[2];
      //int cell_i = (int) ((search_loc_x+0.01)/fac);
      //int cell_j = (int) ((search_loc_y+0.01)/fac);
      //int cell_k = (int) ((search_loc_z+0.01)/fac);

      //int cdim[3] = {levels[2].cdim, levels[2].cdim, levels[2].cdim};
      //int loc_1d = cell_getid(cdim, cell_i, cell_j, cell_k);

      //fprintf(file_swift_density, "cell %d at (%lf , %lf , %lf) has %lf and %lf \n", j, fine->cells[j]->loc[0], fine->cells[j]->loc[1], fine->cells[j]->loc[2], fine->cells[j]->CIC_density, fine->cells[j]->CIC_potential);
      //fprintf(file_swift_density2, "%lf\n", levels[2].cells[j]->CIC_potential);
      //fprintf(file_swift_density2, "%d\n", loc_1d);
    //}
  //}
  //fclose(file_swift_density2);
  //message("Exported the AMR potential data.");
  //sleep(15);

  int export = 0;
  if (export) {
    for (int i=0; i<max_depth-1; i++) {
      struct AMR_levels level = levels[max_depth-i-1];
      message("Summing to level %d", level.depth);
      for (int j=0; j<level.cell_count; j++) {
        if (level.cells[j]->split) {
          level.cells[j]->CIC_potential = 0.;
          for (int k=0; k<8; k++) {
            level.cells[j]->CIC_potential += level.cells[j]->progeny[k]->CIC_potential/8.;
          }
        }
      }
    }

    /* Sum to get the potential at level 1 if a cell is split */
    //for (int i=0; i<levels[1].cell_count; i++) {
      //if (levels[1].cells[i]->split) {
        //levels[1].cells[i]->CIC_potential = 0.;
        //for (int j=0; j<8; j++) {
          //levels[1].cells[i]->CIC_potential += levels[1].cells[i]->progeny[j]->CIC_potential/8.;
        //}
      //}
    //}
    //message("Exporting AMR potential data");
    //FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/AMR_potential_test/converged_alllevels_new_interpolation", "w");
    //for (int j = 0; j < levels[1].cell_count; j++) {
      //if (levels[1].cells[j]->split) {
        //fprintf(file_swift_density, "%lf\n", levels[1].cells[j]->CIC_potential);
      //}
      //else {
        //fprintf(file_swift_density, "%lf\n", 30000.);
      //}
    //}
    //fclose(file_swift_density);
    //message("Exported the AMR potential data.");
    //sleep(5);
  }
}

void modify_tree(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int new_cell_count[max_depth+1]) {
  for (int i=0; i<max_depth+1; i++) {
    new_cell_count[i] = 0;
  }
  for (int i=0; i<max_depth+1; i++) {
    int nr_cells = levels[i].cell_count;
    for (int j=0; j<nr_cells; j++) {
      struct cell *c = levels[i].cells[j];
      if (!(c->split) && c->refine) {
        //if (!c->refine2) {
          //message("Doing cell %d on level %d", j, i);
          //error("Not marked for refinement in original method");
        //}
        space_getcells(s, 8, c->progeny, 0, 0, 0);
        c->split = 1;
        for (int k=0; k<8; k++) {
          struct cell *child = c->progeny[k];
          child->loc[0] = c->loc[0];
          child->loc[1] = c->loc[1];
          child->loc[2] = c->loc[2];

          double width = c->width[0]/2;
          child->width[0] = width;
          child->width[1] = width;
          child->width[2] = width;

          if (k & 4) {
            child->loc[0] += child->width[0];
            c->progeny[k]->neighbours[1] = c->progeny[k-4];
            c->progeny[k-4]->neighbours[0] = c->progeny[k];
          }
          if (k & 2) {
            child->loc[1] += child->width[1];
            c->progeny[k]->neighbours[3] = c->progeny[k-2];
            c->progeny[k-2]->neighbours[2] = c->progeny[k];
          }
          if (k & 1) {
            child->loc[2] += child->width[2];
            c->progeny[k]->neighbours[5] = c->progeny[k-1];
            c->progeny[k-1]->neighbours[4] = c->progeny[k];
          }

          child->parent = c;
        }
        /* Put the particles in the cell */
        if (c->hydro.count != 0 || c->grav.count != 0 || c->stars.count != 0 ||
        c->black_holes.count != 0 || c->sinks.count != 0) {
          divide_particles(c, s);
        }

        /* Also add them to the level */
        levels[i+1].cells = realloc(levels[i+1].cells, ((levels[i+1].cell_count)+8)*sizeof(struct cell *));
        for (int k=0;k<8;k++) {
          levels[i+1].cells[(levels[i+1].cell_count)+k] = c->progeny[k];
        }
        /* Also update cell count for the level */
        levels[i+1].cell_count += 8;
        new_cell_count[i+1] += 8;
        //message("Cell count at level %d was increased and is now %d", levels[i+1].depth,new_cell_count[i+1]);
      }
    }
  }
}

void divide_particles(struct cell *c, struct space *s) {
  //message("Dividing particles");
  /* Allocate buffer, or something like that.. */
  struct cell_buff *restrict buff = NULL;
  struct cell_buff *restrict sbuff = NULL;
  struct cell_buff *restrict bbuff = NULL;
  struct cell_buff *restrict gbuff = NULL;
  struct cell_buff *restrict sink_buff = NULL;

  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;

  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct sink *sinks = c->sinks.parts;

  if (count > 0) {
    if (swift_memalign("tempbuff", (void **)&buff, SWIFT_STRUCT_ALIGNMENT,
                        sizeof(struct cell_buff) * count) != 0)
      error("Failed to allocate temporary indices.");
    for (int k = 0; k < count; k++) {
      buff[k].x[0] = parts[k].x[0];
      buff[k].x[1] = parts[k].x[1];
      buff[k].x[2] = parts[k].x[2];
    }
  }
  if (gcount > 0) {
    if (swift_memalign("tempgbuff", (void **)&gbuff, SWIFT_STRUCT_ALIGNMENT,
                        sizeof(struct cell_buff) * gcount) != 0)
      error("Failed to allocate temporary indices.");
    for (int k = 0; k < gcount; k++) {
      gbuff[k].x[0] = gparts[k].x[0];
      gbuff[k].x[1] = gparts[k].x[1];
      gbuff[k].x[2] = gparts[k].x[2];
    }
  }
  if (scount > 0) {
    if (swift_memalign("tempsbuff", (void **)&sbuff, SWIFT_STRUCT_ALIGNMENT,
                        sizeof(struct cell_buff) * scount) != 0)
      error("Failed to allocate temporary indices.");
    for (int k = 0; k < scount; k++) {
      sbuff[k].x[0] = sparts[k].x[0];
      sbuff[k].x[1] = sparts[k].x[1];
      sbuff[k].x[2] = sparts[k].x[2];
    }
  }
  if (bcount > 0) {
    if (swift_memalign("tempbbuff", (void **)&bbuff, SWIFT_STRUCT_ALIGNMENT,
                        sizeof(struct cell_buff) * bcount) != 0)
      error("Failed to allocate temporary indices.");
    for (int k = 0; k < bcount; k++) {
      bbuff[k].x[0] = bparts[k].x[0];
      bbuff[k].x[1] = bparts[k].x[1];
      bbuff[k].x[2] = bparts[k].x[2];
    }
  }
  if (sink_count > 0) {
    if (swift_memalign("temp_sink_buff", (void **)&sink_buff,
                        SWIFT_STRUCT_ALIGNMENT,
                        sizeof(struct cell_buff) * sink_count) != 0)
      error("Failed to allocate temporary indices.");
    for (int k = 0; k < sink_count; k++) {
      sink_buff[k].x[0] = sinks[k].x[0];
      sink_buff[k].x[1] = sinks[k].x[1];
      sink_buff[k].x[2] = sinks[k].x[2];
    }
  }

  cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts,
              c->black_holes.parts - s->bparts, c->sinks.parts - s->sinks,
              buff, sbuff, bbuff, gbuff, sink_buff);
}

/* Use the ghost slot to mark whether a cell should be refined. Ghost=1 means refine */
void build_refinement_map(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]) {
  for (int i=0; i<max_depth+1; i++) {
    if (levels[max_depth - i].depth < min_depth) continue;
    int nr_cells = levels[max_depth - i].cell_count;
    message("Accessing level %d with cell count %d", levels[max_depth - i].depth, nr_cells);
    for (int j=0; j<nr_cells; j++) {
      //if (levels[max_depth - i].depth==5) message("Checking cell %d", j);
      struct cell *curr_cell = levels[max_depth - i].cells[j];
      if (!curr_cell->split) continue;
        //if (levels[max_depth - i].depth==5) message("Cell is split");
        /* Step 1 */
      for (int k=0;k<8;k++) {
        //if (max_depth-i == 2 && k == 0) message("Checking child %d", k);
        if (curr_cell->progeny[k] == NULL) message("Child %d somehow doesn't exist", k);
        if (curr_cell->progeny[k]->split || curr_cell->progeny[k]->refine) {
          curr_cell->refine = 1;
          //curr_cell->refine2 = 1;
          /* Step 2 */
          //mark_neighbours(s, min_depth, &levels[max_depth - i], curr_cell);
          mark_all_neighbours(s, min_depth, &levels[max_depth -i], curr_cell);
          break;
        }
      }
    }
    int to_refine = 0;
    for (int j=0; j<nr_cells; j++) {
      if (levels[max_depth - i].cells[j]->refine) {
        //message("Marked the cell at (%lf, %lf, %lf) for refinement", levels[max_depth - i].cells[j]->loc[0], levels[max_depth - i].cells[j]->loc[1], levels[max_depth - i].cells[j]->loc[2]);
        to_refine += 1;
      }
    }
    message("At level %d we just marked %d cells for refinement", levels[max_depth - i].depth, to_refine);
  }
}

enum {
    XP = 0,  // x+1
    XM = 1,  // x-1
    YP = 2,  // y+1
    YM = 3,  // y-1
    ZP = 4,  // z+1
    ZM = 5   // z-1
};

void mark_all_neighbours(struct space *s, int min_depth, struct AMR_levels *level, struct cell *curr_cell) {
  //message("Marking neighbours for cell at (%lf, %lf, %lf)", curr_cell->loc[0], curr_cell->loc[1], curr_cell->loc[2]);
  int counter = 0;
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dz = -1; dz <= 1; dz++) {

        struct cell *c = curr_cell;
          if (dx == 0 && dy == 0 && dz == 0) continue;

          if (dx == 1 && c) c = c->neighbours[XP];
          if (dx == -1 && c) c = c->neighbours[XM];

          if (dy == 1 && c) c = c->neighbours[YP];
          if (dy == -1 && c) c = c->neighbours[YM];

          if (dz == 1 && c) c = c->neighbours[ZP];
          if (dz == -1 && c) c = c->neighbours[ZM];

          if (c) {
            if (!c->refine) counter += 1;
            c->refine = 1;
            //message("Marked cell at (%lf, %lf, %lf)", c->loc[0], c->loc[1], c->loc[2]);
          }
          if (counter == 26) break;
      }
    }
  }
}

void mark_neighbours(struct space *s, int min_depth, struct AMR_levels *level, struct cell *curr_cell) {
  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Box wrap the multipole's position */
  const double search_loc_x = box_wrap(curr_cell->loc[0], 0., dim[0]);
  const double search_loc_y = box_wrap(curr_cell->loc[1], 0., dim[1]);
  const double search_loc_z = box_wrap(curr_cell->loc[2], 0., dim[2]);

  int gridsize = level->cdim;
  double box_size = s->dim[0];
  double fac = box_size/gridsize;

  if (level->depth <= min_depth) {
    /* Mark the neighbours using the cellIDs */
    int cdim[3] = {gridsize, gridsize, gridsize};

    int cell_i = (int) ((search_loc_x+0.0001)/fac);
    int i_plus = (cell_i+1)%gridsize;
    int i_min = (cell_i-1>=0) ? (cell_i-1) % cdim[0] : (cell_i-1) % cdim[0] + cdim[0];
    int cell_j = (int) ((search_loc_y+0.0001)/fac);
    int j_plus = (cell_j+1)%gridsize;
    int j_min = (cell_j-1>=0) ? (cell_j-1) % cdim[0] : (cell_j-1) % cdim[0] + cdim[0];
    int cell_k = (int) ((search_loc_z+0.0001)/fac);
    int k_plus = (cell_k+1)%gridsize;
    int k_min = (cell_k-1>=0) ? (cell_k-1) % cdim[0] : (cell_k-1) % cdim[0] + cdim[0];

    //if (level->depth == 2) {
      //message("x-loc = %lf and got set to %d and then imin to %d", search_loc_x, cell_i, i_min);
      //message("y-loc = %lf and got set to %d", search_loc_y, cell_j);
      //message("z-loc = %lf and got set to %d", search_loc_z, cell_k);
      //message("Searching for %d, %d, %d, %d, %d, %d", cell_getid(cdim, i_plus, cell_j, cell_k), cell_getid(cdim, i_min, cell_j, cell_k), cell_getid(cdim, cell_i, j_plus, cell_k), cell_getid(cdim, cell_i, j_min, cell_k), cell_getid(cdim, cell_i, cell_j, k_plus), cell_getid(cdim, cell_i, cell_j, k_min));
    //}

    level->cells[cell_getid(cdim, i_plus, cell_j, cell_k)]->refine = 1;
    level->cells[cell_getid(cdim, i_min, cell_j, cell_k)]->refine = 1;
    level->cells[cell_getid(cdim, cell_i, j_plus, cell_k)]->refine = 1;
    level->cells[cell_getid(cdim, cell_i, j_min, cell_k)]->refine = 1;
    level->cells[cell_getid(cdim, cell_i, cell_j, k_plus)]->refine = 1;
    level->cells[cell_getid(cdim, cell_i, cell_j, k_min)]->refine = 1;

    /* Create extra cells for trilinear interpolation */
    level->cells[cell_getid(cdim, i_plus, j_plus, cell_k)]->refine = 1;
    level->cells[cell_getid(cdim, i_plus, cell_j, k_plus)]->refine = 1;
    level->cells[cell_getid(cdim, cell_i, j_plus, k_plus)]->refine = 1;
    level->cells[cell_getid(cdim, i_plus, j_plus, k_plus)]->refine = 1;

  }

  else {
    for (int i=0; i<6; i++) {
      if (curr_cell->neighbours[i] != NULL) curr_cell->neighbours[i]->refine = 1;
    }

    if (curr_cell->neighbours[0] != NULL) { 
      if (curr_cell->neighbours[0]->neighbours[2] != NULL) { 
        curr_cell->neighbours[0]->neighbours[2]->refine = 1; 
        if (curr_cell->neighbours[0]->neighbours[2]->neighbours[4] != NULL) curr_cell->neighbours[0]->neighbours[2]->neighbours[4]->refine = 1; 
      } 
      if (curr_cell->neighbours[0]->neighbours[4] != NULL) { 
        curr_cell->neighbours[0]->neighbours[4]->refine = 1; 
        if (curr_cell->neighbours[0]->neighbours[4]->neighbours[2] != NULL) curr_cell->neighbours[0]->neighbours[4]->neighbours[2]->refine = 1; 
      } 
    } 
    if (curr_cell->neighbours[2] != NULL) { 
      if (curr_cell->neighbours[2]->neighbours[0] != NULL) { 
        curr_cell->neighbours[2]->neighbours[0]->refine = 1; 
        if (curr_cell->neighbours[2]->neighbours[0]->neighbours[4] != NULL) curr_cell->neighbours[2]->neighbours[0]->neighbours[4]->refine = 1; 
      } 
      if (curr_cell->neighbours[2]->neighbours[4] != NULL) { 
        curr_cell->neighbours[2]->neighbours[4]->refine = 1; 
        if (curr_cell->neighbours[2]->neighbours[4]->neighbours[0] != NULL) curr_cell->neighbours[2]->neighbours[4]->neighbours[0]->refine = 1; 
      } 
    } 
    if (curr_cell->neighbours[4] != NULL) { 
      if (curr_cell->neighbours[4]->neighbours[0] != NULL) { 
        curr_cell->neighbours[4]->neighbours[0]->refine = 1; 
        if (curr_cell->neighbours[4]->neighbours[0]->neighbours[2] != NULL) curr_cell->neighbours[4]->neighbours[0]->neighbours[2]->refine = 1; 
      } 
      if (curr_cell->neighbours[4]->neighbours[2] != NULL) { 
        curr_cell->neighbours[4]->neighbours[2]->refine = 1; 
        if (curr_cell->neighbours[4]->neighbours[2]->neighbours[0] != NULL) curr_cell->neighbours[4]->neighbours[2]->neighbours[0]->refine = 1; 
      } 
    } 
  }
}

int search_neighbours(struct AMR_levels *level, double search_loc[3], double fac) {
  int neighbour = 0;
  //message("Doing nb search");
  /* Do a neighbour search */

  //if (search_loc[0] <2. && search_loc[1]>22. && search_loc[1]<23. && search_loc[2] >40. && search_loc[2] <41.) {
    //message("Searching for (%lf,%lf,%lf)", search_loc[0], search_loc[1], search_loc[2]);
  //}
  //int nr_overlap = 0;
  for (int i=0; i<level->cell_count; i++) {
    double x_loc = level->cells[i]->loc[0];
    double y_loc = level->cells[i]->loc[1];
    double z_loc = level->cells[i]->loc[2];

    double cell_check[6] = {x_loc, x_loc + fac, y_loc, y_loc + fac, z_loc, z_loc + fac};
    double overlap = get_overlap(cell_check, search_loc, level->cells[0]->width);
    //if (cell_check[1] > 2. && cell_check[1] < 3. && cell_check[2]>22. && cell_check[2]<23. && cell_check[4] >40. && cell_check[4] <41.) {
      //message("Comparing with (%lf,%lf,%lf)", cell_check[0], cell_check[2], cell_check[4]);
      //message("The overlap is %lf", overlap);
    //}

    if (overlap > 0.1) { //Set refine flag
      level->cells[i]->refine = 1;
      //nr_overlap++;
      neighbour = 1;
      break;
    }
  }
  //if (nr_overlap>1) message("Found %d overlapping cells", nr_overlap);
  return neighbour;
}

void link_uniform_level(struct AMR_levels *level) {
  int cdim[3] = {level->cdim, level->cdim, level->cdim};
  for (int k=0; k<cdim[2]; k++) {
    int k_plus = (k+1) % cdim[2];
    int k_min = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
    for (int j=0; j<cdim[1]; j++) {
      int j_plus = (j+1) % cdim[1];
      int j_min = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
      for (int i=0; i<cdim[0]; i++) {
        int i_plus = (i+1) % cdim[0];
        int i_min = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
        size_t cid = cell_getid(cdim,i,j,k);

        /* Linking one way*/
        level->cells[cid]->neighbours[0] = level->cells[cell_getid(cdim,i_plus,j,k)];
        level->cells[cid]->neighbours[1] = level->cells[cell_getid(cdim,i_min,j,k)];
        level->cells[cid]->neighbours[2] = level->cells[cell_getid(cdim,i,j_plus,k)];
        level->cells[cid]->neighbours[3] = level->cells[cell_getid(cdim,i,j_min,k)];
        level->cells[cid]->neighbours[4] = level->cells[cell_getid(cdim,i,j,k_plus)];
        level->cells[cid]->neighbours[5] = level->cells[cell_getid(cdim,i,j,k_min)];

        /* Linking the other way */
        level->cells[cell_getid(cdim,i_plus,j,k)]->neighbours[1] = level->cells[cid];
        level->cells[cell_getid(cdim,i_min,j,k)]->neighbours[0] = level->cells[cid];
        level->cells[cell_getid(cdim,i,j_plus,k)]->neighbours[3] = level->cells[cid];
        level->cells[cell_getid(cdim,i,j_min,k)]->neighbours[2] = level->cells[cid];
        level->cells[cell_getid(cdim,i,j,k_plus)]->neighbours[5] = level->cells[cid];
        level->cells[cell_getid(cdim,i,j,k_min)]->neighbours[4] = level->cells[cid];
      }
    }
  }
}

void get_patch_potential(struct space *s, struct AMR_levels *fine, struct AMR_levels *coarse) {
  fine->mean_potential = 0.;
  int coarse_split_count = fine->cell_count/8;
  for (int i=0; i<coarse->cell_count; i++) {
    struct cell *c = coarse->cells[i];
    if (!c->split) continue;
    //if (fine->depth == 3) message("Adding %lf", c->CIC_potential);
    fine->mean_potential += c->CIC_potential;
  }
  fine->mean_potential /= coarse_split_count;
  //if (fine->depth == 3) sleep(15);
}

void get_patch_density(struct space *s, struct AMR_levels *level) {
  double mean_density = 0.0;
  int nr_cells = level->cell_count;
  int N = level->cdim;
  double box_size = s->dim[0];

  // Find the total mass of all particles
  double mass_tot = 0.;
  size_t level_cells = (size_t) N*N*N;
  for (size_t i=0; i<s->nr_gparts; ++i) {
    mass_tot += s->gparts[i].mass;
  }
  mean_density = mass_tot/level_cells;
  //mean_density = sum/nr_cells;
  level->mean_density = mean_density;
  message("The mean density of level %d is %lf and the cell count was %d", level->depth, level->mean_density, nr_cells);
  for (int i = 0; i < nr_cells; i++) {
    if (mean_density >0.) level->cells[i]->CIC_density = (level->cells[i]->CIC_density)/mean_density;
  }

  //sum =0.;
  //for (int i = 0; i < nr_cells; i++) {
    //sum += level->cells[i]->CIC_density;
  //}
  //mean_density = sum/nr_cells;
  //level->mean_density *= 4.*M_PI*nr_cells/(box_size*box_size*box_size);
  double fac = N/box_size;
  level->mean_density *= 4.*M_PI*fac*fac*fac;
}

void perform_multigrid_acceleration(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int current_depth) {
  double tolerance = 10e-5;
  int V_cycles = 0;
  double delta = s->dim[0]/(levels[current_depth].cdim);
  message("The value of delta is %lf", delta);
  int V_max = 5;

  message("The mean density of the current level is %lf", levels[current_depth].mean_density);

  double residual = get_patch_residual(levels[current_depth], delta);
  message("The first residual is %lf", residual);

  while (residual>tolerance && V_cycles < V_max) {
    /* Pre-smoothing for 10 steps */
    for (int i=0; i<10; i++) {
      perform_patch_sweep(&levels[current_depth], delta);
      residual = get_patch_residual(levels[current_depth], delta);
      double msq = 0;
      for (int j=0; j<levels[current_depth].cell_count; j++) {
        msq += fabs(levels[current_depth].cells[j]->CIC_potential);
      }
      msq = msq/(levels[current_depth].cell_count);
      //message("The msq potential of level %d after %d step is %lf", current_depth, i, msq);
      //message("Performed %d pre-smoothing step", i);
      //sleep(20);
    }
    //sleep(10);

    //int nr_cells = levels[1].cell_count;
  //double fac = box_size/levels[1].cdim;
    //if (V_cycles == 5) {
      //message("Exporting AMR density and potential data");
      //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/uniform_check/AMR_16_uniform_middle_result.txt", "w");
      //for (int j = 0; j < nr_cells; j++) {
        //  fprintf(file_swift_density2, "%lf\n", levels[1].cells[j]->CIC_potential);
      //}
      //fclose(file_swift_density2);
      //message("Exported the AMR potential data.");
      //sleep(5);
    //}

    /* Initialise V-cycle recursion */
    transfer_residual_array(levels[current_depth], &levels[current_depth-1],delta);
    to_coarser_patch(levels, delta, current_depth-1, max_depth, current_depth);
    //residual = get_patch_residual(levels[current_depth], delta);

    /* Post-smoothing for 10 steps */
    for (int i=0; i<10; i++) {
      perform_patch_sweep(&levels[current_depth], delta);
    }
    residual = get_patch_residual(levels[current_depth], delta);
    V_cycles += 1;
    message("After %d V_cycle the residual is %lf", V_cycles, residual);
  }

  /* Some extra convergence if the V-cycles weren't enough */
  int counter_post_smoothing = 0;
  while (residual > tolerance) {
    counter_post_smoothing += 1;
    perform_patch_sweep(&levels[current_depth], delta);
    residual = get_patch_residual(levels[current_depth], delta);
    if (counter_post_smoothing % 100 == 0) message("Did %d steps in post-smoothing and the residual is %lf", counter_post_smoothing, residual);
  }
  message("Had to do post-smoothing on level %d for %d steps and the residual is %lf", current_depth, counter_post_smoothing, residual);

  int normalise = 1;
  if (normalise) {

    /* Find the mean potential on the patch */
    double pot_mean = 0.;

    for (int i=0; i<levels[current_depth].cell_count; i++) {
      pot_mean += levels[current_depth].cells[i]->CIC_potential;
    }
    pot_mean /= levels[current_depth].cell_count;

    message("The current mean density is %lf, while it is supposed to be %lf", pot_mean, levels[current_depth].mean_potential);

    /* Subtract the mean potential from the values on the finest level */
    for (int i=0; i<levels[current_depth].cell_count; i++) {
      levels[current_depth].cells[i]->CIC_potential -= (pot_mean - levels[current_depth].mean_potential);
    }
  }

  //Set mean of potential to zero if uniform grid
  if (levels[current_depth].cell_count == levels[current_depth].cdim * levels[current_depth].cdim * levels[current_depth].cdim) { //I.e. secretly uniform
    double sum = 0.;
    for (int i = 0; i < levels[current_depth].cell_count; i++) {
        sum += levels[current_depth].cells[i]->CIC_potential;
    }
    for (int i = 0; i<levels[current_depth].cell_count; i++) {
      levels[current_depth].cells[i]->CIC_potential -= sum/(levels[current_depth].cell_count); 
    }

    double potential_sum = 0.;
    for (int i =0; i<levels[current_depth].cell_count; i++)
      potential_sum += fabs(levels[current_depth].cells[i]->CIC_potential);
    message("The mean absolute value of the potential is %lf", potential_sum/(levels[current_depth].cell_count));
  }

  //for (int i=0; i<levels[current_depth].ghost_count; i++) {
    //struct cell *ghost = levels[current_depth].cells[levels[current_depth].cell_count+i];
    //message("After convergence ghost cell (%lf, %lf, %lf) has value %lf", ghost->loc[0], ghost->loc[1], ghost->loc[2], ghost->CIC_potential);
  //}
  //sleep(20);
 
  //int nr_cells = levels[1].cell_count;
  //message("Exporting AMR density and potential data");
      //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/AMR_16_nooct_ghosts_prolongated_init0.txt", "w");
      //for (int j = 0; j < nr_cells; j++) {
        //fprintf(file_swift_density2, "%lf %.15g %.15g %.15g \n", levels[1].cells[j]->CIC_potential, levels[1].cells[j]->loc[0], levels[1].cells[j]->loc[1], levels[1].cells[j]->loc[2]);
      //}
      //fprintf(file_swift_density2, "%lf\n", 70.);
      //fclose(file_swift_density2);
      //message("Exported the AMR uniform potential data.");
      //sleep(15);
}

void to_coarser_patch(struct AMR_levels *levels, double delta, int current_depth, int max_depth, int active_depth) {
  //message("Converging on level %d", current_depth);
  delta *= 2.0;
  double tolerance = 10e-4;
  double conv_tolerance = 0.9999;
  //double conv_tolerance2 = 0.1;
  int counter = 0;
  int counter_abs=0;
  //double msq = 0;
  //int nr_steps = 10;

  double residual_old = 0.;
  double residual_new = 1.;

  int nr_active = 0;
  for (int i=0;i<levels[current_depth].cell_count;i++) {
    if (levels[current_depth].cells[i]->mask_value > 0) {
      nr_active +=1;
    }
  }
  message("%d active cells for second-order", nr_active);
  if (nr_active == 0) return;

  double msq = 0;
  for (int j=0; j<levels[current_depth].cell_count; j++) {
    if (levels[current_depth].cells[j]->mask_value > 0) {
      msq += fabs(levels[current_depth].cells[j]->conv_residual);
    }
  }
  msq = msq/nr_active;
  message("The initial msq potential of level %d is %lf", current_depth, msq);

  msq = 0;
  for (int j=0; j<levels[current_depth].cell_count; j++) {
    if (levels[current_depth].cells[j]->mask_value > 0) {
      msq += levels[current_depth].cells[j]->CIC_density;
    }
  }
  msq = msq/nr_active;
  message("The initial mean density of level %d is %lf", current_depth, msq);
  if (current_depth == 0) { //On the coarsest grid, so solve the equation exactly 
    while (residual_new > tolerance) {
      counter +=1;
      counter_abs +=1;
      perform_masked_patch_sweep(&levels[current_depth], delta);
      residual_new = get_masked_residual(levels[current_depth], delta, nr_active, counter);
      //msq = get_solution_magnitude(levels[current_depth]);
      msq = 0;
      for (int j=0; j<levels[current_depth].cell_count; j++) {
        if (levels[current_depth].cells[j]->mask_value > 0) {
          msq += fabs(levels[current_depth].cells[j]->conv_residual);
        }
      }
      msq = msq/nr_active;
      //message("The msq potential of level %d after %d step is %lf", current_depth, counter, msq);
      //message("The msq of the current solution is %lf", sqrt(msq/nr_active));
      if (counter_abs>1) {
        //message("The previous residual was %lf and the new one is %lf", residual_old, residual_new);
        if (residual_new/residual_old>conv_tolerance && counter_abs >20) {
          //converge_first_order(&levels[current_depth], delta, nr_steps);
          //message("Finished with first order convergence");
          //counter =0;
          to_finer_patch(levels, current_depth+1, active_depth);
          //sleep(5);
          return;
          //error("Residuals diverging during exact solving at depth %d", current_depth);
        }
      }
      residual_old = residual_new;
    }  
    message("Going back to level %d", current_depth+1);
    to_finer_patch(levels, current_depth+1, active_depth);
  }

  else {
    //message("Calling this block for level %d", current_depth);
    for (int i=0;i<20;i++) {
      perform_masked_patch_sweep(&levels[current_depth], delta);
      residual_new = get_masked_residual(levels[current_depth], delta, nr_active, counter);
      msq = 0;
      for (int j=0; j<levels[current_depth].cell_count; j++) {
        if (levels[current_depth].cells[j]->mask_value > 0) {
          msq += fabs(levels[current_depth].cells[j]->conv_residual);
        }
      }
      msq = msq/nr_active;
      //message("The msq potential of level %d after %d step is %lf", current_depth, i, msq);
      //msq = get_solution_magnitude(levels[current_depth]);
      if (i>0) {
        //message("The previous residual was %lf and the new one is %lf", residual_old, residual_new);
        if (residual_new/residual_old>conv_tolerance && counter_abs>20) {
          //converge_first_order(&levels[current_depth], delta, nr_steps);
          //message("Finished with first order convergence");
          break;
        }
      }
      residual_old = residual_new;
    }

    //sleep(10);
    transfer_coarser_residual(levels[current_depth], &levels[current_depth-1], delta);
    to_coarser_patch(levels, delta, current_depth-1, max_depth, active_depth);

    /* Post-smoothing */
    message("Post-smoothing on level %d", current_depth);
    for (int i=0; i<20; i++) {
      perform_masked_patch_sweep(&levels[current_depth], delta);
      residual_new = get_masked_residual(levels[current_depth], delta, nr_active, counter);
      //msq = get_solution_magnitude(levels[current_depth]);
      if (i>0) {
        //message("The previous residual was %lf and the new one is %lf", residual_old, residual_new);
        if (residual_new/residual_old>conv_tolerance && counter_abs>20) {
          //converge_first_order(&levels[current_depth], delta, nr_steps);
          //message("Did first order convergence");
          //message("The previous residual was %lf and the new one is %lf", residual_old, residual_new);
          //error("Residuals diverging during post-smoothing at depth %d", current_depth);
          //message("Going back to level %d", current_depth+1);
          if (current_depth<max_depth) to_finer_patch(levels, current_depth+1, active_depth);
          //sleep(5);
          return;
        }
      }
      residual_old = residual_new;
    }
    message("Going back to level %d", current_depth+1);
    if (current_depth<max_depth) to_finer_patch(levels, current_depth+1, active_depth);
    //message("Survived");
  }

  //double msq_pot = 0;
  //for (int j=0; j<levels[current_depth].cell_count; j++) {
    //if (levels[current_depth].cells[j]->mask_value > 0) msq_pot += fabs(levels[current_depth].cells[j]->CIC_potential);
  //}
  //msq_pot = msq_pot/(nr_active);
  //message("The msq potential of level %d after %d steps is %lf", current_depth, counter, msq_pot);

  //double mean_dens = 0;
  //for (int j=0; j<levels[current_depth].cell_count; j++) {
    //if (levels[current_depth].cells[j]->mask_value > 0) mean_dens +=levels[current_depth].cells[j]->CIC_density;
  //}
  //mean_dens = mean_dens/(nr_active);
  //message("The mean density/residual of level %d after %d steps is %lf", current_depth, counter, mean_dens);
 
}

void to_finer_patch(struct AMR_levels *levels, int target_depth, int active_depth) {
  int fake=0;
  for (int i=0; i<levels[target_depth].cell_count;i++) {
    //if (target_depth==2) message("the cell count is %d", levels[target_depth].cell_count);
    struct cell *current_cell = levels[target_depth].cells[i];
    
    if (current_cell->mask_value<=0) continue;
    if (levels[target_depth].cells[i]->parent == NULL) {
      if (current_cell->hydro.count == 0 && current_cell->grav.count == 0 && current_cell->stars.count == 0 &&
                current_cell->black_holes.count == 0 && current_cell->sinks.count == 0) fake=1;
      message("Fake value is %d", fake);
    }
    if (current_cell->parent->mask_value<=0) continue;
    if (target_depth<active_depth) levels[target_depth].cells[i]->conv_residual += current_cell->parent->conv_residual;
    else levels[target_depth].cells[i]->CIC_potential += current_cell->parent->conv_residual;
  }
}

void converge_first_order(struct AMR_levels *level, double delta, int nr_steps) {
  double residual;
  for (int i=0; i<nr_steps+100; i++) {
    first_order_sweep(level, delta);
    residual = get_first_order_residual(level, delta);
  }
  message("The residual after %d 1st order steps is %lf", nr_steps+100,residual);
}

double get_first_order_residual(struct AMR_levels *level, double delta) {
  double pot[6];
  double residual = 0.;
  double mask_nb;
  int nr_active=0;

  for (int i=0;i<level->cell_count;i++) {
    struct cell *current_cell = level->cells[i];
    if (current_cell->mask_value<1) continue;
    nr_active +=1;
    for (int j=0;j<6;j++) {
      mask_nb = current_cell->neighbours[j]->mask_value;
      if (mask_nb<1) {
        pot[j] = 0.;
      }
      else pot[j] = current_cell->neighbours[j]->CIC_potential;
    }
    double res = (pot[0] + pot[1] + pot[2] + pot[3] + pot[4] + pot[5] - 6.0*current_cell->CIC_potential)/(delta*delta) - current_cell->CIC_density;
    residual += res*res;
  }
  //message("Found %d active cells doing first order convergence", nr_active);
  return sqrt(residual/(nr_active));

}

void first_order_sweep(struct AMR_levels *level, double delta) {
  double pot[6];
  double mask_nb;

  for (int col=0; col<2; col++){
    for (int i=0; i<level->cell_count;i++) {
      struct cell *current_cell = level->cells[i];
      if (current_cell->mask_value<1) {
        //message("The mask is %lf, so skipped this", current_cell->mask_value);
        continue;
      }
      int location =  (int) ((current_cell->loc[0]+0.01)/delta) + (int) ((current_cell->loc[1]+0.01)/delta) + (int) ((current_cell->loc[2]+0.01)/delta);
      //message("Checking location  (%lf,%lf,%lf) with index (%d,%d,%d)", current_cell->loc[0], current_cell->loc[1], current_cell->loc[2], (int) ((current_cell->loc[0]+0.01)/delta), (int) ((current_cell->loc[1]+0.01)/delta), (int) ((current_cell->loc[2]+0.01)/delta));
      if (location%2 != col) {
        //message("Skipped this cell on location grounds");
        continue;
      }
      for (int j=0;j<6;j++) {
        mask_nb = current_cell->neighbours[j]->mask_value;
        //message("The neighbour has mask %lf", mask_nb);
        if (mask_nb<1) {
          pot[j] = 0.;
          //message("Added lin extrapolated potential %lf", pot[j]);
        }
        else {
          pot[j] = current_cell->neighbours[j]->CIC_potential;
          //message("Added neighbour potential %lf", pot[j]);
        }
      }
      //message("The density is %lf", current_cell->CIC_density);
      level->cells[i]->CIC_potential = ((pot[0] + pot[1] + pot[2] + pot[3] + pot[4] + pot[5] - delta*delta*current_cell->CIC_density)/6.0);
      //if (fabs(level->cells[i]->CIC_potential)>10000.) message("Added %lf", level->cells[i]->CIC_potential);
    }
  }
}

double get_solution_magnitude(struct AMR_levels level) {
  double msq=0;
  for (int i=0;i<level.cell_count;i++) {
    struct cell *current_cell = level.cells[i];
    if (current_cell->mask_value<=0) continue;
    msq += current_cell->CIC_potential*current_cell->CIC_potential;
  }
  return msq;
}

double get_masked_residual(struct AMR_levels level, double delta, int nr_active, int step_nr) {
  int delta_nb[6];
  double pot[6];
  double residual = 0.;
  double mask_nb[6];

  for (int i=0;i<level.cell_count;i++) {
    struct cell *current_cell = level.cells[i];
    if (current_cell->mask_value<=0) continue;
    for (int j=0;j<6;j++) {
      mask_nb[j] = current_cell->neighbours[j]->mask_value;
      if (mask_nb[j]<=0) {
        delta_nb[j] = 0;
      }
      else {
        delta_nb[j] = 1;
        pot[j] = current_cell->neighbours[j]->conv_residual;
      }
    }
    double potential_term = ((double) delta_nb[0]*pot[0] + (double) delta_nb[1]*pot[1] + (double) delta_nb[2]*pot[2] + (double) delta_nb[3]*pot[3] + 
                                (double) delta_nb[4]*pot[4] + (double) delta_nb[5]*pot[5]);
    double ghost_correction = (((((double) 1-delta_nb[0])*mask_nb[0] + ((double) 1-delta_nb[1])*mask_nb[1] + ((double) 1-delta_nb[2])*mask_nb[2] +
                                  ((double) 1-delta_nb[3])*mask_nb[3] + ((double) 1-delta_nb[4])*mask_nb[4] + ((double) 1-delta_nb[5])*mask_nb[5]))/(current_cell->mask_value));
    double res = (potential_term - (6.-ghost_correction)*current_cell->conv_residual)/(delta*delta) - current_cell->CIC_density;
    //double res = (pot[0] + pot[1] + pot[2] + pot[3] + pot[4] + pot[5] - 6.0*current_cell->CIC_potential)/(delta*delta) - current_cell->CIC_density;
    //if (level.depth == 1 && i<30) {
      //message("We just found the residual %lf", res);
    //}
    residual += res*res;
  }
  return sqrt(residual/(nr_active));
}

void perform_masked_patch_sweep(struct AMR_levels *level, double delta) {
  /* Delta is set to 1 if we have a neighbour, and to 0 if a ghost*/
  int delta_nb[6];
  double mask_nb[6];
  double pot[6] = {0,0,0,0,0,0};

  for (int col=0; col<2; col++){
    for (int i=0; i<level->cell_count;i++) {
      struct cell *current_cell = level->cells[i];
      if (current_cell->mask_value<=0) {
        //message("The mask is %lf, so skipped this", current_cell->mask_value);
        continue;
      }
      int location =  (int) ((current_cell->loc[0]+0.0001)/delta) + (int) ((current_cell->loc[1]+0.0001)/delta) + (int) ((current_cell->loc[2]+0.0001)/delta);
      //message("Checking location  (%lf,%lf,%lf) with index (%d,%d,%d)", current_cell->loc[0], current_cell->loc[1], current_cell->loc[2], (int) ((current_cell->loc[0]+0.01)/delta), (int) ((current_cell->loc[1]+0.01)/delta), (int) ((current_cell->loc[2]+0.01)/delta));
      if (location%2 != col) {
        //message("Skipped this cell on location grounds");
        continue;
      }
      for (int j=0;j<6;j++) {
        mask_nb[j] = current_cell->neighbours[j]->mask_value;
        //message("The neighbour has mask %lf", mask_nb);
        if (mask_nb[j]<=0) {
          delta_nb[j] = 0;
          //message("Added lin extrapolated potential %lf", pot[j]);
        }
        else {
          delta_nb[j] = 1;
          pot[j] = current_cell->neighbours[j]->CIC_potential;
          //message("Added neighbour potential %lf", pot[j]);
        }
      }
      //message("The density is %lf", current_cell->CIC_density);
      //if (level->depth == 0) { 
        //message("The density is %lf", current_cell->CIC_density);
        //message("The potentials are (%lf, %lf, %lf, %lf, %lf, %lf)", (double) delta_nb[0]*pot[0], (double) delta_nb[1]*pot[1], (double) delta_nb[2]*pot[2], (double) delta_nb[3]*pot[3], (double) delta_nb[4]*pot[4], (double) delta_nb[5]*pot[5]);
        //message("The old value of the cell is %lf", level->cells[i]->CIC_potential); 
        //message("The delta values are (%d,%d,%d,%d,%d,%d)", delta_nb[0], delta_nb[1], delta_nb[2], delta_nb[3], delta_nb[4], delta_nb[5]);
      //}
      double ghost_correction = (6. - (((double) 1-delta_nb[0])*mask_nb[0] + ((double) 1-delta_nb[1])*mask_nb[1] + ((double) 1-delta_nb[2])*mask_nb[2] +
                                   ((double) 1-delta_nb[3])*mask_nb[3] + ((double) 1-delta_nb[4])*mask_nb[4] + ((double) 1-delta_nb[5])*mask_nb[5])/(current_cell->mask_value));
      double potential_term = ((double) delta_nb[0]*pot[0] + (double) delta_nb[1]*pot[1] + (double) delta_nb[2]*pot[2] + (double) delta_nb[3]*pot[3] + 
                                (double) delta_nb[4]*pot[4] + (double) delta_nb[5]*pot[5] - delta*delta*current_cell->CIC_density);
      level->cells[i]->conv_residual = potential_term/ghost_correction;

      //if (level->depth == 0) {
        //message("The potential term is %lf and the ghost correction %lf", potential_term, ghost_correction);
        //message("The new value of the cell is %lf", level->cells[i]->CIC_potential);
      //}
      //level->cells[i]->CIC_potential = ((pot[0] + pot[1] + pot[2] + pot[3] + pot[4] + pot[5] - delta*delta*current_cell->CIC_density)/6.0);
      //if (fabs(level->cells[i]->CIC_potential)>10000.) message("Added %lf", level->cells[i]->CIC_potential);
    }
  }
}

void transfer_coarser_residual(struct AMR_levels fine, struct AMR_levels *coarse, double delta) {
  message("Transferring the residual to level %d", coarse->depth);
  double sum = 0.;
  int counter = 0;

  int delta_nb[6];
  double pot[6];
  double mask_nb[6];

  /* We don't change anything for the non-split cells, because have been masked. */
  for (int i=0; i<coarse->cell_count;i++) {
    if (coarse->cells[i]->mask_value < 0) continue; //We only restrict if the coarse cell is not masked
    //message("Restricting to unmasked cell %d", i);
    counter += 1;
    coarse->cells[i]->conv_residual = 0.;
    coarse->cells[i]->CIC_density = 0.;
    for (int j=0;j<8;j++) { //Loop over the fine cells and add the residual
      struct cell *child = coarse->cells[i]->progeny[j];
      if (child->mask_value > 0) {
        for (int k=0;k<6;k++) {
          mask_nb[k] = child->neighbours[k]->mask_value;
          if (mask_nb[k]<=0) {
            delta_nb[k] = 0;
          }
          else {
            delta_nb[k] = 1;
            pot[k] = child->neighbours[k]->conv_residual;
          }
        }
        double potential_term = ((double) delta_nb[0]*pot[0] + (double) delta_nb[1]*pot[1] + (double) delta_nb[2]*pot[2] + (double) delta_nb[3]*pot[3] + 
                                (double) delta_nb[4]*pot[4] + (double) delta_nb[5]*pot[5]);
        double ghost_correction = (((((double) 1-delta_nb[0])*mask_nb[0] + ((double) 1-delta_nb[1])*mask_nb[1] + ((double) 1-delta_nb[2])*mask_nb[2] +
                                      + ((double) 1-delta_nb[3])*mask_nb[3] + ((double) 1-delta_nb[4])*mask_nb[4] + ((double) 1-delta_nb[5])*mask_nb[5]))/(child->mask_value));
        /* The below calculates the residual and stores it in the density slot */
        /* Add zero if a child is masked? */
        coarse->cells[i]->CIC_density += ((potential_term - (6.-ghost_correction)*child->conv_residual)/(delta*delta) - child->CIC_density)/8.;
      }
    }
    //Don't forget to multiply by -1
    coarse->cells[i]->CIC_density *= -1.;
  }
  message("The mean of the array being transferred is %lf", sum/counter);
}

void transfer_residual_array(struct AMR_levels fine, struct AMR_levels *coarse, double delta) {
  message("Transferring the residual to level %d", coarse->depth);
  //float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  double mean_density = fine.mean_density;
  double sum = 0.;
  int counter = 0;
  message("Delta when transferring is %lf", delta);

  /* We don't change anything for the non-split cells, because they have been masked. */
  for (int i=0; i<coarse->cell_count;i++) {
    if (coarse->cells[i]->mask_value <= 0) continue; //We only restrict if the coarse cell is not masked
    //message("Restricting to unmasked cell %d", i);
    counter += 1;
    coarse->cells[i]->conv_residual = 0.;
    coarse->cells[i]->CIC_density = 0.;
    for (int j=0;j<8;j++) { //Loop over the fine cells and add the residual
      struct cell *child = coarse->cells[i]->progeny[j];
      if (child->mask_value > 0) {
        /* The below calculates the residual and stores it in the density slot */
        /* Add zero if a child is masked? */
        coarse->cells[i]->CIC_density += (((child->neighbours[0]->CIC_potential + child->neighbours[1]->CIC_potential + 
                  child->neighbours[2]->CIC_potential + child->neighbours[3]->CIC_potential + 
                  child->neighbours[4]->CIC_potential + child->neighbours[5]->CIC_potential - 
                  6.0*child->CIC_potential)/(delta*delta) - mean_density*(child->CIC_density - 1.))/8.0);
      }
    }
    //Don't forget to multiply by -1
    coarse->cells[i]->CIC_density *= -1.;
  }
  message("The mean of the array being transferred is %lf", sum/counter);
}

double get_patch_residual(struct AMR_levels level, double delta) {
  double mean_density = level.mean_density;
  double residual = 0.;
  double mean_res = 0.;

  //message("Exporting AMR nonuniform FFT residual");
  //FILE *file_swift_density2 = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/AMR_16_nonuniform_residual.txt", "w");

  for (int i=0; i<level.cell_count; i++) {
    struct cell *current_cell = level.cells[i];
    if (current_cell->ghost) message("Found a stray ghost :(");
    if (level.depth == 6) {
      //message("Checking cell %d", i);
      if (current_cell->neighbours[0] == NULL) message("Neighbour 0 is NULL");
      if (current_cell->neighbours[1] == NULL) message("Neighbour 1 is NULL");
      if (current_cell->neighbours[2] == NULL) message("Neighbour 2 is NULL");
      if (current_cell->neighbours[3] == NULL) message("Neighbour 3 is NULL");
      if (current_cell->neighbours[4] == NULL) message("Neighbour 4 is NULL");
      if (current_cell->neighbours[5] == NULL) message("Neighbour 5 is NULL");
    }
    double res = ((current_cell->neighbours[0]->CIC_potential + current_cell->neighbours[1]->CIC_potential + 
                  current_cell->neighbours[2]->CIC_potential + current_cell->neighbours[3]->CIC_potential + 
                  current_cell->neighbours[4]->CIC_potential + current_cell->neighbours[5]->CIC_potential - 
                  6.0*current_cell->CIC_potential)/(delta*delta) - mean_density*(current_cell->CIC_density - 1.));
    //if (i<10) {
      //message("Doing cell (%lf, %lf, %lf)", current_cell->loc[0], current_cell->loc[1], current_cell->loc[2]);
      //message("Calculating using %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g %.5g", current_cell->neighbours[0]->CIC_potential, current_cell->neighbours[1]->CIC_potential, current_cell->neighbours[2]->CIC_potential, current_cell->neighbours[3]->CIC_potential, current_cell->neighbours[4]->CIC_potential, current_cell->neighbours[5]->CIC_potential, current_cell->CIC_potential, mean_density, current_cell->CIC_density);
    //}
    //fprintf(file_swift_density2, "%lf %.15g %.15g %.15g \n", res, current_cell->loc[0], current_cell->loc[1], current_cell->loc[2]);
    //message("We just added %lf", temp);
    mean_res += res;
    residual += res*res;
  }
  double final_res = sqrt(residual/(level.cell_count));
  //message("The mean residual is %lf", mean_res/level.cell_count);
  //message("The residual is %lf", final_res);
  //fclose(file_swift_density2);
  //message("Exported the AMR nonuniform residual data.");
  //sleep(15);
  return final_res;
}

void perform_patch_sweep(struct AMR_levels *level, double delta) {
  double mean_density = level->mean_density;

  for (int col=0; col<2; col++){
    for (int i=0; i<level->cell_count;i++) {
      struct cell *current_cell = level->cells[i];
      int location =  (int) ((current_cell->loc[0]+0.0001)/delta) + (int) ((current_cell->loc[1]+0.0001)/delta) + (int) ((current_cell->loc[2]+0.0001)/delta);
      if (location%2 != col) continue;
      if (level->cells[i] == NULL) message("Calculating potential for nonexistent cell");
      for (int j=0;j<6;j++) {
        if (level->cells[i]->neighbours[j] == NULL) message("Neighbour %d nonexistent", j);
        //if (level->cells[i]->neighbours[j]->ghost) message("Neighbour %d of (%lf, %lf, %lf) is a ghost with location (%lf, %lf, %lf) and potential %lf", j, level->cells[i]->loc[0], level->cells[i]->loc[1], level->cells[i]->loc[2],level->cells[i]->neighbours[j]->loc[0], level->cells[i]->neighbours[j]->loc[1], level->cells[i]->neighbours[j]->loc[2], level->cells[i]->neighbours[j]->CIC_potential);
      }
      level->cells[i]->CIC_potential = ((current_cell->neighbours[0]->CIC_potential + current_cell->neighbours[1]->CIC_potential + 
                                          current_cell->neighbours[2]->CIC_potential + current_cell->neighbours[3]->CIC_potential +
                                          current_cell->neighbours[4]->CIC_potential + current_cell->neighbours[5]->CIC_potential - 
                                          delta*delta*mean_density*(current_cell->CIC_density - 1.))/6.0);
      //message("Set cell %d with loc (%lf, %lf, %lf) to %lf", i, level->cells[i]->loc[0], level->cells[i]->loc[1], level->cells[i]->loc[2],level->cells[i]->CIC_potential);
    }
  }
}

void extrapolate_mask_values(struct AMR_levels *coarse) {
  int coarse_count = coarse->cell_count;

  /* Value must remain -1 if a ghost cell, so do nothing in that case */
  for (int i=0; i<coarse_count;i++) {
    if (!(coarse->cells[i]->ghost)) {
      /* Set to zero then sum over children */
      if (!(coarse->cells[i]->split)) {
        //if (i==0) message("Cell zero is not split so mask value set to -1");
        coarse->cells[i]->mask_value = -1.; 
        continue;
      }
      coarse->cells[i]->mask_value = 0.; 
      for (int j=0;j<8;j++) {
        //if (i==0) message("Adding value %lf for cell zero", coarse->cells[i]->progeny[j]->mask_value);
        coarse->cells[i]->mask_value += coarse->cells[i]->progeny[j]->mask_value;
      }
      coarse->cells[i]->mask_value = (coarse->cells[i]->mask_value)/8.;
    }
  }
  if (coarse->depth == 2) message("The mask value of cell zero is %lf", coarse->cells[0]->mask_value);
}

void initialise_link_neighbours(struct AMR_levels *fine, double box_size, int gridsize) {
  message("Linking neighbours on level %d", fine->depth);
  double fac = box_size/gridsize;
  int cells_tot = fine->cell_count + fine->ghost_count;
  for (int i=0;i<cells_tot;i++) {
    if (fine->cells[i]->ghost) continue; //Only create one-way link to ghost cells
    double cell_x = fine->cells[i]->loc[0];
    double cell_y = fine->cells[i]->loc[1];
    double cell_z = fine->cells[i]->loc[2];

    message("Checking neighbours for the cell with (%lf,%lf,%lf)", cell_x,cell_y, cell_z);

    /* Check all six neighbours and link them. Prevent linking if this already happened. */
    double search_loc[3] = {cell_x + fac,0,0};
    if (fine->cells[i]->neighbours[0] == NULL) {
      link_neighbour(fine->cells[i], fine, search_loc, fac, 0);
      if (fine->cells[i]->neighbours[0] == NULL) message("linking failed");
    }
    search_loc[0] = cell_x - fac;
    if (fine->cells[i]->neighbours[1] == NULL) link_neighbour(fine->cells[i], fine, search_loc, fac, 1);

    search_loc[0] = cell_x;
    search_loc[1] = cell_y + fac;
    if (fine->cells[i]->neighbours[2] == NULL) link_neighbour(fine->cells[i], fine, search_loc, fac, 2);
    search_loc[1] = cell_y - fac;
    if (fine->cells[i]->neighbours[3] == NULL) link_neighbour(fine->cells[i], fine, search_loc, fac, 3);

    search_loc[1] = cell_y;
    search_loc[2] = cell_z + fac;
    if (fine->cells[i]->neighbours[4] == NULL) link_neighbour(fine->cells[i], fine, search_loc, fac, 4);
    search_loc[2] = cell_z - fac;
    if (fine->cells[i]->neighbours[5] == NULL) link_neighbour(fine->cells[i], fine, search_loc, fac, 5);
  }
}

void link_neighbour(struct cell *link_cell, struct AMR_levels *fine, double search_loc[3], double fac, int which_neighbour) {
  int cells_tot = fine->cell_count + fine->ghost_count;
  for (int i=0; i<cells_tot; i++) {
    double x_loc = fine->cells[i]->loc[0];
    double y_loc = fine->cells[i]->loc[1];
    double z_loc = fine->cells[i]->loc[2];

    double cell_check[6] = {x_loc, x_loc + fac, y_loc, y_loc + fac, z_loc, z_loc + fac};
    double overlap = get_overlap(cell_check, search_loc, fine->cells[0]->width);
    if (overlap > 0.1) {
      /* Link back and forth */
      link_cell->neighbours[which_neighbour] = fine->cells[i]; 
      if (!(fine->cells[i]->ghost)) {
        int backward_neighbour = (which_neighbour%2 == 0) ? which_neighbour+1 : which_neighbour-1;
        fine->cells[i]->neighbours[backward_neighbour] = link_cell;
      }
      //message("We are searching for cell (%lf,%lf,%lf) and found (%lf,%lf,%lf), yay!", search_loc[0], search_loc[1], search_loc[2], x_loc, y_loc, z_loc);
      break;
    }  
  }
}

void set_patch_guess(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, int nr_cells, int gridsize, double box_size, int min_depth) {
  double fac = box_size/gridsize;
  double fac_coarse = box_size/(gridsize/2);
  int cdim[3] = {fine->cdim, fine->cdim, fine->cdim};
  /* Identify and set the ghost cells on the fine level */
  int neighbour;

  /*Set the progeny of non-split cells to NULL manually*/
  for (int i=0; i<coarse->cell_count; i++) {
    if (coarse->cells[i]->split) continue;
    for (int j=0; j<8; j++) {
      coarse->cells[i]->progeny[j] = NULL;
    }
  }

  /* Check if prolongation is the issue */
  /*if (fine->depth <= min_depth) {
    message("Doing this");
    int cdimh[3] = {fine->cdim, fine->cdim, fine->cdim};
    cdim[0] = coarse->cdim;
    cdim[1] = coarse->cdim;
    cdim[2] = coarse->cdim;
    int N = coarse->cdim;

    for (int i = 0; i<fine->cdim; i++) {
      for (int j = 0; j<fine->cdim; j++) {
        for (int k=0; k<fine->cdim; k++) {
          int iH = i/2;
          int jH = j/2;
          int kH = k/2;
          double dx = (i%2) * 0.5;
          double dy = (j%2) * 0.5;
          double dz = (k%2) * 0.5;

          int iHp = (iH +1)%N;
          int jHp = (jH +1)%N;
          int kHp = (kH +1)%N;

          fine->cells[cell_getid(cdimh,i,j,k)]->CIC_potential = ((1.-dx)*(1.-dy)*(1.-dz)*coarse->cells[cell_getid(cdim,iH,jH,kH)]->CIC_potential + dx*(1.-dy)*(1.-dz)*coarse->cells[cell_getid(cdim,iHp, jH, kH)]->CIC_potential
                                    + (1.-dx)*dy*(1.-dz)*coarse->cells[cell_getid(cdim,iH,jHp,kH)]->CIC_potential + (1.-dx)*(1.-dy)*dz*coarse->cells[cell_getid(cdim,iH,jH,kHp)]->CIC_potential
                                    + dx*dy*(1.-dz)*coarse->cells[cell_getid(cdim,iHp,jHp,kH)]->CIC_potential + dx*(1.-dy)*dz*coarse->cells[cell_getid(cdim,iHp,jH, kHp)]->CIC_potential
                                    + (1.-dx)*dy*dz*coarse->cells[cell_getid(cdim,iH,jHp,kHp)]->CIC_potential+ dx*dy*dz*coarse->cells[cell_getid(cdim,iHp,jHp,kHp)]->CIC_potential);
        }
      }
    }
    return;
  }*/

  message("Going to do this");
  for (int i=0;i<nr_cells;i++) {
    /* Box wrap the multipole's position */
    const double search_loc_x = box_wrap(fine->cells[i]->loc[0], 0., s->dim[0]);
    const double search_loc_y = box_wrap(fine->cells[i]->loc[1], 0., s->dim[1]);
    const double search_loc_z = box_wrap(fine->cells[i]->loc[2], 0., s->dim[2]);

    int cell_i = (int) ((search_loc_x+0.0001)/fac);
    int cell_j = (int) ((search_loc_y+0.0001)/fac);
    int cell_k = (int) ((search_loc_z+0.0001)/fac);

    int i_plus = (cell_i + 1) % cdim[0];
    int i_min = (cell_i-1>=0) ? (cell_i-1) % cdim[0] : (cell_i-1) % cdim[0] + cdim[0];
    int j_plus = (cell_j + 1) % cdim[0];
    int j_min = (cell_j-1>=0) ? (cell_j-1) % cdim[0] : (cell_j-1) % cdim[0] + cdim[0];
    int k_plus = (cell_k + 1) % cdim[0];
    int k_min = (cell_k-1>=0) ? (cell_k-1) % cdim[0] : (cell_k-1) % cdim[0] + cdim[0];

    size_t search_id[6] = {cell_getid(cdim, i_plus, cell_j, cell_k), cell_getid(cdim, i_min, cell_j, cell_k),
                        cell_getid(cdim, cell_i, j_plus, cell_k), cell_getid(cdim, cell_i, j_min, cell_k),
                        cell_getid(cdim, cell_i, cell_j, k_plus), cell_getid(cdim, cell_i, cell_j, k_min)};

    struct cell *parent = fine->cells[i]->parent;

    int neighbour_counter = 0;
    int no_neighbour_counter = 0;
    //message("Checking neighbours for the cell with (%lf,%lf,%lf)", search_loc_x,search_loc_y, search_loc_z);

    //if (cell_i == 0 && cell_j == 0 && cell_k == 2) message("The neighbour %d has location (%lf, %lf, %lf)", 5, fine->cells[i]->neighbours[5]->loc[0], fine->cells[i]->neighbours[5]->loc[1], fine->cells[i]->neighbours[5]->loc[2]);

    for (int j=0; j<6; j++) {
      if (fine->cells[i]->neighbours[j] == NULL) {
        //message("Neighbour %d not connected", j);
        //if (cell_i == 0 && cell_j == 0 && cell_k == 1) message("Neighbour %d not connected", j);
        //message("Searching for neighbour %d for cell %d", j, i);
        int pre_smoothing = 0;
        neighbour = find_neighbour(fine, parent, fine->cells[i], search_id[j], fac, j, pre_smoothing, 0);
        if (cell_i == 0 && cell_j == 0 && cell_k == 1) message("Neighbour found? %d", neighbour);
        no_neighbour_counter += 1;
        //if (fine->depth == 8) message("For cell %d, neighbour %d was not linked", i, j);
      }
      else {
        neighbour = 1;
        neighbour_counter +=1;
        //message("For cell %d, neighbour %d already linked", i, j);
      }
      //message("The neighbour %d value %d for cell %d is", j, neighbour, i);
      if (!neighbour) {
        //if (fine->depth == 8) message("For cell %d, neighbour %d did not exist", i, j);
        //if (search_id[j]<0) message("Search ID %d is %lu and we had entries %d,%d,%d and cdim %d", j, search_id[j], i_min, cell_j, cell_k, cdim[0]);
        int check_parent = 0;
        create_ghost(s, coarse, fine, parent->neighbours[j], search_id[j], fac_coarse, min_depth, j, check_parent, 0);
        fine->cells[i]->neighbours[j] = fine->cells[nr_cells+fine->ghost_count-1];
        int backward_neighbour = (j%2 == 0) ? j+1 : j-1;
        //WAIT DOES THIS LINK A CELL TO ITSELF? Should no longer now
        fine->cells[nr_cells+fine->ghost_count-1]->neighbours[backward_neighbour] = fine->cells[i];
        //message("We just linked the cell at (%lf, %lf, %lf) to the cell at (%lf, %lf, %lf) for neighbour %d", fine->cells[nr_cells+fine->ghost_count-1]->loc[0], fine->cells[nr_cells+fine->ghost_count-1]->loc[1], fine->cells[nr_cells+fine->ghost_count-1]->loc[2], fine->cells[i]->neighbours[j]->loc[0], fine->cells[i]->neighbours[j]->loc[1], fine->cells[i]->neighbours[j]->loc[2], backward_neighbour);
        //sleep(15);
        for (int k=0; k<6; k++) {
          if (fine->cells[nr_cells+fine->ghost_count-1]->neighbours[k] != NULL) {
            //message("Newly made ghost has neighbour %d linked", k);
            //message("This neighbour has location (%lf, %lf, %lf)", fine->cells[nr_cells+fine->ghost_count-1]->neighbours[k]->loc[0], fine->cells[nr_cells+fine->ghost_count-1]->neighbours[k]->loc[1], fine->cells[nr_cells+fine->ghost_count-1]->neighbours[k]->loc[2]);
          }
        }
      }
      //if (!neighbour && fine->cells[i]->neighbours[j] == NULL) message("Failed to create the neighbour %d for cell %d", j, i);
      if (fine->cells[i]->neighbours[j] == NULL) {
      message("Failed to create or link the neighbour %d for cell %d with location (%lf, %lf, %lf)", j, i, fine->cells[i]->loc[0], fine->cells[i]->loc[1], fine->cells[i]->loc[2]);
      }
    }
  }
  message("In total, %d ghosts were created", fine->ghost_count);

  for (int i=0; i<nr_cells; i++) {
    struct cell *c = fine->cells[i];
    for (int j=0; j<6; j++) {
      struct cell *b = c->neighbours[j];
      int dir = (j%2 == 0) ? j : j-1;
      int check[4] = {(dir + 2)%6, (dir + 3)%6, (dir + 4)%6, (dir + 5)%6};
      int nr_neighbours = 4;
      //message("Going to do the diagonal for the cell at (%lf, %lf, %lf) using neighbour at (%lf, %lf, %lf) with ghost value %d", c->loc[0], c->loc[1], c->loc[2], b->loc[0], b->loc[1], b->loc[2], b->ghost);
      check_diagonal1(s, coarse, fine, b, nr_neighbours, check, min_depth, 0);
      //sleep(15);
    }
  }
  message("In total, %d ghosts were created", fine->ghost_count);
  //sleep(15);
  

  /* Check the neighbour assignment */
  /*for (int i=0; i<nr_cells; i++) {
    
    const double search_loc_x = box_wrap(fine->cells[i]->loc[0], 0., s->dim[0]);
    const double search_loc_y = box_wrap(fine->cells[i]->loc[1], 0., s->dim[1]);
    const double search_loc_z = box_wrap(fine->cells[i]->loc[2], 0., s->dim[2]);

    int cell_i = (int) ((search_loc_x+ 0.0001)/fac );
    int cell_j = (int) ((search_loc_y+ 0.0001)/fac );
    int cell_k = (int) ((search_loc_z+ 0.0001)/fac );

    int i_plus = (cell_i+1) % cdim[0];
    int i_min = (cell_i-1>=0) ? (cell_i-1) % cdim[0] : (cell_i-1) % cdim[0] + cdim[0];
    int j_plus = (cell_j+1) % cdim[1];
    int j_min = (cell_j-1>=0) ? (cell_j-1) % cdim[1] : (cell_j-1) % cdim[1] + cdim[1];
    int k_plus = (cell_k+1) % cdim[2];
    int k_min = (cell_k-1>=0) ? (cell_k-1) % cdim[2] : (cell_k-1) % cdim[2] + cdim[2];

    //int id_cell = cell_getid(cdim, cell_i, cell_j, cell_k);

    for (int j=0; j<6; j++) {
      struct cell *cell_check = fine->cells[i]->neighbours[j];
      const double loc_x_nb = box_wrap(cell_check->loc[0], 0., s->dim[0]);
      const double loc_y_nb = box_wrap(cell_check->loc[1], 0., s->dim[1]);
      const double loc_z_nb = box_wrap(cell_check->loc[2], 0., s->dim[2]);

      int cell_i_nb = (int) ((loc_x_nb+ 0.0001)/fac );
      int cell_j_nb = (int) ((loc_y_nb+ 0.0001)/fac );
      int cell_k_nb = (int) ((loc_z_nb+ 0.0001)/fac );

      int id_nb = cell_getid(cdim, cell_i_nb, cell_j_nb, cell_k_nb);
      int id_nb_check;

      if (j==0) id_nb_check = cell_getid(cdim, i_plus, cell_j, cell_k);
      if (j==1) id_nb_check = cell_getid(cdim, i_min, cell_j, cell_k);
      if (j==2) id_nb_check = cell_getid(cdim, cell_i, j_plus, cell_k);
      if (j==3) id_nb_check = cell_getid(cdim, cell_i, j_min, cell_k);
      if (j==4) id_nb_check = cell_getid(cdim, cell_i, cell_j, k_plus);
      if (j==5) id_nb_check = cell_getid(cdim, cell_i, cell_j, k_min);
      //Get id of neighbour by +-1 id_cell
      //Compare
      if (id_nb != id_nb_check) {
        message("Found a wrong neighbour.");
        message("Cell has (%lf,%lf, %lf) and indices (%d,%d,%d)", search_loc_x, search_loc_y, search_loc_z, cell_i, cell_j, cell_k);
        message("Neighbour %d has (%lf,%lf, %lf) and indices (%d,%d,%d)", j, loc_x_nb, loc_y_nb, loc_z_nb, cell_i_nb, cell_j_nb, cell_k_nb);
        //message("Cell has (%lf,%lf, %lf) and nb %d has (%lf, %lf, %lf)", search_loc_x, search_loc_y, search_loc_z, j, loc_x_nb, loc_y_nb, loc_z_nb);
        message("The compared cell IDs are %d and %d", id_nb, id_nb_check);
      }
    }


  }*/


  //int cdimH[3] = {gridsize/2, gridsize/2, gridsize/2}; //For finding cells on the coarse grid

  int trilinear = 1;

  if (trilinear) {
    interpolate_trilinear(coarse, fine);
  }
  else {
    for (int i=0; i<coarse->cell_count + coarse->ghost_count; i++) {
      struct cell *parent = coarse->cells[i];
      if (!parent->split) continue;
      //int cell_x = (int) ((parent->loc[0]+0.01)/fac_coarse);
      //int cell_y = (int) ((parent->loc[1]+0.01)/fac_coarse);
      //int cell_z = (int) ((parent->loc[2]+0.01)/fac_coarse);

      for (int j=0; j<8; j++) {
        struct cell *child = parent->progeny[j];
        child->CIC_potential = 0.;
        //if (!child->ghost) continue;

        int offset[3] = {0,0,0}; //Do we need to interpolate in this direction?
        if (j & 4) offset[0] = 1;
        if (j & 2) offset[1] = 1;
        if (j & 1) offset[2] = 1;

        child->CIC_potential = 1./3. * (((double) offset[0])*(-1./8.*parent->neighbours[1]->CIC_potential + 3./4.*parent->CIC_potential + 3./8.*parent->neighbours[0]->CIC_potential) +
                                ((double) offset[1])*(-1./8.*parent->neighbours[3]->CIC_potential + 3./4.*parent->CIC_potential + 3./8.*parent->neighbours[2]->CIC_potential) +
                                ((double) offset[2])*(-1./8.*parent->neighbours[5]->CIC_potential + 3./4.*parent->CIC_potential + 3./8.*parent->neighbours[4]->CIC_potential) +
                                (double) ((1-offset[0]) + (1-offset[1]) + (1-offset[2])) * parent->CIC_potential);

        /* Set mask value */
        child->mask_value = 1.;
      }
    }
  }
  
  nr_cells += fine->ghost_count;
  message("Done interpolating");

  //if (fine->depth == 2) {
    //message("Exporting AMR density and potential data");
    //FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/optimised_ghost_test/location_onlyghosts_level2", "w");
    //for (int j = 0; j < nr_cells; j++) {
      //if (fine->cells[j]->ghost) {
        //double search_loc_x = fine->cells[j]->loc[0];
        //double search_loc_y = fine->cells[j]->loc[1];
        //double search_loc_z = fine->cells[j]->loc[2];
        //int cell_i = (int) ((search_loc_x+0.01)/fac);
        //int cell_j = (int) ((search_loc_y+0.01)/fac);
        //int cell_k = (int) ((search_loc_z+0.01)/fac);

        //int cdim[3] = {fine->cdim, fine->cdim, fine->cdim};
        //int loc_1d = cell_getid(cdim, cell_i, cell_j, cell_k);

        //fprintf(file_swift_density, "cell %d at (%lf , %lf , %lf) has %lf and %lf \n", j, fine->cells[j]->loc[0], fine->cells[j]->loc[1], fine->cells[j]->loc[2], fine->cells[j]->CIC_density, fine->cells[j]->CIC_potential);
        //fprintf(file_swift_density, "%d\n", loc_1d);
      //}
    //}
    //fclose(file_swift_density);
    //message("Exported the AMR density and potential data.");
    //sleep(5);
  //}
}

void check_diagonal1(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, struct cell *c, int nr_neighbours, int neighbours[nr_neighbours], int min_depth, int double_diag) {
  double fac = s->dim[0]/fine->cdim;
  double fac_coarse = s->dim[0]/coarse->cdim;
  int cdim[3] = {fine->cdim, fine->cdim, fine->cdim};

  const int nbs[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};
  //if (!double_diag) message("Neighbours are %d, %d, %d, %d", neighbours[0], neighbours[1], neighbours[2], neighbours[3]);
  for (int i=0; i<nr_neighbours; i++) {
    int which_neighbour = neighbours[i];
    //message("searching for neighbour %d for cell at (%lf, %lf, %lf)", which_neighbour, c->loc[0], c->loc[1], c->loc[2]);
    //if (c->neighbours[which_neighbour] != NULL) message("Found the neighbour at (%lf, %lf, %lf)", c->neighbours[which_neighbour]->loc[0], c->neighbours[which_neighbour]->loc[1], c->neighbours[which_neighbour]->loc[2]);
    if (c->neighbours[which_neighbour] == NULL) {
      const double search_loc_x = box_wrap(c->loc[0], 0., s->dim[0]);
      const double search_loc_y = box_wrap(c->loc[1], 0., s->dim[1]);
      const double search_loc_z = box_wrap(c->loc[2], 0., s->dim[2]);

      int cell_i = (int) ((search_loc_x+0.0001)/fac);
      int cell_j = (int) ((search_loc_y+0.0001)/fac);
      int cell_k = (int) ((search_loc_z+0.0001)/fac);

      int search_i = (cell_i + nbs[which_neighbour][0] + cdim[0]) % cdim[0];
      int search_j = (cell_j + nbs[which_neighbour][1] + cdim[1]) % cdim[1];
      int search_k = (cell_k + nbs[which_neighbour][2] + cdim[2]) % cdim[2];

      size_t search_id = cell_getid(cdim, search_i, search_j, search_k);
      int check_parent = 1;
      int verbose = 0;
      create_ghost(s, coarse, fine, c->parent, search_id, fac_coarse, min_depth, which_neighbour, check_parent, verbose);

      c->neighbours[which_neighbour] = fine->cells[fine->cell_count+fine->ghost_count-1];
      int backward_neighbour = (which_neighbour%2 == 0) ? which_neighbour+1 : which_neighbour-1;
      fine->cells[fine->cell_count+fine->ghost_count-1]->neighbours[backward_neighbour] = c;
      //message("The new ghost has location (%lf, %lf, %lf)", fine->cells[fine->cell_count+fine->ghost_count-1]->loc[0], fine->cells[fine->cell_count+fine->ghost_count-1]->loc[1], fine->cells[fine->cell_count+fine->ghost_count-1]->loc[2]);
      //message("Parent at location (%lf, %lf, %lf)", fine->cells[fine->cell_count+fine->ghost_count-1]->parent->loc[0], fine->cells[fine->cell_count+fine->ghost_count-1]->parent->loc[1], fine->cells[fine->cell_count+fine->ghost_count-1]->parent->loc[2]);
      for (int k=0; k<6; k++) {
        if (fine->cells[fine->cell_count+fine->ghost_count-1]->neighbours[k] != NULL) {
          //message("Newly made ghost has neighbour %d linked", k);
          //message("This neighbour has location (%lf, %lf, %lf)", fine->cells[fine->cell_count+fine->ghost_count-1]->neighbours[k]->loc[0], fine->cells[fine->cell_count+fine->ghost_count-1]->neighbours[k]->loc[1], fine->cells[fine->cell_count+fine->ghost_count-1]->neighbours[k]->loc[2]);
        }
      }
      //if (!neighbour && fine->cells[i]->neighbours[j] == NULL) message("Failed to create the neighbour %d for cell %d", j, i);
      if (fine->cells[i]->neighbours[which_neighbour] == NULL) {
        message("Failed to create or link the neighbour %d for cell %d with location (%lf, %lf, %lf)", which_neighbour, i, fine->cells[i]->loc[0], fine->cells[i]->loc[1], fine->cells[i]->loc[2]);
      }
    }
    if (double_diag) {
      //message("\n");
      continue;
    }
    struct cell *nb = c->neighbours[which_neighbour];
    int dir = (i%2 == 0) ? i : i-1;

    int check[2] = {neighbours[(dir + 2)%4], neighbours[(dir + 3)%4]};
    //message("Going to search for neighbours %d and %d", check[0], check[1]);
    int nr_neighbours2 = 2;
    //message("Going to do the double diagonal");
    check_diagonal1(s, coarse, fine, nb, nr_neighbours2, check, min_depth, 1);

  }
}

void interpolate_trilinear(struct AMR_levels *coarse, struct AMR_levels *fine) {
  for (int i=0; i<fine->cell_count + fine->ghost_count; i++) {
    struct cell *child = fine->cells[i];
    struct cell *parent = child->parent;
    child->CIC_potential = 0.;

    /* Get the location/offset */
    double offset[3] = {child->loc[0] - parent->loc[0], child->loc[1] - parent->loc[1], child->loc[2] - parent->loc[2]};
    double dx = (offset[0] > 0) ? 0.5 : 0.;
    double dy = (offset[1] > 0) ? 0.5 : 0.;
    double dz = (offset[2] > 0) ? 0.5 : 0.;
    //message("The offsets are (%lf, %lf, %lf) and we calculated (%lf,%lf,%lf)", offset[0], offset[1], offset[2], dx, dy, dz);

    if (parent->neighbours[0] == NULL) message("Neighbour 0 null for cell %d", i);
    if (parent->neighbours[2] == NULL) message("Neighbour 2 null for cell %d", i);
    if (parent->neighbours[4] == NULL) message("Neighbour 4 null for cell %d", i);
    if (parent->neighbours[0]->neighbours[2] == NULL) message("Neighbour 0,2 null for cell %d", i);
    if (parent->neighbours[0]->neighbours[4] == NULL) message("Neighbour 0,4 null for cell %d", i);
    if (parent->neighbours[2]->neighbours[4] == NULL) message("Neighbour 2,4 null for cell %d", i);
    if (parent->neighbours[0]->neighbours[2]->neighbours[4] == NULL) message("Neighbour 0,2,4 null for cell %d", i);
    child->CIC_potential = ((1.-dx)*(1.-dy)*(1.-dz)*parent->CIC_potential + dx*(1.-dy)*(1.-dz)*parent->neighbours[0]->CIC_potential
                              + (1.-dx)*dy*(1.-dz)*parent->neighbours[2]->CIC_potential + (1.-dx)*(1.-dy)*dz*parent->neighbours[4]->CIC_potential
                              + dx*dy*(1.-dz)*parent->neighbours[0]->neighbours[2]->CIC_potential + dx*(1.-dy)*dz*parent->neighbours[0]->neighbours[4]->CIC_potential
                              + (1.-dx)*dy*dz*parent->neighbours[2]->neighbours[4]->CIC_potential+ dx*dy*dz*parent->neighbours[0]->neighbours[2]->neighbours[4]->CIC_potential);

    if (!child->ghost) child->mask_value = 1.;
    //else message("Set ghost at loc (%lf, %lf, %lf) to value %lf", child->loc[0], child->loc[1], child->loc[2], child->CIC_potential);
  }
}

void create_ghost(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, struct cell *parent, size_t search_id, double fac_coarse, int min_depth, int which_neighbour, int check_parent, int verbose) {
  int ghost_count = fine->ghost_count;
  int nr_cells = fine->cell_count;
  int cdim = fine->cdim;
  size_t search_loc[3] = {(search_id/(cdim*cdim))%cdim, (search_id/cdim)%cdim, search_id % cdim};
  double fac = s->dim[0]/cdim;
  fine->cells = realloc(fine->cells, (nr_cells+ghost_count+1)*sizeof(struct cell *));
  
  space_getcells(s, 1, fine->cells, 0, nr_cells + ghost_count, 1);
  if (fine->cells[nr_cells+ghost_count] == NULL) error("Not created");
  struct cell *c = fine->cells[nr_cells + ghost_count];

  c->loc[0] = (double) search_loc[0] * fac; 
  c->loc[1] = (double) search_loc[1] * fac; 
  c->loc[2] = (double) search_loc[2] * fac; 
  c->width[0] = fine->cells[0]->width[0]; 
  c->width[1] = fine->cells[0]->width[1]; 
  c->width[2] = fine->cells[0]->width[2]; 

  c->ghost = 1;
  c->mask_value = -1;
  if (check_parent) { /* In this case the parent is the parent of the cell of which we have now created the neighbour */
    if (c->loc[0] >= parent->loc[0] && c->loc[0] < parent->loc[0] + parent->width[0] && c->loc[1] >= parent->loc[1] && c->loc[1] < parent->loc[1] + parent->width[1] && c->loc[2] >= parent->loc[2] && c->loc[2] < parent->loc[2] + parent->width[2]) {
      c->parent = parent;
    }
    else {
      if (verbose) message("Doing this");
      c->parent = parent->neighbours[which_neighbour];
      parent = parent->neighbours[which_neighbour];
    }
  }
  else c->parent = parent;
  /* parent->split remains 0, but we do add some children. Which child is this? */
  int x_id = search_loc[0]%2;
  int y_id = search_loc[1]%2;
  int z_id = search_loc[2]%2;
  if (parent->progeny[2*2*x_id + 2*y_id + z_id] != NULL) {
    message("We are doing ID %d and the search loc was (%lu, %lu, %lu)", 2*2*x_id + 2*y_id + z_id, search_loc[0], search_loc[1], search_loc[2]);
    message("Created ghost cell at (%lf, %lf, %lf)", (double) search_loc[0] * fac, (double) search_loc[1] * fac, (double) search_loc[2] * fac);
    message("The existing cell had location (%lf, %lf, %lf)", parent->progeny[2*2*x_id + 2*y_id + z_id]->loc[0], parent->progeny[2*2*x_id + 2*y_id + z_id]->loc[1], parent->progeny[2*2*x_id + 2*y_id + z_id]->loc[2]);
    error("Something went wrong");
  }
  parent->progeny[2*2*x_id + 2*y_id + z_id] = c;
  //message("Assigned progeny to parent of ghost");
  /* Link ghost cell to ALL possible neighbours */
  create_link(s, fine, &c, 1, which_neighbour, verbose); //Cell doesn't need to search itself
  fine->ghost_count += 1;
}

/* Link cells to all cells in level. Cells could be level->cells itself, but doesn't need to be 
   nr_cells is the number of cells that need linking. */
void create_link(struct space *s, struct AMR_levels *level, struct cell **cells, int nr_cells, int which_neighbour, int verbose) {
  //message("Going to create link");
  double box_size = s->dim[0];
  double fac = box_size/level->cdim;
  int cdim[3] = {level->cdim, level->cdim, level->cdim};

  for (int i=0; i<nr_cells;i++) {
    int cell_i = (int) ((cells[i]->loc[0]+0.0001)/fac);
    int cell_j = (int) ((cells[i]->loc[1]+0.0001)/fac);
    int cell_k = (int) ((cells[i]->loc[2]+0.0001)/fac);

    int i_plus = (cell_i + 1) % cdim[0];
    int i_min = (cell_i-1>=0) ? (cell_i-1) % cdim[0] : (cell_i-1) % cdim[0] + cdim[0];
    int j_plus = (cell_j + 1) % cdim[0];
    int j_min = (cell_j-1>=0) ? (cell_j-1) % cdim[0] : (cell_j-1) % cdim[0] + cdim[0];
    int k_plus = (cell_k + 1) % cdim[0];
    int k_min = (cell_k-1>=0) ? (cell_k-1) % cdim[0] : (cell_k-1) % cdim[0] + cdim[0];

    int search_id[6] = {cell_getid(cdim, i_plus, cell_j, cell_k), cell_getid(cdim, i_min, cell_j, cell_k),
                        cell_getid(cdim, cell_i, j_plus, cell_k), cell_getid(cdim, cell_i, j_min, cell_k),
                        cell_getid(cdim, cell_i, cell_j, k_plus), cell_getid(cdim, cell_i, cell_j, k_min)};

    for (int j=0; j<6; j++) {
      int pre_smoothing = 0;
      if (cells[i]->neighbours[j] == NULL && j!=which_neighbour) find_neighbour(level, cells[i]->parent, cells[i], search_id[j], fac, j, pre_smoothing, verbose);
    }
  }
}

int find_neighbour(struct AMR_levels *fine, struct cell *parent, struct cell *cell_link, int search_id, double fac, int which_neighbour, int pre_smoothing, int verbose) {
  int neighbour = 0; //Set to one if we do find a neighbour later
  int backward_neighbour = (which_neighbour%2 == 0) ? which_neighbour+1 : which_neighbour-1;
  int cdim[3] = {fine->cdim, fine->cdim, fine->cdim};
  if (verbose) message("Going to find neighbour %d for the cel at (%lf, %lf, %lf)", which_neighbour, cell_link->loc[0], cell_link->loc[1], cell_link->loc[2]);

  /* Check all the non-ghost cells. 
     A ghost cell does have a parent cell, but this parent is not split. But let's say it does have some (not eight) children, all ghosts!
     So parent->split indicates whether we are doing a neighbour search for a ghost or regular cell
     Should we do this, or does this mess up the red-black sweeps? */

  if (parent->split) {
    if (verbose) message("The parent of the cell with ghost value %d is split", cell_link->ghost);
    /* We can go to the parent's neighbour, because we already linked the children of the parent */
    struct cell *search_cand = parent->neighbours[which_neighbour];
    if (search_cand == NULL && pre_smoothing) {
      neighbour = 0;
      return neighbour;
    }
    //else message("Neighbour has children");
    /* Check progeny */
    for (int i=0; i<8; i++) {
      if (search_cand->progeny[i] == NULL) continue;
      int cell_i = (int) ((search_cand->progeny[i]->loc[0]+0.0001)/fac);
      int cell_j = (int) ((search_cand->progeny[i]->loc[1]+0.0001)/fac);
      int cell_k = (int) ((search_cand->progeny[i]->loc[2]+0.0001)/fac);
      int found_id = cell_getid(cdim, cell_i, cell_j, cell_k);

      if (found_id == search_id) {
        neighbour = 1;
        /* Link the cells if this not already happened */
        if (cell_link->neighbours[which_neighbour] == NULL) {
          //message("Creating link");
          cell_link->neighbours[which_neighbour] = search_cand->progeny[i];
          search_cand->progeny[i]->neighbours[backward_neighbour] = cell_link;
        }
        return neighbour;
      }
    }
  }

  /* Check parent cell and neighbour j of the parent */
  else {
    if (verbose) message("Checking progeny of the same parent which has loc (%lf, %lf, %lf)", parent->loc[0], parent->loc[1], parent->loc[2]);
    for (int i=0; i<8; i++) {
      if (parent->progeny[i] == NULL) continue;
      if (verbose) message("Child %d at location (%lf, %lf, %lf) of the same parent exists", i, parent->progeny[i]->loc[0], parent->progeny[i]->loc[1], parent->progeny[i]->loc[2]);
      int cell_i = (int) ((parent->progeny[i]->loc[0]+0.0001)/fac);
      int cell_j = (int) ((parent->progeny[i]->loc[1]+0.0001)/fac);
      int cell_k = (int) ((parent->progeny[i]->loc[2]+0.0001)/fac);
      int found_id = cell_getid(cdim, cell_i, cell_j, cell_k);

      if (found_id == search_id) {
        //message("Found one and going to link");
        neighbour = 1;
        /* Link the cells if this not already happened */
        if (cell_link->neighbours[which_neighbour] == NULL) {
          //message("Creating link");
          cell_link->neighbours[which_neighbour] = parent->progeny[i];
          parent->progeny[i]->neighbours[backward_neighbour] = cell_link;
        }
        return neighbour;
      }
    }
    struct cell *search_cand = parent->neighbours[which_neighbour];
    //message("Trying to search in neighbour");
    if (search_cand == NULL) return neighbour;
    if (verbose) message("Checking the parent neighbour at (%lf, %lf, %lf)", search_cand->loc[0], search_cand->loc[1], search_cand->loc[2]);
    for (int i=0; i<8; i++) {
      //message("Searching neighbour children");
      if (search_cand->progeny[i] == NULL) continue;
      if (verbose) message("Child %d at location (%lf, %lf, %lf) of the parent neighbour %d exists", i, search_cand->progeny[i]->loc[0], search_cand->progeny[i]->loc[1], search_cand->progeny[i]->loc[2], which_neighbour);
      int cell_i = (int) ((search_cand->progeny[i]->loc[0]+0.0001)/fac);
      int cell_j = (int) ((search_cand->progeny[i]->loc[1]+0.0001)/fac);
      int cell_k = (int) ((search_cand->progeny[i]->loc[2]+0.0001)/fac);
      int found_id = cell_getid(cdim, cell_i, cell_j, cell_k);

      if (found_id == search_id) {
        //message("Linking problem");
        neighbour = 1;
        /* Link the cells if this not already happened */
        if (cell_link->neighbours[which_neighbour] == NULL) {
          //message("Creating link");
          cell_link->neighbours[which_neighbour] = search_cand->progeny[i];
          search_cand->progeny[i]->neighbours[backward_neighbour] = cell_link;
        }
        return neighbour;
      }
    }
  }
  return neighbour;
}

double find_coarse_cell(struct AMR_levels *coarse, int cid, int missing, double fac_coarse, int cdimH[3]) {
  double potential=0;
  if (cid - missing<0) {
    for (int i=0;i<cid+1;i++) {
      int x_found = (int) ((coarse->cells[i]->loc[0]+0.0001)/fac_coarse);
      int y_found = (int) ((coarse->cells[i]->loc[1]+0.0001)/fac_coarse);
      int z_found = (int) ((coarse->cells[i]->loc[2]+0.0001)/fac_coarse);
      if ((int) cell_getid(cdimH, x_found, y_found, z_found) == cid) {
        message("We found a match for cid=%d at (%d,%d,%d)", cid, x_found, y_found, z_found);
        /* Return the potential in the cell we found a match for */
        potential = coarse->cells[i]->CIC_potential;
      }
    }
  }
  else{
    for (int i=cid-missing;i<cid+1;i++) {
      int x_found = (int) ((coarse->cells[i]->loc[0]+0.0001)/fac_coarse);
      int y_found = (int) ((coarse->cells[i]->loc[1]+0.0001)/fac_coarse);
      int z_found = (int) ((coarse->cells[i]->loc[2]+0.0001)/fac_coarse);
      if ((int) cell_getid(cdimH, x_found, y_found, z_found) == cid) {
        message("We found a match for cid=%d at (%d,%d,%d)", cid, x_found, y_found, z_found);
        /* Return the potential in the cell we found a match for */
        potential = coarse->cells[i]->CIC_potential;
      }
    }
  }
  return potential;
}

void extract_AMR_patches(struct space *s, struct cell *cells_top, int extract_depth, struct cell ***extracted_cells, int *count) {
  int current_depth = 1;
  for (int i=0; i<s->nr_cells; i++) {
    if (cells_top[i].maxdepth < extract_depth) continue;
    for (int j=0;j<8;j++) {
      if (cells_top[i].progeny[j] == NULL) continue;
      if (current_depth == extract_depth) {
        if (*count == 0) *extracted_cells = malloc(sizeof(struct cell *));
        else *extracted_cells = realloc(*extracted_cells, (*count+1)*sizeof(struct cell *));
        (*extracted_cells)[*count] = cells_top[i].progeny[j];
        *count +=1;
        //message("Assigning the cell at [%lf,%lf,%lf]", cells_top[i].progeny[j]->loc[0],cells_top[i].progeny[j]->loc[1],cells_top[i].progeny[j]->loc[2]);
      }
      else extract_lower(cells_top[i].progeny[j], &current_depth, extract_depth, extracted_cells, count);
    }
  }
}

void extract_lower(struct cell *parent, int *current_depth, int extract_depth, struct cell ***extracted_cells, int *count) {
  *current_depth += 1;
  for (int i=0;i<8;i++) {
    if (parent->progeny[i] == NULL) continue;
    if (*current_depth == extract_depth) {
      if (*count ==0) *extracted_cells = malloc(sizeof(struct cell *));
      else *extracted_cells = realloc(*extracted_cells, (*count+1)*sizeof(struct cell));
      (*extracted_cells)[*count] = parent->progeny[i];
      *count += 1;
      //message("Assigning the cell at [%lf,%lf,%lf]", parent->progeny[i]->loc[0],parent->progeny[i]->loc[1],parent->progeny[i]->loc[2]);
    }
    else extract_lower(parent->progeny[i], current_depth, extract_depth, extracted_cells, count);
  }
  *current_depth -= 1;
}

int perform_uniform_calculation(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]) {
  message("Going to do the uniform calculation");
  const int N_min = s->cdim[0];
  const double box_size = s->dim[0];
  if (min_depth < 0) {
    error("N_levels must be >= 0");
  }
  int array_size = min_depth+1;
  message("Array size is %d", array_size);

  /* Get array of the different grid sizes */
  int n = N_min;
  int level = 0;
  int grid_sizes[array_size];
  while (level< array_size) {
    grid_sizes[level++] = n;
    n *=2;
  }

  //struct cell_basic **cells_uniform[N_levels];
  double *rho[array_size];
  double *pot[array_size];
  double mean_density[array_size];

  for (int i = 0; i < array_size; i++) {
    rho[i] = NULL;
    pot[i] = NULL;
    mean_density[i] = 0.0;
  }

  /* Initialise arrays for density calculation. */
  for (int i=0; i<array_size; i++) {
    message("Doing i=%d",i);
    //cells_uniform[i] = calloc()
    rho[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    pot[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));

    if (!rho[i]) error("Uninitialised value");
    if (!pot[i]) error("Uninitialised value");

    /* Pass density information from cells to uniform grid */
    for (int j=0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++) {
      rho[i][j] = levels[i].cells[j]->CIC_density;
    }

    /*if (i==1) {
      struct cic_mapper_data data; 
      data.N = grid_sizes[i];
      data.rho = rho[i];
      data.fac = grid_sizes[i]/box_size;
      data.dim[0] = s->dim[0];
      data.dim[1] = s->dim[1];
      data.dim[2] = s->dim[2];
      int cdim[3] = {grid_sizes[i], grid_sizes[i], grid_sizes[i]};
      get_pm_potential(&data, grid_sizes[i], box_size, &s->e->threadpool, cdim);
    }*/
    message("Getting mean density for i=%d", i);
    mean_density[i] = get_mean_density(rho[i], grid_sizes[i]);
  }

  /* Set initial guess on the coarsest grid */
  int cdim[3] = {N_min, N_min, N_min};
  set_initial_guess(pot[0], cdim);

  apply_GS(rho[0], pot[0], cdim, mean_density[0], box_size);
  message("Got through level zero");

  /* Loop to get solutions on higher levels */
  for (int i = 1; i<array_size; i++) {
    message("Going to the level with grid size %d \n", grid_sizes[i]);
    prolongate_solution(pot[i-1], pot[i], grid_sizes[i-1], grid_sizes[i]); //After this function is completed, potential_array corresponds the solution prolongated to one level finer
    double potential_sum = 0;
    for (int j =0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++)
      potential_sum += fabs(pot[i][j]);
    message("The mean absolute value of the potential is %lf", potential_sum/(grid_sizes[i]*grid_sizes[i]*grid_sizes[i]));
    
    //Export the potential array to a .txt file to be read in Python
    //FILE *fptr;
    //fptr = fopen("/data1/vandervlugt/PythonFiles/smoothing_test/AMR_interpolated_potential_smooth_32.txt", "w");
    //for (int j=0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++) {
      //fprintf(fptr, "%lf\n", pot[i][j]);
    //}
    //fclose(fptr);
    //message("Done");
    //if (i==2) break;
    int N = grid_sizes[i];
    cdim[0] = N;
    cdim[1] = N;
    cdim[2] = N;
    //message("Going to apply multigrid");
    apply_multigrid(rho[i], pot[i], cdim, mean_density[i], box_size, N_min, N, 5);
  }

  /* Now pass the information calculated on the grid back to the cells */
  for (int i=0; i<array_size; i++) {
    for (int j=0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++) {
      levels[i].cells[j]->CIC_potential = pot[i][j];
    }
  }

  for (int i=0; i<array_size; i++) {
    free(rho[i]);
    free(pot[i]);
  }

  return grid_sizes[min_depth];
}

void potential_to_cells(struct cell *cells_top, double *pot, int grid_size, double grid_top, double fac, int level) {
  //double cdim[3] = {grid_size,grid_size,grid_size};
  double pot_assign;
  int current_level = 0;
  //We should be able to find the correct top-level cell in one go
  for (int i=0;i<grid_size*grid_size*grid_size;i++) {
    pot_assign = pot[i];
    double x_loc = (double) ((int) (i/(grid_size*grid_size))) * fac;
    double y_loc = (double) ((int) (i/grid_size) % grid_size) *fac;
    double z_loc = (double) (i % grid_size) * fac;
    //message("Assigning the potential %lf for location [%lf,%lf,%lf]", pot_assign, x_loc,y_loc,z_loc);
    double pot_cell[6] = {x_loc, x_loc + fac, y_loc, y_loc + fac, z_loc, z_loc + fac};
    for (int j=0;j<grid_top*grid_top*grid_top;j++) {
      double overlap = get_overlap(pot_cell, cells_top[j].loc, cells_top[j].width);
      if (overlap < 10e-2) continue;
      //message("At top level we found the overlap %lf for j=%d", overlap,j);
      to_lower_level(&cells_top[j], pot_cell, pot_assign, level, &current_level); //We never assign to the highest level, so immediately continue
    }
  }
}

void to_lower_level(struct cell *cell, double pot_cell[6], double pot, int level, int *current_level) {
  double overlap;
  *current_level += 1;
  for (int i=0;i<8;i++) {
    overlap = get_overlap(pot_cell,cell->progeny[i]->loc, cell->progeny[i]->width);
    if (overlap<10e-2) continue;
    if (*current_level == level) {
      //message("We found the overlap %lf for i=%d", overlap,i);
      cell->progeny[i]->CIC_potential = pot;
      //message("Assigned %lf at level %d", cell->progeny[i]->CIC_potential, *current_level);
    }
    else to_lower_level(cell->progeny[i], pot_cell, pot, level, current_level);
  }
  *current_level -= 1;
}

void sort_lower_level(struct cell *parent, double *rho, double fac, int max_level, int *current_level, int cdim[3]) {
  *current_level +=1;
  //message("Went to level %d", *current_level);
  for (int i=0;i<8;i++) { //Working with uniform grids, so all eight children must exist
    if (parent->progeny[i] == NULL) error("This cell must exist");
    if (*current_level == max_level) {
      int cell_x = (int) ((parent->progeny[i]->loc[0]+0.0001)/fac);
      int cell_y = (int) ((parent->progeny[i]->loc[1]+0.0001)/fac);
      int cell_z = (int) ((parent->progeny[i]->loc[2]+0.0001)/fac);
      size_t cid = cell_getid(cdim, cell_x, cell_y, cell_z);
      //message("Assigning cid [%d,%d,%d] for location [%lf,%lf,%lf]", cell_x, cell_y, cell_z, (parent->progeny[i]->loc[0]+0.1)/fac,(parent->progeny[i]->loc[1]+0.1)/fac, (parent->progeny[i]->loc[2]+0.1)/fac);
      if (rho[cid]!= 0.) error("We got the wrong location. Cell was already assigned");
      rho[cid] = parent->progeny[i]->CIC_density;
    }
    else {
      sort_lower_level(parent->progeny[i], rho, fac, max_level, current_level, cdim);
    }
  }
  *current_level -= 1;
}

void test_density_assignment(struct cell *cells_top, int nr_cells, const double boxsize, int nr_gparts, struct space *s, int level_check) {
  struct gpart *gparts = s->gparts;
  for (int i=0; i<nr_gparts;i++) {
    struct gpart* gp = &gparts[i];
    if (gp->x[0] > 137){
      message("Doing this function for [%lf,%lf,%lf]", gp->x[0],gp->x[1], gp->x[2]);
      get_CIC_results(gp, s->cells_top, s);

      initialise_tree_search(gp,&(cells_top[0]),cells_top,boxsize,nr_cells, level_check, 0);
    }
  }
}

void get_CIC_results(struct gpart *gp, struct cell *cells_top, struct space *s) {
  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int N = 2*s->cdim[0];
  const double box_size = s->dim[0];
  const double fac = N/box_size;
  const double cell_width = box_size/N;

  /* Box wrap the multipole's position */
  double pos_x = box_wrap(gp->x[0], 0., dim[0]);
  double pos_y = box_wrap(gp->x[1], 0., dim[1]);
  double pos_z = box_wrap(gp->x[2], 0., dim[2]);

  /* Workout the CIC coefficients */
  int i2 = (int)(fac * pos_x);
  if (i2 >= N) i2 = N - 1;
  double dx = fac * pos_x - i2;
  double tx = 1. - dx;

  int j = (int)(fac * pos_y);
  if (j >= N) j = N - 1;
  double dy = fac * pos_y - j;
  double ty = 1. - dy;

  int k = (int)(fac * pos_z);
  if (k >= N) k = N - 1;
  double dz = fac * pos_z - k;
  double tz = 1. - dz;

  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*tx * ty * tz, cell_width*(double)i2,cell_width*(double)j,cell_width*(double)k);
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*tx * ty * dz, cell_width*(double)i2,cell_width*(double)j,cell_width*(double)(k+1));
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*tx * dy * tz, cell_width*(double)i2,cell_width*(double)(j+1),cell_width*(double)k);
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*tx * dy * dz, cell_width*(double)i2,cell_width*(double)(j+1),cell_width*(double)(k+1));
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*dx * ty * tz, cell_width*(double)(i2+1),cell_width*(double)j,cell_width*(double)k);
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*dx * ty * dz, cell_width*(double)(i2+1),cell_width*(double)j,cell_width*(double)(k+1));
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*dx * dy * tz, cell_width*(double)(i2+1),cell_width*(double)(j+1),cell_width*(double)k);
  message("Assigning %lf to the cell at [%lf,%lf,%lf]", gp->mass*dx * dy * dz, cell_width*(double)(i2+1),cell_width*(double)(j+1),cell_width*(double)(k+1));

}

struct tree_data {
  struct cell* cells;
  int nr_cells;
  int level_check;
  double box_size;
};

void assign_densities(struct cell *cells_top, struct threadpool *tp, int nr_cells, const double boxsize, int level_check) {
  message("Assign densities called");
  for (int i=0;i<nr_cells;i++){
    cells_top[i].CIC_density = 0.;
  }

  int split;
  int level = 0;
  for (int i=0;i<nr_cells; ++i) {
    split = cells_top[i].split;
    //cells_top[i].CIC_density = 0.;
    if (split && level_check!=0) { 
      //message("Going to do this");
      for (int j=0;j<8;j++) {
        check_lower_level(cells_top[i].progeny[j], cells_top, boxsize, nr_cells, &level, level_check);
      }
    }
    else {
      //message("Initialising tree search at level %d", level);
      //message("At least one top level cell is not split");
      int nr_gparts = cells_top[i].grav.count;
      for (int j=0; j<nr_gparts; j++) {
        initialise_tree_search(&(cells_top[i].grav.parts[j]), &(cells_top[i]), cells_top, boxsize, nr_cells, level_check, 0);
        //for (int k=0; k<6; k++) {
          //if (tot_assigned[k]>41.564) error("Too large density assigned at level %d!", k);
        //}
      }
    }
  }

  //message("assigning the parent densities");
  //for (int i=0; i<nr_cells; i++) {
    //int level2 = 0;
    //assign_parent_densities(&cells_top[i], &level2);
    //message("The mass of top-level cell %d is %lf", i, cells_top[i].CIC_density);
  //}

  /* Export information at a certain level for debugging */
  //export_grid_data(cells_top, nr_cells);

  //message("Exporting AMR density data");
  //FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/particle_acc_test/AMR_pot_to_cells/density_test_AMR8_testtest", "w");
  //for (int j = 0; j < nr_cells; j++) {
    //fprintf(file_swift_density, "%lf\n", cells_top[j].CIC_density);
  //}
  //fclose(file_swift_density);
  //message("Exported the AMR density data.");
  //sleep(5);
}

void cells_top_mapper(void* map_data, const int num, void* extra) {
  message("Function gets called. We are processing %d cells", num);
  struct tree_data* data = (struct tree_data*)extra;

  /* Unpack the array */
  //struct cell *cells = data->cells;
  const int nr_cells = data->nr_cells;

  /* Unpack the conditions for the loop */
  const int level_check = data->level_check;
  const double boxsize = data->box_size;

  int split;
  int level = 0;

  /* Pointer to the chunk to be processed */
  struct cell* cells_top = (struct cell*)map_data;

  for (int i=0;i<num; ++i) {
    split = cells_top[i].split;
    //cells_top[i].CIC_density = 0.;
    if (split && level_check!=0) {
      for (int j=0;j<8;j++) {
        check_lower_level(cells_top[i].progeny[j], cells_top, boxsize, nr_cells, &level, level_check);
      }
    }
    else {
      //message("Initialising tree search at level %d", level);
      //message("At least one top level cell is not split");
      int nr_gparts = cells_top[i].grav.count;
      for (int j=0; j<nr_gparts; j++) {
        initialise_tree_search(&(cells_top[i].grav.parts[j]), &(cells_top[i]), cells_top, boxsize, nr_cells, level_check, 0);
      }
    }
  }
}

void export_grid_data(struct cell *cells_top, int nr_cells) {
  //message("Exporting information about grid level 2");
  FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/AMR_test/level_6/AMR_1", "w");

  for (int i = 0; i<nr_cells; i++) {
    if (cells_top[i].split) {
      for (int j=0;j<8;j++) {
        //if (cells_top[i].progeny[j]->split) {
          //for (int k=0;k<8;k++) fprintf(file_swift_density, "loc=[%lf,%lf, %lf]\t\t%lf\n", cells_top[i].progeny[j]->progeny[k]->loc[0], cells_top[i].progeny[j]->progeny[k]->loc[1], cells_top[i].progeny[j]->progeny[k]->loc[2],cells_top[i].progeny[j]->progeny[k]->CIC_density);
          //fprintf(file_swift_density, "loc=[%lf,%lf, %lf]\t\t%lf\n", cells_top[i].progeny[j]->loc[0], cells_top[i].progeny[j]->loc[1], cells_top[i].progeny[j]->loc[2],cells_top[i].progeny[j]->CIC_density);
          fprintf(file_swift_density, "%lf\n", cells_top[i].progeny[j]->CIC_density);
        //}
        //else message("Somehow, this is not split!");
      }
    }
  }
  fclose(file_swift_density);
  message("Exported the grid info");
  sleep(7);
}

void assign_parent_densities(struct cell *parent, int *level) {
  //message("function getting called");
  *level += 1;
  //message("We are on level %d", *level);
  int split = parent->split;
  if (!split) {
    //message("Doing this");
    *level -=1;
    return; //In case we have a leaf cell
  }
  for (int i=0;i<8; i++) {
    assign_parent_densities(parent->progeny[i], level);
  }
  double mass_sum=0.;
  for (int i=0;i<8;i++) {
    mass_sum += parent->progeny[i]->CIC_density;
  }
  //message("Assigning %lf", mass_sum);
  parent->CIC_density = mass_sum;
  *level -= 1;
}

void initialise_tree_search(struct gpart *part, struct cell *home_cell, struct cell *cells_top, const double boxsize, int nr_topcells, int level_check, int verbose) {
  //message("Initialising the tree search");
  double temp=0.;
  /* Box wrap the multipole's position */
  double pos_x = box_wrap(part->x[0], 0., boxsize);
  double pos_y = box_wrap(part->x[1], 0., boxsize);
  double pos_z = box_wrap(part->x[2], 0., boxsize);

  double x[3] = {pos_x, pos_y, pos_z};
  double width[3] = {home_cell->width[0], home_cell->width[1], home_cell->width[2]}; //Width of constant density cube by which we represent the particle
  //message("We are assigning the width %lf,%lf,%lf", width[0],width[1],width[2]);
  double pbox[6] = {x[0], x[0] + width[0], x[1],  x[1] + width[1], x[2], x[2] + width[2]};
  double pbox_shift[6];

  double x_shift[2] = {-boxsize, 0.};
  double y_shift[2] = {-boxsize, 0.};
  double z_shift[2] = {-boxsize, 0.};

  for (int i=0; i<2;i++) {
    for (int j=0; j<2;j++) {
      for (int k=0; k<2;k++) {
        pbox_shift[0] = pbox[0]+x_shift[i];
        pbox_shift[1] = pbox[1]+x_shift[i];
        pbox_shift[2] = pbox[2]+y_shift[j];
        pbox_shift[3] = pbox[3]+y_shift[j];
        pbox_shift[4] = pbox[4]+z_shift[k];
        pbox_shift[5] = pbox[5]+z_shift[k];
        if (pbox_shift[1]<0. || pbox_shift[3] <0. || pbox_shift[5]<0.) continue;
        //if (verbose) message("Initialising for i=%d, j=%d, k=%d", i, j, k);
        search_tree(part, pbox_shift, cells_top, width[0], nr_topcells, &temp, level_check, 0, verbose, 0);
      }
    }
  }
  //if (verbose) message("Exiting initialise search tree");
}

double get_overlap(double pbox[6], double cloc[3], double cwidth[3]) {
  double cbox[6] = {cloc[0], cloc[0]+cwidth[0], cloc[1], cloc[1]+cwidth[1], cloc[2], cloc[2]+cwidth[2]};
  double ox = fmax(0.0, fmin(pbox[1], cbox[1]) - fmax(pbox[0], cbox[0]));
  //if (fmin(pbox[1], cbox[1]) - fmax(pbox[0], cbox[0]) > 0) {
    //message("locations are pbox: [%lf, %lf, %lf, %lf, %lf, %lf] ", pbox[0],pbox[1], pbox[2],pbox[3],pbox[4],pbox[5]);
    //message("and cbox: [%lf, %lf, %lf, %lf, %lf, %lf]", cbox[0],cbox[1], cbox[2],cbox[3],cbox[4],cbox[5]);
    //message("ox =%lf", ox);
  //}
  //message("The y values are pbox = [%lf, %lf] and cbox = [%lf,%lf]", pbox[2],pbox[3],cbox[2],cbox[3]);
  double oy = fmax(0.0, fmin(pbox[3], cbox[3]) - fmax(pbox[2], cbox[2]));
  //if (fmin(pbox[3], cbox[3]) - fmax(pbox[2], cbox[2]) > 0) {
    //message("locations are pbox: [%lf, %lf] and cbox: [%lf, %lf]", pbox[2], pbox[3],cbox[2], cbox[3]);
    //message("oy =%lf", oy);
  //}
  double oz = fmax(0.0, fmin(pbox[5], cbox[5]) - fmax(pbox[4], cbox[4]));
  //if (fmin(pbox[5], cbox[5]) - fmax(pbox[4], cbox[4]) > 0) {
    //message("locations are pbox: [%lf, %lf] and cbox: [%lf, %lf]", pbox[4], pbox[5],cbox[4], cbox[5]);
    //message("oz =%lf", oz);
  //}
  //if (ox*oy*oz >0) message("Returning nonzero %lf", ox*oy*oz);
  return ox * oy * oz;
}

void search_tree(struct gpart *part, double pbox[6], struct cell *cells_top, double width, int nr_topcells, double *acc, int level_check, int reverse, int verbose, int calc_acc) {
  double overlap = 0.;
  int level =0;

  for (int i=0; i<nr_topcells;i++) {
    overlap = get_overlap(pbox,cells_top[i].loc, cells_top[i].width);
    if (overlap >0. && reverse && level_check == 3) message("Overlap was %lf between pbox with edge (%lf, %lf, %lf) and top cell at (%lf, %lf, %lf)", overlap, pbox[0], pbox[2], pbox[4], cells_top[i].loc[0], cells_top[i].loc[1], cells_top[i].loc[2]);
    if (overlap == 0.) continue;
    if (overlap<0.) error("Overlap smaller than zero detected..");
    if (!reverse) {
      //if (verbose) message("Assigning to top cell %d", i);
      atomic_add_d(&(cells_top[i].CIC_density), part->mass*overlap/(width*width*width));
    }
    //message("Succesfully added a value");
    //cells_top[i].CIC_density += part->mass*overlap/(width*width*width);
    if (cells_top[i].split && level_check!=0) {
    //if (cells_top[i].split && level <2) {
      /* Go a level lower */
      //if (reverse) message("Going to recurse");
      assign_lower_level(part, pbox, &cells_top[i], width, &level, acc, level_check, reverse, verbose, calc_acc);
    }
    else if (reverse) {
      //message("Doing this");
      //message("Just assigned %lf", cells_top[i].CIC_potential * overlap/(width*width*width));
      if (calc_acc) {
        acc[0] += cells_top[i].CIC_acc[0] * overlap/(width*width*width);
        acc[1] += cells_top[i].CIC_acc[1] * overlap/(width*width*width);
        acc[2] += cells_top[i].CIC_acc[2] * overlap/(width*width*width);
      }
      else *acc += cells_top[i].CIC_potential * overlap/(width*width*width);
    }
  }
  //if (verbose) message("Exiting search_tree");
}

void assign_lower_level(struct gpart *part, double pbox[6], struct cell *parent, double width, int *level, double *acc, int level_check, int reverse, int verbose, int calc_acc) {
  double overlap = 0.;
  *level +=1;
  for (int i=0;i<8;i++) {
    if (parent->progeny[i] == NULL) {
      error("This cell should have been created");
    }
    overlap = get_overlap(pbox, parent->progeny[i]->loc, parent->progeny[i]->width);
    if (overlap >0. && reverse && level_check == 3) message("Overlap was %lf between pbox with edge (%lf, %lf, %lf) and child at (%lf, %lf, %lf)", overlap, pbox[0], pbox[2], pbox[4], parent->progeny[i]->loc[0], parent->progeny[i]->loc[1], parent->progeny[i]->loc[2]);
    if (overlap == 0.) continue;
    if (overlap<0.) error("Overlap smaller than zero detected..");
    //if (*level>5) message("Assigning at level %d", *level);
    if (!reverse) {
      //if (verbose) message("Assigning to progeny %d at level %d", i, *level);
      atomic_add_d(&(parent->progeny[i]->CIC_density), part->mass * overlap/(width*width*width));
    }
    //parent->progeny[i]->CIC_density += part->mass * overlap/(width*width*width);
    //if (verbose) message("The split value is %d", parent->progeny[i]->split);
    if (parent->progeny[i]->split && *level<level_check) {
    //if (parent->progeny[i]->split && *level<3) {
      //message("The overlap was nonzero. Going to lower level");
      //if (reverse) message("Going to recurse further");
      assign_lower_level(part, pbox, parent->progeny[i], width, level, acc, level_check, reverse, verbose, calc_acc);
    }
    else if (reverse) {
      //message("Going to assign");
      if (calc_acc) {
        acc[0] += parent->progeny[i]->CIC_acc[0] * overlap/(width*width*width);
        acc[1] += parent->progeny[i]->CIC_acc[1] * overlap/(width*width*width);
        acc[2] += parent->progeny[i]->CIC_acc[2] * overlap/(width*width*width);
      }
      else *acc += parent->progeny[i]->CIC_potential * overlap/(width*width*width);
      if (level_check == 3) message("Just assigned %lf", parent->progeny[i]->CIC_potential * overlap/(width*width*width));
    }
    //else {
      //message("The overlap was nonzero. Assigning mass");
      //if (*level>1) message("Assigning at level %d", *level);
      //parent->progeny[i]->CIC_density += part->mass * overlap/(width*width*width);
      //message("Assigning %lf to the cell at location [%lf,%lf,%lf]", part->mass * overlap/(width*width*width), parent->progeny[i]->loc[0], parent->progeny[i]->loc[1], parent->progeny[i]->loc[2]);
    //}
  }
  *level -=1;
  //if (verbose) message("Going back to level %d", *level);
}

void check_lower_level(struct cell *parent, struct cell *cells_top, const double boxsize, int nr_topcells, int *level, int level_check) {
  //message("check_lower_level called");
  *level += 1;
  int verbose = 0;
  //if (*level>4) message("We are on level %d", *level);
  if (parent->split && *level <level_check) {
  //if (parent->split && *level<2) {
  for (int i=0;i<8;i++) {
      if (parent->progeny[i] == NULL) {
        error("On level %d this cell %d should exist", *level, i);
      }
      else if (parent->progeny[i]->hydro.count == 0 && parent->progeny[i]->grav.count == 0 && parent->progeny[i]->stars.count == 0 &&
                parent->progeny[i]->black_holes.count == 0 && parent->progeny[i]->sinks.count == 0) {
        //message("Found fake daughter");
        /* We continue as there are by definition no particles in this fake daughter cell */
        continue;
      }
      else check_lower_level(parent->progeny[i], cells_top, boxsize, nr_topcells, level, level_check);
    }
  }
  else { 
    //if (*level == 4) message("Initialising tree search at level %d", *level);
    //if (*level>5) message("Finding particles at level %d", *level);
    int nr_gparts = parent->grav.count;
    for (int j=0;j<nr_gparts;j++) {
      //double tot_assigned[6] = {0.};
      if (*level == 4) {
        verbose = 1;
        //message("Doing particle number %d", j);
      }
      initialise_tree_search(&(parent->grav.parts[j]), parent, cells_top, boxsize, nr_topcells, level_check, verbose);
      //for (int k=0; k<6; k++) {
        //if (tot_assigned[k]>41.564) error("Too large density assigned at level %d!", k);
      //}
    }
  }
  *level -= 1;
  //message("Going back to level %d", *level);
}

void destroy_daughters(struct space *s, struct cell *cells_top, int nr_cells) {
  int split;
  int daughter_list = 0;
  for (int i = 0; i<nr_cells; i++) {
    split = cells_top[i].split;
    if (split) {
      level_down(s, &(cells_top[i]), &daughter_list);
    }
  }
  message("We destroyed %d daughter cells", daughter_list);
}

void level_down(struct space *s, struct cell *parent, int *length) {
  for (int i= 0; i<8; i++) {
    if (parent->split){
      if (parent->progeny[i] == NULL) error("This cell should exist");
      else if (parent->progeny[i]->hydro.count == 0 && parent->progeny[i]->grav.count == 0 && parent->progeny[i]->stars.count == 0 &&
                parent->progeny[i]->black_holes.count == 0 && parent->progeny[i]->sinks.count == 0) {
        space_recycle(s, parent->progeny[i]);
        parent->progeny[i] = NULL;
        //message("Destroyed a daughter cell");
        *length += 1;
      }
      else {
        level_down(s, parent->progeny[i], length);
      }
    }
  }
}

void check_progeny(struct space *s, struct cell *parent, int *length, int *level) {
  *level += 1;
  //message("On level %d", *level);
  struct cell *child = NULL;
  for (int i = 0; i<8; i++) {
    child = parent->progeny[i];
    if (child == NULL) {
      space_getcells(s, 1, parent->progeny, 0, i, 0);
      parent->progeny[i]->loc[0] = parent->loc[0];
      parent->progeny[i]->loc[1] = parent->loc[1];
      parent->progeny[i]->loc[2] = parent->loc[2];

      double width = parent->width[0]/2;
      if (i==1 || i==3 || i==5 || i==7) parent->progeny[i]->loc[2] += width;
      if (i==2 || i==3 || i==6 || i==7) parent->progeny[i]->loc[1] += width;
      if (i==4 || i==5 || i==6 || i==7) parent->progeny[i]->loc[0] += width;

      parent->progeny[i]->width[0] = parent->width[0]/2;
      parent->progeny[i]->width[1] = parent->width[1]/2;
      parent->progeny[i]->width[2] = parent->width[2]/2;

      parent->progeny[i]->parent = parent;

      //for (int j= 0; j<8; j++)
        //parent->progeny[i]->progeny[j] = NULL;
      
      //parent->progeny[i]->parent = parent;
      //parent->progeny[i]->top = parent->top;
      //parent->progeny[i]->CIC_density = 0.;
      //parent->progeny[i]->grav.count = 0.;
      //parent->progeny[i]->grav.parts = NULL;

      /* Flag that we need to destroy this cell later */
      //parent->progeny[i]->fake_daughter = 1;
      //parent->progeny[i]->split = 0;
      
      /* Also add to the list of AMR cells */
      *length +=1;
      //construct_daughter(parent, i, length);
      if (parent->progeny[i] == NULL) message("The child is still the null pointer...");
    }
    else if (child->split) {
      if (child->width[0] != child->width[1] || 
        child->width[0] != child->width[2]){
          message("The cell dimensions are not equal!");
        }
      child->CIC_density =0.;
      check_progeny(s, child, length, level);
      *level -= 1;
      //message("Back to level %d", *level);
    }
    else {
      //*length += 1;
      child->CIC_density =0.;
      //message("Added one to the list. Length is now %d. And the i-index is %d", *length, i);
    }
  }
}

void construct_daughter(struct cell *parent, int i, int *length) {

if (swift_memalign("extra.daughter", (void **)&parent->progeny[i], cell_align,
                      sizeof(struct cell)) != 0)
    error("Failed to allocate daughter cell.");
  bzero(parent->progeny[i], sizeof(struct cell));

  /*! The cell location on the grid (corner nearest to the origin). */
  //for (int j =0; j<i; j++) {
    //message("The location of cell %d was %lf, %lf, and %lf", j, parent->progeny[j]->loc[0], parent->progeny[j]->loc[1], parent->progeny[j]->loc[2]);
  //}
  //message("The location of the parent cell was: %lf, %lf, and %lf", parent->loc[0], parent->loc[1], parent->loc[2]);
  parent->progeny[i]->loc[0] = parent->loc[0];
  parent->progeny[i]->loc[1] = parent->loc[1];
  parent->progeny[i]->loc[2] = parent->loc[2];

  double width = parent->width[0]/2;
  if (i==1 || i==3 || i==5 || i==7) parent->progeny[i]->loc[2] += width;
  if (i==2 || i==3 || i==6 || i==7) parent->progeny[i]->loc[1] += width;
  if (i==4 || i==5 || i==6 || i==7) parent->progeny[i]->loc[0] += width;

  parent->progeny[i]->width[0] = parent->width[0]/2;
  parent->progeny[i]->width[1] = parent->width[1]/2;
  parent->progeny[i]->width[2] = parent->width[2]/2;

  for (int j= 0; j<8; j++)
    parent->progeny[i]->progeny[j] = NULL;
  
  parent->progeny[i]->parent = parent;
  parent->progeny[i]->top = parent->top;
  parent->progeny[i]->CIC_density = 0.;
  parent->progeny[i]->grav.count = 0.;
  parent->progeny[i]->grav.parts = NULL;

  /* Flag that we need to destroy this cell later */
  //parent->progeny[i]->fake_daughter = 1;
  parent->progeny[i]->split = 0;
  
  /* Also add to the list of AMR cells */
  *length +=1;
  //message("Added a NULL daughter to the list. Length is now %d. And the i-index is %d", *length, i);
}

void space_apply_FMG(struct space *s, struct engine *e) {
  const int N_max = 64;
  const int N_min = 16;
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int cdim_max[3] = {N_max, N_max, N_max};

  if (dim[0] != dim[1] || dim[0] != dim[2]) {
    message("Noncubic domain for density calculation. Did not fill density array");
    return;
  }

  /* Get array of the different grid sizes */
  int n = N_min;
  int level = 0;
  int grid_sizes[32];
  while (n<=N_max) {
    grid_sizes[level++] = n;
    n *=2;
  }

  const int N_levels = level;
  double *rho[N_levels];
  double *pot[N_levels];

  double mean_density[N_levels];

  struct cic_mapper_data data;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];

  /* Initialise arrays for density calculation. */
  for (int i=0; i<N_levels; i++) {
    rho[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    pot[i] = calloc(grid_sizes[i]*grid_sizes[i]*grid_sizes[i], sizeof(double));
    
    /* Assign densities to all arrays */ 
    data.N = grid_sizes[i];
    //message("The number of cells is %d", data.N);
    data.rho = rho[i];
    data.fac = grid_sizes[i]/box_size;
    gpart_to_mesh_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);

    /* Get pm potential for verification of FMG method */
    if (i == N_levels-1) {
      get_pm_potential(&data, N_max, box_size, &e->threadpool, cdim_max);
    }

    /* Get mean density for every level. Multiplying by the right conversion factor happens later. */
    //if (i==N_levels-1) {
      //message("Exporting FMG density data");
      //FILE *file_swift_density = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/FMG_check/FMG_density_64.txt", "w");
      //for (int j = 0; j < N_max*N_max*N_max; j++) {
        //fprintf(file_swift_density, "%lf\n", rho[i][j]);
      //}
      //fclose(file_swift_density);
      //message("Exported the FMG density data.");
      //sleep(5);
    //}
    mean_density[i] = get_mean_density(data.rho, grid_sizes[i]);
  }

  if (N_min%2 != 0){
    message("Uneven domain for density calculation. Did not apply Gauss-Seidel");
    return;
  }  
  
  /* Set initial guess on the coarsest grid */
  int cdim[3] = {N_min, N_min, N_min};
  set_initial_guess(pot[0], cdim);

  apply_GS(rho[0], pot[0], cdim, mean_density[0], box_size);

  /* Loop to get solutions on higher levels */
  for (int i = 1; i<N_levels; i++) {
    message("Going to the level with grid size %d \n", grid_sizes[i]);
    prolongate_solution(pot[i-1], pot[i], grid_sizes[i-1], grid_sizes[i]); //After this function is completed, potential_array corresponds the solution prolongated to one level finer
    //Export the potential array to a .txt file to be read in Python
    //FILE *fptr;
    //fptr = fopen("/data1/vandervlugt/PythonFiles/smoothing_test/FMG_interpolated_potential_32.txt", "w");
    //for (int j=0; j<grid_sizes[i]*grid_sizes[i]*grid_sizes[i]; j++) {
      //fprintf(fptr, "%lf\n", pot[i][j]);
    //}
    //fclose(fptr);
    
    int N = grid_sizes[i];
    cdim[0] = N;
    cdim[1] = N;
    cdim[2] = N;
    apply_multigrid(rho[i], pot[i], cdim, mean_density[i], box_size, N_min, N, 3);
  }

  FILE *fptr;
  fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/FMG_check/FMG_64_cells.txt", "w");
  for (int j=0; j<grid_sizes[N_levels-1]*grid_sizes[N_levels-1]*grid_sizes[N_levels-1]; j++) {
    fprintf(fptr, "%lf\n", pot[N_levels-1][j]);
  }
  fclose(fptr);

  int step = e->step;

  /* Now calculate the accelerations from the potential */
  get_accelerations(&data, pot[N_levels-1], &e->threadpool, s, N_max, step);

  for (int i=0; i<N_levels; i++) {
    free(rho[i]);
    free(pot[i]);
  }

}

void get_accelerations(struct cic_mapper_data* data, double *pot, struct threadpool* tp, struct space *s,const int N, int step) {
  double box_size = s->dim[0];

  /* Gather the mesh shared information to be used by the threads */
  //data.cells = s->cells_top; //Liever niet
  data->rho = NULL;
  data->potential = pot;
  data->fac = N / box_size;
  data->dim[0] = s->dim[0];
  data->dim[1] = s->dim[1];
  data->dim[2] = s->dim[2];
  data->const_G = s->e->physical_constants->const_newton_G;

  mesh_to_gpart_CIC_mapper(s->gparts, s->nr_gparts, (void*)data);

  //threadpool_map(tp, mesh_to_gpart_CIC_mapper, s->gparts, s->nr_gparts,
                   //sizeof(struct gpart), threadpool_auto_chunk_size,
                   //(void*)data);

  /* Als het goed is hebben we nu de accelerations! */
  /* Stored in gp->a_grav_mesh[0] *= const_G; */

  //struct gpart* gp = &(s->gparts)[0];
  //message("We have a_grav_mesh = %lf", gp->a_grav_mesh[0]);

  int num = s->nr_gparts;
  struct gpart* gp;
  //double acc;

  message("Exporting FMG acceleration data");
  FILE *file_swift_acc = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/FMG_check/FMG_64_particles", "w");
  for (int j = 0; j < num; j++) {
    gp = &(s->gparts)[j];
    //acc = gp->a_grav_mesh[0];
    fprintf(file_swift_acc, "%lf %.15g %.15g %.15g \n", gp->potential_mesh, gp->x[0], gp->x[1], gp->x[2]);
  }
  fclose(file_swift_acc);
  message("Exported the FMG acceleration data.");
  sleep(5);

  //message("Exporting FMG data");
  //FILE *file_swift_acc = fopen("/data1/vandervlugt/PythonFiles/acceleration_test/FMG_accx_new", "w");
  //for (int i=0; i<num; i++) {
      //gp = &(s->gparts)[i];
      //acc = gp->a_grav_mesh[0];
      //fprintf(file_swift_acc, "%lf\n", acc);
    //}
  //fclose(file_swift_acc);

  //for (int i = 0; i<num; i++) {
    //gp = &(s->gparts)[i];
    //gp->a_grav_mesh[0] = 0.;
    //gp->a_grav_mesh[1] = 0.;
    //gp->a_grav_mesh[2] = 0.;

    //FILE *fptry;
    //fptry = fopen("/data1/vandervlugt/PythonFiles/acceleration_test/FMG_y.txt", "w");
    //for (int i=0; i<num; i++) {
      //gp = &(s->gparts)[i];
      //acc = gp->a_grav_mesh[1];
      //fprintf(fptry, "%lf\n", acc);
    //}
    //fclose(fptry);

    //FILE *fptrz;
    //fptrz = fopen("/data1/vandervlugt/PythonFiles/acceleration_test/FMG_z.txt", "w");
    //for (int i=0; i<num; i++) {
      //gp = &(s->gparts)[i];
      //acc = gp->a_grav_mesh[2];
      //fprintf(fptrz, "%lf\n", acc);
    //}
    //fclose(fptrz);
  //}
}

void prolongate_solution(double *pot_coarse, double *pot_fine, const int N, const int N_double) {

  int cdimh[3] = {N_double, N_double, N_double};
  int cdim[3] = {N,N,N};

  for (int i = 0; i<N_double; i++) {
    for (int j = 0; j<N_double; j++) {
      for (int k=0; k<N_double; k++) {
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

void space_get_density(struct space *s, struct swift_params *params, struct engine *e, int use_multigrid) {
  const int N = 128;
  int cdim[3] = {N,N,N};
  const double box_size = s->dim[0];
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  if (dim[0] != dim[1] || dim[0] != dim[2]) {
    message("Noncubic domain for density calculation. Did not fill density array");
    return;
  }

  double *density_array = NULL;
  /* Allocate the memory for the density array */
  density_array = (double*)calloc(N * N * N, sizeof(double));
  if (density_array == NULL)
    error("Error allocating memory for the density mesh.");
  memuse_log_allocation("mesh.density", density_array, 1,
                        sizeof(double) * N * N * N);
                        
  struct cic_mapper_data data; 
  message("The number of cells is %d", N);
  data.N = N;
  data.rho = density_array;
  data.fac = N/box_size;
  data.dim[0] = dim[0];
  data.dim[1] = dim[1];
  data.dim[2] = dim[2];

  gpart_to_mesh_CIC_mapper(s->gparts, s->nr_gparts, (void*)&data);
  //get_pm_potential(&data, N, box_size, &e->threadpool, cdim);  //PM method to get the potential is implemented here

  //double spacing = box_size/N;

  //message("Exporting regular density data");
  //FILE *file_regular_density = fopen("/data1/vandervlugt/PythonFiles/particle_acc_test/PM_pot_to_cells/density_CIC_8", "w");

  //for (int i=0; i<N*N*N; i++) {
    //fprintf(file_regular_density, "%lf\n", data.rho[i]);
  //}
  
  /*
  for (int i=0;i<N;i++) {
    for (int j=0; j<N;j++) {
      for (int k=0;k<N;k++) {
        int cid = cell_getid(cdim,i,j,k);
        double loc_x = spacing*(double)i;
        if (i==4) message("calculated locx=%lf", loc_x);
        double loc_y = spacing*(double)j;
        double loc_z = spacing*(double)k;
        fprintf(file_regular_density, "loc=[%lf,%lf, %lf]\t\t%lf\n", loc_x,loc_y,loc_z,data.rho[cid]);
      }
    }
  }
  */
  //fclose(file_regular_density);
  //message("Exported the regular density data.");
  //sleep(5);

  //if (N == s->cdim[0]){
    //density_to_cells(&data, s->cells_top, s->nr_cells, s->cdim); 
  //}
  
  double *potential_array = NULL;
  //Allocate the memory for the potential array 
  potential_array = (double*)malloc(sizeof(double) * N * N * N);
  if (potential_array == NULL){
    error("Error allocating memory for the potential mesh.");
  }
  memuse_log_allocation("mesh.potential", potential_array, 1, sizeof(double)*N*N*N);
  for (int i = 0; i<N*N*N; i++)
      potential_array[i] = 0.0;

  if (N%2 != 0){
    message("Uneven domain for density calculation. Did not apply Gauss-Seidel");
    return;
  }
 
  double mean_density = get_mean_density(density_array, N);
  double sum = 0.0;
  double rms = 0.0;
  for (int i = 0; i < N*N*N; i++) {
      rms += fabs(density_array[i]-1);
      sum += density_array[i]-1;
  }
  set_initial_guess(potential_array, cdim);

  if (use_multigrid) {
    apply_multigrid(density_array, potential_array, cdim, mean_density, box_size, 16, 128, 5);
  } else {
    apply_GS(density_array, potential_array, cdim, mean_density, box_size);
  }

  message("Exporting PM FFT data");
  FILE *fptr;
  fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/FMG_check/MG_128_cells.txt", "w");
  for (int i=0; i<N*N*N; i++) {
    fprintf(fptr, "%lf\n", potential_array[i]);
  }
  fclose(fptr); 

  free(density_array);
  free(potential_array);
}

void apply_multigrid(double *density, double *pot, int cdim[3], double mean_density, double box_size, int N_stop, int N_start, int V_max) {
  //message("Applying the multigrid method...");
  int N = cdim[0];
  double *residual_array = NULL;
  /* Allocate the memory for the potential array */
  residual_array = (double*)malloc(sizeof(double) * N * N * N);
  if (residual_array == NULL){
    error("Error allocating memory for the residual array.");
  }
  memuse_log_allocation("residual.array", residual_array, 1, sizeof(double)*N*N*N);

  double delta = box_size/N;
  double residual; 
  double fac2 = N/box_size;
  mean_density *= 4.*M_PI*fac2*fac2*fac2;
  get_residual(pot, density, cdim, mean_density, &residual, delta);
  message("The first residual is %lf", residual);
  double tolerance = 10e-6; //Choose reasonable value here
  int counter = 0;
  //const int fine_steps = 10; //Arbitrary for now
  double sum = 0;

  int V_cycles=0;
  V_max = 3;
  //const int N_stop = 16; //Dimension of coarsest grid, on which to solve the equation exactly. 
  //const int N_start = 128; //Dimension of the finest grid

  while (residual > tolerance && V_cycles < V_max) {
    //message("Performing V-cycle %d", V_cycles);
    /* Pre-smoothing for 10 steps */
    for (int i2=0; i2<10; i2++) {
      //message("Performing red-black sweep %d", i);
      perform_red_black_sweep(pot, density, cdim, mean_density, delta); 
      /*if (V_cycles == 0 && i2==0) {
        for (int col=0; col<2; col++){
          for (int k=0; k<cdim[2]; k++){
            for (int j=0; j<cdim[1]; j++){
              for (int i=0; i<cdim[0]; i++) {
                if ((i+j+k)%2 != col) continue;
                if (cell_getid(cdim,i,j,k)<50) message("Set cell %lu with loc (%lf, %lf, %lf) to %lf", cell_getid(cdim,i,j,k), (double) i * delta, (double) j * delta, (double) k * delta,pot[cell_getid(cdim,i,j,k)]);
              }
            }
          }
        }
      }*/
    }
    //if (V_cycles == 5) {
      //message("Exporting PM FFT data");
      //FILE *fptr;
      //fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/uniform_check/FMG_16_middle_result.txt", "w");
      //for (int i=0; i<cdim[0]*cdim[0]*cdim[0]; i++) {
        //fprintf(fptr, "%lf\n", pot[i]);
      //}
      //fclose(fptr); 
      //message("Exported post pre-smoothing");
    //}
    get_residual(pot, density, cdim, mean_density, &residual, delta);
    //message("After 10 pre-smoothing steps and %d V-cycles, the residual of the 'normal' Gauss-Seidel is %lf", V_cycles, residual);
    get_residual_array(pot, density, cdim, mean_density, residual_array, delta);
    /* Get coarse-grid correction */
    solve_coarser_problem(pot, residual_array, cdim, delta, N_stop, N_start);
    get_residual(pot, density, cdim, mean_density, &residual, delta);
    /* Post-smoothing for 10 steps */
    for (int i=0; i<10; i++) {
        perform_red_black_sweep(pot, density, cdim, mean_density, delta); 
    }
    get_residual(pot, density, cdim, mean_density, &residual, delta);
    //message("After 10 post-smoothing steps and %d V-cycles, the residual of the 'normal' Gauss-Seidel is %lf", V_cycles, residual);    
    V_cycles +=1;
  }

  //message("The total number of V-cycles was %d", V_cycles);

  //sum = 0.0;
  //for (int i=0; i<N*N*N; i++) {
    //sum += residual_array[i];
  //}
  //message("The first index of the original residual array is %lf", residual_array[0]);
  
  /* Post-smoothing until convergence */
  while (residual > tolerance) {
    perform_red_black_sweep(pot, density, cdim, mean_density, delta);
    get_residual(pot, density, cdim, mean_density, &residual, delta);
    counter +=1;
    if (counter % 100 ==0){
      sum = 0;
      for (int i = 0; i < N*N*N; i++) {
          sum += pot[i];
      }
      //message("Performed %d steps and the residual is %lf", counter, residual);
    }
  }
  get_residual(pot, density, cdim, mean_density, &residual, delta);
  message("The total number of steps on the finest grid is %d and the residual is %lf", counter, residual);
  free(residual_array);

  //Set mean of potential to zero
  sum = 0.;
  for (int i = 0; i < N*N*N; i++) {
      sum += pot[i];
  }
  for (int i = 0; i<N*N*N; i++) {
    pot[i] -= sum/(N*N*N); 
  }

  double potential_sum = 0.;
  for (int i =0; i<N*N*N; i++)
    potential_sum += fabs(pot[i]);
  message("The mean absolute value of the potential is %lf", potential_sum/(N*N*N));

  //if (N==16) {
  //Export the potential array to a .txt file to be read in Python
    //FILE *fptr;
    //fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/uniform_check/FMG_16_check2.txt", "w");
    //for (int i=0; i<N*N*N; i++) {
      //fprintf(fptr, "%lf\n", pot[i]);
    //}
    //fprintf(fptr, "%lf\n", 50.);
    //fclose(fptr);
  //}
}

void apply_GS(double *density, double *pot, int cdim[3], double mean_density, double box_size) {
  message("Applying Gauss-Seidel...");
  //message("The mean of the density is %lf", mean_density);
  int N = cdim[0];
  double delta = box_size/N;
  double residual; 
  mean_density *= 4.*M_PI*N*N*N/(box_size*box_size*box_size);
  get_residual(pot, density, cdim, mean_density, &residual, delta);
  double tolerance = 10e-6; //Choose reasonable value here
  int counter = 0;
  double sum = 0;

  //sum = 0;
  //for (int i = 0; i < N*N*N; i++) {
      //sum += pot[i];
  //}
  //message("Performed %d steps and the residual is %lf", counter, residual);
  //message("The mean of the potential is %lf", sum/(N*N*N));

  while (residual >= tolerance) {
    perform_red_black_sweep(pot, density, cdim, mean_density, delta);
    get_residual(pot, density, cdim, mean_density, &residual, delta);
    counter +=1;
    if (counter % 50 ==0){
      sum = 0;
      for (int i = 0; i < N*N*N; i++) {
          sum += pot[i];
      }
      //message("Performed %d steps and the residual is %lf", counter, residual);
      //message("The mean of the potential is %lf", sum/(N*N*N));
    }
  }
  message("We needed %d steps to converge", counter);

  //Set mean of potential to zero
  sum = 0;
  for (int i = 0; i < N*N*N; i++) {
      sum += pot[i];
  }
  for (int i = 0; i<N*N*N; i++) {
    pot[i] -= sum/(N*N*N);
  }

  double potential_sum = 0;
  for (int i =0; i<N*N*N; i++)
    potential_sum += fabs(pot[i]);
  //message("The mean absolute value of the potential is %lf", potential_sum/(N*N*N));

  //get_residual(pot, density, cdim, mean_density, &residual, delta);
  //message("The final residual is %lf", residual);

  //Export the potential array to a .txt file to be read in Python
  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/density_test_plots/GS_64.txt", "w");
  //for (int i=0; i<N*N*N; i++) {
    //fprintf(fptr, "%lf\n", pot[i]);
  //}
  //fclose(fptr);
  
}

void solve_coarser_problem(double *pot, double *residual, int cdim[3], double delta, const int N_stop, const int N_start) {
  //message("The value of delta is %lf", delta);
  int N = cdim[0];
  delta = delta*2.0;
  double *coarser_residual = NULL;
  /* Allocate the memory for the potential array */
  coarser_residual = (double*)malloc(sizeof(double) * N * N * N);
  if (coarser_residual == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.grid", coarser_residual, 1, sizeof(double)*N*N*N);

  double *coarser_solution = NULL; //The coarser solution is the error obtained from solving on the coarser array. So we suppose it is zero. 
  /* Allocate the memory for the potential array */
  coarser_solution = (double*)malloc(sizeof(double) * N * N * N);
  if (coarser_solution == NULL){
    error("Error allocating memory for the coarser-grid array.");
  }
  memuse_log_allocation("coarser.solution", coarser_solution, 1, sizeof(double)*N*N*N);
  for (int i =0; i<N*N*N; i++)
    coarser_solution[i] = 0.0;

  double *new_residual_array = NULL; //The coarser solution is the error obtained from solving on the coarser array. So we suppose it is zero. 
  /* Allocate the memory for the potential array */
  new_residual_array = (double*)malloc(sizeof(double) * N * N * N);
  if (new_residual_array == NULL){
    error("Error allocating memory for the new coarser residual.");
  }
  memuse_log_allocation("coarser.newresidual", new_residual_array, 1, sizeof(double)*N*N*N);

  double sum = 0.0;
  //for (int i=0; i<N*N*N; i++) {
    //sum += residual[i];
  //}
  //message("The sum of the original residual array is %lf", sum);

  restrict_problem(coarser_residual, residual, cdim, delta); //Transfer residual to the coarser grid
  N = N/2;
  int cdimH[] = {cdim[0]/2, cdim[1]/2, cdim[2]/2}; 

  double multiplier = 1.0; 
  double coarser_residual_abs = 1.0; //Zodat we het iig 1x doen (testen!) 
  double tolerance = 10e-4; //Beetje random
  int counter = 0;

  //Multiply the residual array on H by -1 in order to solve the right equation
  for (int i =0; i<N*N*N; i++) 
    coarser_residual[i] = -1.0*coarser_residual[i];
  //sum = 0.0;
  //for (int i =0; i<N*N*N; i++)
    //sum += coarser_residual[i];
  //message("The sum of the coarser residual array is %lf", sum);

  double rhs_check = 0.;

  if (N==N_stop) {
    //message("We are now on the coarsest grid. Solving the residual equation exactly.");
    while (coarser_residual_abs >= tolerance) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
      get_residual(coarser_solution, coarser_residual, cdimH, multiplier, &coarser_residual_abs, delta);
      counter +=1;
      if (counter % 20 ==0){
        sum = 0.;
        rhs_check = 0.;
        for (int i = 0; i < N*N*N; i++) {
            sum += fabs(coarser_solution[i]);
            rhs_check += coarser_residual[i];
        }
        //message("Performed %d steps and the residual is %lf", counter, coarser_residual_abs);
        //message("The rms of the approximated solution is %lf", sum/(N*N*N));
        //message("The mean of the rhs is %lf", rhs_check);
      }
    }
    get_residual(coarser_solution, coarser_residual, cdimH, multiplier, &coarser_residual_abs, delta);
    message("The total number of steps on the coarsest grid is %d and the residual is %lf", counter, coarser_residual_abs);
    prolongate_problem(coarser_solution, pot, cdimH);
  }

  else { 
    //message("Not yet on the coarsest grid. The current dimension is %d", N);
    counter = 0;
    for (int i=0; i<50; i++) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
      get_residual(coarser_solution, coarser_residual, cdimH, multiplier, &coarser_residual_abs, delta);
      counter +=1;
      if (counter%10==0) {
        sum = 0.;
        rhs_check = 0.;
        for (int j = 0; j < N*N*N; j++) {
            sum += coarser_solution[j];
            rhs_check += coarser_residual[j];
        }
        //message("Performed %d steps and the residual is %lf", counter, coarser_residual_abs);
        //message("The rms of the approximated solution is %lf", sum/(N*N*N));
        //message("The mean of the rhs is %lf", rhs_check);
      }
    }
    //message("We are passing the values cdim=%d, mean_density=%lf, delta=%lf", cdim[0], mean_density, delta);
    get_residual(coarser_solution, coarser_residual, cdimH, multiplier, &coarser_residual_abs, delta);
    //message("After 10 steps, the residual of the 'normal' Gauss-Seidel is %lf", residual);
    get_residual_array(coarser_solution, coarser_residual, cdimH, multiplier, new_residual_array, delta); //Now we have a residual array instead of just a number

    /* Get coarse-grid correction */
    solve_coarser_problem(coarser_solution, new_residual_array, cdimH, delta, N_stop, N_start);

    /* Post-smoothing */
    for (int i=0; i<50; i++) {
      perform_red_black_sweep(coarser_solution, coarser_residual, cdimH, multiplier, delta); 
    }
    prolongate_problem(coarser_solution, pot, cdimH);
  }

  /* The coarser-grid correction has now been added to the finer-grid solution for the potential, so discard used arrays. */
  free(coarser_residual);
  free(coarser_solution);
  //message("Working back up again. The dimension is now %d", N);
}

void prolongate_problem(double *coarser_solution, double *pot, int cdim[3]) {
  int cdimh[] = {cdim[0]*2, cdim[1]*2, cdim[2]*2};
  for (int i=0; i<cdim[0]; i++) {
    for (int j=0; j<cdim[0]; j++) {
      for (int k=0; k<cdim[0]; k++) {
        const size_t cid = cell_getid(cdim,i,j,k);
        pot[cell_getid(cdimh,2*i,2*j,2*k)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i+1,2*j,2*k)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i,2*j+1,2*k)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i,2*j,2*k+1)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i+1,2*j+1,2*k)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i+1,2*j,2*k+1)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i,2*j+1,2*k+1)] += coarser_solution[cid];
        pot[cell_getid(cdimh,2*i+1,2*j+1,2*k+1)] += coarser_solution[cid];
      }
    }
  }
  /* Prepare for the next layer. */
  cdim[0] *= 2.;
  cdim[1] *= 2.;
  cdim[2] *= 2.;
}

void restrict_problem(double *H_array, double *residual, int cdim[3], double delta) {
  int cdimH[] = {cdim[0]/2, cdim[1]/2, cdim[2]/2};
  for (int i=0; i<cdimH[0]; i++) {
    for (int j=0; j<cdimH[0]; j++) {
      for (int k=0; k<cdimH[0]; k++) {
        //Might be useful to have this information in a tree or something
        H_array[cell_getid(cdimH,i,j,k)] = (residual[cell_getid(cdim,2*i,2*j,2*k)] + residual[cell_getid(cdim,2*i+1,2*j,2*k)]
                                            +residual[cell_getid(cdim,2*i,2*j+1,2*k)] + residual[cell_getid(cdim,2*i,2*j,2*k+1)]
                                            +residual[cell_getid(cdim,2*i+1,2*j+1,2*k)] + residual[cell_getid(cdim,2*i+1,2*j,2*k+1)]
                                            +residual[cell_getid(cdim,2*i,2*j+1,2*k+1)] + residual[cell_getid(cdim,2*i+1,2*j+1,2*k+1)])/8.0;
      }
    }
  }
}

/* Niet heel netjes om hier twee functies voor te hebben, laten we dit later aanpassen */
void get_residual_array(double *pot, double *dens, int cdim[3], double multiplier, double *residual, double delta) {
  /* Make sure array gets set to zero */
  const int N = cdim[0];
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  for (int i = 0; i<N*N*N; i++)
      residual[i] = 0.0;
  //double m2 = (multiplier==1.0) ? 1.0 : (double) delta*cdim[0]/3.0; //Only extra multiplier if we are on the finest grid
  //double m2 = (double) (delta*cdim[0]/3.0); 
  double m2 = 1.0;
  for (int k=0; k<cdim[2]; k++){
    int k_plus = (k+1) % cdim[2];
    int k_min = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
    for (int j=0; j<cdim[1]; j++){
      int j_plus = (j+1) % cdim[1];
      int j_min = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
      for (int i=0; i<cdim[0]; i++) {
        int i_plus = (i+1) % cdim[0];
        int i_min = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
        residual[cell_getid(cdim,i,j,k)] = ((pot[cell_getid(cdim,i_min,j,k)]+ pot[cell_getid(cdim,i_plus,j,k)]+
                                                pot[cell_getid(cdim,i,j_min,k)]+pot[cell_getid(cdim,i,j_plus,k)]+
                                                pot[cell_getid(cdim,i,j,k_min)]+pot[cell_getid(cdim,i,j,k_plus)]
                                                - 6.0*pot[cell_getid(cdim,i,j,k)])/(delta*delta) - m2*multiplier*(dens[cell_getid(cdim,i,j,k)]-rhs_correction));
        }
    }
  }
}

void get_residual(double *pot, double *dens, int cdim[3], double multiplier, double *residual, double delta){
  *residual = 0.0;
  //double m2 = (multiplier==1.0) ? 1.0 : (double) (delta*cdim[0]/3.0); //Only extra multiplier if we are on the finest grid
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  //double m2 = (double) (delta*cdim[0]/3.0);
  double m2 = 1.0;
  for (int k=0; k<cdim[2]; k++){
    int k_plus = (k+1) % cdim[2];
    int k_min = (k-1>=0) ? (k-1) % cdim[2] : (k-1) % cdim[2] + cdim[2];
    for (int j=0; j<cdim[1]; j++){
      int j_plus = (j+1) % cdim[1];
      int j_min = (j-1>=0) ? (j-1) % cdim[1] : (j-1) % cdim[1] + cdim[1];
      for (int i=0; i<cdim[0]; i++) {
        int i_plus = (i+1) % cdim[0];
        int i_min = (i-1>=0) ? (i-1) % cdim[0] : (i-1) % cdim[0] + cdim[0];
        double res = ((pot[cell_getid(cdim,i_min,j,k)]+ pot[cell_getid(cdim,i_plus,j,k)]+
                                                pot[cell_getid(cdim,i,j_min,k)]+pot[cell_getid(cdim,i,j_plus,k)]+
                                                pot[cell_getid(cdim,i,j,k_min)]+pot[cell_getid(cdim,i,j,k_plus)]
                                                - 6.0*pot[cell_getid(cdim,i,j,k)])/(delta*delta) - m2*multiplier*(dens[cell_getid(cdim,i,j,k)]-rhs_correction));
      *residual += res*res;
      }
    }
  }
  int N = cdim[0];
  *residual = sqrt(*residual/(N*N*N));
}

void perform_red_black_sweep(double *pot, double *dens, int cdim[3], double multiplier, double delta) {
  //double m2 = (multiplier==1.0) ? 1.0 : (double) delta*cdim[0]/3.0; //Only extra multiplier if we are on the finest grid
  float rhs_correction = (multiplier==1.0) ? 0.0 : 1.0;
  //double m2 = (double) (delta*cdim[0]/3.0);
  double m2 = 1;
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
                                              delta*delta*m2*multiplier*(dens[cell_getid(cdim,i,j,k)]-rhs_correction))/6.0);
          
        }
      }
    }
  }
}

//First semi-uniform, then adds a random value between 0 and 200 to every cell
void set_initial_guess(double *potential_array, int cdim[3]) {
  const int N = cdim[0];
  for (int i = 0; i<N*N*N; i++) {
    potential_array[i] = -1000;
  }
  for (int k=0; k<cdim[2]; k++){
    for (int j=0; j<cdim[1]; j++){
      for (int i=0; i<cdim[0]; i++) {
        if ((i+j+k)%2 != 0)
          continue;
      const size_t cid = cell_getid(cdim, i, j, k);
      potential_array[cid] += 2000;
      }  
    }
  } 
  srand(time(NULL));
  for (int i = 0; i<N*N*N; i++){
    potential_array[i] += rand() % (200);
  }
}

void density_to_cells(struct cic_mapper_data* data, struct cell* top_cells, int nr_cells, int cdim[3]) {
  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/test_pm_solution/density_array16.txt", "w");
  for (int i=0;i<cdim[0];i++) {
    for (int j=0;j<cdim[1];j++) {
      for (int k=0;k<cdim[2];k++) {
        size_t cid = cell_getid(cdim,i,j,k);
        top_cells[cid].CIC_density = data->rho[cid];
      }
    }
  }
  
  //fclose(fptr); 
}

void get_pm_potential(struct cic_mapper_data* data, const int N, const double box_size, struct threadpool* tp, int cdim[3]) {
  const int N_half = N / 2;
  double* rho_copy = malloc(N*N*N *sizeof(double));
  memcpy(rho_copy, data->rho, N*N*N* sizeof(double));
  if (rho_copy == NULL) error("Error allocating memory for density mesh");

  //double* density = malloc(N*N*N *sizeof(double));
  //memcpy(density, data->rho, N*N*N* sizeof(double));

  //double mean_density = get_mean_density(density, N);

  //double delta = box_size/N; 

  //double sum = 0.0;  
  //for (int i = 0; i < N*N*N; i++) {
      //sum += rho_copy[i];
  //}
  //message("The mean density is %lf", sum/(N*N*N));

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
  message("The msq of the potential is %lf", sum/(N*N*N));

  //mesh_to_gpart_CIC_mapper()

  //double residual = 0.;
  //mean_density *= 4.*M_PI*N*N*N/(box_size*box_size*box_size);
  //message("The mean density passed is %lf", mean_density);

  //get_residual(rho_copy, density, cdim, mean_density, &residual, delta);
  //message("The residual for the pm solution is %lf", residual);

  //Export the potential array to a .txt file to be read in Python
  //int cdim[3] = {N,N,N};
  //double fac = box_size/N;
  
  //message("Exporting PM FFT data");
  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/new_AMR_tests/missing_cell_check/FFT_16_noncropped.txt", "w");
  //for (int k=0; k<N; k++) {
    //for (int j=0; j<N; j++) {
      //for (int i=0; i<N; i++) {
        //if (i>=2 || j>=2 || k>=2) fprintf(fptr, "%.15g %.15g %.15g %.15g \n", rho_copy[cell_getid(cdim, i,j,k)], (double) i*fac, (double) j*fac, (double) k*fac);
        //fprintf(fptr, "%.15g %.15g %.15g %.15g \n", rho_copy[cell_getid(cdim, i,j,k)], (double) i*fac, (double) j*fac, (double) k*fac);
      //}
    //}
  //}
  //fclose(fptr); 
  //sleep(10);

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);
  free(rho_copy);
  //free(density);
  sleep(10);
}

double get_mean_density(double *dens, const int N) {
  double sum = 0.0;
  double mean_density = 0.0;

  //if (N == 64) {
  //FILE *fptr;
  //fptr = fopen("/data1/vandervlugt/PythonFiles/FMG_test/density_64.txt", "w");
  //for (int i=0; i<N*N*N; i++) {
    //fprintf(fptr, "%lf\n", dens[i]);
  //}
  //fclose(fptr);
  //}

  // Find the sum of all elements
  for (int i = 0; i < N*N*N; i++) {
      sum += dens[i];
  }
  mean_density = sum/(N*N*N);
  for (int i = 0; i < N*N*N; i++) {
      dens[i] = dens[i] / mean_density;
  }

  message("The mean density is %lf", mean_density);
  return mean_density;
}