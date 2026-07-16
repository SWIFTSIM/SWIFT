/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk)
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

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/**
 * @brief A frontier of cells waiting to be split.
 *
 * A frontier is a flat array of cells, at a single depth, which need to be
 * considered for splitting.
 */
struct space_split_frontier {

  /*! Cells to process at this BFS level. */
  struct cell **cells;

  /*! Number of valid entries in cells. */
  int count;

  /*! Number of slots allocated in cells. */
  int capacity;
};

/**
 * @brief Grow a frontier so it can hold at least "needed" pointers.
 *
 * If the frontier is already large enough, this does nothing. Otherwise, it
 * will allocate a new array of pointers, copy the old pointers into it, free
 * the old array and attach the new array to the frontier.
 *
 * @param frontier The frontier to resize.
 * @param needed Required capacity (in pointers).
 */
static void space_split_frontier_ensure_capacity(
    struct space_split_frontier *frontier, const int needed) {

  /* Nothing to do if we already have enough space. */
  if (needed <= frontier->capacity) return;

  /* Get the array for the new cell pointers. */
  struct cell **new_cells = NULL;
  if (swift_memalign("split_frontier", (void **)&new_cells,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct cell *) * needed) != 0) {
    error("Failed to (re)allocate BFS frontier of size %d", needed);
  }

  /* Copy the old cell pointers into the new array. */
  for (int i = 0; i < frontier->count; i++) new_cells[i] = frontier->cells[i];

  /* Free the old array. */
  if (frontier->cells != NULL) swift_free("split_frontier", frontier->cells);

  /* Attach the new array to the frontier. */
  frontier->cells = new_cells;
  frontier->capacity = needed;
}

/**
 * @brief Append a cell to a next frontier.
 *
 * If the frontier is full, it will be expanded to hold 2x the number of cells.
 *
 * @param frontier The frontier being produced.
 * @param c The cell to append.
 */
__attribute__((always_inline)) static INLINE void space_split_frontier_append(
    struct space_split_frontier *frontier, struct cell *c) {

  /* Do we need to expand the frontier's array of cell pointers? */
  if (frontier->count == frontier->capacity) {

    /* Frontier is full, so we need to grow it. If the capacity is non-zero, we
     * double it. Otherwise, we start with a capacity of 256 (though the
     * capacity by the time this is called is almost certainly non-zero, better
     * to be safe than sorry). */
    const int new_capacity =
        frontier->capacity == 0 ? 256 : 2 * frontier->capacity;
    space_split_frontier_ensure_capacity(frontier, new_capacity);
  }

  /* Append the cell to the frontier. */
  frontier->cells[frontier->count++] = c;
}

/**
 * @brief Merge thread-local next-frontier buffers into a flat frontier.
 *
 * Each thread accumulates the next frontier in its own local buffer to avoid
 * contention. This function merges those local buffers into a single flat
 * frontier for the next round of processing.
 *
 * @param next_frontier The flat frontier produced for the next round.
 * @param local_frontiers The per-thread frontier buffers to concatenate.
 * @param nr_threads Number of thread-local frontier buffers.
 */
static void space_split_merge_local_frontiers(
    struct space_split_frontier *next_frontier,
    struct space_split_frontier *local_frontiers, const int nr_threads) {

  /* Count the total number of cells in all the local frontiers. */
  int needed = 0;
  for (int tid = 0; tid < nr_threads; tid++) {
    needed += local_frontiers[tid].count;
  }

  /* Ensure the next frontier is big enough to hold "needed" cells. */
  space_split_frontier_ensure_capacity(next_frontier, needed);
  next_frontier->count = needed;

  /* Copy the local frontiers into the next frontier. */
  int offset = 0;
  for (int tid = 0; tid < nr_threads; tid++) {

    /* Pointer to the local frontier for this thread. */
    struct space_split_frontier *local = &local_frontiers[tid];

    /* Copy the local frontier into the next frontier. */
    for (int i = 0; i < local->count; i++) {
      next_frontier->cells[offset + i] = local->cells[i];
    }

    /* Update the offset and reset the local frontier count. */
    offset += local->count;
    local->count = 0;
  }
}

/**
 * @brief Recursively split a cell, appending subtrees to a frontier once the
 *        DFS depth budget is exhausted.
 *
 * If buff, sbuff, bbuff, gbuff and sink_buff are all NULL, fresh local
 * buffers are allocated and populated for this DFS chunk exactly as would be
 * done at the top level for the DFS approach.
 *
 * @param s The #space.
 * @param next The frontier being produced for the next BFS level.
 * @param c The cell to process.
 * @param buff A buffer for particle sorting, should be of size at least
 *        c->hydro.count or NULL.
 * @param sbuff A buffer for particle sorting, should be of size at least
 *        c->stars.count or NULL.
 * @param bbuff A buffer for particle sorting, should be of size at least
 *        c->black_holes.count or NULL.
 * @param gbuff A buffer for particle sorting, should be of size at least
 *        c->grav.count or NULL.
 * @param sink_buff A buffer for particle sorting, should be of size at least
 *        c->sinks.count or NULL.
 * @param tpid The threadpool tid of the calling worker.
 * @param depth_remaining Number of additional DFS levels to descend before
 *        enqueuing survivors (space_dfs_levels_per_frontier - 1 on entry
 *        from the mapper).
 */
static void space_split_recursive(
    struct space *s, struct space_split_frontier *next, struct cell *c,
    struct cell_buff *restrict buff, struct cell_buff *restrict sbuff,
    struct cell_buff *restrict bbuff, struct cell_buff *restrict gbuff,
    struct cell_buff *restrict sink_buff, const short int tpid,
    const int depth_remaining) {

  /* Unpack cell information. */
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int with_self_gravity = s->with_self_gravity;
  const int depth = c->depth;

  /* Top-level cells inherit the worker's tpid (preserves the legacy
   * behaviour where a top cell's progeny share its tpid). */
  if (c->depth == 0) c->tpid = tpid;

  /* Do we need to allocate the sorting buffers? If so, do it here. */
  const int allocate_buffer = (buff == NULL && gbuff == NULL && sbuff == NULL &&
                               bbuff == NULL && sink_buff == NULL);
  if (allocate_buffer) {
    space_split_allocate_and_fill_buffers(c, &buff, &gbuff, &sbuff, &bbuff,
                                          &sink_buff);
  }

  /* If the depth is too large, we have a problem and should stop. */
  if (depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting cells, aborting. This is "
        "most likely due to having too many particles at the exact same "
        "position, making the construction of a tree impossible.",
        space_cell_maxdepth);
  }

  /* Split or let it be? */
  if ((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize))) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny, tpid);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->hydro.count = 0;
      cp->grav.count = 0;
      cp->stars.count = 0;
      cp->sinks.count = 0;
      cp->black_holes.count = 0;
      cp->hydro.count_total = 0;
      cp->grav.count_total = 0;
      cp->sinks.count_total = 0;
      cp->stars.count_total = 0;
      cp->black_holes.count_total = 0;
      cp->hydro.ti_old_part = c->hydro.ti_old_part;
      cp->grav.ti_old_part = c->grav.ti_old_part;
      cp->grav.ti_old_multipole = c->grav.ti_old_multipole;
      cp->stars.ti_old_part = c->stars.ti_old_part;
      cp->sinks.ti_old_part = c->sinks.ti_old_part;
      cp->black_holes.ti_old_part = c->black_holes.ti_old_part;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->width[0] = c->width[0] / 2;
      cp->width[1] = c->width[1] / 2;
      cp->width[2] = c->width[2] / 2;
      cp->dmin = c->dmin / 2;
      cp->h_min_allowed = cp->dmin * 0.5 * (1. / kernel_gamma);
      cp->h_max_allowed = cp->dmin * (1. / kernel_gamma);
      if (k & 4) cp->loc[0] += cp->width[0];
      if (k & 2) cp->loc[1] += cp->width[1];
      if (k & 1) cp->loc[2] += cp->width[2];
      cp->depth = c->depth + 1;
      cp->split = 0;
      cp->hydro.h_max = 0.f;
      cp->hydro.h_max_active = 0.f;
      cp->hydro.dx_max_part = 0.f;
      cp->hydro.dx_max_sort = 0.f;
      cp->stars.h_max = 0.f;
      cp->stars.h_max_active = 0.f;
      cp->stars.dx_max_part = 0.f;
      cp->stars.dx_max_sort = 0.f;
      cp->sinks.h_max = 0.f;
      cp->sinks.h_max_active = 0.f;
      cp->sinks.dx_max_part = 0.f;
      cp->black_holes.h_max = 0.f;
      cp->black_holes.h_max_active = 0.f;
      cp->black_holes.dx_max_part = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->top = c->top;
      cp->super = NULL;
      cp->hydro.super = NULL;
      cp->grav.super = NULL;
      cp->flags = 0;
      star_formation_logger_init(&cp->stars.sfh);
#ifdef WITH_MPI
      cp->mpi.tag = -1;
#endif  // WITH_MPI
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
      cell_assign_cell_index(cp, c);
#endif
    }

    /* Split the cell's particle data. */
    cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts,
               c->black_holes.parts - s->bparts, c->sinks.parts - s->sinks,
               buff, sbuff, bbuff, gbuff, sink_buff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff, *progeny_bbuff = bbuff,
                     *progeny_sink_buff = sink_buff;

    for (int k = 0; k < 8; k++) {
      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->hydro.count == 0 && cp->grav.count == 0 && cp->stars.count == 0 &&
          cp->black_holes.count == 0 && cp->sinks.count == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {

        /* We need to keep splitting, if we are are not at a frontier limit
         * (i.e. if we have not finished a DFS portion of the BFS splitting)
         * then we just recurse. Otherwise, append to the next frontier. */
        if (depth_remaining > 0) {
          space_split_recursive(s, next, cp, progeny_buff, progeny_sbuff,
                                progeny_bbuff, progeny_gbuff, progeny_sink_buff,
                                tpid, depth_remaining - 1);
        } else {
          space_split_frontier_append(next, cp);
        }

        /* Update the pointers in the buffers */
        progeny_buff += cp->hydro.count;
        progeny_gbuff += cp->grav.count;
        progeny_sbuff += cp->stars.count;
        progeny_bbuff += cp->black_holes.count;
        progeny_sink_buff += cp->sinks.count;
      }
    }

  } /* Split or let it be? */

  /* Otherwise, this cell remains a leaf: finalise it directly from its own
   * particles instead of aggregating from progeny. */
  else {

    /* Finalise the leaf cell. */
    space_split_finalise_leaf(s, c);
  }

  /* Clean up. */
  if (allocate_buffer) {
    if (buff != NULL) swift_free("tempbuff", buff);
    if (gbuff != NULL) swift_free("tempgbuff", gbuff);
    if (sbuff != NULL) swift_free("tempsbuff", sbuff);
    if (bbuff != NULL) swift_free("tempbbuff", bbuff);
    if (sink_buff != NULL) swift_free("temp_sink_buff", sink_buff);
  }
}

/**
 * @brief Recursively aggregate progeny properties onto parent cells.
 *
 * Once we have completed the BFS splitting phase we can propagate cell level
 * properties from the leaves up to the root. This is a depth-first traversal of
 * the tree. If we are running with self-gravity, we also build the multipoles
 * of each cell from the bottom up (M2M).
 *
 * @param s The #space the cell lives in.
 * @param c The cell to aggregate.
 */
static void space_split_aggregate_recursive(struct space *s, struct cell *c) {

  /* Leaves were already finalised during the BFS split phase. */
  if (!c->split) return;

  /* Aggregate progeny first. */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL)
      space_split_aggregate_recursive(s, c->progeny[k]);
  }

  /* Reset to the identity of the max/min reduction that
   * space_split_accumulate_props() performs below, one surviving child at a
   * time. */
  c->hydro.h_max = 0.f;
  c->hydro.h_max_active = 0.f;
  c->stars.h_max = 0.f;
  c->stars.h_max_active = 0.f;
  c->black_holes.h_max = 0.f;
  c->black_holes.h_max_active = 0.f;
  c->sinks.h_max = 0.f;
  c->sinks.h_max_active = 0.f;
  c->hydro.ti_end_min = max_nr_timesteps;
  c->hydro.ti_beg_max = 0;
  c->rt.ti_rt_end_min = max_nr_timesteps;
  c->rt.ti_rt_beg_max = 0;
  c->rt.ti_rt_min_step_size = max_nr_timesteps;
  c->grav.ti_end_min = max_nr_timesteps;
  c->grav.ti_beg_max = 0;
  c->stars.ti_end_min = max_nr_timesteps;
  c->stars.ti_beg_max = 0;
  c->sinks.ti_end_min = max_nr_timesteps;
  c->sinks.ti_beg_max = 0;
  c->black_holes.ti_end_min = max_nr_timesteps;
  c->black_holes.ti_beg_max = 0;
  int maxdepth = 0;

  /* Loop over progeny and accumulate their properties into this cell. */
  for (int k = 0; k < 8; k++) {
    const struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;

    /* Update the cell-wide properties. */
    space_split_accumulate_props(c, cp);

    /* Update the maximum depth. */
    maxdepth = max(maxdepth, cp->maxdepth);
  }

  /* Deal with the multipole */
  if (s->with_self_gravity) space_split_populate_multipole(c);

  /* The per-family summary fields were already written directly into c by
   * space_split_accumulate_props() as each child was folded in. */
  c->maxdepth = maxdepth;

  /* No runner owns this cell yet. We assign those during scheduling. */
  c->owner = -1;
}

/**
 * @brief #threadpool mapper function for one BFS level.
 *
 * This mapper will extract a cell from the current frontier and then call the
 * recursive splitting function on it. If the DFS depth limit is reached, the
 * surviving children will be appended to the thread-local next frontier buffer
 * for splitting in a later round. This gives a hybrid BFS/DFS splitting
 * strategy that makes the most of thread local particle sorting buffers while
 * still allowing for parallelism across expensive subtrees.
 *
 * @param map_data Pointer to the start of a chunk of cell pointers.
 * @param num_cells Number of cells in this chunk.
 * @param extra_data Pointer to an array containing the #space and the
 *        thread-local next frontier buffers.
 */
static void space_split_frontier_mapper(void *map_data, int num_cells,
                                        void *extra_data) {

  /* Unpack the inputs. */
  void **mapper_data = (void **)extra_data;
  struct space *s = (struct space *)mapper_data[0];
  struct cell **cells = (struct cell **)map_data;

  /* Threadpool id of current thread. */
  const short int tpid = threadpool_gettid();

  /* Get the local frontier we will append to once we reach the DFS depth limit
   * for this frontier. */
  struct space_split_frontier *local_frontiers =
      (struct space_split_frontier *)mapper_data[1];
  struct space_split_frontier *next = &local_frontiers[tpid];

  /* Loop over the cells in this chunk of the frontier. */
  for (int i = 0; i < num_cells; i++) {

    /* Split this cell recursively. */
    space_split_recursive(s, next, cells[i], NULL, NULL, NULL, NULL, NULL, tpid,
                          space_dfs_levels_per_frontier - 1);
  }
}

/**
 * @brief #threadpool mapper function for accumulating cell properties after
 * splitting.
 *
 * This function will loop over all top level cells and call
 * #space_split_aggregate_recursive on each one, which will recurse all the way
 * to the leaves and accumulate cell properties and popualte multipoles from the
 * bottom up.
 *
 * @param map_data Pointer to the start of a chunk of cell indices.
 * @param num_cells Number of indices in this chunk.
 * @param extra_data Pointer to the parent #space.
 */
static void space_split_aggregate_mapper(void *map_data, int num_cells,
                                         void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Initialise some global information about the top-level m-poles */
  float min_a_grav = FLT_MAX;
  float max_softening = 0.f;
  float max_mpole_power[SELF_GRAVITY_MULTIPOLE_ORDER + 1] = {0.f};

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {

    /* Get this cell and aggregate its progeny's properties into it. */
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    space_split_aggregate_recursive(s, c);

    /* If we are running with self-gravity, collect the global min/max values
     * of the multipole properties. */
    if (s->with_self_gravity) {
      min_a_grav =
          min(min_a_grav, c->grav.multipole->m_pole.min_old_a_grav_norm);
      max_softening =
          max(max_softening, c->grav.multipole->m_pole.max_softening);
      for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
        max_mpole_power[n] =
            max(max_mpole_power[n], c->grav.multipole->m_pole.power[n]);
    }

    /* The deepest cell in this top tree determines the contribution to
     * the global max depth. */
    atomic_max(&s->maxdepth, c->maxdepth);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
  }
#endif

  /* Update the global min/max values of the multipole properties. */
  atomic_min_f(&s->min_a_grav, min_a_grav);
  atomic_max_f(&s->max_softening, max_softening);
  for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
    atomic_max_f(&s->max_mpole_power[n], max_mpole_power[n]);
}

/**
 * @brief Set up the frontier buffers for BFS splitting.
 *
 * Allocates the frontiers (this_level and next_level), one for the current BFS
 * level and one to populate for the next pass, plus the local_frontiers array
 * of per-thread buffers that each worker appends surviving children to
 * (avoiding contention) before they are merged into a flat frontier for the
 * next round.
 *
 * @param this_level Pointer to the frontier pointer for the current BFS
 *        level, allocated in place.
 * @param next_level Pointer to the frontier pointer for the next BFS level,
 *        allocated in place.
 * @param local_frontiers Pointer to the thread-local frontier buffers,
 *        allocated in place.
 * @param nr_threads Number of threads in the #threadpool.
 */
static void space_split_init_frontiers(
    struct space_split_frontier **this_level,
    struct space_split_frontier **next_level,
    struct space_split_frontier **local_frontiers, const int nr_threads) {

  /* Allocate the 2 frontiers, one for the current BFS level and one to
   * populate for the next pass. */
  *this_level = (struct space_split_frontier *)swift_malloc(
      "split_frontier", sizeof(struct space_split_frontier));
  *next_level = (struct space_split_frontier *)swift_malloc(
      "split_frontier", sizeof(struct space_split_frontier));
  if (*this_level == NULL || *next_level == NULL)
    error("Failed to allocate BFS frontier level buffers.");
  **this_level = (struct space_split_frontier){NULL, 0, 0};
  **next_level = (struct space_split_frontier){NULL, 0, 0};

  /* Allocate the thread-local frontier buffers. It's into these that we
   * will append the surviving children on each thread to avoid contention. Once
   * all threads have finished, we will merge these into a single flat frontier
   * for the next round. */
  *local_frontiers = (struct space_split_frontier *)swift_malloc(
      "split_frontier", sizeof(struct space_split_frontier) * nr_threads);

  /* Make sure the thread-local frontier buffers are safe and zeroed. */
  if (*local_frontiers == NULL)
    error("Failed to allocate thread-local BFS frontiers.");
  bzero(*local_frontiers, sizeof(struct space_split_frontier) * nr_threads);
}

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * This function employs a hybrid BFS/DFS splitting strategy to allow for
 * parallelism across expensive subtrees while still making the most of
 * thread-local particle sorting buffers.
 *
 * Unlike the DFS splitting strategy this process is done in 2 phases:
 *
 * 1. Hybrid split phase: The initial frontier is just the local top-level
 * cells. Each round, the current frontier is distributed over threads and
 * each of its cells is recursed into with a depth-first search (DFS) of
 * space_dfs_levels_per_frontier levels. Once that depth is reached, any
 * surviving children are appended to a thread-local next frontier buffer. After
 * all threads have finished, the local frontiers are merged into a single flat
 * frontier for the next round. This process then repeats until the frontier is
 * empty, at which point all cells have been split.
 *
 * 2. Aggregation phase: Once all cells have been split, we need to derive
 * cell properties from the particles at the leaves and then propagate these
 * up the tree to the top. This includes smoothing length and timestep
 * information. Additionally, if we are running with self-gravity, we also need
 * to build the multipoles of each cell from the bottom up (M2M).
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split_bfs_frontiers(struct space *s, int verbose) {

  /* Nothing to do if we have no local cells with particles. */
  const int nr_local_cells = s->nr_local_cells_with_particles;
  if (nr_local_cells == 0) return;

  /* Get the threadpool and number of threads. */
  struct threadpool *tp = &s->e->threadpool;
  const int nr_threads = tp->num_threads;

  /* Prepare the frontier buffers. These will hold the current, next and
   * thread-local collections of cells for splitting. */
  ticks tic = getticks();
  struct space_split_frontier *this_level;
  struct space_split_frontier *next_level;
  struct space_split_frontier *local_frontiers;
  space_split_init_frontiers(&this_level, &next_level, &local_frontiers,
                             nr_threads);

  /* Both mapper passes share the same #space and thread-local frontier
   * buffers. */
  void *mapper_data[2] = {s, local_frontiers};

  if (verbose)
    message("Frontier setup: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Hybrid BFS/DFS splitting phase */

  /* Seed the first frontier with the local top-level cells themselves. */
  tic = getticks();
  space_split_frontier_ensure_capacity(this_level, nr_local_cells);
  for (int i = 0; i < nr_local_cells; i++) {
    this_level->cells[i] = &s->cells_top[s->local_cells_with_particles_top[i]];
  }
  this_level->count = nr_local_cells;

  /* Begin at frontier level 0. */
  int frontier_level = 0;

  /* Keep splitting successive frontiers of cells in parallel (each cell
   * recursing DFS-style for space_dfs_levels_per_frontier levels before any
   * surviving children are deferred to the next frontier) until nothing is left
   * to split. */
  while (this_level->count > 0) {

    /* Split this frontier's cells in parallel. */
    threadpool_map(tp, space_split_frontier_mapper, this_level->cells,
                   this_level->count, sizeof(struct cell *),
                   threadpool_auto_chunk_size, mapper_data);

    /* Merge the thread-local next frontier buffers into a single flat
     * frontier for the next round. */
    space_split_merge_local_frontiers(next_level, local_frontiers, nr_threads);

    /* Swap this and next frontier for the next round. */
    struct space_split_frontier *tmp = this_level;
    this_level = next_level;
    next_level = tmp;
    frontier_level++;

    /* Report the number of cells in this frontier. */
    if (verbose && this_level->count > 0) {
      message("Frontier level %d has %d cells.", frontier_level,
              this_level->count);
    }
  }

  if (verbose)
    message("BFS split loop: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Clean up the frontier buffers. */
  if (this_level->cells != NULL)
    swift_free("split_frontier", this_level->cells);
  if (next_level->cells != NULL)
    swift_free("split_frontier", next_level->cells);
  swift_free("split_frontier", this_level);
  swift_free("split_frontier", next_level);
  for (int tid = 0; tid < nr_threads; tid++) {
    if (local_frontiers[tid].cells != NULL)
      swift_free("split_frontier", local_frontiers[tid].cells);
  }
  swift_free("split_frontier", local_frontiers);

  /* Aggregation phase */

  /* Map over all local top-level cells with particles and aggregate their
   * progeny's properties into them. This is a depth-first traversal of the
   * tree, so we will recurse all the way to the leaves and accumulate cell
   * properties and populate multipoles from the bottom up (if running with
   * self-gravity). */
  tic = getticks();
  threadpool_map(tp, space_split_aggregate_mapper,
                 s->local_cells_with_particles_top, nr_local_cells, sizeof(int),
                 threadpool_auto_chunk_size, s);

  if (verbose)
    message("Aggregate pass: %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
