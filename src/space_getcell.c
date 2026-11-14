/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Some standard headers. */
#include <stdint.h>
#include <stdlib.h>

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "threadpool.h"


/**
 * @brief List of the chunks of memory allocated for one buffer of unused
 * sub-cells or multipoles.
 *
 * Note: The chunks are never freed.
 */
struct space_sub_chunks {

  /*! The start of each chunk. */
  void **chunks;

  /*! Number of chunks recorded. */
  int count;

  /*! Allocated size of the chunks array. */
  int size;
};

/**
 * @brief One chunk to re-link and the chunk that follows it in its buffer.
 */
struct space_relink_item {
  void *chunk;
  void *next;
};

/*! Chunks of the buffers of sub-cells, one list per buffer. */
static struct space_sub_chunks *space_cell_chunks = NULL;

/*! Chunks of the buffers of multipoles, one list per buffer. */
static struct space_sub_chunks *space_multipole_chunks = NULL;

/*! Number of buffers (i.e. of entries in the two lists). */
static int space_sub_chunks_nr_buffers = 0;

/*! Initial size of the lists */
static const int space_sub_chunks_initial_size = 64;

/**
 * @brief Release the lists of chunks of the buffers of sub-cells and
 * multipoles.
 *
 * The chunks themselves are not freed.
 */
void space_cell_sub_chunks_clean(void) {

  for (int k = 0; k < space_sub_chunks_nr_buffers; k++) {
    free(space_cell_chunks[k].chunks);
    free(space_multipole_chunks[k].chunks);
  }
  free(space_cell_chunks);
  free(space_multipole_chunks);
  space_cell_chunks = NULL;
  space_multipole_chunks = NULL;
  space_sub_chunks_nr_buffers = 0;
}

/**
 * @brief Prepare the empty lists of chunks of the buffers of sub-cells and
 * multipoles.
 *
 * Note: Must be called whenever the buffers (#space cells_sub and
 * multipoles_sub) are allocated empty, i.e. in engine_config().
 *
 * @param nr_buffers The number of buffers.
 */
void space_cell_sub_chunks_init(const int nr_buffers) {

  space_cell_sub_chunks_clean();
  space_cell_chunks = (struct space_sub_chunks *)calloc(
      nr_buffers, sizeof(struct space_sub_chunks));
  space_multipole_chunks = (struct space_sub_chunks *)calloc(
      nr_buffers, sizeof(struct space_sub_chunks));
  if (space_cell_chunks == NULL || space_multipole_chunks == NULL)
    error("Failed to allocate the list of sub-cell chunks.");
  space_sub_chunks_nr_buffers = nr_buffers;
}

/**
 * @brief Record a newly allocated chunk of one buffer.
 *
 * @param lists The lists (cells or multipoles).
 * @param tpid The buffer the chunk belongs to.
 * @param chunk The start of the chunk.
 */
static void space_sub_chunks_record(struct space_sub_chunks *lists,
                                    const short int tpid, void *chunk) {

  /* No list: return */
  if (lists == NULL) return;

  /* Get the list of our thread */
  if (tpid < 0 || tpid >= space_sub_chunks_nr_buffers)
    error("Invalid buffer index %d (%d buffers).", tpid,
          space_sub_chunks_nr_buffers);

  struct space_sub_chunks *l = &lists[tpid];

  /* Does the list need to grow? */
  if (l->count == l->size) {
    l->size = l->size > 0 ? 2 * l->size : space_sub_chunks_initial_size;
    l->chunks = (void **)realloc(l->chunks, l->size * sizeof(void *));
    if (l->chunks == NULL) error("Failed to grow the list of chunks.");
  }

  /* Add the chunk to the list */
  l->chunks[l->count++] = chunk;
}

/**
 * @brief Comparison function sorting pointers by increasing address.
 */
static int space_sub_chunks_cmp(const void *a, const void *b) {
  const uintptr_t pa = (uintptr_t)(*(void *const *)a);
  const uintptr_t pb = (uintptr_t)(*(void *const *)b);
  return (pa > pb) - (pa < pb);
}

/**
 * @brief #threadpool mapper function re-linking the cells of chunks of
 * sub-cells: the cells of a chunk are linked in address order and the last
 * one is linked to the first cell of the next chunk of its buffer.
 *
 * @param map_data Pointer towards the #space_relink_item of the chunks.
 * @param num_elements The number of chunks to treat.
 * @param extra_data Unused.
 */
static void space_relink_cells_mapper(void *map_data, int num_elements,
                                      void *extra_data) {

  struct space_relink_item *items = (struct space_relink_item *)map_data;
  for (int i = 0; i < num_elements; i++) {
    struct cell *c = (struct cell *)items[i].chunk;
    for (int k = 0; k < space_cellallocchunk - 1; k++) c[k].next = &c[k + 1];
    c[space_cellallocchunk - 1].next = (struct cell *)items[i].next;
  }
}

/**
 * @brief #threadpool mapper function re-linking the multipoles of chunks of
 * multipoles: the multipoles of a chunk are linked in address order and the
 * last one is linked to the first multipole of the next chunk of its buffer.
 *
 * @param map_data Pointer towards the #space_relink_item of the chunks.
 * @param num_elements The number of chunks to treat.
 * @param extra_data Unused.
 */
static void space_relink_multipoles_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct space_relink_item *items = (struct space_relink_item *)map_data;
  for (int i = 0; i < num_elements; i++) {
    struct gravity_tensors *m = (struct gravity_tensors *)items[i].chunk;
    for (int k = 0; k < space_cellallocchunk - 1; k++) m[k].next = &m[k + 1];
    m[space_cellallocchunk - 1].next = (struct gravity_tensors *)items[i].next;
  }
}

/**
 * @brief Re-link one kind of buffers (cells or multipoles) in address order.
 *
 * Each buffer is rebuilt from the chunks it allocated, sorted by address and
 * chained one after the other. The previous content of the buffers is
 * discarded, so all their elements must be free.
 *
 * @param s The #space.
 * @param lists The lists of chunks of the buffers (cells or multipoles).
 * @param buffers The heads of the buffers (#space cells_sub or
 * multipoles_sub), overwritten with the first chunk of each buffer (or
 * @c NULL if the buffer has no chunk).
 * @param fun The mapper re-linking the elements of the chunks.
 *
 * @return The number of chunks re-linked.
 */
static int space_relink_buffers(struct space *s,
				struct space_sub_chunks *lists,
				void **buffers, threadpool_map_function fun) {

  int nr_chunks = 0;
  for (int k = 0; k < space_sub_chunks_nr_buffers; k++)
    nr_chunks += lists[k].count;

  /* Nothing to do? */
  if (nr_chunks == 0) return 0;

  struct space_relink_item *items = (struct space_relink_item *)malloc(
      nr_chunks * sizeof(struct space_relink_item));
  if (items == NULL) error("Failed to allocate the list of chunks.");

  int n = 0;
  for (int k = 0; k < space_sub_chunks_nr_buffers; k++) {
    struct space_sub_chunks *l = &lists[k];
    qsort(l->chunks, l->count, sizeof(void *), space_sub_chunks_cmp);
    for (int i = 0; i < l->count; i++) {
      items[n].chunk = l->chunks[i];
      items[n].next = (i + 1 < l->count) ? l->chunks[i + 1] : NULL;
      n++;
    }
    buffers[k] = (l->count > 0) ? l->chunks[0] : NULL;
  }

  threadpool_map(&s->e->threadpool, fun, items, nr_chunks,
                 sizeof(struct space_relink_item), threadpool_auto_chunk_size,
                 NULL);
  free(items);
  return nr_chunks;
}

/**
 * @brief Re-order the buffers of unused sub-cells and multipoles by memory
 * address.
 *
 * Cells come back to the buffers in whatever order the trees are freed in,
 * and possibly to the buffer of another thread than the one that allocated
 * them, so after many rebuilds consecutive cells of a buffer, and hence the
 * progeny of a cell in a freshly split tree, are scattered over all the
 * memory of the pools. Once *all* the sub-cells are back in the buffers,
 * the buffers together hold exactly all the chunks ever allocated, so we can
 * rebuild every buffer from the chunks it allocated itself, chained chunk by
 * chunk in address order. The next trees are then built from contiguous
 * memory first touched by the thread using it.
 *
 * @param s The #space.
 */
void space_relink_cell_sub_buffers(struct space *s) {

  if (space_cell_chunks == NULL) return;

  /* Safety check: Are all the sub-cells back in the buffers? They should be,
   * as this is called once all the trees have been recycled. */
  if (s->tot_cells != s->nr_cells) {
#ifdef SWIFT_DEBUG_CHECKS
    error("%d sub-cells still in use after recycling all the trees.",
          s->tot_cells - s->nr_cells);
#endif
    if (s->e->verbose)
      message("%d sub-cells still in use, not re-ordering the buffers.",
              s->tot_cells - s->nr_cells);
    return;
  }

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  long long nr_free_cells = 0, nr_free_multipoles = 0;
  for (int k = 0; k < space_sub_chunks_nr_buffers; k++) {
    for (struct cell *c = s->cells_sub[k]; c != NULL; c = c->next)
      nr_free_cells++;
    if (s->with_self_gravity)
      for (struct gravity_tensors *m = s->multipoles_sub[k]; m != NULL;
           m = m->next)
        nr_free_multipoles++;
  }
#endif

  const int nr_cell_chunks =
      space_relink_buffers(s, space_cell_chunks, (void **)s->cells_sub,
                           space_relink_cells_mapper);
  const int nr_multipole_chunks =
      s->with_self_gravity
          ? space_relink_buffers(s, space_multipole_chunks,
                                 (void **)s->multipoles_sub,
                                 space_relink_multipoles_mapper)
          : 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (nr_free_cells != (long long)nr_cell_chunks * space_cellallocchunk)
    error("Found %lld free cells in the buffers but %d chunks of %d cells.",
          nr_free_cells, nr_cell_chunks, space_cellallocchunk);
  if (s->with_self_gravity &&
      nr_free_multipoles !=
          (long long)nr_multipole_chunks * space_cellallocchunk)
    error(
        "Found %lld free multipoles in the buffers but %d chunks of %d "
        "multipoles.",
        nr_free_multipoles, nr_multipole_chunks, space_cellallocchunk);
#endif

  if (s->e->verbose)
    message("Re-linked %d chunks of cells and %d of multipoles, took %.3f %s.",
            nr_cell_chunks, nr_multipole_chunks,
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

/**
 * @brief Get a new empty (sub-)#cell.
 *
 * The cells (and their multipoles) are taken from the head of the buffer of
 * thread @c tpid. If that buffer is empty, a new chunk of memory is allocated
 * (and recorded for space_relink_cell_sub_buffers()) and the cells are picked
 * from there.
 *
 * The buffers are not locked: a thread may only use its own buffer, i.e.
 * pass its threadpool id as @c tpid from within a #threadpool mapper. Outside
 * of a mapper, the main thread may use any buffer.
 *
 * @param s The #space.
 * @param nr_cells Number of #cell to pick up.
 * @param cells Array of @c nr_cells #cell pointers in which to store the
 *        new cells.
 * @param tpid ID of the threadpool thread whose buffer (cells_sub and
 *        multipoles_sub) is used.
 */
void space_getcells(struct space *s, int nr_cells, struct cell **cells,
                    const short int tpid) {

  /* For each requested cell... */
  for (int j = 0; j < nr_cells; j++) {

    /* Is the cell buffer empty? */
    if (s->cells_sub[tpid] == NULL) {
      if (swift_memalign("cells_sub", (void **)&s->cells_sub[tpid], cell_align,
                         space_cellallocchunk * sizeof(struct cell)) != 0)
        error("Failed to allocate more cells.");

      /* This allocation is never correctly freed, that is ok. */
      swift_ignore_leak(s->cells_sub[tpid]);

      /* Record the chunk such that the buffer can be re-ordered later. */
      space_sub_chunks_record(space_cell_chunks, tpid, s->cells_sub[tpid]);

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

      /* Record the chunk such that the buffer can be re-ordered later. */
      space_sub_chunks_record(space_multipole_chunks, tpid,
                              s->multipoles_sub[tpid]);

      /* Constructed a linked list */
      for (int k = 0; k < space_cellallocchunk - 1; k++)
        s->multipoles_sub[tpid][k].next = &s->multipoles_sub[tpid][k + 1];
      s->multipoles_sub[tpid][space_cellallocchunk - 1].next = NULL;
    }

    /* Pick off the next cell. */
    cells[j] = s->cells_sub[tpid];
    s->cells_sub[tpid] = cells[j]->next;

    /* Hook the multipole */
    if (s->with_self_gravity) {
      cells[j]->grav.multipole = s->multipoles_sub[tpid];
      s->multipoles_sub[tpid] = cells[j]->grav.multipole->next;
    }
  }

  /* Count the new cells. */
  atomic_add(&s->tot_cells, nr_cells);

  /* Init some things in the cell we just got. */
  for (int j = 0; j < nr_cells; j++) {
    cell_free_hydro_sorts(cells[j]);
    cell_free_stars_sorts(cells[j]);

    struct gravity_tensors *temp = cells[j]->grav.multipole;
    bzero(cells[j], sizeof(struct cell));
    cells[j]->grav.multipole = temp;
    cells[j]->nodeID = -1;
    cells[j]->tpid = tpid;
    if (lock_init(&cells[j]->hydro.lock) != 0 ||
        lock_init(&cells[j]->hydro.extra_sort_lock) != 0 ||
        lock_init(&cells[j]->grav.plock) != 0 ||
        lock_init(&cells[j]->grav.mlock) != 0 ||
        lock_init(&cells[j]->stars.lock) != 0 ||
        lock_init(&cells[j]->sinks.lock) != 0 ||
        lock_init(&cells[j]->sinks.sink_formation_lock) != 0 ||
        lock_init(&cells[j]->black_holes.lock) != 0 ||
        lock_init(&cells[j]->stars.star_formation_lock) != 0 ||
        lock_init(&cells[j]->grav.star_formation_lock) != 0)
      error("Failed to initialize cell spinlocks.");
  }
}
