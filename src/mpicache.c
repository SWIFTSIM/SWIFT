/* This file is part of SWIFT.
 * Copyright (c) 2020 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/**
 *  @file mpicache.c
 *  @brief file of routines to cache MPI task messages
 */
/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

/* Local defines. */
#include "mpicache.h"
#include "scheduler.h"

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"
#include "memuse_rnodes.h"
#include "task.h"

/* Size of the cache entry buffer increment. */
#define CACHE_INC_SIZE 5000

/* Bit shift to accomodate all the bits of the maximum subtype. */
int mpicache_subtype_shift = 0;

/* Comparator function to sort by node, subtype, tag order. */
#ifdef WITH_MPI
static int mpicache_entry_cmp(const void *a, const void *b) {
  const struct mpicache_entry *e1 = (const struct mpicache_entry *)a;
  const struct mpicache_entry *e2 = (const struct mpicache_entry *)b;

  int comp = e1->node - e2->node;
  if (comp < 0) return -1;
  if (comp > 0) return 1;

  /* Same node, check subtype. */
  comp = e1->task->subtype - e2->task->subtype;
  if (comp < 0) return -1;
  if (comp > 0) return 1;

  /* Same subtype, check tag. */
  return e1->task->flags - e2->task->flags;
}
#endif

/**
 * @brief reallocate the cache entries to make more space available.
 */
static void mpicache_reallocate(struct mpicache *cache, int ind) {

  if (ind == 0) {

    /* Need to perform initialization. */
    if ((cache->entries = (struct mpicache_entry *)malloc(
             CACHE_INC_SIZE * sizeof(struct mpicache_entry))) == NULL) {
      error("Failed to allocate MPI cache.");
    }

    /* Last action, releases threads. */
    cache->entries_size = CACHE_INC_SIZE;

  } else {
    struct mpicache_entry *new_entries;
    if ((new_entries = (struct mpicache_entry *)malloc(
             sizeof(struct mpicache_entry) *
             (cache->entries_size + CACHE_INC_SIZE))) == NULL) {
      error("Failed to re-allocate MPI cache.");
    }

    /* Wait for all writes to the old buffer to complete. */
    while (cache->entries_done < cache->entries_size)
      ;

    /* Copy to new buffer. */
    memcpy(new_entries, cache->entries,
           sizeof(struct mpicache_entry) * cache->entries_size);
    free(cache->entries);
    cache->entries = new_entries;

    /* Last action, releases waiting threads. */
    atomic_add(&cache->entries_size, CACHE_INC_SIZE);
  }
}

/**
 * @brief Initialize an mpicache for use.
 *
 * @return the new mpicache
 */
struct mpicache *mpicache_init() {

  /* Initialise the bitshifts needed to create the compact key. */
  if (mpicache_subtype_shift == 0) {
    mpicache_subtype_shift =
        (sizeof(int) * CHAR_BIT) - __builtin_clz(task_subtype_count);
  }

  /* Create an initial cache. */
  struct mpicache *cache =
      (struct mpicache *)calloc(1, sizeof(struct mpicache));
  return cache;
}

/**
 * @brief Add a new MPI task to the cache.
 *
 * @param cache the #mpicache
 * @param send_node the MPI rank that the message originates from.
 * @param t the #task struct.
 */
void mpicache_add(struct mpicache *cache, int node, struct task *t) {

  /* Append this to the cache. */
  size_t ind = atomic_inc(&(cache->nr_entries));
  if (ind == cache->entries_size) mpicache_reallocate(cache, ind);

  /* Wait if needed while the reallocation occurs. */
  while (ind > cache->entries_size)
    ;

  /* Derive the size in bytes, may not be set yet. */
#ifdef WITH_MPI
  t->size = scheduler_mpi_size(t);
#endif

  /* And store. */
  cache->entries[ind].node = node;
  cache->entries[ind].task = t;
  atomic_inc(&(cache->entries_done));

  return;
}

/**
 * @brief Free any data associated with the cache.
 *
 * @param cache the #mpicache
 */
void mpicache_destroy(struct mpicache *cache) {
  free(cache->entries);
  free(cache);
}

/**
 * @brief Apply the cache.
 *
 * The expected results are the updating of the MPI tasks so that their
 * offsets fields match the MPI window used and the expected sizes of the
 * windows are updated in the cache.
 *
 * @param cache the #mpicache.
 * @param nr_ranks number of MPI ranks. 
 */
void mpicache_apply(struct mpicache *cache, int nr_ranks, const char *prefix) {

#ifdef WITH_MPI
  /* Do nothing if we have nothing. */
  if (cache->nr_entries == 0) {
    return;
  }

  /* First job is to sort the entries to gives us the entries in
   * node|subtask|tag order. Within each node section the subtask|tag should
   * be the same on the send and recv sides, that is necessary so that the
   * offsets match. */
  qsort(cache->entries, cache->nr_entries, sizeof(struct mpicache_entry),
        mpicache_entry_cmp);

  /* Lists of offsets for each node in the subtype window. Keep this square
   * for convenience. */
  cache->window_node_offsets = (size_t *)malloc(task_subtype_count * nr_ranks * sizeof(size_t));
  for (int k = 0; k < task_subtype_count * nr_ranks; k++)
    cache->window_node_offsets[k] = LLONG_MAX;
  
  size_t offset[task_subtype_count] = {0};
  for (size_t k = 0; k < cache->nr_entries; k++) {

    struct task *task =  cache->entries[k].task;
    int node = cache->entries[k].node;
    int subtype = task->subtype;

#define INDEX2(nx, x, y) (nx * y + x)
    if (cache->window_node_offsets[INDEX2(task_subtype_count, subtype, node)] == LLONG_MAX) {

      /* Offset for this node and subtype not seen, so this is the first. */
      cache->window_node_offsets[INDEX2(task_subtype_count, subtype, node)] = offset[subtype];
      //message("node %d offset for subtype %d = %zu", node, subtype, offset[subtype]);

      /* And now first in this subtype once more. */
      offset[subtype] = 0;
    }

    /* Window sizes are in bytes. */
    cache->window_sizes[subtype] += task->size + scheduler_osmpi_tobytes(scheduler_osmpi_header_size);

    /* Offsets are in blocks. */
    task->offset = offset[subtype];
    offset[subtype] += scheduler_osmpi_toblocks(task->size) + scheduler_osmpi_header_size;

    //message("%s %d applied task %d subtype %d at %zd tag %lld node %d size %zu blocks %zu node offset %zu",
    //        prefix, task->ci->cellID, task->type, subtype, task->offset,
    //        task->flags, node, task->size,
    //        scheduler_osmpi_toblocks(task->size) + scheduler_osmpi_header_size,
    //        cache->window_node_offsets[INDEX2(task_subtype_count, subtype, node)]);
  }


#endif /* WITH_MPI */
}
