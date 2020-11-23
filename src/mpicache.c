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

/* Comparator function to sort in node, subtype, tag order. */
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
 * @param node the MPI rank that the message originates.
 * @param t the #task struct.
 */
void mpicache_add(struct mpicache *cache, int node, struct task *t) {

  /* Append this to the cache. */
  size_t ind = atomic_inc(&(cache->nr_entries));
  if (ind == cache->entries_size) mpicache_reallocate(cache, ind);

  /* Wait if needed while the reallocation occurs. */
  while (ind > cache->entries_size)
    ;

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
 */
void mpicache_apply(struct mpicache *cache) {

#ifdef WITH_MPI

  /* First job is to sort the entries to gives us the entries in
   * node|subtask|tag order. Within each node section the subtask|tag should
   * be the same on the send and recv sides, that is necessary so that the
   * offsets match. */
  qsort(cache->entries, cache->nr_entries, sizeof(struct mpicache_entry),
        mpicache_entry_cmp);

  /* Now we go through the entries and generate the offsets. */
  int node = -1;
  int subtype = -1;
  cache->window_sizes = calloc(26, sizeof(int));  // 26 *task_subtype_count?
  cache->window_nodes = calloc(26, sizeof(int));
  cache->window_subtypes = calloc(26, sizeof(int));
  cache->window_size = 26;
  cache->nr_windows = 0;
  for (size_t k = 0; k < cache->nr_entries; k++) {

    /* New node started, so we loop until the start of the next one. */
    node = cache->entries[k].node;
    for (; k < cache->nr_entries; k++) {
      if (cache->entries[k].node != node) break;

      /* New subtype started, so we start a new set of offsets. */
      subtype = cache->entries[k].task->subtype;

      if (cache->nr_windows == cache->window_size) {
        cache->window_size += 26;
        cache->window_sizes =
            realloc(cache->window_sizes, cache->window_size * sizeof(int));
        cache->window_nodes =
            realloc(cache->window_nodes, cache->window_size * sizeof(int));
        cache->window_subtypes =
            realloc(cache->window_subtypes, cache->window_size * sizeof(int));
      }
      cache->window_sizes[cache->nr_windows] = 0;
      cache->window_nodes[cache->nr_windows] = node;
      cache->window_subtypes[cache->nr_windows] = subtype;

      size_t offset = 0;
      for (; k < cache->nr_entries; k++) {
        if (cache->entries[k].task->subtype != subtype) break;

        /* Offsets are in osmpi blocks, but window sizes are in bytes. */
        cache->window_sizes[cache->nr_windows] +=
            cache->entries[k].task->size +
            scheduler_osmpi_tobytes(scheduler_osmpi_header_size);
        cache->entries[k].task->offset = offset;

        /* Make room for this message and the control header next loop. */
        offset += scheduler_osmpi_toblocks(cache->entries[k].task->size) +
                  scheduler_osmpi_header_size;
      }
      cache->nr_windows++;
    }
  }
#endif /* WITH_MPI */
}
