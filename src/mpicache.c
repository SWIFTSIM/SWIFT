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
 *  @brief file of routines to cache MPI task messages that arrive without a
 *         ready recv task.
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

/* Local includes. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"
#include "memuse_rnodes.h"
#include "task.h"

#define CHECKS 1

/* A single cache entry. */
struct mpicache_entry {
  size_t size;
  void *data;
#ifdef CHECKS
  int rank;
  int subtype;
  int tag;
#endif
};

/* Mapping a key into an 8 byte pointer. */
union key {
  size_t keyval;
  uint8_t ptr[sizeof(size_t)];
};

/* Bit shift to accomodate all the bits of the maximum rank id. */
static int rank_shift = 0;

/* Bit shift to accomodate all the bits of the maximum subtype. */
static int subtype_shift = 0;

/**
 * @brief Convert ranks, subtype and tag into a unique compact key.
 *
 * @param rank the MPI rank that the message came from.
 * @param subtype the task subtype of the message.
 * @param tag the message tag.
 *
 * @result the key for this information.
 * 
 * XXX do we also need size, if so need a word twice as large.
 */
size_t mpicache_makekey(int rank, int subtype, int tag) {

  /* Compact all these values into a single size_t word. */
  size_t result = subtype | engine_rank << subtype_shift |
                  rank << (subtype_shift + rank_shift) |
                  tag << (subtype_shift * 2 + rank_shift);
  return result;
}

/**
 * @brief Initialize an mpicache for use.
 *
 * @param nr_ranks the number of MPI ranks.
 *
 * @return the new mpicache
 */
struct mpicache *mpicache_init(int nr_ranks) {

  /* Create the root of the cache. */
  struct mpicache *cache = (struct mpicache *)
    calloc(1, sizeof(struct mpicache));
  cache->root.value = -1;

  /* Initialise the bitshifts needed to create the compact key. */
  if (rank_shift == 0) {
    rank_shift = (sizeof(int) * CHAR_BIT) - __builtin_clz(nr_ranks);
    subtype_shift = (sizeof(int) * CHAR_BIT) - __builtin_clz(task_subtype_count);
  }

  return cache;
}

/**
 * @brief Add a new MPI messages to the cache.
 *
 * @param cache the #mpicache
 * @param rank the MPI rank that the message came from.
 * @param subtype the task subtype of the message.
 * @param tag the message tag.
 * @param size size of message in bytes.
 * @param data the message data.
 */
void mpicache_add(struct mpicache *cache, int rank, int subtype, int tag,
                  size_t size, void *data) {

  union key key;
  key.keyval = mpicache_makekey(rank, subtype, tag);

  /* Check if this record is already present. */
  struct memuse_rnode *child =
    memuse_rnode_find_child(&cache->root, 0, key.ptr, sizeof(size_t));
  if (child != NULL && child->value != -1) {

    /* Already present, this is an error. */
    error("Attempt to add existing MPI message to cache");
  }

  struct mpicache_entry *entry = (struct mpicache_entry *)
    calloc(1, sizeof(struct mpicache_entry));
  entry->data = data;
  entry->size = size;
#ifdef CHECKS
  entry->rank = rank;
  entry->subtype = subtype;
  entry->tag = tag;
#endif

  memuse_rnode_insert_child(&cache->root, 0, key.ptr,
                            sizeof(size_t), (int64_t)entry);
  return;
}

/**
 * @brief Fetch an MPI messages from the cache, if it exists. Once found the
 *        message is released.
 *
 * @param cache the #mpicache
 * @param rank the MPI rank that the message is expected from
 * @param subtype the task subtype of the message.
 * @param tag the message tag.
 * @param size size of message in bytes if found, zero it not.
 * @param data the message data.
 */
void mpicache_fetch(struct mpicache *cache, int rank, int subtype, int tag,
                    size_t *size, void **data) {

  union key key;
  key.keyval = mpicache_makekey(rank, subtype, tag);

  /* Check if this record is already present. */
  struct memuse_rnode *child =
    memuse_rnode_find_child(&cache->root, 0, key.ptr, sizeof(size_t));
  if (child != NULL && child->value != -1) {

    struct mpicache_entry *entry = (struct mpicache_entry *) child->value;
    *size = entry->size;
    *data = entry->data;
#ifdef CHECKS
    if (entry->rank != rank || entry->subtype != subtype || entry->tag != tag)
      error("Fetched the wrong entry");
#endif

    free(entry);
    child->value = -1;
  } else {
    *size = 0;
    *data = NULL;
  }
  return;
}

/**
 * @brief Free any data associated with the cache.
 *
 * @param cache the #mpicache
 */
void mpicache_destroy(struct mpicache *cache) {
  free(cache);
}

