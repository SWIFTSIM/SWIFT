/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 James Willis (james.s.willis@durham.ac.uk)
 *                    Pedro Gonnet (pedro.gonnet@gmail.com)
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
/*
 * Generic hashmap manipulation functions
 *
 * Originally by Elliot C Back -
 * http://elliottback.com/wp/hashmap-implementation-in-c/
 *
 * Modified by Pete Warden to fix a serious performance problem, support strings
 * as keys and removed thread synchronization - http://petewarden.typepad.com
 */
#ifndef SWIFT_HASHMAP_UNION_FIND_H
#define SWIFT_HASHMAP_UNION_FIND_H

/* Some standard headers. */
#include <stdbool.h>
#include <stddef.h>

/* Local headers. */
#include "align.h"

// Type used for chunk bitmasks.
typedef size_t hashmap_uf_mask_t;

#define HASHMAP_UF_BITS_PER_MASK ((int)sizeof(hashmap_uf_mask_t) * 8)

// Type used for the hashmap_uf keys (must have a valid '==' operation).
#ifndef hashmap_uf_key_t
#define hashmap_uf_key_t size_t
#endif

// Type used for the hashmap_uf values (must have a valid '==' operation).
#ifndef hashmap_uf_value_t
typedef struct _hashmap_uf_struct {
  long long value_st;
} hashmap_uf_struct_t;
#define hashmap_uf_value_t hashmap_uf_struct_t
#endif

/* We need to keep keys and values */
typedef struct _hashmap_uf_element {
  hashmap_uf_key_t key;
  hashmap_uf_value_t value;
} hashmap_uf_element_t;

/* Make sure a chunk fits in a given size. */
#define HASHMAP_UF_TARGET_CHUNK_BYTES (80 * 1024)
#define HASHMAP_UF_BITS_PER_ELEMENT ((int)sizeof(hashmap_uf_element_t) * 8 + 1)
#define HASHMAP_UF_ELEMENTS_PER_CHUNK \
  ((HASHMAP_UF_TARGET_CHUNK_BYTES * 8) / HASHMAP_UF_BITS_PER_ELEMENT)
#define HASHMAP_UF_MASKS_PER_CHUNK                                  \
  ((HASHMAP_UF_ELEMENTS_PER_CHUNK + HASHMAP_UF_BITS_PER_MASK - 1) / \
   HASHMAP_UF_BITS_PER_MASK)

#define HASHMAP_UF_ALLOCS_INITIAL_SIZE (80 * 1024)
#define HASHMAP_UF_ALLOC_SIZE_FRACTION (0.1)

#define HASHMAP_UF_MAX_CHAIN_LENGTH (HASHMAP_UF_ELEMENTS_PER_CHUNK / 8)
#ifndef HASHMAP_UF_DEBUG_OUTPUT
#define HASHMAP_UF_DEBUG_OUTPUT (0)
#endif  // HASHMAP_UF_DEBUG_OUTPUT

/* A chunk of hashmap_uf_element, with the corresponding bitmask. */
typedef struct _hashmap_uf_chunk {
  union {
    hashmap_uf_mask_t masks[HASHMAP_UF_MASKS_PER_CHUNK];
    void *next;
  };
  hashmap_uf_element_t data[HASHMAP_UF_ELEMENTS_PER_CHUNK];
} SWIFT_STRUCT_ALIGN hashmap_uf_chunk_t;

/* A hashmap_uf has some maximum size and current size,
 * as well as the data to hold. */
typedef struct _hashmap_uf {
  size_t table_size;
  size_t size;
  size_t nr_chunks;
  hashmap_uf_chunk_t *
      *chunks;  // Pointer to chunks in use, but not densely populated.
  hashmap_uf_chunk_t
      *graveyard;  // Pointer to allocated, but currently unused chunks.

  void **allocs;        // Pointers to allocated blocks of chunks.
  size_t allocs_size;   // Size of the allocs array.
  size_t allocs_count;  // Number of elements in the allocs array.

#if HASHMAP_UF_DEBUG_OUTPUT
  /* Chain lengths, used for debugging only. */
  size_t chain_length_counts[HASHMAP_MAX_CHAIN_LENGTH];
#endif
} hashmap_uf_t;

/**
 * Pointer to a function that can take a key, a pointer to a value, and a
 * void pointer extra data payload.
 */
typedef void (*hashmap_uf_mapper_t)(hashmap_uf_key_t, hashmap_uf_value_t *,
                                    void *);

/**
 * @brief Initialize a hashmap.
 */
void hashmap_uf_init(hashmap_uf_t *m);

/**
 * @brief Re-size the hashmap.
 *
 * Note that the hashmap size does not necessarily correspond to its
 * capacity, since it will grow if too many collisions occur. As a rule
 * of thumb, allocate twice as many elements as you think you will need.
 *
 * @param m The hasmmap to grow.
 * @param new_size New table size. If zero, the current size will be increase
 *                 by a fixed rate.
 */
void hashmap_uf_grow(hashmap_uf_t *m, size_t new_size);

/**
 * @brief Add a key/value pair to the hashmap, overwriting whatever was
 * previously there.
 */
extern void hashmap_uf_put(hashmap_uf_t *m, hashmap_uf_key_t key,
                           hashmap_uf_value_t value);

/**
 * @brief Get the value for a given key. If no value exists a new one will be
 * created.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_uf_value_t *hashmap_uf_get(hashmap_uf_t *m,
                                          hashmap_uf_key_t key);

/**
 * @brief Get the value for a given key. If no value exists a new one will be
 * created. Return a flag indicating whether a new element has been added.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_uf_value_t *hashmap_uf_get_new(hashmap_uf_t *m,
                                              hashmap_uf_key_t key,
                                              int *created_new_element);

/**
 * @brief Look for the given key and return a pointer to its value or NULL if
 * it is not in the hashmap.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_uf_value_t *hashmap_uf_lookup(hashmap_uf_t *m,
                                             hashmap_uf_key_t key);

/**
 * @brief Iterate the function parameter over each element in the hashmap.
 *
 * The function `f` takes three arguments, the first and second are the element
 * key and a pointer to the correspondig value, respectively, while the third
 * is the `void *data` argument.
 */
extern void hashmap_uf_iterate(hashmap_uf_t *m, hashmap_uf_mapper_t f,
                               void *data);

/**
 * @brief De-allocate memory associated with this hashmap, clears all the
 * entries.
 *
 * After a call to `hashmap_free`, the hashmap cna be re-initialized with
 * `hashmap_init`.
 */
extern void hashmap_uf_free(hashmap_uf_t *m);

/**
 * Get the current size of a hashmap
 */
extern size_t hashmap_uf_size(hashmap_uf_t *m);

/**
 * @brief Print all sorts of stats on the given hashmap.
 */
void hashmap_uf_print_stats(hashmap_uf_t *m);

#endif /* SWIFT_HASHMAP_UNION_FIND_H */
