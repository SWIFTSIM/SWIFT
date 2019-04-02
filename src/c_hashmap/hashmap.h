/*
 * Generic hashmap manipulation functions
 *
 * Originally by Elliot C Back -
 * http://elliottback.com/wp/hashmap-implementation-in-c/
 *
 * Modified by Pete Warden to fix a serious performance problem, support strings
 * as keys and removed thread synchronization - http://petewarden.typepad.com
 */
#ifndef __HASHMAP_H__
#define __HASHMAP_H__

/* Some standard headers. */
#include <stdbool.h>
#include <stddef.h>

// Type used for chunk bitmasks.
typedef size_t hashmap_mask_t;

#define HASHMAP_BITS_PER_MASK ((int)sizeof(hashmap_mask_t) * 8)

// Type used for the hashmap keys (must have a valid '==' operation).
#ifndef hashmap_key_t
#define hashmap_key_t size_t
#endif

// Type used for the hashmap values (must have a valid '==' operation).
#ifndef hashmap_value_t
typedef union _hashmap_union {
  size_t value_st;
  double value_dbl;
} hashmap_union_t;
#define hashmap_value_t hashmap_union_t
#endif

/* We need to keep keys and values */
typedef struct _hashmap_element {
  hashmap_key_t key;
  hashmap_value_t value;
} hashmap_element_t;

/* Make sure a chunk fits in a given size. */
#define HASHMAP_TARGET_CHUNK_BYTES (4 * 1024)
#define HASHMAP_BITS_PER_ELEMENT ((int)sizeof(hashmap_element_t) * 8 + 1)
#define HASHMAP_ELEMENTS_PER_CHUNK \
  ((HASHMAP_TARGET_CHUNK_BYTES * 8) / HASHMAP_BITS_PER_ELEMENT)
#define HASHMAP_MASKS_PER_CHUNK                               \
  ((HASHMAP_ELEMENTS_PER_CHUNK + HASHMAP_BITS_PER_MASK - 1) / \
   HASHMAP_BITS_PER_MASK)

#define HASHMAP_ALLOCS_INITIAL_SIZE (8)
#define HASHMAP_ALLOC_SIZE_FRACTION (0.1)

#define HASHMAP_MAX_CHAIN_LENGTH (HASHMAP_ELEMENTS_PER_CHUNK / 8)
#ifndef HASHMAP_DEBUG_OUTPUT
#define HASHMAP_DEBUG_OUTPUT (0)
#endif  // HASHMAP_DEBUG_OUTPUT

/* A chunk of hashmap_element, with the corresponding bitmask. */
typedef struct _hashmap_chunk {
  union {
    hashmap_mask_t masks[HASHMAP_MASKS_PER_CHUNK];
    void *next;
  };
  hashmap_element_t data[HASHMAP_ELEMENTS_PER_CHUNK];
} hashmap_chunk_t;

/* A hashmap has some maximum size and current size,
 * as well as the data to hold. */
typedef struct _hashmap {
  size_t table_size;
  size_t size;
  size_t nr_chunks;
  hashmap_chunk_t *
      *chunks;  // Pointer to chunks in use, but not densely populated.
  hashmap_chunk_t
      *graveyard;  // Pointer to allocated, but currently unused chunks.

  void **allocs;        // Pointers to allocated blocks of chunks.
  size_t allocs_size;   // Size of the allocs array.
  size_t allocs_count;  // Number of elements in the allocs array.

#if HASHMAP_DEBUG_OUTPUT
  /* Chain lengths, used for debugging only. */
  size_t chain_length_counts[HASHMAP_MAX_CHAIN_LENGTH];
#endif
} hashmap_t;

/**
 * Pointer to a function that can take a key, a pointer to a value, and a
 * void pointer extra data payload.
 */
typedef void (*hashmap_mapper_t)(hashmap_key_t, hashmap_value_t *, void *);

/**
 * @brief Initialize a hashmap.
 */
void hashmap_init(hashmap_t *m);

/**
 * @brief Add a key/value pair to the hashmap, overwriting whatever was
 * previously there.
 */
extern void hashmap_put(hashmap_t *m, hashmap_key_t key, hashmap_value_t value);

/**
 * @brief Get the value for a given key. If no value exists a new one will be
 * created.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_value_t *hashmap_get(hashmap_t *m, hashmap_key_t key);

/**
 * @brief Get the value for a given key. If no value exists a new one will be
 * created. Return a flag indicating whether a new element has been added.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_value_t *hashmap_get_new(hashmap_t *m, hashmap_key_t key, int *created_new_element);

/**
 * @brief Look for the given key and return a pointer to its value or NULL if
 * it is not in the hashmap.
 *
 * Note that the returned pointer is volatile and will be invalidated if the
 * hashmap is re-hashed!
 */
extern hashmap_value_t *hashmap_lookup(hashmap_t *m, hashmap_key_t key);

/**
 * @brief Iterate the function parameter over each element in the hashmap.
 *
 * The function `f` takes three arguments, the first and second are the element
 * key and a pointer to the correspondig value, respectively, while the third
 * is the `void *data` argument.
 */
extern void hashmap_iterate(hashmap_t *m, hashmap_mapper_t f, void *data);

/**
 * @brief De-allocate memory associated with this hashmap, clears all the
 * entries.
 *
 * After a call to `hashmap_free`, the hashmap cna be re-initialized with
 * `hashmap_init`.
 */
extern void hashmap_free(hashmap_t *m);

/**
 * Get the current size of a hashmap
 */
extern size_t hashmap_size(hashmap_t *m);

/**
 * @brief Print all sorts of stats on the given hashmap.
 */
void hashmap_print_stats(hashmap_t *m);

#endif /*__HASHMAP_H__*/
