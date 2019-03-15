/*
 * Generic hashmap manipulation functions
 *
 * Originally by Elliot C Back - http://elliottback.com/wp/hashmap-implementation-in-c/
 *
 * Modified by Pete Warden to fix a serious performance problem, support strings as keys
 * and removed thread synchronization - http://petewarden.typepad.com
 */
#ifndef __HASHMAP_H__
#define __HASHMAP_H__

/* Some standard headers. */
#include <stddef.h>
#include <stdbool.h>

// Type used for chunk bitmasks.
typedef size_t hashmap_mask_t;

#define HASHMAP_BITS_PER_MASK ((int)sizeof(hashmap_mask_t) * 8)
#define HASHMAP_MASKS_PER_CHUNK 4
#define HASHMAP_ELEMENTS_PER_CHUNK (HASHMAP_BITS_PER_MASK * HASHMAP_MASKS_PER_CHUNK)
#define CHUNKS_PER_ALLOC 8

/* We need to keep keys and values */
typedef struct _hashmap_element{
	size_t key;
	size_t value;
} hashmap_element_t;

/* A chunk of hashmap_element, with the corresponding bitmask. */
typedef struct _hashmap_chunk {
	union {
		hashmap_mask_t masks[HASHMAP_MASKS_PER_CHUNK];
		void *next;
	};
	hashmap_element_t data[HASHMAP_ELEMENTS_PER_CHUNK];
} hashmap_chunk_t;

struct _hashmap_alloc {
  hashmap_chunk chunks[CHUNKS_PER_ALLOC];
  void *next;
} hashmap_alloc_t;

/* A hashmap has some maximum size and current size,
 * as well as the data to hold. */
typedef struct _hashmap_map{
	int table_size;
	int size;
	int num_chunks;
	hashmap_chunk_t **chunks;    // Pointer to chunks in use, but not densely populated.
	hashmap_chunk_t *graveyard;  // Pointer to allocated, but currently unused chunks.
    hashmap_alloc_t *allocs;	// Pointer to the allocated chunks of chunks, needed for cleanup.
} hashmap_t;

/*
 * any_t is a pointer.  This allows you to put arbitrary structures in
 * the hashmap.
 */
typedef void *any_t;

/*
 * PFany is a pointer to a function that can take two any_t arguments
 * and return an integer. Returns status code..
 */
typedef int (*PFany)(any_t, any_t);

/*
 * map_t is a pointer to an internally maintained data structure.
 * Clients of this package do not need to know how hashmaps are
 * represented.  They see and manipulate only map_t's.
 */
typedef any_t map_t;

/**
 * @brief Initialize a hashmap.
 */
void hashmap_init(hashmap_t *m);

/*
 * Iteratively call f with argument (item, data) for
 * each element data in the hashmap. The function must
 * return a map status code. If it returns anything other
 * than MAP_OK the traversal is terminated. f must
 * not reenter any hashmap functions, or deadlock may arise.
 */
extern int hashmap_iterate(map_t in, PFany f, any_t item);

/*
 * Add an element to the hashmap. Return MAP_OK or MAP_OMEM.
 */
extern void hashmap_put(hashmap_t *m, size_t key, size_t value);

/*
 * Get an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern size_t* hashmap_get(hashmap_t *m, size_t key);
extern size_t* hashmap_lookup(hashmap_t *m, size_t key);

/*
 * Remove an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern int hashmap_remove(hashmap_t *m, size_t key);

/*
 * Get any element. Return MAP_OK or MAP_MISSING.
 * remove - should the element be removed from the hashmap
 */
extern int hashmap_get_one(map_t in, any_t *arg, int remove);

/*
 * Free the hashmap
 */
extern void hashmap_free(hashmap_t *m);

/*
 * Get the current size of a hashmap
 */
extern int hashmap_length(hashmap_t *m);

#endif /*__HASHMAP_H__*/
