/*
 * Generic map implementation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../error.h"
#include "hashmap.h"

#define INITIAL_SIZE (1024)
#define MAX_CHAIN_LENGTH (8)
#define HASHMAP_GROWTH_FACTOR (2)

/**
 * brief: Pre-allocate a number of chunks for the graveyard.
 */
void hashmap_allocate_chunks(hashmap_t *m) {
  /* Allocate a fresh set of chunks. */
  hashmap_alloc_t *alloc;
  if ((alloc = (hashmap_alloc_t *)calloc(1, sizeof(hashmap_alloc_t))) == NULL) {
    error("Unable to allocate chunks.");
  }

  /* Hook up the alloc, so that we can clean it up later. */
  alloc->next = m->allocs;
  m->allocs = alloc;

  /* Link the chunks together. */
  for (int k = 0; k < CHUNKS_PER_ALLOC - 1; k++) {
    alloc->chunks[k].next = &alloc->chunks[k + 1];
  }

  /* Last chunk points to current graveyard. */
  alloc->chunks[CHUNKS_PER_ALLOC - 1].next = m->graveyard;

  /* Graveyard points to first new chunk. */
  m->graveyard = &alloc->chunks[0];
}

void hashmap_init(hashmap_t *m) {
  /* Allocate the first (empty) list of chunks. */
  m->nr_chunks = (INITIAL_SIZE + HASHMAP_ELEMENTS_PER_CHUNK - 1) / HASHMAP_ELEMENTS_PER_CHUNK;
  if ((m->chunks = (hashmap_chunk_t **)calloc(
           m->nr_chunks, sizeof(hashmap_chunk_t *))) == NULL) {
    error("Unable to allocate hashmap chunks.");
  }

  /* The graveyard is currently empty. */
  m->graveyard = NULL;
  m->allocs = NULL;

  /* Set initial sizes. */
  m->table_size = m->nr_chunks * HASHMAP_ELEMENTS_PER_CHUNK;
  m->size = 0;

  /* Inform the men. */
  message("Created hash table of size: %ld each element is %ld bytes. Allocated %d empty chunks.",
          INITIAL_SIZE * sizeof(hashmap_element_t), sizeof(hashmap_element_t), m->nr_chunks);
}

/**
 * @brief Put a used chunk back into the recycling bin.
 */
void hashmap_release_chunk(hashmap_t *m, hashmap_chunk_t *chunk) {
  /* Clear all the chunk's data. */
  memset(chunk, 0, sizeof(hashmap_chunk_t));

  /* Hook it up with the other stiffs in the graveyard. */
  chunk->next = m->graveyard;
  m->graveyard = chunk;
}

/**
 * @brief Return a new chunk, either recycled or freshly allocated.
 */
hashmap_chunk_t *hashmap_get_chunk(hashmap_t *m) {
  if (m->graveyard == NULL) {
    hashmap_allocate_chunks(m);
  }

  hashmap_chunk_t *res = m->graveyard;
  m->graveyard = res->next;
  res->next = NULL;

  return res;
}

/**
 * @brief Looks for the given key and retuns a pointer to the corresponding element.
 *
 * The returned element is either the one that already existed in the hashmap, or a
 * newly-reseverd element initialized to zero.
 *
 * If the hashmap is full, NULL is returned.
 *
 * We use `rand_r` as a hashing function. The key is first hashed to obtain an initial
 * global position. If there is a collision, the hashing function is re-applied to the
 * key to obtain a new offset *within the same bucket*. This is repeated for at most
 * MAX_CHAIN_LENGTH steps, at which point insertion fails.
 */
hashmap_element_t *hashmap_find(hashmap_t *m, size_t key, int create_new) {
  /* If full, return immediately */
  if (m->size >= (m->table_size / 2)) return NULL;

  /* We will use rand_r as our hash function. */
  unsigned int curr = (unsigned int)key;

  /* Get offsets to the entry, its chunk, it's mask, etc. */
  const int offset = rand_r(&curr) % m->table_size;
  const int chunk_offset = offset / HASHMAP_ELEMENTS_PER_CHUNK;
  int offset_in_chunk = offset - chunk_offset * HASHMAP_ELEMENTS_PER_CHUNK;

  /* Allocate the chunk if needed. */
  if (m->chunks[chunk_offset] == NULL) {
    /* Quit here if we don't want to create a new entry. */
    if (!create_new) return NULL;

    /* Get a new chunk for this offset. */
    m->chunks[chunk_offset] = hashmap_get_chunk(m);
  }
  hashmap_chunk_t *chunk = m->chunks[chunk_offset];

  /* Linear probing (well, not really, but whatever). */
  for (int i = 0; i < MAX_CHAIN_LENGTH; i++) {
    /* Compute the offsets within the masks of this chunk. */
    const int mask_offset = offset_in_chunk / HASHMAP_BITS_PER_MASK;
    const int offset_in_mask =
        offset_in_chunk - mask_offset * HASHMAP_BITS_PER_MASK;

    /* Is the offset empty? */
    hashmap_mask_t search_mask = ((hashmap_mask_t)1) << offset_in_mask;
    if (!(chunk->masks[mask_offset] & search_mask)) {
      /* Quit here if we don't want to create a new element. */
      if (!create_new) return NULL;

      /* Mark this element as taken and increase the size counter. */
      chunk->masks[mask_offset] |= search_mask;
      m->size += 1;

      /* Set the key. */
      chunk->data[offset_in_chunk].key = key;

      /* Return a pointer to the new element. */
      return &chunk->data[offset_in_chunk];
    }

    /* Does the offset by chance contain the key we are looking for? */
    else if (chunk->data[offset_in_chunk].key == key) {
      return &chunk->data[offset_in_chunk];
    }

    /* None of the above, so this is a collision. Re-hash, but within the same
       chunk. I guess this is Hopscotch Hashing? */
    else {
      offset_in_chunk = rand_r(&curr) % HASHMAP_ELEMENTS_PER_CHUNK;
    }
  }

  /* We lucked out, so return nothing. */
  return NULL;
}

/**
 * @brief Grows the hashmap and rehashes all the elements
 */
void hashmap_grow(hashmap_t *m) {
  /* Hold on to the old data. */
  const int old_table_size = m->table_size;
  hashmap_chunk_t **old_chunks = m->chunks;

  /* Re-allocate the chunk array. */
  m->table_size = HASHMAP_GROWTH_FACTOR * m->table_size;
  m->nr_chunks = m->table_size / HASHMAP_ELEMENTS_PER_CHUNK;
  if ((m->chunks = (hashmap_chunk_t **)calloc(
           m->nr_chunks, sizeof(hashmap_chunk_t *))) == NULL) {
    error("Unable to allocate hashmap chunks.");
  }

  /* Iterate over the chunks and add their entries to the new table. */
  for (int cid = 0; cid < old_table_size / HASHMAP_ELEMENTS_PER_CHUNK; cid++) {
    /* Skip empty chunks. */
    hashmap_chunk_t *chunk = old_chunks[cid];
    if (!chunk) continue;
    
    /* Loop over the masks in this chunk. */
    for (int mid = 0; mid < HASHMAP_MASKS_PER_CHUNK; mid++) {
      /* Skip empty masks. */
      if (chunk->masks[mid] == 0) continue;

      /* Loop over the mask entries. */
      for (int eid = 0; eid < HASHMAP_BITS_PER_MASK; eid++) {
        hashmap_mask_t element_mask = ((hashmap_mask_t)1) << eid;
        if (chunk->masks[mid] & element_mask) {
          hashmap_element_t *element =
              &chunk->data[mid * HASHMAP_BITS_PER_MASK + eid];

          /* Copy the element over to the new hashmap. */
          hashmap_element_t *new_element = hashmap_find(m, element->key, /*create_new=*/1);
          if (!new_element) {
            /* TODO(pedro): Deal with this type of failure more elegantly. */
            error("Failed to re-hash element.");
          }
          new_element->value = element->value;
        }
      }
    }

    /* We're through with this chunk, recycle it. */
    hashmap_release_chunk(m, chunk);
  }

  /* Free the old list of chunks. */
  free(old_chunks);
}

void hashmap_put(hashmap_t *m, size_t key, size_t value) {
  /* Try to find an element for the given key. */
  hashmap_element_t *element = hashmap_find(m, key, /*create_new=*/1);

  /* Loop around, trying to find our place in the world. */
  while (!element) {
    hashmap_grow(m);
    element = hashmap_find(m, key, /*create_new=*/1);
  }

  /* Set the value. */
  element->value = value;
}

size_t* hashmap_get(hashmap_t *m, size_t key) {
  /* Look for the given key. */
  hashmap_element_t *element = hashmap_find(m, key, /*create_new=*/1);
  while (!element) {
    hashmap_grow(m);
    element = hashmap_find(m, key, /*create_new=*/1);
  }
  return &element->value;
}

size_t* hashmap_lookup(hashmap_t *m, size_t key) {
  hashmap_element_t *element = hashmap_find(m, key, /*create_new=*/0);
  return element ? &element->value : NULL;
}

void hashmap_iterate(hashmap_t *m, hashmap_mapper_t f, void *data) {
  /* Loop over the chunks. */
  for (int cid = 0; cid < m->nr_chunks; cid++) {
    hashmap_chunk_t *chunk = m->chunks[cid];
    if (!chunk) continue;

    /* Loop over the masks. */
    for (int mid = 0; mid < HASHMAP_MASKS_PER_CHUNK; mid++) {
      hashmap_mask_t mask = chunk->masks[mid];
      if (!mask) continue;

      /* Loop over each element in the mask. */
      for (int eid = 0; 
           eid < HASHMAP_BITS_PER_MASK && 
               mid * HASHMAP_BITS_PER_MASK + eid < HASHMAP_ELEMENTS_PER_CHUNK; 
           eid++) {

        /* If the element exists, call the function on it. */
        if (mask & (((hashmap_mask_t)1) << eid)) {
          hashmap_element_t *element = &chunk->data[mid * HASHMAP_BITS_PER_MASK + eid];
          f(element->key, &element->value, data);
        }
      }
    }
  }
}

void hashmap_free(hashmap_t *m) {
  /* Free the list of active chunks. Note that the actual chunks will be freed
     as part of the allocs below. */
  free(m->chunks);

  /* Re-set some pointers and values, just in case. */
  m->chunks = NULL;
  m->graveyard = NULL;
  m->size = 0;
  m->table_size = 0;
  m->nr_chunks = 0;
 
  /* Free the chunk allocs. */ 
  while (m->allocs) {
    hashmap_alloc_t *tmp = m->allocs;
    m->allocs = tmp->next;
    free(tmp);
  }
}

int hashmap_size(hashmap_t *m) {
  if (m != NULL)
    return m->size;
  else
    return 0;
}
