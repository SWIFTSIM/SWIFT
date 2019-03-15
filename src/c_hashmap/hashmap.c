/*
 * Generic map implementation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../error.h"
#include "hashmap.h"

#define INITIAL_SIZE (4096)
#define MAX_CHAIN_LENGTH (8)
#define CHUNKS_PER_ALLOC 8
#define HASHMAP_GROWTH_FACTOR 2

/**
 * brief: Pre-allocate a number of chunks for the graveyard.
 */
void hashmap_allocate_chunks(hashmap_t *m, int num_chunks) {
  /* Allocate a fresh set of chunks. */
  hashmap_chunk_t *chunks;
  if ((chunks = (hashmap_chunk_t *)calloc(num_chunks,
                                          sizeof(hashmap_chunk_t))) == NULL) {
    error("Unable to allocate chunks.");
  }

  /* Link the chunks together. */
  for (int k = 0; k < num_chunks - 1; k++) {
    chunks[k].next = &chunks[k + 1];
  }

  /* Last chunk points to current graveyard. */
  chunks[num_chunks - 1].next = m->graveyard;

  /* Graveyard points to first new chunk. */
  m->graveyard = chunks;
}

/**
 * @brief Initialize a hashmap.
 */
void hashmap_init(hashmap_t *m) {
  /* Allocate the first (empty) list of chunks. */
  const int nr_chunks = INITIAL_SIZE / HASHMAP_ELEMENTS_PER_CHUNK;
  if ((m->chunks = (hashmap_chunk_t **)calloc(
           nr_chunks, sizeof(hashmap_chunk_t *))) == NULL) {
    error("Unable to allocate hashmap chunks.");
  }

  /* The graveyard is currently empty. */
  m->graveyard = NULL;

  /* Set initial sizes. */
  m->table_size = INITIAL_SIZE;
  m->size = 0;

  /* Inform the men. */
  message("Created hash table of size: %ld each element is %ld bytes. Allocated %d empty chunks.",
          INITIAL_SIZE * sizeof(hashmap_element_t), sizeof(hashmap_element_t), nr_chunks);
}

/**
 * @brief Put a used chunk back into the recycling bin.
 */
void hashmap_release_chunk(hashmap_t *m, hashmap_chunk_t *chunk) {
  /* Clear all the chunk's data. */
  memset(chunk, 0, sizeof(hashmap_chunk_t));

  /* Hook it up with the other stiffs in the graveyard. */
  chunk->next = m->graveyard;
  m->graveyard = chunk->next;
}

/**
 * @brief Return a new chunk, either recycled or freshly allocated.
 */
hashmap_chunk_t *hashmap_get_chunk(hashmap_t *m) {
  if (m->graveyard == NULL) {
    hashmap_allocate_chunks(m, CHUNKS_PER_ALLOC);
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

  /* Linear probing (well, not really). */
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
  const int nr_chunks = m->table_size / HASHMAP_ELEMENTS_PER_CHUNK;
  if ((m->chunks = (hashmap_chunk_t **)calloc(
           nr_chunks, sizeof(hashmap_chunk_t *))) == NULL) {
    error("Unable to allocate hashmap chunks.");
  }

  /* Iterate over the chunks and add their entries to the new table. */
  for (int cid = 0; cid < old_table_size / HASHMAP_ELEMENTS_PER_CHUNK; cid++) {
    /* Skip empty chunks. */
    if (old_chunks[cid] == NULL) continue;

    hashmap_chunk_t *chunk = old_chunks[cid];

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
}

/**
 * @brief Add a key/value pair to the hashmap, overwriting whatever was previously there.
 */
void hashmap_put(hashmap_t *m, size_t key, size_t value) {
    hashmap_element_t *element = hashmap_find(m, key, /*create_new=*/1);

  /* Loop around, trying to find our place in the world. */
  while (!element) {
    hashmap_grow(m);
    element = hashmap_find(m, key, /*create_new=*/1);
  }

  /* Set the value. */
  element->value = value;
}

/**
 * @brief Get the value for a given key. If no value exists and create_new is true,
 * a new one will be created.
 *
 * Note that the returned pointer is volatile and will be invalidated if the hashmap
 * is re-hashed!
 */
size_t* hashmap_get(hashmap_t *m, size_t key, int create_new) {
  /* Look for the given key. */
  hashmap_element_t *element = hashmap_find(m, key, create_new);

  /* If we weren't asked to create a new entry, return NULL if the key was not found. */
  if (!create_new) {
    return element ? &element->value : NULL;
  }

  /* Otherwise, re-hash and repeat. */
  else {
    while (!element) {
      hashmap_grow(m);
      element = hashmap_find(m, key, /*create_new=*/1);
    }
    return &element->value;
  }
}

/*
 * Iterate the function parameter over each element in the hashmap.  The
 * additional any_t argument is passed to the function as its first
 * argument and the hashmap element is the second.
 */
int hashmap_iterate(map_t in, PFany f, any_t item) {
//  int i;
//
//  /* Cast the hashmap */
//  hashmap_t *m = (hashmap_t *)in;
//
//  /* On empty hashmap, return immediately */
//  if (hashmap_length(m) <= 0) return MAP_MISSING;
//
//  /* Linear probing */
//  for (i = 0; i < m->table_size; i++)
//    if (m->chunks[i / 64].in_use & (1 << (i % 64))) {
//      any_t data = (any_t) & (m->chunks[i / 64].data[i % 64]);
//      int status = f(item, data);
//      if (status != MAP_OK) {
//        return status;
//      }
//    }
//
  return MAP_OK;
}

/*
 * Remove an element with that key from the map
 */
//int hashmap_remove(map_t in, size_t key) {
//  int i;
//  int curr;
//  hashmap_map *m;
//
//  /* Cast the hashmap */
//  m = (hashmap_map *)in;
//
//  /* Find key */
//  curr = hashmap_hash_int(m, key);
//
//  /* Linear probing, if necessary */
//  for (i = 0; i < MAX_CHAIN_LENGTH; i++) {
//    int in_use = m->chunks[curr / 64].in_use & (1 << (curr % 64));
//    if (in_use == 1) {
//      if (m->chunks[curr / 64].data[curr % 64].key == key) {
//        /* Blank out the fields */
//        m->chunks[curr / 64].in_use &= ~(1 << (curr % 64));
//        m->chunks[curr / 64].data[curr % 64].key = 0;
//        m->chunks[curr / 64].data[curr % 64].group_size = 0;
//
//        /* Reduce the size */
//        m->size--;
//        return MAP_OK;
//      }
//    }
//    curr = (curr + 1) % m->table_size;
//  }
//
//  /* Data not found */
//  return MAP_MISSING;
//}

/* Deallocate the hashmap */
void hashmap_free(hashmap_t *m) {
  free(m->chunks);
 
  /* Free chunks */ 
  hashmap_chunk_t *head = m->graveyard;
  hashmap_chunk_t *tmp;
  while(head != NULL) {
  
    tmp = head;
    head = head->next;
    free(tmp);

  }
}

/* Return the length of the hashmap */
int hashmap_length(hashmap_t *m) {
  if (m != NULL)
    return m->size;
  else
    return 0;
}
