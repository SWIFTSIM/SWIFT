#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdio.h>

#include "particle_buffer.h"

/**
 * @brief Initialize a particle buffer.
 *
 * Stores a sequence of data objects of size element_size,
 * allowing new elements to be appended by multiple threads
 * simultaneously.
 *
 * Objects are stored in a linked list of blocks and new blocks
 * are allocated as needed.
 *
 * @param buffer The #particle_buffer
 * @param element_size Size of a single element
 * @param elements_per_block Number of elements to store in each block
 *
 */
void particle_buffer_init(struct particle_buffer *buffer, size_t element_size,
                          size_t elements_per_block) {
  
  buffer->element_size = element_size;
  buffer->elements_per_block = elements_per_block;
  buffer->first_block = malloc(sizeof(struct particle_buffer_block));
  buffer->first_block->num_elements = 0;
  buffer->first_block->data = malloc(elements_per_block*element_size);
  buffer->first_block->next = NULL;
  buffer->last_block = buffer->first_block;
  lock_init(&buffer->lock);
}


/**
 * @brief Deallocate a particle buffer.
 *
 * @param buffer The #particle_buffer
 *
 */
void particle_buffer_free(struct particle_buffer *buffer) {
  
  struct particle_buffer_block *block = buffer->first_block;
  while(block) {
    struct particle_buffer_block *next = block->next;
    free(block->data);
    free(block);
    block = next;
  }
  buffer->first_block = NULL;
  buffer->last_block = NULL;
  lock_destroy(&buffer->lock);
}


/**
 * @brief Append an element to a particle buffer.
 *
 * @param buffer The #particle_buffer
 * @param data The element to append
 *
 */
void particle_buffer_append(struct particle_buffer *buffer, void *data) {
 
  const size_t element_size = buffer->element_size;
  const size_t elements_per_block = buffer->elements_per_block;
 
  while(1) {

    /* Find the current block */
    struct particle_buffer_block *block = buffer->last_block;

    /* Find next available index in current block */
    size_t index = __atomic_fetch_add(&block->num_elements, 1, __ATOMIC_SEQ_CST);

    /* Increment total number of elements */
    __atomic_fetch_add(&buffer->total_num_elements, 1, __ATOMIC_SEQ_CST);

    if(index < buffer->elements_per_block) {
      /* We reserved a valid index, so copy the data */
      memcpy(block->data+index*buffer->element_size, data, buffer->element_size);
      return;
    } else {
      /* No space left, so we need to allocate a new block */
      lock_lock(&buffer->lock);
      /* Check no-one else already did it before we got the lock */
      if(!block->next) {
        /* Allocate and initialize the new block */
        struct particle_buffer_block *new_block = malloc(sizeof(struct particle_buffer_block));
        char *new_block_data = malloc(elements_per_block*element_size);
        __atomic_store_n(&new_block->data, new_block_data, __ATOMIC_SEQ_CST);
        __atomic_store_n(&new_block->num_elements, 0, __ATOMIC_SEQ_CST);
        __atomic_store_n(&new_block->next, NULL, __ATOMIC_SEQ_CST);
        __atomic_store_n(&block->next, new_block, __ATOMIC_SEQ_CST);
        __atomic_thread_fence(__ATOMIC_SEQ_CST);
        /* After this store other threads will write to the new block,
           so all initialization must complete before this. */
        __atomic_store_n(&buffer->last_block, new_block, __ATOMIC_SEQ_CST);
      }
      lock_unlock(&buffer->lock);
      /* Now we have space, will try again */
    }
  }
}


/**
 * @brief Iterate over data blocks in particle buffer.
 *
 * @param buffer The #particle_buffer
 * @param block Initially null, returns pointer to next data block
 * @param num_elements Returns number of elements in this block
 * @param data Returns pointer to data in this block
 *
 */
void particle_buffer_iterate(struct particle_buffer *buffer,
                             struct particle_buffer_block **block,
                             size_t *num_elements, void **data) {

  if(!*block) {
    *block = buffer->first_block;
  } else {
    *block = (*block)->next;
  }

  if(*block) {
    *data = (*block)->data;
    *num_elements = (*block)->num_elements;
    if(*num_elements > buffer->elements_per_block)
      *num_elements = buffer->elements_per_block;
  } else {
    *data = NULL;
    *num_elements = 0;
  }

}
