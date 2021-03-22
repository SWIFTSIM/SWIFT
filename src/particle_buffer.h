#include "lock.h"

#ifndef SWIFT_PARTICLE_BUFFER_H
#define SWIFT_PARTICLE_BUFFER_H

#define PARTICLE_BUFFER_NAME_LENGTH 100

struct particle_buffer_block {
  size_t num_elements;
  char *data;
  struct particle_buffer_block *next;
};

struct particle_buffer {
  size_t element_size;
  size_t elements_per_block;
  struct particle_buffer_block *first_block;
  struct particle_buffer_block *last_block;
  swift_lock_type lock;
  char name[PARTICLE_BUFFER_NAME_LENGTH];
};


void particle_buffer_init(struct particle_buffer *buffer, size_t element_size,
                          size_t elements_per_block, char *name);

void particle_buffer_free(struct particle_buffer *buffer);

void particle_buffer_empty(struct particle_buffer *buffer);

void particle_buffer_append(struct particle_buffer *buffer, void *data);

void particle_buffer_iterate(struct particle_buffer *buffer,
                             struct particle_buffer_block **block,
                             size_t *num_elements, void **data);

size_t particle_buffer_num_elements(struct particle_buffer *buffer);

#endif /* SWIFT_PARTICLE_BUFFER_H */
