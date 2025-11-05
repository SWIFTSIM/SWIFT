/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#include "lock.h"
#include "threadpool.h"

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

size_t particle_buffer_memory_use(struct particle_buffer *buffer);

#endif /* SWIFT_PARTICLE_BUFFER_H */
