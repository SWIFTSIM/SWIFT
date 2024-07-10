/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <limits.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "space_unique_id.h"

/* Local headers. */
#include "engine.h"
#include "lock.h"
#include "space.h"

/**
 * @brief Update the unique id structure by requesting a
 * new batch if required.
 *
 * @param s The #space.
 */
void space_update_unique_id(struct space *s) {
  /* Do we need unique IDs? */
  if (!star_formation_need_unique_id && !sink_need_unique_id) {
    return;
  }

  int require_new_batch = s->unique_id.next_batch.current == 0;

#ifdef WITH_MPI
  const struct engine *e = s->e;

  /* Check if the other ranks need a batch. */
  int *all_requires = (int *)malloc(sizeof(int) * e->nr_nodes);

  /* Do the communication */
  MPI_Allgather(&require_new_batch, 1, MPI_INT, all_requires, 1, MPI_INT,
                MPI_COMM_WORLD);

  /* Compute the position of this rank batch and the position of
     the next free batch. */
  int local_index = 0;
  int total_shift = 0;
  for (int i = 0; i < e->nr_nodes; i++) {
    total_shift += all_requires[i];
    if (i < e->nodeID) {
      local_index += all_requires[i];
    }
  }

  /* Free the allocated resources. */
  free(all_requires);

#else

  int local_index = 0;
  int total_shift = require_new_batch;

#endif  // WITH_MPI

  /* Compute the size of each batch. */
  const long long batch_size = (space_extra_parts + space_extra_sparts +
                                space_extra_gparts + space_extra_bparts) *
                               s->nr_cells;

  /* Get a new batch. */
  if (require_new_batch) {
    /* First check against an overflow. */
    const long long local_shift = local_index * batch_size;
    if (s->unique_id.global_next_id > LLONG_MAX - (local_shift + batch_size)) {
      error("Overflow for the unique IDs.");
    }
    /* Now assign it. */
    s->unique_id.next_batch.current = s->unique_id.global_next_id + local_shift;
    s->unique_id.next_batch.max =
        s->unique_id.global_next_id + local_shift + batch_size;
  }

  /* Shift the position of the next available batch. */
  const long long shift = total_shift * batch_size;
  if (s->unique_id.global_next_id > LLONG_MAX - shift) {
    error("Overflow for the unique IDs.");
  }
  s->unique_id.global_next_id += shift;
}

/**
 * @brief Get a new unique ID.
 *
 * @param s the #space.
 *
 * @return The new unique ID
 */
long long space_get_new_unique_id(struct space *s) {
  /* Do we need unique IDs? */
  if (!star_formation_need_unique_id && !sink_need_unique_id) {
    error("The scheme selected does not seem to use unique ID.");
  }

  /* Get the lock. */
  lock_lock(&s->unique_id.lock);

  /* Get the current available id. */
  const long long id = s->unique_id.current_batch.current;

  /* Update the counter. */
  s->unique_id.current_batch.current++;

  /* Check if everything is fine */
  if (s->unique_id.current_batch.current > s->unique_id.current_batch.max) {
    error("Failed to get a new ID");
  }

  /* Check if need to move to the next batch. */
  else if (s->unique_id.current_batch.current ==
           s->unique_id.current_batch.max) {

    /* Check if the next batch is already used */
    if (s->unique_id.next_batch.current == 0) {
      error("Failed to obtain a new unique ID.");
    }

    s->unique_id.current_batch = s->unique_id.next_batch;

    /* Reset the next batch. */
    s->unique_id.next_batch.current = 0;
    s->unique_id.next_batch.max = 0;
  }

  /* Release the lock. */
  if (lock_unlock(&s->unique_id.lock) != 0) {
    error("Impossible to unlock the unique id.");
  }

  return id;
}

/**
 * @brief Initialize the computation of unique IDs.
 *
 * @param s The #space.
 * @param nr_nodes The number of MPI ranks.
 */
void space_init_unique_id(struct space *s, int nr_nodes) {
  /* Do we need unique IDs? */
  if (!star_formation_need_unique_id && !sink_need_unique_id) {
    return;
  }

  /* Set the counter to 0. */
  s->unique_id.global_next_id = 0;

  /* Check the parts for the max id. */
  for (size_t i = 0; i < s->nr_parts; i++) {
    s->unique_id.global_next_id =
        max(s->unique_id.global_next_id, s->parts[i].id);
  }

  /* Check the gparts for the max id. */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    s->unique_id.global_next_id =
        max(s->unique_id.global_next_id, s->gparts[i].id_or_neg_offset);
  }

  /* Check the sparts for the max id. */
  for (size_t i = 0; i < s->nr_sparts; i++) {
    s->unique_id.global_next_id =
        max(s->unique_id.global_next_id, s->sparts[i].id);
  }

  /* Check the bparts for the max id. */
  for (size_t i = 0; i < s->nr_bparts; i++) {
    s->unique_id.global_next_id =
        max(s->unique_id.global_next_id, s->bparts[i].id);
  }

#ifdef WITH_MPI
  /* Find the global max. */
  MPI_Allreduce(MPI_IN_PLACE, &s->unique_id.global_next_id, 1, MPI_LONG_LONG,
                MPI_MAX, MPI_COMM_WORLD);
#endif  // WITH_MPI

  /* Get the first unique id. */
  if (s->unique_id.global_next_id == LLONG_MAX) {
    error("Overflow for the unique id.");
  }
  s->unique_id.global_next_id++;

  /* Compute the size of each batch. */
  const long long batch_size = (space_extra_parts + space_extra_sparts +
                                space_extra_gparts + space_extra_bparts) *
                               s->nr_cells;

  /* Compute the initial batchs (each rank has 2 batchs). */
  if (s->unique_id.global_next_id > LLONG_MAX - 2 * engine_rank * batch_size) {
    error("Overflow for the unique id.");
  }
  const long long init =
      s->unique_id.global_next_id + 2 * engine_rank * batch_size;

  /* Set the batchs and check for overflows. */
  s->unique_id.current_batch.current = init;

  if (init > LLONG_MAX - batch_size) {
    error("Overflow for the unique id.");
  }
  s->unique_id.current_batch.max = init + batch_size;
  s->unique_id.next_batch.current = s->unique_id.current_batch.max;

  if (s->unique_id.next_batch.current > LLONG_MAX - batch_size) {
    error("Overflow for the unique id.");
  }
  s->unique_id.next_batch.max = s->unique_id.next_batch.current + batch_size;

  /* Update the next available id */
  if (s->unique_id.global_next_id > LLONG_MAX - 2 * batch_size * nr_nodes) {
    error("Overflow for the unique id.");
  }
  s->unique_id.global_next_id += 2 * batch_size * nr_nodes;

  /* Initialize the lock. */
  if (lock_init(&s->unique_id.lock) != 0)
    error("Failed to init spinlock for the unique ids.");
}
