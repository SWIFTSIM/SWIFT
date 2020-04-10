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

#ifndef SWIFT_SPACE_UNIQUE_ID_H
#define SWIFT_SPACE_UNIQUE_ID_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "lock.h"

/* Predefine the space structure */
struct space;

/**
 * @brief Batch of unique IDs for particle creation.
 */
struct batch {
  /*! Current free unique id */
  long long current;

  /*! Maximal unique id in this batch  (not included) */
  long long max;
};

/*! Structure dealing with the computation of a unique ID */
struct unique_id {
  /*! Current batch of unique ids */
  struct batch current_batch;

  /*! Next batch of unique ids */
  struct batch next_batch;

  /* Global next slot available */
  long long global_next_id;

  /* Lock for the unique ids */
  swift_lock_type lock;
};

void space_update_unique_id(struct space *s);
long long space_get_new_unique_id(struct space *s);
void space_init_unique_id(struct space *s, int nr_nodes);

#endif  // SWIFT_SPACE_UNIQUE_ID_H
