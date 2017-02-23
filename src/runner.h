/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_H
#define SWIFT_RUNNER_H

#include "cache.h"

extern const double runner_shift[13][3];
extern const char runner_flip[27];

struct cell;
struct engine;

/**
 * @brief A struct representing a runner's thread and its data.
 */
struct runner {

  /*! The id of this thread. */
  int id;

  /*! The actual thread which it is running. */
  pthread_t thread;

  /*! The queue to use to get tasks. */
  int cpuid, qid;

  /*! The engine owing this runner. */
  struct engine *e;

  /*! The particle cache of this runner. */
  struct cache par_cache;
};

/* Function prototypes. */
void runner_do_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_extra_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_sort(struct runner *r, struct cell *c, int flag, int clock);
void runner_do_drift_particles(struct runner *r, struct cell *c, int timer);
void runner_do_kick1(struct runner *r, struct cell *c, int timer);
void runner_do_kick2(struct runner *r, struct cell *c, int timer);
void runner_do_end_force(struct runner *r, struct cell *c, int timer);
void runner_do_init(struct runner *r, struct cell *c, int timer);
void runner_do_cooling(struct runner *r, struct cell *c, int timer);
void runner_do_grav_external(struct runner *r, struct cell *c, int timer);
void *runner_main(void *data);
void runner_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data);
void runner_do_drift_mapper(void *map_data, int num_elements, void *extra_data);

#endif /* SWIFT_RUNNER_H */
