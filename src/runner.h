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

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "cache.h"
#include "gravity_cache.h"

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

  /*! The particle gravity_cache of cell ci. */
  struct gravity_cache ci_gravity_cache;

  /*! The particle gravity_cache of cell cj. */
  struct gravity_cache cj_gravity_cache;

#ifdef WITH_VECTORIZATION

  /*! The particle cache of cell ci. */
  struct cache ci_cache;

  /*! The particle cache of cell cj. */
  struct cache cj_cache;
#endif
};

/* Function prototypes. */
void runner_do_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_extra_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_hydro_sort(struct runner *r, struct cell *c, int flag,
                          int cleanup, int clock);
void runner_do_stars_sort(struct runner *r, struct cell *c, int flag,
                          int cleanup, int clock);
void runner_do_drift_part(struct runner *r, struct cell *c, int timer);
void runner_do_drift_gpart(struct runner *r, struct cell *c, int timer);
void runner_do_drift_spart(struct runner *r, struct cell *c, int timer);
void runner_do_kick1(struct runner *r, struct cell *c, int timer);
void runner_do_kick2(struct runner *r, struct cell *c, int timer);
void runner_do_end_hydro_force(struct runner *r, struct cell *c, int timer);
void runner_do_init(struct runner *r, struct cell *c, int timer);
void runner_do_cooling(struct runner *r, struct cell *c, int timer);
void runner_do_grav_external(struct runner *r, struct cell *c, int timer);
void runner_do_grav_fft(struct runner *r, int timer);
void runner_do_logger(struct runner *r, struct cell *c, int timer);
void runner_do_fof_self(struct runner *r, struct cell *c, int timer);
void runner_do_fof_pair(struct runner *r, struct cell *ci, struct cell *cj,
                        int timer);
void *runner_main(void *data);
void runner_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data);
void runner_do_drift_all_mapper(void *map_data, int num_elements,
                                void *extra_data);

#endif /* SWIFT_RUNNER_H */
