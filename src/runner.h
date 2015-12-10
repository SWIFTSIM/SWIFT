/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Some standard headers. */
#include <pthread.h>

/* Includes. */
#include "cell.h"
#include "inline.h"

/* Forward-declare the engine type to avoid cyclic header dependencies. */
struct engine;

/* Some constants/flags. */
#define runner_prefetch 0

/* SID stuff. */
extern const char runner_flip[];

/* Counters. */
enum runner_counters {
  runner_counter_swap = 0,
  runner_counter_stall,
  runner_counter_steal_stall,
  runner_counter_steal_empty,
  runner_counter_keep,
  runner_counter_iact,
  runner_counter_count,
};
extern int runner_counter[runner_counter_count];

/* Counter macros. */
#ifdef COUNTER
#define COUNT(c) (__sync_add_and_fetch(&runner_counter[c], 1))
#else
#define COUNT(c)
#endif

/* Histogram functions. */
#define runner_hist_a 1.0
#define runner_hist_b 100.0
#define runner_hist_N 99
extern long long int runner_hist_bins[runner_hist_N];
#define runner_hist_hit(x)                                                   \
  __sync_add_and_fetch(                                                      \
      &runner_hist_bins[(int)fmax(                                           \
          0.0, fmin(runner_hist_N - 1, ((x) - runner_hist_a) /               \
                                           (runner_hist_b - runner_hist_a) * \
                                           runner_hist_N))],                 \
      1)

/* A struct representing a runner's thread and its data. */
struct runner {

  /* The id of this thread. */
  int id;

  /* The thread which it is running. */
  pthread_t thread;

  /* The queue to use to get tasks. */
  int cpuid, qid;

  /* The underlying runner. */
  struct engine *e;
};

/* Function prototypes. */
void runner_doghost(struct runner *r, struct cell *c);
void runner_dopair_density(struct runner *r, struct cell *ci, struct cell *cj);
void runner_doself_density(struct runner *r, struct cell *c);
void runner_dosub_density(struct runner *r, struct cell *ci, struct cell *cj,
                          int flags);
void runner_dosort(struct runner *r, struct cell *c, int flag, int clock);
void runner_dogsort(struct runner *r, struct cell *c, int flag, int clock);
void runner_dokick(struct runner *r, struct cell *c, int timer);
void runner_dodrift(struct runner *r, struct cell *c, int timer);
void runner_doinit(struct runner *r, struct cell *c);
void *runner_main(void *data);

#endif /* SWIFT_RUNNER_H */
