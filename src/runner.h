/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* Local headers. */
#include "cache.h"
#include "gravity_cache.h"

struct cell;
struct engine;
struct task;

/* Unique identifier of loop types */
#define TASK_LOOP_DENSITY 0
#define TASK_LOOP_GRADIENT 1
#define TASK_LOOP_FORCE 2
#define TASK_LOOP_LIMITER 3
#define TASK_LOOP_FEEDBACK 4
#define TASK_LOOP_SWALLOW 5
#define TASK_LOOP_SINK_SWALLOW 6
#define TASK_LOOP_SINK_DO_SINK_SWALLOW 7
#define TASK_LOOP_SINK_DO_GAS_SWALLOW 8
#define TASK_LOOP_STARS_PREP1 9
#define TASK_LOOP_STARS_PREP2 10
#define TASK_LOOP_RT_GRADIENT 11
#define TASK_LOOP_RT_TRANSPORT 12

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

  /*! Time this runner was active during the last engine_launch. */
  ticks active_time;

#ifdef WITH_VECTORIZATION

  /*! The particle cache of cell ci. */
  struct cache ci_cache;

  /*! The particle cache of cell cj. */
  struct cache cj_cache;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /*! Pointer to the task this runner is currently performing */
  const struct task *t;
#endif
};

/* Function prototypes. */
void runner_do_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_extra_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_stars_ghost(struct runner *r, struct cell *c, int timer);
void runner_do_black_holes_density_ghost(struct runner *r, struct cell *c,
                                         int timer);
void runner_do_black_holes_swallow_ghost(struct runner *r, struct cell *c,
                                         int timer);
void runner_do_init_grav(struct runner *r, struct cell *c, int timer);
void runner_do_hydro_sort(struct runner *r, struct cell *c, int flag,
                          int cleanup, int rt_requests_sort, int clock);
void runner_do_stars_sort(struct runner *r, struct cell *c, int flag,
                          int cleanup, int clock);
void runner_do_all_hydro_sort(struct runner *r, struct cell *c);
void runner_do_all_stars_sort(struct runner *r, struct cell *c);
void runner_do_drift_part(struct runner *r, struct cell *c, int timer);
void runner_do_drift_gpart(struct runner *r, struct cell *c, int timer);
void runner_do_drift_spart(struct runner *r, struct cell *c, int timer);
void runner_do_drift_sink(struct runner *r, struct cell *c, int timer);
void runner_do_drift_bpart(struct runner *r, struct cell *c, int timer);
void runner_do_kick1(struct runner *r, struct cell *c, int timer);
void runner_do_kick2(struct runner *r, struct cell *c, int timer);
void runner_do_timestep(struct runner *r, struct cell *c, int timer);
void runner_do_timestep_collect(struct runner *r, struct cell *c, int timer);
void runner_do_end_hydro_force(struct runner *r, struct cell *c, int timer);
void runner_do_end_grav_force(struct runner *r, struct cell *c, int timer);
void runner_do_init(struct runner *r, struct cell *c, int timer);
void runner_do_cooling(struct runner *r, struct cell *c, int timer);
void runner_do_limiter(struct runner *r, struct cell *c, int force, int timer);
void runner_do_sync(struct runner *r, struct cell *c, int force, int timer);
void runner_do_grav_mesh(struct runner *r, struct cell *c, int timer);
void runner_do_grav_external(struct runner *r, struct cell *c, int timer);
void runner_do_grav_fft(struct runner *r, int timer);
void runner_do_csds(struct runner *r, struct cell *c, int timer);
void runner_do_fof_search_self(struct runner *r, struct cell *c, int timer);
void runner_do_fof_search_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer);
void runner_do_fof_attach_self(struct runner *r, struct cell *c, int timer);
void runner_do_fof_attach_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer);
void runner_do_rt_ghost1(struct runner *r, struct cell *c, int timer);
void runner_do_rt_ghost2(struct runner *r, struct cell *c, int timer);
void runner_do_rt_tchem(struct runner *r, struct cell *c, int timer);
void runner_do_gas_swallow_self(struct runner *r, struct cell *c, int timer);
void runner_do_bh_swallow_self(struct runner *r, struct cell *c, int timer);
void runner_do_gas_swallow_pair(struct runner *r, struct cell *ci,
                                struct cell *cj, int timer);
void runner_do_bh_swallow_pair(struct runner *r, struct cell *ci,
                               struct cell *cj, int timer);
void runner_do_star_formation(struct runner *r, struct cell *c, int timer);
void runner_do_star_formation_sink(struct runner *r, struct cell *c, int timer);
void runner_do_sink_formation(struct runner *r, struct cell *c);
void runner_do_stars_resort(struct runner *r, struct cell *c, const int timer);

void runner_do_recv_gpart(struct runner *r, struct cell *c, int timer);
void runner_do_recv_part(struct runner *r, struct cell *c, int clear_sorts,
                         int timer);
void runner_do_recv_spart(struct runner *r, struct cell *c, int clear_sorts,
                          int timer);
void runner_do_recv_bpart(struct runner *r, struct cell *c, int clear_sorts,
                          int timer);
void runner_do_pack_limiter(struct runner *r, struct cell *c, void **buffer,
                            const int timer);
void runner_do_unpack_limiter(struct runner *r, struct cell *c, void *buffer,
                              const int timer);
void runner_do_neutrino_weighting(struct runner *r, struct cell *c, int timer);
void runner_do_rt_advance_cell_time(struct runner *r, struct cell *c,
                                    int timer);
void runner_do_collect_rt_times(struct runner *r, struct cell *c,
                                const int timer);
void *runner_main(void *data);

ticks runner_get_active_time(const struct runner *restrict r);
void runner_reset_active_time(struct runner *restrict r);

#endif /* SWIFT_RUNNER_H */
