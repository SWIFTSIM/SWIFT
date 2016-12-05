/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
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
#ifndef SWIFT_TIMER_H
#define SWIFT_TIMER_H

/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Some standard headers. */
#include <pthread.h>
#include <stdio.h>

/* Includes. */
#include "clocks.h"
#include "cooling_struct.h"
#include "parser.h"
#include "partition.h"
#include "potential.h"
#include "runner.h"
#include "scheduler.h"
#include "sourceterms_struct.h"
#include "space.h"
#include "task.h"
#include "units.h"

struct timer {

  /* File pointers for timing info. */
  FILE *file_engine_collect_timesteps;
  FILE *file_engine_drift;
  FILE *file_engine_rebuild;
  FILE *file_scheduler_reweight;
  FILE *file_scheduler_clear_waits;
  FILE *file_scheduler_re_wait;
  FILE *file_scheduler_enqueue;
  FILE *file_engine_stats;
  FILE *file_engine_launch;
  FILE *file_space_rebuild;
  FILE *file_engine_maketasks;
  FILE *file_engine_marktasks;
  FILE *file_space_regrid;
  FILE *file_space_parts_sort;
  FILE *file_space_split;
  FILE *file_space_parts_get_cell_id;
  FILE *file_space_count_parts;

  ticks collect_timesteps_time;
  ticks drift_time;
  ticks rebuild_time;
  ticks reweight_time;   
  ticks clear_waits_time;   
  ticks re_wait_time;   
  ticks enqueue_time;   
  ticks stats_time;
  ticks launch_time;
  ticks space_rebuild_time;
  ticks engine_maketasks_time;
  ticks engine_marktasks_time;
  ticks space_regrid_time;
  ticks space_parts_sort_time;
  ticks space_split_time;
  ticks space_parts_get_cell_id_time;
  ticks space_count_parts_time;
};

void timer_reset_timers(struct timer *profiler);
void timer_write_timing_info_header(struct engine *e, char *fileName, char *functionName, FILE **file);
void timer_write_all_timing_info_headers(struct engine *e, struct timer *profiler);
void timer_write_timing_info(struct engine *e, ticks time, FILE **file);
void timer_write_all_timing_info(struct engine *e, struct timer *profiler);

#endif /* SWIFT_TIMER_H */
