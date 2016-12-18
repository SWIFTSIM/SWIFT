/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 James S. Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_PROFILER_H
#define SWIFT_PROFILER_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "engine.h"

/* Profiler that holds file pointers and time taken in functions. */
struct profiler {

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

  /* Time taken in functions. */
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

/* Function prototypes. */
void profiler_reset_timers(struct profiler *profiler);
void profiler_write_all_timing_info_headers(const struct engine *e,
                                            struct profiler *profiler);
void profiler_write_all_timing_info(const struct engine *e,
                                    struct profiler *profiler);
void profiler_close_files(struct profiler *profiler);

#endif /* SWIFT_PROFILER_H */
