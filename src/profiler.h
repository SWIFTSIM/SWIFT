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

/* Enumerator to be used as an index into the timers and files array. To add an
 * extra timer extend this list, before the profiler_length value.*/
enum profiler_types {
  profiler_engine_collect_timesteps = 0,
  profiler_engine_drift,
  profiler_engine_rebuild,
  profiler_scheduler_reweight,
  profiler_scheduler_clear_waits,
  profiler_scheduler_re_wait,
  profiler_scheduler_enqueue,
  profiler_engine_stats,
  profiler_engine_launch,
  profiler_space_rebuild,
  profiler_engine_maketasks,
  profiler_engine_marktasks,
  profiler_space_regrid,
  profiler_space_parts_sort,
  profiler_space_split,
  profiler_space_parts_get_cell_id,
  profiler_space_count_parts,
  profiler_length
};

/* Profiler that holds file pointers and time taken in functions. */
struct profiler {
  FILE *files[profiler_length];
  ticks times[profiler_length];
};

/* Function prototypes. */
void profiler_reset_timers(struct profiler *profiler);
void profiler_write_all_timing_info_headers(const struct engine *e,
                                            struct profiler *profiler);
void profiler_write_all_timing_info(const struct engine *e,
                                    struct profiler *profiler);
void profiler_close_files(struct profiler *profiler);

#endif /* SWIFT_PROFILER_H */
