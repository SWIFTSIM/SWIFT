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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <sched.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef HAVE_LIBNUMA
#include <numa.h>
#endif

/* This object's header. */
#include "profiler.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "hydro.h"
#include "minmax.h"
#include "parallel_io.h"
#include "part.h"
#include "partition.h"
#include "proxy.h"
#include "runner.h"
#include "serial_io.h"
#include "single_io.h"
#include "statistics.h"
#include "timers.h"
#include "tools.h"
#include "units.h"
#include "version.h"

void profiler_reset_timers(struct profiler *profiler) {

  profiler->collect_timesteps_time = 0;
  profiler->drift_time = 0;
  profiler->rebuild_time = 0;
  profiler->reweight_time = 0;   
  profiler->clear_waits_time = 0;   
  profiler->re_wait_time = 0;   
  profiler->enqueue_time = 0;   
  profiler->stats_time = 0;
  profiler->launch_time = 0;
  profiler->space_rebuild_time = 0;
  profiler->engine_maketasks_time = 0;
  profiler->engine_marktasks_time = 0;
  profiler->space_regrid_time = 0;
  profiler->space_parts_sort_time = 0;
  profiler->space_split_time = 0;
  profiler->space_parts_get_cell_id_time = 0;
  profiler->space_count_parts_time = 0;
}

void profiler_write_timing_info_header(struct engine *e, char *fileName, char *functionName, FILE **file) {

  char fullFileName[200] = "";
  sprintf(fullFileName + strlen(fullFileName), "%s_%d.txt", fileName, e->nr_nodes * e->nr_threads);

  *file = fopen(fullFileName, "w");
  fprintf(*file,
      "# Host: %s\n# Branch: %s\n# Revision: %s\n# Compiler: %s, "
      "Version: %s \n# "
      "Number of threads: %d\n# Number of MPI ranks: %d\n# Hydrodynamic "
      "scheme: %s\n# Hydrodynamic kernel: %s\n# No. of neighbours: %.2f "
      "+/- %.2f\n# Eta: %f\n"
      "# %6s %14s %14s %10s %10s %16s [%s]\n", 
      hostname(), functionName, git_revision(), compiler_name(),
      compiler_version(), e->nr_threads, e->nr_nodes, SPH_IMPLEMENTATION,
      kernel_name, e->hydro_properties->target_neighbours,
      e->hydro_properties->delta_neighbours,
      e->hydro_properties->eta_neighbours,
      "Step","Time", "Time-step", "Updates", "g-Updates", "Wall-clock time",
      clocks_getunit());

  fflush(*file);
}

void profiler_write_all_timing_info_headers(struct engine *e, struct profiler *profiler) {

  profiler_write_timing_info_header(e,"enginecollecttimesteps","engine_collect_timesteps",&profiler->file_engine_collect_timesteps);
  profiler_write_timing_info_header(e,"enginedrift","engine_drift",&profiler->file_engine_drift);
  profiler_write_timing_info_header(e,"enginerebuild","engine_rebuild",&profiler->file_engine_rebuild);
  profiler_write_timing_info_header(e,"schedulerreweight","scheduler_reweight",&profiler->file_scheduler_reweight);
  profiler_write_timing_info_header(e,"schedulerclearwaits","scheduler_clear_waits",&profiler->file_scheduler_clear_waits);
  profiler_write_timing_info_header(e,"schedulerrewait","scheduler_rewait",&profiler->file_scheduler_re_wait);
  profiler_write_timing_info_header(e,"schedulerenqueue","scheduler_enqueue",&profiler->file_scheduler_enqueue);
  profiler_write_timing_info_header(e,"engineprintstats","engine_print_stats",&profiler->file_engine_stats);
  profiler_write_timing_info_header(e,"enginelaunch","engine_launch",&profiler->file_engine_launch);
  profiler_write_timing_info_header(e,"spacerebuild","space_rebuild",&profiler->file_space_rebuild);
  profiler_write_timing_info_header(e,"enginemaketasks","engine_maketasks",&profiler->file_engine_maketasks);
  profiler_write_timing_info_header(e,"enginemarktasks","engine_marktasks",&profiler->file_engine_marktasks);
  profiler_write_timing_info_header(e,"spaceregrid","space_regrid",&profiler->file_space_regrid);
  profiler_write_timing_info_header(e,"spacepartssort","space_parts_sort",&profiler->file_space_parts_sort);
  profiler_write_timing_info_header(e,"spacesplit","space_split",&profiler->file_space_split);
  profiler_write_timing_info_header(e,"spacegetcellid","space_get_cell_id",&profiler->file_space_parts_get_cell_id);
  profiler_write_timing_info_header(e,"spacecountparts","space_count_parts",&profiler->file_space_count_parts);
}

void profiler_write_timing_info(struct engine *e, ticks time, FILE **file) {

  fprintf(*file, "  %6d %14e %14e %10zu %10zu %21.3f\n", e->step,
      e->time, e->timeStep, e->updates, e->g_updates, clocks_from_ticks(time));
  fflush(*file);
}

void profiler_write_all_timing_info(struct engine *e, struct profiler *profiler) {

  profiler_write_timing_info(e,profiler->drift_time,&profiler->file_engine_drift);
  profiler_write_timing_info(e,profiler->rebuild_time,&profiler->file_engine_rebuild);
  profiler_write_timing_info(e,profiler->reweight_time,&profiler->file_scheduler_reweight);
  profiler_write_timing_info(e,profiler->clear_waits_time,&profiler->file_scheduler_clear_waits);
  profiler_write_timing_info(e,profiler->re_wait_time,&profiler->file_scheduler_re_wait);
  profiler_write_timing_info(e,profiler->enqueue_time,&profiler->file_scheduler_enqueue);
  profiler_write_timing_info(e,profiler->stats_time,&profiler->file_engine_stats);
  profiler_write_timing_info(e,profiler->launch_time,&profiler->file_engine_launch);
  profiler_write_timing_info(e,profiler->space_rebuild_time,&profiler->file_space_rebuild);
  profiler_write_timing_info(e,profiler->engine_maketasks_time,&profiler->file_engine_maketasks);
  profiler_write_timing_info(e,profiler->engine_marktasks_time,&profiler->file_engine_marktasks);
  profiler_write_timing_info(e,profiler->space_regrid_time,&profiler->file_space_regrid);
  profiler_write_timing_info(e,profiler->space_parts_sort_time,&profiler->file_space_parts_sort);
  profiler_write_timing_info(e,profiler->space_split_time,&profiler->file_space_split);
  profiler_write_timing_info(e,profiler->space_parts_get_cell_id_time,&profiler->file_space_parts_get_cell_id);
  profiler_write_timing_info(e,profiler->space_count_parts_time,&profiler->file_space_count_parts);

  profiler_reset_timers(profiler); 
}
