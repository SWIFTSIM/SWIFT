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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>
#include <string.h>

/* This object's header. */
#include "profiler.h"

/* Local includes */
#include "clocks.h"
#include "hydro.h"
#include "version.h"

/* Array to store the list of file names. Order must match profiler_types
 * enumerator and profiler_func_names. */
const char *profiler_file_names[profiler_length] = {"enginecollecttimesteps",
                                                    "enginedrift",
                                                    "enginerebuild",
                                                    "schedulerreweight",
                                                    "schedulerclearwaits",
                                                    "schedulerrewait",
                                                    "schedulerenqueue",
                                                    "engineprintstats",
                                                    "enginelaunch",
                                                    "spacerebuild",
                                                    "enginemaketasks",
                                                    "enginemarktasks",
                                                    "spaceregrid",
                                                    "spacepartssort",
                                                    "spacesplit",
                                                    "spacegetcellid",
                                                    "spacecountparts"};

/* Array to store the list of function names. Order must match profiler_types
 * enumerator and profiler_file_names. */
const char *profiler_func_names[profiler_length] = {"engine_collect_timesteps",
                                                    "engine_drift",
                                                    "engine_rebuild",
                                                    "scheduler_reweight",
                                                    "scheduler_clear_waits",
                                                    "scheduler_rewait",
                                                    "scheduler_enqueue",
                                                    "engine_print_stats",
                                                    "engine_launch",
                                                    "space_rebuild",
                                                    "engine_maketasks",
                                                    "engine_marktasks",
                                                    "space_regrid",
                                                    "space_parts_sort",
                                                    "space_split",
                                                    "space_get_cell_id",
                                                    "space_count_parts"};

/**
 * @brief Resets all timers.
 *
 * @param profiler #profiler object that holds file pointers and
 * function timers.
 */
void profiler_reset_timers(struct profiler *profiler) {
  /* Iterate over times array and reset values. */
  for (int i = 0; i < profiler_length; i++) profiler->times[i] = 0;
}

/**
 * @brief Opens an output file and populates the header.
 *
 * @param e #engine object to get various properties.
 * @param fileName name of file to be written to.
 * @param functionName name of function that is being timed.
 * @param file (return) pointer used to open output file.
 */
void profiler_write_timing_info_header(const struct engine *e,
                                       const char *fileName,
                                       const char *functionName, FILE **file) {

  /* Create the file name in the format: "fileName_(no. of threads)" */
  char fullFileName[200] = "";
  sprintf(fullFileName + strlen(fullFileName), "%s_%d.txt", fileName,
          e->nr_nodes * e->nr_threads);

  /* Open the file and write the header. */
  *file = fopen(fullFileName, "w");
  fprintf(*file,
          "# Host: %s\n# Branch: %s\n# Revision: %s\n# Compiler: %s, "
          "Version: %s \n# "
          "Number of threads: %d\n# Number of MPI ranks: %d\n# Hydrodynamic "
          "scheme: %s\n# Hydrodynamic kernel: %s\n# No. of neighbours: %.2f "
          "+/- %.2f\n# Eta: %f\n"
          "# %6s %14s %14s %10s %10s %10s %16s [%s]\n",
          hostname(), functionName, git_revision(), compiler_name(),
          compiler_version(), e->nr_threads, e->nr_nodes, SPH_IMPLEMENTATION,
          kernel_name, e->hydro_properties->target_neighbours,
          e->hydro_properties->delta_neighbours,
          e->hydro_properties->eta_neighbours, "Step", "Time", "Time-step",
          "Updates", "g-Updates", "s-Updates", "Wall-clock time",
          clocks_getunit());

  fflush(*file);
}

/**
 * @brief Writes the headers for all output files. Should be called once at the
 * start of the simulation, it could be called in engine_init() for example.
 *
 * @param e #engine object to get various properties.
 * @param profiler #profiler object that holds file pointers and
 * function timers.
 */
void profiler_write_all_timing_info_headers(const struct engine *e,
                                            struct profiler *profiler) {
  /* Iterate over files array and write file headers. */
  for (int i = 0; i < profiler_length; i++) {
    profiler_write_timing_info_header(
        e, profiler_file_names[i], profiler_func_names[i], &profiler->files[i]);
  }
}

/**
 * @brief Writes timing info to the output file.
 *
 * @param e #engine object to get various properties.
 * @param time Time in ticks to be written to the output file.
 * @param file pointer used to open output file.
 */
void profiler_write_timing_info(const struct engine *e, ticks time,
                                FILE *file) {

  fprintf(file, "  %6d %14e %14e %10lld %10lld %10lld %21.3f\n", e->step,
          e->time, e->time_step, e->updates, e->g_updates, e->s_updates,
          clocks_from_ticks(time));
  fflush(file);
}

/**
 * @brief Writes timing info to all output files. Should be called at the end of
 * every time step, in engine_step() for example.
 *
 * @param e #engine object to get various properties.
 * @param profiler #profiler object that holds file pointers and
 * function timers.
 */
void profiler_write_all_timing_info(const struct engine *e,
                                    struct profiler *profiler) {

  /* Iterate over times array and print timing info to files. */
  for (int i = 0; i < profiler_length; i++) {
    profiler_write_timing_info(e, profiler->times[i], profiler->files[i]);
  }
  /* Reset timers. */
  profiler_reset_timers(profiler);
}

/**
 * @brief Closes all output files, should be called at the end of the
 * simulation.
 *
 * @param profiler #profiler object that holds file pointers and
 * function timers.
 */
void profiler_close_files(struct profiler *profiler) {

  /* Iterate over files array and close files. */
  for (int i = 0; i < profiler_length; i++) fclose(profiler->files[i]);
}
