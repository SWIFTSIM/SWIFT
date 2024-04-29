/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* This object's header. */
#include "timers.h"

/* Some standard headers. */
#include <stdio.h>

/* Local includes. */
#include "clocks.h"
#include "error.h"

/* The timers. */
ticks timers[timer_count];

/* Timer names. */
const char* timers_names[timer_count] = {
    "none",
    "prepare",
    "init",
    "init_grav",
    "drift_part",
    "drift_gpart",
    "drift_spart",
    "drift_bpart",
    "kick1",
    "kick2",
    "timestep",
    "end_hydro_force",
    "end_grav_force",
    "dosort",
    "doself_density",
    "doself_gradient",
    "doself_force",
    "doself_limiter",
    "doself_stars_density",
    "doself_stars_feedback",
    "doself_bh_density",
    "doself_bh_swallow",
    "doself_bh_feedback",
    "doself_grav_pp",
    "doself_sink_swallow",
    "dopair_density",
    "dopair_gradient",
    "dopair_force",
    "dopair_limiter",
    "dopair_stars_density",
    "dopair_stars_feedback",
    "dopair_bh_density",
    "dopair_bh_swallow",
    "dopair_bh_feedback",
    "dopair_grav_mm",
    "dopair_grav_pp",
    "dopair_sink_swallow",
    "dograv_external",
    "dograv_down",
    "dograv_mesh",
    "dograv_top_level",
    "dograv_long_range",
    "dosub_self_density",
    "dosub_self_gradient",
    "dosub_self_force",
    "dosub_self_limiter",
    "dosub_self_stars_density",
    "dosub_self_stars_feedback",
    "dosub_self_bh_density",
    "dosub_self_bh_swallow",
    "dosub_self_bh_feedback",
    "dosub_self_grav",
    "dosub_self_sink_swallow",
    "dosub_pair_density",
    "dosub_pair_gradient",
    "dosub_pair_force",
    "dosub_pair_limiter",
    "dosub_pair_stars_density",
    "dosub_pair_stars_feedback",
    "dosub_pair_bh_density",
    "dosub_pair_bh_swallow",
    "dosub_pair_bh_feedback",
    "dosub_pair_grav",
    "dosub_pair_sink_swallow",
    "doself_subset",
    "dopair_subset",
    "dopair_subset_naive",
    "dosub_subset",
    "do_ghost",
    "do_extra_ghost",
    "do_stars_ghost",
    "do_black_holes_ghost",
    "dorecv_part",
    "dorecv_gpart",
    "dorecv_spart",
    "dorecv_bpart",
    "do_limiter",
    "do_cooling",
    "do_star_formation",
    "do_star_evol",
    "gettask",
    "qget",
    "qsteal",
    "locktree",
    "runners",
    "step",
    "csds",
    "do_stars_sort",
    "do_stars_resort",
    "fof_self",
    "fof_pair",
    "drift_sink",
    "rt_ghost1",
    "rt_ghost2",
    "doself_rt_gradient",
    "dopair_rt_gradient",
    "dosub_self_rt_gradient",
    "dosub_pair_rt_gradient",
    "doself_rt_transport",
    "dopair_rt_transport",
    "dosub_self_rt_transport",
    "dosub_pair_rt_transport",
    "rt_tchem",
    "rt_advance_cell_time",
    "rt_collect_times",
    "do_sync",
    "neutrino_weighting",
};

/* File to store the timers */
static FILE* timers_file;

/**
 * @brief Re-set all the timers.
 *
 */
void timers_reset_all(void) {

  for (int k = 0; k < timer_count; k++) timers[k] = 0;
}

/**
 * @brief Outputs all the timers to the timers dump file.
 *
 * @param step The current step.
 */
void timers_print(int step) {
  fprintf(timers_file, "%d\t", step);
  for (int k = 0; k < timer_count; k++)
    fprintf(timers_file, "%25.3f ", clocks_from_ticks(timers[k]));
  fprintf(timers_file, "\n");
  fflush(timers_file);
}

/**
 * @brief Opens the file to contain the timers info and print a header
 *
 * @param rank The MPI rank of the file.
 */
void timers_open_file(int rank) {

  char buff[100];
  sprintf(buff, "timers_%d.txt", rank);
  timers_file = fopen(buff, "w");
  if (timers_file == NULL) error("Could not create file '%s'.", buff);

  fprintf(timers_file, "# timers: \n# step |");
  for (int k = 0; k < timer_count; k++)
    fprintf(timers_file, "%25s ", timers_names[k]);
  fprintf(timers_file, "\n");
}

/**
 * @brief Close the file containing the timer info.
 */
void timers_close_file(void) { fclose(timers_file); }
