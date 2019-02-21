/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
    "kick1",
    "kick2",
    "timestep",
    "endforce",
    "dosort",
    "doself_density",
    "doself_gradient",
    "doself_force",
    "doself_grav_pp",
    "dopair_density",
    "dopair_gradient",
    "dopair_force",
    "dopair_grav_mm",
    "dopair_grav_pp",
    "dograv_external",
    "dograv_down",
    "dograv_mesh",
    "dograv_top_level",
    "dograv_long_range",
    "dosub_self_density",
    "dosub_self_gradient",
    "dosub_self_force",
    "dosub_self_grav",
    "dosub_pair_density",
    "dosub_pair_gradient",
    "dosub_pair_force",
    "dosub_pair_grav",
    "doself_subset",
    "dopair_subset",
    "dopair_subset_naive",
    "dosub_subset",
    "do_ghost",
    "do_extra_ghost",
    "dorecv_part",
    "dorecv_gpart",
    "dorecv_spart",
    "do_cooling",
    "do_star_formation",
    "gettask",
    "qget",
    "qsteal",
    "locktree",
    "runners",
    "step",
    "do_stars_ghost",
    "logger",
    "do_stars_sort",
    "fof_self",
    "fof_pair",
};

/* File to store the timers */
static FILE* timers_file;

/**
 * @brief Re-set the timers.
 *
 * @param mask A bitmask of the timers to re-set.
 *
 * To reset all timers, use the mask #timers_mask_all.
 */
void timers_reset(unsigned long long mask) {

  /* Loop over the timers and set the masked ones to zero. */
  for (int k = 0; k < timer_count; k++)
    if (mask & (1ull << k)) timers[k] = 0;
}

/**
 * @brief Re-set all the timers.
 *
 */
void timers_reset_all(void) { timers_reset(timers_mask_all); }

/**
 * @brief Outputs all the timers to the timers dump file.
 *
 * @param step The current step.
 */
void timers_print(int step) {
  fprintf(timers_file, "%d\t", step);
  for (int k = 0; k < timer_count; k++)
    fprintf(timers_file, "%18.3f ", clocks_from_ticks(timers[k]));
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

  fprintf(timers_file, "# timers: \n# step | ");
  for (int k = 0; k < timer_count; k++)
    fprintf(timers_file, "%18s ", timers_names[k]);
  fprintf(timers_file, "\n");
}

/**
 * @brief Close the file containing the timer info.
 */
void timers_close_file(void) { fclose(timers_file); }
