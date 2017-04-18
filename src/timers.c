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
#include "../config.h"

/* This object's header. */
#include "timers.h"

/* The timers. */
ticks timers[timer_count];

/* Timer names. */
char *timers_names[timer_count] = {
    "none",
    "prepare",
    "init",
    "init_grav",
    "drift",
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
    "dopair_grav_pm",
    "dopair_grav_mm",
    "dopair_grav_pp",
    "dograv_external",
    "dograv_down",
    "dograv_long_range",
    "dosource",
    "dosub_self_density",
    "dosub_self_gradient",
    "dosub_self_force",
    "dosub_self_grav",
    "dosub_pair_density",
    "dosub_pair_gradient",
    "dosub_pair_force",
    "dosub_pair_grav",
    "dopair_subset",
    "do_ghost",
    "do_extra_ghost",
    "dorecv_part",
    "dorecv_gpart",
    "dorecv_spart",
    "gettask",
    "qget",
    "qsteal",
    "runners",
    "step",
    "do_cooling",
};

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
