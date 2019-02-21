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
#ifndef SWIFT_TIMERS_H
#define SWIFT_TIMERS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "atomic.h"
#include "cycle.h"
#include "inline.h"

/**
 * @brief The timers themselves.
 *
 * If you modify this list, be sure to change timers_names in timers.c as
 * well!
 **/
enum {
  timer_none = 0,
  timer_prepare,
  timer_init,
  timer_init_grav,
  timer_drift_part,
  timer_drift_gpart,
  timer_kick1,
  timer_kick2,
  timer_timestep,
  timer_endforce,
  timer_dosort,
  timer_doself_density,
  timer_doself_gradient,
  timer_doself_force,
  timer_doself_grav_pp,
  timer_dopair_density,
  timer_dopair_gradient,
  timer_dopair_force,
  timer_dopair_grav_mm,
  timer_dopair_grav_pp,
  timer_dograv_external,
  timer_dograv_down,
  timer_dograv_mesh,
  timer_dograv_top_level,
  timer_dograv_long_range,
  timer_dosub_self_density,
  timer_dosub_self_gradient,
  timer_dosub_self_force,
  timer_dosub_self_grav,
  timer_dosub_pair_density,
  timer_dosub_pair_gradient,
  timer_dosub_pair_force,
  timer_dosub_pair_grav,
  timer_doself_subset,
  timer_dopair_subset,
  timer_dopair_subset_naive,
  timer_dosub_subset,
  timer_do_ghost,
  timer_do_extra_ghost,
  timer_dorecv_part,
  timer_dorecv_gpart,
  timer_dorecv_spart,
  timer_do_cooling,
  timer_do_star_formation,
  timer_gettask,
  timer_qget,
  timer_qsteal,
  timer_locktree,
  timer_runners,
  timer_step,
  timer_dostars_ghost,
  timer_logger,
  timer_do_stars_sort,
  timer_fof_self,
  timer_fof_pair,
  timer_count,
};

/* The timers. */
extern ticks timers[timer_count];

/* The timer names. */
extern const char *timers_names[];

/* Mask for all timers. */
#define timers_mask_all ((1ull << timer_count) - 1)

/* Define the timer macros. */
#ifdef SWIFT_USE_TIMERS
#define TIMER_TIC const ticks tic = getticks();
#define TIMER_TOC(t) timers_toc(t, tic)
#define TIMER_TIC2 const ticks tic2 = getticks();
#define TIMER_TOC2(t) timers_toc(t, tic2)
INLINE static ticks timers_toc(unsigned int t, ticks tic) {
  const ticks d = (getticks() - tic);
  atomic_add(&timers[t], d);
  return d;
}
#else
#define TIMER_TIC
#define TIMER_TOC(t) (void)0
#define TIMER_TIC2
#define TIMER_TOC2(t) (void)0
#endif

/* Function prototypes. */
void timers_reset_all(void);
void timers_reset(unsigned long long mask);
void timers_open_file(int rank);
void timers_close_file(void);
void timers_print(int step);

#endif /* SWIFT_TIMERS_H */
