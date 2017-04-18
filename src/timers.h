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

/* Includes. */
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
  timer_drift,
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
  timer_dopair_grav_pm,
  timer_dopair_grav_mm,
  timer_dopair_grav_pp,
  timer_dograv_external,
  timer_dograv_down,
  timer_dograv_long_range,
  timer_dosource,
  timer_dosub_self_density,
  timer_dosub_self_gradient,
  timer_dosub_self_force,
  timer_dosub_self_grav,
  timer_dosub_pair_density,
  timer_dosub_pair_gradient,
  timer_dosub_pair_force,
  timer_dosub_pair_grav,
  timer_dopair_subset,
  timer_do_ghost,
  timer_do_extra_ghost,
  timer_dorecv_part,
  timer_dorecv_gpart,
  timer_dorecv_spart,
  timer_gettask,
  timer_qget,
  timer_qsteal,
  timer_runners,
  timer_step,
  timer_do_cooling,
  timer_count,
};

/* The timers. */
extern ticks timers[timer_count];

/* The timer names. */
extern char *timers_names[];

/* Mask for all timers. */
#define timers_mask_all ((1ull << timer_count) - 1)

/* Define the timer macros. */
#ifdef TIMER
#define TIMER_TIC_ND tic = getticks();
#define TIMER_TIC2_ND ticks tic2 = getticks();
#define TIMER_TIC ticks tic = getticks();
#define TIMER_TOC(t) timers_toc(t, tic)
#define TIMER_TIC2 ticks tic2 = getticks();
#define TIMER_TOC2(t) timers_toc(t, tic2)
INLINE static ticks timers_toc(unsigned int t, ticks tic) {
  ticks d = (getticks() - tic);
  atomic_add(&timers[t], d);
  return d;
}
#else
#define TIMER_TIC
#define TIMER_TOC(t)
#define TIMER_TIC2
#define TIMER_TOC2(t)
#endif

/* Function prototypes. */
void timers_reset(unsigned long long mask);

#endif /* SWIFT_TIMERS_H */
