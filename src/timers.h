/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "cycle.h"
#include "inline.h"

/* The timers themselves. */
enum {
  timer_none = 0,
  timer_prepare,
  timer_init,
  timer_drift,
  timer_kick,
  timer_dosort,
  timer_doself_density,
  timer_doself_force,
  timer_doself_grav,
  timer_dopair_density,
  timer_dopair_force,
  timer_dopair_grav,
  timer_dograv_external,
  timer_dosub_density,
  timer_dosub_force,
  timer_dosub_grav,
  timer_dopair_subset,
  timer_doghost,
  timer_gettask,
  timer_qget,
  timer_qsteal,
  timer_runners,
  timer_step,
  timer_count,
};

/* The timers. */
extern ticks timers[timer_count];

/* Mask for all timers. */
#define timers_mask_all ((1 << timer_count) - 1)

/* Define the timer macros. */
#ifdef TIMER_VERBOSE
#ifndef TIMER
#define TIMER
#endif
#endif
#ifdef TIMER
#define TIMER_TIC_ND tic = getticks();
#define TIMER_TIC2_ND ticks tic2 = getticks();
#define TIMER_TIC ticks tic = getticks();
#define TIMER_TOC(t) timers_toc(t, tic)
#define TIMER_TIC2 ticks tic2 = getticks();
#define TIMER_TOC2(t) timers_toc(t, tic2)
INLINE static ticks timers_toc(int t, ticks tic) {
  ticks d = (getticks() - tic);
  __sync_add_and_fetch(&timers[t], d);
  return d;
}
#else
#define TIMER_TIC
#define TIMER_TOC(t)
#define TIMER_TIC2
#define TIMER_TOC2(t)
#endif

/* Function prototypes. */
void timers_reset(unsigned int mask);

#endif /* SWIFT_TIMERS_H */
