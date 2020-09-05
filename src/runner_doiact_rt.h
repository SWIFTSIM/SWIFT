/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOSELF1_RT(f) PASTE(runner_doself_rt, f)
#define DOSELF1_RT _DOSELF1_RT(FUNCTION)

#define _DOPAIR1_RT(f) PASTE(runner_dopair_rt, f)
#define DOPAIR1_RT _DOPAIR1_RT(FUNCTION)

#define _DOPAIR1_NONSYM_RT(f) PASTE(runner_dopair_nonsym_rt, f)
#define DOPAIR1_NONSYM_RT _DOPAIR1_NONSYM_RT(FUNCTION)

#define _DOSELF1_BRANCH_RT(f) PASTE(runner_doself_branch_rt, f)
#define DOSELF1_BRANCH_RT _DOSELF1_BRANCH_RT(FUNCTION)

#define _DOPAIR1_BRANCH_RT(f) PASTE(runner_dopair_branch_rt, f)
#define DOPAIR1_BRANCH_RT _DOPAIR1_BRANCH_RT(FUNCTION)

#define _DOSUB_PAIR1_RT(f) PASTE(runner_dosub_pair_rt, f)
#define DOSUB_PAIR1_RT _DOSUB_PAIR1_RT(FUNCTION)

#define _DOSUB_SELF1_RT(f) PASTE(runner_dosub_self_rt, f)
#define DOSUB_SELF1_RT _DOSUB_SELF1_RT(FUNCTION)

#define _TIMER_DOSELF_RT(f) PASTE(timer_doself_rt, f)
#define TIMER_DOSELF_RT _TIMER_DOSELF_RT(FUNCTION)

#define _TIMER_DOPAIR_RT(f) PASTE(timer_dopair_rt, f)
#define TIMER_DOPAIR_RT _TIMER_DOPAIR_RT(FUNCTION)

#define _TIMER_DOSUB_SELF_RT(f) PASTE(timer_dosub_self_rt, f)
#define TIMER_DOSUB_SELF_RT _TIMER_DOSUB_SELF_RT(FUNCTION)

#define _TIMER_DOSUB_PAIR_RT(f) PASTE(timer_dosub_pair_rt, f)
#define TIMER_DOSUB_PAIR_RT _TIMER_DOSUB_PAIR_RT(FUNCTION)

#define _IACT_RT(f) PASTE(runner_iact_rt, f)
#define IACT_RT _IACT_RT(FUNCTION)

void DOSELF1_BRANCH_RT(struct runner *r, struct cell *c, int timer);
void DOPAIR1_BRANCH_RT(struct runner *r, struct cell *ci, struct cell *cj,
                       int timer);

void DOSUB_SELF1_RT(struct runner *r, struct cell *ci, int timer);
void DOSUB_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj,
                    int timer);
