/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOPAIR1_BRANCH(f) PASTE(runner_dopair1_branch, f)
#define DOPAIR1_BRANCH _DOPAIR1_BRANCH(FUNCTION)

#define _DOPAIR1(f) PASTE(runner_dopair1, f)
#define DOPAIR1 _DOPAIR1(FUNCTION)

#define _DOPAIR2_BRANCH(f) PASTE(runner_dopair2_branch, f)
#define DOPAIR2_BRANCH _DOPAIR2_BRANCH(FUNCTION)

#define _DOPAIR2(f) PASTE(runner_dopair2, f)
#define DOPAIR2 _DOPAIR2(FUNCTION)

#define _DOPAIR1_NAIVE(f) PASTE(runner_dopair1_naive, f)
#define DOPAIR1_NAIVE _DOPAIR1_NAIVE(FUNCTION)

#define _DOPAIR2_NAIVE(f) PASTE(runner_dopair2_naive, f)
#define DOPAIR2_NAIVE _DOPAIR2_NAIVE(FUNCTION)

#define _DOSELF1_NAIVE(f) PASTE(runner_doself1_naive, f)
#define DOSELF1_NAIVE _DOSELF1_NAIVE(FUNCTION)

#define _DOSELF2_NAIVE(f) PASTE(runner_doself2_naive, f)
#define DOSELF2_NAIVE _DOSELF2_NAIVE(FUNCTION)

#define _DOSELF1_BRANCH(f) PASTE(runner_doself1_branch, f)
#define DOSELF1_BRANCH _DOSELF1_BRANCH(FUNCTION)

#define _DOSELF1(f) PASTE(runner_doself1, f)
#define DOSELF1 _DOSELF1(FUNCTION)

#define _DOSELF2_BRANCH(f) PASTE(runner_doself2_branch, f)
#define DOSELF2_BRANCH _DOSELF2_BRANCH(FUNCTION)

#define _DOSELF2(f) PASTE(runner_doself2, f)
#define DOSELF2 _DOSELF2(FUNCTION)

#define _DOSUB_SELF1(f) PASTE(runner_dosub_self1, f)
#define DOSUB_SELF1 _DOSUB_SELF1(FUNCTION)

#define _DOSUB_PAIR1(f) PASTE(runner_dosub_pair1, f)
#define DOSUB_PAIR1 _DOSUB_PAIR1(FUNCTION)

#define _DOSUB_SELF2(f) PASTE(runner_dosub_self2, f)
#define DOSUB_SELF2 _DOSUB_SELF2(FUNCTION)

#define _DOSUB_PAIR2(f) PASTE(runner_dosub_pair2, f)
#define DOSUB_PAIR2 _DOSUB_PAIR2(FUNCTION)

#define _IACT_NONSYM(f) PASTE(runner_iact_nonsym, f)
#define IACT_NONSYM _IACT_NONSYM(FUNCTION)

#define _IACT(f) PASTE(runner_iact, f)
#define IACT _IACT(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

#define _TIMER_DOSUB_SELF(f) PASTE(timer_dosub_self, f)
#define TIMER_DOSUB_SELF _TIMER_DOSUB_SELF(FUNCTION)

#define _TIMER_DOSUB_PAIR(f) PASTE(timer_dosub_pair, f)
#define TIMER_DOSUB_PAIR _TIMER_DOSUB_PAIR(FUNCTION)

void DOSELF1_BRANCH(struct runner *r, struct cell *c);

void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1(struct runner *r, struct cell *ci, int gettimer);

void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer);
