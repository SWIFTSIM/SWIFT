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

#ifndef SWIFT_RUNNER_DOIACT_DARK_MATTER_H
#define SWIFT_RUNNER_DOIACT_DARK_MATTER_H

#define PASTE(x, y) x##_##y

#define _DOSELF1_NAIVE(f) PASTE(runner_doself1_naive, f)
#define DOSELF1_NAIVE _DOSELF1_NAIVE(FUNCTION)

#define _DOSELF1(f) PASTE(runner_doself1, f)
#define DOSELF1 _DOSELF1(FUNCTION)

#define _DOSELF1_BRANCH(f) PASTE(runner_doself1_branch, f)
#define DOSELF1_BRANCH _DOSELF1_BRANCH(FUNCTION)

#define _DOPAIR1(f) PASTE(runner_dopair1, f)
#define DOPAIR1 _DOPAIR1(FUNCTION)

#define _DOPAIR1_NAIVE(f) PASTE(runner_dopair1_naive, f)
#define DOPAIR1_NAIVE _DOPAIR1_NAIVE(FUNCTION)

#define _DOPAIR1_BRANCH(f) PASTE(runner_dopair1_branch, f)
#define DOPAIR1_BRANCH _DOPAIR1_BRANCH(FUNCTION)

#define _DOSUB_SELF1(f) PASTE(runner_dosub_self1, f)
#define DOSUB_SELF1 _DOSUB_SELF1(FUNCTION)

#define _DOSUB_PAIR1(f) PASTE(runner_dosub_pair1, f)
#define DOSUB_PAIR1 _DOSUB_PAIR1(FUNCTION)

#define _DOSELF_SUBSET(f) PASTE(runner_doself_subset, f)
#define DOSELF_SUBSET _DOSELF_SUBSET(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE(f) PASTE(runner_dopair_subset_naive, f)
#define DOPAIR_SUBSET_NAIVE _DOPAIR_SUBSET_NAIVE(FUNCTION)

#define _DOPAIR_SUBSET_BRANCH(f) PASTE(runner_dopair_subset_branch, f)
#define DOPAIR_SUBSET_BRANCH _DOPAIR_SUBSET_BRANCH(FUNCTION)

#define _DOSUB_SUBSET(f) PASTE(runner_dosub_subset, f)
#define DOSUB_SUBSET _DOSUB_SUBSET(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

#define _TIMER_DOSUB_SELF(f) PASTE(timer_dosub_self, f)
#define TIMER_DOSUB_SELF _TIMER_DOSUB_SELF(FUNCTION)

#define _TIMER_DOSUB_PAIR(f) PASTE(timer_dosub_pair, f)
#define TIMER_DOSUB_PAIR _TIMER_DOSUB_PAIR(FUNCTION)

#define _TIMER_DOSELF_SUBSET(f) PASTE(timer_doself_subset, f)
#define TIMER_DOSELF_SUBSET _TIMER_DOSELF_SUBSET(FUNCTION)

#define _TIMER_DOPAIR_SUBSET(f) PASTE(timer_dopair_subset, f)
#define TIMER_DOPAIR_SUBSET _TIMER_DOPAIR_SUBSET(FUNCTION)

void DOSELF1_BRANCH(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj);
void DOSUB_SELF1(struct runner *r, struct cell *ci);
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj);
void DOSELF_SUBSET(struct runner *r, struct cell *restrict ci,
                          struct dmpart *restrict dmparts, int *restrict ind,
                          int count);

void DOPAIR_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct dmpart *restrict dmparts_i, int *restrict ind,
                          struct cell *restrict cj, int count);

void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct dmpart *dmparts,
                  int *ind, int count, struct cell *cj);

#endif
