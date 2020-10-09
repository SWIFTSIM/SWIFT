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

#define _DOSELF2_NAIVE(f) PASTE(runner_doself2_naive, f)
#define DOSELF2_NAIVE _DOSELF2_NAIVE(FUNCTION)

#define _DOSELF2_BRANCH(f) PASTE(runner_doself2_branch, f)
#define DOSELF2_BRANCH _DOSELF2_BRANCH(FUNCTION)

#define _DOSELF2(f) PASTE(runner_doself2, f)
#define DOSELF2 _DOSELF2(FUNCTION)

#define _DOPAIR2_NAIVE(f) PASTE(runner_dopair2_naive, f)
#define DOPAIR2_NAIVE _DOPAIR2_NAIVE(FUNCTION)

#define _DOPAIR2_BRANCH(f) PASTE(runner_dopair2_branch, f)
#define DOPAIR2_BRANCH _DOPAIR2_BRANCH(FUNCTION)

#define _DOPAIR2(f) PASTE(runner_dopair2, f)
#define DOPAIR2 _DOPAIR2(FUNCTION)

#define _DOSUB_SELF2(f) PASTE(runner_dosub_self2, f)
#define DOSUB_SELF2 _DOSUB_SELF2(FUNCTION)

#define _DOSUB_PAIR2(f) PASTE(runner_dosub_pair2, f)
#define DOSUB_PAIR2 _DOSUB_PAIR2(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

void DOSELF2_BRANCH(struct runner *r, struct cell *c);

void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF2(struct runner *r, struct cell *ci);

void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj);


#endif
