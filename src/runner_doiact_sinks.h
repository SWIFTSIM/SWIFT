/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#define _DOSELF1_SINK(f) PASTE(runner_doself_sinks, f)
#define DOSELF1_SINK _DOSELF1_SINK(FUNCTION)

#define _DO_SYM_PAIR1_SINK(f) PASTE(runner_do_sym_pair_sinks, f)
#define DO_SYM_PAIR1_SINK _DO_SYM_PAIR1_SINK(FUNCTION)

#define _DO_NONSYM_PAIR1_SINK_NAIVE(f) PASTE(runner_do_nonsym_pair_sinks_naive, f)
#define DO_NONSYM_PAIR1_SINK_NAIVE _DO_NONSYM_PAIR1_SINK_NAIVE(FUNCTION)

#define _DOPAIR1_SINK_NAIVE(f) PASTE(runner_dopair_sinks_naive, f)
#define DOPAIR1_SINK_NAIVE _DOPAIR1_SINK_NAIVE(FUNCTION)

#define _DOSELF1_BRANCH_SINK(f) PASTE(runner_doself_branch_sinks, f)
#define DOSELF1_BRANCH_SINK _DOSELF1_BRANCH_SINK(FUNCTION)

#define _DOPAIR1_BRANCH_SINK(f) PASTE(runner_dopair_branch_sinks, f)
#define DOPAIR1_BRANCH_SINK _DOPAIR1_BRANCH_SINK(FUNCTION)

#define _DOSUB_PAIR1_SINK(f) PASTE(runner_dosub_pair_sinks, f)
#define DOSUB_PAIR1_SINK _DOSUB_PAIR1_SINK(FUNCTION)

#define _DOSUB_SELF1_SINK(f) PASTE(runner_dosub_self_sinks, f)
#define DOSUB_SELF1_SINK _DOSUB_SELF1_SINK(FUNCTION)

#define _TIMER_DOSELF_SINK(f) PASTE(timer_doself_sinks, f)
#define TIMER_DOSELF_SINK _TIMER_DOSELF_SINK(FUNCTION)

#define _TIMER_DOPAIR_SINK(f) PASTE(timer_dopair_sinks, f)
#define TIMER_DOPAIR_SINK _TIMER_DOPAIR_SINK(FUNCTION)

#define _TIMER_DOSUB_SELF_SINK(f) PASTE(timer_dosub_self_sinks, f)
#define TIMER_DOSUB_SELF_SINK _TIMER_DOSUB_SELF_SINK(FUNCTION)

#define _TIMER_DOSUB_PAIR_SINK(f) PASTE(timer_dosub_pair_sinks, f)
#define TIMER_DOSUB_PAIR_SINK _TIMER_DOSUB_PAIR_SINK(FUNCTION)

#define _IACT_SINK_GAS(f) PASTE(runner_iact_nonsym_sinks_gas, f)
#define IACT_SINK_GAS _IACT_SINK_GAS(FUNCTION)

#define _IACT_SINK_SINK(f) PASTE(runner_iact_nonsym_sinks_sink, f)
#define IACT_SINK_SINK _IACT_SINK_SINK(FUNCTION)

/* The sinks don't have the "interaction subset" stuff that the BHs do, but someday? */

// #define _DOPAIR1_SUBSET_SINK(f) PASTE(runner_dopair_subset_sinks, f)
// #define DOPAIR1_SUBSET_SINK _DOPAIR1_SUBSET_SINK(FUNCTION)

// #define _DOPAIR1_SUBSET_SINK_NAIVE(f) PASTE(runner_dopair_subset_sinks_naive, f)
// #define DOPAIR1_SUBSET_SINK_NAIVE _DOPAIR1_SUBSET_SINK_NAIVE(FUNCTION)

// #define _DOSELF1_SUBSET_SINK(f) PASTE(runner_doself_subset_sinks, f)
// #define DOSELF1_SUBSET_SINK _DOSELF1_SUBSET_SINK(FUNCTION)

// #define _DOSELF1_SUBSET_BRANCH_SINK(f) PASTE(runner_doself_subset_branch_sinks, f)
// #define DOSELF1_SUBSET_BRANCH_SINK _DOSELF1_SUBSET_BRANCH_SINK(FUNCTION)

// #define _DOPAIR1_SUBSET_BRANCH_SINK(f) PASTE(runner_dopair_subset_branch_sinks, f)
// #define DOPAIR1_SUBSET_BRANCH_SINK _DOPAIR1_SUBSET_BRANCH_SINK(FUNCTION)

// #define _DOSUB_SUBSET_SINK(f) PASTE(runner_dosub_subset_sinks, f)
// #define DOSUB_SUBSET_SINK _DOSUB_SUBSET_SINK(FUNCTION)

// void DOSELF1_SUBSET_BRANCH_SINK(struct runner *r, struct cell *restrict ci,
//                               struct bpart *restrict bparts, int *restrict ind,
//                               const int bcount);

// void DOPAIR1_SUBSET_BRANCH_SINK(struct runner *r, struct cell *restrict ci,
//                               struct bpart *restrict bparts_i,
//                               int *restrict ind, int const bcount,
//                               struct cell *restrict cj);

// void DOSUB_SUBSET_SINK(struct runner *r, struct cell *ci, struct bpart *bparts,
//                      int *ind, const int bcount, struct cell *cj, int gettimer);


void DOSELF1_BRANCH_SINK(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH_SINK(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1_SINK(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_PAIR1_SINK(struct runner *r, struct cell *ci, struct cell *cj,
                    int gettimer);

