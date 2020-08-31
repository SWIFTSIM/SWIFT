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
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOSELF1_DM(f) PASTE(runner_doself_dm, f)
#define DOSELF1_DM _DOSELF1_DM(FUNCTION)

#define _DO_SYM_PAIR1_DM(f) PASTE(runner_do_sym_pair_dm, f)
#define DO_SYM_PAIR1_DM _DO_SYM_PAIR1_DM(FUNCTION)

#define _DO_NONSYM_PAIR1_DM_NAIVE(f) PASTE(runner_do_nonsym_pair_dm_naive, f)
#define DO_NONSYM_PAIR1_DM_NAIVE _DO_NONSYM_PAIR1_DM_NAIVE(FUNCTION)

#define _DOPAIR1_DM_NAIVE(f) PASTE(runner_dopair_dm_naive, f)
#define DOPAIR1_DM_NAIVE _DOPAIR1_DM_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_DM(f) PASTE(runner_dopair_subset_dm, f)
#define DOPAIR1_SUBSET_DM _DOPAIR1_SUBSET_DM(FUNCTION)

#define _DOPAIR1_SUBSET_DM_NAIVE(f) PASTE(runner_dopair_subset_dm_naive, f)
#define DOPAIR1_SUBSET_DM_NAIVE _DOPAIR1_SUBSET_DM_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_DM(f) PASTE(runner_doself_subset_dm, f)
#define DOSELF1_SUBSET_DM _DOSELF1_SUBSET_DM(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_DM(f) PASTE(runner_doself_subset_branch_dm, f)
#define DOSELF1_SUBSET_BRANCH_DM _DOSELF1_SUBSET_BRANCH_DM(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_DM(f) PASTE(runner_dopair_subset_branch_dm, f)
#define DOPAIR1_SUBSET_BRANCH_DM _DOPAIR1_SUBSET_BRANCH_DM(FUNCTION)

#define _DOSUB_SUBSET_DM(f) PASTE(runner_dosub_subset_dm, f)
#define DOSUB_SUBSET_DM _DOSUB_SUBSET_DM(FUNCTION)

#define _DOSELF1_BRANCH_DM(f) PASTE(runner_doself_branch_dm, f)
#define DOSELF1_BRANCH_DM _DOSELF1_BRANCH_DM(FUNCTION)

#define _DOPAIR1_BRANCH_DM(f) PASTE(runner_dopair_branch_dm, f)
#define DOPAIR1_BRANCH_DM _DOPAIR1_BRANCH_DM(FUNCTION)

#define _DOSUB_PAIR1_DM(f) PASTE(runner_dosub_pair_dm, f)
#define DOSUB_PAIR1_DM _DOSUB_PAIR1_DM(FUNCTION)

#define _DOSUB_SELF1_DM(f) PASTE(runner_dosub_self_dm, f)
#define DOSUB_SELF1_DM _DOSUB_SELF1_DM(FUNCTION)

#define _TIMER_DOSELF_DM(f) PASTE(timer_doself_dm, f)
#define TIMER_DOSELF_DM _TIMER_DOSELF_DM(FUNCTION)

#define _TIMER_DOPAIR_DM(f) PASTE(timer_dopair_dm, f)
#define TIMER_DOPAIR_DM _TIMER_DOPAIR_DM(FUNCTION)

#define _TIMER_DOSUB_SELF_DM(f) PASTE(timer_dosub_self_dm, f)
#define TIMER_DOSUB_SELF_DM _TIMER_DOSUB_SELF_DM(FUNCTION)

#define _TIMER_DOSUB_PAIR_DM(f) PASTE(timer_dosub_pair_dm, f)
#define TIMER_DOSUB_PAIR_DM _TIMER_DOSUB_PAIR_DM(FUNCTION)

void DOSELF1_BRANCH_DM(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH_DM(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1_DM(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_PAIR1_DM(struct runner *r, struct cell *ci, struct cell *cj,
                    int gettimer);

void DOSELF1_SUBSET_BRANCH_DM(struct runner *r, struct cell *restrict ci,
                              struct dmpart *restrict dmparts, int *restrict ind,
                              const int dmcount);
void DOPAIR1_SUBSET_BRANCH_DM(struct runner *r, struct cell *restrict ci,
                              struct dmpart *restrict dmparts_i,
                              int *restrict ind, int const dmcount,
                              struct cell *restrict cj);

void DOSUB_SUBSET_DM(struct runner *r, struct cell *ci, struct dmpart *dmparts,
                     int *ind, const int dmcount, struct cell *cj, int gettimer);
