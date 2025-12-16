/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Katy Proctor (katy.proctor@fysik.su.se)
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

#define _DOSELF1_SIDM(f) PASTE(runner_doself1_sidm, f)
#define DOSELF1_SIDM _DOSELF1_SIDM(FUNCTION)

#define _DOPAIR1_SIDM_NAIVE(f) PASTE(runner_dopair1_sidm_naive, f)
#define DOPAIR1_SIDM_NAIVE _DOPAIR1_SIDM_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_SIDM(f) PASTE(runner_dopair1_subset_sidm, f)
#define DOPAIR1_SUBSET_SIDM _DOPAIR1_SUBSET_SIDM(FUNCTION)

#define _DOPAIR1_SUBSET_SIDM_NAIVE(f) PASTE(runner_dopair1_subset_sidm_naive, f)
#define DOPAIR1_SUBSET_SIDM_NAIVE _DOPAIR1_SUBSET_SIDM_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_SIDM(f) PASTE(runner_doself1_subset_sidm, f)
#define DOSELF1_SUBSET_SIDM _DOSELF1_SUBSET_SIDM(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_SIDM(f) \
  PASTE(runner_doself1_subset_branch_sidm, f)
#define DOSELF1_SUBSET_BRANCH_SIDM _DOSELF1_SUBSET_BRANCH_SIDM(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_SIDM(f) \
  PASTE(runner_dopair1_subset_branch_sidm, f)
#define DOPAIR1_SUBSET_BRANCH_SIDM _DOPAIR1_SUBSET_BRANCH_SIDM(FUNCTION)

#define _DOSELF1_BRANCH_SIDM(f) PASTE(runner_doself1_branch_sidm, f)
#define DOSELF1_BRANCH_SIDM _DOSELF1_BRANCH_SIDM(FUNCTION)

#define _DOPAIR1_BRANCH_SIDM(f) PASTE(runner_dopair1_branch_sidm, f)
#define DOPAIR1_BRANCH_SIDM _DOPAIR1_BRANCH_SIDM(FUNCTION)

#define _DOSUB_PAIR1_SIDM(f) PASTE(runner_dosub_pair_sidm, f)
#define DOSUB_PAIR1_SIDM _DOSUB_PAIR1_SIDM(FUNCTION)

#define _DOSUB_SELF1_SIDM(f) PASTE(runner_dosub_self_sidm, f)
#define DOSUB_SELF1_SIDM _DOSUB_SELF1_SIDM(FUNCTION)

#define _DOSUB_SELF1_SUBSET_SIDM(f) PASTE(runner_dosub_self_subset_sidm, f)
#define DOSUB_SELF1_SUBSET_SIDM _DOSUB_SELF1_SUBSET_SIDM(FUNCTION)

#define _DOSUB_PAIR1_SUBSET_SIDM(f) PASTE(runner_dosub_pair_subset_sidm, f)
#define DOSUB_PAIR1_SUBSET_SIDM _DOSUB_PAIR1_SUBSET_SIDM(FUNCTION)

#define _TIMER_DOSELF_SIDM(f) PASTE(timer_doself_sidm, f)
#define TIMER_DOSELF_SIDM _TIMER_DOSELF_SIDM(FUNCTION)

#define _TIMER_DOPAIR_SIDM(f) PASTE(timer_dopair_sidm, f)
#define TIMER_DOPAIR_SIDM _TIMER_DOPAIR_SIDM(FUNCTION)

#define _TIMER_DOSUB_SELF_SIDM(f) PASTE(timer_dosub_self_sidm, f)
#define TIMER_DOSUB_SELF_SIDM _TIMER_DOSUB_SELF_SIDM(FUNCTION)

#define _TIMER_DOSUB_PAIR_SIDM(f) PASTE(timer_dosub_pair_sidm, f)
#define TIMER_DOSUB_PAIR_SIDM _TIMER_DOSUB_PAIR_SIDM(FUNCTION)

#define _FIND_SUB_SIDM(f) PASTE(runner_find_sub_sidm, f)
#define FIND_SUB_SIDM _FIND_SUB_SIDM(FUNCTION)

#define _IACT_NONSYM_SIDM(f) PASTE(runner_iact_nonsym_sidm, f)
#define IACT_NONSYM_SIDM _IACT_NONSYM_SIDM(FUNCTION)

#define _IACT_SIDM(f) PASTE(runner_iact_sidm, f)
#define IACT_SIDM _IACT_SIDM(FUNCTION)

void DOSELF1_BRANCH_SIDM(struct runner *r, const struct cell *c,
                         const int limit_min_h, const int limit_max_h);
void DOPAIR1_BRANCH_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                         const int limit_min_h, const int limit_max_h);

void DOSUB_SELF1_SIDM(struct runner *r, struct cell *c, int recurse_below_h_max,
                      const int gettimer);
void DOSUB_PAIR1_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                      int recurse_below_h_max, const int gettimer);

void DOSUB_SELF1_SUBSET_SIDM(struct runner *r, struct cell *ci,
                             struct sipart *siparts, const int *ind,
                             const int sicount, const int gettimer);
void DOSUB_PAIR1_SUBSET_SIDM(struct runner *r, struct cell *ci,
                             struct sipart *siparts, const int *ind,
                             const int sicount, struct cell *cj,
                             const int gettimer);

void DOSELF1_SUBSET_BRANCH_SIDM(struct runner *r, const struct cell *ci,
                                struct sipart *restrict siparts, const int *ind,
                                const int sicount);
void DOPAIR1_SUBSET_BRANCH_SIDM(struct runner *r,
                                const struct cell *restrict ci,
                                struct sipart *restrict siparts_i,
                                const int *ind, const int sicount,
                                struct cell *restrict cj);
