/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
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

#define _DOSELF1_BH(f) PASTE(runner_doself_bh, f)
#define DOSELF1_BH _DOSELF1_BH(FUNCTION)

#define _DO_SYM_PAIR1_BH(f) PASTE(runner_do_sym_pair_bh, f)
#define DO_SYM_PAIR1_BH _DO_SYM_PAIR1_BH(FUNCTION)

#define _DO_NONSYM_PAIR1_BH_NAIVE(f) PASTE(runner_do_nonsym_pair_bh_naive, f)
#define DO_NONSYM_PAIR1_BH_NAIVE _DO_NONSYM_PAIR1_BH_NAIVE(FUNCTION)

#define _DOPAIR1_BH_NAIVE(f) PASTE(runner_dopair_bh_naive, f)
#define DOPAIR1_BH_NAIVE _DOPAIR1_BH_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_BH(f) PASTE(runner_dopair_subset_bh, f)
#define DOPAIR1_SUBSET_BH _DOPAIR1_SUBSET_BH(FUNCTION)

#define _DOPAIR1_SUBSET_BH_NAIVE(f) PASTE(runner_dopair_subset_bh_naive, f)
#define DOPAIR1_SUBSET_BH_NAIVE _DOPAIR1_SUBSET_BH_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_BH(f) PASTE(runner_doself_subset_bh, f)
#define DOSELF1_SUBSET_BH _DOSELF1_SUBSET_BH(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_BH(f) PASTE(runner_doself_subset_branch_bh, f)
#define DOSELF1_SUBSET_BRANCH_BH _DOSELF1_SUBSET_BRANCH_BH(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_BH(f) PASTE(runner_dopair_subset_branch_bh, f)
#define DOPAIR1_SUBSET_BRANCH_BH _DOPAIR1_SUBSET_BRANCH_BH(FUNCTION)

#define _DOSUB_SELF_SUBSET_BH(f) PASTE(runner_dosub_self_subset_bh, f)
#define DOSUB_SELF_SUBSET_BH _DOSUB_SELF_SUBSET_BH(FUNCTION)

#define _DOSUB_PAIR_SUBSET_BH(f) PASTE(runner_dosub_pair_subset_bh, f)
#define DOSUB_PAIR_SUBSET_BH _DOSUB_PAIR_SUBSET_BH(FUNCTION)

#define _FIND_SUB_BH(f) PASTE(runner_find_sub_bh, f)
#define FIND_SUB_BH _FIND_SUB_BH(FUNCTION)
  
#define _DOSELF1_BRANCH_BH(f) PASTE(runner_doself_branch_bh, f)
#define DOSELF1_BRANCH_BH _DOSELF1_BRANCH_BH(FUNCTION)

#define _DOPAIR1_BRANCH_BH(f) PASTE(runner_dopair_branch_bh, f)
#define DOPAIR1_BRANCH_BH _DOPAIR1_BRANCH_BH(FUNCTION)

#define _DOSUB_PAIR1_BH(f) PASTE(runner_dosub_pair_bh, f)
#define DOSUB_PAIR1_BH _DOSUB_PAIR1_BH(FUNCTION)

#define _DOSUB_SELF1_BH(f) PASTE(runner_dosub_self_bh, f)
#define DOSUB_SELF1_BH _DOSUB_SELF1_BH(FUNCTION)

#define _TIMER_DOSELF_BH(f) PASTE(timer_doself_bh, f)
#define TIMER_DOSELF_BH _TIMER_DOSELF_BH(FUNCTION)

#define _TIMER_DOPAIR_BH(f) PASTE(timer_dopair_bh, f)
#define TIMER_DOPAIR_BH _TIMER_DOPAIR_BH(FUNCTION)

#define _TIMER_DOSUB_SELF_BH(f) PASTE(timer_dosub_self_bh, f)
#define TIMER_DOSUB_SELF_BH _TIMER_DOSUB_SELF_BH(FUNCTION)

#define _TIMER_DOSUB_PAIR_BH(f) PASTE(timer_dosub_pair_bh, f)
#define TIMER_DOSUB_PAIR_BH _TIMER_DOSUB_PAIR_BH(FUNCTION)

#define _IACT_BH_GAS(f) PASTE(runner_iact_nonsym_bh_gas, f)
#define IACT_BH_GAS _IACT_BH_GAS(FUNCTION)

#define _IACT_BH_BH(f) PASTE(runner_iact_nonsym_bh_bh, f)
#define IACT_BH_BH _IACT_BH_BH(FUNCTION)

void DOSELF1_BRANCH_BH(struct runner *r, struct cell *c, const int limit_min,
                       const int limit_max);
void DOPAIR1_BRANCH_BH(struct runner *r, struct cell *ci, struct cell *cj,
                       const int limit_min, const int limit_max);

void DOSUB_SELF1_BH(struct runner *r, struct cell *ci, int recurse_below_h_max,
                    const int gettimer);
void DOSUB_PAIR1_BH(struct runner *r, struct cell *ci, struct cell *cj,
                    int recurse_below_h_max, const int gettimer);

void DOSUB_SUBSET_BH(struct runner *r, struct cell *ci, struct bpart *bparts,
                     int *ind, const int bcount, struct cell *cj, int gettimer);

void DOSUB_SELF_SUBSET_BH(struct runner *r, struct cell *ci,
			  struct bpart *bparts, const int *ind,
			  const int bcount, const int gettimer);
void DOSUB_PAIR_SUBSET_BH(struct runner *r, struct cell *ci,
			  struct bpart *bparts, const int *ind,
			  const int bcount, struct cell *cj,
			  const int gettimer);

