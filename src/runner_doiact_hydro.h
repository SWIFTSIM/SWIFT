/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#define _DOPAIR_SUBSET(f) PASTE(runner_dopair_subset, f)
#define DOPAIR_SUBSET _DOPAIR_SUBSET(FUNCTION)

#define _DOPAIR_SUBSET_BRANCH(f) PASTE(runner_dopair_subset_branch, f)
#define DOPAIR_SUBSET_BRANCH _DOPAIR_SUBSET_BRANCH(FUNCTION)

#define _DOPAIR_SUBSET_NOSORT(f) PASTE(runner_dopair_subset_nosort, f)
#define DOPAIR_SUBSET_NOSORT _DOPAIR_SUBSET_NOSORT(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE(f) PASTE(runner_dopair_subset_naive, f)
#define DOPAIR_SUBSET_NAIVE _DOPAIR_SUBSET_NAIVE(FUNCTION)

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

#define _DOSELF_SUBSET(f) PASTE(runner_doself_subset, f)
#define DOSELF_SUBSET _DOSELF_SUBSET(FUNCTION)

#define _DOSELF_SUBSET_BRANCH(f) PASTE(runner_doself_subset_branch, f)
#define DOSELF_SUBSET_BRANCH _DOSELF_SUBSET_BRANCH(FUNCTION)

#define _DOSUB_SELF1(f) PASTE(runner_dosub_self1, f)
#define DOSUB_SELF1 _DOSUB_SELF1(FUNCTION)

#define _DOSUB_PAIR1(f) PASTE(runner_dosub_pair1, f)
#define DOSUB_PAIR1 _DOSUB_PAIR1(FUNCTION)

#define _DOSUB_SELF2(f) PASTE(runner_dosub_self2, f)
#define DOSUB_SELF2 _DOSUB_SELF2(FUNCTION)

#define _DOSUB_PAIR2(f) PASTE(runner_dosub_pair2, f)
#define DOSUB_PAIR2 _DOSUB_PAIR2(FUNCTION)

#define _DOSUB_SUBSET(f) PASTE(runner_dosub_subset, f)
#define DOSUB_SUBSET _DOSUB_SUBSET(FUNCTION)

#define _IACT_NONSYM(f) PASTE(runner_iact_nonsym, f)
#define IACT_NONSYM _IACT_NONSYM(FUNCTION)

#define _IACT(f) PASTE(runner_iact, f)
#define IACT _IACT(FUNCTION)

#if ((FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY) ||  \
     (FUNCTION_TASK_LOOP == TASK_LOOP_GRADIENT) || \
     (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE))

#define _IACT_NONSYM_MHD(f) PASTE(runner_iact_nonsym_mhd, f)
#define IACT_NONSYM_MHD _IACT_NONSYM_MHD(FUNCTION)

#define _IACT_MHD(f) PASTE(runner_iact_mhd, f)
#define IACT_MHD _IACT_MHD(FUNCTION)

#define GET_MU0() \
  const double mu_0 = e->physical_constants->const_vacuum_permeability;

#else

#define IACT_NONSYM_MHD(x1, x2, x3, x4, x5, x6, x7, x8, x9) \
  {}
#define IACT_MHD(x1, x2, x3, x4, x5, x6, x7, x8, x9) \
  {}
#define GET_MU0() \
  {}
#endif

#if (FUNCTION_TASK_LOOP == TASK_LOOP_RT_GRADIENT) || \
    (FUNCTION_TASK_LOOP == TASK_LOOP_RT_TRANSPORT)
/* RT specific function calls */
#define PART_IS_ACTIVE part_is_rt_active
#define CELL_IS_ACTIVE cell_is_rt_active
#define CELL_ARE_PART_DRIFTED cell_are_part_drifted_rt_sub_cycle
#else
/* default hydro behaviour. */
#define PART_IS_ACTIVE part_is_active
#define CELL_IS_ACTIVE cell_is_active_hydro
#define CELL_ARE_PART_DRIFTED cell_are_part_drifted
/* when running with RT subcycling, we can have RT active
 * particles in a normal swift step that aren't drifted to
 * the current time, so we don't do those checks there. */
#define DO_DRIFT_DEBUG_CHECKS 1
#endif

#define _IACT_NONSYM_VEC(f) PASTE(runner_iact_nonsym_vec, f)
#define IACT_NONSYM_VEC _IACT_NONSYM_VEC(FUNCTION)

#define _IACT_VEC(f) PASTE(runner_iact_vec, f)
#define IACT_VEC _IACT_VEC(FUNCTION)

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
void DOSELF2_BRANCH(struct runner *r, struct cell *c);

void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj);
void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_SELF2(struct runner *r, struct cell *ci, int gettimer);

void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer);
void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer);

void DOSELF_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct part *restrict parts, int *restrict ind,
                          int count);

void DOPAIR_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct part *restrict parts_i, int *restrict ind,
                          int count, struct cell *restrict cj);

void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                  int *ind, int count, struct cell *cj, int gettimer);
