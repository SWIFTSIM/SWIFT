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

#define _DOPAIR1_BRANCH_HYDRO_SINKS(f) PASTE(runner_dopair1_branch_hydro_sinks, f)
#define DOPAIR1_BRANCH_HYDRO_SINKS _DOPAIR1_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOPAIR1_HYDRO_SINKS(f) PASTE(runner_dopair1_hydro_sinks, f)
#define DOPAIR1_HYDRO_SINKS _DOPAIR1_HYDRO_SINKS(FUNCTION)

#define _DOPAIR2_BRANCH_HYDRO_SINKS(f) PASTE(runner_dopair2_branch_hydro_sinks, f)
#define DOPAIR2_BRANCH_HYDRO_SINKS _DOPAIR2_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOPAIR2_HYDRO_SINKS(f) PASTE(runner_dopair2_hydro_sinks, f)
#define DOPAIR2_HYDRO_SINKS _DOPAIR2_HYDRO_SINKS(FUNCTION)

#define _DOPAIR_SUBSET_HYDRO_SINKS(f) PASTE(runner_dopair_subset_hydro_sinks, f)
#define DOPAIR_SUBSET_HYDRO_SINKS _DOPAIR_SUBSET_HYDRO_SINKS(FUNCTION)

#define _DOPAIR_SUBSET_BRANCH_HYDRO_SINKS(f) PASTE(runner_dopair_subset_branch_hydro_sinks, f)
#define DOPAIR_SUBSET_BRANCH_HYDRO_SINKS _DOPAIR_SUBSET_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOPAIR_SUBSET_NOSORT_HYDRO_SINKS(f) PASTE(runner_dopair_subset_nosort_hydro_sinks, f)
#define DOPAIR_SUBSET_NOSORT_HYDRO_SINKS _DOPAIR_SUBSET_NOSORT_HYDRO_SINKS(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE_HYDRO_SINKS(f) PASTE(runner_dopair_subset_naive_hydro_sinks, f)
#define DOPAIR_SUBSET_NAIVE_HYDRO_SINKS _DOPAIR_SUBSET_NAIVE_HYDRO_SINKS(FUNCTION)

#define _DOPAIR1_NAIVE_HYDRO_SINKS(f) PASTE(runner_dopair1_naive_hydro_sinks, f)
#define DOPAIR1_NAIVE_HYDRO_SINKS _DOPAIR1_NAIVE_HYDRO_SINKS(FUNCTION)

#define _DOPAIR2_NAIVE_HYDRO_SINKS(f) PASTE(runner_dopair2_naive_hydro_sinks, f)
#define DOPAIR2_NAIVE_HYDRO_SINKS _DOPAIR2_NAIVE_HYDRO_SINKS(FUNCTION)

#define _DOSELF1_NAIVE_HYDRO_SINKS(f) PASTE(runner_doself1_naive_hydro_sinks, f)
#define DOSELF1_NAIVE_HYDRO_SINKS _DOSELF1_NAIVE_HYDRO_SINKS(FUNCTION)

#define _DOSELF2_NAIVE_HYDRO_SINKS(f) PASTE(runner_doself2_naive_hydro_sinks, f)
#define DOSELF2_NAIVE_HYDRO_SINKS _DOSELF2_NAIVE_HYDRO_SINKS(FUNCTION)

#define _DOSELF1_BRANCH_HYDRO_SINKS(f) PASTE(runner_doself1_branch_hydro_sinks, f)
#define DOSELF1_BRANCH_HYDRO_SINKS _DOSELF1_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOSELF1_HYDRO_SINKS(f) PASTE(runner_doself1_hydro_sinks, f)
#define DOSELF1_HYDRO_SINKS _DOSELF1_HYDRO_SINKS(FUNCTION)

#define _DOSELF2_BRANCH_HYDRO_SINKS(f) PASTE(runner_doself2_branch_hydro_sinks, f)
#define DOSELF2_BRANCH_HYDRO_SINKS _DOSELF2_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOSELF2_HYDRO_SINKS(f) PASTE(runner_doself2_hydro_sinks, f)
#define DOSELF2_HYDRO_SINKS _DOSELF2_HYDRO_SINKS(FUNCTION)

#define _DOSELF_SUBSET_HYDRO_SINKS(f) PASTE(runner_doself_subset_hydro_sinks, f)
#define DOSELF_SUBSET_HYDRO_SINKS _DOSELF_SUBSET_HYDRO_SINKS(FUNCTION)

#define _DOSELF_SUBSET_BRANCH_HYDRO_SINKS(f) PASTE(runner_doself_subset_branch_hydro_sinks, f)
#define DOSELF_SUBSET_BRANCH_HYDRO_SINKS _DOSELF_SUBSET_BRANCH_HYDRO_SINKS(FUNCTION)

#define _DOSUB_SELF1_HYDRO_SINKS(f) PASTE(runner_dosub_self1_hydro_sinks, f)
#define DOSUB_SELF1_HYDRO_SINKS _DOSUB_SELF1_HYDRO_SINKS(FUNCTION)

#define _DOSUB_PAIR1_HYDRO_SINKS(f) PASTE(runner_dosub_pair1_hydro_sinks, f)
#define DOSUB_PAIR1_HYDRO_SINKS _DOSUB_PAIR1_HYDRO_SINKS(FUNCTION)

#define _DOSUB_SELF2_HYDRO_SINKS(f) PASTE(runner_dosub_self2_hydro_sinks, f)
#define DOSUB_SELF2_HYDRO_SINKS _DOSUB_SELF2_HYDRO_SINKS(FUNCTION)

#define _DOSUB_PAIR2_HYDRO_SINKS(f) PASTE(runner_dosub_pair2_hydro_sinks, f)
#define DOSUB_PAIR2_HYDRO_SINKS _DOSUB_PAIR2_HYDRO_SINKS(FUNCTION)

#define _DOSUB_SELF_SUBSET_HYDRO_SINKS(f) PASTE(runner_dosub_self_subset_hydro_sinks, f)
#define DOSUB_SELF_SUBSET_HYDRO_SINKS _DOSUB_SELF_SUBSET_HYDRO_SINKS(FUNCTION)

#define _DOSUB_PAIR_SUBSET_HYDRO_SINKS(f) PASTE(runner_dosub_pair_subset_hydro_sinks, f)
#define DOSUB_PAIR_SUBSET_HYDRO_SINKS _DOSUB_PAIR_SUBSET_HYDRO_SINKS(FUNCTION)

#define _FIND_SUB_HYDRO_SINKS(f) PASTE(runner_find_sub_hydro_sinks, f)
#define FIND_SUB_HYDRO_SINKS _FIND_SUB_HYDRO_SINKS(FUNCTION)

#define _IACT_NONSYM_HYDRO_SINKS(f) PASTE(runner_iact_nonsym_hydro_sinks, f)
#define IACT_NONSYM_HYDRO_SINKS _IACT_NONSYM_HYDRO_SINKS(FUNCTION)

#define _IACT_HYDRO_SINKS(f) PASTE(runner_iact_hydro_sinks, f)
#define IACT_HYDRO_SINKS _IACT_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOSELF_HYDRO_SINKS(f) PASTE(timer_doself_hydro_sinks, f)
#define TIMER_DOSELF_HYDRO_SINKS _TIMER_DOSELF_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOPAIR_HYDRO_SINKS(f) PASTE(timer_dopair_hydro_sinks, f)
#define TIMER_DOPAIR_HYDRO_SINKS _TIMER_DOPAIR_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOSUB_SELF_HYDRO_SINKS(f) PASTE(timer_dosub_self_hydro_sinks, f)
#define TIMER_DOSUB_SELF_HYDRO_SINKS _TIMER_DOSUB_SELF_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOSUB_PAIR_HYDRO_SINKS(f) PASTE(timer_dosub_pair_hydro_sinks, f)
#define TIMER_DOSUB_PAIR_HYDRO_SINKS _TIMER_DOSUB_PAIR_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOSELF_SUBSET_HYDRO_SINKS(f) PASTE(timer_doself_subset_hydro_sinks, f)
#define TIMER_DOSELF_SUBSET_HYDRO_SINKS _TIMER_DOSELF_SUBSET_HYDRO_SINKS(FUNCTION)

#define _TIMER_DOPAIR_SUBSET_HYDRO_SINKS(f) PASTE(timer_dopair_subset_hydro_sinks, f)
#define TIMER_DOPAIR_SUBSET_HYDRO_SINKS _TIMER_DOPAIR_SUBSET_HYDRO_SINKS(FUNCTION)

/* Active/drifted status macros for the sink formation loop. */
#define PART_IS_ACTIVE part_is_active
#define CELL_IS_ACTIVE cell_is_active_hydro
#define CELL_ARE_PART_DRIFTED cell_are_part_drifted
#define DO_DRIFT_DEBUG_CHECKS 1
#define GET_MU0() \
  {               \
  }

void DOSELF1_BRANCH_HYDRO_SINKS(struct runner *r, const struct cell *c,
                                 const float r_cut);
void DOSELF2_BRANCH_HYDRO_SINKS(struct runner *r, const struct cell *c,
                                 const float r_cut);

void DOPAIR1_BRANCH_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                 struct cell *cj, const float r_cut);
void DOPAIR2_BRANCH_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                 struct cell *cj, const float r_cut);

void DOSUB_SELF1_HYDRO_SINKS(struct runner *r, struct cell *c,
                              const float r_cut, const int gettimer);
void DOSUB_SELF2_HYDRO_SINKS(struct runner *r, struct cell *c,
                              const float r_cut, const int gettimer);

void DOSUB_PAIR1_HYDRO_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                              const float r_cut, const int gettimer);
void DOSUB_PAIR2_HYDRO_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                              const float r_cut, const int gettimer);

void DOSELF_SUBSET_BRANCH_HYDRO_SINKS(struct runner *r, const struct cell *ci,
                                       struct part *restrict parts,
                                       const int *ind, const int count,
                                       const float r_cut);

void DOPAIR_SUBSET_BRANCH_HYDRO_SINKS(struct runner *r,
                                       const struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *ind, const int count,
                                       struct cell *restrict cj,
                                       const float r_cut);

void DOSUB_PAIR_SUBSET_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                    struct part *parts, const int *ind,
                                    const int count, struct cell *cj,
                                    const float r_cut, const int gettimer);

void DOSUB_SELF_SUBSET_HYDRO_SINKS(struct runner *r, struct cell *ci,
                                    struct part *parts, const int *ind,
                                    const int count, const float r_cut,
                                    const int gettimer);
