/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#define _DOSELF1_STARS(f) PASTE(runner_doself_stars, f)
#define DOSELF1_STARS _DOSELF1_STARS(FUNCTION)

#define _DO_SYM_PAIR1_STARS(f) PASTE(runner_do_sym_pair_stars, f)
#define DO_SYM_PAIR1_STARS _DO_SYM_PAIR1_STARS(FUNCTION)

#define _DO_NONSYM_PAIR1_STARS_NAIVE(f) \
  PASTE(runner_do_nonsym_pair_stars_naive, f)
#define DO_NONSYM_PAIR1_STARS_NAIVE _DO_NONSYM_PAIR1_STARS_NAIVE(FUNCTION)

#define _DOPAIR1_STARS_NAIVE(f) PASTE(runner_dopair_stars_naive, f)
#define DOPAIR1_STARS_NAIVE _DOPAIR1_STARS_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_STARS(f) PASTE(runner_dopair_subset_stars, f)
#define DOPAIR1_SUBSET_STARS _DOPAIR1_SUBSET_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_STARS_NAIVE(f) \
  PASTE(runner_dopair_subset_stars_naive, f)
#define DOPAIR1_SUBSET_STARS_NAIVE _DOPAIR1_SUBSET_STARS_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_STARS(f) PASTE(runner_doself_subset_stars, f)
#define DOSELF1_SUBSET_STARS _DOSELF1_SUBSET_STARS(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_doself_subset_branch_stars, f)
#define DOSELF1_SUBSET_BRANCH_STARS _DOSELF1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_dopair_subset_branch_stars, f)
#define DOPAIR1_SUBSET_BRANCH_STARS _DOPAIR1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOSUB_SUBSET_STARS(f) PASTE(runner_dosub_subset_stars, f)
#define DOSUB_SUBSET_STARS _DOSUB_SUBSET_STARS(FUNCTION)

#define _DOSELF1_BRANCH_STARS(f) PASTE(runner_doself_branch_stars, f)
#define DOSELF1_BRANCH_STARS _DOSELF1_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_BRANCH_STARS(f) PASTE(runner_dopair_branch_stars, f)
#define DOPAIR1_BRANCH_STARS _DOPAIR1_BRANCH_STARS(FUNCTION)

#define _DOSUB_PAIR1_STARS(f) PASTE(runner_dosub_pair_stars, f)
#define DOSUB_PAIR1_STARS _DOSUB_PAIR1_STARS(FUNCTION)

#define _DOSUB_SELF1_STARS(f) PASTE(runner_dosub_self_stars, f)
#define DOSUB_SELF1_STARS _DOSUB_SELF1_STARS(FUNCTION)

#define _TIMER_DOSELF_STARS(f) PASTE(timer_doself_stars, f)
#define TIMER_DOSELF_STARS _TIMER_DOSELF_STARS(FUNCTION)

#define _TIMER_DOPAIR_STARS(f) PASTE(timer_dopair_stars, f)
#define TIMER_DOPAIR_STARS _TIMER_DOPAIR_STARS(FUNCTION)

#define _TIMER_DOSUB_SELF_STARS(f) PASTE(timer_dosub_self_stars, f)
#define TIMER_DOSUB_SELF_STARS _TIMER_DOSUB_SELF_STARS(FUNCTION)

#define _TIMER_DOSUB_PAIR_STARS(f) PASTE(timer_dosub_pair_stars, f)
#define TIMER_DOSUB_PAIR_STARS _TIMER_DOSUB_PAIR_STARS(FUNCTION)

#define _IACT_STARS(f) PASTE(runner_iact_nonsym_stars, f)
#define IACT_STARS _IACT_STARS(FUNCTION)

void DOSELF1_BRANCH_STARS(struct runner *r, struct cell *c);
void DOPAIR1_BRANCH_STARS(struct runner *r, struct cell *ci, struct cell *cj);

void DOSUB_SELF1_STARS(struct runner *r, struct cell *ci, int gettimer);
void DOSUB_PAIR1_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                       int gettimer);

void DOSELF1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts,
                                 int *restrict ind, int scount);

void DOPAIR1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts_i,
                                 int *restrict ind, int scount,
                                 struct cell *restrict cj);

void DOSUB_SUBSET_STARS(struct runner *r, struct cell *ci, struct spart *sparts,
                        int *ind, int scount, struct cell *cj, int gettimer);
