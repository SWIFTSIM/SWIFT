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

#define _DOPAIR_SUBSET_DM_NAIVE(f) PASTE(runner_dopair_subset_dm_naive, f)
#define DOPAIR_SUBSET_DM_NAIVE _DOPAIR_SUBSET_DM_NAIVE(FUNCTION)

#define _DOPAIR1_DM_NAIVE(f) PASTE(runner_dopair1_dm_naive, f)
#define DOPAIR1_DM_NAIVE _DOPAIR1_DM_NAIVE(FUNCTION)

#define _DOPAIR2_DM_NAIVE(f) PASTE(runner_dopair2_dm_naive, f)
#define DOPAIR2_DM_NAIVE _DOPAIR2_DM_NAIVE(FUNCTION)

#define _DOSELF1_DM_NAIVE(f) PASTE(runner_doself1_dm_naive, f)
#define DOSELF1_DM_NAIVE _DOSELF1_DM_NAIVE(FUNCTION)

#define _DOSELF2_DM_NAIVE(f) PASTE(runner_doself2_dm_naive, f)
#define DOSELF2_DM_NAIVE _DOSELF2_DM_NAIVE(FUNCTION)

#define _DOSELF_SUBSET_DM_NAIVE(f) PASTE(runner_doself_subset_dm_naive, f)
#define DOSELF_SUBSET_DM_NAIVE _DOSELF_SUBSET_DM_NAIVE(FUNCTION)

#define _IACT_NONSYM_DM(f) PASTE(runner_iact_nonsym_dm, f)
#define IACT_NONSYM_DM _IACT_NONSYM_DM(FUNCTION)

#define _IACT_DM(f) PASTE(runner_iact_dm, f)
#define IACT_DM _IACT_DM(FUNCTION)

#define _TIMER_DOSELF_DM(f) PASTE(timer_doself_dm, f)
#define TIMER_DOSELF_DM _TIMER_DOSELF_DM(FUNCTION)

#define _TIMER_DOPAIR_DM(f) PASTE(timer_dopair_dm, f)
#define TIMER_DOPAIR_DM _TIMER_DOPAIR_DM(FUNCTION)

#define _TIMER_DOSELF_SUBSET_DM(f) PASTE(timer_doself_subset_dm, f)
#define TIMER_DOSELF_SUBSET_DM _TIMER_DOSELF_SUBSET_DM(FUNCTION)

#define _TIMER_DOPAIR_SUBSET_DM(f) PASTE(timer_dopair_subset_dm, f)
#define TIMER_DOPAIR_SUBSET_DM _TIMER_DOPAIR_SUBSET_DM(FUNCTION)

