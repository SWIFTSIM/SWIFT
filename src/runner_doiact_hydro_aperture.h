/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Before including this file, define FUNCTION (e.g. "prep_sink_formation"),
   which expands to the full set of interaction function names used by the
   fixed-aperture gas-gas neighbour loop:

     runner_dopair1_branch_hydro_aperture_FUNCTION
     runner_dopair1_hydro_aperture_FUNCTION
     runner_dopair1_naive_hydro_aperture_FUNCTION
     runner_doself1_branch_hydro_aperture_FUNCTION
     runner_doself1_hydro_aperture_FUNCTION
     runner_doself1_naive_hydro_aperture_FUNCTION
     runner_dosub_self1_hydro_aperture_FUNCTION
     runner_dosub_pair1_hydro_aperture_FUNCTION

   Only the non-symmetric _1_ variants exist.  The loop uses a pure gather
   pattern (each active particle accumulates from neighbours), so no symmetric
   _2_ variants or subset-reiteration functions are needed. */

#define PASTE(x, y) x##_##y

#define _DOPAIR1_BRANCH_HYDRO_APERTURE(f) \
  PASTE(runner_dopair1_branch_hydro_aperture, f)
#define DOPAIR1_BRANCH_HYDRO_APERTURE _DOPAIR1_BRANCH_HYDRO_APERTURE(FUNCTION)

#define _DOPAIR1_HYDRO_APERTURE(f) PASTE(runner_dopair1_hydro_aperture, f)
#define DOPAIR1_HYDRO_APERTURE _DOPAIR1_HYDRO_APERTURE(FUNCTION)

#define _DOPAIR1_NAIVE_HYDRO_APERTURE(f) \
  PASTE(runner_dopair1_naive_hydro_aperture, f)
#define DOPAIR1_NAIVE_HYDRO_APERTURE _DOPAIR1_NAIVE_HYDRO_APERTURE(FUNCTION)

#define _DOSELF1_BRANCH_HYDRO_APERTURE(f) \
  PASTE(runner_doself1_branch_hydro_aperture, f)
#define DOSELF1_BRANCH_HYDRO_APERTURE _DOSELF1_BRANCH_HYDRO_APERTURE(FUNCTION)

#define _DOSELF1_HYDRO_APERTURE(f) PASTE(runner_doself1_hydro_aperture, f)
#define DOSELF1_HYDRO_APERTURE _DOSELF1_HYDRO_APERTURE(FUNCTION)

#define _DOSELF1_NAIVE_HYDRO_APERTURE(f) \
  PASTE(runner_doself1_naive_hydro_aperture, f)
#define DOSELF1_NAIVE_HYDRO_APERTURE _DOSELF1_NAIVE_HYDRO_APERTURE(FUNCTION)

#define _DOSUB_SELF1_HYDRO_APERTURE(f) \
  PASTE(runner_dosub_self1_hydro_aperture, f)
#define DOSUB_SELF1_HYDRO_APERTURE _DOSUB_SELF1_HYDRO_APERTURE(FUNCTION)

#define _DOSUB_PAIR1_HYDRO_APERTURE(f) \
  PASTE(runner_dosub_pair1_hydro_aperture, f)
#define DOSUB_PAIR1_HYDRO_APERTURE _DOSUB_PAIR1_HYDRO_APERTURE(FUNCTION)

#define _IACT_NONSYM_HYDRO_APERTURE(f) \
  PASTE(runner_iact_nonsym_hydro_aperture, f)
#define IACT_NONSYM_HYDRO_APERTURE _IACT_NONSYM_HYDRO_APERTURE(FUNCTION)

#define _IACT_HYDRO_APERTURE(f) PASTE(runner_iact_hydro_aperture, f)
#define IACT_HYDRO_APERTURE _IACT_HYDRO_APERTURE(FUNCTION)

#define _TIMER_DOSELF_HYDRO_APERTURE(f) \
  PASTE(timer_doself_hydro_aperture, f)
#define TIMER_DOSELF_HYDRO_APERTURE _TIMER_DOSELF_HYDRO_APERTURE(FUNCTION)

#define _TIMER_DOPAIR_HYDRO_APERTURE(f) \
  PASTE(timer_dopair_hydro_aperture, f)
#define TIMER_DOPAIR_HYDRO_APERTURE _TIMER_DOPAIR_HYDRO_APERTURE(FUNCTION)

#define _TIMER_DOSUB_SELF_HYDRO_APERTURE(f) \
  PASTE(timer_dosub_self_hydro_aperture, f)
#define TIMER_DOSUB_SELF_HYDRO_APERTURE \
  _TIMER_DOSUB_SELF_HYDRO_APERTURE(FUNCTION)

#define _TIMER_DOSUB_PAIR_HYDRO_APERTURE(f) \
  PASTE(timer_dosub_pair_hydro_aperture, f)
#define TIMER_DOSUB_PAIR_HYDRO_APERTURE \
  _TIMER_DOSUB_PAIR_HYDRO_APERTURE(FUNCTION)

/* Active/drifted status macros for the fixed-aperture gas-gas loop. */
#define PART_IS_ACTIVE part_is_active
#define CELL_IS_ACTIVE cell_is_active_hydro
#define CELL_ARE_PART_DRIFTED cell_are_part_drifted
#define DO_DRIFT_DEBUG_CHECKS 1
#define GET_MU0() \
  {               \
  }

void DOSELF1_BRANCH_HYDRO_APERTURE(struct runner *r, const struct cell *c,
                                   const float r_cut);

void DOPAIR1_BRANCH_HYDRO_APERTURE(struct runner *r, struct cell *ci,
                                   struct cell *cj, const float r_cut);

void DOSUB_SELF1_HYDRO_APERTURE(struct runner *r, struct cell *c,
                                const float r_cut, const int gettimer);

void DOSUB_PAIR1_HYDRO_APERTURE(struct runner *r, struct cell *ci,
                                struct cell *cj, const float r_cut,
                                const int gettimer);
