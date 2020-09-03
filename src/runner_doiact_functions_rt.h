/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#include "runner_doiact_rt.h"


/**
 * @brief TODO
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_RT(struct runner *r, struct cell *c, int timer){
  TIMER_TIC; 
  message("message from DOSELF1_RT");
  if (timer) TIMER_TOC(TIMER_DOSELF_RT);
}

/**
 * @brief TODO
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_SYM_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer){
  TIMER_TIC; 
  message("message from DOPAIR1_RT");
  if (timer) TIMER_TOC(TIMER_DOPAIR_RT);
}

/**
 * @brief Determine which version of DOSELF1_RT needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_BRANCH_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

/**
 * @brief Determine which version of DOPAIR1_RT needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_BRANCH_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
  DOPAIR1_SYM_RT(r, ci, cj, timer);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

/**
 * @brief Compute grouped sub-cell interactions for pair tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer) {
  DOPAIR1_SYM_RT(r, ci, cj, timer);
}
