/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_DO_CELLS_H
#define SWIFT_RT_DO_CELLS_H

/* Local includes. */
#include "active.h"

/**
 * @file rt_do_cells.h
 * @brief This file contains the checks for whether radiative transfer
 * (injection) tasks should be activated for given cells.
 * The functions differ from the checks in `src/active.h` in that not only
 * time-step data is being checked, but more properties as well. So for
 * consistency, they get their own file.
 */

/**
 * @brief Does a cell contain particles that should do RT this step?
 * This function is for a self-type interaction, where we need a cell
 * to have active hydro particles and star particles in any state.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int rt_should_do_cell(
    const struct cell *c, const struct engine *e) {

  return ((cell_is_active_hydro(c, e) && (c->hydro.count > 0)) &&
          (c->stars.count > 0));
}

/**
 * @brief Does a cell contain particles that should do RT this step?
 * This function is for a pair-type interaction, where we take stars from
 * cell ci and hydro particles from cell cj.
 *
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int rt_should_do_cell_pair(
    const struct cell *ci, const struct cell *cj, const struct engine *e) {

  return (cell_is_active_hydro(cj, e) && (cj->hydro.count > 0) &&
          (ci->stars.count > 0));
}

/**
 * @brief Do we need to unskip this cell's RT (injection) related tasks?
 * For unskipping (recursively), don't check about the star's count here.
 * Pair tasks might need to be checked for as well first. This way, we can
 * be more restrictive for the check for pair type tasks and include the
 * check for star count > 0 there on a pair by pair basis.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell needs to activate tasks
 */
__attribute__((always_inline)) INLINE static int rt_should_do_unskip_cell(
    const struct cell *c, const struct engine *e) {

  return (cell_is_active_hydro(c, e) && (c->hydro.count > 0));
}

#endif /* SWIFT_RT_DO_CELLS_H */
