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
#ifndef SWIFT_RT_ACTIVE_H
#define SWIFT_RT_ACTIVE_H

/* Local includes. */
#include "active.h"

/**
 * @file rt_active.h
 * @brief This file contains the checks for whether radiative transfer
 * (injection) tasks should be activated for given cells, parts, or sparts.
 * The functions differ from the checks in `src/active.h` in that not only
 * time-step data is being checked, but more properties as well. So for
 * consistency, they get their own file. Finally, the functions are gathered
 * in one RT related file to concentrate all #ifdef macro shenanigans in a
 * single place as far as possible.
 */

/**
 * @brief Pick whether we need to check whether an spart is active depending
 * on which version of injection we are using. This is done here to minimize
 * #ifdef macros throughout this file.
 *
 * Returns 1 if spart is active.
 *
 * @param sp star particle
 * @param e the engine
 */
__attribute__((always_inline)) INLINE static int rt_is_spart_active_in_loop(
    struct spart *restrict sp, const struct engine *e) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  return 1; /* ignore stellar activity when gas pulls radiation from stars */
#else
  return spart_is_active(sp, e);
#endif
}

/**
 * @brief Pick whether we need to check whether a part is active depending
 * on which version of injection we are using. This is done here to minimize
 * #ifdef macros throughout this file.
 *
 * Returns 1 if spart is active.
 *
 * @param p the part
 * @param e the engine
 */
__attribute__((always_inline)) INLINE static int rt_is_part_active_in_loop(
    struct part *restrict p, const struct engine *e) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  return part_is_active(p, e);
#else
  return 1; /* ignore hydro activity when stars push radiation onto gas */
#endif
}

/**
 * @brief Does a cell contain particles that should do RT this step?
 * This function is for a self-type interaction, where we need a cell
 * to have active hydro particles and star particles in any state.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int rt_should_iact_cell(
    const struct cell *c, const struct engine *e) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  return ((cell_is_active_hydro(c, e) && (c->hydro.count > 0)) &&
          (c->stars.count > 0));
#else
  return ((cell_is_active_stars(c, e) && (c->stars.count > 0)) &&
          (c->hydro.count > 0));
#endif
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
__attribute__((always_inline)) INLINE static int rt_should_iact_cell_pair(
    const struct cell *ci, const struct cell *cj, const struct engine *e) {

#ifdef RT_HYDRO_CONTROLLED_INJECTION
  return (cell_is_active_hydro(cj, e) && (cj->hydro.count > 0) &&
          (ci->stars.count > 0));
#else
  return (cell_is_active_stars(ci, e) && (ci->stars.count > 0) &&
          (cj->hydro.count > 0));
#endif
}

/**
 * @brief Do we need to unskip this cell's RT (injection) related tasks?
 * For unskipping (recursively), don't check about the star's count here:
 * Pair-type interactions don't require stars in every cell. This way, we
 * can include the check for star count > 0 there on a pair by pair basis.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell needs to activate tasks
 */
__attribute__((always_inline)) INLINE static int rt_should_do_unskip_cell(
    const struct cell *c, const struct engine *e) {
  /* whether it's hydro controlled or not, we need to check for hydro
   * activity at the top level. We also need to check for star activity
   * so we can activate the rt_in implicit tasks to catch dependencies
   * before the injection and not be overwritten by work in star density
   * ghosts. */
  return ((cell_is_active_hydro(c, e) && (c->hydro.count > 0)) ||
          cell_is_active_stars(c, e));
}

#endif /* SWIFT_RT_ACTIVE_H */
