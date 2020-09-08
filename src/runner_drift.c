/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"

/**
 * @brief Drift all part in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_part(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_part(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_part);
}

/**
 * @brief Drift all gpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_gpart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_gpart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_gpart);
}

/**
 * @brief Drift all spart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_spart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_spart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_spart);
}

/**
 * @brief Drift all bpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_bpart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_bpart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_bpart);
}

/**
 * @brief Drift all dmpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_dmpart(struct runner *r, struct cell *c, int timer) {
    
    TIMER_TIC;
    
    cell_drift_dmpart(c, r->e, 0);
    
    if (timer) TIMER_TOC(timer_drift_dmpart);
}
