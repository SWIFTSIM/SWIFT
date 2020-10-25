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
#ifndef SWIFT_CELL_HYDRO_H
#define SWIFT_CELL_HYDRO_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "timeline.h"

/**
 * @brief Hydro-related cell variables.
 */
struct cell_hydro {

  /*! Pointer to the #part data. */
  struct part *parts;

  /*! Pointer to the #xpart data. */
  struct xpart *xparts;

  /*! Pointer for the sorted indices. */
  struct sort_entry *sort;

  /*! Super cell, i.e. the highest-level parent cell that has a hydro
   * pair/self tasks */
  struct cell *super;

  /*! The task computing this cell's sorts. */
  struct task *sorts;

  /*! The drift task for parts */
  struct task *drift;

  /*! Linked list of the tasks computing this cell's hydro density. */
  struct link *density;

  /* Linked list of the tasks computing this cell's hydro gradients. */
  struct link *gradient;

  /*! Linked list of the tasks computing this cell's hydro forces. */
  struct link *force;

  /*! Linked list of the tasks computing this cell's limiter. */
  struct link *limiter;

  /*! Dependency implicit task for the ghost  (in->ghost->out)*/
  struct task *ghost_in;

  /*! Dependency implicit task for the ghost  (in->ghost->out)*/
  struct task *ghost_out;

  /*! The ghost task itself */
  struct task *ghost;

  /*! The extra ghost task for complex hydro schemes */
  struct task *extra_ghost;

  /*! The task to end the force calculation */
  struct task *end_force;

  /*! Dependency implicit task for cooling (in->cooling->out) */
  struct task *cooling_in;

  /*! Dependency implicit task for cooling (in->cooling->out) */
  struct task *cooling_out;

  /*! Task for cooling */
  struct task *cooling;

  /*! Task for star formation */
  struct task *star_formation;

  /*! Task for sink formation */
  struct task *sink_formation;

  /*! Task for sorting the stars again after a SF event */
  struct task *stars_resort;

  /*! Radiative transfer ghost in task */
  struct task *rt_in;

  /*! Radiative transfer ghost out task */
  struct task *rt_out;

  /*! Radiative transfer ghost1 task (finishes up injection) */
  struct task *rt_ghost1;

  /*! Task for self/pair injection step of radiative transfer */
  struct link *rt_inject;

  /*! Last (integer) time the cell's part were drifted forward in time. */
  integertime_t ti_old_part;

  /*! Minimum end of (integer) time step in this cell for hydro tasks. */
  integertime_t ti_end_min;

  /*! Maximum end of (integer) time step in this cell for hydro tasks. */
  integertime_t ti_end_max;

  /*! Maximum beginning of (integer) time step in this cell for hydro tasks.
   */
  integertime_t ti_beg_max;

  /*! Spin lock for various uses (#part case). */
  swift_lock_type lock;

  /*! Max smoothing length in this cell. */
  float h_max;

  /*! Max smoothing length of active particles in this cell. */
  float h_max_active;

  /*! Maximum part movement in this cell since last construction. */
  float dx_max_part;

  /*! Maximum particle movement in this cell since the last sort. */
  float dx_max_sort;

  /*! Values of h_max before the drifts, used for sub-cell tasks. */
  float h_max_old;

  /*! Values of dx_max before the drifts, used for sub-cell tasks. */
  float dx_max_part_old;

  /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
  float dx_max_sort_old;

  /*! Nr of #part in this cell. */
  int count;

  /*! Nr of #part this cell can hold after addition of new #part. */
  int count_total;

  /*! Number of #part updated in this cell. */
  int updated;

  /*! Is the #part data of this cell being used in a sub-cell? */
  int hold;

  /*! Bit mask of sort directions that will be needed in the next timestep. */
  uint16_t requires_sorts;

  /*! Bit mask of sorts that need to be computed for this cell. */
  uint16_t do_sort;

  /*! Bit-mask indicating the sorted directions */
  uint16_t sorted;

  /*! Bit-mask indicating the sorted directions */
  uint16_t sort_allocated;

#ifdef SWIFT_DEBUG_CHECKS

  /*! Last (integer) time the cell's sort arrays were updated. */
  integertime_t ti_sort;

#endif
};

#endif /* SWIFT_CELL_HYDRO_H */
