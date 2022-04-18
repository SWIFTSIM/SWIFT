/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CELL_GRAV_H
#define SWIFT_CELL_GRAV_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "timeline.h"

/**
 * @brief Gravity-related cell variables.
 */

struct cell_grav {

  /*! Pointer to the #gpart data. */
  struct gpart *parts;

  /*! Pointer to the #spart data at rebuild time. */
  struct gpart *parts_rebuild;

  /*! This cell's multipole. */
  struct gravity_tensors *multipole;

  /*! Super cell, i.e. the highest-level parent cell that has a grav pair/self
   * tasks */
  struct cell *super;

  /*! The drift task for gparts */
  struct task *drift;

  /*! Implicit task (going up- and down the tree) for the #gpart drifts */
  struct task *drift_out;

  /*! Linked list of the tasks computing this cell's gravity forces. */
  struct link *grav;

  /*! Linked list of the tasks computing this cell's gravity M-M forces. */
  struct link *mm;

  /*! The multipole initialistation task */
  struct task *init;

  /*! Implicit task for the gravity initialisation */
  struct task *init_out;

  /*! Task computing long range non-periodic gravity interactions */
  struct task *long_range;

  /*! Implicit task for the down propagation */
  struct task *down_in;

  /*! Task propagating the multipole to the particles */
  struct task *down;

  /*! The task to end the force calculation */
  struct task *end_force;

  /*! Task for weighting neutrino particles */
  struct task *neutrino_weight;

  /*! Minimum end of (integer) time step in this cell for gravity tasks. */
  integertime_t ti_end_min;

  /*! Maximum beginning of (integer) time step in this cell for gravity tasks.
   */
  integertime_t ti_beg_max;

  /*! Last (integer) time the cell's gpart were drifted forward in time. */
  integertime_t ti_old_part;

  /*! Last (integer) time the cell's multipole was drifted forward in time. */
  integertime_t ti_old_multipole;

  /*! Spin lock for various uses (#gpart case). */
  swift_lock_type plock;

  /*! Spin lock for various uses (#multipole case). */
  swift_lock_type mlock;

  /*! Spin lock for star formation use. */
  swift_lock_type star_formation_lock;

  /*! Nr of #gpart in this cell. */
  int count;

  /*! Nr of #gpart this cell can hold after addition of new #gpart. */
  int count_total;

  /*! Number of #gpart updated in this cell. */
  int updated;

  /*! Is the #gpart data of this cell being used in a sub-cell? */
  int phold;

  /*! Is the #multipole data of this cell being used in a sub-cell? */
  int mhold;

  /*! Number of M-M tasks that are associated with this cell. */
  short int nr_mm_tasks;
};

#endif /* SWIFT_CELL_GRAV_H */
