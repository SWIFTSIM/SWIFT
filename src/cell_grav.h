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
#include <config.h>

/* Local includes. */
#include "gravity.h"
#include "lock.h"
#include "timeline.h"

/**
 * @brief Gravity-related cell variables.
 */

struct cell_grav {

  union {

    /*! Pointer to the #gpart data. */
    struct gpart *parts;

    /*! or #gpart_foreign data. */
    struct gpart_foreign *parts_foreign;

    /*! or #gpart_fof_foreign data. */
    struct gpart_fof_foreign *parts_fof_foreign;
  };

  union {

    /*! Pointer to the #gpart data at rebuild time. */
    struct gpart *parts_rebuild;

    /*! Pointer to the #gpart_foreign data at rebuild time. */
    struct gpart_foreign *parts_foreign_rebuild;

    /*! Pointer to the #gpart_fof_foreign data at rebuild time. */
    struct gpart_fof_foreign *parts_fof_foreign_rebuild;
  };

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

#ifdef SWIFT_DEBUG_CHECKS
  /*! ti_old_part on entry to the last cell_drift_gpart() call on this cell,
   * captured before any drift happens. Diagnostic for the undrifted-gpart
   * pack check: if this already equalled ti_current, the drift's own guard
   * skipped the cell entirely. */
  integertime_t ti_old_part_on_entry;

  /*! Number of gparts actually visited by the last direct (non-recursive)
   * drift of this cell, or -1 if the direct drift branch wasn't taken
   * (cell was split, or the drift guard skipped it). */
  int count_drifted;

  /*! Whether force was already true when the last cell_drift_gpart() call
   * reached this cell (i.e. cell_flag_do_grav_drift was set directly on
   * it, not just inherited/cell_flag_do_grav_sub_drift on an ancestor). If
   * 0 for a split cell, only the sub-path with its own flags got drifted;
   * siblings outside that path were skipped despite ti_old_part being
   * stamped as current. */
  int drift_force_on_entry;
#endif

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

#ifdef WITH_MPI
  /*! #count at the last gpart data pack, before particle creation can change it.
   * Sent to the foreign copy as #count_valid via task_subtype_grav_counts. */
  int count_packed;
#endif

  /*! For a foreign cell, the sender's #count_packed at the last data
   * delivery. Lags #count on spawn steps; use this, not #count, to bound
   * any read of #parts_foreign. */
  int count_valid;

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
