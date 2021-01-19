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
#ifndef SWIFT_CELL_SINKS_H
#define SWIFT_CELL_SINKS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "timeline.h"

/**
 * @brief Sinks-related cell variables.
 */
struct cell_sinks {

  /* If we are not using sinks, compact as much of the unecessary variables
     into an anonymous union to save memory in the cell structure. */
#ifdef SINK_NONE
  union {
#endif

    /*! Pointer to the #sink data. */
    struct sink *parts;

    /*! Linked list of the tasks computing this cell's sink formation checks. */
    struct link *compute_formation;

    /*! The drift task for sinks */
    struct task *drift;

    /*! Implicit tasks marking the entry of the sink block of tasks */
    struct task *sink_in;

    /*! Implicit tasks marking the exit of the sink block of tasks */
    struct task *sink_out;

    /*! Last (integer) time the cell's sink were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Maximum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_max;

    /*! Spin lock for sink formation use. */
    swift_lock_type sink_formation_lock;

    /*! Nr of #sink this cell can hold after addition of new one. */
    int count_total;

    /*! Max cut off radius of active particles in this cell. */
    float r_cut_max_active;

    /*! Values of r_cut_max before the drifts, used for sub-cell tasks. */
    float r_cut_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

#ifdef SINK_NONE
  };
#endif

  /*! Minimum end of (integer) time step in this cell for sink tasks. */
  integertime_t ti_end_min;

  /*! Maximum beginning of (integer) time step in this cell for sink
   * tasks. */
  integertime_t ti_beg_max;

  /*! Spin lock for various uses (#sink case). */
  swift_lock_type lock;

  /*! Max cut off radius in this cell. */
  float r_cut_max;

  /*! Number of #sink updated in this cell. */
  int updated;

  /*! Is the #sink data of this cell being used in a sub-cell? */
  int hold;

  /*! Nr of #sink in this cell. */
  int count;
};

#endif /* SWIFT_CELL_SINKS_H */
