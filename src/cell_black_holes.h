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
#ifndef SWIFT_CELL_BLACK_HOLES_H
#define SWIFT_CELL_BLACK_HOLES_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "lock.h"
#include "timeline.h"

/**
 * @brief BHs-related cell variables.
 */
struct cell_black_holes {

  /* If we are not using BHs, compact as much of the unecessary variables
     into an anonymous union to save memory in the cell structure. */
#ifdef BLACK_HOLES_NONE
  union {
#endif

    /*! Pointer to the #bpart data. */
    struct bpart *parts;

    /*! The drift task for bparts */
    struct task *drift;

    /*! Implicit tasks marking the entry of the BH physics block of tasks */
    struct task *black_holes_in;

    /*! Implicit tasks marking the exit of the BH physics block of tasks */
    struct task *black_holes_out;

    /*! The black hole ghost task itself */
    struct task *density_ghost;

    /*! The ghost tasks related to BH swallowing */
    struct task *swallow_ghost_1;
    struct task *swallow_ghost_2;
    struct task *swallow_ghost_3;

    /*! Linked list of the tasks computing this cell's BH density. */
    struct link *density;

    /*! Linked list of the tasks computing this cell's BH swallowing and
     * merging. */
    struct link *swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_gas_swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_bh_swallow;

    /*! Linked list of the tasks computing this cell's BH feedback. */
    struct link *feedback;

    /*! Last (integer) time the cell's bpart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Nr of #bpart this cell can hold after addition of new #bpart. */
    int count_total;

    /*! Max smoothing length of active particles in this cell. */
    float h_max_active;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

#ifdef BLACK_HOLES_NONE
  };
#endif

  /*! Maximum end of (integer) time step in this cell for black tasks. */
  integertime_t ti_end_min;

  /*! Maximum beginning of (integer) time step in this cell for black hole
   * tasks. */
  integertime_t ti_beg_max;

  /*! Spin lock for various uses (#bpart case). */
  swift_lock_type lock;

  /*! Max smoothing length in this cell. */
  float h_max;

  /*! Is the #bpart data of this cell being used in a sub-cell? */
  int hold;

  /*! Number of #bpart updated in this cell. */
  int updated;

  /*! Nr of #bpart in this cell. */
  int count;
};

#endif /* SWIFT_CELL_BLACK_HOLES_H */
