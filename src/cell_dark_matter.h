/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Camila Correa (camila.correa@uva.nl)
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
#ifndef SWIFT_CELL_DARK_MATTER_H
#define SWIFT_CELL_DARK_MATTER_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "lock.h"
#include "dark_matter_logger_struct.h"
#include "timeline.h"

/**
 * @brief Dark matter-related cell variables.
 */

struct cell_dark_matter {

/* If we are not using self-interacting dark matter, compact as much of the unnecessary variables
   into an anonymous union to save memory in the cell structure. */
#ifdef SIDM_NONE
union {
#endif

    /*! Pointer to the #dmpart data. */
    struct dmpart *parts;

    /*! Pointer to the #dmpart data at rebuild time. */
    struct dmpart *parts_rebuild;

    /*! Super cell, i.e. the highest-level parent cell that has a grav pair/self
     * tasks */
    struct cell *super;

    /*! Pointer for the sorted indices. */
    struct sort_entry *sort;

    /*! The task computing this cell's sorts. */
    struct task *sorts;

    /*! Linked list of the tasks computing this cell's dm self-interactions. */
    struct link *sidm;

    /*! The drift task for dmparts */
    struct task *drift;

    /*! The dark matter ghost task itself */
    struct task *ghost;

    /*! Linked list of the tasks computing this cell's dark matter density. */
    struct link *density;

    /*! kick due to DM-DM interactions */
    struct task *sidm_kick;

    /*! The task to synchronize the time-step of inactive particles hit by
     * feedback */
    struct task *timestep_sync;

    /*! Last (integer) time the cell's spart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Spin lock for various uses (#dmpart case). */
    swift_lock_type lock;

    /*! Nr of #dmpart in this cell. */
    int count;

    /*! Nr of #dmpart this cell can hold after addition of new #dmpart. */
    int count_total;

    /*! Max smoothing length in this cell. */
    float h_max;

    /*! Max smoothing length of active particles in this cell. */
    float h_max_active;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Maximum particle movement in this cell since the last sort. */
    float dx_max_sort;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

    /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
    float dx_max_sort_old;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sort_allocated;

    /*! Maximum end of (integer) time step in this cell for star tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for star tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for star tasks.
     */
    integertime_t ti_beg_max;

    /*! Number of #spart updated in this cell. */
    int updated;

    /*! Is the #dmpart data of this cell being used in a sub-cell? */
    int hold;

    /*! Bit mask of sort directions that will be needed in the next timestep. */
    uint16_t requires_sorts;

    /*! Bit mask of sorts that need to be computed for this cell. */
    uint16_t do_sort;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sorted;

    /*! SIDM history struct */
    struct sidm_history sh;

#ifdef SWIFT_DEBUG_CHECKS
    /*! Last (integer) time the cell's sort arrays were updated. */
    integertime_t ti_sort;
#endif

#ifdef SIDM_NONE
    };
#endif

};

#endif /* SWIFT_CELL_DARK_MATTER_H */