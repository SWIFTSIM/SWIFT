/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_CELL_RT_H
#define SWIFT_CELL_RT_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "timeline.h"

/**
 * @brief Radiative transfer related cell variables.
 */
struct cell_rt {

  /* If we are not using RT, compact as much of the unecessary variables
     into an anonymous union to save memory in the cell structure. */
#ifdef RT_NONE
  union {
#endif

    /*! Radiative transfer ghost in task */
    struct task *rt_in;

    /*! Radiative transfer ghost1 task (finishes up injection) */
    struct task *rt_ghost1;

    /*! Task for self/pair gradient step of radiative transfer */
    struct link *rt_gradient;

    /*! Radiative transfer ghost2 task */
    struct task *rt_ghost2;

    /*! Task for self/pair transport step of radiative transfer */
    struct link *rt_transport;

    /*! Radiative transfer transport out task */
    struct task *rt_transport_out;

    /*! Radiative transfer thermochemistry task */
    struct task *rt_tchem;

    /*! Radiative transfer cell time advancement task */
    struct task *rt_advance_cell_time;

    /*! Sort a cell after a recv rt gradients */
    struct task *rt_sorts;

    /*! Collect the cell times from the super to the top level */
    struct task *rt_collect_times;

    /*! Radiative transfer ghost out task */
    struct task *rt_out;

    /*! Bit mask of sorts that need to be computed for this cell.
     * Needed to be able to skip sorting undrifted cells. */
    uint16_t do_sort;

#ifdef RT_NONE
  };
#endif

#ifdef SWIFT_RT_DEBUG_CHECKS
  /*! has rt_advance_cell_time run on this cell? */
  int advanced_time;
#endif

  /*! Minimum end of (integer) time step in this cell for RT tasks. */
  integertime_t ti_rt_end_min;

  /*! Maximum beginning of (integer) time step in this cell for RT tasks. */
  integertime_t ti_rt_beg_max;

  /*! Minimum (integer) time step size in this cell for RT tasks. */
  integertime_t ti_rt_min_step_size;

  /*! Number of #part updated for RT in this cell */
  int updated;
};

#endif /* SWIFT_CELL_RT_H */
