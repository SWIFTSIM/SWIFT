/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_COLLECTGROUP_H
#define SWIFT_COLLECTGROUP_H

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stddef.h>

/* Local headers. */
#include "timeline.h"

/* Forward declaration of engine struct (to avoid cyclic include). */
struct engine;

/* A collection of global quantities that can be processed at the same time. */
struct collectgroup1 {

  /* Number of particles updated */
  long long updated, g_updated, s_updated;

  /* Number of particles inhibited */
  long long inhibited, g_inhibited, s_inhibited;

  /* Times for the time-step */
  integertime_t ti_hydro_end_min, ti_hydro_end_max, ti_hydro_beg_max;
  integertime_t ti_gravity_end_min, ti_gravity_end_max, ti_gravity_beg_max;

  /* Force the engine to rebuild? */
  int forcerebuild;

  /* Totals of cells and tasks. */
  long long total_nr_cells;
  long long total_nr_tasks;

  /* Maximum value of actual tasks per cell across all ranks. */
  float tasks_per_cell_max;
};

void collectgroup_init(void);
void collectgroup1_apply(struct collectgroup1 *grp1, struct engine *e);
void collectgroup1_init(
    struct collectgroup1 *grp1, size_t updated, size_t g_updated,
    size_t s_updated, size_t inhibited, size_t g_inhibited, size_t s_inhibited,
    integertime_t ti_hydro_end_min, integertime_t ti_hydro_end_max,
    integertime_t ti_hydro_beg_max, integertime_t ti_gravity_end_min,
    integertime_t ti_gravity_end_max, integertime_t ti_gravity_beg_max,
    int forcerebuild, long long total_nr_cells, long long total_nr_tasks,
    float tasks_per_cell);
void collectgroup1_reduce(struct collectgroup1 *grp1);

#endif /* SWIFT_COLLECTGROUP_H */
