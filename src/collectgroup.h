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
  size_t updates, g_updates, s_updates;

  /* Times for the time-step */
  integertime_t ti_end_min, ti_end_max, ti_beg_max;

  /* Force the engine to rebuild? */
  int forcerebuild;
};

void collectgroup_init();
void collectgroup1_apply(struct collectgroup1 *grp1, struct engine *e);
void collectgroup1_init(struct collectgroup1 *grp1, size_t updates,
                        size_t g_updates, size_t s_updates,
                        integertime_t ti_end_min,
                        integertime_t ti_end_max,
                        integertime_t ti_beg_max,
                        int forcerebuild);
void collectgroup1_reduce(struct collectgroup1 *grp1);

#endif /* SWIFT_COLLECTGROUP_H */
