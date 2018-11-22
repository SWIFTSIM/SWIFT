/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 STFC (Author email aidan.chalk@stfc.ac.uk)
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
#ifndef SWIFT_NUMA_HELPER_H
#define SWIFT_NUMA_HELPER_H

#include <stdint.h>
#include "cell.h"

void swiftnuma_cell_move_hydro_parts_threadpool_map(void *map_data, int num_elements, void *extra_data);
void swiftnuma_cell_move_hydro_parts(struct cell *c, int32_t node, int32_t verbose);

#endif
