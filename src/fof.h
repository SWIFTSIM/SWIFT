/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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

#ifndef SWIFT_FOF_H
#define SWIFT_FOF_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "cell.h"
#include "space.h"

/* Function prototypes. */
void fof_search_serial(struct space *s);
void fof_search_cell(struct space *s, struct cell *c);
void fof_search_pair_cells(struct space *s, struct cell *ci, struct cell *cj);
void fof_search_tree_serial(struct space *s);
void fof_dump_group_data(char *out_file, const size_t nr_gparts, int *group_id, int *num_in_groups);

#endif /* SWIFT_FOF_H */
