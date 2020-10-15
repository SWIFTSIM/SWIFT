/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 James Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_FOF6D_H
#define SWIFT_FOF6D_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "align.h"
#include "parser.h"
#include "fof.h"

/* Avoid cyclic inclusions */
struct gpart;
struct space;
struct engine;
struct unit_system;
struct phys_const;
struct cosmology;

/* Function prototypes. */
void fof6d_calc_vel_disp(struct fof_props *props, struct space *s, const size_t num_parts_in_groups); 
void fof6d_split_groups(struct fof_props *props, struct space *s, const size_t num_parts_in_groups, const double *v_disp, const size_t *part_index); 

#endif /* SWIFT_FOF6D_H */
