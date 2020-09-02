/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

#include "logger_gravity.h"

/* Define the size of all the fields. */
#define member_size(type, member) sizeof(((type *)0)->member)

const int gravity_logger_field_size[gravity_logger_field_count] = {
    member_size(struct gpart, x),                 // coordinates
    member_size(struct gpart, v_full),            // velocities
    member_size(struct gpart, a_grav),            // accelerations
    member_size(struct gpart, mass),              // massses
    member_size(struct gpart, id_or_neg_offset),  // IDs
};

int gravity_logger_local_to_global[gravity_logger_field_count];
