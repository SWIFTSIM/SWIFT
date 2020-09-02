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

#include "logger_stars.h"

/* Define the size of all the fields. */
#define member_size(type, member) sizeof(((type *)0)->member)

const int stars_logger_field_size[stars_logger_field_count] = {
    member_size(struct spart, x),       // coordinates
    member_size(struct spart, v),       // velocities
    member_size(struct gpart, a_grav),  // accelerations -> stored inside gparts
    member_size(struct spart, mass),    // massses
    member_size(struct spart, h),       // Smoothing Length
    member_size(struct spart, id),      // IDs
};

int stars_logger_local_to_global[stars_logger_field_count];
