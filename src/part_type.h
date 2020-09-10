/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_PART_TYPES_H
#define SWIFT_PART_TYPES_H

/**
 * @brief The different types of particles a #gpart can link to.
 *
 * Note we use the historical values from Gadget for these fields.
 */
enum part_type {
  swift_type_gas = 0,
  swift_type_dark_matter = 1,
  swift_type_dark_matter_background = 2,
  swift_type_sink = 3,
  swift_type_stars = 4,
  swift_type_black_hole = 5,
  swift_type_count
} __attribute__((packed));

extern const char* part_type_names[];

#endif /* SWIFT_PART_TYPES_H */
