/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_STAR_FORMATION_STRUCT_H
#define SWIFT_NONE_STAR_FORMATION_STRUCT_H

/* Do we need unique IDs (only useful when spawning
   new particles, conversion gas->stars does not need unique IDs) */
#define star_formation_need_unique_id 0

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {};

/**
 * @brief Star-formation-related properties stored in the particle data.
 */
struct star_formation_part_data {};

/**
 * @brief Star-formation-related properties stored in the star particle
 * data.
 */
struct star_formation_spart_data {};

#endif /* SWIFT_NONE_STAR_FORMATION_STRUCT_H */
