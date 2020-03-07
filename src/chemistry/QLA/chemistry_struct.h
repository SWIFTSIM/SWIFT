/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_STRUCT_QLA_H
#define SWIFT_CHEMISTRY_STRUCT_QLA_H

/**
 * @file src/chemistry/none/chemistry_struct.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/**
 * @brief The individual elements traced in the model.
 */
enum chemistry_element { chemistry_element_count = 0 };

/**
 * @brief Global chemical abundance information.
 *
 * Nothing here.
 */
struct chemistry_global_data {};

/**
 * @brief Chemistry properties carried by the #part.
 *
 * Nothing here.
 */
struct chemistry_part_data {};

/**
 * @brief Chemistry properties carried by the #bpart.
 *
 * Nothing here.
 */
struct chemistry_bpart_data {};

#endif /* SWIFT_CHEMISTRY_STRUCT_QLA_H */
