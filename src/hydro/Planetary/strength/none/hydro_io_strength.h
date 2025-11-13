/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#ifndef SWIFT_PLANETARY_STRENGTH_IO_NONE_H
#define SWIFT_PLANETARY_STRENGTH_IO_NONE_H

/**
 * @file Planetary/strength/none/hydro_io_strength.h
 * @brief Planetary implementation of SPH with default material strength
 */

#include "adiabatic_index.h"
#include "hydro.h"
#include "hydro_parameters.h"
#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle strength fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void hydro_read_particles_strength(struct part* parts,
                                                 struct io_props* list,
                                                 int* num_fields) {
    
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
}

/**
 * @brief Specifies which particle strength fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void hydro_write_particles_strength(const struct part* parts,
                                         const struct xpart* xparts,
                                         struct io_props* list,
                                         int* num_fields) {}

#endif /* SWIFT_PLANETARY_STRENGTH_IO_NONE_H */