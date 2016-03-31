/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_PHYSICAL_CONSTANTS_H
#define SWIFT_PHYSICAL_CONSTANTS_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "units.h"

/* physical constants in in defined programme units */
struct phys_const {
  double newton_gravity;
};

/**
 * @brief Converts physical constants to the internal unit system
 *
 * @param us The current internal system of units.
 * @param internal_const The physical constants to initialize.
 */
void initPhysicalConstants(struct UnitSystem* us,
                           struct phys_const* internal_const);

#endif /* SWIFT_PHYSICAL_CONSTANTS_H */
