/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H
#define SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H

#include "hydro_properties.h"
#include "inline.h"

/**
 * @brief Properties of the DF model.
 */
struct df_props {};

/**
 * @brief Initialize the global properties of the df scheme.
 *
 * Nothing to do here for the no feedback model.
 *
 * @param dfp The #df_props.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param hydro_props The already read-in properties of the hydro scheme.
 * @param cosmo The cosmological model.
 */
INLINE static void df_props_init(struct df_props *dfp,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct hydro_props *hydro_props,
                                       const struct cosmology *cosmo) {}

#endif /* SWIFT_BASIC_DYNAMICAL_FRICTION_PROPERTIES_H */
