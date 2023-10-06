/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_EQUATION_OF_STATE_H
#define SWIFT_EQUATION_OF_STATE_H

/**
 * @file src/equation_of_state.h
 * @brief Defines the equation of state of the gas we simulate in the form of
 * relations between thermodynamic quantities. These are later used internally
 * by all hydro schemes
 */

/* Config parameters. */
#include <config.h>

/* Import the right functions */
#if defined(EOS_IDEAL_GAS)
#include "./equation_of_state/ideal_gas/equation_of_state.h"
#elif defined(EOS_ISOTHERMAL_GAS)
#include "./equation_of_state/isothermal/equation_of_state.h"
#elif defined(EOS_PLANETARY)
#include "./equation_of_state/planetary/equation_of_state.h"
#elif defined(EOS_BAROTROPIC_GAS)
#include "./equation_of_state/barotropic/equation_of_state.h"
#else
#error "Invalid choice of equation of state"
#endif

#endif /* SWIFT_EQUATION_OF_STATE_H */
