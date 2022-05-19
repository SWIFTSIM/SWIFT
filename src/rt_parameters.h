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
#ifndef SWIFT_RT_PARAMETERS_H
#define SWIFT_RT_PARAMETERS_H

/**
 * @file src/rt_parameters.h
 * @brief Defines radiative transfer parameters to make them
 * available during hydro task like interactions without explicitly
 * passing them as function parameters
 * */

/* Config parameters. */
#include "../config.h"

/* Import the right functions */
#if defined(RT_DEBUG)
#include "./rt/debug/rt_parameters.h"
#elif defined(RT_GEAR)
#include "./rt/GEAR/rt_parameters.h"
#elif defined(RT_SPHM1RT)
#include "./rt/SPHM1RT/rt_parameters.h"
#elif defined(RT_NONE)
#include "./rt/none/rt_parameters.h"
#else
#error "Invalid choice of radiative transfer scheme"
#endif

#endif /* SWIFT_RT_PARAMETERS_H */
