/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

/* This object's header. */
#include "rt_parameters.h"

/* Some parameters related to the Radiative Transfer
 * (temporary ugly solution as a global variable) */

#ifdef __APPLE__
/*
 * The clang compiler and linker on OSX incorrectly optimize
 * out the global object before the final linking stage, which
 * leads to a compilation error.
 * The fake initialisation below forces the compiler to keep the
 * instance and pass it to the linker stage.
 */
#if defined(RT_GEAR)
struct rt_parameters rt_params = {.reduced_speed_of_light = 1.f,
                                  .reduced_speed_of_light_inverse = 1.f};
#else
struct rt_parameters rt_params;
#endif

#else  /* i.e. not __APPLE__ */
struct rt_parameters rt_params;
#endif /* __APPLE__ */
