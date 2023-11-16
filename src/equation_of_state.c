/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "equation_of_state.h"

/* Equation of state for the physics model
 * (temporary ugly solution as a global variable) */

#ifdef __APPLE__
/*
 * The clang compiler and linker on OSX incorrectly optimize
 * out the eos global object before the final linking stage, which
 * leads to a compilation error.
 * The fake initialisation below forces the compiler to keep the
 * instance and pass it to the linker stage.
 */
#if defined(EOS_PLANETARY)
struct eos_parameters eos = {.Til_iron.rho_0 = -1.f};
#elif defined(EOS_ISOTHERMAL_GAS)
struct eos_parameters eos = {.isothermal_internal_energy = -1.};
#elif defined(EOS_BAROTROPIC_GAS)
struct eos_parameters eos = {.vacuum_sound_speed2 = -1.,
                             .inverse_core_density = -1.};
#else
struct eos_parameters eos;
#endif

#else  /* i.e. not __APPLE__ */
struct eos_parameters eos;
#endif /* __APPLE__ */
