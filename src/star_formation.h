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
#ifndef SWIFT_STAR_FORMATION_H
#define SWIFT_STAR_FORMATION_H

/**
 * @file src/star_formation.h
 * @brief Branches between the different star formation recipies.
 */

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "hydro.h"

/* Import the right star formation law definition */
#if defined(STAR_FORMATION_NONE)
#define swift_star_formation_model_creates_stars 1
#include "./star_formation/none/star_formation.h"
#elif defined(STAR_FORMATION_QLA)
#define swift_star_formation_model_creates_stars 0
#include "./star_formation/QLA/star_formation.h"
#elif defined(STAR_FORMATION_EAGLE)
#define swift_star_formation_model_creates_stars 1
#include "./star_formation/EAGLE/star_formation.h"
#elif defined(STAR_FORMATION_GEAR)
#define swift_star_formation_model_creates_stars 1
#include "./star_formation/GEAR/star_formation.h"
#elif defined(STAR_FORMATION_SIMBA)
#define swift_star_formation_model_creates_stars 1
#include "./star_formation/SIMBA/star_formation.h"
#else
#error "Invalid choice of star formation law"
#endif

/* General functions defined in the source file */
void starformation_init(struct swift_params* parameter_file,
                        const struct phys_const* phys_const,
                        const struct unit_system* us,
                        const struct hydro_props* hydro_props,
                        struct star_formation* starform);

void starformation_print(const struct star_formation* starform);

/* Dump store */
void starformation_struct_dump(const struct star_formation* starform,
                               FILE* stream);

void starformation_struct_restore(const struct star_formation* starform,
                                  FILE* stream);

#endif /* SWIFT_STAR_FORMATION_H */
