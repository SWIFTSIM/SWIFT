/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_H
#define SWIFT_CHEMISTRY_H

/**
 * @file src/chemistry.h
 * @brief Branches between the different chemistry functions.
 */

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "chemistry_struct.h"

/* Import the right chemistry definition */
#if defined(CHEMISTRY_NONE)
#include "./chemistry/none/chemistry.h"
#include "./chemistry/none/chemistry_iact.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/GEAR/chemistry.h"
#include "./chemistry/GEAR/chemistry_iact.h"
#elif defined(CHEMISTRY_GEAR_DIFFUSION)
#include "./chemistry/GEAR_DIFFUSION/chemistry.h"
#include "./chemistry/GEAR_DIFFUSION/chemistry_iact.h"
#elif defined(CHEMISTRY_AGORA)
#include "./chemistry/AGORA/chemistry.h"
#include "./chemistry/AGORA/chemistry_iact.h"
#elif defined(CHEMISTRY_QLA)
#include "./chemistry/QLA/chemistry.h"
#include "./chemistry/QLA/chemistry_iact.h"
#elif defined(CHEMISTRY_EAGLE)
#include "./chemistry/EAGLE/chemistry.h"
#include "./chemistry/EAGLE/chemistry_iact.h"
#else
#error "Invalid choice of chemistry function."
#endif

/* Common functions */
void chemistry_init(struct swift_params* parameter_file,
                    const struct unit_system* us,
                    const struct phys_const* phys_const,
                    struct chemistry_global_data* data);

void chemistry_print(const struct chemistry_global_data* data);

/* Dump/restore. */
void chemistry_struct_dump(const struct chemistry_global_data* chemistry,
                           FILE* stream);
void chemistry_struct_restore(const struct chemistry_global_data* chemistry,
                              FILE* stream);

#endif /* SWIFT_CHEMISTRY_H */
