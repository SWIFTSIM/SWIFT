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
#ifndef SWIFT_COOLING_H
#define SWIFT_COOLING_H

/**
 * @file src/cooling.h
 * @brief Branches between the different cooling functions.
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "inline.h"
#include "parser.h"
#include "physical_constants.h"
#include "restart.h"
#include "space.h"
#include "units.h"

/* Import the right cooling definition */
#if defined(COOLING_NONE)
#include "./cooling/none/cooling.h"
#elif defined(COOLING_CONST_DU)
#include "./cooling/const_du/cooling.h"
#elif defined(COOLING_CONST_LAMBDA)
#include "./cooling/const_lambda/cooling.h"
#elif defined(COOLING_COMPTON)
#include "./cooling/Compton/cooling.h"
#elif defined(COOLING_GRACKLE)
#include "./cooling/grackle/cooling.h"
#elif defined(COOLING_GEAJAY)
#include "./cooling/geajay/cooling.h"
#elif defined(COOLING_QLA)
#include "./cooling/QLA/cooling.h"
#elif defined(COOLING_QLA_EAGLE)
#include "./cooling/QLA_EAGLE/cooling.h"
#elif defined(COOLING_EAGLE)
#include "./cooling/EAGLE/cooling.h"
#elif defined(COOLING_PS2020)
#include "./cooling/PS2020/cooling.h"
#elif defined(COOLING_TREECOOL)
#include "./cooling/TREECOOL/cooling.h"
#else
#error "Invalid choice of cooling function."
#endif

/* A module defines COOLING_MISSING_INIT_MESSAGE in its header to name its own
 * parameter block and keys in cooling_get_missing_init_message(). */
#ifndef COOLING_MISSING_INIT_MESSAGE
#define COOLING_MISSING_INIT_MESSAGE                                       \
  "A feedback channel that needs the cooling module initialised "          \
  "(e.g. GEARFeedback:with_interstellar_radiation_field) is enabled, but " \
  "neither --cooling nor --temperature was passed. Pass --cooling or "     \
  "--temperature."
#endif

/**
 * @brief Error message for a feedback channel that needs the cooling module
 * initialised when neither --cooling nor --temperature was passed.
 *
 * @return The message text, with the module's own parameter block and keys.
 */
__attribute__((always_inline)) INLINE static const char *
cooling_get_missing_init_message(void) {
  return COOLING_MISSING_INIT_MESSAGE;
}

/* Common functions */
void cooling_init(struct swift_params *parameter_file,
                  const struct unit_system *us,
                  const struct phys_const *phys_const,
                  const struct hydro_props *hydro_props,
                  struct cooling_function_data *cooling);

void cooling_print(const struct cooling_function_data *cooling);

/* Dump/restore. */
void cooling_struct_dump(const struct cooling_function_data *cooling,
                         FILE *stream);
void cooling_struct_restore(struct cooling_function_data *cooling, FILE *stream,
                            const struct cosmology *cosmo);

#endif /* SWIFT_COOLING_H */
