/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_FEEDBACK_GEAR_RADIATION_SELECTION_H
#define SWIFT_FEEDBACK_GEAR_RADIATION_SELECTION_H

/**
 * @file src/feedback/GEAR/radiation_selection.h
 * @brief The parts of the GEAR subgrid radiation compiled into this build
 * (./configure --with-subgrid-radiation) and the reading of their run-time
 * switches.
 *
 * Standalone, so that the cooling and stars modules can read a switch too.
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "error.h"
#include "inline.h"
#include "parser.h"

/* Standard includes */
#include <string.h>

#ifdef GEAR_SUBGRID_RADIATION_PRESSURE
#define RADIATION_COMPILED_PRESSURE 1
#else
#define RADIATION_COMPILED_PRESSURE 0
#endif

#ifdef GEAR_SUBGRID_RADIATION_HII
#define RADIATION_COMPILED_HII 1
#else
#define RADIATION_COMPILED_HII 0
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF
#define RADIATION_COMPILED_ISRF 1
#else
#define RADIATION_COMPILED_ISRF 0
#endif

/*! Bit of a part in #radiation_selection_absent_mask(). */
#define RADIATION_SELECTION_BIT_PRESSURE 1
#define RADIATION_SELECTION_BIT_HII 2
#define RADIATION_SELECTION_BIT_ISRF 4

/**
 * @brief The parts of the subgrid radiation NOT compiled into this build.
 *
 * @return Sum of the RADIATION_SELECTION_BIT_* of the absent parts; 0 for
 * --with-subgrid-radiation=all.
 */
__attribute__((always_inline)) INLINE static int
radiation_selection_absent_mask(void) {
  return (RADIATION_COMPILED_PRESSURE ? 0 : RADIATION_SELECTION_BIT_PRESSURE) |
         (RADIATION_COMPILED_HII ? 0 : RADIATION_SELECTION_BIT_HII) |
         (RADIATION_COMPILED_ISRF ? 0 : RADIATION_SELECTION_BIT_ISRF);
}

/**
 * @brief Read the on/off switch of a part of the subgrid radiation.
 *
 * @param params The parsed parameters.
 * @param name The switch, e.g. "GEARFeedback:with_photoionization".
 * @param compiled Is the part compiled into this build?
 * @return The value of the switch (default 0), or 0 if the part is absent.
 */
__attribute__((always_inline)) INLINE static int radiation_selection_get_switch(
    struct swift_params *params, const char *name, const int compiled) {
  if (!compiled) return 0;
  return parser_get_opt_param_int(params, name, 0);
}

/**
 * @brief Warn once if the parameter file sets keys of an absent part.
 *
 * A key belongs to the part if its name starts with one of @p prefixes. The
 * keys are ignored, never an error.
 *
 * @param params The parsed parameters.
 * @param part The part's configure name (rp, hii or isrf).
 * @param prefixes The key prefixes of the part, NULL-terminated.
 */
__attribute__((always_inline)) INLINE static void
radiation_selection_warn_absent_part_keys(const struct swift_params *params,
                                          const char *part,
                                          const char *const *prefixes) {
  int count = 0;
  const char *first = NULL;
  for (int i = 0; i < params->paramCount; i++) {
    for (int k = 0; prefixes[k] != NULL; k++) {
      if (strncmp(params->data[i].name, prefixes[k], strlen(prefixes[k])) ==
          0) {
        if (first == NULL) first = params->data[i].name;
        count++;
        break;
      }
    }
  }
  if (count > 0 && engine_rank == 0)
    warning(
        "%d parameter(s) of the subgrid radiation part '%s' (first: %s) are "
        "ignored: this build has no '%s' (./configure "
        "--with-subgrid-radiation).",
        count, part, first, part);
}

/**
 * @brief Warn once per absent part whose keys the parameter file sets.
 *
 * @param params The parsed parameters.
 */
__attribute__((always_inline)) INLINE static void
radiation_selection_warn_absent_parts(const struct swift_params *params) {
  if (!RADIATION_COMPILED_PRESSURE) {
    const char *const keys[] = {"GEARFeedback:with_radiation_pressure",
                                "GEARFeedback:radiation_pressure_", NULL};
    radiation_selection_warn_absent_part_keys(params, "rp", keys);
  }
  if (!RADIATION_COMPILED_HII) {
    const char *const keys[] = {"GEARFeedback:with_photoionization",
                                "GEARFeedback:HII_", "Stars:HII_", NULL};
    radiation_selection_warn_absent_part_keys(params, "hii", keys);
  }
  if (!RADIATION_COMPILED_ISRF) {
    const char *const keys[] = {
        "GEARFeedback:with_interstellar_radiation_field", "GEARFeedback:ISRF_",
        NULL};
    radiation_selection_warn_absent_part_keys(params, "isrf", keys);
  }
}

#endif /* SWIFT_FEEDBACK_GEAR_RADIATION_SELECTION_H */
