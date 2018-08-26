/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "const.h"
#include "hydro.h"
#include "parser.h"
#include "units.h"

/* This object's header. */
#include "sourceterms.h"

/**
 * @brief Initialises the sourceterms
 *
 * @param parameter_file The parsed parameter file
 * @param us The current internal system of units
 * @param source the structure that has all the source term properties
 */
void sourceterms_init(struct swift_params *parameter_file,
                      struct unit_system *us, struct sourceterms *source) {
#ifdef SOURCETERMS_SN_FEEDBACK
  supernova_init(parameter_file, us, source);
#endif /* SOURCETERMS_SN_FEEDBACK */
};

/**
 * @brief Prints the properties of the source terms to stdout
 * @param source the structure that has all the source term properties
 */
void sourceterms_print(struct sourceterms *source) {
#ifdef SOURCETERMS_NONE
  error(" no sourceterms defined yet you ran with -F");
#ifdef SOURCETERMS_SN_FEEDBACK
#error "can't have sourceterms when defined SOURCETERMS_NONE"
#endif
#endif
#ifdef SOURCETERMS_SN_FEEDBACK
  supernova_print(source);
#endif /* SOURCETERMS_SN_FEEDBACK */
};

/**
 * @brief Write a sourceterms struct to the given FILE as a stream of bytes.
 *
 * @param sourceterms the struct
 * @param stream the file stream
 */
void sourceterms_struct_dump(const struct sourceterms *sourceterms,
                             FILE *stream) {
  restart_write_blocks((void *)sourceterms, sizeof(struct sourceterms), 1,
                       stream, "sourceterms", "sourceterms");
}

/**
 * @brief Restore a sourceterms struct from the given FILE as a stream of
 * bytes.
 *
 * @param sourceterms the struct
 * @param stream the file stream
 */
void sourceterms_struct_restore(const struct sourceterms *sourceterms,
                                FILE *stream) {
  restart_read_blocks((void *)sourceterms, sizeof(struct sourceterms), 1,
                      stream, NULL, "sourceterms");
}
