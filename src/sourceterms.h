/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_SOURCETERMS_H
#define SWIFT_SOURCETERMS_H

/**
 * @file src/sourceterms.h
 * @brief Branches between the different sourceterms functions.
 */

#include "./const.h"
#include "runner.h"

#ifdef SOURCETERMS_SN_FEEDBACK
#include "sourceterms/sn_feedback/sn_feedback_struct.h"
#endif

/* So far only one model here */
struct sourceterms {
#ifdef SOURCETERMS_SN_FEEDBACK
  struct supernova_struct supernova;
#endif
};
#ifdef SOURCETERMS_SN_FEEDBACK
#include "sourceterms/sn_feedback/sn_feedback.h"
#endif

void sourceterms_init(const struct swift_params* parameter_file,
                      struct unit_system* us, struct sourceterms* source);
void sourceterms_print(struct sourceterms* source);

/* Dump/restore. */
void sourceterms_struct_dump(const struct sourceterms* source, FILE* stream);
void sourceterms_struct_restore(const struct sourceterms* source, FILE* stream);

/**
 * @brief Routines related to source terms
 * @param cell_min: corner of cell to test
 * @param cell_width: width of cell to test
 * @param sourceterms: properties of source terms to test
 * @param dimen: dimensionality of the problem
 *
 * This routine tests whether a source term should be applied to this cell
 * return: 1 if yes, return: 0 if no
 */

__attribute__((always_inline)) INLINE static int sourceterms_test_cell(
    const double cell_min[], const double cell_width[],
    struct sourceterms* sourceterms, const int dimen) {
#ifdef SOURCETERMS_SN_FEEDBACK
  return supernova_feedback_test_cell(cell_min, cell_width, sourceterms, dimen);
#endif
  return 0;
};

__attribute__((always_inline)) INLINE static void sourceterms_apply(
    struct runner* r, struct sourceterms* sourceterms, struct cell* c) {
#ifdef SOURCETERMS_SN_FEEDBACK
  supernova_feedback_apply(r, sourceterms, c);
#endif
};
#endif /*  SWIFT_SOURCETERMS_H */
