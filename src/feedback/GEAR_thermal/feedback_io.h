/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FEEDBACK_IO_GEAR_H
#define SWIFT_FEEDBACK_IO_GEAR_H

#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from the ICs.
 *
 * Reads "FUVSpecificEnergy"/"LWSpecificEnergy" (singular: the codebase's
 * own IC-read convention, e.g. hydro's "Density"/"SmoothingLength", vs. the
 * plural "Densities"/"SmoothingLengths" used for the matching snapshot
 * output) as OPTIONAL fields into #feedback_part_data.u_FUV/u_LW. The
 * snapshot output for the same quantity
 * (#convert_part_u_FUV/convert_part_u_LW, src/tracers/GEAR/tracers_io.h)
 * deliberately uses the plural "FUVSpecificEnergies"/"LWSpecificEnergies"
 * instead, so a snapshot cannot be fed back in as an IC unmodified -- this
 * is a validation/testing tool, not a normal production IC input, and is
 * not meant to make that round-trip easy. It lets a test set an arbitrary,
 * analytically-known initial LW/FUV field shape (a value range across
 * particles, a pulse, a step) and watch only Grackle's chemistry, or only
 * the LW/FUV propagation PDE, evolve it, decoupled from the star and
 * injection machinery. It only makes physical sense in a run with no star,
 * though this reader does not itself enforce that.
 *
 * An IC without these fields is unaffected: #radiation_first_init_part no
 * longer zeroes #feedback_part_data.u_FUV/u_LW (and seeds u_FUV_prev/
 * u_LW_prev from them, not from 0.f, so a supplied value also survives the
 * very first propagation update when `GEARFeedback:LW_FUV_propagation` is
 * on) so that a supplied value survives first-init, but every #part is
 * bzero'd before this read runs (single_io.c/parallel_io.c/serial_io.c),
 * so a missing field still leaves exactly 0.f, matching pre-existing
 * behaviour.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_particles(struct part *parts,
                                          struct io_props *list) {

  list[0] = io_make_input_field("FUVSpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.u_FUV);
  list[1] = io_make_input_field("LWSpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.u_LW);

  return 2;
}

#endif /* SWIFT_FEEDBACK_IO_GEAR_H */
