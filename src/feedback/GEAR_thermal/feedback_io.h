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
 * Reads "PESpecificEnergy"/"LWSpecificEnergy" (singular: the codebase's
 * own IC-read convention, e.g. hydro's "Density"/"SmoothingLength", vs. the
 * plural "Densities"/"SmoothingLengths" used for the matching snapshot
 * output) as OPTIONAL fields into #feedback_isrf_moment_data.u.
 * The snapshot output for the same quantity
 * (#convert_part_u_PE/convert_part_u_LW, src/tracers/GEAR/tracers_io.h)
 * deliberately uses the plural "PESpecificEnergies"/"LWSpecificEnergies"
 * instead, so a snapshot cannot be fed back in as an IC unmodified. This
 * is a validation/testing tool, not a normal production IC input, and is
 * not meant to make that round-trip easy. It lets a test set an arbitrary,
 * analytically-known initial LW/PE field shape (a value range across
 * particles, a pulse, a step) and watch only Grackle's chemistry, or only
 * the LW/PE propagation PDE, evolve it, decoupled from the star and
 * injection machinery. It only makes physical sense in a run with no star,
 * though this reader does not itself enforce that.
 *
 * An IC without these fields is unaffected: #radiation_first_init_part no
 * longer zeroes #feedback_isrf_moment_data.u (and seeds
 * #feedback_isrf_moment_data.u_prev from it, not from 0.f, so a supplied value
 * also survives the very first propagation update when
 * `GEARFeedback:ISRF_propagation` is on) so that a supplied value survives
 * first-init, but every #part is bzero'd before this read runs
 * (single_io.c/parallel_io.c/serial_io.c), so a missing field still leaves
 * exactly 0.f, matching pre-existing behaviour.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_particles(struct part *parts,
                                          struct io_props *list) {

  /* Both are PHYSICAL, mass-specific quantities: no scale-factor exponent
     beyond what UNIT_CONV_ENERGY_PER_UNIT_MASS implies, and an input field
     carries no a-exponent slot at all, so an IC value is taken verbatim. */

  list[0] = io_make_input_field("PESpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_PE].u);
  list[1] = io_make_input_field("LWSpecificEnergy", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                feedback_data.isrf_moment[ISRF_MOMENT_LW].u);

  /* LWPhotonSpecificEnergy = 0 together with a nonzero LWSpecificEnergy is
     not representable after first-init: radiation_first_init_part()
     overwrites it with the LW value, for either sign, since a seeded LW
     field with no attribution is a field at the reference photon energy by
     definition. An IC author who wants no photon moment carried must also
     zero LWSpecificEnergy. */
  list[2] =
      io_make_input_field("LWPhotonSpecificEnergy", FLOAT, 1, OPTIONAL,
                          UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                          feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u);

  return 3;
}

#endif /* SWIFT_FEEDBACK_IO_GEAR_H */
