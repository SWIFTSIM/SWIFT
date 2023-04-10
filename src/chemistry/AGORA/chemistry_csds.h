/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_AGORA_CHEMISTRY_CSDS_H
#define SWIFT_CHEMISTRY_AGORA_CHEMISTRY_CSDS_H

#include "csds_io.h"

#ifdef WITH_CSDS

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_chemistry_define_fields_parts(
    struct csds_field *fields) {
  return 0;
}

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_chemistry_define_fields_sparts(
    struct csds_field *fields) {
  return 0;
}

#endif  // WITH_CSDS
#endif  // SWIFT_CHEMISTRY_AGORA_CHEMISTRY_CSDS_H
