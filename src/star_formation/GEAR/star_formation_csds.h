/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STAR_FORMATION_GEAR_STAR_FORMATION_CSDS_H
#define SWIFT_STAR_FORMATION_GEAR_STAR_FORMATION_CSDS_H

/* Other Includes */
#include "csds_io.h"

#ifdef WITH_CSDS

/**
 * @brief Group all the fields for the star formation together.
 *
 * @param p The #spart
 * @param e The #engine
 * @param buffer Allocated buffer for writing the particle.
 *
 * @return Buffer after the bits written.
 */
INLINE static void *csds_star_formation_convert(const struct spart *sp,
                                                const struct engine *e,
                                                void *buffer) {
  /* Write the terms into the buffer */
  float *out = (float *)buffer;
  out[0] = sp->sf_data.birth_density;
  out[1] = sp->sf_data.birth_mass;
  long long *id = (long long *)(out + 2);
  *id = sp->sf_data.progenitor_id;

  return id + 1;
}

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_star_formation_define_fields(struct csds_field *fields) {

  /* Write all the fields together. */
  csds_define_field_from_function_stars(fields[0], "GEARStarFormation",
                                        csds_star_formation_convert,
                                        2 * sizeof(float) + sizeof(long long));

  return 1;
}

#endif  // WITH_CSDS
#endif  // SWIFT_STAR_FORMATION_NONE_STAR_FORMATION_CSDS_H
