/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_NONE_CHEMISTRY_CSDS_H
#define SWIFT_CHEMISTRY_NONE_CHEMISTRY_CSDS_H

/* Other Includes */
#include "csds_io.h"
#include "hydro.h"

#ifdef WITH_CSDS

/**
 * @brief Write the metallicities.
 *
 * @param p The #part
 * @param xp Its related #xpart
 * @param e The #engine
 * @param buffer Allocated buffer for writing the particle.
 *
 * @return Buffer after the bits written.
 */
INLINE static void *csds_chemistry_write_part(const struct part *p,
                                              const struct xpart *xp,
                                              const struct engine *e,
                                              void *buffer) {

  /* Add the smoothed metals */
  const size_t size = sizeof(p->chemistry_data.smoothed_metal_mass_fraction);
  memcpy(buffer, p->chemistry_data.smoothed_metal_mass_fraction, size);
  buffer += size;

  /* Add the metal mass */
  double *metals = buffer;
  const float m = hydro_get_mass(p);
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    metals[i] = p->chemistry_data.metal_mass[i] / m;
  }

  return metals + GEAR_CHEMISTRY_ELEMENT_COUNT;
}

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_chemistry_define_fields_parts(
    struct csds_field *fields) {

  /* Write the metallicities and non smoothed metallicities together. */
  csds_define_field_from_function_hydro(
      fields[0], "GEARChemistryParts", csds_chemistry_write_part,
      2 * GEAR_CHEMISTRY_ELEMENT_COUNT * sizeof(double));

  return 1;
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

  csds_define_hydro_standard_field(fields[0], "GEARChemistrySparts",
                                   struct spart,
                                   chemistry_data.metal_mass_fraction,
                                   /* saving_xpart */ 0);

  return 1;
}

#endif  // WITH_CSDS
#endif  // SWIFT_CHEMISTRY_NONE_CHEMISTRY_CSDS_H
