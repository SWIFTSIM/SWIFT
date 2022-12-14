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
#ifndef SWIFT_MULTISOFTENING_GRAVITY_CSDS_H
#define SWIFT_MULTISOFTENING_GRAVITY_CSDS_H

#include "csds_io.h"
#include "gravity_part.h"

#ifdef WITH_CSDS

/**
 * @brief Compute the acceleration and writes it.
 *
 * @param gp The #gpart
 * @param e The #engine
 * @param buffer Allocated buffer for writing the particle.
 *
 * @return Buffer after the bits written.
 */
INLINE static void *csds_gravity_convert_acc(const struct gpart *gp,
                                             const struct engine *e,
                                             void *buffer) {
  /* Compute the acceleration due to hydro and gravity */
  float *acc = (float *)buffer;
  acc[0] = gp->a_grav[0] + gp->a_grav_mesh[0];
  acc[1] = gp->a_grav[1] + gp->a_grav_mesh[1];
  acc[2] = gp->a_grav[2] + gp->a_grav_mesh[2];

  return acc + 3;
}

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_gravity_define_fields(struct csds_field *fields) {

  /* Positions */
  csds_define_standard_field(fields[0], "Coordinates", struct gpart, x);

  /* Velocities */
  csds_define_standard_field(fields[1], "Velocities", struct gpart, v_full);

  /* Accelerations */
  struct gpart p;
  csds_define_field_from_function_gravity(
      fields[2], "Accelerations", csds_gravity_convert_acc, sizeof(p.a_grav));

  /* Masses */
  csds_define_standard_field(fields[3], "Masses", struct gpart, mass);

  /* Particle IDs */
  csds_define_standard_field(fields[4], "ParticleIDs", struct gpart,
                             id_or_neg_offset);

  return 5;
}
#endif  // WITH_CSDS
#endif  // SWIFT_MULTISOFTENING_GRAVITY_CSDS_H
