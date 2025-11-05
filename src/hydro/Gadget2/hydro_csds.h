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
#ifndef SWIFT_GADGET2_HYDRO_CSDS_H
#define SWIFT_GADGET2_HYDRO_CSDS_H

/* Other Includes */
#include "csds_io.h"
#include "hydro.h"

#ifdef WITH_CSDS

/**
 * @brief Compute the acceleration and writes it.
 *
 * @param p The #part
 * @param xp Its related #xpart
 * @param e The #engine
 * @param buffer Allocated buffer for writing the particle.
 *
 * @return Buffer after the bits written.
 */
INLINE static void *csds_hydro_convert_acc(const struct part *p,
                                           const struct xpart *xp,
                                           const struct engine *e,
                                           void *buffer) {
  /* Compute the acceleration due to hydro and gravity */
  float *acc = (float *)buffer;

  /* The hydro and gravity do not have the same factors */
  /* Convert everything into gravity acceleration */
  const double factor =
      e->cosmology->a_factor_hydro_accel / e->cosmology->a_factor_grav_accel;

  acc[0] = p->a_hydro[0] * factor;
  acc[1] = p->a_hydro[1] * factor;
  acc[2] = p->a_hydro[2] * factor;
  if (p->gpart) {
    acc[0] += p->gpart->a_grav[0] + p->gpart->a_grav_mesh[0];
    acc[1] += p->gpart->a_grav[1] + p->gpart->a_grav_mesh[1];
    acc[2] += p->gpart->a_grav[2] + p->gpart->a_grav_mesh[2];
  }

  return acc + 3;
}

/**
 * @brief Defines the fields to write in the CSDS.
 *
 * @param fields (output) The list of fields to write (already allocated).
 *
 * @return The number of fields.
 */
INLINE static int csds_hydro_define_fields(struct csds_field *fields) {

  /* Positions */
  csds_define_hydro_standard_field(fields[0], "Coordinates", struct part, x,
                                   /* saving_xpart */ 0);

  /* Velocities */
  csds_define_hydro_standard_field(fields[1], "Velocities", struct part, v,
                                   /* saving_xpart */ 0);

  /* Accelerations */
  struct part p;
  csds_define_field_from_function_hydro(
      fields[2], "Accelerations", csds_hydro_convert_acc, sizeof(p.a_hydro));

  /* Masses */
  csds_define_hydro_standard_field(fields[3], "Masses", struct part, mass,
                                   /* saving_xpart */ 0);

  /* Smoothing lengths */
  csds_define_hydro_standard_field(fields[4], "SmoothingLengths", struct part,
                                   h, /* saving_xpart */ 0);

  /* Entropies */
  csds_define_hydro_standard_field(fields[5], "Entropies", struct part, entropy,
                                   /* saving_xpart */ 0);

  /* Particle IDs */
  csds_define_hydro_standard_field(fields[6], "ParticleIDs", struct part, id,
                                   /* saving_xpart */ 0);

  /* Densities */
  csds_define_hydro_standard_field(fields[7], "Densities", struct part, rho,
                                   /* saving_xpart */ 0);

  return 8;
}

#endif  // WITH_CSDS
#endif  // SWIFT_GADGET2_HYDRO_CSDS_H
