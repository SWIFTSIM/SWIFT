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
#ifndef SWIFT_SPHENIX_HYDRO_CSDS_H
#define SWIFT_SPHENIX_HYDRO_CSDS_H

/* Other Includes */
#include "hydro.h"
#include "csds_io.h"

#ifdef WITH_CSDS

/*
 * List of all possible mask.
 * Outside the module, only hydro_csds_field_count is used.
 */
enum hydro_csds_fields {
  hydro_csds_field_coordinates = 0,
  hydro_csds_field_velocities,
  hydro_csds_field_accelerations,
  hydro_csds_field_masses,
  hydro_csds_field_smoothing_lengths,
  hydro_csds_field_internal_energies,
  hydro_csds_field_particle_ids,
  hydro_csds_field_densities,
  /* Due to the low number of masks available,
   * we group the fields together: pressure, entopries,
   * viscosity and diffusion parameters, laplacian of the energy
   * and the velocity divergence along its derivative.
   */
  hydro_csds_field_secondary_fields,
  hydro_csds_field_count,
};

/* Name of each possible mask. */
extern const char *hydro_csds_field_names[hydro_csds_field_count];

/**
 * @brief Initialize the csds.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_csds_fields!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int hydro_csds_writer_populate_mask_data(
    struct mask_data *mask_data) {
  mask_data[hydro_csds_field_coordinates] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_coordinates],
      3 * sizeof(double));

  mask_data[hydro_csds_field_velocities] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_velocities],
      3 * sizeof(float));

  mask_data[hydro_csds_field_accelerations] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_accelerations],
      3 * sizeof(float));

  mask_data[hydro_csds_field_masses] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_masses], sizeof(float));

  mask_data[hydro_csds_field_smoothing_lengths] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_smoothing_lengths],
      sizeof(float));

  mask_data[hydro_csds_field_internal_energies] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_internal_energies],
      sizeof(float));

  mask_data[hydro_csds_field_particle_ids] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_particle_ids],
      sizeof(long long));

  mask_data[hydro_csds_field_densities] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_densities], sizeof(float));

  const size_t size_secondary = 7 * sizeof(float);
  mask_data[hydro_csds_field_secondary_fields] = csds_create_mask_entry(
      hydro_csds_field_names[hydro_csds_field_secondary_fields],
      size_secondary);

  return hydro_csds_field_count;
}

/**
 * @brief Generates the mask and compute the size of the record.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_csds_fields!
 *
 * @param masks The list of masks (same order than in #hydro_csds_init).
 * @param part The #part that will be written.
 * @param xpart The #xpart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void hydro_csds_compute_size_and_mask(
    const struct mask_data *masks, const struct part *part,
    const struct xpart *xpart, const int write_all, size_t *buffer_size,
    unsigned int *mask) {

  /* Here you can decide your own writing logic */

  /* Add the coordinates. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_coordinates],
                                    buffer_size);

  /* Add the velocities. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_velocities],
                                    buffer_size);

  /* Add the accelerations. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_accelerations],
                                    buffer_size);

  /* Add the masses. */
  *mask |=
      csds_add_field_to_mask(masks[hydro_csds_field_masses], buffer_size);

  /* Add the smoothing lengths. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_smoothing_lengths],
                                    buffer_size);

  /* Add the energies. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_internal_energies],
                                    buffer_size);

  /* Add the ID. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_particle_ids],
                                    buffer_size);

  /* Add the density. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_densities],
                                    buffer_size);

  /* Add the secondary fields. */
  *mask |= csds_add_field_to_mask(masks[hydro_csds_field_secondary_fields],
                                    buffer_size);
}

/**
 * @brief Write a particle to the csds.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_csds_fields!
 *
 * @param masks The list of masks (same order than in #hydro_csds_init).
 * @param p The #part to write.
 * @param xp The #xpart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *hydro_csds_write_particle(
    const struct mask_data *mask_data, const struct part *p,
    const struct xpart *xp, unsigned int *mask, char *buff) {

  /* Write the coordinate. */
  if (csds_should_write_field(mask_data[hydro_csds_field_coordinates],
                                mask)) {
    memcpy(buff, p->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Write the velocity. */
  if (csds_should_write_field(mask_data[hydro_csds_field_velocities],
                                mask)) {
    memcpy(buff, p->v, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the acceleration. */
  if (csds_should_write_field(mask_data[hydro_csds_field_accelerations],
                                mask)) {

    /* Compute the acceleration due to hydro and gravity */
    float *acc = (float *)buff;
    acc[0] = p->a_hydro[0];
    acc[1] = p->a_hydro[1];
    acc[2] = p->a_hydro[2];
    if (p->gpart) {
      acc[0] += p->gpart->a_grav[0];
      acc[1] += p->gpart->a_grav[1];
      acc[2] += p->gpart->a_grav[2];
    }

    memcpy(buff, acc, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Write the mass. */
  if (csds_should_write_field(mask_data[hydro_csds_field_masses], mask)) {
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the smoothing length. */
  if (csds_should_write_field(mask_data[hydro_csds_field_smoothing_lengths],
                                mask)) {
    memcpy(buff, &p->h, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the energy. */
  if (csds_should_write_field(mask_data[hydro_csds_field_internal_energies],
                                mask)) {
    memcpy(buff, &p->u, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the Id. */
  if (csds_should_write_field(mask_data[hydro_csds_field_particle_ids],
                                mask)) {
    memcpy(buff, &p->id, sizeof(long long));
    buff += sizeof(long long);
  }

  /* Write the density. */
  if (csds_should_write_field(mask_data[hydro_csds_field_densities],
                                mask)) {
    memcpy(buff, &p->rho, sizeof(float));
    buff += sizeof(float);
  }

  /* Write the secondary fields. */
  if (csds_should_write_field(mask_data[hydro_csds_field_secondary_fields],
                                mask)) {
    const float secondary[7] = {
        hydro_get_comoving_entropy(p, xp),
        hydro_get_comoving_pressure(p),
        p->viscosity.alpha * p->force.balsara,
        p->diffusion.alpha,
        p->diffusion.laplace_u,
        p->viscosity.div_v,
        p->viscosity.div_v_dt,
    };

    memcpy(buff, &secondary, sizeof(secondary));
    buff += sizeof(secondary);
  }

  return buff;
}

#endif  // WITH_CSDS
#endif  // SWIFT_SPHENIX_HYDRO_CSDS_H
