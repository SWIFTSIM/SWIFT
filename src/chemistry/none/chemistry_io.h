/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_IO_NONE_H
#define SWIFT_CHEMISTRY_IO_NONE_H

#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_read_particles(struct part* parts,
                                           struct io_props* list) {

  /* update list according to hydro_io */

  /* Return the number of fields to read */
  return 0;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_particles(const struct part* parts,
                                            const struct xpart* xparts,
                                            struct io_props* list,
                                            const int with_cosmology) {

  /* update list according to hydro_io */

  /* Return the number of fields to write */
  return 0;
}

/**
 * @brief Specifies which sparticle fields to write to a dataset
 *
 * @param sparts The sparticle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_sparticles(const struct spart* sparts,
                                             struct io_props* list) {

  /* update list according to hydro_io */

  /* Return the number of fields to write */
  return 0;
}

/**
 * @brief Specifies which bparticle fields to write to a dataset
 *
 * @param bparts The bparticle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_bparticles(const struct bpart* bparts,
                                             struct io_props* list) {

  /* update list according to hydro_io */

  /* Return the number of fields to write */
  return 0;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of chemistry to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The #engine.
 */
INLINE static void chemistry_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                           const struct engine* e) {

  io_write_attribute_s(h_grp, "Chemistry Model", "None");
}
#endif

#endif /* SWIFT_CHEMISTRY_IO_NONE_H */
