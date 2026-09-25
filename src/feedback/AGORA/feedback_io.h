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
#ifndef SWIFT_FEEDBACK_IO_AGORA_H
#define SWIFT_FEEDBACK_IO_AGORA_H

#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from the ICs. Nothing to do
 * here.
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_particles(struct part *parts,
                                          struct io_props *list) {
  return 0;
}

/**
 * @brief Specifies which particle fields to write to a dataset. Nothing to
 * do here.
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int feedback_write_particles(const struct part *parts,
                                           const struct xpart *xparts,
                                           struct io_props *list,
                                           const int with_cosmology) {
  return 0;
}

/**
 * @brief Specifies which star particle fields to write to a dataset.
 * Nothing to do here.
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology switched on?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int feedback_write_sparticles(const struct spart *sparts,
                                            struct io_props *list,
                                            const int with_cosmology) {
  return 0;
}

/**
 * @brief Specifies which star particle fields to read from the ICs or a
 * restart file. Nothing to do here.
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int feedback_read_sparticles(struct spart *sparts,
                                           struct io_props *list) {
  return 0;
}

#endif /* SWIFT_FEEDBACK_IO_AGORA_H */
