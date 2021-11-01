/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Willem Elbers (willem.h.elbers@durham.ac.uk)
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
#ifndef SWIFT_NONE_NEUTRINO_IO_H
#define SWIFT_NONE_NEUTRINO_IO_H

#include "../../gravity.h"

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param gparts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int neutrino_write_particles(
    const struct gpart* gparts, struct io_props* list) {

  return 0;
}

#endif /* SWIFT_NONE_NEUTRINO_IO_H */
