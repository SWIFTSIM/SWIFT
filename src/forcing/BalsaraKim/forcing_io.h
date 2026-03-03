/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Maarten Elion (elion@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_BALSARAKIM_IO_H
#define SWIFT_FORCING_BALSARAKIM_IO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "forcing.h"
#include "engine.h"
#include "io_properties.h"

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extended data particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
__attribute__((always_inline)) INLINE static int forcing_write_particles(
    const struct part *parts, const struct xpart *xparts,
    struct io_props *list) {

  list[0] = io_make_physical_output_field(
      "ForcingInjectedEnergies", FLOAT, 1, UNIT_CONV_ENERGY, 
      0.f, xparts, forcing_data.forcing_injected_energy,
      /*convertable_to_comoving=*/0,
      "Physical energies injected by the forcing scheme");

  return 1;
}

#endif /* SWIFT_FORCING_BALSARAKIM_IO_H */