/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_HYDRO_SPACE_H
#define SWIFT_HYDRO_SPACE_H

#include "../config.h"

struct space;

/**
 * @brief Extra space information that is needed for some hydro schemes.
 */
#ifdef SHADOWFAX_SPH
struct hydro_space {
  /*! Anchor of the simulation space. */
  double anchor[3];

  /*! Side lengths of the simulation space. */
  double side[3];
};
#else
struct hydro_space {};
#endif

void hydro_space_init(struct hydro_space *hs, const struct space *s);

#endif /* SWIFT_HYDRO_SPACE_H */
