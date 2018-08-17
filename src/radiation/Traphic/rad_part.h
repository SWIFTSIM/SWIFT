/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_TRAPHIC_RAD_PART_H
#define SWIFT_TRAPHIC_RAD_PART_H

/* Some standard headers. */
#include <stdlib.h>

/**
 * @brief Particle fields for the radiation particles.
 */
struct rpart {

  /*! Particle position. */
  double x[3];

  /*! Particle smoothing length. */
  float h;

  /*! Particle density. */
  double density;

  /*! Particle hydrogen neutral fraction. */
  double hydrogen_neutral_fraction;

  /*! Particle received ionising luminosity. */
  double ionising_luminosity;

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_TRAPHIC_RAD_PART_H */
