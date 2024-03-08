/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ADAPTIVE_SOFTENING_STRUCT_H
#define SWIFT_ADAPTIVE_SOFTENING_STRUCT_H

/* Config parameters. */
#include <config.h>

#ifdef ADAPTIVE_SOFTENING

/**
 * @brief Particle-carried fields for the adaptive softening scheme.
 */
struct adaptive_softening_part_data {

  /*! Correction term for energy conservation */
  float zeta;
};

#else

/**
 * @brief Particle-carried fields for the adaptive softening scheme.
 */
struct adaptive_softening_part_data {};

#endif

#endif /* SWIFT_ADAPTIVE_SOFTENING_STRUCT_H */
