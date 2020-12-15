/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Josh Borrow (josh@joshborrow.com)
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

#ifndef SWIFT_PARTICLE_SPLITTING_STRUCT_H
#define SWIFT_PARTICLE_SPLITTING_STRUCT_H

#include <stdint.h>

/**
 * @brief  Data stored in the xpart, spart, or bpart
 *         to track the split history of each particle.
 */
struct particle_splitting_data {
  /*! Particle ID of the progenitor */
  long long progenitor_id;

  /*! Binary tree used to show the outcome of splitting events */
  long long split_tree;

  /*! Number of times this particle has been split. */
  uint8_t split_count;
};

#endif /* SWIFT_PARTICLE_SPLITTING_STRUCT_H */