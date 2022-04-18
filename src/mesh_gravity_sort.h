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
#ifndef SWIFT_MESH_GRAVITY_SORT_H
#define SWIFT_MESH_GRAVITY_SORT_H

/* Config parameters. */
#include "../config.h"

/* Standard includes */
#include <stdlib.h>
#include <string.h>

struct threadpool;

/**
 * @brief Store contributions to the mesh as (index, mass) pairs
 */
struct mesh_key_value_rho {
  size_t key;
  double value;
};

/**
 * @brief Store contributions to the mesh as (cell, index, mass) tuples
 */
struct mesh_key_value_pot {
  size_t cell_index;
  size_t key;
  double value;
};

void bucket_sort_mesh_key_value_rho(const struct mesh_key_value_rho *array_in,
                                    const size_t count, const int N,
                                    struct threadpool *tp,
                                    struct mesh_key_value_rho *array_out,
                                    size_t *bucket_offsets);

void bucket_sort_mesh_key_value_pot(const struct mesh_key_value_pot *array_in,
                                    const size_t count, const int N,
                                    struct threadpool *tp,
                                    struct mesh_key_value_pot *array_out,
                                    size_t *bucket_offsets);

void bucket_sort_mesh_key_value_pot_index(
    const struct mesh_key_value_pot *array_in, const size_t count, const int N,
    struct threadpool *tp, struct mesh_key_value_pot *array_out);

#endif /* SWIFT_MESH_GRAVITY_SORT_H */
