/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT developers.
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
#ifndef SWIFT_BOUNDARY_CONDITIONS_H
#define SWIFT_BOUNDARY_CONDITIONS_H

#include "periodic.h"

struct swift_params;

enum boundary_condition_type {
  boundary_condition_open,
  boundary_condition_periodic,
  boundary_condition_reflective
};

/* Bits for the six faces of the simulation box. */
enum boundary_wall {
  boundary_wall_x_min = 1 << 0,
  boundary_wall_x_max = 1 << 1,
  boundary_wall_y_min = 1 << 2,
  boundary_wall_y_max = 1 << 3,
  boundary_wall_z_min = 1 << 4,
  boundary_wall_z_max = 1 << 5
};

struct boundary_conditions {
  /*! Box boundary policy. */
  enum boundary_condition_type type;

  /*! Reflective faces; unlisted faces are open in reflective runs. */
  unsigned int walls;
};

void boundary_conditions_init(struct boundary_conditions *bc,
                              struct swift_params *params);

/** Normalize a coordinate when the particle is assigned a new cell. */
static inline double boundary_conditions_wrap_periodic(const double x,
                                                       const double dim) {
  const double wrapped = box_wrap(x, 0., dim);
  /* A rounded value on the upper face belongs to the lower one. */
  return wrapped == dim ? 0. : wrapped;
}

#endif /* SWIFT_BOUNDARY_CONDITIONS_H */
