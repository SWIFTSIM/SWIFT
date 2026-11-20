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

#include "minmax.h"
#include "periodic.h"

#include <math.h>

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

/* How a reflective wall acts on the fluid. */
enum boundary_wall_model {
  /*! Particles crossing a wall are bounced back; kernels are truncated. */
  boundary_wall_model_bounce,
  /*! As bounce, plus a layer of mirror particles outside the wall. */
  boundary_wall_model_mirror
};

struct boundary_conditions {
  /*! Box boundary policy. */
  enum boundary_condition_type type;

  /*! Reflective faces; unlisted faces are open in reflective runs. */
  unsigned int walls;

  /*! Treatment of the fluid at the reflective faces. */
  enum boundary_wall_model wall_model;
};

void boundary_conditions_init(struct boundary_conditions *bc,
                              struct swift_params *params);

/**
 * @brief Normalize a coordinate when the particle is assigned a new cell.
 *
 * Periodic boxes wrap. Other boxes never do: the boundary phase has already
 * reflected or removed anything outside, so only a particle exactly on the
 * upper face is left, and it is moved just inside.
 */
static inline double boundary_conditions_wrap(const double x, const double dim,
                                              const int periodic) {
  if (!periodic) return min(max(x, 0.), nextafter(dim, 0.));
  const double wrapped = box_wrap(x, 0., dim);
  /* A rounded value on the upper face belongs to the lower one. */
  return wrapped == dim ? 0. : wrapped;
}

/**
 * @brief Number of layers of mirror cells outside one face of the box.
 *
 * @param bc The boundary conditions.
 * @param axis The axis normal to the face (0, 1 or 2).
 * @param upper 0 for the lower face along that axis, 1 for the upper one.
 */
static inline int boundary_conditions_mirror_pad(
    const struct boundary_conditions *bc, const int axis, const int upper) {
  return bc->wall_model == boundary_wall_model_mirror &&
         ((bc->walls >> (2 * axis + upper)) & 1u);
}

/**
 * @brief Index along one axis of the cell reflected in a top-level cell.
 *
 * Cells inside the box map to themselves; a mirror cell maps to the cell at
 * the same distance on the other side of the wall.
 *
 * @param i Index of the cell along this axis of the (padded) grid.
 * @param cdim Number of cells along this axis, including the padding.
 * @param pad_lo Number of layers of mirror cells below the box.
 * @param pad_hi Number of layers of mirror cells above the box.
 */
static inline int boundary_conditions_mirror_index(const int i, const int cdim,
                                                   const int pad_lo,
                                                   const int pad_hi) {
  if (i < pad_lo) return 2 * pad_lo - 1 - i;
  if (i >= cdim - pad_hi) return 2 * (cdim - pad_hi) - 1 - i;
  return i;
}

#endif /* SWIFT_BOUNDARY_CONDITIONS_H */
