/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_MESH_GRAVITY_H
#define SWIFT_MESH_GRAVITY_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "gravity_properties.h"

/* Forward declarations */
struct engine;
struct gpart;

/**
 * @brief Data structure for the long-range periodic forces using a mesh
 */
struct pm_mesh {

  /*! Side-length of the mesh */
  int N;

  /*! Conversion factor between box and mesh size */
  double cell_fac;

  /*! (Comoving) side-length of the volume covered by the mesh */
  double box_size;

  /*! Scale over which we smooth the forces */
  double a_smooth;

  /*! Potential field */
  double *potential;
};

void pm_mesh_init(struct pm_mesh *mesh, const struct gravity_props *props,
                  double box_size);
void pm_mesh_compute_potential(struct pm_mesh *mesh, const struct engine *e);
void pm_mesh_interpolate_forces(const struct pm_mesh *mesh,
                                const struct engine *e, struct gpart *gparts,
                                int gcount);
void pm_mesh_clean(struct pm_mesh *mesh);

#endif /* SWIFT_MESH_GRAVITY_H */
