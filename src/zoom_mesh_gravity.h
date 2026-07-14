/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk).
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
#ifndef SWIFT_ZOOM_MESH_GRAVITY_H
#define SWIFT_ZOOM_MESH_GRAVITY_H

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <stdio.h>

/* Forward declarations. */
struct cell;
struct gravity_props;
struct space;
struct swift_params;
struct threadpool;

/**
 * @brief Data structure for the high-resolution zoom gravity mesh.
 */
struct zoom_pm_mesh {

  /*! Is the zoom mesh enabled? */
  int enabled;

  /*! Side-length of the mesh along each axis. */
  int N[3];

  /*! Lower corner of the mesh in box coordinates. */
  double loc[3];

  /*! Side-length of the mesh region along each axis. */
  double dim[3];

  /*! Conversion factor between box and mesh coordinates. */
  double cell_fac[3];

  /*! Scale over which we smooth the forces. */
  double r_s;

  /*! Inverse of the scale over which we smooth the forces. */
  double r_s_inv;

  /*! Distance beyond which tree forces are handled by the zoom mesh. */
  double r_cut_max;

  /*! Distance below which tree forces are retained. */
  double r_cut_min;

  /*! Full zero-padded mesh used as density and potential field. */
  double *potential_global;
};

/* Prototypes. */
int zoom_mesh_can_use_mesh(const struct zoom_pm_mesh *mesh,
                           const struct space *s, const struct cell *ci,
                           const struct cell *cj);
int zoom_mesh_can_use_mesh_between_rebuilds(const struct zoom_pm_mesh *mesh,
                                            const struct space *s,
                                            const struct cell *ci,
                                            const struct cell *cj);
void zoom_mesh_init(struct zoom_pm_mesh *mesh, struct swift_params *params,
                    const struct gravity_props *props, const struct space *s,
                    const int verbose);
void zoom_mesh_clean(struct zoom_pm_mesh *mesh);
void zoom_mesh_compute_potential(struct zoom_pm_mesh *mesh,
                                 const struct space *s, struct threadpool *tp,
                                 int verbose);
int zoom_mesh_cell_is_covered(const struct zoom_pm_mesh *mesh,
                              const struct cell *c, const int use_max_dx);
int zoom_mesh_cells_are_covered(const struct zoom_pm_mesh *mesh,
                                const struct cell *ci, const struct cell *cj,
                                const int use_max_dx);
void zoom_mesh_struct_dump(const struct zoom_pm_mesh *mesh, FILE *stream);
void zoom_mesh_struct_restore(struct zoom_pm_mesh *mesh, FILE *stream);

#endif /* SWIFT_ZOOM_MESH_GRAVITY_H */
