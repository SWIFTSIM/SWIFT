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
#include "timeline.h"

/* Forward declarations */
struct engine;
struct space;
struct gpart;
struct threadpool;

/**
 * @brief Data structure for the long-range periodic forces using a mesh
 */
struct pm_mesh {

  /*! Is the calculation using periodic BCs? */
  int periodic;

  /*! The number of threads used by the FFTW library */
  int nr_threads;

  /*! Side-length of the mesh */
  int N;

  /*! Integer time-step end of the mesh force for the last step */
  integertime_t ti_end_mesh_last;

  /*! Integer time-step beginning of the mesh force for the last step */
  integertime_t ti_beg_mesh_last;

  /*! Integer time-step end of the mesh force for the next step*/
  integertime_t ti_end_mesh_next;

  /*! Integer time-step beginning of the mesh force for the next step */
  integertime_t ti_beg_mesh_next;

  /*! Conversion factor between box and mesh size */
  double cell_fac;

  /*! (Comoving) side-length of the box along the three axis */
  double dim[3];

  /*! Scale over which we smooth the forces */
  double r_s;

  /*! Inverse of the scale over which we smooth the forces */
  double r_s_inv;

  /*! Distance beyond which tree forces are neglected */
  double r_cut_max;

  /*! Distance below which tree forces are Newtonian */
  double r_cut_min;

  /*! Potential field */
  double *potential;
};

void pm_mesh_init(struct pm_mesh *mesh, const struct gravity_props *props,
                  const double dim[3], int nr_threads);
void pm_mesh_init_no_mesh(struct pm_mesh *mesh, double dim[3]);
void pm_mesh_compute_potential(struct pm_mesh *mesh, const struct space *s,
                               struct threadpool *tp, int verbose);
void pm_mesh_clean(struct pm_mesh *mesh);

void pm_mesh_allocate(struct pm_mesh *mesh);
void pm_mesh_free(struct pm_mesh *mesh);

/* Dump/restore. */
void pm_mesh_struct_dump(const struct pm_mesh *p, FILE *stream);
void pm_mesh_struct_restore(struct pm_mesh *p, FILE *stream);

#endif /* SWIFT_MESH_GRAVITY_H */
