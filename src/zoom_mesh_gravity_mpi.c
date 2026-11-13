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
 *******************************************************************************/

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "zoom_mesh_gravity_mpi.h"

/* Local includes. */
#include "active.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "exchange_structs.h"
#include "lock.h"
#include "mesh_gravity_patch.h"
#include "mesh_gravity_sort.h"
#include "neutrino.h"
#include "part.h"
#include "periodic.h"
#include "row_major_id.h"
#include "space.h"
#include "threadpool.h"
#include "zoom_mesh_gravity.h"

/**
 * @brief Initialize a mesh patch to cover a cell on the zoom mesh
 *
 * Zoom mesh equivalent of pm_mesh_patch_init(). The patch coordinates are
 * indices into the zero-padded zoom mesh and the particles are not wrapped.
 *
 * @param patch A pointer to the mesh patch
 * @param mesh The #zoom_pm_mesh
 * @param cell The cell which the mesh should cover
 * @param boundary_size Size of the boundary layer to include
 */
static void zoom_mesh_patch_init(struct pm_mesh_patch *patch,
                                 const struct zoom_pm_mesh *mesh,
                                 const struct cell *cell,
                                 const int boundary_size) {

  const int gcount = cell->grav.count;
  const struct gpart *gparts = cell->grav.parts;
  const double fac = mesh->cell_fac;

  /* The FFTs are performed on the zero-padded mesh */
  patch->N = 2 * mesh->N;
  patch->fac = fac;

  /* No wrapping on the isolated mesh, positions are relative to its corner */
  for (int i = 0; i < 3; i++) {
    patch->wrap_min[i] = 0.;
    patch->wrap_max[i] = mesh->dim[i];
  }

  /* Find the extent of the particle distribution in the cell */
  double pos_min[3];
  double pos_max[3];
  for (int i = 0; i < 3; i++) {
    pos_min[i] = patch->wrap_max[i];
    pos_max[i] = patch->wrap_min[i];
  }
  for (int ipart = 0; ipart < gcount; ipart++) {

    const struct gpart *gp = &gparts[ipart];

    if (gp->time_bin == time_bin_inhibited) continue;

    for (int i = 0; i < 3; i++) {
      const double pos = gp->x[i] - mesh->loc[i];
      if (pos < pos_min[i]) pos_min[i] = pos;
      if (pos > pos_max[i]) pos_max[i] = pos;
    }
  }

  /* Determine the integer size and coordinates of the mesh */
  int num_cells = 1;
  for (int i = 0; i < 3; i++) {
    patch->mesh_min[i] = floor(pos_min[i] * fac) - boundary_size;
    /* CIC interpolation requires one extra element in the positive direction */
    patch->mesh_max[i] = floor(pos_max[i] * fac) + boundary_size + 1;
    patch->mesh_size[i] = patch->mesh_max[i] - patch->mesh_min[i] + 1;
    num_cells *= patch->mesh_size[i];
  }

  /* Allocate the mesh */
  if (swift_memalign("mesh_patch", (void **)&patch->mesh, SWIFT_CACHE_ALIGNMENT,
                     num_cells * sizeof(double)) != 0) {
    error("Failed to allocate array for mesh patch!");
  }
}

/**
 * @brief Accumulate contributions from cell to the zoom mesh density field
 *
 * Zoom mesh equivalent of accumulate_cell_to_local_patch().
 *
 * Allocates a temporary mesh which covers the top level cell and
 * accumulates mass contributions to this mesh.
 *
 * @param mesh The #zoom_pm_mesh
 * @param cell The #cell containing the particles.
 * @param patch The local mesh patch
 * @param nu_model Struct with neutrino constants
 *
 */
void zoom_mesh_accumulate_cell_to_local_patch(
    const struct zoom_pm_mesh *mesh, const struct cell *cell,
    struct pm_mesh_patch *patch, const struct neutrino_model *nu_model) {

  /* If the cell is empty, then there's nothing to do
     (and the code to find the extent of the cell would fail) */
  if (cell->grav.count == 0) return;

  const double fac = mesh->cell_fac;

  /* Initialise the local mesh patch */
  zoom_mesh_patch_init(patch, mesh, cell, /*boundary_size=*/1);
  pm_mesh_patch_zero(patch);

  /* Loop over particles in this cell */
  for (int ipart = 0; ipart < cell->grav.count; ipart++) {

    const struct gpart *gp = &(cell->grav.parts[ipart]);

    /* Position relative to the corner of the zoom mesh */
    const double pos_x = gp->x[0] - mesh->loc[0];
    const double pos_y = gp->x[1] - mesh->loc[1];
    const double pos_z = gp->x[2] - mesh->loc[2];

    /* Workout the CIC coefficients */
    int i = (int)floor(fac * pos_x);
    const double dx = fac * pos_x - i;
    const double tx = 1. - dx;

    int j = (int)floor(fac * pos_y);
    const double dy = fac * pos_y - j;
    const double ty = 1. - dy;

    int k = (int)floor(fac * pos_z);
    const double dz = fac * pos_z - k;
    const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
    if (i < 0 || i >= mesh->N) error("Invalid gpart position in x");
    if (j < 0 || j >= mesh->N) error("Invalid gpart position in y");
    if (k < 0 || k >= mesh->N) error("Invalid gpart position in z");
#endif

    /* Get coordinates within the mesh patch */
    const int ii = i - patch->mesh_min[0];
    const int jj = j - patch->mesh_min[1];
    const int kk = k - patch->mesh_min[2];

    /* Compute weight (for neutrino delta-f weighting) */
    double weight = 1.0;
    if (gp->type == swift_type_neutrino) {
      gpart_neutrino_weight_mesh_only(gp, nu_model, &weight);
    }

    /* Accumulate contributions to the local mesh patch */
    const double mass = gp->mass;
    const double value = mass * weight;
    pm_mesh_patch_CIC_set(patch, ii, jj, kk, tx, ty, tz, dx, dy, dz, value);
  }
}
