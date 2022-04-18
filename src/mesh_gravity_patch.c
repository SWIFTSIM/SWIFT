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

/* Config parameters. */
#include "../config.h"

#include <math.h>

/* This object's header. */
#include "mesh_gravity_patch.h"

/* Local includes. */
#include "cell.h"
#include "error.h"
#include "row_major_id.h"

/**
 * @brief Initialize a mesh patch to cover a cell
 *
 * @param patch A pointer to the mesh patch
 * @param cell The cell which the mesh should cover
 * @param fac Inverse of the FFT mesh size
 * @param dim Size of the full volume in each dimension
 * @param boundary_size Size of the boundary layer to include
 */
void pm_mesh_patch_init(struct pm_mesh_patch *patch, const struct cell *cell,
                        const int N, const double fac, const double dim[3],
                        const int boundary_size) {

  const int gcount = cell->grav.count;
  const struct gpart *gparts = cell->grav.parts;

  patch->N = N;
  patch->fac = fac;

  /* Will need to wrap particles to position nearest the cell centre */
  double wrap_min[3];
  double wrap_max[3];
  for (int i = 0; i < 3; i++) {
    wrap_min[i] = cell->loc[i] + 0.5 * cell->width[i] - 0.5 * dim[i];
    wrap_max[i] = cell->loc[i] + 0.5 * cell->width[i] + 0.5 * dim[i];
  }

  /* Store the wraps */
  for (int i = 0; i < 3; i++) {
    patch->wrap_min[i] = wrap_min[i];
    patch->wrap_max[i] = wrap_max[i];
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
      const double pos_wrap = box_wrap(gp->x[i], wrap_min[i], wrap_max[i]);
      if (pos_wrap < pos_min[i]) pos_min[i] = pos_wrap;
      if (pos_wrap > pos_max[i]) pos_max[i] = pos_wrap;
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
                     num_cells * sizeof(double)) != 0)
    error("Failed to allocate array for mesh patch!");
}

/**
 * @brief Set all values in a mesh patch to zero
 *
 * @param patch A pointer to the mesh patch
 */
void pm_mesh_patch_zero(struct pm_mesh_patch *patch) {

  const int num =
      patch->mesh_size[0] * patch->mesh_size[1] * patch->mesh_size[2];
  memset(patch->mesh, 0, num * sizeof(double));
}

/**
 * @brief Free the memory associated with a mesh patch.
 *
 * @param patch A pointer to the mesh patch
 */
void pm_mesh_patch_clean(struct pm_mesh_patch *patch) {

  swift_free("mesh_patch", patch->mesh);

  /* Zero everything and give a silly mesh size to help debugging */
  memset(patch, 0, sizeof(struct pm_mesh_patch));
  patch->N = -1;
}
