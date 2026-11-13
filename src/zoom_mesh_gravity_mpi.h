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
#ifndef SWIFT_ZOOM_MESH_GRAVITY_MPI_H
#define SWIFT_ZOOM_MESH_GRAVITY_MPI_H

/* Config parameters. */
#include <config.h>

/* Forward declarations */
struct cell;
struct zoom_pm_mesh;
struct pm_mesh_patch;
struct neutrino_model;

void zoom_mesh_accumulate_cell_to_local_patch(
    const struct zoom_pm_mesh *mesh, const struct cell *cell,
    struct pm_mesh_patch *patch, const struct neutrino_model *nu_model);

#endif /* SWIFT_ZOOM_MESH_GRAVITY_MPI_H */
