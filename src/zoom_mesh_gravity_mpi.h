/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will J. Roper (w.roper@sussex.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *******************************************************************************/
#ifndef SWIFT_ZOOM_MESH_GRAVITY_MPI_H
#define SWIFT_ZOOM_MESH_GRAVITY_MPI_H

struct space;
struct threadpool;
struct zoom_pm_mesh;

void zoom_mesh_compute_potential_distributed(struct zoom_pm_mesh *mesh,
                                             const struct space *s,
                                             struct threadpool *tp,
                                             int verbose);

#endif /* SWIFT_ZOOM_MESH_GRAVITY_MPI_H */
