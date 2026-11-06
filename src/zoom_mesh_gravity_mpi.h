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

struct space;
struct threadpool;
struct zoom_pm_mesh;

/**
 * @brief Compute the isolated zoom correction with a distributed FFT.
 *
 * @param mesh The #zoom_pm_mesh with its cached Fourier kernel.
 * @param s The #space containing locally owned gravity particles.
 * @param tp The #threadpool used for deposition and interpolation.
 * @param verbose Are we talkative?
 */
void zoom_mesh_compute_potential_distributed(struct zoom_pm_mesh *mesh,
                                             const struct space *s,
                                             struct threadpool *tp,
                                             int verbose);

#endif /* SWIFT_ZOOM_MESH_GRAVITY_MPI_H */
