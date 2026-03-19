/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will Roper (w.roper@sussex.ac.uk)
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

/* Config */
#include <config.h>

/* Standard includes. */
#include <stdlib.h>
#include <string.h>

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "partition.h"
#include "space.h"
#include "zoom.h"

/*  Grid support */
/*  ============ */

/**
 *  @brief Partition the space using a simple grid.
 *
 *  Each cell grid is partitioned independently.
 */
void partition_zoom_grid(struct partition *initial_partition, int nr_nodes,
                         struct space *s) {

  if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                      initial_partition->grid[2])
    error("Grid size does not match number of nodes.");

  for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {
    struct cell *c = &s->zoom_props->zoom_cells_top[k];
    int ind[3];
    for (int j = 0; j < 3; j++) {
      ind[j] = (c->loc[j] - s->zoom_props->region_lower_bounds[j]) /
               s->zoom_props->dim[j] * initial_partition->grid[j];
    }
    c->nodeID = ind[0] + initial_partition->grid[0] *
                             (ind[1] + initial_partition->grid[1] * ind[2]);
  }

  for (int k = 0; k < s->zoom_props->nr_bkg_cells; k++) {
    struct cell *c = &s->zoom_props->bkg_cells_top[k];
    int ind[3];
    for (int j = 0; j < 3; j++) {
      ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
    }
    c->nodeID = ind[0] + initial_partition->grid[0] *
                             (ind[1] + initial_partition->grid[1] * ind[2]);
  }
}

/*  Vector support */
/*  ============== */

#ifdef WITH_MPI
/**
 * @brief Partition the space using a vectorised list of sample positions.
 */
void partition_zoom_vector(int nr_nodes, struct space *s) {

  int *samplecells = NULL;
  if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
    error("Failed to allocate samplecells");

  if (s->e->nodeID == 0) {
    pick_vector(s->zoom_props->cdim, nr_nodes, samplecells);
  }

  int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample zoom cells.");

  split_vector(s->zoom_props->zoom_cells_top, s->zoom_props->cdim, nr_nodes,
               samplecells);

  bzero(samplecells, sizeof(int) * nr_nodes * 3);

  if (s->e->nodeID == 0) {
    pick_vector(s->cdim, nr_nodes, samplecells);
  }

  res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample background cells.");

  split_vector(s->zoom_props->bkg_cells_top, s->cdim, nr_nodes, samplecells);

  free(samplecells);
}
#endif
