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
 *  @param initial_partition The initial partition settings.
 *  @param nr_nodes The number of MPI ranks.
 *  @param s The #space containing the zoom and background top-level cells.
 */
void partition_zoom_grid(struct partition *initial_partition, int nr_nodes,
                         struct space *s) {

  /* Ensure we have a compatible grid size. */
  if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                      initial_partition->grid[2]) {
    error("Grid size does not match number of nodes.");
  }

  /* Apply the grid partitioning to the zoom top-level cells. */
  for (int cid = 0; cid < s->zoom_props->nr_zoom_cells; cid++) {

    /* Get the cell and compute its (i,j,k) index in the grid. */
    struct cell *c = &s->zoom_props->zoom_cells_top[cid];
    int i = (c->loc[0] - s->zoom_props->region_lower_bounds[0]) /
            s->zoom_props->dim[0] * initial_partition->grid[0];
    int j = (c->loc[1] - s->zoom_props->region_lower_bounds[1]) /
            s->zoom_props->dim[1] * initial_partition->grid[1];
    int k = (c->loc[2] - s->zoom_props->region_lower_bounds[2]) /
            s->zoom_props->dim[2] * initial_partition->grid[2];

    /* Compute the nodeID from the (i,j,k) index. */
    c->nodeID =
        i + initial_partition->grid[0] * (j + initial_partition->grid[1] * k);
  }

  /* Now the same for the background top-level cells. */
  for (int cid = 0; cid < s->zoom_props->nr_bkg_cells; cid++) {

    /* Get the cell and compute its (i,j,k) index in the grid. */
    struct cell *c = &s->zoom_props->bkg_cells_top[cid];
    int i = c->loc[0] / s->dim[0] * initial_partition->grid[0];
    int j = c->loc[1] / s->dim[1] * initial_partition->grid[1];
    int k = c->loc[2] / s->dim[2] * initial_partition->grid[2];

    /* Compute the nodeID from the (i,j,k) index. */
    c->nodeID =
        i + initial_partition->grid[0] * (j + initial_partition->grid[1] * k);
  }
}

/*  Vector support */
/*  ============== */

#ifdef WITH_MPI
/**
 * @brief Partition the space using a vectorised list of sample positions.
 *
 * The zoom and background top-level grids are vectorised independently using
 * their own grid geometry.
 *
 * @param nr_nodes The number of MPI ranks.
 * @param s The #space containing the zoom and background top-level cells.
 */
void partition_zoom_vector(int nr_nodes, struct space *s) {

  /* Allocate a buffer for the sample cell indices. We need 3 integers per node
   * to store the (i,j,k) index of the sample cell. */
  int *samplecells = NULL;
  if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
    error("Failed to allocate samplecells");

  /* First distribute the zoom top-level cells. */
  if (s->e->nodeID == 0) {
    pick_vector(s->zoom_props->cdim, nr_nodes, samplecells);
  }

  /* Broadcast the sample cell indices to all ranks. */
  int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample zoom cells.");

  /* Now split the zoom top-level cells according to the sample cell indices. */
  split_vector(s->zoom_props->zoom_cells_top, s->zoom_props->cdim, nr_nodes,
               samplecells);

  /* Clear the sample cell buffer ready for the background top-level cells. */
  bzero(samplecells, sizeof(int) * nr_nodes * 3);

  /* Then reuse the same buffer for the background top-level cells. */
  if (s->e->nodeID == 0) {
    pick_vector(s->cdim, nr_nodes, samplecells);
  }

  /* Broadcast the sample cell indices to all ranks. */
  res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample background cells.");

  /* Now split the background top-level cells according to the sample cell
   * indices. */
  split_vector(s->zoom_props->bkg_cells_top, s->cdim, nr_nodes, samplecells);

  free(samplecells);
}
#endif

/*  Void cell support */
/*  ================= */

/**
 *  @brief Partition the void cells.
 *
 *  For now we simply place all void cells on all ranks. In the future we may
 *  want to only assign void cells to be local if they have local zoom
 *  progeny, but this is a bit complex downstream for now.
 *
 *  @param s The #space containing the void top-level cells.
 *  @param nodeID The local MPI rank.
 */
void zoom_partition_voids(struct space *s, int nodeID) {

  for (int k = 0; k < s->zoom_props->nr_void_cells; k++) {
    s->cells_top[s->zoom_props->void_cell_indices[k]].nodeID = nodeID;
  }
}
