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
#include "partition.h"
#include "space.h"
#include "zoom.h"

/*  Grid support */
/*  ============ */

/**
 *  @brief Partition the space using a simple grid.
 *
 *  Using the sample positions as seeds pick cells that are geometrically
 *  closest and apply the partition to the space.
 *
 *  This is the function used for a volume with a zoom region.
 *
 *  Each cell grid will be partitioned onto the grid separately. This will
 *  not result in anything close to an optimal partitioning, but it will
 *  suffice given the relative cost between the cell grids.
 *
 *  @param initial_partition The initial partition data.
 *  @param nr_nodes The number of nodes.
 *  @param s The #space.
 */
void partition_zoom_grid(struct partition *initial_partition, int nr_nodes,
                         struct space *s) {

  /* If we've got the wrong number of nodes, fail. */
  if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                      initial_partition->grid[2])
    error("Grid size does not match number of nodes.");

  /* Run through the zoom cells and set their nodeID. */
  for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {

    /* Get the cell. */
    struct cell *c = &s->zoom_props->zoom_cells_top[k];

    /* Compute it's ijk grid index. */
    int ind[3];
    for (int j = 0; j < 3; j++) {
      ind[j] = (c->loc[j] - s->zoom_props->region_lower_bounds[j]) /
               s->zoom_props->dim[j] * initial_partition->grid[j];
    }

    /* Convert the grid index into the node ID. */
    c->nodeID = ind[0] + initial_partition->grid[0] *
                             (ind[1] + initial_partition->grid[1] * ind[2]);
  }

  /* Run through the background cells and set their nodeID. */
  for (int k = 0; k < s->zoom_props->nr_bkg_cells; k++) {

    /* Get the cell. */
    struct cell *c = &s->zoom_props->bkg_cells_top[k];

    /* Compute it's ijk grid index. */
    int ind[3];
    for (int j = 0; j < 3; j++) {
      ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
    }

    /* Convert the grid index into the node ID. */
    c->nodeID = ind[0] + initial_partition->grid[0] *
                             (ind[1] + initial_partition->grid[1] * ind[2]);
  }

  /* Run through the buffer cells and set their nodeID. (We don't always
   * have these) */
  for (int k = 0; k < s->zoom_props->nr_buffer_cells; k++) {

    /* Get the cell. */
    struct cell *c = &s->zoom_props->buffer_cells_top[k];

    /* Compute it's ijk grid index. */
    int ind[3];
    for (int j = 0; j < 3; j++) {
      ind[j] =
          (c->loc[j] - s->zoom_props->buffer_lower_bounds[j]) /
          (s->zoom_props->buffer_width[j] * s->zoom_props->buffer_cdim[j]) *
          initial_partition->grid[j];
    }

    /* Convert the grid index into the node ID. */
    c->nodeID = ind[0] + initial_partition->grid[0] *
                             (ind[1] + initial_partition->grid[1] * ind[2]);
  }
}

#ifdef WITH_MPI
/**
 * @brief Partition the space using a vectorised list of sample positions.
 *
 * Like the grid approach this will apply the vector independently to each
 * cell grid. This will not result in anything close to an optimal partitioning,
 * but it will suffice given the relative cost between the cell grids.
 *
 * @param nr_nodes The number of nodes.
 * @param s The #space.
 */
void partition_zoom_vector(int nr_nodes, struct space *s) {

  /* Vectorised selection, guaranteed to work for samples less than the
   * number of cells, but not very clumpy in the selection of regions. */
  int *samplecells = NULL;
  if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
    error("Failed to allocate samplecells");

  /* Do the zoom cells... */
  if (s->e->nodeID == 0) {
    pick_vector(s->zoom_props->cdim, nr_nodes, samplecells);
  }

  /* Share the samplecells around all the nodes. */
  int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample zoom cells.");

  /* And apply to our zoom cells */
  split_vector(s->zoom_props->zoom_cells_top, s->zoom_props->cdim, nr_nodes,
               samplecells);

  /* Clear out the samplecells. */
  bzero(samplecells, sizeof(int) * nr_nodes * 3);

  /* Do the buffer cells if we have any... */
  if (s->zoom_props->with_buffer_cells) {
    if (s->e->nodeID == 0) {
      pick_vector(s->zoom_props->buffer_cdim, nr_nodes, samplecells);
    }

    /* Share the samplecells around all the nodes. */
    res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the partition sample buffer cells.");

    /* And apply to our buffer cells */
    split_vector(s->zoom_props->buffer_cells_top, s->zoom_props->buffer_cdim,
                 nr_nodes, samplecells);

    /* Clear out the samplecells. */
    bzero(samplecells, sizeof(int) * nr_nodes * 3);
  }

  /* Finally do the background cells... */
  if (s->e->nodeID == 0) {
    pick_vector(s->cdim, nr_nodes, samplecells);
  }

  /* Share the samplecells around all the nodes. */
  res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the partition sample background cells.");

  /* And apply to our background cells */
  split_vector(s->zoom_props->bkg_cells_top, s->cdim, nr_nodes, samplecells);

  /* We're done */
  free(samplecells);
}
#endif

/*  Void cell support */
/*  ================ */

/**
 *  @brief Partition the void cells
 *
 *  If a void cell has one or more zoom cell leaves on this rank, then it
 *  is local.
 *
 *  @param s The #space.
 *  @param nodeID The node ID.
 */
void zoom_partition_voids(struct space *s, int nodeID) {

  /* Run through the zoom cells. */
  for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {

    /* Get the cell. */
    struct cell *c = &s->zoom_props->zoom_cells_top[k];

    /* Skip foreign cells. */
    if (c->nodeID != nodeID) continue;

    /* Find the top level void cell index. */
    int vid;
    if (s->zoom_props->with_buffer_cells) {
      vid = cell_getid_offset(
          s->zoom_props->buffer_cdim, s->zoom_props->buffer_cell_offset,
          c->loc[0] + (c->width[0] / 2), c->loc[1] + (c->width[1] / 2),
          c->loc[2] + (c->width[2] / 2));
    } else {
      vid = cell_getid_offset(s->cdim, s->zoom_props->bkg_cell_offset,
                              c->loc[0] + (c->width[0] / 2),
                              c->loc[1] + (c->width[1] / 2),
                              c->loc[2] + (c->width[2] / 2));
    }

    /* Get the void cell. */
    struct cell *vc = &s->cells_top[vid];

    /* Set the void cell's nodeID. */
    vc->nodeID = nodeID;
  }
}
