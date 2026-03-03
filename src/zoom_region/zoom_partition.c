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

/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif

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
}

/*  Vector support */
/*  ============ */

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

  /* All void cells are local for now. */
  for (int k = 0; k < s->zoom_props->nr_void_cells; k++) {
    s->cells_top[s->zoom_props->void_cell_indices[k]].nodeID = nodeID;
  }

  // /* Initialise the void cell nodeIDs to -1 (non-local). */
  // for (int k = 0; k < s->zoom_props->nr_void_cells; k++) {
  //   s->cells_top[s->zoom_props->void_cell_indices[k]].nodeID = -1;
  // }
  //
  // /* Run through the local zoom cells and mark their containing void cells as
  //  * local. */
  // for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {
  //
  //   /* Get the cell. */
  //   struct cell *c = &s->cells_top[k];
  //
  //   /* Skip foreign cells. */
  //   if (c->nodeID != nodeID) continue;
  //
  //   /* Find the top level void cell index. */
  //   double x = c->loc[0] + (c->width[0] / 2);
  //   double y = c->loc[1] + (c->width[1] / 2);
  //   double z = c->loc[2] + (c->width[2] / 2);
  //   int vid =
  //       cell_getid_offset(s->cdim, s->zoom_props->bkg_cell_offset, x, y, z);
  //
  //   /* Get the void cell. */
  //   struct cell *vc = &s->cells_top[vid];
  //
  //   /* Set the void cell's nodeID. */
  //   vc->nodeID = nodeID;
  // }
}

/*  Metis/Parmetis support */
/*  ====================== */

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/**
 * @brief Count the number of edges each cell has in the adjacency graph.
 *
 * An "edge" is a relationship between two cells that are adjacent and could
 * have interactions between them (i.e. gravity pair tasks etc.). For a
 * background cell, this is the number of immediate background cell neighbours
 * plus the zoom cells nested within any neighbouring void cells. Symmetrically,
 * for a zoom cell, this is the number of immediate zoom cell neighbours plus
 * the background cells that neighbour the void cell containing this zoom cell.
 *
 * We need to count these edges to allocate the adjacency list for
 * Metis/ParMetis.
 *
 * @param s The #space.
 * @param verbose Are we talking?
 */
int zoom_partition_count_vertex_edges(struct space *s, int periodic,
                                      int *cell_edge_offsets) {

  int *zoom_cdim = s->zoom_props->cdim;
  int *bkg_cdim = s->zoom_props->bkg_cdim;
  struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

  /* Temporary array to track edge counts per cell during counting */
  int *edge_counts = (int *)calloc(s->nr_cells, sizeof(int));
  if (edge_counts == NULL) error("Failed to allocate temporary edge counts");

  /* Loop over zoom cells and count the number of edges */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        /* Get the cell index. */
        int cid = cell_getid(zoom_cdim, i, j, k);

        /* Get the cell. */
        struct cell *ci = &zoom_cells[cid];

        /* Loop over the immediate neighbours in this grid, note that zoom
         * cells are never periodic. */
        for (int ii = -1; ii <= 1; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;
          for (int jj = -1; jj <= 1; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;
            for (int kk = -1; kk <= 1; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              /* Get the cell index. */
              int cjd = cell_getid(zoom_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We found an edge - count for both cells */
              edge_counts[cid]++;
              edge_counts[cjd]++;
            }
          }
        }

        /* Find the i,j,k position of this cell in the background. */
        int bkg_i = (ci->loc[0] + (ci->width[0] / 2)) * s->iwidth[0];
        int bkg_j = (ci->loc[1] + (ci->width[1] / 2)) * s->iwidth[1];
        int bkg_k = (ci->loc[2] + (ci->width[2] / 2)) * s->iwidth[2];

        /* Loop over the background neighbours. */
        for (int bkg_ii = -1; bkg_ii <= 1; bkg_ii++) {
          int bkg_iii = bkg_i + bkg_ii;
          if (!periodic && (bkg_iii < 0 || bkg_iii >= bkg_cdim[0])) continue;
          for (int bkg_jj = -1; bkg_jj <= 1; bkg_jj++) {
            int bkg_jjj = bkg_j + bkg_jj;
            if (!periodic && (bkg_jjj < 0 || bkg_jjj >= bkg_cdim[1])) continue;
            for (int bkg_kk = -1; bkg_kk <= 1; bkg_kk++) {
              int bkg_kkk = bkg_k + bkg_kk;
              if (!periodic && (bkg_kkk < 0 || bkg_kkk >= bkg_cdim[2]))
                continue;

              /* Skip the void cell containing this zoom cell */
              if (bkg_i == bkg_iii && bkg_j == bkg_jjj && bkg_k == bkg_kkk)
                continue;

              /* Get the cell index */
              int bkg_cjd = cell_getid(bkg_cdim, bkg_iii, bkg_jjj, bkg_kkk);

              /* We have a background neighbour - count for both cells */
              edge_counts[cid]++;
              edge_counts[bkg_cjd]++;
            }
          }
        }
      }
    }
  }

  /* Loop over all background cells and count the number of edges. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the cell. */
        int cid = cell_getid(bkg_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid (we've done the zoom
         * edges above) */
        for (int ii = -1; ii <= 1; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= bkg_cdim[0])) continue;
          iii = (iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int jj = -1; jj <= 1; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= bkg_cdim[1])) continue;
            jjj = (jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int kk = -1; kk <= 1; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= bkg_cdim[2])) continue;
              kkk = (kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Get the cell index */
              int cjd = cell_getid(bkg_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We have a neighbour - count for both cells */
              edge_counts[cid]++;
              edge_counts[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build cumulative offsets from edge counts */
  cell_edge_offsets[0] = 0;
  for (int cid = 0; cid < s->nr_cells; cid++) {
    cell_edge_offsets[cid + 1] = cell_edge_offsets[cid] + edge_counts[cid];
  }

  /* Clean up */
  free(edge_counts);

  /* Return total number of edges */
  return cell_edge_offsets[s->nr_cells];
}

void zoom_partition_graph_init(struct space *s, int periodic, idx_t *adjncy,
                               int *nadjcny, idx_t *xadj, int *nxadj,
                               const int *cell_edge_offsets) {

  int *zoom_cdim = s->zoom_props->cdim;
  int *bkg_cdim = s->zoom_props->bkg_cdim;
  struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

  /* Track current position in adjncy for each cell using temporary array */
  int *edge_pos = (int *)malloc(s->nr_cells * sizeof(int));
  if (edge_pos == NULL) error("Failed to allocate temporary edge positions");

  /* Initialize positions from offsets */
  for (int i = 0; i < s->nr_cells; i++) {
    edge_pos[i] = cell_edge_offsets[i];
  }

  /* Loop over zoom cells and populate adjacency array */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        /* Get the cell index. */
        int ci_id = cell_getid(zoom_cdim, i, j, k);

        /* Get the cell. */
        struct cell *ci = &zoom_cells[ci_id];

        /* Loop over the immediate neighbours in this grid, note that zoom
         * cells are never periodic. */
        for (int ii = -1; ii <= 1; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;
          for (int jj = -1; jj <= 1; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;
            for (int kk = -1; kk <= 1; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              /* Get the cell index. */
              int cjd = cell_getid(zoom_cdim, iii, jjj, kkk);

              /* Avoid duplicates - same logic as counting */
              if (ci_id >= cjd) continue;

              /* We found an edge - record in both directions */
              if (adjncy != NULL) {
                adjncy[edge_pos[ci_id]] = cjd;
                adjncy[edge_pos[cjd]] = ci_id;
              }
              edge_pos[ci_id]++;
              edge_pos[cjd]++;
            }
          }
        }

        /* Find the i,j,k position of this cell in the background. */
        int bkg_i = (ci->loc[0] + (ci->width[0] / 2)) * s->iwidth[0];
        int bkg_j = (ci->loc[1] + (ci->width[1] / 2)) * s->iwidth[1];
        int bkg_k = (ci->loc[2] + (ci->width[2] / 2)) * s->iwidth[2];

        /* Loop over the background neighbours. */
        for (int bkg_ii = -1; bkg_ii <= 1; bkg_ii++) {
          int bkg_iii = bkg_i + bkg_ii;
          if (!periodic && (bkg_iii < 0 || bkg_iii >= bkg_cdim[0])) continue;
          for (int bkg_jj = -1; bkg_jj <= 1; bkg_jj++) {
            int bkg_jjj = bkg_j + bkg_jj;
            if (!periodic && (bkg_jjj < 0 || bkg_jjj >= bkg_cdim[1])) continue;
            for (int bkg_kk = -1; bkg_kk <= 1; bkg_kk++) {
              int bkg_kkk = bkg_k + bkg_kk;
              if (!periodic && (bkg_kkk < 0 || bkg_kkk >= bkg_cdim[2]))
                continue;

              /* Skip the void cell containing this zoom cell */
              if (bkg_i == bkg_iii && bkg_j == bkg_jjj && bkg_k == bkg_kkk)
                continue;

              /* Get the cell index */
              int bkg_cjd = cell_getid(bkg_cdim, bkg_iii, bkg_jjj, bkg_kkk);

              /* We have a background neighbour - record in both directions */
              if (adjncy != NULL) {
                adjncy[edge_pos[ci_id]] = bkg_cjd;
                adjncy[edge_pos[bkg_cjd]] = ci_id;
              }
              edge_pos[ci_id]++;
              edge_pos[bkg_cjd]++;
            }
          }
        }
      }
    }
  }

  /* Loop over all background cells and populate adjacency array. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the cell. */
        int ci_id = cell_getid(bkg_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid (we've done the zoom
         * edges above) */
        for (int ii = -1; ii <= 1; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= bkg_cdim[0])) continue;
          iii = (iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int jj = -1; jj <= 1; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= bkg_cdim[1])) continue;
            jjj = (jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int kk = -1; kk <= 1; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= bkg_cdim[2])) continue;
              kkk = (kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Get the cell index */
              int cjd = cell_getid(bkg_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (ci_id >= cjd) continue;

              /* We have a neighbour - record in both directions */
              if (adjncy != NULL) {
                adjncy[edge_pos[ci_id]] = cjd;
                adjncy[edge_pos[cjd]] = ci_id;
              }
              edge_pos[ci_id]++;
              edge_pos[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build xadj array from cell_edge_offsets if needed */
  if (xadj != NULL) {
    for (int i = 0; i <= s->nr_cells; i++) {
      xadj[i] = cell_edge_offsets[i];
    }
  }

  /* Clean up */
  free(edge_pos);

  /* Set the output counters */
  if (nadjcny != NULL) *nadjcny = cell_edge_offsets[s->nr_cells];
  if (nxadj != NULL) *nxadj = s->nr_cells;
}
#endif
