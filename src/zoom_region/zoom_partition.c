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
#include "clocks.h"
#include "engine.h"
#include "error.h"
#include "memuse.h"
#include "partition.h"
#include "sort_part.h"
#include "space.h"
#include "tools.h"
#include "zoom.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

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

/*  Metis/Parmetis support */
/*  ====================== */

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/**
 * @brief Return the directional edge weight for a neighbour offset.
 *
 * The uniform-grid implementation normalised these weights by the fixed
 * 26-neighbour stencil. In the zoom graph the emitted CSR row length is no
 * longer fixed, so we normalise by the actual number of outgoing edges for each
 * cell to avoid over-weighting interface cells.
 *
 * @param count The number of particles in the source cell.
 * @param cid The source cell index.
 * @param cell_edge_offsets The cumulative edge offsets for each cell, i.e.
 * cell_edge_offsets[cid+1] - cell_edge_offsets[cid] gives the number of
 * outgoing edges for this cell.
 * @param di The i-offset of the neighbour within the local adjacent cells
 * around cid (-1, 0, or 1).
 * @param dj The j-offset of the neighbour within the local adjacent cells
 * around cid (-1, 0, or 1).
 * @param dk The k-offset of the neighbour within the local adjacent cells
 * around cid (-1, 0, or 1).
 *
 * @return The edge weight for the neighbour at the given offset, normalised by
 * the number of outgoing edges for the source cell.
 */
__attribute__((always_inline)) INLINE static double zoom_partition_edge_weight(
    const double count, const int cid, const int *cell_edge_offsets,
    const int di, const int dj, const int dk) {
  const int sid = ((di + 1) * 9 + (dj + 1) * 3 + (dk + 1));
  const int nr_edges = cell_edge_offsets[cid + 1] - cell_edge_offsets[cid];
  return count * sid_scale[sortlistID[sid]] / nr_edges;
}

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
 * @param periodic Whether the background grid is periodic.
 * @param cell_edge_offsets Output array of cumulative edge offsets.
 */
int zoom_partition_count_vertex_edges(struct space *s, int periodic,
                                      int *cell_edge_offsets) {

  int *zoom_cdim = s->zoom_props->cdim;
  int *bkg_cdim = s->zoom_props->bkg_cdim;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
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
          if (periodic) bkg_iii = (bkg_iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int bkg_jj = -1; bkg_jj <= 1; bkg_jj++) {
            int bkg_jjj = bkg_j + bkg_jj;
            if (!periodic && (bkg_jjj < 0 || bkg_jjj >= bkg_cdim[1])) continue;
            if (periodic) bkg_jjj = (bkg_jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int bkg_kk = -1; bkg_kk <= 1; bkg_kk++) {
              int bkg_kkk = bkg_k + bkg_kk;
              if (!periodic && (bkg_kkk < 0 || bkg_kkk >= bkg_cdim[2]))
                continue;
              if (periodic) bkg_kkk = (bkg_kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Skip the void cell containing this zoom cell */
              if (bkg_i == bkg_iii && bkg_j == bkg_jjj && bkg_k == bkg_kkk)
                continue;

              /* Get the cell index */
              int bkg_cjd = cell_getid_offset(bkg_cdim, bkg_offset, bkg_iii,
                                              bkg_jjj, bkg_kkk);

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
        int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);

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
              int cjd = cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);

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

/**
 * @brief Make edge weights from the accumulated particle sizes per cell for the
 * zoom adjacency graph.
 *
 * This preserves the directional weighting used by the legacy 26-neighbour
 * graph, but applies it to the zoom-specific CSR graph by using the geometric
 * offset associated with each zoom-zoom, zoom-background, and background-
 * background edge.
 *
 * @param s the space containing the cells.
 * @param counts the number of bytes in particles per cell.
 * @param edges weights for the graph edges in CSR order.
 * @param cell_edge_offsets array[nr_cells+1] with cumulative edge offsets.
 */
void zoom_partition_sizes_to_edges(struct space *s, double *counts,
                                   double *edges,
                                   const int *cell_edge_offsets) {

  int *zoom_cdim = s->zoom_props->cdim;
  int *bkg_cdim = s->zoom_props->bkg_cdim;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
  struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

  /* Clear the output edge weights to be repopulated */
  const int nedges = cell_edge_offsets[s->nr_cells];
  bzero(edges, sizeof(double) * nedges);

  /* Temporary array to track the current edge position for each cell. */
  int *edge_pos = (int *)malloc(s->nr_cells * sizeof(int));
  if (edge_pos == NULL) error("Failed to allocate temporary edge positions");

  /* Initialise the edge positions from the cumulative edge offsets. */
  for (int i = 0; i < s->nr_cells; i++) {
    edge_pos[i] = cell_edge_offsets[i];
  }

  /* Loop over zoom cells and assign the edge weights. Note that zoom cells
   * are never periodic */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        /* Get the cell. */
        const int cid = cell_getid(zoom_cdim, i, j, k);
        struct cell *ci = &zoom_cells[cid];
        /* Loop over the immediate neighbours in the zoom grid. */
        for (int ii = -1; ii <= 1; ii++) {
          const int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;
          for (int jj = -1; jj <= 1; jj++) {
            const int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;
            for (int kk = -1; kk <= 1; kk++) {
              const int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              /* Get the neighbour index */
              const int cjd = cell_getid(zoom_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We have a neighbour: symmetrically assign the weights. */
              edges[edge_pos[cid]] = zoom_partition_edge_weight(
                  counts[cid], cid, cell_edge_offsets, ii, jj, kk);
              edges[edge_pos[cjd]] = zoom_partition_edge_weight(
                  counts[cjd], cjd, cell_edge_offsets, -ii, -jj, -kk);
              edge_pos[cid]++;
              edge_pos[cjd]++;
            }
          }
        }

        /* Get the top level bkg cell indices */
        const int bkg_i = (ci->loc[0] + (ci->width[0] / 2)) * s->iwidth[0];
        const int bkg_j = (ci->loc[1] + (ci->width[1] / 2)) * s->iwidth[1];
        const int bkg_k = (ci->loc[2] + (ci->width[2] / 2)) * s->iwidth[2];

        /* Loop over the neighbouring background cells. */
        for (int bkg_ii = -1; bkg_ii <= 1; bkg_ii++) {
          int bkg_iii = bkg_i + bkg_ii;
          if (!s->periodic && (bkg_iii < 0 || bkg_iii >= bkg_cdim[0])) continue;
          if (s->periodic) bkg_iii = (bkg_iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int bkg_jj = -1; bkg_jj <= 1; bkg_jj++) {
            int bkg_jjj = bkg_j + bkg_jj;
            if (!s->periodic && (bkg_jjj < 0 || bkg_jjj >= bkg_cdim[1]))
              continue;
            if (s->periodic) bkg_jjj = (bkg_jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int bkg_kk = -1; bkg_kk <= 1; bkg_kk++) {
              int bkg_kkk = bkg_k + bkg_kk;
              if (!s->periodic && (bkg_kkk < 0 || bkg_kkk >= bkg_cdim[2]))
                continue;
              if (s->periodic) bkg_kkk = (bkg_kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Skip the void cell containing this zoom cell. */
              if (bkg_i == bkg_iii && bkg_j == bkg_jjj && bkg_k == bkg_kkk)
                continue;

              /* Get the bkg neighbour index */
              const int bkg_cjd = cell_getid_offset(bkg_cdim, bkg_offset,
                                                    bkg_iii, bkg_jjj, bkg_kkk);

              /* We have a neighbour: symmetrically assign the weights. */
              edges[edge_pos[cid]] = zoom_partition_edge_weight(
                  counts[cid], cid, cell_edge_offsets, bkg_ii, bkg_jj, bkg_kk);
              edges[edge_pos[bkg_cjd]] = zoom_partition_edge_weight(
                  counts[bkg_cjd], bkg_cjd, cell_edge_offsets, -bkg_ii, -bkg_jj,
                  -bkg_kk);
              edge_pos[cid]++;
              edge_pos[bkg_cjd]++;
            }
          }
        }
      }
    }
  }

  /* Loop over background cells and assign the remaining edge weights. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the cell index */
        const int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);
        /* Loop over the immediate neighbours in the background grid. */
        for (int ii = -1; ii <= 1; ii++) {
          int iii = i + ii;
          if (!s->periodic && (iii < 0 || iii >= bkg_cdim[0])) continue;
          iii = (iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int jj = -1; jj <= 1; jj++) {
            int jjj = j + jj;
            if (!s->periodic && (jjj < 0 || jjj >= bkg_cdim[1])) continue;
            jjj = (jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int kk = -1; kk <= 1; kk++) {
              int kkk = k + kk;
              if (!s->periodic && (kkk < 0 || kkk >= bkg_cdim[2])) continue;
              kkk = (kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Get the neighbour index */
              const int cjd =
                  cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We have a neighbour: symmetrically assign the weights. */
              edges[edge_pos[cid]] = zoom_partition_edge_weight(
                  counts[cid], cid, cell_edge_offsets, ii, jj, kk);
              edges[edge_pos[cjd]] = zoom_partition_edge_weight(
                  counts[cjd], cjd, cell_edge_offsets, -ii, -jj, -kk);
              edge_pos[cid]++;
              edge_pos[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Clean up. */
  free(edge_pos);
}

/**
 * @brief Fill the CSR adjacency arrays for the zoom-specific cell graph.
 *
 * Builds the combined zoom/background adjacency structure using the same edge
 * ordering as zoom_partition_count_vertex_edges() and
 * zoom_partition_sizes_to_edges().
 *
 * @param s The #space.
 * @param periodic Whether the background grid is periodic.
 * @param adjncy The adjacency array to fill.
 * @param nadjcny The number of adjacency entries written.
 * @param xadj The CSR row-offset array to fill.
 * @param nxadj The number of xadj entries written.
 * @param cell_edge_offsets Input array of cumulative edge offsets.
 */
void zoom_partition_graph_init(struct space *s, int periodic, idx_t *adjncy,
                               int *nadjcny, idx_t *xadj, int *nxadj,
                               const int *cell_edge_offsets) {

  int *zoom_cdim = s->zoom_props->cdim;
  int *bkg_cdim = s->zoom_props->bkg_cdim;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
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
          if (periodic) bkg_iii = (bkg_iii + bkg_cdim[0]) % bkg_cdim[0];
          for (int bkg_jj = -1; bkg_jj <= 1; bkg_jj++) {
            int bkg_jjj = bkg_j + bkg_jj;
            if (!periodic && (bkg_jjj < 0 || bkg_jjj >= bkg_cdim[1])) continue;
            if (periodic) bkg_jjj = (bkg_jjj + bkg_cdim[1]) % bkg_cdim[1];
            for (int bkg_kk = -1; bkg_kk <= 1; bkg_kk++) {
              int bkg_kkk = bkg_k + bkg_kk;
              if (!periodic && (bkg_kkk < 0 || bkg_kkk >= bkg_cdim[2]))
                continue;
              if (periodic) bkg_kkk = (bkg_kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Skip the void cell containing this zoom cell */
              if (bkg_i == bkg_iii && bkg_j == bkg_jjj && bkg_k == bkg_kkk)
                continue;

              /* Get the cell index */
              int bkg_cjd = cell_getid_offset(bkg_cdim, bkg_offset, bkg_iii,
                                              bkg_jjj, bkg_kkk);

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
        int ci_id = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);

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
              int cjd = cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);

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

/*  Separate zoom/background subgraph partitioning */
/*  ============================================== */

/**
 * @brief Edge weight for a subgraph edge (zoom or background).
 *
 * Mirrors zoom_partition_edge_weight() but takes the subgraph-local cumulative
 * edge offsets so it can be used while building each subgraph in isolation.
 *
 * @param count The number of particles in the source cell.
 * @param sub_cid The source cell index within the subgraph.
 * @param sub_offsets Subgraph cumulative edge offsets array.
 * @param di The i-offset of the neighbour (-1, 0, or 1).
 * @param dj The j-offset of the neighbour (-1, 0, or 1).
 * @param dk The k-offset of the neighbour (-1, 0, or 1).
 *
 * @return The edge weight, normalised by the number of outgoing edges.
 */
__attribute__((always_inline)) INLINE static double
zoom_partition_subgraph_edge_weight(const double count, const int sub_cid,
                                    const int *sub_offsets, const int di,
                                    const int dj, const int dk) {
  const int sid = ((di + 1) * 9 + (dj + 1) * 3 + (dk + 1));
  const int nr_edges = sub_offsets[sub_cid + 1] - sub_offsets[sub_cid];
  return count * sid_scale[sortlistID[sid]] / nr_edges;
}

/**
 * @brief Count the edges of the zoom-only subgraph.
 *
 * Considers only zoom-zoom adjacencies in the zoom top-level grid. The zoom
 * grid is never periodic.
 *
 * @param s The #space.
 * @param sub_offsets Output array[nr_zoom_cells+1] of cumulative edge offsets.
 *
 * @return Total number of edges in the zoom subgraph.
 */
static int zoom_partition_count_zoom_subgraph_edges(struct space *s,
                                                    int *sub_offsets) {

  const int *zoom_cdim = s->zoom_props->cdim;
  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Temporary array to track edge counts per cell during counting. */
  int *edge_counts = (int *)calloc(nr_zoom_cells, sizeof(int));
  if (edge_counts == NULL)
    error("Failed to allocate temporary zoom subgraph edge counts");

  /* Loop over zoom cells and count zoom-zoom edges. */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        /* Get the cell index. */
        const int cid = cell_getid(zoom_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid, note that zoom
         * cells are never periodic. */
        for (int ii = -1; ii <= 1; ii++) {
          const int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;
          for (int jj = -1; jj <= 1; jj++) {
            const int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;
            for (int kk = -1; kk <= 1; kk++) {
              const int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              /* Get the neighbour index. */
              const int cjd = cell_getid(zoom_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We found an edge - count for both cells. */
              edge_counts[cid]++;
              edge_counts[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build cumulative offsets from edge counts. */
  sub_offsets[0] = 0;
  for (int cid = 0; cid < nr_zoom_cells; cid++) {
    sub_offsets[cid + 1] = sub_offsets[cid] + edge_counts[cid];
  }

  /* Clean up. */
  free(edge_counts);

  return sub_offsets[nr_zoom_cells];
}

/**
 * @brief Count the edges of the background-only subgraph.
 *
 * Considers only background-background adjacencies in the background top-level
 * grid. Voids are part of the background grid and are included in the
 * subgraph; they will be reassigned to ranks post-partition by
 * zoom_partition_voids().
 *
 * @param s The #space.
 * @param periodic Whether the background grid is periodic.
 * @param sub_offsets Output array[nr_bkg_cells+1] of cumulative edge offsets.
 *
 * @return Total number of edges in the background subgraph.
 */
static int zoom_partition_count_bkg_subgraph_edges(struct space *s, int periodic,
                                                   int *sub_offsets) {

  const int *bkg_cdim = s->zoom_props->bkg_cdim;
  const int nr_bkg_cells = s->zoom_props->nr_bkg_cells;

  /* Temporary array to track edge counts per cell during counting. */
  int *edge_counts = (int *)calloc(nr_bkg_cells, sizeof(int));
  if (edge_counts == NULL)
    error("Failed to allocate temporary background subgraph edge counts");

  /* Loop over background cells and count bkg-bkg edges. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the local (subgraph) cell index. */
        const int cid = cell_getid(bkg_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid. */
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

              /* Get the neighbour index. */
              const int cjd = cell_getid(bkg_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We found an edge - count for both cells. */
              edge_counts[cid]++;
              edge_counts[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build cumulative offsets from edge counts. */
  sub_offsets[0] = 0;
  for (int cid = 0; cid < nr_bkg_cells; cid++) {
    sub_offsets[cid + 1] = sub_offsets[cid] + edge_counts[cid];
  }

  /* Clean up. */
  free(edge_counts);

  return sub_offsets[nr_bkg_cells];
}

/**
 * @brief Fill the CSR adjacency arrays and edge weights for the zoom subgraph.
 *
 * Builds xadj and adjncy for the zoom-only subgraph and, when edge weights are
 * requested, fills the per-edge weights derived from the per-cell vertex
 * weights using the same directional weighting scheme as the combined-graph
 * implementation.
 *
 * @param s The #space.
 * @param vertexw Per-cell vertex weights indexed across the full space, or
 *        NULL if no edge weights are needed.
 * @param adjncy Output adjacency array (subgraph CSR).
 * @param xadj Output row-offset array of size nr_zoom_cells+1.
 * @param edgew Output edge weights array, or NULL if vertexw is NULL.
 * @param sub_offsets Cumulative edge offsets from
 *        zoom_partition_count_zoom_subgraph_edges().
 */
static void zoom_partition_zoom_subgraph_init(struct space *s, double *vertexw,
                                              idx_t *adjncy, idx_t *xadj,
                                              double *edgew,
                                              const int *sub_offsets) {

  const int *zoom_cdim = s->zoom_props->cdim;
  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Track current adjncy position for each cell. */
  int *edge_pos = (int *)malloc(nr_zoom_cells * sizeof(int));
  if (edge_pos == NULL)
    error("Failed to allocate temporary zoom subgraph edge positions");

  /* Initialise positions from the cumulative offsets. */
  for (int i = 0; i < nr_zoom_cells; i++) edge_pos[i] = sub_offsets[i];

  /* Clear the output edge weights to be repopulated. */
  if (edgew != NULL) bzero(edgew, sizeof(double) * sub_offsets[nr_zoom_cells]);

  /* Loop over zoom cells and populate adjacency + edge weights. */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        /* Get the cell index. */
        const int cid = cell_getid(zoom_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid (never periodic). */
        for (int ii = -1; ii <= 1; ii++) {
          const int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;
          for (int jj = -1; jj <= 1; jj++) {
            const int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;
            for (int kk = -1; kk <= 1; kk++) {
              const int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              /* Get the neighbour index. */
              const int cjd = cell_getid(zoom_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We found an edge - record both directions. */
              if (adjncy != NULL) {
                adjncy[edge_pos[cid]] = cjd;
                adjncy[edge_pos[cjd]] = cid;
              }

              /* Symmetrically assign the edge weights. The vertex weights are
               * indexed in the full-space frame, but for the zoom subgraph
               * the zoom-cell index equals its full-space index. */
              if (edgew != NULL && vertexw != NULL) {
                edgew[edge_pos[cid]] = zoom_partition_subgraph_edge_weight(
                    vertexw[cid], cid, sub_offsets, ii, jj, kk);
                edgew[edge_pos[cjd]] = zoom_partition_subgraph_edge_weight(
                    vertexw[cjd], cjd, sub_offsets, -ii, -jj, -kk);
              }

              edge_pos[cid]++;
              edge_pos[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build xadj from sub_offsets. */
  if (xadj != NULL) {
    for (int i = 0; i <= nr_zoom_cells; i++) xadj[i] = sub_offsets[i];
  }

  /* Clean up. */
  free(edge_pos);
}

/**
 * @brief Fill the CSR adjacency arrays and edge weights for the bkg subgraph.
 *
 * As zoom_partition_zoom_subgraph_init() but for the background subgraph.
 * Vertex weights are indexed in the full-space frame and are accessed using
 * the bkg_cell_offset to convert subgraph indices to full-space indices.
 *
 * @param s The #space.
 * @param periodic Whether the background grid is periodic.
 * @param vertexw Per-cell vertex weights indexed across the full space, or
 *        NULL if no edge weights are needed.
 * @param adjncy Output adjacency array (subgraph CSR).
 * @param xadj Output row-offset array of size nr_bkg_cells+1.
 * @param edgew Output edge weights array, or NULL if vertexw is NULL.
 * @param sub_offsets Cumulative edge offsets from
 *        zoom_partition_count_bkg_subgraph_edges().
 */
static void zoom_partition_bkg_subgraph_init(struct space *s, int periodic,
                                             double *vertexw, idx_t *adjncy,
                                             idx_t *xadj, double *edgew,
                                             const int *sub_offsets) {

  const int *bkg_cdim = s->zoom_props->bkg_cdim;
  const int nr_bkg_cells = s->zoom_props->nr_bkg_cells;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;

  /* Track current adjncy position for each cell. */
  int *edge_pos = (int *)malloc(nr_bkg_cells * sizeof(int));
  if (edge_pos == NULL)
    error("Failed to allocate temporary bkg subgraph edge positions");

  /* Initialise positions from the cumulative offsets. */
  for (int i = 0; i < nr_bkg_cells; i++) edge_pos[i] = sub_offsets[i];

  /* Clear the output edge weights to be repopulated. */
  if (edgew != NULL) bzero(edgew, sizeof(double) * sub_offsets[nr_bkg_cells]);

  /* Loop over background cells and populate adjacency + edge weights. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the local (subgraph) cell index. */
        const int cid = cell_getid(bkg_cdim, i, j, k);

        /* Loop over the immediate neighbours in this grid. */
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

              /* Get the neighbour index. */
              const int cjd = cell_getid(bkg_cdim, iii, jjj, kkk);

              /* Avoid duplicates. */
              if (cid >= cjd) continue;

              /* We found an edge - record both directions. */
              if (adjncy != NULL) {
                adjncy[edge_pos[cid]] = cjd;
                adjncy[edge_pos[cjd]] = cid;
              }

              /* Symmetrically assign the edge weights. Vertex weights are in
               * the full-space frame so we offset to find the bkg entry. */
              if (edgew != NULL && vertexw != NULL) {
                edgew[edge_pos[cid]] = zoom_partition_subgraph_edge_weight(
                    vertexw[cid + bkg_offset], cid, sub_offsets, ii, jj, kk);
                edgew[edge_pos[cjd]] = zoom_partition_subgraph_edge_weight(
                    vertexw[cjd + bkg_offset], cjd, sub_offsets, -ii, -jj, -kk);
              }

              edge_pos[cid]++;
              edge_pos[cjd]++;
            }
          }
        }
      }
    }
  }

  /* Build xadj from sub_offsets. */
  if (xadj != NULL) {
    for (int i = 0; i <= nr_bkg_cells; i++) xadj[i] = sub_offsets[i];
  }

  /* Clean up. */
  free(edge_pos);
}

/**
 * @brief Run METIS on a single subgraph and write its partitioning.
 *
 * This is the core METIS driver shared by zoom_partition_pick_metis() between
 * the zoom and background subgraph passes. Only rank 0 computes the
 * partitioning; the result is broadcast to all ranks. Mirrors the structure
 * of partition_pick_metis() but operates on a pre-built subgraph CSR.
 *
 * @param nodeID The local MPI rank.
 * @param nregions The number of partitions to produce.
 * @param nverts The number of vertices in the subgraph.
 * @param nedges The number of edges in the subgraph.
 * @param vertexw Per-vertex weights in the subgraph index space, or NULL.
 * @param edgew Per-edge weights in CSR order, or NULL.
 * @param xadj Subgraph CSR row offsets, size nverts + 1.
 * @param adjncy Subgraph CSR adjacency, size nedges.
 * @param sub_part Output partition array, size nverts.
 */
static void zoom_partition_metis_subgraph(int nodeID, int nregions, int nverts,
                                          int nedges, double *vertexw,
                                          double *edgew, idx_t *xadj,
                                          idx_t *adjncy, int *sub_part) {

  /* Single-rank shortcut. */
  if (nregions == 1) {
    for (int i = 0; i < nverts; i++) sub_part[i] = 0;
    return;
  }

  /* Only rank 0 runs METIS, the result is then broadcast. */
  if (nodeID == 0) {

    idx_t *weights_v = NULL;
    if (vertexw != NULL) {
      if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
        error("Failed to allocate subgraph vertex weights array");
      for (int k = 0; k < nverts; k++) {
        weights_v[k] = (vertexw[k] > 1) ? (idx_t)vertexw[k] : 0;
      }
    }

    idx_t *weights_e = NULL;
    if (edgew != NULL) {
      if ((weights_e = (idx_t *)swift_malloc("partition_edge_weights",
                                             sizeof(idx_t) * nedges)) == NULL)
        error("Failed to allocate subgraph edge weights array");
      for (int k = 0; k < nedges; k++) {
        weights_e[k] = (edgew[k] > 1) ? (idx_t)edgew[k] : 1;
      }
    }

    idx_t *regionid;
    if ((regionid = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate subgraph regionid array");

    /* Set the METIS options. */
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_NCUTS] = 10;
    options[METIS_OPTION_NITER] = 20;

    /* Call METIS. */
    idx_t one = 1;
    idx_t idx_nverts = nverts;
    idx_t idx_nregions = nregions;
    idx_t objval;

    if (METIS_PartGraphKway(&idx_nverts, &one, xadj, adjncy, weights_v, NULL,
                            weights_e, &idx_nregions, NULL, NULL, options,
                            &objval, regionid) != METIS_OK)
      error("Call to METIS_PartGraphKway failed.");

    /* Copy back, checking for sanity. */
    for (int k = 0; k < nverts; k++) {
      if (regionid[k] < 0 || regionid[k] >= nregions)
        error("Got bad nodeID %" PRIDX " for subgraph cell %i.", regionid[k], k);
      sub_part[k] = regionid[k];
    }

    /* Clean up. */
    if (weights_v != NULL) free(weights_v);
    if (weights_e != NULL) swift_free("partition_edge_weights", weights_e);
    free(regionid);
  }

  /* Calculations done, everyone gets a copy. */
  int res = MPI_Bcast(sub_part, nverts, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS)
    mpi_error(res, "Failed to broadcast subgraph celllist");
}

/**
 * @brief Partition the zoom and background grids independently with METIS.
 *
 * The zoom and background top-level grids are partitioned as two independent
 * subgraphs. Only intra-grid edges are considered; cross-grid (zoom-bkg)
 * coupling is dropped, which decouples the optimal balance of the two grids
 * and prevents the combined-graph partitioner from concentrating the dense
 * zoom region on a single rank.
 *
 * Voids are bkg cells in the index space and are partitioned with the
 * background subgraph; their final assignment is overwritten by
 * zoom_partition_voids() during partition_split_metis().
 *
 * @param nodeID The local MPI rank.
 * @param s The #space containing the zoom and background top-level cells.
 * @param nregions The number of partitions to produce.
 * @param vertexw Per-cell vertex weights, sizeof number of cells, or NULL.
 * @param edgew Combined-graph edge weights from the caller. Ignored by the
 *        separate-subgraph path; kept in the signature so the calling
 *        convention matches partition_pick_metis().
 * @param celllist Output partition array, size s->nr_cells, indexed in the
 *        full-space frame.
 */
void zoom_partition_pick_metis(int nodeID, struct space *s, int nregions,
                               double *vertexw, double *edgew, int *celllist) {

  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;
  const int nr_bkg_cells = s->zoom_props->nr_bkg_cells;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;

  /* The combined-graph edge weights (if any) are not used by the separate
   * subgraph path; we recompute per-subgraph edge weights from vertexw. */
  (void)edgew;

  /* --- Zoom subgraph --- */

  int *zoom_sub_offsets =
      (int *)malloc(sizeof(int) * (nr_zoom_cells + 1));
  if (zoom_sub_offsets == NULL)
    error("Failed to allocate zoom subgraph offsets");

  const int zoom_nedges =
      zoom_partition_count_zoom_subgraph_edges(s, zoom_sub_offsets);

  idx_t *zoom_xadj =
      (idx_t *)malloc(sizeof(idx_t) * (nr_zoom_cells + 1));
  if (zoom_xadj == NULL) error("Failed to allocate zoom subgraph xadj");

  idx_t *zoom_adjncy = (idx_t *)swift_malloc(
      "partition_adjncy", sizeof(idx_t) * zoom_nedges);
  if (zoom_adjncy == NULL) error("Failed to allocate zoom subgraph adjncy");

  double *zoom_edgew = NULL;
  if (vertexw != NULL) {
    if ((zoom_edgew = (double *)swift_malloc(
             "partition_edge_weights", sizeof(double) * zoom_nedges)) == NULL)
      error("Failed to allocate zoom subgraph edge weights");
  }

  zoom_partition_zoom_subgraph_init(s, vertexw, zoom_adjncy, zoom_xadj,
                                    zoom_edgew, zoom_sub_offsets);

  int *zoom_part = (int *)malloc(sizeof(int) * nr_zoom_cells);
  if (zoom_part == NULL) error("Failed to allocate zoom subgraph celllist");

  zoom_partition_metis_subgraph(nodeID, nregions, nr_zoom_cells, zoom_nedges,
                                vertexw, zoom_edgew, zoom_xadj, zoom_adjncy,
                                zoom_part);

  free(zoom_sub_offsets);
  free(zoom_xadj);
  swift_free("partition_adjncy", zoom_adjncy);
  if (zoom_edgew != NULL)
    swift_free("partition_edge_weights", zoom_edgew);

  /* --- Background subgraph --- */

  int *bkg_sub_offsets = (int *)malloc(sizeof(int) * (nr_bkg_cells + 1));
  if (bkg_sub_offsets == NULL)
    error("Failed to allocate bkg subgraph offsets");

  const int bkg_nedges =
      zoom_partition_count_bkg_subgraph_edges(s, s->periodic, bkg_sub_offsets);

  idx_t *bkg_xadj = (idx_t *)malloc(sizeof(idx_t) * (nr_bkg_cells + 1));
  if (bkg_xadj == NULL) error("Failed to allocate bkg subgraph xadj");

  idx_t *bkg_adjncy = (idx_t *)swift_malloc("partition_adjncy",
                                            sizeof(idx_t) * bkg_nedges);
  if (bkg_adjncy == NULL) error("Failed to allocate bkg subgraph adjncy");

  /* The bkg subgraph wants vertex weights in the bkg index frame. */
  double *bkg_vertexw = NULL;
  if (vertexw != NULL) {
    if ((bkg_vertexw = (double *)malloc(sizeof(double) * nr_bkg_cells)) == NULL)
      error("Failed to allocate bkg subgraph vertex weights");
    for (int k = 0; k < nr_bkg_cells; k++)
      bkg_vertexw[k] = vertexw[k + bkg_offset];
  }

  double *bkg_edgew = NULL;
  if (vertexw != NULL) {
    if ((bkg_edgew = (double *)swift_malloc(
             "partition_edge_weights", sizeof(double) * bkg_nedges)) == NULL)
      error("Failed to allocate bkg subgraph edge weights");
  }

  zoom_partition_bkg_subgraph_init(s, s->periodic, vertexw, bkg_adjncy,
                                   bkg_xadj, bkg_edgew, bkg_sub_offsets);

  int *bkg_part = (int *)malloc(sizeof(int) * nr_bkg_cells);
  if (bkg_part == NULL) error("Failed to allocate bkg subgraph celllist");

  zoom_partition_metis_subgraph(nodeID, nregions, nr_bkg_cells, bkg_nedges,
                                bkg_vertexw, bkg_edgew, bkg_xadj, bkg_adjncy,
                                bkg_part);

  free(bkg_sub_offsets);
  free(bkg_xadj);
  swift_free("partition_adjncy", bkg_adjncy);
  if (bkg_vertexw != NULL) free(bkg_vertexw);
  if (bkg_edgew != NULL) swift_free("partition_edge_weights", bkg_edgew);

  /* --- Stitch the two partitions into the full-space celllist. --- */

  for (int k = 0; k < nr_zoom_cells; k++) celllist[k] = zoom_part[k];
  for (int k = 0; k < nr_bkg_cells; k++)
    celllist[k + bkg_offset] = bkg_part[k];

  free(zoom_part);
  free(bkg_part);
}

#ifdef HAVE_PARMETIS

/* Forward declaration of partition_metis.c helper used below. */
extern void permute_regions(int *newlist, int *oldlist, int nregions,
                            int ncells, int *permcelllist);

/**
 * @brief Run ParMETIS on a single subgraph and write its partitioning.
 *
 * Mirrors partition_pick_parmetis() but operates on a pre-built subgraph CSR.
 * Rank 0 holds the full subgraph and scatters per-rank slabs; all ranks then
 * call ParMETIS_V3_*; rank 0 gathers regionids, runs permute_regions() if a
 * prior partition is available, and finally MPI_Bcasts the result.
 *
 * @param nodeID The local MPI rank.
 * @param nregions The number of partitions to produce.
 * @param nverts Total number of vertices in the subgraph.
 * @param nedges Total number of edges in the subgraph (size of full_adjncy).
 * @param vertexw Per-vertex weights in the subgraph index space, or NULL.
 * @param edgew Per-edge weights in CSR order, or NULL.
 * @param full_xadj Full subgraph CSR row offsets, size nverts + 1.
 * @param full_adjncy Full subgraph CSR adjacency, size nedges.
 * @param refine Whether to refine an existing partition (sub_celllist input).
 * @param adaptive Whether to use adaptive repartition (only with refine).
 * @param itr Inter-process communication / data redistribution ratio.
 * @param sub_celllist In/out partition array, size nverts. On entry holds the
 *        prior partition (when refine=1); on exit holds the new partition.
 */
static void zoom_partition_parmetis_subgraph(
    int nodeID, int nregions, int nverts, int nedges, double *vertexw,
    double *edgew, idx_t *full_xadj, idx_t *full_adjncy, int refine,
    int adaptive, float itr, int *sub_celllist) {

  int res;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  /* Single-rank shortcut. */
  if (nregions == 1) {
    for (int i = 0; i < nverts; i++) sub_celllist[i] = 0;
    return;
  }

  /* Build vtxdist: equal partitioning of vertices across ranks. */
  idx_t *vtxdist;
  if ((vtxdist = (idx_t *)malloc(sizeof(idx_t) * (nregions + 1))) == NULL)
    error("Failed to allocate vtxdist buffer.");

  if (nodeID == 0) {
    vtxdist[0] = 0;
    int k = nverts;
    for (int i = 0; i < nregions; i++) {
      int l = k / (nregions - i);
      vtxdist[i + 1] = vtxdist[i] + l;
      k -= l;
    }
    res = MPI_Bcast((void *)vtxdist, nregions + 1, IDX_T, 0, comm);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast vtxdist.");
  } else {
    res = MPI_Bcast((void *)vtxdist, nregions + 1, IDX_T, 0, comm);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast vtxdist.");
  }

  /* Local sizes. */
  int nverts_local = vtxdist[nodeID + 1] - vtxdist[nodeID];
  int nedge_local =
      full_xadj[vtxdist[nodeID + 1]] - full_xadj[vtxdist[nodeID]];

  /* On non-root ranks the sliced full_xadj/adjncy are NULL; we need the local
   * pieces in zero-based form for ParMETIS. */
  idx_t *xadj = NULL;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (nverts_local + 1))) == NULL)
    error("Failed to allocate local xadj buffer.");

  idx_t *adjncy = NULL;
  if ((adjncy = (idx_t *)swift_malloc(
           "partition_adjncy_local", sizeof(idx_t) * nedge_local)) == NULL)
    error("Failed to allocate local adjncy array.");

  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts_local)) == NULL)
      error("Failed to allocate local vertex weights array");

  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)swift_malloc(
             "partition_edge_weights_local", sizeof(idx_t) * nedge_local)) ==
        NULL)
      error("Failed to allocate local edge weights array");

  idx_t *regionid = NULL;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * (nverts_local + 1))) == NULL)
    error("Failed to allocate regionid array");

  /* Async request bookkeeping. */
  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 5 * nregions)) ==
      NULL)
    error("Failed to allocate MPI request list.");
  for (int k = 0; k < 5 * nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  MPI_Status *stats;
  if ((stats = (MPI_Status *)malloc(sizeof(MPI_Status) * 5 * nregions)) == NULL)
    error("Failed to allocate MPI status list.");

  if (nodeID == 0) {

    /* Build a 0-based xadj per rank into a packed buffer. */
    idx_t *full_xadj_packed = NULL;
    if ((full_xadj_packed = (idx_t *)malloc(
             sizeof(idx_t) * (nverts + nregions + 1))) == NULL)
      error("Failed to allocate full_xadj_packed buffer.");

    for (int rank = 0, i = 0, j = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];
      int offset = full_xadj[j];
      for (int k = 0; k < nvt; k++) {
        full_xadj_packed[i] = full_xadj[j] - offset;
        j++;
        i++;
      }
      full_xadj_packed[i] = full_xadj[j] - offset;
      i++;
    }

    /* Cast/clamp weights into idx_t in-place buffers. */
    idx_t *full_weights_v = NULL;
    if (vertexw != NULL) {
      if ((full_weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
        error("Failed to allocate full vertex weights array");
      for (int k = 0; k < nverts; k++) {
        full_weights_v[k] = (vertexw[k] > 1) ? (idx_t)vertexw[k] : 0;
      }
    }

    idx_t *full_weights_e = NULL;
    if (edgew != NULL) {
      if ((full_weights_e = (idx_t *)swift_malloc(
               "partition_edge_weights_full", sizeof(idx_t) * nedges)) == NULL)
        error("Failed to allocate full edge weights array");
      for (int k = 0; k < nedges; k++) {
        full_weights_e[k] = (edgew[k] > 1) ? (idx_t)edgew[k] : 1;
      }
    }

    idx_t *full_regionid = NULL;
    if (refine) {
      if ((full_regionid = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
        error("Failed to allocate full regionid array");
      for (int k = 0; k < nverts; k++) full_regionid[k] = sub_celllist[k];
    }

    /* Send slabs to other ranks; copy own slab locally. */
    for (int rank = 0, j1 = 0, j2 = 0, j3 = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];
      int nedge = full_xadj[vtxdist[rank + 1]] - full_xadj[vtxdist[rank]];

      if (rank == 0) {
        memcpy(xadj, &full_xadj_packed[j1], sizeof(idx_t) * (nvt + 1));
        memcpy(adjncy, &full_adjncy[j2], sizeof(idx_t) * nedge);
        if (weights_e != NULL)
          memcpy(weights_e, &full_weights_e[j2], sizeof(idx_t) * nedge);
        if (weights_v != NULL)
          memcpy(weights_v, &full_weights_v[j3], sizeof(idx_t) * nvt);
        if (refine) memcpy(regionid, full_regionid, sizeof(idx_t) * nvt);
      } else {
        res = MPI_Isend(&full_xadj_packed[j1], nvt + 1, IDX_T, rank, 0, comm,
                        &reqs[5 * rank + 0]);
        if (res == MPI_SUCCESS)
          res = MPI_Isend(&full_adjncy[j2], nedge, IDX_T, rank, 1, comm,
                          &reqs[5 * rank + 1]);
        if (res == MPI_SUCCESS && weights_e != NULL)
          res = MPI_Isend(&full_weights_e[j2], nedge, IDX_T, rank, 2, comm,
                          &reqs[5 * rank + 2]);
        if (res == MPI_SUCCESS && weights_v != NULL)
          res = MPI_Isend(&full_weights_v[j3], nvt, IDX_T, rank, 3, comm,
                          &reqs[5 * rank + 3]);
        if (refine && res == MPI_SUCCESS)
          res = MPI_Isend(&full_regionid[j3], nvt, IDX_T, rank, 4, comm,
                          &reqs[5 * rank + 4]);
        if (res != MPI_SUCCESS) mpi_error(res, "Failed to send graph data");
      }
      j1 += nvt + 1;
      j2 += nedge;
      j3 += nvt;
    }

    int result;
    if ((result = MPI_Waitall(5 * nregions, reqs, stats)) != MPI_SUCCESS)
      error("Failed during waitall sending subgraph data.");

    if (full_weights_v != NULL) free(full_weights_v);
    if (full_weights_e != NULL)
      swift_free("partition_edge_weights_full", full_weights_e);
    free(full_xadj_packed);
    if (full_regionid != NULL) free(full_regionid);

  } else {

    /* Receive xadj first to confirm local edge count. */
    res = MPI_Recv(xadj, nverts_local + 1, IDX_T, 0, 0, comm,
                   MPI_STATUS_IGNORE);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive xadj data");

    const int nedge = xadj[nverts_local];
    if (nedge != nedge_local)
      error("Inconsistent local edge count (xadj=%d, offsets=%d)", nedge,
            nedge_local);

    res = MPI_Irecv(adjncy, nedge, IDX_T, 0, 1, comm, &reqs[0]);
    if (res == MPI_SUCCESS && weights_e != NULL)
      res = MPI_Irecv(weights_e, nedge, IDX_T, 0, 2, comm, &reqs[1]);
    if (res == MPI_SUCCESS && weights_v != NULL)
      res = MPI_Irecv(weights_v, nverts_local, IDX_T, 0, 3, comm, &reqs[2]);
    if (refine && res == MPI_SUCCESS)
      res +=
          MPI_Irecv((void *)regionid, nverts_local, IDX_T, 0, 4, comm,
                    &reqs[3]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive graph data");

    const int nreq = 1 + (weights_e != NULL) + (weights_v != NULL) + refine;
    int result;
    if ((result = MPI_Waitall(nreq, reqs, stats)) != MPI_SUCCESS)
      error("Failed during waitall receiving subgraph data.");
  }

  /* tpwgts: balanced. */
  real_t *tpwgts;
  if ((tpwgts = (real_t *)malloc(sizeof(real_t) * nregions)) == NULL)
    error("Failed to allocate tpwgts array");
  for (int i = 0; i < nregions; i++) tpwgts[i] = 1.0 / (real_t)nregions;

  idx_t options[4];
  options[0] = 1;
  options[1] = 0;

  idx_t edgecut;
  idx_t ncon = 1;
  idx_t nparts = nregions;
  idx_t numflag = 0;
  idx_t wgtflag = 0;
  if (edgew != NULL) wgtflag += 1;
  if (vertexw != NULL) wgtflag += 2;

  real_t ubvec[1];
  ubvec[0] = 1.001;

  if (refine) {
    options[3] = PARMETIS_PSR_UNCOUPLED;
    options[2] = clocks_random_seed();

    if (adaptive) {
      real_t itr_real_t = itr;
      if (ParMETIS_V3_AdaptiveRepart(
              vtxdist, xadj, adjncy, weights_v, NULL, weights_e, &wgtflag,
              &numflag, &ncon, &nparts, tpwgts, ubvec, &itr_real_t, options,
              &edgecut, regionid, &comm) != METIS_OK)
        error("Call to ParMETIS_V3_AdaptiveRepart failed.");
    } else {
      if (ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                                 &wgtflag, &numflag, &ncon, &nparts, tpwgts,
                                 ubvec, options, &edgecut, regionid,
                                 &comm) != METIS_OK)
        error("Call to ParMETIS_V3_RefineKway failed.");
    }
  } else {
    /* Best-of-N PartKway like the combined-graph driver. */
    idx_t best_edgecut = 0;
    idx_t *best_regionid = NULL;
    if ((best_regionid = (idx_t *)malloc(
             sizeof(idx_t) * (nverts_local + 1))) == NULL)
      error("Failed to allocate best_regionid array");

    for (int i = 0; i < 10; i++) {
      options[2] = clocks_random_seed();

      if (ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                               &wgtflag, &numflag, &ncon, &nparts, tpwgts,
                               ubvec, options, &edgecut, regionid,
                               &comm) != METIS_OK)
        error("Call to ParMETIS_V3_PartKway failed.");

      if (i == 0 || (best_edgecut > edgecut)) {
        best_edgecut = edgecut;
        memcpy(best_regionid, regionid, sizeof(idx_t) * (nverts_local + 1));
      }
    }
    memcpy(regionid, best_regionid, sizeof(idx_t) * (nverts_local + 1));
    free(best_regionid);
  }

  /* Gather regionid arrays back to rank 0. */
  for (int k = 0; k < nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  if (nodeID != 0) {
    res = MPI_Isend(regionid, vtxdist[nodeID + 1] - vtxdist[nodeID], IDX_T, 0,
                    1, comm, &reqs[0]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to send new regionids");
    int err;
    if ((err = MPI_Wait(reqs, stats)) != MPI_SUCCESS)
      mpi_error(err, "Failed during wait sending regionids.");
  } else {
    idx_t *remoteids = NULL;
    if ((remoteids = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate remoteids buffer");

    int nvt = vtxdist[1] - vtxdist[0];
    memcpy(remoteids, regionid, sizeof(idx_t) * nvt);

    for (int rank = 1, j = nvt; rank < nregions; rank++) {
      nvt = vtxdist[rank + 1] - vtxdist[rank];
      res = MPI_Irecv((void *)&remoteids[j], nvt, IDX_T, rank, 1, comm,
                      &reqs[rank]);
      if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive new regionids");
      j += nvt;
    }

    int err;
    if ((err = MPI_Waitall(nregions, reqs, stats)) != MPI_SUCCESS)
      error("Failed during waitall receiving regionid data.");

    int *newcelllist = NULL;
    if ((newcelllist = (int *)malloc(sizeof(int) * nverts)) == NULL)
      error("Failed to allocate new celllist");
    for (int k = 0; k < nverts; k++) newcelllist[k] = remoteids[k];
    free(remoteids);

    int bad = 0;
    for (int k = 0; k < nverts; k++) {
      if (newcelllist[k] < 0 || newcelllist[k] >= nregions) {
        message("Got bad nodeID %d for subgraph cell %i.", newcelllist[k], k);
        bad++;
      }
    }
    if (bad) error("Bad node IDs located");

    /* Permute against the input partition (when refining or when caller passed
     * a non-zero prior partition). */
    int permute = 1;
    if (!refine) {
      int nsum = 0;
      for (int k = 0; k < nverts; k++) nsum += sub_celllist[k];
      if (nsum == 0) permute = 0;
    }

    if (permute) {
      int *permcelllist = NULL;
      if ((permcelllist = (int *)malloc(sizeof(int) * nverts)) == NULL)
        error("Failed to allocate perm celllist array");
      permute_regions(newcelllist, sub_celllist, nregions, nverts,
                      permcelllist);
      memcpy(sub_celllist, permcelllist, sizeof(int) * nverts);
      free(permcelllist);
    } else {
      memcpy(sub_celllist, newcelllist, sizeof(int) * nverts);
    }
    free(newcelllist);
  }

  /* Broadcast final subgraph partition. */
  res = MPI_Bcast(sub_celllist, nverts, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast subgraph celllist");

  free(reqs);
  free(stats);
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL)
    swift_free("partition_edge_weights_local", weights_e);
  free(vtxdist);
  free(tpwgts);
  free(xadj);
  swift_free("partition_adjncy_local", adjncy);
  free(regionid);
  MPI_Comm_free(&comm);
}

/**
 * @brief Partition the zoom and background grids independently with ParMETIS.
 *
 * As zoom_partition_pick_metis() but driving ParMETIS. The two top-level
 * subgraphs (zoom and background) are partitioned independently and stitched
 * into the full-space celllist on exit.
 *
 * @param nodeID The local MPI rank.
 * @param s The #space containing the zoom and background top-level cells.
 * @param nregions The number of partitions to produce.
 * @param vertexw Per-cell vertex weights, sizeof number of cells, or NULL.
 * @param edgew Combined-graph edge weights from the caller, ignored here.
 * @param refine Whether to refine an existing partition.
 * @param adaptive Whether to use adaptive repartition (only with refine).
 * @param itr Inter-process communication / data redistribution ratio.
 * @param celllist In/out partition array, size s->nr_cells, indexed in the
 *        full-space frame.
 */
void zoom_partition_pick_parmetis(int nodeID, struct space *s, int nregions,
                                  double *vertexw, double *edgew, int refine,
                                  int adaptive, float itr, int *celllist) {

  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;
  const int nr_bkg_cells = s->zoom_props->nr_bkg_cells;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;

  /* The combined-graph edge weights are not used; we recompute per-subgraph
   * weights from vertexw inside each subgraph initialiser. */
  (void)edgew;

  /* --- Zoom subgraph --- */

  int *zoom_sub_offsets =
      (int *)malloc(sizeof(int) * (nr_zoom_cells + 1));
  if (zoom_sub_offsets == NULL)
    error("Failed to allocate zoom subgraph offsets");

  const int zoom_nedges =
      zoom_partition_count_zoom_subgraph_edges(s, zoom_sub_offsets);

  idx_t *zoom_xadj =
      (idx_t *)malloc(sizeof(idx_t) * (nr_zoom_cells + 1));
  if (zoom_xadj == NULL) error("Failed to allocate zoom subgraph xadj");

  idx_t *zoom_adjncy = (idx_t *)swift_malloc(
      "partition_adjncy", sizeof(idx_t) * zoom_nedges);
  if (zoom_adjncy == NULL) error("Failed to allocate zoom subgraph adjncy");

  double *zoom_edgew = NULL;
  if (vertexw != NULL) {
    if ((zoom_edgew = (double *)swift_malloc(
             "partition_edge_weights", sizeof(double) * zoom_nedges)) == NULL)
      error("Failed to allocate zoom subgraph edge weights");
  }

  zoom_partition_zoom_subgraph_init(s, vertexw, zoom_adjncy, zoom_xadj,
                                    zoom_edgew, zoom_sub_offsets);

  int *zoom_part = (int *)malloc(sizeof(int) * nr_zoom_cells);
  if (zoom_part == NULL) error("Failed to allocate zoom subgraph celllist");

  /* Seed sub_celllist from caller's input (for refine and permute). */
  for (int k = 0; k < nr_zoom_cells; k++) zoom_part[k] = celllist[k];

  zoom_partition_parmetis_subgraph(nodeID, nregions, nr_zoom_cells,
                                   zoom_nedges, vertexw, zoom_edgew,
                                   zoom_xadj, zoom_adjncy, refine, adaptive,
                                   itr, zoom_part);

  free(zoom_sub_offsets);
  free(zoom_xadj);
  swift_free("partition_adjncy", zoom_adjncy);
  if (zoom_edgew != NULL)
    swift_free("partition_edge_weights", zoom_edgew);

  /* --- Background subgraph --- */

  int *bkg_sub_offsets = (int *)malloc(sizeof(int) * (nr_bkg_cells + 1));
  if (bkg_sub_offsets == NULL)
    error("Failed to allocate bkg subgraph offsets");

  const int bkg_nedges =
      zoom_partition_count_bkg_subgraph_edges(s, s->periodic, bkg_sub_offsets);

  idx_t *bkg_xadj = (idx_t *)malloc(sizeof(idx_t) * (nr_bkg_cells + 1));
  if (bkg_xadj == NULL) error("Failed to allocate bkg subgraph xadj");

  idx_t *bkg_adjncy = (idx_t *)swift_malloc(
      "partition_adjncy", sizeof(idx_t) * bkg_nedges);
  if (bkg_adjncy == NULL) error("Failed to allocate bkg subgraph adjncy");

  double *bkg_vertexw = NULL;
  if (vertexw != NULL) {
    if ((bkg_vertexw = (double *)malloc(sizeof(double) * nr_bkg_cells)) == NULL)
      error("Failed to allocate bkg subgraph vertex weights");
    for (int k = 0; k < nr_bkg_cells; k++)
      bkg_vertexw[k] = vertexw[k + bkg_offset];
  }

  double *bkg_edgew = NULL;
  if (vertexw != NULL) {
    if ((bkg_edgew = (double *)swift_malloc(
             "partition_edge_weights", sizeof(double) * bkg_nedges)) == NULL)
      error("Failed to allocate bkg subgraph edge weights");
  }

  zoom_partition_bkg_subgraph_init(s, s->periodic, vertexw, bkg_adjncy,
                                   bkg_xadj, bkg_edgew, bkg_sub_offsets);

  int *bkg_part = (int *)malloc(sizeof(int) * nr_bkg_cells);
  if (bkg_part == NULL) error("Failed to allocate bkg subgraph celllist");

  /* Seed sub_celllist from caller's input (for refine and permute). */
  for (int k = 0; k < nr_bkg_cells; k++)
    bkg_part[k] = celllist[k + bkg_offset];

  zoom_partition_parmetis_subgraph(nodeID, nregions, nr_bkg_cells, bkg_nedges,
                                   bkg_vertexw, bkg_edgew, bkg_xadj,
                                   bkg_adjncy, refine, adaptive, itr,
                                   bkg_part);

  free(bkg_sub_offsets);
  free(bkg_xadj);
  swift_free("partition_adjncy", bkg_adjncy);
  if (bkg_vertexw != NULL) free(bkg_vertexw);
  if (bkg_edgew != NULL) swift_free("partition_edge_weights", bkg_edgew);

  /* --- Stitch the two partitions into the full-space celllist. --- */

  for (int k = 0; k < nr_zoom_cells; k++) celllist[k] = zoom_part[k];
  for (int k = 0; k < nr_bkg_cells; k++)
    celllist[k + bkg_offset] = bkg_part[k];

  free(zoom_part);
  free(bkg_part);
}

#endif /* HAVE_PARMETIS */
#endif

/* Void cell support */
/* ================= */

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
