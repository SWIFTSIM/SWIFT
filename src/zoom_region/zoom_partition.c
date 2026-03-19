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
