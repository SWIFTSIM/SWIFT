/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "partition.h"
#include "proxy.h"
#include "space.h"
#include "zoom_region.h"

#include <float.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Check if all regions have been assigned a node in the
 *        cells of a space.
 *
 * @param s the space containing the cells to check.
 * @param nregions number of regions expected.
 * @param verbose if true report the missing regions.
 * @return true if all regions have been found, false otherwise.
 */
static int check_complete(struct space *s, int verbose, int nregions) {

  int *present = NULL;
  if ((present = (int *)malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate present array");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    if (s->cells_top[i].nodeID <= nregions)
      present[s->cells_top[i].nodeID]++;
    else
      message("Bad nodeID: s->cells_top[%d].nodeID = %d", i,
              s->cells_top[i].nodeID);
  }
  for (int i = 0; i < nregions; i++) {
    if (!present[i]) {
      failed = 1;
      if (verbose) message("Region %d is not present in partition", i);
    }
  }
  free(present);
  return (!failed);
}

/**
 * @brief Partition a space of cells based on another space of cells.
 *
 * The two spaces are expected to be at different cell sizes, so what we'd
 * like to do is assign the second space to geometrically closest nodes
 * of the first, with the effect of minimizing particle movement when
 * rebuilding the second space from the first.
 *
 * Since two spaces cannot exist simultaneously the old space is actually
 * required in a decomposed state. These are the old cells sizes and counts
 * per dimension, along with a list of the old nodeIDs. The old nodeIDs are
 * indexed by the cellid (see cell_getid()), so should be stored that way.
 *
 * On exit the new space cells will have their nodeIDs assigned.
 *
 * @param oldh the cell dimensions of old space.
 * @param oldcdim number of cells per dimension in old space.
 * @param oldnodeIDs the nodeIDs of cells in the old space, indexed by old
 *cellid.
 * @param s the space to be partitioned.
 *
 * @return 1 if the new space contains nodeIDs from all nodes, 0 otherwise.
 */
int partition_space_to_space_zoom(double *oldh, double *oldcdim,
                                  double *oldzoomh, double *oldzoomcdim,
                                  int *oldnodeIDs, struct space *s) {

  /* Define the old tl_cell_offset */
  const int old_bkg_cell_offset =
      oldzoomcdim[0] * oldzoomcdim[1] * oldzoomcdim[2];

  /* Loop over all the new zoom cells. */
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {

        /* Scale indices to old cell space. */
        const int ii = rint(i * s->zoom_props->iwidth[0] * oldzoomh[0]);
        const int jj = rint(j * s->zoom_props->iwidth[1] * oldzoomh[1]);
        const int kk = rint(k * s->zoom_props->iwidth[2] * oldzoomh[2]);

        const int cid = cell_getid(s->zoom_props->cdim, i, j, k);
        const int oldcid = cell_getid(oldzoomcdim, ii, jj, kk);
        s->cells_top[cid].nodeID = oldnodeIDs[oldcid];
      }
    }
  }

  /* Loop over all the new cells. */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Scale indices to old cell space. */
        const int ii = rint(i * s->iwidth[0] * oldh[0]);
        const int jj = rint(j * s->iwidth[1] * oldh[1]);
        const int kk = rint(k * s->iwidth[2] * oldh[2]);

        const int cid =
            cell_getid(s->cdim, i, j, k) + s->zoom_props->tl_cell_offset;
        const int oldcid =
            cell_getid(oldcdim, ii, jj, kk) + old_bkg_cell_offset;
        s->cells_top[cid].nodeID = oldnodeIDs[oldcid];
      }
    }
  }

  /* Check we have all nodeIDs present in the resample. */
  return check_complete(s, 1, s->e->nr_nodes);
}

/**
 * @brief A genric looping function to handle duplicated looping methods done
 *        for edge counting, adjancency, and edges weighting.
 *
 * The function will do the correct operation based on what is passed.
 *
 * @param s the space of cells.
 * @param adjncy the adjncy array to fill.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 * @param counts the number of bytes in particles per cell.
 * @param edges weights for the edges of these regions.
 */
void edge_loop(int *cdim, int cell_type, struct space *s,
               idx_t *adjncy, idx_t *xadj, double *counts, double *edges,
               int iedge) {
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  /* Get some useful constants. */
  const int periodic = s->periodic;
  const int bkg_cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const int buffer_cdim[3] = {s->zoom_props->buffer_cdim[0],
                              s->zoom_props->buffer_cdim[1],
                              s->zoom_props->buffer_cdim[2]};
  const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                            s->zoom_props->cdim[2]};
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;
  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;
  const double  buffer_bounds[6] = {
    s->zoom_props->buffer_bounds[0], s->zoom_props->buffer_bounds[1],
    s->zoom_props->buffer_bounds[2], s->zoom_props->buffer_bounds[3],
    s->zoom_props->buffer_bounds[4], s->zoom_props->buffer_bounds[5]};
  struct cell *restrict ci;
  struct cell *restrict zoom_ci;
  struct cell *restrict buffer_ci;
  struct cell *restrict cj;
  int top_i, top_j, top_k;

  /* Define the cell type */
  enum tl_cell_types tl_cell_type;
  if (cell_type == 0) {
    tl_cell_type = tl_cell;
  } else if (cell_type == 3) {
    tl_cell_type = zoom_tl_cell;
  } else {
    tl_cell_type = buffer_tl_cell;
  }

  /* Loop over the provided cells and find their edges. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
          
        /* Get the cell index. */
        const size_t cid = cell_getid(cdim, i, j, k);

        /* Include the correct offset. */
        if (tl_cell_type == tl_cell) {
          cid += bkg_cell_offset;
        } else if (tl_cell_type == buffer_tl_cell) {
          cid += buffer_cell_offset;
        }

        /* Get the cell. */
        ci = &s->cells_top[cid];

#ifdef SWIFT_DEBUG_CHECKS
        if (xadj != NULL) {
             
          /* Ensure the previous cell has found enough edges. */
          if ((cid > 0) &&
              ((iedge - xadj[cid - 1]) != s->cells_top[cid - 1].nr_vertex_edges))
            error("Found too few edges (nedges=%d, c->nr_vertex_edges=%d)",
                  iedge - xadj[cid - 1], s->cells_top[cid - 1].nr_vertex_edges);
          
        }
#endif

        /* If given set METIS xadj. */
        if (xadj != NULL) {
          xadj[cid] = iedge;
          
          /* Set edges start pointer for this cell. */
          ci->edges_start = iedge;
        }

        /* Loop over a shell of cells with the same type. */
        for (int ii = i - 1; ii <= i + 1; ii++) {
          if ((tl_cell_type == zoom_tl_cell ||
               tl_cell_type == buffer_tl_cell || !periodic) &&
              (ii < 0 || ii >= cdim[0])) continue;
          for (int jj = j - 1; jj <= j + 1; jj++) {
            if ((tl_cell_type == zoom_tl_cell ||
                 tl_cell_type == buffer_tl_cell || !periodic) &&
                (jj < 0 || jj >= cdim[1])) continue;
            for (int kk = k - 1; kk <= k + 1; kk++) {
              if ((tl_cell_type == zoom_tl_cell ||
                   tl_cell_type == buffer_tl_cell || !periodic) &&
                  (kk < 0 || kk >= cdim[2])) continue;

              /* Apply periodic BC (not harmful if not using periodic BC) */
              const int iii = (ii + cdim[0]) % cdim[0];
              const int jjj = (jj + cdim[1]) % cdim[1];
              const int kkk = (kk + cdim[2]) %cdim[2];
                
              /* Get cell index. */
              size_t cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Include the correct offset. */
              if (tl_cell_type == tl_cell) {
                cjd += bkg_cell_offset;
              } else if (tl_cell_type == buffer_tl_cell) {
                cjd += buffer_cell_offset;
              }

              /* Skip self. */
              if (cid == cjd) continue;

              /* If we have an edge do the operation. */
              if (is_edge) {
                
                /* Handle size_to_edges case */
                if (edges != NULL) {
                  /* Store this edge. */
                  edges[iedge] = counts[zoom_cjd];
                  iedge++;
                }
                
                /* Handle graph_init case */
                else if (adjncy != NULL) {
                  adjncy[iedge] = cbd;
                  iedge++;
                }

                /* Handle find_vertex_edges case */
                else if (adjcny == NULL && counts == NULL) {
                  /* If not self record an edge. */
                  c->nr_vertex_edges++;
                }
                
                else {
                  error("Incompatible arguments supplied!");
                }
              }
              
            } /* neighbour k loop */
          } /* neighbour j loop */
        } /* neighbour i loop */

        /* Loop over a shell of cells with the same type and find all relations
         * between levels. */
        for (int ii = i - 1; ii <= i + 1; ii++) {
          if ((tl_cell_type == zoom_tl_cell ||
               tl_cell_type == buffer_tl_cell || !periodic) &&
              (ii < 0 || ii >= cdim[0])) continue;
          for (int jj = j - 1; jj <= j + 1; jj++) {
            if ((tl_cell_type == zoom_tl_cell ||
                 tl_cell_type == buffer_tl_cell || !periodic) &&
                (jj < 0 || jj >= cdim[1])) continue;
            for (int kk = k - 1; kk <= k + 1; kk++) {
              if ((tl_cell_type == zoom_tl_cell ||
                   tl_cell_type == buffer_tl_cell || !periodic) &&
                  (kk < 0 || kk >= cdim[2])) continue;

              /* Apply periodic BC (not harmful if not using periodic BC) */
              const int iii = (ii + cdim[0]) % cdim[0];
              const int jjj = (jj + cdim[1]) % cdim[1];
              const int kkk = (kk + cdim[2]) %cdim[2];
              
              /* Get cell index. */
              size_t cjd = cell_getid(cdim, iii, jjj, kkk);
              
              /* Include the correct offset. */
              if (tl_cell_type == tl_cell) {
                cjd += bkg_cell_offset;
              } else if (tl_cell_type == buffer_tl_cell) {
                cjd += buffer_cell_offset;
              }

              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Loop over zoom cells and find the edges due to nesting.  */
              for (int zoom_i = 0; zoom_i < zoom_cdim[0]; zoom_i++) {
                for (int zoom_j = 0; zoom_j < zoom_cdim[1]; zoom_j++) {
                  for (int zoom_k = 0; zoom_k < zoom_cdim[2]; zoom_k++) {

                    /* Flag for whether we have an edge. */
                    int is_edge = 0; 

                    /* Early skip if we don't need zoom cells. */
                    if (tl_cell_type == zoom_tl_cell ||
                        (tl_cell_type == tl_cell ||
                         s->zoom_props->with_buffer_cells))
                      continue;

                    /* Get the cell index. */
                    const size_t zoom_cid =
                      cell_getid(zoom_cdim, zoom_i, zoom_j, zoom_k);

                    /* Get the cell. */
                    zoom_ci = &s->cells_top[zoom_cid];

                    /* Handle background cell case */
                    if (tl_cell_type == tl_cell) {
                      
                      /* Get the (i,j,k) location of the top-level cell in the lower
                       * resolution grid we are in. */
                      top_i = zoom_ci->loc[0] * s->iwidth[0];
                      top_j = zoom_ci->loc[1] * s->iwidth[1];
                      top_k = zoom_ci->loc[2] * s->iwidth[2];

                      /* Get cell index. */
                      const size_t zoom_top =
                        cell_getid(bkg_cdim, top_i, top_j, top_k) +
                        bkg_cell_offset;

                      /* Do we have an edge? */
                      if (cjd == zoom_top)
                        is_edge = 1;
                    }
                    
                    /* Handle buffer cell case. */
                    else if (tl_cell_type == buffer_tl_cell) {

                      /* Get the (i,j,k) location of the top-level cell in the
                       * lower resolution grid we are in. */
                      top_i = (zoom_ci->loc[0] - buffer_bounds[0]) *
                        s->zoom_props->buffer_iwidth[0];
                      top_j = (zoom_ci->loc[1] - buffer_bounds[2]) *
                        s->zoom_props->buffer_iwidth[1];
                      top_k = (zoom_ci->loc[2] - buffer_bounds[4]) *
                        s->zoom_props->buffer_iwidth[2];

                      /* Get cell index. */
                      const size_t zoom_top =
                        cell_getid(buffer_cdim, top_i, top_j, top_k) +
                        buffer_cell_offset;

                      /* Do we have an edge? */
                      if (cjd == zoom_top)
                        is_edge = 1;
                    }

                    /* If we have an edge do the operation. */
                    if (is_edge) {

                      /* Handle size_to_edges case */
                      if (edges != NULL) {
                        /* Store this edge. */
                        edges[iedge] = counts[zoom_cjd];
                        iedge++;
                      }

                      /* Handle graph_init case */
                      else if (adjncy != NULL) {
                        adjncy[iedge] = cbd;
                        iedge++;
                      }

                      /* Handle find_vertex_edges case */
                      else if (adjcny == NULL && counts == NULL) {
                        /* If not self record an edge. */
                        c->nr_vertex_edges++;
                      }

                      else {
                        error("Incompatible arguments supplied!");
                      }
                    }

                  } /* zoom k loop */
                } /* zoom j loop */
              } /* zoom i loop */

              /* Loop over background cells and find the edges due to nesting. */
              for (int bkg_i = 0; bkg_i < bkg_cdim[0]; bkg_i++) {
                for (int bkg_j = 0; bkg_j < bkg_cdim[1]; bkg_j++) {
                  for (int bkg_k = 0; bkg_k < bkg_cdim[2]; bkg_k++) {

                    /* Flag for whether we have an edge. */
                    int is_edge = 0; 

                    /* Early skip if we don't need background cells. */
                    if (tl_cell_type == tl_cell ||
                        (tl_cell_type == zoom_tl_cell ||
                         s->zoom_props->with_buffer_cells))
                      continue;

                    /* Get the cell index. */
                    const size_t bkg_cid =
                      cell_getid(bkg_cdim, bkg_i, bkg_j, bkg_k) +
                      bkg_cell_offset;

                    /* Get the (i,j,k) location of the top-level cell in the lower
                     * resolution grid we are in. */
                    top_i = cj->loc[0] * s->iwidth[0];
                    top_j = cj->loc[1] * s->iwidth[1];
                    top_k = cj->loc[2] * s->iwidth[2];

                    /* Get cell index. */
                    const size_t bkg_top =
                      cell_getid(bkg_cdim, top_i, top_j, top_k) +
                      bkg_cell_offset;

                    /* Do we have an edge? */
                    if (bkg_cid == bkg_top)
                      is_edge = 1;

                    /* If we have an edge do the operation. */
                    if (is_edge) {

                      /* Handle size_to_edges case */
                      if (edges != NULL) {
                        /* Store this edge. */
                        edges[iedge] = counts[zoom_cjd];
                        iedge++;
                      }

                      /* Handle graph_init case */
                      else if (adjncy != NULL) {
                        adjncy[iedge] = cbd;
                        iedge++;
                      }

                      /* Handle find_vertex_edges case */
                      else if (adjcny == NULL && counts == NULL) {
                        /* If not self record an edge. */
                        c->nr_vertex_edges++;
                      }

                      else {
                        error("Incompatible arguments supplied!");
                      }
                    }

                  } /* background k loop */
                } /* background j loop */
              } /* background i loop */

              /* Loop over buffer cells and find the edges due to nesting. */
              for (int buff_i = 0; buff_i < buff_cdim[0]; buff_i++) {
                for (int buff_j = 0; buff_j < buff_cdim[1]; buff_j++) {
                  for (int buff_k = 0; buff_k < buff_cdim[2]; buff_k++) {

                    /* Flag for whether we have an edge. */
                    int is_edge = 0; 

                    /* Early skip if we don't need buffer cells. */
                    if (tl_cell_type == buffer_tl_cell) continue;

                    /* Get the cell index. */
                    const size_t buffer_cid =
                      cell_getid(bkg_cdim, buff_i, buff_j, buff_k) +
                      buffer_cell_offset;
                    
                    /* Get the cell. */
                    buffer_ci = &s->cells_top[buffer_cid];

                    /* Handle zoom cell case */
                    if (tl_cell_type == zoom_tl_cell) {

                      /* Get the (i,j,k) location of the top-level cell in the
                       * lower resolution grid we are in. */
                      top_i = (cj->loc[0] - buffer_bounds[0]) *
                        s->zoom_props->buffer_iwidth[0];
                      top_j = (cj->loc[1] - buffer_bounds[2]) *
                        s->zoom_props->buffer_iwidth[1];
                      top_k = (cj->loc[2] - buffer_bounds[4]) *
                        s->zoom_props->buffer_iwidth[2];

                      /* Get cell index. */
                      const size_t buffer_top =
                        cell_getid(bkg_cdim, top_i, top_j, top_k) +
                        bkg_cell_offset;

                    /* Do we have an edge? */
                    if (buffer_cid == buffer_top)
                      is_edge = 1;
                    }

                    /* Handle background cell case. */
                    else if (tl_cell_type == tl_cell) {

                      /* Get the (i,j,k) location of the top-level cell in the lower
                       * resolution grid we are in. */
                      top_i = buffer_ci->loc[0] * s->iwidth[0];
                      top_j = buffer_ci->loc[1] * s->iwidth[1];
                      top_k = buffer_ci->loc[2] * s->iwidth[2];

                      /* Get cell index. */
                      const size_t buffer_top =
                        cell_getid(bkg_cdim, top_i, top_j, top_k) +
                        bkg_cell_offset;

                      /* Do we have an edge? */
                      if (cjd == buffer_top)
                        is_edge = 1;
                    }

                    /* If we have an edge do the operation. */
                    /* If we have an edge do the operation. */
                    if (is_edge) {

                      /* Handle size_to_edges case */
                      if (edges != NULL) {
                        /* Store this edge. */
                        edges[iedge] = counts[zoom_cjd];
                        iedge++;
                      }

                      /* Handle graph_init case */
                      else if (adjncy != NULL) {
                        adjncy[iedge] = cbd;
                        iedge++;
                      }

                      /* Handle find_vertex_edges case */
                      else if (adjcny == NULL && counts == NULL) {
                        /* If not self record an edge. */
                        c->nr_vertex_edges++;
                      }

                      else {
                        error("Incompatible arguments supplied!");
                      }
                    }

                  } /* buffer k loop */
                } /* buffer j loop */
              } /* buffer i loop */
              
            } /* neighbour k loop */
          } /* neighbour j loop */
        } /* neighbour i loop */
        
        if (adjcny == NULL && counts == NULL) {
          /* Include this edge count in the total. */
          s->zoom_props->nr_edges += c->nr_vertex_edges;
        }
      } /* k loop */
    } /* j loop */
  } /* i loop */
#endif
}
