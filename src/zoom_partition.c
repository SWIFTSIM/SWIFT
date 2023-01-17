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
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
void edge_loop(const int *cdim, int offset, struct space *s,
               idx_t *adjncy, idx_t *xadj, double *counts, double *edges,
               int *iedge) {

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
  const double  buffer_bounds[6] = {
    s->zoom_props->buffer_bounds[0], s->zoom_props->buffer_bounds[1],
    s->zoom_props->buffer_bounds[2], s->zoom_props->buffer_bounds[3],
    s->zoom_props->buffer_bounds[4], s->zoom_props->buffer_bounds[5]};
  struct cell *restrict ci;
  struct cell *restrict cj;
  struct cell *restrict zoom_cj;
  struct cell *restrict bkg_cj;
  struct cell *restrict buffer_cj;
  int top_i, top_j, top_k;

  /* Define the cell type */
  enum tl_cell_types tl_cell_type;
  if (offset == 0) {
    tl_cell_type = zoom_tl_cell;
  } else if (offset == bkg_cell_offset) {
    tl_cell_type = tl_cell;
  } else {
    tl_cell_type = buffer_tl_cell;
  }

  /* Loop over the provided cells and find their edges. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
          
        /* Get the cell index. */
        const size_t cid = cell_getid(cdim, i, j, k) + offset;

        /* Get the cell. */
        ci = &s->cells_top[cid];

#ifdef SWIFT_DEBUG_CHECKS
        if (xadj != NULL) {
             
          /* Ensure the previous cell has found enough edges. */
          if ((cid > 0) &&
              ((*iedge - xadj[cid - 1]) != s->cells_top[cid - 1].nr_vertex_edges))
            error("Found too few edges (nedges=%ld, c->nr_vertex_edges=%d)",
                  *iedge - xadj[cid - 1], s->cells_top[cid - 1].nr_vertex_edges);
          
        }
#endif

        /* If given set METIS xadj. */
        if (xadj != NULL) {
          xadj[cid] = *iedge;
          
          /* Set edges start pointer for this cell. */
          ci->edges_start = *iedge;
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
              const size_t cjd = cell_getid(cdim, iii, jjj, kkk) + offset;
              
              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Skip self. */
              if (cid == cjd) continue;
                
              /* Handle size_to_edges case */
              if (edges != NULL) {
                /* Store this edge. */
                edges[*iedge] = counts[cjd];
                (*iedge)++;
              }
                
              /* Handle graph_init case */
              else if (adjncy != NULL) {
                adjncy[*iedge] = cjd;
                (*iedge)++;
              }

              /* Handle find_vertex_edges case */
              else {
                /* If not self record an edge. */
                ci->nr_vertex_edges++;
                (*iedge)++;
              }
              
            } /* neighbour k loop */
          } /* neighbour j loop */
        } /* neighbour i loop */

        /* Handle cells above the current cell in the heirarchy. */

        /* The cell is below the background level. */
        if ((!s->zoom_props->with_buffer_cells &&
             tl_cell_type == zoom_tl_cell) || tl_cell_type == buffer_tl_cell) {

          /* Get the (i,j,k) location of the background cell above this one. */
          top_i = ci->loc[0] * s->iwidth[0];
          top_j = ci->loc[1] * s->iwidth[1];
          top_k = ci->loc[2] * s->iwidth[2];

          /* Loop over a shell of background cells */
          for (int ii = top_i - 1; ii <= top_i + 1; ii++) {
            if (!periodic && (ii < 0 || ii >= cdim[0])) continue;
            for (int jj = top_j - 1; jj <= top_j + 1; jj++) {
              if (!periodic && (jj < 0 || jj >= cdim[1])) continue;
              for (int kk = top_k - 1; kk <= top_k + 1; kk++) {
                if (!periodic && (kk < 0 || kk >= cdim[2])) continue;

                /* Apply periodic BC (not harmful if not using periodic BC) */
                const int iii = (ii + cdim[0]) % cdim[0];
                const int jjj = (jj + cdim[1]) % cdim[1];
                const int kkk = (kk + cdim[2]) % cdim[2];

                /* Get the cell index. */
                const size_t bkg_cjd =
                  cell_getid(bkg_cdim, iii, jjj, kkk) + bkg_cell_offset;

                /* Get the cell. */
                bkg_cj = &s->cells_top[bkg_cjd];

                /* if (bkg_cj->tl_cell_type == void_tl_cell || */
                /*     bkg_cj->tl_cell_type == void_tl_cell_neighbour) */
                /* continue; */

                /* Handle size_to_edges case */
                if (edges != NULL) {
                  /* Store this edge. */
                  edges[*iedge] = counts[bkg_cjd];
                  (*iedge)++;
                }

                /* Handle graph_init case */
                else if (adjncy != NULL) {
                  adjncy[*iedge] = bkg_cjd;
                  (*iedge)++;
                }

                /* Handle find_vertex_edges case */
                else {
                  /* If not self record an edge. */
                  ci->nr_vertex_edges++;
                  (*iedge)++;
                }

              } /* background k loop */
            } /* background j loop */
          } /* background i loop */
        } /* Below background level */

        /* Below the buffer cells level */
        if (tl_cell_type == zoom_tl_cell) {
          
          /* Get the (i,j,k) location of the top-level cell in the
           * lower resolution grid we are in. */
          top_i = (ci->loc[0] - buffer_bounds[0]) *
            s->zoom_props->buffer_iwidth[0];
          top_j = (ci->loc[1] - buffer_bounds[2]) *
            s->zoom_props->buffer_iwidth[1];
          top_k = (ci->loc[2] - buffer_bounds[4]) *
            s->zoom_props->buffer_iwidth[2];

          /* Loop over a shell of background cells */
          for (int ii = top_i - 1; ii <= top_i + 1; ii++) {
            if (ii < 0 || ii >= cdim[0]) continue;
            for (int jj = top_j - 1; jj <= top_j + 1; jj++) {
              if (jj < 0 || jj >= cdim[1]) continue;
              for (int kk = top_k - 1; kk <= top_k + 1; kk++) {
                if (kk < 0 || kk >= cdim[2]) continue;

                /* Get the cell index. */
                const size_t buff_cjd =
                  cell_getid(buffer_cdim, ii, jj, kk) + buffer_cell_offset;

                 /* Get the cell. */
                buffer_cj = &s->cells_top[buff_cjd];

                /* if (buffer_cj->tl_cell_type == void_tl_cell) */
                /*   continue; */

                /* Handle size_to_edges case */
                if (edges != NULL) {
                  /* Store this edge. */
                  edges[*iedge] = counts[buff_cjd];
                  (*iedge)++;
                }

                /* Handle graph_init case */
                else if (adjncy != NULL) {
                  adjncy[*iedge] = buff_cjd;
                  (*iedge)++;
                }

                /* Handle find_vertex_edges case */
                else {
                  /* If not self record an edge. */
                  ci->nr_vertex_edges++;
                  (*iedge)++;
                }

              } /* buffer k loop */
            } /* buffer j loop */
          } /* buffer i loop */      
        } /* Below buffer cells */

        /* Handle cells above the current cell in the heirarchy.  */

        /* Cells above the zoom cells */
        if ((!s->zoom_props->with_buffer_cells &&
             tl_cell_type == tl_cell) || tl_cell_type == buffer_tl_cell) {
          
          /* Loop over zoom cells and find the edges due to nesting.  */
          for (int zoom_i = 0; zoom_i < zoom_cdim[0]; zoom_i++) {
            for (int zoom_j = 0; zoom_j < zoom_cdim[1]; zoom_j++) {
              for (int zoom_k = 0; zoom_k < zoom_cdim[2]; zoom_k++) {

                /* Get the cell index. */
                const size_t zoom_cjd =
                  cell_getid(zoom_cdim, zoom_i, zoom_j, zoom_k);

                /* Get the cell. */
                zoom_cj = &s->cells_top[zoom_cjd];

                /* Is this cell inside ci? */
                if (zoom_cj->loc[0] >= ci->loc[0] &&
                    zoom_cj->loc[0] <= (ci->loc[0] + ci->width[0]) &&
                    zoom_cj->loc[1] >= ci->loc[1] &&
                    zoom_cj->loc[1] <= (ci->loc[1] + ci->width[1]) &&
                    zoom_cj->loc[2] >= ci->loc[2] &&
                    zoom_cj->loc[2] <= (ci->loc[2] + ci->width[2])) {

                  /* Handle size_to_edges case */
                  if (edges != NULL) {
                    /* Store this edge. */
                    edges[*iedge] = counts[zoom_cjd];
                    (*iedge)++;
                  }
                  
                  /* Handle graph_init case */
                  else if (adjncy != NULL) {
                    adjncy[*iedge] = zoom_cjd;
                    (*iedge)++;
                  }
                  
                  /* Handle find_vertex_edges case */
                  else {
                    /* If not self record an edge. */
                    ci->nr_vertex_edges++;
                    (*iedge)++;
                  }
                  
                } /* cj is inside ci */
              } /* zoom k loop */
            } /* zoom j loop */
          } /* zoom i loop */
        } /* Cell is above the zoom cells */

        /* Cells above the buffer cells */
        if (tl_cell_type == tl_cell) {
          
          /* Loop over buffer cells and find the edges due to nesting.  */
          for (int buffer_i = 0; buffer_i < buffer_cdim[0]; buffer_i++) {
            for (int buffer_j = 0; buffer_j < buffer_cdim[1]; buffer_j++) {
              for (int buffer_k = 0; buffer_k < buffer_cdim[2]; buffer_k++) {

                /* Get the cell index. */
                const size_t buffer_cjd =
                  cell_getid(buffer_cdim, buffer_i, buffer_j, buffer_k) +
                  buffer_cell_offset;

                /* Get the cell. */
                buffer_cj = &s->cells_top[buffer_cjd];

                /* if (buffer_cj->tl_cell_type == void_tl_cell) */
                /*   continue; */

                /* Is this cell inside ci? */
                if (buffer_cj->loc[0] >= ci->loc[0] &&
                    buffer_cj->loc[0] <= (ci->loc[0] + ci->width[0]) &&
                    buffer_cj->loc[1] >= ci->loc[1] &&
                    buffer_cj->loc[1] <= (ci->loc[1] + ci->width[1]) &&
                    buffer_cj->loc[2] >= ci->loc[2] &&
                    buffer_cj->loc[2] <= (ci->loc[2] + ci->width[2])) {

                  /* Handle size_to_edges case */
                  if (edges != NULL) {
                    /* Store this edge. */
                    edges[*iedge] = counts[buffer_cjd];
                    (*iedge)++;
                  }
                  
                  /* Handle graph_init case */
                  else if (adjncy != NULL) {
                    adjncy[*iedge] = buffer_cjd;
                    (*iedge)++;
                  }
                  
                  /* Handle find_vertex_edges case */
                  else {
                    /* If not self record an edge. */
                    ci->nr_vertex_edges++;
                    (*iedge)++;
                  }
                  
                } /* cj is inside ci */
              } /* buffer k loop */
            } /* buffer j loop */
          } /* buffer i loop */
        } /* Cell is above the buffer cells */
      } /* k loop */
    } /* j loop */
  } /* i loop */
}
#endif

