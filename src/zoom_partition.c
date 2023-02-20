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
            if (!periodic && (ii < 0 || ii >= bkg_cdim[0])) continue;
            for (int jj = top_j - 1; jj <= top_j + 1; jj++) {
              if (!periodic && (jj < 0 || jj >= bkg_cdim[1])) continue;
              for (int kk = top_k - 1; kk <= top_k + 1; kk++) {
                if (!periodic && (kk < 0 || kk >= bkg_cdim[2])) continue;

                /* Apply periodic BC (not harmful if not using periodic BC) */
                const int iii = (ii + bkg_cdim[0]) % bkg_cdim[0];
                const int jjj = (jj + bkg_cdim[1]) % bkg_cdim[1];
                const int kkk = (kk + bkg_cdim[2]) % bkg_cdim[2];
                /* const int iii = top_i; */
                /* const int jjj = top_j; */
                /* const int kkk = top_k; */

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
            if (ii < 0 || ii >= buffer_cdim[0]) continue;
            for (int jj = top_j - 1; jj <= top_j + 1; jj++) {
              if (jj < 0 || jj >= buffer_cdim[1]) continue;
              for (int kk = top_k - 1; kk <= top_k + 1; kk++) {
                if (kk < 0 || kk >= buffer_cdim[2]) continue;

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

                /* What cell is above this zoom cell? */
                if (tl_cell_type == tl_cell) {
                  top_i = zoom_cj->loc[0] * s->iwidth[0];
                  top_j = zoom_cj->loc[1] * s->iwidth[1];
                  top_k = zoom_cj->loc[2] * s->iwidth[2];

                } else {

                  top_i = (zoom_cj->loc[0] - buffer_bounds[0]) *
                    s->zoom_props->buffer_iwidth[0];
                  top_j = (zoom_cj->loc[1] - buffer_bounds[2]) *
                    s->zoom_props->buffer_iwidth[1];
                  top_k = (zoom_cj->loc[2] - buffer_bounds[4]) *
                    s->zoom_props->buffer_iwidth[2];
                }

                /* Get the cell index. */
                 const size_t top_cid =
                   cell_getid(cdim, top_i, top_j, top_k) + offset;

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

                      /* Apply periodic BC (not harmful if not using periodic
                       * BC) */
                      const int iii = (ii + cdim[0]) % cdim[0];
                      const int jjj = (jj + cdim[1]) % cdim[1];
                      const int kkk = (kk + cdim[2]) % cdim[2];
                
                      /* Get cell index. */
                      const size_t top_cjd =
                        cell_getid(cdim, iii, jjj, kkk) + offset;

                      /* Is this cell inside ci? */
                      if (top_cid == top_cjd) {                  
                        
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
                        
                      } /* neighbour loop i */
                    } /* neighbour loop j */
                  } /* neighbour loop k */
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

                /* What cell is above this buffer cell? */
                top_i = buffer_cj->loc[0] * s->iwidth[0];
                top_j = buffer_cj->loc[1] * s->iwidth[1];
                top_k = buffer_cj->loc[2] * s->iwidth[2];

                /* Get the cell index. */
                const size_t top_cid =
                  cell_getid(cdim, top_i, top_j, top_k) + offset;

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

                      /* Apply periodic BC (not harmful if not using periodic
                       * BC) */
                      const int iii = (ii + cdim[0]) % cdim[0];
                      const int jjj = (jj + cdim[1]) % cdim[1];
                      const int kkk = (kk + cdim[2]) %cdim[2];
                
                      /* Get cell index. */
                      const size_t top_cjd =
                        cell_getid(cdim, iii, jjj, kkk) + offset;

                      /* Is this cell inside ci? */
                      if (top_cid == top_cjd) {                   

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
                        
                      } /* neighbour loop i */
                    } /* neighbour loop j */
                  } /* neighbour loop k */
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


/**
 * @brief Create and fill the proxies for the natural cells.
 *
 * @param e The #engine.
 */
void engine_makeproxies_natural_cells(struct engine *e) {
#ifdef WITH_MPI
  
  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the zoom domain */
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;

  /* Some info about the domain */
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[bkg_cell_offset].width[0],
                                cells[bkg_cell_offset].width[1],
                                cells[bkg_cell_offset].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Compute how many cells away we need to walk */
  int delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max * theta_crit_inv;
    delta_cells = (int)(distance / cells[bkg_cell_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (periodic) {
    if (delta_cells >= cdim[0] / 2) {
      if (cdim[0] % 2 == 0) {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2 - 1;
      } else {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2;
      }
    }
  } else {
    if (delta_cells > cdim[0]) {
      delta_m = cdim[0];
      delta_p = cdim[0];
    }
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k) + bkg_cell_offset;
        struct cell *ci = &cells[cid];

        /* Skip the void cell. */
        if (ci->tl_cell_type == void_tl_cell ||
            ci->tl_cell_type == void_tl_cell_neighbour) continue;

        /* Loop over every other cell within (Manhattan) range delta */
        for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

          /* Escape if non-periodic and beyond range */
          if (!periodic && (ii < 0 || ii >= cdim[0])) continue;

          for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

            /* Escape if non-periodic and beyond range */
            if (!periodic && (jj < 0 || jj >= cdim[1])) continue;
            
            for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

              /* Escape if non-periodic and beyond range */
              if (!periodic && (kk < 0 || kk >= cdim[2])) continue;

              /* Apply periodic BC (not harmful if not using periodic BC) */
              const int iii = (ii + cdim[0]) % cdim[0];
              const int jjj = (jj + cdim[1]) % cdim[1];
              const int kkk = (kk + cdim[2]) % cdim[2];

              /* Get the second cell */
              const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_cell_offset;
              struct cell *cj = &cells[cjd];

              /* Skip the void cell. */
              if (cj->tl_cell_type == void_tl_cell ||
                  cj->tl_cell_type == void_tl_cell_neighbour) continue;

              /* Avoid duplicates, and completely local and foreign
               * pairs */
              if (cid >= cjd ||
                  (ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;
              
              int proxy_type = 0;

              /* In the hydro case, only care about direct neighbours */
              if (with_hydro) {

                // MATTHIEU: to do: Write a better expression for the
                // non-periodic case.

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1)))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* First just add the direct neighbours. Then look for
                   some further out if the opening angle demands it */

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1))) {

                  proxy_type |= (int)proxy_cell_type_gravity;
                } else {

                  /* We don't have multipoles yet (or their CoMs) so we will
                     have to cook up something based on cell locations only. We
                     hence need a lower limit on the distance that the CoMs in
                     those cells could have and an upper limit on the distance
                     of the furthest particle in the multipole from its CoM.
                     We then can decide whether we are too close for an M2L
                     interaction and hence require a proxy as this pair of cells
                     cannot rely on just an M2L calculation. */

                  /* Minimal distance between any two points in the cells */
                  const double min_dist_CoM2 = cell_min_dist2_same_size(
                      &cells[cid], &cells[cjd], periodic, dim);

                  /* Are we beyond the distance where the truncated forces are 0
                   * but not too far such that M2L can be used? */
                  if (periodic) {

                    if ((min_dist_CoM2 < max_mesh_dist2) &&
                        !(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2))
                      proxy_type |= (int)proxy_cell_type_gravity;

                  } else {

                    if (!(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2)) {
                      proxy_type |= (int)proxy_cell_type_gravity;
                    }
                  }
                }
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cj->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cj->nodeID);

                  /* Store the information */
                  e->proxy_ind[cj->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                ci->mpi.sendto |= (1ULL << proxy_id);
              }

              /* Same for the symmetric case? */
              if (cj->nodeID == nodeID && ci->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[ci->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             ci->nodeID);

                  /* Store the information */
                  e->proxy_ind[ci->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cj->mpi.sendto |= (1ULL << proxy_id);
              }
            }
          }
        }
      }
    }
  }
#endif /* WITH_MPI */
}

/**
 * @brief Create and fill the proxies for the zoom cells.
 *
 * @param e The #engine.
 */
void engine_makeproxies_buffer_cells(struct engine *e) {
#ifdef WITH_MPI
  
  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

    /* Some info about the zoom domain */
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

  /* Some info about the domain */
  const int cdim[3] = {s->zoom_props->buffer_cdim[0],
                       s->zoom_props->buffer_cdim[1],
                       s->zoom_props->buffer_cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[buffer_cell_offset].width[0],
                                cells[buffer_cell_offset].width[1],
                                cells[buffer_cell_offset].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Compute how many cells away we need to walk */
  int delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max * theta_crit_inv;
    delta_cells = (int)(distance / cells[buffer_cell_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k) + buffer_cell_offset;
        struct cell *ci = &cells[cid];

        /* Skip the void cell. */
        if (ci->tl_cell_type == void_tl_cell) continue;

        /* Loop over every other cell within (Manhattan) range delta */
        for (int ii = i - delta_m; ii <= i + delta_p; ii++) {
          
          /* Buffer cells are never periodic. */
          if (ii < 0 || ii >= cdim[0]) continue;

          for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

            /* Buffer cells are never periodic. */
            if (jj < 0 || jj >= cdim[1]) continue;

            for (int kk = k - delta_m; kk <= k + delta_p; kk++) {
              
              /* Buffer cells are never periodic. */
              if (kk < 0 || kk >= cdim[2]) continue;

              /* Get the second cell */
              const int cjd = cell_getid(cdim, ii, jj, kk) + buffer_cell_offset;
              struct cell *cj = &cells[cjd];

              /* Skip the void cell. */
              if (cj->tl_cell_type == void_tl_cell) continue;
              
              /* Avoid duplicates, and completely local and
               * foreign pairs */
              if (cid >= cjd ||
                  (ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              int proxy_type = 0;

              /* In the hydro case, only care about direct neighbours */
              if (with_hydro) {

                /* Check for direct neighbours without periodic BC */
                if ((abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
                    (abs(k - kk) <= 1))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* First just add the direct neighbours. Then look for
                   some further out if the opening angle demands it */

                /* Check for direct neighbours without periodic BC */
                if ((abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
                    (abs(k - kk) <= 1)) {

                  proxy_type |= (int)proxy_cell_type_gravity;

                } else {

                  /* We don't have multipoles yet (or their CoMs) so we will
                     have to cook up something based on cell locations only. We
                     hence need a lower limit on the distance that the CoMs in
                     those cells could have and an upper limit on the distance
                     of the furthest particle in the multipole from its CoM.
                     We then can decide whether we are too close for an M2L
                     interaction and hence require a proxy as this pair of cells
                     cannot rely on just an M2L calculation. */

                  /* Minimal distance between any two points in the cells */
                  const double min_dist_CoM2 = cell_min_dist2_same_size(
                      &cells[cid], &cells[cjd], periodic, dim);

                  /* Are we beyond the distance where the truncated forces are 0
                   * but not too far such that M2L can be used? */
                  if (periodic) {

                    if ((min_dist_CoM2 < max_mesh_dist2) &&
                        !(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2))
                      proxy_type |= (int)proxy_cell_type_gravity;

                  } else {

                    if (!(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2)) {
                      proxy_type |= (int)proxy_cell_type_gravity;
                    }
                  }
                }
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cj->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cj->nodeID);

                  /* Store the information */
                  e->proxy_ind[cj->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                ci->mpi.sendto |= (1ULL << proxy_id);
              }

              /* Same for the symmetric case? */
              if (cj->nodeID == nodeID && ci->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[ci->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             ci->nodeID);

                  /* Store the information */
                  e->proxy_ind[ci->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cj->mpi.sendto |= (1ULL << proxy_id);
              }
            }
          }
        }
      }
    }
  }
#endif
}

/**
 * @brief Create and fill the proxies for the zoom cells.
 *
 * @param e The #engine.
 */
void engine_makeproxies_zoom_cells(struct engine *e) {
#ifdef WITH_MPI
  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const int cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                       s->zoom_props->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Compute how many cells away we need to walk */
  int delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max * theta_crit_inv;
    delta_cells = (int)(distance / cells[0].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);
        struct cell *ci = &cells[cid];
        
        /* Loop over every other cell within (Manhattan) range delta */
        for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

          /* Zoom cells are never periodic, exit if beyond zoom region */
          if (ii < 0 || ii >= cdim[0]) continue;

          for (int jj = j - delta_m; jj <= j + delta_p; jj++) {
            
            /* Zoom cells are never periodic, exit if beyond zoom region */
            if (jj < 0 || jj >= cdim[1]) continue;
            
            for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

              /* Zoom cells are never periodic, exit if beyond zoom region */
              if (kk < 0 || kk >= cdim[2]) continue;
          
              /* Get the second cell */
              const int cjd = cell_getid(cdim, ii, jj, kk);
              struct cell *cj = &cells[cjd];

              /* Avoid duplicates, and completely local and
               * foreign pairs */
              if (cid >= cjd || 
                  (ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              int proxy_type = 0;

              /* In the hydro case, only care about direct neighbours */
              if (with_hydro) {

                /* Check for direct neighbours without periodic BC */
                if ((abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
                    (abs(k - kk) <= 1))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* First just add the direct neighbours. Then look for
                   some further out if the opening angle demands it */

                /* Check for direct neighbours without periodic BC */
                if ((abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
                    (abs(k - kk) <= 1)) {

                  proxy_type |= (int)proxy_cell_type_gravity;

                } else {

                  /* We don't have multipoles yet (or their CoMs) so we will
                     have to cook up something based on cell locations only. We
                     hence need a lower limit on the distance that the CoMs in
                     those cells could have and an upper limit on the distance
                     of the furthest particle in the multipole from its CoM.
                     We then can decide whether we are too close for an M2L
                     interaction and hence require a proxy as this pair of cells
                     cannot rely on just an M2L calculation. */

                  /* Minimal distance between any two points in the cells */
                  const double min_dist_CoM2 = cell_min_dist2_same_size(
                      &cells[cid], &cells[cjd], periodic, dim);

                  /* Are we beyond the distance where the truncated forces are 0
                   * but not too far such that M2L can be used? */
                  if (periodic) {

                    if ((min_dist_CoM2 < max_mesh_dist2) &&
                        !(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2))
                      proxy_type |= (int)proxy_cell_type_gravity;

                  } else {

                    if (!(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2)) {
                      proxy_type |= (int)proxy_cell_type_gravity;
                    }
                  }
                }
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cj->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cj->nodeID);

                  /* Store the information */
                  e->proxy_ind[cj->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                ci->mpi.sendto |= (1ULL << proxy_id);
              }

              /* Same for the symmetric case? */
              if (cj->nodeID == nodeID && ci->nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[ci->nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             ci->nodeID);

                  /* Store the information */
                  e->proxy_ind[ci->nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cj->mpi.sendto |= (1ULL << proxy_id);
              }
            }
          }
        }
      }
    }
  }
#endif
}

/**
 * @brief Create and fill the proxies for relations between zoom and bkg cell
 *        grids.
 *
 * NOTE: If we have a buffer region there are no relations between these grids.
 *
 * This is done "lazily" by just making proxies for all neighbour cells
 * as these are defined to be within the gravity criterion.
 *
 * @param e The #engine.
 */
void engine_makeproxies_between_zoom_bkg(struct engine *e) {
#ifdef WITH_MPI

  /* Useful local information */
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Gravity information */
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Get the neighbouring background cells and their offset. */
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;
  const int nr_neighbours = s->zoom_props->nr_neighbour_cells;
  const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[neighbour_cells[0]].width[0],
                                cells[neighbour_cells[0]].width[1],
                                cells[neighbour_cells[0]].width[2]};

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Loop over each zoom cell in the space. */
  for (int cid = 0; cid < bkg_cell_offset; cid++) {

    /* Get the cell. */
    struct cell *ci = &s->cells_top[cid];
    
    /* Loop over every neighbouring background cells */
    for (int k = 0; k < nr_neighbours; k++) {

      /* Get the cell */
      int cjd = neighbour_cells[k];
      struct cell *cj = &s->cells_top[cjd];

      /* Skip the void cell and neighbour background cells. */
      if (cj->tl_cell_type == void_tl_cell ||
          cj->tl_cell_type == void_tl_cell_neighbour) continue;

      /* Early abort (both same node) */
      if (ci->nodeID == nodeID && cj->nodeID == nodeID) continue;

      /* Early abort (both foreign node) */
      if (ci->nodeID != nodeID && cj->nodeID != nodeID) continue;

      int proxy_type = 0;

      /* Minimal distance between any two points in the cells */
      const double min_dist_CoM2 =
        cell_min_dist2_diff_size(&cells[cid], &cells[cjd], periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0
       * but not too far such that M2L can be used? */
      if (periodic) {

        if ((min_dist_CoM2 < max_mesh_dist2) &&
            !(4. * r_max * r_max <
              theta_crit * theta_crit * min_dist_CoM2))
          proxy_type |= (int)proxy_cell_type_gravity;
        
      } else {

        if (!(4. * r_max * r_max <
              theta_crit * theta_crit * min_dist_CoM2)) {
          proxy_type |= (int)proxy_cell_type_gravity;
        }
      }

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Add to proxies? */
      if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[cj->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");

          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     cj->nodeID);
          
          /* Store the information */
          e->proxy_ind[cj->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);
        
        /* Store info about where to send the cell */
        ci->mpi.sendto |= (1ULL << proxy_id);
      }
      
      /* Same for the symmetric case? */
      if (cj->nodeID == nodeID && ci->nodeID != nodeID) {
        
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[ci->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");
          
          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     ci->nodeID);
          
          /* Store the information */
          e->proxy_ind[ci->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);
        
        /* Store info about where to send the cell */
        cj->mpi.sendto |= (1ULL << proxy_id);
      }
    }
  }
#endif
}

/**
 * @brief Create and fill the proxies for relations between buffer and bkg cell
 *        grids.
 *
 * @param e The #engine.
 */
void engine_makeproxies_between_buffer_bkg(struct engine *e) {
#ifdef WITH_MPI

  /* Useful local information */
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Gravity information */
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Get the cell offsets. */
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[bkg_cell_offset].width[0],
                                cells[bkg_cell_offset].width[1],
                                cells[bkg_cell_offset].width[2]};

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Loop over each buffer cell in the space. */
  for (int cid = buffer_cell_offset; cid < s->nr_cells; cid++) {

    /* Get the cell. */
    struct cell *ci = &s->cells_top[cid];

    /* Skip the void cells. */
    if (ci->tl_cell_type == void_tl_cell) continue;
    
    /* Loop over every neighbouring background cells */
    for (int cjd = bkg_cell_offset; cjd < buffer_cell_offset; cjd++) {

      /* Get the cell */
      struct cell *cj = &s->cells_top[cjd];
      
      /* Skip the void cell and neighbour background cells. */
      if (cj->tl_cell_type == void_tl_cell ||
          cj->tl_cell_type == void_tl_cell_neighbour) continue;

      /* Early abort (both same node) */
      if (ci->nodeID == nodeID && cj->nodeID == nodeID) continue;

      /* Early abort (both foreign node) */
      if (ci->nodeID != nodeID && cj->nodeID != nodeID) continue;

      int proxy_type = 0;

      /* Minimal distance between any two points in the cells */
      const double min_dist_CoM2 =
        cell_min_dist2_diff_size(&cells[cid], &cells[cjd], periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0
       * but not too far such that M2L can be used? */
      if (periodic) {

        if ((min_dist_CoM2 < max_mesh_dist2) &&
            !(4. * r_max * r_max <
              theta_crit * theta_crit * min_dist_CoM2))
          proxy_type |= (int)proxy_cell_type_gravity;
        
      } else {

        if (!(4. * r_max * r_max <
              theta_crit * theta_crit * min_dist_CoM2)) {
          proxy_type |= (int)proxy_cell_type_gravity;
        }
      }

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Add to proxies? */
      if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[cj->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");

          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     cj->nodeID);
          
          /* Store the information */
          e->proxy_ind[cj->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);
        
        /* Store info about where to send the cell */
        ci->mpi.sendto |= (1ULL << proxy_id);
      }
      
      /* Same for the symmetric case? */
      if (cj->nodeID == nodeID && ci->nodeID != nodeID) {
        
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[ci->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");
          
          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     ci->nodeID);
          
          /* Store the information */
          e->proxy_ind[ci->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);
        
        /* Store info about where to send the cell */
        cj->mpi.sendto |= (1ULL << proxy_id);
      }
    }
  }
#endif
}

/**
 * @brief Create and fill the proxies including the zoom region.
 *
 * This replaces the function in engine_proxy when running with a zoom region.
 *
 * @param e The #engine.
 */
void engine_makeproxies_with_zoom_region(struct engine *e) {
#ifdef WITH_MPI

  /* Let's time this */
  const ticks tic = getticks();

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Useful local information */
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Gravity information */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;

  /* Set up some width and distance variables. */
  double r_diag2, r_diag, r_max;

  /* Loop over each zoom cell in the space. */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Get the cell. */
    struct cell *ci = &s->cells_top[cid];

    /* Skip the void cells. */
    if (ci->tl_cell_type == void_tl_cell ||
        ci->tl_cell_type == void_tl_cell_neighbour) continue;

    /* Get the ijk coordinates */
    int i, j, k;
    if (cid < s->zoom_props->nr_zoom_cells) {
      i = cid / (s->zoom_props->cdim[1] * s->zoom_props->cdim[2]);
      j = (cid / s->zoom_props->cdim[2]) % s->zoom_props->cdim[1];
      k = cid % s->zoom_props->cdim[2];
    } else if (cid < s->zoom_props->nr_zoom_cells +
               s->zoom_props->nr_bkg_cells) {
      i = (cid - s->zoom_props->tl_cell_offset) / (s->cdim[1] * s->cdim[2]);
      j = ((cid - s->zoom_props->tl_cell_offset) / s->cdim[2]) % s->cdim[1];
      k = (cid - s->zoom_props->tl_cell_offset) % s->cdim[2];
    } else {
      i = (cid - s->zoom_props->buffer_cell_offset) /
        (s->zoom_props->buffer_cdim[1] * s->zoom_props->buffer_cdim[2]);
      j =
        ((cid - s->zoom_props->buffer_cell_offset) /
         s->zoom_props->buffer_cdim[2]) % s->zoom_props->buffer_cdim[1];
      k = (cid - s->zoom_props->buffer_cell_offset) %
        s->zoom_props->buffer_cdim[2];
    }

    /* Distance between centre of the cell and corners */
    r_diag2 =
      ci->width[0] * ci->width[0] +
      ci->width[1] * ci->width[1] +
      ci->width[2] * ci->width[2];
    r_diag = 0.5 * sqrt(r_diag2);

    /* Maximal distance from shifted CoM to any corner */
    r_max = 2 * r_diag;

    /* Loop over every other cell avoiding duplicates. */
    for (int cjd = 0; cjd < cid; cjd++) {

      /* Get the cell. */
      struct cell *cj = &cells[cjd];

      /* Avoid completely local and foreign pairs */
      if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
          (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      /* Skip the void cells. */
      if (cj->tl_cell_type == void_tl_cell ||
          cj->tl_cell_type == void_tl_cell_neighbour) continue;

      /* Get the ijk coordinates */
      int ii, jj, kk;
      if (cjd < s->zoom_props->nr_zoom_cells) {
        ii = cjd / (s->zoom_props->cdim[1] * s->zoom_props->cdim[2]);
        jj = (cjd / s->zoom_props->cdim[2]) % s->zoom_props->cdim[1];
        kk = cjd % s->zoom_props->cdim[2];
      } else if (cjd < s->zoom_props->nr_zoom_cells +
                 s->zoom_props->nr_bkg_cells) {
        ii = (cjd - s->zoom_props->tl_cell_offset) / (s->cdim[1] * s->cdim[2]);
        jj = ((cjd - s->zoom_props->tl_cell_offset) / s->cdim[2]) % s->cdim[1];
        kk = (cjd - s->zoom_props->tl_cell_offset) % s->cdim[2];
      } else {
        ii = (cjd - s->zoom_props->buffer_cell_offset) /
          (s->zoom_props->buffer_cdim[1] * s->zoom_props->buffer_cdim[2]);
        jj =
          ((cjd - s->zoom_props->buffer_cell_offset) /
           s->zoom_props->buffer_cdim[2]) % s->zoom_props->buffer_cdim[1];
        kk = (cjd - s->zoom_props->buffer_cell_offset) %
          s->zoom_props->buffer_cdim[2];
      }

      /* Compute the neighbours distance between centre of the cell
       * and corners */
      r_diag2 =
        cj->width[0] * cj->width[0] +
        cj->width[1] * cj->width[1] +
        cj->width[2] * cj->width[2];
      r_diag = 0.5 * sqrt(r_diag2);

      /* Use the largest r_max. */
      r_max = max(r_max, 2 * r_diag);

      int proxy_type = 0;

      /* In the hydro case, only care about direct neighbours */
      if (with_hydro && ci->tl_cell_type == zoom_tl_cell &&
          ci->tl_cell_type == cj->tl_cell_type) {

        /* Check for direct neighbours without periodic BC */
        if ((abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
            (abs(k - kk) <= 1))
          proxy_type |= (int)proxy_cell_type_hydro;
      }

      /* In the gravity case, check distances using the MAC. */
      if (with_gravity) {
        
        /* First just add the direct neighbours. Then look for
           some further out if the opening angle demands it */
        
        /* Check for direct neighbours without periodic BC */
        if (ci->tl_cell_type == cj->tl_cell_type &&
            (abs(i - ii) <= 1) && (abs(j - jj) <= 1) &&
            (abs(k - kk) <= 1)) {
          
          proxy_type |= (int)proxy_cell_type_gravity;
          
        } else {

          /* We don't have multipoles yet (or their CoMs) so we will
             have to cook up something based on cell locations only. We
             hence need a lower limit on the distance that the CoMs in
             those cells could have and an upper limit on the distance
             of the furthest particle in the multipole from its CoM.
             We then can decide whether we are too close for an M2L
             interaction and hence require a proxy as this pair of cells
             cannot rely on just an M2L calculation. */

          /* Minimal distance between any two points in the cells */
          const double min_dist_CoM2 =
            cell_min_dist2_same_size(&cells[cid], &cells[cjd], periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0
           * but not too far such that M2L can be used? */
          if (periodic) {
            
            if ((min_dist_CoM2 < max_mesh_dist2) &&
                !(4. * r_max * r_max <
                  theta_crit * theta_crit * min_dist_CoM2))
              proxy_type |= (int)proxy_cell_type_gravity;
            
          } else {
            
            if (!(4. * r_max * r_max <
                  theta_crit * theta_crit * min_dist_CoM2)) {
              proxy_type |= (int)proxy_cell_type_gravity;
            }
          }
        }
      }

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;
      
      /* Add to proxies? */
      if (ci->nodeID == nodeID && cj->nodeID != nodeID) {
        
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[cj->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");
          
          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     cj->nodeID);
          
          /* Store the information */
          e->proxy_ind[cj->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);
        
        /* Store info about where to send the cell */
        ci->mpi.sendto |= (1ULL << proxy_id);
      }
      
      /* Same for the symmetric case? */
      if (cj->nodeID == nodeID && ci->nodeID != nodeID) {
        
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[ci->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");
          
          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     ci->nodeID);
          
          /* Store the information */
          e->proxy_ind[ci->nodeID] = e->nr_proxies;
          proxy_id = e->nr_proxies;
          e->nr_proxies += 1;
          
          /* Check the maximal proxy limit */
          if ((size_t)proxy_id > 8 * sizeof(long long))
            error(
                  "Created more than %zd proxies. cell.mpi.sendto will "
                  "overflow.",
                  8 * sizeof(long long));
        }
        
        /* Add the cell to the proxy */
        proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
        proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);
        
        /* Store info about where to send the cell */
        cj->mpi.sendto |= (1ULL << proxy_id);
      }
    }
  }
  
  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
