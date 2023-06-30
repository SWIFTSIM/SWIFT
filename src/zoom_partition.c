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

  /* Define the old bkg_cell_offset */
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
            cell_getid(s->cdim, i, j, k) + s->zoom_props->bkg_cell_offset;
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
 *        for edge counting, adjancency, and edges weighting. This function is
 *        used when wedges are used for the background deomposition.
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
void wedge_edge_loop(const int *cdim, int offset, struct space *s,
                     idx_t *adjncy, idx_t *xadj, double *counts, double *edges,
                     int *iedge) {

  /* The number of slices in theta. */
  int theta_nslices = s->zoom_props->theta_nslices;
  int phi_nslices = s->zoom_props->phi_nslices;

  /* The number of zoom cells. */
  int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Declare some variables. */
  struct cell *restrict ci;
  struct cell *restrict cj;

  /* Loop over the provided cells and find their edges. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
          
        /* Get the cell index. */
        const size_t cid = cell_getid(cdim, i, j, k);

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
          if (ii < 0 || ii >= cdim[0]) continue;
          for (int jj = j - 1; jj <= j + 1; jj++) {
            if (jj < 0 || jj >= cdim[1]) continue;
            for (int kk = k - 1; kk <= k + 1; kk++) {
              if (kk < 0 || kk >= cdim[2]) continue;

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

        /* /\* Which wedge is this zoom cell in? *\/ */
        /* int wedge_ind = get_wedge_index(s, ci); */

        /* Now loop over the wedges. */
        for (int ii = 0; ii < theta_nslices; ii++) {
          for (int jj = 0; jj < phi_nslices; jj++) {

            /* Find the wedge index. */
            const int jwedge_ind = ii * phi_nslices + jj;

            /* Handle size_to_edges case */
            if (edges != NULL) {
              /* Store this edge. */
              edges[*iedge] = counts[nr_zoom_cells + jwedge_ind];
              (*iedge)++;
            }
                
            /* Handle graph_init case */
            else if (adjncy != NULL) {
              adjncy[*iedge] = nr_zoom_cells + jwedge_ind;
              (*iedge)++;
            }

            /* Handle find_vertex_edges case */
            else {
              /* Record an edge. */
              ci->nr_vertex_edges++;
              (*iedge)++;
            }
          }
        }
      }
    }
  }

  /* Now loop over the wedges. */
  for (int i = 0; i < theta_nslices; i++) {
    for (int j = 0; j < phi_nslices; j++) {

      /* Find the wedge index. */
      const int iwedge_ind = i * phi_nslices + j;

      /* Define the current "cell" index. */
      int cid = nr_zoom_cells + iwedge_ind;

      
#ifdef SWIFT_DEBUG_CHECKS
        if (xadj != NULL) {
             
          /* Ensure the previous cell has found enough edges. */
          if ((iwedge_ind > 0) && ((*iedge - xadj[cid - 1]) != s->zoom_props->nr_wedge_edges[iwedge_ind - 1]))
            error("Found too few edges (nedges=%ld, wedge->nr_vertex_edges=%d)",
                  *iedge - xadj[cid - 1], s->zoom_props->nr_wedge_edges[iwedge_ind - 1]);
          
        }
#endif

      /* If given set METIS xadj. */
      if (xadj != NULL) {
        xadj[cid] = *iedge;
        
        /* Set edges start pointer for this wedge. */
        s->zoom_props->wedge_edges_start[iwedge_ind] = *iedge;
      }

      /* Loop over neighbouring cells */
      for (int ii = i - 1; ii <= i + 1; ii++) {
        for (int jj = j - 1; jj <= j + 1; jj++) {
          
          /* Wrap the indices around the sphere. */
          const int iii = (ii + theta_nslices) % theta_nslices;
          const int jjj = (jj + phi_nslices) % phi_nslices;

          /* Find the wedge index. */
          const int jwedge_ind = iii * phi_nslices + jjj;

          /* Skip self. */
          if (iwedge_ind == jwedge_ind) continue;

          /* Handle size_to_edges case */
          if (edges != NULL) {
            /* Store this edge. */
            edges[*iedge] = counts[nr_zoom_cells + jwedge_ind];
            (*iedge)++;
          }
          
          /* Handle graph_init case */
          else if (adjncy != NULL) {
            adjncy[*iedge] = nr_zoom_cells + jwedge_ind;
            (*iedge)++;
          }
          
          /* Handle find_vertex_edges case */
          else {
            /* Record an edge. */
            s->zoom_props->nr_wedge_edges[iwedge_ind]++;
            (*iedge)++;
          }
        }
      }

      /* Now find the zoom cell edges for this wedge. */
      for (int zoom_ii = 0; zoom_ii < cdim[0]; zoom_ii++) {
        for (int zoom_jj = 0; zoom_jj < cdim[1]; zoom_jj++) {
          for (int zoom_kk = 0; zoom_kk < cdim[2]; zoom_kk++) {
            
            /* Get cell ID. */
            const int cjd = cell_getid(cdim, zoom_ii, zoom_jj, zoom_kk);
            
            /* Get the cell */
            cj = &s->cells_top[cjd];
            
            /* /\* Get the wedge index of this cell. *\/ */
            /* int jwedge_ind = get_wedge_index(s, cj); */
            
            /* /\* Skip if not in this wedge. *\/ */
            /* if (iwedge_ind != jwedge_ind) continue; */
            
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
              /* Record an edge. */
              s->zoom_props->nr_wedge_edges[iwedge_ind]++;
              (*iedge)++;
            } 
          }
        }
      }
    }
  }
}
#endif

/**
 * @brief A genric looping function to handle duplicated looping methods done
 *        for edge counting, adjancency, and edges weighting. This function is
 *        used when treating each individual grid individually.
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
void simple_edge_loop(const int *cdim, int offset, struct space *s,
                      idx_t *adjncy, idx_t *xadj, double *counts, double *edges,
                      int *iedge) {

  /* Declare some variables. */
  struct cell *restrict ci;
  struct cell *restrict cj;

  /* Loop over the provided cells and find their edges. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
          
        /* Get the cell index. */
        const size_t cid = cell_getid(cdim, i, j, k) + offset;

        /* Get the vertex index of this cell. */
        const int vid = cid - offset;

        /* Get the cell. */
        ci = &s->cells_top[cid];

#ifdef SWIFT_DEBUG_CHECKS
        if (xadj != NULL) {
             
          /* Ensure the previous cell has found enough edges. */
          if ((cid > 0) &&
              ((*iedge - xadj[vid - 1]) != s->cells_top[vid - 1].nr_vertex_edges))
            error("Found too few edges (nedges=%ld, c->nr_vertex_edges=%d)",
                  *iedge - xadj[vid - 1], s->cells_top[vid - 1].nr_vertex_edges);
          
        }
#endif

        /* If given set METIS xadj. */
        if (xadj != NULL) {
          xadj[vid] = *iedge;
          
          /* Set edges start pointer for this cell. */
          ci->edges_start = *iedge;
        }

        /* Loop over a shell of cells with the same type. */
        for (int ii = i - 1; ii <= i + 1; ii++) {
          if (ii < 0 || ii >= cdim[0]) continue;
          for (int jj = j - 1; jj <= j + 1; jj++) {
            if (jj < 0 || jj >= cdim[1]) continue;
            for (int kk = k - 1; kk <= k + 1; kk++) {
              if (kk < 0 || kk >= cdim[2]) continue;

              /* Apply periodic BC (not harmful if not using periodic BC) */
              const int iii = (ii + cdim[0]) % cdim[0];
              const int jjj = (jj + cdim[1]) % cdim[1];
              const int kkk = (kk + cdim[2]) %cdim[2];
                
              /* Get cell index. */
              const size_t cjd = cell_getid(cdim, iii, jjj, kkk) + offset;

              /* Get the vertex index of this cell. */
              const int vjd = cjd - offset;
              
              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Skip self. */
              if (cid == cjd) continue;
                
              /* Handle size_to_edges case */
              if (edges != NULL) {
                /* Store this edge. */
                edges[*iedge] = counts[vjd];
                (*iedge)++;
              }
                
              /* Handle graph_init case */
              else if (adjncy != NULL) {
                adjncy[*iedge] = vjd;
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
      } /* i loop */
    } /* j loop */
  } /* k loop */
}
#endif

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

  /* If we are running with wedges call that edge loop function. */
  if (s->zoom_props->use_bkg_wedges) {
    wedge_edge_loop(cdim, offset, s, adjncy, xadj, counts, edges, iedge);
    return;
  }
  
  /* Are we doing each grid separately? */
  else if (s->zoom_props->separate_decomps) {
        simple_edge_loop(cdim, offset, s, adjncy, xadj, counts, edges, iedge);
    return;
  }

  /* Some info about the domain */
  const int periodic = s->periodic;
  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Set up cell offsets. */
  const int bkg_offset = s->zoom_props->bkg_cell_offset ;
  const int buff_offset = s->zoom_props->buffer_cell_offset;

  /* Define the looping bounds. */
  int delta_m = 1;
  int delta_p = 1;

  /* Define cell variables we will need. */
  struct cell *restrict ci;
  struct cell *restrict cj;

  /* Get region boundaries. */
  const double  buffer_bounds[6] = {
    s->zoom_props->buffer_bounds[0], s->zoom_props->buffer_bounds[1],
    s->zoom_props->buffer_bounds[2], s->zoom_props->buffer_bounds[3],
    s->zoom_props->buffer_bounds[4], s->zoom_props->buffer_bounds[5]};

  /* ======================= Start with zoom cells ======================= */

  /* Get cdim. */
  cdim = s->zoom_props->cdim;

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell index. */
        const size_t cid = cell_getid(cdim, i, j, k);

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
        
        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;

              /* Get cell index. */
              const size_t cjd = cell_getid(cdim, iii, jjj, kkk);
              
              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Skip self. */
              if (cid == cjd) continue;

              /* Will we need a task here? Uses geometric criterion. */
              if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
                continue;
                
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
            }
          }
        }

        /* And now loop over neighbour background/buffer cells. */

        /* Get the neighbouring background cells. */
        const int nr_neighbours = s->zoom_props->nr_neighbour_cells;
        const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

        /* Now loop over the neighbouring background cells.  */
        for (int n = 0; n < nr_neighbours; n++) {

          /* Get the cell index of this neighbour. */
          int cjd = neighbour_cells[n];

          /* Handle on the neighbouring background cell. */
          cj = &s->cells_top[cjd];
                
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
        }
      }
    }
  }

  /* ======================= Now buffer cells ======================= */

  if (s->zoom_props->with_buffer_cells) {

    /* Get cdim. */
    cdim = s->zoom_props->buffer_cdim;

    /* Loop over each cell in the space. */
    for (int i = 0; i < cdim[0]; i++) {
      for (int j = 0; j < cdim[1]; j++) {
        for (int k = 0; k < cdim[2]; k++) {

          /* Get the cell ID. */
          const int cid = cell_getid(cdim, i, j, k) + buff_offset;

          /* Get the cell. */
          ci = &s->cells_top[cid];

          /* Skip the void cell. */
          if (ci->subtype == void_cell) continue;

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

          /* Loop over all its neighbours in range. */
          for (int ii = -delta_m; ii <= delta_p; ii++) {
            int iii = i + ii;
            if (iii < 0 || iii >= cdim[0]) continue;
            for (int jj = -delta_m; jj <= delta_p; jj++) {
              int jjj = j + jj;
              if (jjj < 0 || jjj >= cdim[1]) continue;
              for (int kk = -delta_m; kk <= delta_p; kk++) {
                int kkk = k + kk;
                if (kkk < 0 || kkk >= cdim[2]) continue;

                /* Get the cell ID. */
                const int cjd = cell_getid(cdim, iii, jjj, kkk) + buff_offset;

                /* Get the cell. */
                cj = &s->cells_top[cjd];

                /* Skip self. */
                if (cid == cjd) continue;

                if (cj->subtype == void_cell) continue;

                /* Will we need a task here? Uses geometric criterion. */
                if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
                  continue;
                
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
              }
            }
          }

          /* If this is a neighbour cell it must we need to find edges with zoom
           * cells. */
          if (ci->subtype == neighbour) {

            /* Loop over the zoom cells. */
            for (int cjd = 0; cjd < nr_zoom_cells; cjd++) {

              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Get the (i,j,k) location of the zoom cell in the buffer grid. */
              int ii =
                (cj->loc[0] - buffer_bounds[0]) * s->zoom_props->buffer_iwidth[0];
              int jj =
                (cj->loc[1] - buffer_bounds[2]) * s->zoom_props->buffer_iwidth[1];
              int kk =
                (cj->loc[2] - buffer_bounds[4]) * s->zoom_props->buffer_iwidth[2];
              
              /* Will we need a task here? Uses geometric criterion. */
              if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
                continue;

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
            }
          }

          /* Now we need to check which background cells for edges. */
          for (int cjd = bkg_offset; cjd < buff_offset; cjd++) {
            
            /* Get the cell. */
            cj = &s->cells_top[cjd];

            if (cj->subtype == empty) continue;

            /* Get the (i,j,k) location of the buffer cell in the background
             * grid. */
            int ii = ci->loc[0] * s->iwidth[0];
            int jj = ci->loc[1] * s->iwidth[1];
            int kk = ci->loc[2] * s->iwidth[2];

            /* Get the (i,j,k) location of the background cell. */
            int iii = cj->loc[0] * s->iwidth[0];
            int jjj = cj->loc[1] * s->iwidth[1];
            int kkk = cj->loc[2] * s->iwidth[2];
            
            /* Will we need a task here? Uses geometric criterion. */
            if (abs(ii - iii) <= 1 && abs(jj - jjj) <= 1 && abs(kk - kkk) <= 1)
              continue;
            
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
          }
        }
      }
    }
  }

  /* ======================= Now background cells ======================= */

  /* Get cdim. */
  cdim = s->cdim;

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k) + bkg_offset;

        /* Get the cell. */
        ci = &s->cells_top[cid];

        /* Skip the void cell. */
        if (ci->subtype == void_cell || ci->subtype == empty) continue;

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
        
        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_offset;

              /* Early abort  */
              if (cid == cjd) continue;

              /* Skip the void cell. */
              if (cj->subtype == void_cell || cj->subtype == empty) continue;
              
              /* Get the cell. */
              cj = &s->cells_top[cjd];

              /* Will we need a task here? Uses geometric criterion. */
              if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
                continue;
              
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
            }
          }
        }

        /* If this is a neighbour cell it must we need to find edges with zoom
         * cells. */
        if (ci->subtype == neighbour) {

          /* Loop over the zoom cells. */
          for (int cjd = 0; cjd < nr_zoom_cells; cjd++) {

            /* Get the cell. */
            cj = &s->cells_top[cjd];

            /* Get the (i,j,k) location of the zoom cell in the background
             * grid. */
            int ii = (cj->loc[0] + (cj->width[0] / 2)) * s->iwidth[0];
            int jj = (cj->loc[1] + (cj->width[1] / 2)) * s->iwidth[1];
            int kk = (cj->loc[2] + (cj->width[2] / 2)) * s->iwidth[2];
              
            /* Will we need a task here? Uses geometric criterion. */
            if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
              continue;

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
          }
        }

        /* If we have buffer cells we need to find edges with buffer cells. */
        if (s->zoom_props->with_buffer_cells) {

          /* Loop over the zoom cells. */
          for (int cjd = bkg_offset; cjd < s->nr_cells; cjd++) {

            /* Get the cell. */
            cj = &s->cells_top[cjd];

            /* Skip the void cell */
            if (cj->subtype == void_cell) continue;

            /* Get the (i,j,k) location of the buffer cell in the background
             * grid. */
            int ii = (cj->loc[0] + (cj->width[0] / 2)) * s->iwidth[0];
            int jj = (cj->loc[1] + (cj->width[1] / 2)) * s->iwidth[1];
            int kk = (cj->loc[2] + (cj->width[2] / 2)) * s->iwidth[2];
              
            /* Will we need a task here? Uses geometric criterion. */
            if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
              continue;

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
          }
        }
      }
    }
  }
}
#endif

#ifdef WITH_MPI
/**
 * @brief What type of proxy do we need between these cells?
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param i integer coordinate of ci.
 * @param j integer coordinate of ci.
 * @param k integer coordinate of ci.
 * @param ii integer coordinate of cj.
 * @param jj integer coordinate of cj.
 * @param kk integer coordinate of cj.
 * @param r_max The minimum separation of ci and cj.
 */
int find_proxy_type(struct cell *ci, struct cell *cj, struct engine *e,
                    int i, int j, int k, int ii, int jj, int kk,
                    double r_max, const double *dim, const int periodic) {

  /* Gravity information */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;
  
  int proxy_type = 0;

  /* In the hydro case, only care about direct neighbours of the same
   * type. */
  if (with_hydro && ci->type == zoom && ci->type == cj->type) {

    /* Check for direct neighbours without periodic BC */
    if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
      proxy_type |= (int)proxy_cell_type_hydro;
  }

  /* In the gravity case, check distances using the MAC. */
  if (with_gravity) {

    /* We don't have multipoles yet (or their CoMs) so we will
       have to cook up something based on cell locations only. We
       hence need a lower limit on the distance that the CoMs in
       those cells could have and an upper limit on the distance
       of the furthest particle in the multipole from its CoM.
       We then can decide whether we are too close for an M2L
       interaction and hence require a proxy as this pair of cells
       cannot rely on just an M2L calculation. */
    
    /* Minimal distance between any two points in the cells. */
    const double min_dist_CoM2 = cell_min_dist2(ci, cj, periodic, dim);
    
    /* Are we beyond the distance where the truncated forces are 0
     * but not too far such that M2L can be used? */
    if (periodic) {
      
      if ((min_dist_CoM2 < max_mesh_dist2) &&
          !(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2))
        proxy_type |= (int)proxy_cell_type_gravity;
      
    } else {
      
      if (!(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2)) {
        proxy_type |= (int)proxy_cell_type_gravity;
      }
    }
  }

  return proxy_type;
}

/**
 * @brief What type of proxy do we need between these cells?
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param e The #engine.
 * @param proxies The proxies themselves.
 * @param nodeID What rank is this?
 * @param proxy_type What sort of proxy is this?
 */
void add_proxy(struct cell *ci, struct cell *cj, struct engine *e,
               struct proxy *proxies, const int nodeID, const int proxy_type) {

  /* Add to proxies? */
  if (ci->nodeID == nodeID && cj->nodeID != nodeID) {
        
    /* Do we already have a relationship with this node? */
    int proxy_id = e->proxy_ind[cj->nodeID];
    if (proxy_id < 0) {
      if (e->nr_proxies == engine_maxproxies)
        error("Maximum number of proxies exceeded.");
      
      /* Ok, start a new proxy for this pair of nodes */
      proxy_init(&proxies[e->nr_proxies], nodeID,
                 cj->nodeID);
          
      /* Store the information */
      e->proxy_ind[cj->nodeID] = e->nr_proxies;
      proxy_id = e->nr_proxies;
      e->nr_proxies += 1;
          
      /* Check the maximal proxy limit */
      if ((size_t)proxy_id > 8 * sizeof(long long))
        error("Created more than %zd proxies. cell.mpi.sendto will "
              "overflow.",
              8 * sizeof(long long));
    }
        
    /* Add the cell to the proxy */
    proxy_addcell_in(&proxies[proxy_id], cj, proxy_type);
    proxy_addcell_out(&proxies[proxy_id], ci, proxy_type);
        
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
        error("Created more than %zd proxies. cell.mpi.sendto will "
              "overflow.",
              8 * sizeof(long long));
    }
    
    /* Add the cell to the proxy */
    proxy_addcell_in(&proxies[proxy_id], ci, proxy_type);
    proxy_addcell_out(&proxies[proxy_id], cj, proxy_type);
    
    /* Store info about where to send the cell */
    cj->mpi.sendto |= (1ULL << proxy_id);
  }
}
#endif

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
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;

  /* Set up cell offsets. */
  const int bkg_offset = s->zoom_props->bkg_cell_offset ;
  const int buff_offset = s->zoom_props->buffer_cell_offset;

  /* Set up some width and distance variables. */
  double r_diag2, r_diag, r_max_zoom, r_max_buff, r_max_bkg;

  /* Declare looping variables. */
  int delta_cells, delta_m, delta_p, *cdim;

  /* Calculate r_max for each level. */
  
  /* Distance between centre of the cell and corners */
  r_diag2 = ((cells[0].width[0] * cells[0].width[0]) +
             (cells[0].width[1] * cells[0].width[1]) +
             (cells[0].width[2] * cells[0].width[2]));
  r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  r_max_zoom = 2 * r_diag;
  
  /* Distance between centre of the cell and corners */
  r_diag2 = ((cells[bkg_offset].width[0] * cells[bkg_offset].width[0]) +
             (cells[bkg_offset].width[1] * cells[bkg_offset].width[1]) +
             (cells[bkg_offset].width[2] * cells[bkg_offset].width[2]));
  r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  r_max_bkg = 2 * r_diag;

  /* Do we have buffer cells? */
  if (s->zoom_props->with_buffer_cells) {
    
    /* Distance between centre of the cell and corners */
    r_diag2 = ((cells[buff_offset].width[0] * cells[buff_offset].width[0]) +
               (cells[buff_offset].width[1] * cells[buff_offset].width[1]) +
               (cells[buff_offset].width[2] * cells[buff_offset].width[2]));
    r_diag = 0.5 * sqrt(r_diag2);
    
    /* Maximal distance from shifted CoM to any corner */
    r_max_buff = 2 * r_diag;
  } else {
    r_max_buff = 0;
  }

  /* ======================= Start with zoom cells ======================= */

  /* Get cdim. */
  cdim = s->zoom_props->cdim;

  /* Compute how many cells away we need to walk */
  delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max_zoom / theta_crit;
    delta_cells = (int)(distance / cells[0].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  delta_m = delta_cells;
  delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d top-level zoom cells away (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Get the cell. */
        struct cell *ci = &cells[cid];
        
        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Early abort  */
              if (cid >= cjd) continue;
              
              /* Get the cell. */
              struct cell *cj = &cells[cjd];

              /* Avoid completely local and foreign pairs */
              if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              /* What type of proxy do we need?
               * (proxy_cell_type_none if no proxy needed). */
              int proxy_type  = find_proxy_type(ci, cj, e, i, j, k,
                                                iii, jjj, kkk,
                                                r_max_zoom, dim, periodic);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;
              
              /* Make the proxies. */
              add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
            }
          }
        }
      }
    }
  }

  /* ======================= Now buffer cells ======================= */

  if (s->zoom_props->with_buffer_cells) {

    /* Get cdim. */
    cdim = s->zoom_props->buffer_cdim;

    /* Compute how many cells away we need to walk */
    delta_cells = 1; /*hydro case */

    /* Gravity needs to take the opening angle into account */
    if (with_gravity) {
      const double distance = 2. * r_max_buff / theta_crit;
      delta_cells = (int)(distance / cells[buff_offset].dmin) + 1;
    }

    /* Turn this into upper and lower bounds for loops */
    delta_m = delta_cells;
    delta_p = delta_cells;

    /* Special case where every cell is in range of every other one */
    if (delta_cells > cdim[0]) {
      delta_m = cdim[0];
      delta_p = cdim[0];
    }

    /* Let's be verbose about this choice */
    if (e->verbose)
      message(
          "Looking for proxies up to %d top-level buffer cells away "
          "(delta_m=%d delta_p=%d)", delta_cells, delta_m, delta_p);

    /* Loop over each cell in the space. */
    for (int i = 0; i < cdim[0]; i++) {
      for (int j = 0; j < cdim[1]; j++) {
        for (int k = 0; k < cdim[2]; k++) {

          /* Get the cell ID. */
          const int cid = cell_getid(cdim, i, j, k) + buff_offset;

          /* Get the cell. */
          struct cell *ci = &cells[cid];

          /* Loop over all its neighbours in range. */
          for (int ii = -delta_m; ii <= delta_p; ii++) {
            int iii = i + ii;
            if (iii < 0 || iii >= cdim[0]) continue;
            for (int jj = -delta_m; jj <= delta_p; jj++) {
              int jjj = j + jj;
              if (jjj < 0 || jjj >= cdim[1]) continue;
              for (int kk = -delta_m; kk <= delta_p; kk++) {
                int kkk = k + kk;
                if (kkk < 0 || kkk >= cdim[2]) continue;

                /* Get the cell ID. */
                const int cjd = cell_getid(cdim, iii, jjj, kkk) + buff_offset;

                /* Early abort  */
                if (cid >= cjd) continue;

                /* Get the cell. */
                struct cell *cj = &cells[cjd];

                /* Skip void cells, we handle these below. */
                if (ci->subtype == void_cell || cj->subtype == void_cell)
                  continue;

                /* Avoid completely local and foreign pairs */
                if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                    (ci->nodeID != nodeID && cj->nodeID != nodeID))
                  continue;

                /* What type of proxy do we need?
                 * (proxy_cell_type_none if no proxy needed). */
                int proxy_type  = find_proxy_type(ci, cj, e, i, j, k, iii, jjj,
                                                  kkk, r_max_buff, dim,
                                                  periodic);

                /* Abort if not in range at all */
                if (proxy_type == proxy_cell_type_none) continue;

                /* Make the proxies. */
                add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
              }
            }
          }
          
          /* If this is a void cell we need to loop over zoom cells. */
          if (ci->subtype == void_cell) {

            int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

            /* Loop over zoom cells. */
            for (int zoom_cid = 0; zoom_cid < nr_zoom_cells; zoom_cid++) {

              /* Get the cell. */
              struct cell *zoom_ci = &cells[zoom_cid];

              /* Loop over all its neighbours in range. */
              for (int ii = -delta_m; ii <= delta_p; ii++) {
                int iii = i + ii;
                if (iii < 0 || iii >= cdim[0]) continue;
                for (int jj = -delta_m; jj <= delta_p; jj++) {
                  int jjj = j + jj;
                  if (jjj < 0 || jjj >= cdim[1]) continue;
                  for (int kk = -delta_m; kk <= delta_p; kk++) {
                    int kkk = k + kk;
                    if (kkk < 0 || kkk >= cdim[2]) continue;
                    
                    /* Get the cell ID. */
                    const int cjd =
                      cell_getid(cdim, iii, jjj, kkk) + buff_offset;
                  
                    /* Get the cell. */
                    struct cell *cj = &cells[cjd];

                    /* Skip itself. */
                    if (cj->subtype == void_cell)
                      continue;

                    /* Avoid completely local and foreign pairs */
                    if ((zoom_ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                        (zoom_ci->nodeID != nodeID && cj->nodeID != nodeID))
                      continue;

                    /* What type of proxy do we need?
                     * (proxy_cell_type_none if no proxy needed). */
                    int proxy_type  = find_proxy_type(zoom_ci, cj, e, i, j, k,
                                                      iii, jjj, kkk,
                                                      r_max_buff, dim,
                                                      periodic);
                    
                    /* Abort if not in range at all */
                    if (proxy_type == proxy_cell_type_none) continue;
                    
                    /* Make the proxies. */
                    add_proxy(zoom_ci, cj, e, proxies, nodeID, proxy_type);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /* ======================= Now background cells ======================= */

  /* Get cdim. */
  cdim = s->cdim;

  /* Compute how many cells away we need to walk */
  delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max_bkg / theta_crit;
    delta_cells = (int)(distance / cells[bkg_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  delta_m = delta_cells;
  delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
    }
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d top-level background cells away "
        "(delta_m=%d delta_p=%d)", delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k) + bkg_offset;

        /* Get the cell. */
        struct cell *ci = &cells[cid];
        
        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_offset;

              /* Early abort  */
              if (cid >= cjd) continue;
              
              /* Get the cell. */
              struct cell *cj = &cells[cjd];

              /* Skip void cells, we handle these below, and empty cells. */
              if (ci->subtype == void_cell || ci->subtype == empty ||
                  cj->subtype == void_cell || cj->subtype == empty)
                continue;

              /* Avoid completely local and foreign pairs */
              if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              /* What type of proxy do we need?
               * (proxy_cell_type_none if no proxy needed). */
              int proxy_type  = find_proxy_type(ci, cj, e, i, j, k,
                                                iii, jjj, kkk,
                                                r_max_bkg, dim, periodic);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;
              
              /* Make the proxies. */
              add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
            }
          }
        }
        
        /* If this is a void cell we need to loop over zoom cells. */
        if (ci->subtype == void_cell) {

          int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

          /* Loop over zoom cells. */
          for (int zoom_cid = 0; zoom_cid < nr_zoom_cells; zoom_cid++) {

            /* Get the cell. */
            struct cell *zoom_ci = &cells[zoom_cid];

            /* Loop over all the background neighbours in range. */
            for (int ii = -delta_m; ii <= delta_p; ii++) {
              int iii = i + ii;
              if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
              iii = (iii + cdim[0]) % cdim[0];
              for (int jj = -delta_m; jj <= delta_p; jj++) {
                int jjj = j + jj;
                if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
                jjj = (jjj + cdim[1]) % cdim[1];
                for (int kk = -delta_m; kk <= delta_p; kk++) {
                  int kkk = k + kk;
                  if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
                  kkk = (kkk + cdim[2]) % cdim[2];
                  
                  /* Get the cell ID. */
                  const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_offset;
                  
                  /* Get the cell. */
                  struct cell *cj = &cells[cjd];

                  /* Skip itself. */
                  if (cj->subtype == void_cell)
                    continue;

                  /* Avoid completely local and foreign pairs */
                  if ((zoom_ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                      (zoom_ci->nodeID != nodeID && cj->nodeID != nodeID))
                    continue;

                  /* What type of proxy do we need?
                   * (proxy_cell_type_none if no proxy needed). */
                  int proxy_type  = find_proxy_type(zoom_ci, cj, e, i, j, k,
                                                    iii, jjj, kkk,
                                                    r_max_bkg, dim, periodic);
                  
                  /* Abort if not in range at all */
                  if (proxy_type == proxy_cell_type_none) continue;
                  
                  /* Make the proxies. */
                  add_proxy(zoom_ci, cj, e, proxies, nodeID, proxy_type);
                }
              }
            }
          }
        }

        /* If this cell contains buffer cells we need to loop over them. */
        if (ci->subtype == empty) {

          /* Loop over zoom cells. */
          for (int buff_cid = buff_offset; buff_cid < s->nr_cells; buff_cid++) {

            /* Get the cell. */
            struct cell *buff_ci = &cells[buff_cid];

            /* Skip void cells. */
            if (buff_ci->subtype == void_cell)
              continue;

            /* Loop over all the background neighbours in range. */
            for (int ii = -delta_m; ii <= delta_p; ii++) {
              int iii = i + ii;
              if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
              iii = (iii + cdim[0]) % cdim[0];
              for (int jj = -delta_m; jj <= delta_p; jj++) {
                int jjj = j + jj;
                if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
                jjj = (jjj + cdim[1]) % cdim[1];
                for (int kk = -delta_m; kk <= delta_p; kk++) {
                  int kkk = k + kk;
                  if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
                  kkk = (kkk + cdim[2]) % cdim[2];
                  
                  /* Get the cell ID. */
                  const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_offset;
                  
                  /* Get the cell. */
                  struct cell *cj = &cells[cjd];

                  /* Skip itself. */
                  if (cj->subtype == empty)
                    continue;

                  /* Avoid completely local and foreign pairs */
                  if ((buff_ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                      (buff_ci->nodeID != nodeID && cj->nodeID != nodeID))
                    continue;

                  /* What type of proxy do we need?
                   * (proxy_cell_type_none if no proxy needed). */
                  int proxy_type  = find_proxy_type(buff_ci, cj, e, i, j, k,
                                                    iii, jjj, kkk,
                                                    r_max_bkg, dim, periodic);
                  
                  /* Abort if not in range at all */
                  if (proxy_type == proxy_cell_type_none) continue;
                  
                  /* Make the proxies. */
                  add_proxy(buff_ci, cj, e, proxies, nodeID, proxy_type);
                }
              }
            }
          }
        }
      }
    }
  }
    
  /* Be clear about the time */
  if (e->verbose) {
    message("Number of proxies: %d", e->nr_proxies);
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Make edge weights from the accumulated particle sizes per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of bytes in particles per cell.
 * @param edges weights for the edges of these regions. Should be 26 * counts.
 */
void sizes_to_edges_zoom(struct space *s, double *counts, double *edges,
                         int offset, int *cdim) {

  /* Get some useful constants. */
  if (cdim == NULL)
    cdim = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
            s->zoom_props->cdim[2]};
  int iedge = 0;

  /* Find adjacency arrays for zoom cells. */
  edge_loop(cdim, offset, s, /*adjncy*/ NULL, /*xadj*/ NULL, counts, edges,
            &iedge);
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Fill the adjncy array defining the graph of cells in a space.
 *
 * See the ParMETIS and METIS manuals if you want to understand this
 * format. The cell graph consists of all nodes as vertices with edges as the
 * connections to all neighbours, so we have 26 per vertex for periodic
 * boundary, fewer than 26 on the space edges when non-periodic. Note you will
 * also need an xadj array, for METIS that would be:
 *
 *   xadj[0] = 0;
 *   for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
 *
 * but each rank needs a different xadj when using ParMETIS (each segment
 * should be rezeroed).
 *
 * @param s the space of cells.
 * @param periodic whether to assume a periodic space (fixed 26 edges).
 * @param weights_e the edge weights for the cells, if used. On input
 *                  assumed to be ordered with a fixed 26 edges per cell, so
 *                  will need reordering for non-periodic spaces.
 * @param adjncy the adjncy array to fill, must be of size 26 * the number of
 *               cells in the space.
 * @param nadjcny number of adjncy elements used, can be less if not periodic.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 * @param nxadj the number of xadj element used.
 * @param nverts the number of vertices.
 * @param offset the offset into the cell grid.
 * @param cdim the cdim of the current grid (only used when doing grids
 *                                           separately).
 */
void graph_init_zoom(struct space *s, int periodic, idx_t *weights_e,
                     idx_t *adjncy, int *nadjcny, idx_t *xadj,
                     int *nxadj, int nverts, int offset, int *cdim) {

  /* Loop over all cells in the space. */
  *nadjcny = 0;

  int iedge = 0;

  /* Find adjacency arrays for zoom cells. */
  edge_loop(cdim, offset, s, adjncy, xadj, /*counts*/ NULL, /*edges*/ NULL,
            &iedge);

  /* Set the number of adjacncy entries. */
  *nadjcny = iedge;

  /* If given set METIS xadj. */
  if (xadj != NULL) {
    xadj[nverts] = iedge;
    *nxadj = nverts;
  }

#ifdef SWIFT_DEBUG_CHECKS

  /* How many edges and vertices do we have? */
  int nverts = s->zoom_props->nr_zoom_cells + s->zoom_props->nwedges;
  int nedges = s->zoom_props->nr_edges;
  
  /* Check our adjcncy array. */
  for (int i = 0; i < nedges; i++) {
    if (adjncy[i] < 0 || adjncy[i] >= nverts)
      error("Vertex found with an incompatible adjncy (edge=%d, "
            "adjncy[edge]=%ld)", i, adjncy[i]);
  }
  /* Check our xadj array. */
  int max_edges = 0;
  int max_edge = -1;
  for (int i = 1; i < nverts; i++) {
    if ((xadj[i] - xadj[i - 1]) > max_edges) {
      max_edges = xadj[i] - xadj[i - 1];
      max_edge = i;
    }
    if (xadj[i] < 0 || xadj[i] > nedges)
      error("Vertex found with an incompatible xadj (vertex=%d, "
            "xadj[vertex]=%ld)", i, xadj[i]);
  }
  message("Max number of edges for a vertex is %d (vertex=%d)", max_edges,
          max_edge);
#endif
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Apply METIS cell-list partitioning to a cell structure.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
int get_wedge_index(struct space *s, struct cell *c) {

  if (s->zoom_props->use_bkg_wedges) {

  /* The number of slices in theta. */
  int theta_nslices = s->zoom_props->theta_nslices;
  int phi_nslices = s->zoom_props->phi_nslices;

  /* Calculate the size of a slice in theta and phi. */
  double theta_width = s->zoom_props->theta_width;
  double phi_width = s->zoom_props->phi_width;
  
  /* Center cell coordinates. */
  double dx = c->loc[0] - (s->dim[0] / 2) + c->width[0] / 2;
  double dy = c->loc[1] - (s->dim[1] / 2) + c->width[1] / 2;
  double dz = c->loc[2] - (s->dim[2] / 2) + c->width[2] / 2;

  /* Handle the central cell, just put it in wedge 0, there won't
   * be particles here anyway. */
  if (dx < (c->width[0] / 2) &&
      dy < (c->width[1] / 2) &&
      dz < (c->width[2] / 2)) {
    return 0;
  }

  /* Calculate the spherical version of these coordinates. */
  double r = sqrt(dx * dx + dy * dy + dz * dz);
  double theta = atan2(dy, dx) + M_PI;
  double phi = acos(dz / r);
    
  /* Find this wedge index. */
  int phi_ind = ((int)floor(phi / phi_width) + phi_nslices) % phi_nslices;
  int theta_ind =
    ((int)floor(theta / theta_width) + theta_nslices) % theta_nslices;
  return theta_ind * phi_nslices + phi_ind;
  } else {
    return cell_getid_pos(s,
                          c->loc[0] + c->width[0] / 2,
                          c->loc[1] + c->width[1] / 2,
                          c->loc[2] + c->width[2] / 2)
      - s->zoom_props->nr_zoom_cells;
  }
  
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Apply METIS cell-list partitioning to a cell structure.
 *
 * This version of the function handles wedges in the background.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
void split_metis_wedges(struct space *s, int nregions, int *celllist) {

  /* Get the cells array. */
  struct cell *cells = s->cells_top;
  
  /* Get how many cells we are dealing with. */
  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Setup cell counters to track the type of cells on a rank. */
  int zoom_cell_counts = 0;
  int bkg_cell_counts = 0;
  int buff_cell_counts = 0;

  /* First do the zoom cells. */
  for (int cid = 0; cid < nr_zoom_cells; cid++) {
    cells[cid].nodeID = celllist[cid];

    if (s->cells_top[cid].nodeID == s->e->nodeID)
      zoom_cell_counts++;
    
  }

  /* Now we need to loop over all the background cells and assign them based
   * on the predominant rank of zoom cells in their slice. */
    
  /* Define variables for selection */
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

  /* Loop over natural cells. Decomp these into radial slices. */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;

        /* Get the cell */
        struct cell *c = &cells[cid];

        /* Get the wedge index of this cell. */
        int wedge_ind = get_wedge_index(s, c);
        
        /* Store the rank. */
        s->cells_top[cid].nodeID = celllist[nr_zoom_cells + wedge_ind];

        if (s->cells_top[cid].nodeID == s->e->nodeID)
          bkg_cell_counts++;
        
      }
    }
  }

  /* Loop over buffer cells  */
  for (int i = 0; i < s->zoom_props->buffer_cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->buffer_cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->buffer_cdim[2]; k++) {

        /* Get cell ID. */
        const int cid =
          cell_getid(s->zoom_props->buffer_cdim, i, j, k) + buffer_cell_offset;

        /* Get the cell */
        struct cell *c = &cells[cid];

        /* Get the wedge index of this cell. */
        int wedge_ind = get_wedge_index(s, c);
        
        /* Store the rank. */
        s->cells_top[cid].nodeID = celllist[nr_zoom_cells + wedge_ind];

        if (s->cells_top[cid].nodeID == s->e->nodeID)
          buff_cell_counts++;
      }
    }
  }

  /* Report How many cells we have on a rank of each type. */
  if (s->e->verbose) {
    if (s->zoom_props->with_buffer_cells)
      message("Rank %d has %d zoom cells, %d buffer cells and %d bkg cells.",
              s->e->nodeID, zoom_cell_counts, buff_cell_counts,
              bkg_cell_counts);
    else
      message("Rank %d has %d zoom cells and %d bkg cells.", s->e->nodeID,
              zoom_cell_counts, bkg_cell_counts);
  }

  /* To check or visualise the partition dump all the cells. */
  /*if (engine_rank == 0) dumpCellRanks("metis_partition", s->cells_top,
                                      s->nr_cells);*/
}

/**
 * @brief Apply METIS cell-list partitioning to a cell structure.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
void split_metis_zoom(struct space *s, int nregions, int *celllist, int ncells,
                      int offset) {

  /* Are we doing a special decomp? */
  if (s->zoom_props->use_bkg_wedges) {
    split_metis_wedges(s, nregions, celllist);
    return;
  }

  /* Get the cells array. */
  struct cell *cells = s->cells_top;

  /* First do the zoom cells. */
  for (int cid = offset; cid < offset + ncells; cid++)
    cells[cid].nodeID = celllist[cid - offset];

  /* To check or visualise the partition dump all the cells. */
  /*if (engine_rank == 0) dumpCellRanks("metis_partition", s->cells_top,
                                      s->nr_cells);*/
}
#endif
