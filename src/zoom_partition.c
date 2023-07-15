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
#include <config.h>

/* Standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/* Include int min and max values. Define these limits in C++ as well. */
#define __STDC_LIMIT_MACROS
#include <stdint.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "gravity_properties.h"
#include "partition.h"
#include "restart.h"
#include "space.h"
#include "threadpool.h"
#include "tools.h"
#include "zoom_region.h"

/*
 * Repartition fixed costs per type/subtype. These are determined from the
 * statistics output produced when running with task debugging enabled.
 */
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
static double repartition_costs[task_type_count][task_subtype_count];
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

/*  Vectorisation support */
/*  ===================== */

#if defined(WITH_MPI)
/**
 *  @brief Pick a number of cell positions from a vectorised list.
 *
 *  Vectorise the cell space and pick positions in it for the number of
 *  expected regions using a single step. Vectorisation is guaranteed
 *  to work, providing there are more cells than regions.
 *
 *  @param s the space.
 *  @param nregions the number of regions
 *  @param samplecells the list of sample cell positions, size of 3*nregions
 */
static void pick_vector(struct space *s, int *cdim, int nregions,
                        int *samplecells) {

  /* Get length of space and divide up. */
  int length = cdim[0] * cdim[1] * cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }

  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
        if (n == 0 && l < nregions) {
          samplecells[m++] = i;
          samplecells[m++] = j;
          samplecells[m++] = k;
          l++;
        }
        n++;
        if (n == step) n = 0;
      }
    }
  }
}
#endif

#if defined(WITH_MPI)
/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometrically
 * closest and apply the partition to the space.
 */
static void split_vector(struct space *s, int *cdim, int nregions,
                         int *samplecells, int offset) {

  int n = 0;
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells_top[n++ + offset].nodeID = select;
      }
    }
  }
}
#endif

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
          
        /* Get the vertex index of this cell. */
        const size_t vid = cell_getid(cdim, i, j, k);

        /* Get the cell index. */
        const size_t cid = vid + offset;

        /* Get the cell. */
        ci = &s->cells_top[cid];

#ifdef SWIFT_DEBUG_CHECKS
        if (xadj != NULL) {
             
          /* Ensure the previous cell has found enough edges. */
          if ((vid > 0) &&
              ((*iedge - xadj[vid - 1]) != s->cells_top[cid - 1].nr_vertex_edges))
            error("Found too few edges (cid=%ld nedges=%ld, c->nr_vertex_edges=%d)",
                  cid, *iedge - xadj[vid - 1],
                  s->cells_top[cid - 1].nr_vertex_edges);
          
        }
#endif

        /* If given set METIS xadj. */
        if (xadj != NULL) {
          xadj[vid] = *iedge;
          
          /* Set edges start pointer for this cell. */
          ci->edges_start = *iedge;
        }

        /* Are we wrapping? */
        int periodic = 0;

        /* Loop over a shell of cells with the same type. */
        for (int ii = i - 1; ii <= i + 1; ii++) {
          if (!periodic && (ii < 0 || ii >= cdim[0])) continue;
          for (int jj = j - 1; jj <= j + 1; jj++) {
            if (!periodic && (jj < 0 || jj >= cdim[1])) continue;
            for (int kk = k - 1; kk <= k + 1; kk++) {
              if (!periodic && (kk < 0 || kk >= cdim[2])) continue;

              /* Apply periodic BC (not harmful if not using periodic BC) */
              const int iii = (ii + cdim[0]) % cdim[0];
              const int jjj = (jj + cdim[1]) % cdim[1];
              const int kkk = (kk + cdim[2]) % cdim[2];
                
              /* Get the vertex index of this cell. */
              const size_t vjd = cell_getid(cdim, iii, jjj, kkk);
              
              /* Get the cell index. */
              const size_t cjd = vjd + offset;
              
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
  
  /* Otherwise, loop over the cell grid we've been handed */
  else {
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
    cdim = s->zoom_props->cdim;
  
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
  
  /* Check our adjcncy array. */
  for (int i = 0; i < iedge; i++) {
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
    if (xadj[i] < 0 || xadj[i] > iedge)
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

  /* Do the current set of cells. */
  for (int cid = offset; cid < offset + ncells; cid++)
    cells[cid].nodeID = celllist[cid - offset];

  /* To check or visualise the partition dump all the cells. */
  /*if (engine_rank == 0) dumpCellRanks("metis_partition", s->cells_top,
                                      s->nr_cells);*/
}
#endif

/* ========================================================================= */

#ifdef WITH_MPI
/**
 * @brief Partition the into radial slices.
 *
 * This simply slices the box into wedges along the x-y plane.
 */
static void split_radial_wedges(struct space *s, int nregions,
                                double *weights_v, int nslices,
                                int nwedges) {

  double r, theta, phi;
    
  /* Define variables for selection */
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

  /* Calculate the size of a radial slice. */
  float slice_width = 2 * M_PI / nslices;

  /* Set up an array to store slice weights. */
  double tot_weight = 0;
  double *slice_weights;
  if ((slice_weights = (double *)malloc(sizeof(double) * nwedges)) == NULL)
    error("Failed to allocate slice_weights buffer.");
  bzero(slice_weights, sizeof(double) * nwedges);

  /* Get the weight of each slice*/

  /* Loop over zoom  */
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->zoom_props->cdim, i, j, k);

        /* Center cell coordinates. */
        int ii = i - (s->zoom_props->cdim[0] / 2);
        int jj = j - (s->zoom_props->cdim[1] / 2);
        int kk = k - (s->zoom_props->cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        slice_weights[wedge_ind] += weights_v[cid];
        tot_weight += weights_v[cid];
      }
    }
  }

  /* Loop over natural cells. Decomp these into radial slices. */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;

        /* Center cell coordinates. */
        int ii = i - (s->cdim[0] / 2);
        int jj = j - (s->cdim[1] / 2);
        int kk = k - (s->cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        slice_weights[wedge_ind] += weights_v[cid];
        tot_weight += weights_v[cid];
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

        /* Center cell coordinates. */
        int ii = i - (s->zoom_props->buffer_cdim[0] / 2);
        int jj = j - (s->zoom_props->buffer_cdim[1] / 2);
        int kk = k - (s->zoom_props->buffer_cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        slice_weights[wedge_ind] += weights_v[cid];
        tot_weight += weights_v[cid];
      }
    }
  }

  /* What would a perfectly distributed weight look like? */
  double split_weight = tot_weight / nregions;

  /* Set up an array dictating where each slice ends up. */
  int *slicelist;
  double *region_weights;
  if ((slicelist = (int *)malloc(sizeof(int) * nwedges)) == NULL)
    error("Failed to allocate slicelist");
  if ((region_weights = (double *)malloc(sizeof(double) * nregions)) == NULL)
    error("Failed to allocate region_weights buffer.");
  bzero(region_weights, sizeof(double) * nregions);

  /* Lets distribute these slices. */
  int select = 0;
  for (int islice = 0; islice < nwedges; islice++) {

    /* Assign this slice and include its weight. */
    slicelist[islice] = select;
    region_weights[select] += slice_weights[islice];

    /* Have we filled this region/rank? */
    if (region_weights[select] > split_weight && select < nregions - 1)
      select++;
  }

  /* Now lets tell each cell where it is. */
  
  /* Loop over zoom  */
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->zoom_props->cdim, i, j, k);

        /* Center cell coordinates. */
        int ii = i - (s->zoom_props->cdim[0] / 2);
        int jj = j - (s->zoom_props->cdim[1] / 2);
        int kk = k - (s->zoom_props->cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        s->cells_top[cid].nodeID = slicelist[wedge_ind];
      }
    }
  }

  /* Loop over natural cells. Decomp these into radial slices. */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;

        /* Center cell coordinates. */
        int ii = i - (s->cdim[0] / 2);
        int jj = j - (s->cdim[1] / 2);
        int kk = k - (s->cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        s->cells_top[cid].nodeID = slicelist[wedge_ind];
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

        /* Center cell coordinates. */
        int ii = i - (s->zoom_props->buffer_cdim[0] / 2);
        int jj = j - (s->zoom_props->buffer_cdim[1] / 2);
        int kk = k - (s->zoom_props->buffer_cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        int phi_ind = phi / slice_width / 2;
        int theta_ind = theta / slice_width;
        int wedge_ind = theta_ind * nslices + phi_ind;
        s->cells_top[cid].nodeID = slicelist[wedge_ind];
      }
    }
  }

  free(slice_weights);
  free(slicelist);
  free(region_weights);

}
#endif


/* METIS/ParMETIS support (optional)
 * =================================
 *
 * METIS/ParMETIS partitions using a multi-level k-way scheme. We support
 * using this in a unweighted scheme, which works well and seems to be
 * guaranteed, and a weighted by the number of particles scheme.
 *
 * Repartitioning is based on ParMETIS and uses weights determined from the
 * estimated costs that a cells tasks will take or the relative time bins of
 * the cells next updates.
 */
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/* qsort support. */
struct indexval {
  int index;
  int count;
  int old_val;
  int new_val;
};
static int indexvalcmp(const void *p1, const void *p2) {
  const struct indexval *iv1 = (const struct indexval *)p1;
  const struct indexval *iv2 = (const struct indexval *)p2;
  return iv2->count - iv1->count;
}

/**
 * @brief Check if there is a permutation of the region indices of our cells
 *        that will reduce the amount of particle movement and return it.
 *
 * @param newlist the new list of regions for our cells.
 * @param oldlist the old list of regions for our cells.
 * @param nregions the number of regions.
 * @param ncells the number of cells.
 * @param permlist the permutation of the newlist.
 */
void permute_regions_zoom(int *newlist, int *oldlist, int nregions, int ncells,
                     int *permlist) {

  /* We want a solution in which the current region assignments of the cells
   * are preserved when possible, to avoid unneccesary particle movement.  So
   * create a 2d-array of counts of cells that are common to all pairs of old
   * and new lists. Each element of the array has a count of cells and an
   * unique index so we can sort into decreasing counts.
   */
  int indmax = nregions * nregions;
  struct indexval *ivs = NULL;
  if ((ivs = (struct indexval *)malloc(sizeof(struct indexval) * indmax)) ==
      NULL)
    error("Failed to allocate ivs structs");
  bzero(ivs, sizeof(struct indexval) * indmax);

  for (int k = 0; k < ncells; k++) {
    int index = newlist[k] + nregions * oldlist[k];
    ivs[index].count++;
    ivs[index].index = index;
    ivs[index].old_val = oldlist[k];
    ivs[index].new_val = newlist[k];
  }
  qsort(ivs, indmax, sizeof(struct indexval), indexvalcmp);

  /* Go through the ivs using the largest counts first, these are the
   * regions with the most cells in common, old partition to new. If not
   * returning the permutation, avoid the associated work. */
  int *oldmap = NULL;
  int *newmap = NULL;
  oldmap = permlist; /* Reuse this */
  if ((newmap = (int *)malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate newmap array");

  for (int k = 0; k < nregions; k++) {
    oldmap[k] = -1;
    newmap[k] = -1;
  }

  for (int k = 0; k < indmax; k++) {

    /* Stop when all regions with common cells have been considered. */
    if (ivs[k].count == 0) break;

    /* Store old and new IDs, if not already used. */
    if (newmap[ivs[k].new_val] == -1 && oldmap[ivs[k].old_val] == -1) {
      newmap[ivs[k].new_val] = ivs[k].old_val;
      oldmap[ivs[k].old_val] = ivs[k].new_val;
    }
  }

  /* Handle any regions that did not get selected by picking an unused rank
   * from oldmap and assigning to newmap. */
  int spare = 0;
  for (int k = 0; k < nregions; k++) {
    if (newmap[k] == -1) {
      for (int j = spare; j < nregions; j++) {
        if (oldmap[j] == -1) {
          newmap[k] = j;
          oldmap[j] = j;
          spare = j;
          break;
        }
      }
    }
  }

  /* Permute the newlist into this order. */
  for (int k = 0; k < ncells; k++) {
    permlist[k] = newmap[newlist[k]];
  }
  free(newmap);
  free(ivs);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
/**
 * @brief Partition the given space into a number of connected regions using
 *        ParMETIS.
 *
 * Split the space using PARMETIS to derive a partitions using the
 * given edge and vertex weights. If no weights are given then an
 * unweighted partition is performed. If refine is set then an existing
 * partition is assumed to be present from the last call to this routine
 * in the celllist argument, that will get a refined partition, not a new
 * one.
 *
 * Assumes MPI is up and running and the number of ranks is the same as the
 * number of regions.
 *
 * @param nodeID our nodeID.
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param ncells the number of vertices in the graph.
 * @param nedges the total number of edges in the graph.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights. Need to be in the range of idx_t.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array. Need to be in the range of
 *        idx_t.
 * @param refine whether to refine an existing partition, or create a new one.
 * @param adaptive whether to use an adaptive reparitition of an existing
 *        partition or simple refinement. Adaptive repartition is controlled
 *        by the itr parameter.
 * @param itr the ratio of inter-process communication time to data
 *            redistribution time. Used to weight repartitioning edge cuts
 *            when refine and adaptive are true.
 * @param celllist on exit this contains the ids of the selected regions,
 *        size of number of cells. If refine is 1, then this should contain
 *        the old partition on entry.
 * @param cell_offset the offset into the cell grid.
 * @param cdim the cdim of the current grid (only used when doing grids
 *                                           separately).
 */
static void pick_parmetis(int nodeID, struct space *s, int nregions,
                          int ncells, int nedges, double *vertexw,
                          double *edgew, int refine, int adaptive,
                          float itr, int *celllist, int cell_offset,
                          int *cdim) {

  int res;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  
  /* Nothing much to do if only using a single MPI rank. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* We all get one of these with the same content. It defines the ranges of
   * vertices that are found on each rank. This contiguity constraint seems to
   * stop efficient local processing, since our cell distributions do not
   * meet this requirement. That means the graph and related information needs
   * to be all brought to one node and redistributed for processing in
   * approproiate batches. */
  idx_t *vtxdist;
  if ((vtxdist = (idx_t *)malloc(sizeof(idx_t) * (nregions + 1))) == NULL)
    error("Failed to allocate vtxdist buffer.");

  if (nodeID == 0) {

    /* Construct vtxdist and send it to all ranks. Each rank gets an equal
     * number of vertices. */
    vtxdist[0] = 0;
    int k = ncells;
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

  /* Number of cells on this node and space for the expected arrays. */
  int nverts = vtxdist[nodeID + 1] - vtxdist[nodeID];

  /* We need to count how many edges are on this rank in the zoom case. */
  int nr_my_edges = 0;
  for (int cid = vtxdist[nodeID]; cid < vtxdist[nodeID + 1]; cid++) {
    if (cid < s->zoom_props->nr_zoom_cells || s->zoom_props->separate_decomps)
      nr_my_edges += s->cells_top[cid].nr_vertex_edges;
    else
      nr_my_edges +=
        s->zoom_props->nr_wedge_edges[cid - s->zoom_props->nr_zoom_cells];
  }

  idx_t *xadj = NULL;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
    error("Failed to allocate xadj buffer.");

  idx_t *adjncy = NULL;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * nr_my_edges)) == NULL)
    error("Failed to allocate adjncy array (nr_local_edges=%d).", nr_my_edges);

  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate vertex weights array");

  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)malloc(nr_my_edges * sizeof(idx_t))) == NULL)
      error("Failed to allocate edge weights array");

  idx_t *regionid = NULL;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
    error("Failed to allocate regionid array");

  /* Prepare MPI requests for the asynchronous communications */
  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 5 * nregions)) ==
      NULL)
    error("Failed to allocate MPI request list.");
  for (int k = 0; k < 5 * nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  MPI_Status *stats;
  if ((stats = (MPI_Status *)malloc(sizeof(MPI_Status) * 5 * nregions)) == NULL)
    error("Failed to allocate MPI status list.");

  /* Only use one rank to organize everything. */
  if (nodeID == 0) {

    /* Space for largest lists. */
    idx_t *full_xadj = NULL;
    if ((full_xadj =
             (idx_t *)malloc(sizeof(idx_t) * (ncells + nregions + 1))) == NULL)
      error("Failed to allocate full xadj buffer.");
    idx_t *std_xadj = NULL;
    if ((std_xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
      error("Failed to allocate std xadj buffer.");
    idx_t *full_adjncy = NULL;
    if ((full_adjncy = (idx_t *)malloc(sizeof(idx_t) * nedges)) == NULL)
      error("Failed to allocate full adjncy array.");
    idx_t *full_weights_v = NULL;
    if (weights_v != NULL)
      if ((full_weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate full vertex weights array");
    idx_t *full_weights_e = NULL;
    if (weights_e != NULL)
      if ((full_weights_e = (idx_t *)malloc(nedges * sizeof(idx_t))) == NULL)
        error("Failed to allocate full edge weights array");

    idx_t *full_regionid = NULL;
    if (refine) {
      if ((full_regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate full regionid array");
    }

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          full_weights_v[k] = vertexw[k];
        } else {
          full_weights_v[k] = 0;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < ncells; k++) {
        if ((idx_t)vertexw[k] < 0) {
          message("Input vertex weight out of range: %ld", (long)vertexw[k]);
          failed++;
        }
        if (full_weights_v[k] < 0) {
          message("Used vertex weight  out of range: %" PRIDX,
                  full_weights_v[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d vertex weights are out of range", failed);
#endif
    }

    /* Init the edges weights array. */
    if (edgew != NULL) {
      for (int k = 0; k < nedges; k++) {
        if (edgew[k] > 1) {
          full_weights_e[k] = edgew[k];
        } else {
          full_weights_e[k] = 1;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < nedges; k++) {

        if ((idx_t)edgew[k] < 0) {
          message("Input edge weight out of range: %ld", (long)edgew[k]);
          failed++;
        }
        if (full_weights_e[k] < 1) {
          message("Used edge weight out of range: %" PRIDX, full_weights_e[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d edge weights are out of range", failed);
#endif
    }

    message("Got to init");

    /* Define the cell graph. Keeping the edge weights association. */
    int nadjcny = 0;
    int nxadj = 0;
    graph_init_zoom(s, s->periodic, full_weights_e, full_adjncy, &nadjcny,
                    std_xadj, &nxadj, ncells, cell_offset, cdim);
    
    /* Dump graphs to disk files for testing. */
    /*dumpMETISGraph("parmetis_graph", ncells, 1, std_xadj, full_adjncy,
                   full_weights_v, NULL, full_weights_e); */

    /* xadj is set for each rank, different to serial version in that each
     * rank starts with 0, so we need to re-offset. */
    for (int rank = 0, i = 0, j = 0; rank < nregions; rank++) {

      /* Number of vertices for this rank. */
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      /* Each xadj section starts at 0 and terminates like a complete one. */
      int offset = std_xadj[j];
      for (int k = 0; k < nvt; k++) {
        full_xadj[i] = std_xadj[j] - offset;
        j++;
        i++;
      }
      full_xadj[i] = std_xadj[j] - offset;
      i++;
    }

    /* Send ranges to the other ranks and keep our own. */
    for (int rank = 0, j1 = 0, j2 = 0, j3 = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];
      int nedge = std_xadj[vtxdist[rank + 1]] - std_xadj[vtxdist[rank]];

      if (refine)
        for (int i = 0; i < nvt; i++) full_regionid[j3 + i] = celllist[j3 + i];

      if (rank == 0) {
        memcpy(xadj, &full_xadj[j1], sizeof(idx_t) * (nvt + 1));
        memcpy(adjncy, &full_adjncy[j2], sizeof(idx_t) * nedge);
        if (weights_e != NULL)
          memcpy(weights_e, &full_weights_e[j2], sizeof(idx_t) * nedge);
        if (weights_v != NULL)
          memcpy(weights_v, &full_weights_v[j3], sizeof(idx_t) * nvt);
        if (refine) memcpy(regionid, full_regionid, sizeof(idx_t) * nvt);

      } else {
        res = MPI_Isend(&full_xadj[j1], nvt + 1, IDX_T, rank, 0, comm,
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

      /* Note we send 26 edges, but only increment by the correct number. */
      j2 += nedge;
      j3 += nvt;
    }

    /* Wait for all sends to complete. */
    int result;
    if ((result = MPI_Waitall(5 * nregions, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5 * nregions; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &result);
        message("send request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall sending repartition data.");
    }

    /* Clean up. */
    if (weights_v != NULL) free(full_weights_v);
    if (weights_e != NULL) free(full_weights_e);
    free(full_xadj);
    free(std_xadj);
    free(full_adjncy);
    if (refine) free(full_regionid);

  } else {

    /* Receive stuff from rank 0. */
    res = MPI_Irecv(xadj, nverts + 1, IDX_T, 0, 0, comm, &reqs[0]);
    if (res == MPI_SUCCESS)
      res = MPI_Irecv(adjncy, nr_my_edges, IDX_T, 0, 1, comm, &reqs[1]);
    if (res == MPI_SUCCESS && weights_e != NULL)
      res = MPI_Irecv(weights_e, nr_my_edges, IDX_T, 0, 2, comm, &reqs[2]);
    if (res == MPI_SUCCESS && weights_v != NULL)
      res = MPI_Irecv(weights_v, nverts, IDX_T, 0, 3, comm, &reqs[3]);
    if (refine && res == MPI_SUCCESS)
      res += MPI_Irecv((void *)regionid, nverts, IDX_T, 0, 4, comm, &reqs[4]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive graph data");

    /* Wait for all recvs to complete. */
    int result;
    if ((result = MPI_Waitall(5, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &result);
        message("recv request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall receiving repartition data.");
    }
  }

  /* Set up the tpwgts array. This is just 1/nregions. */
  real_t *tpwgts;
  if ((tpwgts = (real_t *)malloc(sizeof(real_t) * nregions)) == NULL)
    error("Failed to allocate tpwgts array");
  for (int i = 0; i < nregions; i++) tpwgts[i] = 1.0 / (real_t)nregions;

  /* Common parameters. */
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

  message("Lets enter parmetis");

  if (refine) {
    /* Refine an existing partition, uncouple as we do not have the cells
     * present on their expected ranks. */
    options[3] = PARMETIS_PSR_UNCOUPLED;

    /* Seed for randoms. */
    options[2] = clocks_random_seed();

    /* Choice is whether to use an adaptive repartition or a simple
     * refinement. */
    if (adaptive) {

      /* Balance between cuts and movement. */
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

    /* Create a new partition. Use a number of guesses as that is similar to
     * the way that serial METIS works (serial METIS usually gives the best
     * quality partitions). */
    idx_t best_edgecut = 0;
    idx_t *best_regionid = NULL;
    if ((best_regionid = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
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
        memcpy(best_regionid, regionid, sizeof(idx_t) * (nverts + 1));
      }
    }

    /* Keep the best edgecut. */
    memcpy(regionid, best_regionid, sizeof(idx_t) * (nverts + 1));
    free(best_regionid);
  }

  message("Got out of parmetis");

  /* Need to gather all the regionid arrays from the ranks. */
  for (int k = 0; k < nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  if (nodeID != 0) {

    /* Send our regions to node 0. */
    res = MPI_Isend(regionid, vtxdist[nodeID + 1] - vtxdist[nodeID], IDX_T, 0,
                    1, comm, &reqs[0]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to send new regionids");

    /* Wait for send to complete. */
    int err;
    if ((err = MPI_Wait(reqs, stats)) != MPI_SUCCESS) {
      mpi_error(err, "Failed during wait sending regionids.");
    }

  } else {

    /* Node 0 */
    idx_t *remoteids = NULL;
    if ((remoteids = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate remoteids buffer");

    int nvt = vtxdist[1] - vtxdist[0];
    memcpy(remoteids, regionid, sizeof(idx_t) * nvt);

    /* Receive from other ranks. */
    for (int rank = 1, j = nvt; rank < nregions; rank++) {
      nvt = vtxdist[rank + 1] - vtxdist[rank];
      res = MPI_Irecv((void *)&remoteids[j], nvt, IDX_T, rank, 1, comm,
                      &reqs[rank]);
      if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive new regionids");
      j += nvt;
    }

    int err;
    if ((err = MPI_Waitall(nregions, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &err);
        message("recv request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall receiving regionid data.");
    }

    /* Copy: idx_t -> int. */
    int *newcelllist = NULL;
    if ((newcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
      error("Failed to allocate new celllist");
    for (int k = 0; k < ncells; k++) newcelllist[k] = remoteids[k];
    free(remoteids);

    /* Check that the region ids are all good. */
    int bad = 0;
    for (int k = 0; k < ncells; k++) {
      if (newcelllist[k] < 0 || newcelllist[k] >= nregions) {
        message("Got bad nodeID %d for cell %i.", newcelllist[k], k);
        bad++;
      }
    }
    if (bad) error("Bad node IDs located");

    /* Now check the similarity to the old partition and permute if necessary.
     * Checks show that refinement can return a permutation of the partition,
     * we need to check that and correct as necessary. */
    int permute = 1;
    if (!refine && s->zoom_props->separate_decomps) {

      /* No old partition was given, so we need to construct the existing
       * partition from the cells, if one existed. */
      int nsum = 0;
      for (int i = 0; i < ncells; i++) {
        celllist[i] = s->cells_top[cell_offset + i].nodeID;
        nsum += celllist[i];
      }

      /* If no previous partition then all nodeIDs will be set to 0. */
      if (nsum == 0) permute = 0;
    }

    else if (!refine) {

      /* No old partition was given, so we need to construct the existing
       * partition from the cells, if one existed. Unlike the periodic case
       * we need to zero the array to account for the possibility of
       * unpopulated wedges which need to be zeroed. */
      bzero(celllist, sizeof(int) * ncells);
      int nsum = 0;
      for (int i = 0; i < s->nr_cells; i++) {

        /* If this is a zoom cell then the vertex index == cell index. */
        if (i < s->zoom_props->nr_zoom_cells) {
          celllist[i] = s->cells_top[i].nodeID;
          nsum += celllist[i];
        }

        /* Otherwise, we need to find the wedge index. */
        else {
          int iwedge =
            get_wedge_index(s, &s->cells_top[i]) + s->zoom_props->nr_zoom_cells;
          celllist[iwedge] = s->cells_top[i].nodeID;
          nsum += celllist[iwedge];
        }
      }

      /* If no previous partition then all nodeIDs will be set to 0. */
      if (nsum == 0) permute = 0;
    }

    /* Ensure the celllist is valid. */
    for (int k = 0; k < ncells; k++) {
      if (celllist[k] < 0 || celllist[k] >= nregions) {
        message("Got bad nodeID %d for cell %i.", celllist[k], k);
        bad++;
      }
    }
    if (bad) error("Bad node IDs located (refine=%d)", refine);

    if (permute) {
      
      int *permcelllist = NULL;
      if ((permcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
        error("Failed to allocate perm celllist array");
      permute_regions_zoom(newcelllist, celllist, nregions, ncells,
                           permcelllist);

      /* And keep. */
      memcpy(celllist, permcelllist, sizeof(int) * ncells);
      free(permcelllist);

    } else {
      memcpy(celllist, newcelllist, sizeof(int) * ncells);
    }
    free(newcelllist);
  }

  message("Completed permutation");

  /* And everyone gets a copy. */
  res = MPI_Bcast(celllist, ncells, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");

  /* Clean up. */
  free(reqs);
  free(stats);
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(vtxdist);
  free(tpwgts);
  free(xadj);
  free(adjncy);
  free(regionid);
}
#endif


#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the given edge and
 * vertex weights. If no weights are given then an unweighted partition is
 * performed.
 *
 * @param nodeID the rank of our node.
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param nverts the number of vertices in the graph.
 * @param nedges the total number of edges in the graph.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights. Need to be in the range of idx_t.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array. Need to be in the range of
 *        idx_t.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 * @param offset the offset into the cell grid.
 * @param cdim the cdim of the current grid (only used when doing grids
 *                                           separately).
 */
static void pick_metis(int nodeID, struct space *s, int nregions, int nverts,
                       int nedges,  double *vertexw, double *edgew,
                       int *celllist, int offset, int *cdim) {

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (nregions == 1) {
    for (int i = 0; i < nverts; i++) celllist[i] = 0;
    return;
  }

  /* Only one node needs to calculate this. */
  if (nodeID == 0) {

    /* Allocate adjacency and weights arrays . */
    idx_t *xadj;
    if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
      error("Failed to allocate xadj buffer.");
    idx_t *adjncy;
    if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * nedges)) == NULL)
      error("Failed to allocate adjncy array.");
    idx_t *weights_v = NULL;
    if (vertexw != NULL)
      if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
        error("Failed to allocate vertex weights array");
    idx_t *weights_e = NULL;
    if (edgew != NULL)
      if ((weights_e = (idx_t *)malloc(nedges * sizeof(idx_t))) == NULL)
        error("Failed to allocate edge weights array");
    idx_t *regionid;
    if ((regionid = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate regionid array");

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < nverts; k++) {
        if (vertexw[k] > 1) {
          weights_v[k] = vertexw[k];
        } else {
          weights_v[k] = 0;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < nverts; k++) {
        if ((idx_t)vertexw[k] < 0) {
          message("Input vertex weight out of range: %ld", (long)vertexw[k]);
          failed++;
        }
        if (weights_v[k] < 0) {
          message("Used vertex weight  out of range: %" PRIDX, weights_v[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d vertex weights are out of range", failed);
#endif
    }

    /* Init the edges weights array. */

    if (edgew != NULL) {
      for (int k = 0; k < nedges; k++) {
        if (edgew[k] > 1) {
          weights_e[k] = edgew[k];
        } else {
          weights_e[k] = 1;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < nedges; k++) {

        if ((idx_t)edgew[k] < 0) {
          message("Input edge weight out of range: %ld", (long)edgew[k]);
          failed++;
        }
        if (weights_e[k] < 1) {
          message("Used edge weight out of range: %" PRIDX, weights_e[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d edge weights are out of range", failed);
#endif
    }

    /* Define the cell graph. Keeping the edge weights association. */
    int nadjcny = 0;
    int nxadj = 0;
    graph_init_zoom(s, s->periodic, weights_e, adjncy, &nadjcny, xadj, &nxadj,
               nverts, offset, cdim);

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

    /* Dump graph in METIS format */
    /* dumpMETISGraph("metis_graph", idx_nverts, one, xadj, adjncy, weights_v, */
    /*                NULL, weights_e); */

    if (METIS_PartGraphKway(&idx_nverts, &one, xadj, adjncy, weights_v, NULL,
                            weights_e, &idx_nregions, NULL, NULL, options,
                            &objval, regionid) != METIS_OK)
      error("Call to METIS_PartGraphKway failed.");

    /* Check that the regionids are ok. */
    for (int k = 0; k < nverts; k++) {
      if (regionid[k] < 0 || regionid[k] >= nregions)
        error("Got bad nodeID %" PRIDX " for cell %i.", regionid[k], k);

      /* And keep. */
      celllist[k] = regionid[k];
    }

    /* Clean up. */
    if (weights_v != NULL) free(weights_v);
    if (weights_e != NULL) free(weights_e);
    free(xadj);
    free(adjncy);
    free(regionid);
  }

  /* Calculations all done, now everyone gets a copy. */
  int res = MPI_Bcast(celllist, nverts, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");
}
#endif


#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/* Helper struct for partition_gather weights. */
struct weights_mapper_data {
  double *weights_e;
  double *weights_v;
  idx_t *inds;
  int eweights;
  int nodeID;
  int timebins;
  int vweights;
  int nedges;
  int nr_cells;
  int use_ticks;
  struct cell *cells;
};

/* #ifdef SWIFT_DEBUG_CHECKS */
/* static void check_weights(struct task *tasks, int nr_tasks, */
/*                           struct weights_mapper_data *weights_data, */
/*                           double *weights_v, double *weights_e); */
/* #endif */

/**
 * @brief Threadpool mapper function to gather cell edge and vertex weights
 *        from the associated tasks.
 *
 * @param map_data part of the data to process in this mapper.
 * @param num_elements the number of data elements to process.
 * @param extra_data additional data for the mapper context.
 */
void partition_gather_weights_zoom(void *map_data, int num_elements,
                                   void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  struct weights_mapper_data *mydata = (struct weights_mapper_data *)extra_data;

  double *weights_e = mydata->weights_e;
  double *weights_v = mydata->weights_v;
  idx_t *inds = mydata->inds;
  int eweights = mydata->eweights;
  int nodeID = mydata->nodeID;
  int nedges = mydata->nedges;
  int timebins = mydata->timebins;
  int vweights = mydata->vweights;
  int use_ticks = mydata->use_ticks;

  struct cell *cells = mydata->cells;

  /* Loop over the tasks... */
  for (int i = 0; i < num_elements; i++) {
    struct task *t = &tasks[i];

    /* Skip un-interesting tasks. */
    if (t->type == task_type_send || t->type == task_type_recv ||
        t->type == task_type_csds || t->implicit || t->ci == NULL)
      continue;

    /* Get weight for this task. Either based on fixed costs or task timings. */
    double w = 0.0;
    if (use_ticks) {
      w = (double)t->toc - (double)t->tic;
    } else {
      w = repartition_costs[t->type][t->subtype];
    }
    if (w <= 0.0) continue;

    /* Get the top-level cells involved. */
    struct cell *ci, *cj;
    for (ci = t->ci; ci->parent != NULL; ci = ci->parent)
      ;
    if (t->cj != NULL)
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent)
        ;
    else
      cj = NULL;

    /* Get the cell IDs. */
    int cid = ci - cells;

    /* Different weights for different tasks. */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_ghost || t->type == task_type_extra_ghost ||
        t->type == task_type_kick1 || t->type == task_type_kick2 ||
        t->type == task_type_end_hydro_force ||
        t->type == task_type_end_grav_force || t->type == task_type_cooling ||
        t->type == task_type_star_formation || t->type == task_type_timestep ||
        t->type == task_type_init_grav || t->type == task_type_grav_down ||
        t->type == task_type_grav_long_range ||
        t->type == task_type_grav_long_range_bkg) {

      /* Particle updates add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);
    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub_self && cj == NULL &&
              ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);

    }

    /* Pair? */
    else if (t->type == task_type_pair || (t->type == task_type_sub_pair)) {

      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (vweights) atomic_add_d(&weights_v[cid], w);

      }

      /* Distinct cells. */
      else {

        /* Index of the jth cell. */
        int cjd = cj - cells;
        
        /* Local cells add weight to vertices. */
        if (vweights && ci->nodeID == nodeID) {
          atomic_add_d(&weights_v[cid], 0.5 * w);
          if (cj->nodeID == nodeID) atomic_add_d(&weights_v[cjd], 0.5 * w);
        }

        if (eweights) {

          /* Find indices of ci/cj neighbours. Note with gravity these cells may
           * not be neighbours, in that case we ignore any edge weight for that
           * pair. */
          int ik = -1;
          for (int k = ci->edges_start; k < nedges; k++) {
            if (inds[k] == cjd) {
              ik = k;
              break;
            }
          }

          /* cj */
          int jk = -1;
          for (int k = cj->edges_start; k < nedges; k++) {
            if (inds[k] == cid) {
              jk = k;
              break;
            }
          }

          if (ik != -1 && jk != -1) {

            if (timebins) {
              /* Add weights to edge for all cells based on the expected
               * interaction time (calculated as the time to the last expected
               * time) as we want to avoid having active cells on the edges, so
               * we cut for that. Note that weight is added to the local and
               * remote cells, as we want to keep both away from any cuts, this
               * can overflow int, so take care. */
              int dti = num_time_bins - get_time_bin(ci->hydro.ti_end_min);
              int dtj = num_time_bins - get_time_bin(cj->hydro.ti_end_min);
              double dt = (double)(1 << dti) + (double)(1 << dtj);
              atomic_add_d(&weights_e[ik], dt);
              atomic_add_d(&weights_e[jk], dt);

            } else {

              /* Add weights from task costs to the edge. */
              atomic_add_d(&weights_e[ik], w);
              atomic_add_d(&weights_e[jk], w);
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Repartition the zoom cells amongst the nodes using weights based on
 *        the memory use of particles in the cells.
 *
 * @param repartition the partition struct of the local engine.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 */
void repart_memory_metis_zoom(struct repartition *repartition, int nodeID,
                              int nr_nodes, struct space *s) {

  /* Total number of cells. */
  int ncells = s->zoom_props->nr_zoom_cells + s->zoom_props->nwedges;

  /* Get the particle weights in all cells. */
  double *cell_weights;
  if ((cell_weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
    error("Failed to allocate cell_weights buffer.");
  accumulate_sizes(s, s->e->verbose, cell_weights);

  /* Space for counts of particle memory use per cell. */
  double *weights = NULL;
  if ((weights = (double *)malloc(sizeof(double) * ncells)) == NULL)
    error("Failed to allocate cell weights buffer.");

  /* Check each particle and accumulate the sizes per cell. */
  for (int cid = 0; cid < s->zoom_props->nr_zoom_cells; cid++)
        weights[cid] = cell_weights[cid];

  /* Get the wedge weights. */
  for (int cid = s->zoom_props->nr_zoom_cells; cid < s->nr_cells; cid++) {

    /* Get the cell. */
    struct cell *c = &s->cells_top[cid];
      
    /* Find this wedge index. */
    int wedge_ind = get_wedge_index(s, c);

    /* Add this weight. */
    weights[s->zoom_props->nr_zoom_cells + wedge_ind] += cell_weights[cid];
  }

  /* Allocate cell list for the partition. If not already done. */
#ifdef HAVE_PARMETIS
  int refine = 1;
#endif
  if (repartition->ncelllist != ncells) {
#ifdef HAVE_PARMETIS
    refine = 0;
#endif
    free(repartition->celllist);
    repartition->ncelllist = 0;
    if ((repartition->celllist = (int *)malloc(sizeof(int) * ncells)) ==
        NULL)
      error("Failed to allocate celllist");
    repartition->ncelllist = ncells;
  }

  /* We need to rescale the sum of the weights so that the sum is
   * less than IDX_MAX, that is the range of idx_t. */
  double sum = 0.0;
  for (int k = 0; k < ncells; k++) sum += weights[k];
  if (sum > (double)IDX_MAX) {
    double scale = (double)(IDX_MAX - 1000) / sum;
    for (int k = 0; k < ncells; k++) weights[k] *= scale;
  }

  /* And repartition. */
#ifdef HAVE_PARMETIS
  if (repartition->usemetis) {
    pick_metis(nodeID, s, nr_nodes, ncells, 0, weights, NULL,
               repartition->celllist, 0, NULL);
  } else {
    pick_parmetis(nodeID, s, nr_nodes, ncells, 0, weights, NULL, refine,
                  repartition->adaptive, repartition->itr,
                  repartition->celllist, 0, NULL);
  }
#else
  pick_metis(nodeID, s, nr_nodes, ncells, 0, weights, NULL, repartition->celllist,
             0, NULL);
#endif

  /* Check that all cells have good values. All nodes have same copy, so just
   * check on one. */
  if (nodeID == 0) {
    for (int k = 0; k < ncells; k++)
      if (repartition->celllist[k] < 0 || repartition->celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", repartition->celllist[k], k);
  }

  /* Check that the zoom partition is complete and all nodes have some cells. */
  int present[nr_nodes];
  int failed = 0;
  for (int i = 0; i < nr_nodes; i++) present[i] = 0;
  for (int i = 0; i < ncells; i++) present[repartition->celllist[i]]++;
  for (int i = 0; i < nr_nodes; i++) {
    if (!present[i]) {
      failed = 1;
      if (nodeID == 0) message("Node %d is not present after repartition", i);
    }
  }

  /* If zoom partition failed continue with the current one, but make this
   * clear. */
  if (failed) {
    if (nodeID == 0)
      message(
          "WARNING: repartition has failed, continuing with the current"
          " partition, load balance will not be optimal");
    for (int k = 0; k < ncells; k++)
      repartition->celllist[k] = s->cells_top[k].nodeID;
  }

  /* And apply to our cells */
  split_metis_zoom(s, nr_nodes, repartition->celllist, ncells, 0);

  free(cell_weights);
}

/**
 * @brief Repartition the zoom cells amongst the nodes using weights of
 *        various kinds.
 *
 * @param vweights whether vertex weights will be used.
 * @param eweights whether weights will be used.
 * @param timebins use timebins as the edge weights.
 * @param repartition the partition struct of the local engine.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
static void repart_edge_metis_zoom(int vweights, int eweights, int timebins,
                                   struct repartition *repartition, int nodeID,
                                   int nr_nodes, struct space *s,
                                   struct task *tasks, int nr_tasks) {

  /* Before doing anything else we need to do the zoom cells. */
  
  /* Define the number of vertices */
  int nverts = s->zoom_props->nr_zoom_cells;
  
  /* Define the number of edges we have to handle. */
  int nedges = 0;
  for (int cid = 0; cid < nverts; cid++) {
    nedges += s->cells_top[cid].nr_vertex_edges;
  }

  /* Get the cells */
  struct cell *cells = s->cells_top;
  
  /* Allocate and fill the adjncy indexing array defining the graph of
   * cells. */
  idx_t *inds;
  if ((inds = (idx_t *)malloc(sizeof(idx_t) * nedges)) == NULL)
    error("Failed to allocate the inds array");
  int nadjcny = 0;
  int nxadj = 0;
  graph_init_zoom(s, 1 /* periodic */, NULL /* no edge weights */, inds,
                  &nadjcny,  NULL /* no xadj needed */, &nxadj, nverts, 0,
                  s->zoom_props->cdim);

  /* Allocate and init weights. */
  double *weights_v = NULL;
  double *weights_e = NULL;
  if (vweights) {
    if ((weights_v = (double *)malloc(sizeof(double) * nverts)) == NULL)
      error("Failed to allocate vertex weights arrays.");
    bzero(weights_v, sizeof(double) * nverts);
  }
  if (eweights) {
    if ((weights_e = (double *)malloc(sizeof(double) * nedges)) == NULL)
      error("Failed to allocate edge weights arrays.");
    bzero(weights_e, sizeof(double) * nedges);
  }

  /* Gather weights. */
  struct weights_mapper_data weights_data;

  weights_data.cells = cells;
  weights_data.eweights = eweights;
  weights_data.inds = inds;
  weights_data.nodeID = nodeID;
  weights_data.nedges = nedges;
  weights_data.nr_cells = nverts;
  weights_data.timebins = timebins;
  weights_data.vweights = vweights;
  weights_data.weights_e = weights_e;
  weights_data.weights_v = weights_v;
  weights_data.use_ticks = repartition->use_ticks;

  ticks tic = getticks();

  threadpool_map(&s->e->threadpool, partition_gather_weights_zoom, tasks,
                 nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                 &weights_data);
  if (s->e->verbose)
    message("weight mapper took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

/* #ifdef SWIFT_DEBUG_CHECKS */
/*   check_weights(tasks, nr_tasks, &weights_data, weights_v, weights_e); */
/* #endif */

  /* Merge the weights arrays across all nodes. */
  int res;
  if (vweights) {
    res = MPI_Allreduce(MPI_IN_PLACE, weights_v, nverts, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to allreduce vertex weights.");
  }

  if (eweights) {
    res = MPI_Allreduce(MPI_IN_PLACE, weights_e, nedges,
                        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to allreduce edge weights.");
  }

  /* Allocate cell list for the partition. If not already done. */
#ifdef HAVE_PARMETIS
  int refine = 1;
#endif
  if (repartition->ncelllist != nverts) {
#ifdef HAVE_PARMETIS
    refine = 0;
#endif
    free(repartition->celllist);
    repartition->ncelllist = 0;
    if ((repartition->celllist = (int *)malloc(sizeof(int) * nverts)) == NULL)
      error("Failed to allocate celllist");
    repartition->ncelllist = nverts;
  }

  /* We need to rescale the sum of the weights so that the sums of the two
   * types of weights are less than IDX_MAX, that is the range of idx_t.  */
  double vsum = 0.0;
  if (vweights)
    for (int k = 0; k < nverts; k++) vsum += weights_v[k];
  double esum = 0.0;
  if (eweights)
    for (int k = 0; k < nedges; k++) esum += weights_e[k];

  /* Do the scaling, if needed, keeping both weights in proportion. */
  double vscale = 1.0;
  double escale = 1.0;
  if (vweights && eweights) {
    if (vsum > esum) {
      if (vsum > (double)IDX_MAX) {
        vscale = (double)(IDX_MAX - 10000) / vsum;
        escale = vscale;
      }
    } else {
      if (esum > (double)IDX_MAX) {
        escale = (double)(IDX_MAX - 10000) / esum;
        vscale = escale;
      }
    }
  } else if (vweights) {
    if (vsum > (double)IDX_MAX) {
      vscale = (double)(IDX_MAX - 10000) / vsum;
    }
  } else if (eweights) {
    if (esum > (double)IDX_MAX) {
      escale = (double)(IDX_MAX - 10000) / esum;
    }
  }

  if (vweights && vscale != 1.0) {
    vsum = 0.0;
    for (int k = 0; k < nverts; k++) {
      weights_v[k] *= vscale;
      vsum += weights_v[k];
    }
    vscale = 1.0;
  }
  if (eweights && escale != 1.0) {
    esum = 0.0;
    for (int k = 0; k < nedges; k++) {
      weights_e[k] *= escale;
      esum += weights_e[k];
    }
    escale = 1.0;
  }

  /* Balance edges and vertices when the edge weights are timebins, as these
   * have no reason to have equivalent scales, we use an equipartition. */
  if (timebins && eweights) {

    /* Make sums the same. */
    if (vsum > esum) {
      escale = vsum / esum;
      for (int k = 0; k < nedges; k++) weights_e[k] *= escale;
    } else {
      vscale = esum / vsum;
      for (int k = 0; k < nverts; k++) weights_v[k] *= vscale;
    }
  }

  /* And repartition/ partition, using both weights or not as requested. */
#ifdef HAVE_PARMETIS
  if (repartition->usemetis) {
    pick_metis(nodeID, s, nr_nodes, nverts, nedges, weights_v, weights_e,
               repartition->celllist, 0, NULL);
  } else {
    pick_parmetis(nodeID, s, nr_nodes, nverts, nedges, weights_v, weights_e,
                  refine, repartition->adaptive, repartition->itr,
                  repartition->celllist, 0, NULL);
  }
#else
  pick_metis(nodeID, s, nr_nodes, nverts, nedges, weights_v, weights_e,
             repartition->celllist, 0, NULL);
#endif

  /* Check that all cells have good values. All nodes have same copy, so just
   * check on one. */
  if (nodeID == 0) {
    for (int k = 0; k < nverts; k++)
      if (repartition->celllist[k] < 0 || repartition->celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", repartition->celllist[k], k);
  }

  /* Check that the partition is complete and all nodes have some work. */
  int present[nr_nodes];
  int failed = 0;
  for (int i = 0; i < nr_nodes; i++) present[i] = 0;
  for (int i = 0; i < nverts; i++) present[repartition->celllist[i]]++;
  for (int i = 0; i < nr_nodes; i++) {
    if (!present[i]) {
      failed = 1;
      if (nodeID == 0) message("Node %d is not present after repartition", i);
    }
  }

  /* If partition failed continue with the current one, but make this clear. */
  if (failed) {
    if (nodeID == 0)
      message(
          "WARNING: repartition has failed, continuing with the current"
          " partition, load balance will not be optimal");
    for (int k = 0; k < nverts; k++)
      repartition->celllist[k] = cells[k].nodeID;
  }

  /* And apply to our cells */
  split_metis_zoom(s, nr_nodes, repartition->celllist, nverts, 0);

  /* Clean up. */
  free(inds);
  if (vweights) free(weights_v);
  if (eweights) free(weights_e);
}
#endif

/**
 * @brief Repartition the space using the given repartition type.
 *
 * Note that at the end of this process all the cells will be re-distributed
 * across the nodes, but the particles themselves will not be.
 *
 * @param reparttype #repartition struct
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
void partition_repartition_zoom(struct repartition *reparttype, int nodeID,
                                int nr_nodes, struct space *s,
                                struct task *tasks,
                                int nr_tasks) {

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  ticks tic = getticks();

  /* Before doing anything else we need to do the zoom cells. */
  if (reparttype->type == REPART_METIS_VERTEX_EDGE_COSTS) {
    repart_edge_metis_zoom(1, 1, 0, reparttype, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype->type == REPART_METIS_EDGE_COSTS) {
    repart_edge_metis_zoom(0, 1, 0, reparttype, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COSTS_TIMEBINS) {

    /* We can't use timecosts with a zoom decomp. */
    error("Repartition type 'timecosts' is incompatible with a zoom region");

  } else if (reparttype->type == REPART_METIS_VERTEX_COUNTS) {
    repart_memory_metis_zoom(reparttype, nodeID, nr_nodes, s);

  } else if (reparttype->type == REPART_NONE) {
    /* Doing nothing. */

  } else {
    error("Impossible repartition type");
  }

  /* Now handle the background in the desired way */
  if (s->zoom_props->separate_decomps) {

    /* Here we will decompose each level separately. */
    
  } else if (s->zoom_props->use_bkg_wedges) {

    /* Here we will decompose the background first into wedges and then
     * decompose those wedges onto each rank. */
      
  } else {

    /* Put all background cells on a single node */
    for (int cid = s->zoom_props->nr_zoom_cells; cid < s->nr_cells; cid++) {
      s->cells_top[cid].nodeID = 0;
    }
    
  }

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with METIS or ParMETIS support.");
#endif
}

/**
 * @brief Initial partition of space cells.
 *
 * Cells are assigned to a node on the basis of various schemes, all of which
 * should attempt to distribute them in geometrically close regions to
 * minimise the movement of particles.
 *
 * Note that the partition type is a suggestion and will be ignored if that
 * scheme fails. In that case we fallback to a vectorised scheme, that is
 * guaranteed to work provided we have more cells than nodes.
 *
 * @param initial_partition the type of partitioning to try.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells.
 */
void partition_initial_partition_zoom(struct partition *initial_partition,
                                 int nodeID, int nr_nodes, struct space *s) {
  ticks tic = getticks();

  /* Geometric grid partitioning. */
  if (initial_partition->type == INITPART_GRID) {
      error("Grid is not compatible with a zoom region!");
  } else if (initial_partition->type == INITPART_RADIAL) {
#if defined(WITH_MPI)

    /* How many wedges do we have? Start by treating each cell as an area on the
     * spheres surface. */
    int nwedges = s->zoom_props->nwedges;
    int nslices = sqrt(nwedges);
    nwedges = nslices * nslices;
    
    /* Particles sizes per cell, which will be used as weights. */
    double *weights_v = NULL;
    if ((weights_v = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
      error("Failed to allocate weights_v buffer.");

    /* Check each particle and accumulate the sizes per cell. */
    accumulate_sizes(s, s->e->verbose, weights_v);

    /* Do a simple radial wedge decomposition. */
    split_radial_wedges(s, nr_nodes, weights_v, nslices, nwedges);

    /* The radial technique shouldn't fail, but lets be safe. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("Grid initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
      return;
    }
    
#endif

  } else if (s->with_zoom_region &&
             (initial_partition->type == INITPART_METIS_WEIGHT ||
              initial_partition->type == INITPART_METIS_WEIGHT_EDGE ||
              initial_partition->type == INITPART_METIS_NOWEIGHT)) {
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
    /* Simple k-way partition selected by METIS using cell particle
     * counts as weights or not. Should be best when starting with a
     * inhomogeneous dist.
     *
     * There are multiple ways of treating the background cells but we will
     * always have Metis handle the zoom cells.
     */
    
    /* First get the particle weights for all cells. */
    double *cell_weights;
    if ((cell_weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate cell_weights buffer.");
    accumulate_sizes(s, s->e->verbose, cell_weights);

    /* Declare variables we'll need at all levels */
    int *cdim;
    int nverts, nedges, offset;
    
    /* ======================= Zoom Cells ======================= */

    message("Partitioning Zoom cells...");
      
    /* Define the number of vertices */
    nverts = s->zoom_props->nr_zoom_cells;

    /* Get the cell's offset. */
    offset = 0;

    /* Define the number of edges we have to handle. */
    nedges = 0;
    for (int cid = 0; cid < nverts; cid++) {
      nedges += s->cells_top[cid].nr_vertex_edges;
    }

    /* Get this levels cdim. */
    cdim = s->zoom_props->cdim;

    double *zoom_weights_v = NULL;
    double *zoom_weights_e = NULL;
    double sum = 0.0;
    if (initial_partition->type == INITPART_METIS_WEIGHT) {

      /* Particles sizes per cell or wedge, which will be used as weights. */
      if ((zoom_weights_v = (double *)
           malloc(sizeof(double) * nverts)) == NULL)
        error("Failed to allocate zoom_weights_v buffer.");
      bzero(zoom_weights_v, nverts * sizeof(double));

      /* Get the zoom cell weights. */
      for (int cid = offset; cid < offset + nverts; cid++) {
        zoom_weights_v[cid - offset] = cell_weights[cid];
        sum += cell_weights[cid];
      }

      /* Keep the sum of particles across all ranks in the range of IDX_MAX. */
      if (sum > (double)(IDX_MAX - 10000)) {
        double vscale = (double)(IDX_MAX - 10000) / sum;
        for (int k = 0; k < nverts; k++) zoom_weights_v[k] *= vscale;
      }

    } else if (initial_partition->type == INITPART_METIS_WEIGHT_EDGE) {
      
      /* Particle sizes also counted towards the edges. */
      if ((zoom_weights_v = (double *)
           malloc(sizeof(double) * nverts)) == NULL)
        error("Failed to allocate zoom_weights_v buffer.");
      bzero(zoom_weights_v, sizeof(double) * nverts);
      if ((zoom_weights_e = (double *)malloc(sizeof(double) * nedges)) == NULL)
        error("Failed to allocate zoom_weights_e buffer.");
      bzero(zoom_weights_e, sizeof(double) * nedges);

      /* Get the zoom cell weights. */
      for (int cid = offset; cid < offset + nverts; cid++) {
        zoom_weights_v[cid - offset] = cell_weights[cid];
        sum += cell_weights[cid];
      }

      /* Keep the sum of particles across all ranks in the range of IDX_MAX. */
      if (sum > (double)(IDX_MAX - 10000)) {
        double vscale = (double)(IDX_MAX - 10000) / sum;
        for (int k = 0; k < nverts; k++) zoom_weights_v[k] *= vscale;
      }

      /* Spread these into edge weights. */
      sizes_to_edges_zoom(s, zoom_weights_v, zoom_weights_e, offset, cdim);
    }

#ifdef SWIFT_DEBUG_CHECKS
    for (int i = 0; i < nverts; i++) {
      if (!(zoom_weights_v[i] >= 0))
        error("Found zero weighted cell. (i=%d, zoom_weights_e[i]=%.2f)",
              i, zoom_weights_v[i]);
    }
    if (zoom_weights_e != NULL) {
      for (int i = 0; i < nedges; i++) {
        if (!(zoom_weights_e[i] >= 0))
          error("Found zero weighted edge. (i=%d, zoom_weights_e[i]=%.2f)", i,
                zoom_weights_e[i]);
      }
    }
#endif

    /* Do the calculation. */
    int *zoom_celllist = NULL;
    if ((zoom_celllist = (int *)malloc(sizeof(int) * nverts)) == NULL)
      error("Failed to allocate zoom_celllist");
#ifdef HAVE_PARMETIS
    if (initial_partition->usemetis) {
      pick_metis(nodeID, s, nr_nodes, nverts, nedges, zoom_weights_v,
                 zoom_weights_e, zoom_celllist, offset, cdim);
    } else {
      pick_parmetis(nodeID, s, nr_nodes, nverts, nedges, zoom_weights_v,
                    zoom_weights_e, 0, 0, 0.0f, zoom_celllist, offset, cdim);
    }
#else
    pick_metis(nodeID, s, nr_nodes, nverts, nedges, zoom_weights_v,
               zoom_weights_e, zoom_celllist, offset, cdim);
#endif

    /* And apply to our cells */
    split_metis_zoom(s, nr_nodes, zoom_celllist, nverts, offset);

    free(zoom_celllist);

    message("Completed partitioning zoom cells with %d vertices and %d edges",
            nverts, nedges);

    /* Now handle the background in the desired way */
    if (s->zoom_props->separate_decomps) {

      /* Here we will decompose each level separately. */
      
    } else if (s->zoom_props->use_bkg_wedges) {

      /* Here we will decompose the background first into wedges and then
       * decompose those wedges onto each rank. */
      
    } else {

      /* Put all background cells on a single node */
      for (int cid = s->zoom_props->nr_zoom_cells; cid < s->nr_cells; cid++) {
        s->cells_top[cid].nodeID = 0;
      }
      
    }

    /* It's not known if this can fail, but check for this before
     * proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("METIS initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
    }

#else
    error("SWIFT was not compiled with METIS or ParMETIS support");
#endif

  } else if (initial_partition->type == INITPART_VECTORIZE) {

#if defined(WITH_MPI)
    /* Vectorised selection, guaranteed to work for samples less than the
     * number of cells, but not very clumpy in the selection of regions. */
    int *samplecells = NULL;
    if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
      error("Failed to allocate samplecells");

    if (nodeID == 0) {
      pick_vector(s, s->cdim, nr_nodes, samplecells);
    }

    /* Share the samplecells around all the nodes. */
    int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the partition sample cells.");

    /* Do the zoom cells if we are running with them */
    if (s->with_zoom_region) {

      /* With a zoom region we must apply the background offset */
      split_vector(s, s->cdim, nr_nodes, samplecells,
                   s->zoom_props->bkg_cell_offset);
      free(samplecells);

      int *zoom_samplecells = NULL;
      if ((zoom_samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) ==
          NULL)
        error("Failed to allocate zoom_samplecells");

      if (nodeID == 0) {
        pick_vector(s, s->zoom_props->cdim, nr_nodes, zoom_samplecells);
      }

      /* Share the zoom_samplecells around all the nodes. */
      res =
          MPI_Bcast(zoom_samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
      if (res != MPI_SUCCESS)
        mpi_error(res, "Failed to bcast the partition sample cells.");

      /* And apply to our zoom cells */
      split_vector(s, s->zoom_props->cdim, nr_nodes, zoom_samplecells, 0);
      free(zoom_samplecells);
    }

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure everyone agrees how many cells they should have. */

  /* Set up array to hold cell counts. */
  int ncells_on_rank[nr_nodes];
  for (int rank = 0; rank < nr_nodes; rank++)
    ncells_on_rank[rank] = 0;

  /* Count cells on each rank. */
  for (int ind = 0; ind < s->nr_cells; ind++)
    ncells_on_rank[s->cells_top[ind].nodeID]++;

  /* Set up array to hold everyone's cell counts. */
  int rank_cell_counts[nr_nodes * nr_nodes];
  for (int ind = 0; ind < nr_nodes * nr_nodes; ind++)
    rank_cell_counts[ind] = 0;
  for (int irank = 0; irank < nr_nodes; irank++) {
    if (irank == nodeID)
      for (int jrank = 0; jrank < nr_nodes; jrank++)
        rank_cell_counts[irank * nr_nodes + jrank] = ncells_on_rank[jrank];
  }
  
  /* Tell everyone what we've found. */
  for (int rank = 0; rank < nr_nodes; rank++) {
    res =
      MPI_Bcast(&rank_cell_counts[rank], nr_nodes, MPI_INT, rank, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the cell counts of rank %d.", rank);
  }

  /* Let's check we all agree. */
  for (int jrank = 0; jrank < nr_nodes; jrank++) {
    for (int irank = 0; irank < nr_nodes; irank++) {
      if (rank_cell_counts[irank * nr_nodes + jrank] != ncells_on_rank[jrank])
        error("Rank %d disagrees with rank %d about how many cells it "
              "should have (rank %d = %d, rank %d = %d)", nodeID, irank,
              nodeID, irank,
              ncells_on_rank[jrank],
              rank_cell_counts[irank * nr_nodes + jrank]);
    }
  }
  
#endif
#else
    error("SWIFT was not compiled with MPI support");
#endif /* WITH_MPI */
  }

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
