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
               int *iedge, int nslices, float slice_width,
               double *slice_weights) {

  /* Declare some variables. */
  struct cell *restrict ci;
  struct cell *restrict cj;
  int phi_ind, theta_ind;
  double r, theta, phi;

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

        /* Which wedge is this cell in? */

        /* Get the centred ijk coordinates. */
        int ii = i - (s->zoom_props->cdim[0] / 2);
        int jj = j - (s->zoom_props->cdim[1] / 2);
        int kk = k - (s->zoom_props->cdim[2] / 2);

        /* Calculate the spherical version of these coordinates. */
        r = sqrt(ii * ii + jj * jj + kk * kk);
        theta = atan2(jj, ii) + M_PI;
        phi = acos(kk / r);

        /* Add this cells weight. */
        phi_ind = phi / slice_width / 2;
        theta_ind = theta / slice_width;
        int iwedge_ind = phi_ind * nslices + theta_ind;

        /* Loop over cells neighbouring the zoom region. */
        for (int n = 0; n < s->zoom_props->nr_neighbour_cells; n++) {

          /* Get the cell. */
          int cjd = s->zoom_props->neighbour_cells_top[n];
          cj = &s->cells_top[cjd];

          /* Get the centred coordinates. */
          double xx = cj->loc[0] - (s->dim[0] / 2);
          double yy = cj->loc[1] - (s->dim[1] / 2);
          double zz = cj->loc[2] - (s->dim[2] / 2);

          /* Calculate the spherical version of these coordinates. */
          r = sqrt(xx * xx + yy * yy + zz * zz);
          theta = atan2(yy, xx) + M_PI;
          phi = acos(zz / r);

          /* Add this cells weight. */
          phi_ind = phi / slice_width / 2;
          theta_ind = theta / slice_width;
          int jwedge_ind = phi_ind * nslices + theta_ind;

          /* Is this the same wedge as ci? */
          if (iwedge_ind == jwedge_ind) {

            /* Handle size_to_edges case */
            if (edges != NULL) {
              /* Store this edge. */
              edges[*iedge] = counts[cjd];
              (*iedge)++;
            }
                
            /* Handle graph_init case */
            else if (adjncy != NULL) {
              adjncy[*iedge] = jwedge_ind + s->zoom_props->nr_zoom_cells;
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
        
      } /* k loop */
    } /* j loop */
  } /* i loop */

  /* Now we need to consider the wedges. */

  /* Loop over each surface element. */
  for (int i = 0; i < nslices; i++) {
    for (int j = 0; j < nslices; j++) {

      /* What is the index of this wedge? */
      int cid = j * nslices + i;

      /* Loop over neighbouring elements. */
      for (int ii = i - 1; ii <= i + 1; ii++) {
        for (int jj = j - 1; jj <= i + 1; jj++) {

          /* Wrap around the sphere. */
          const int iii = (ii + nslices) % nslices;
          const int jjj = (jj + nslices) % nslices;

          /* What is the index of this wedge? */
          int cjd = jjj * nslices + iii;

          /* Skip self. */
          if (cid == cjd) continue;
          
          /* Handle size_to_edges case */
          if (edges != NULL) {
            /* Store this edge. */
            edges[*iedge] = slice_weights[cjd + s->zoom_props->nr_zoom_cells];
            (*iedge)++;
          }
                
          /* Handle graph_init case */
          else if (adjncy != NULL) {
            adjncy[*iedge] = cjd + s->zoom_props->nr_zoom_cells;
            (*iedge)++;
          }

          /* Handle find_vertex_edges case */
          else {
            /* If not self record an edge. */
            (*iedge)++;
          }
        }
      }
    }
  }
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

  /* threadpool_map(&e->threadpool, engine_makeproxies_zoom_cells_mapper, NULL, */
  /*                s->zoom_props->nr_zoom_cells, 1, threadpool_auto_chunk_size, */
  /*                e); */
  /* threadpool_map(&e->threadpool, engine_makeproxies_natural_cells_mapper, NULL, */
  /*                s->zoom_props->nr_bkg_cells, 1, threadpool_auto_chunk_size, */
  /*                e); */
  /* threadpool_map(&e->threadpool, engine_makeproxies_between_zoom_bkg_mapper, NULL, */
  /*                s->zoom_props->nr_zoom_cells, 1, threadpool_auto_chunk_size, */
  /*                e); */
  /* if (e->s->zoom_props->with_buffer_cells) { */
  /*   threadpool_map(&e->threadpool, engine_makeproxies_buffer_cells_mapper, NULL, */
  /*                  s->zoom_props->nr_buffer_cells, 1, */
  /*                  threadpool_auto_chunk_size, e); */
  /*   threadpool_map(&e->threadpool, engine_makeproxies_between_buffer_bkg_mapper, */
  /*                  NULL, s->zoom_props->nr_buffer_cells, 1, */
  /*                  threadpool_auto_chunk_size, e); */
  /* } */

  engine_makeproxies_zoom_cells(e);
  engine_makeproxies_buffer_cells(e);
  engine_makeproxies_natural_cells(e);
  engine_makeproxies_between_zoom_bkg(e);
  if (e->s->zoom_props->with_buffer_cells)
    engine_makeproxies_between_buffer_bkg(e);

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Partition the into radial slices.
 *
 * This simply slices the box into wedges along the x-y plane.
 */
void split_bkg_radial_wedges(struct space *s, int nregions,
                             double *slice_weights, double *cell_weights,
                             int nslices, int nwedges, float slice_width) {

  double r, theta, phi;
    
  /* Define variables for selection */
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

  /* Get the weight of each slice*/

  /* Loop over zoom cells and store the weight. */
  for (int i = 0; i < s->zoom_props->cdim[0]; i++) {
    for (int j = 0; j < s->zoom_props->cdim[1]; j++) {
      for (int k = 0; k < s->zoom_props->cdim[2]; k++) {

        /* Get cell ID. */
        const int cid = cell_getid(s->zoom_props->cdim, i, j, k);

        /* Store this cells weight. */
        slice_weights[cid] = cell_weights[cid];
        
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
        int wedge_ind = phi_ind * nslices + theta_ind;
        slice_weights[wedge_ind + s->zoom_props->nr_zoom_cells] +=
          cell_weights[cid];
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
        int wedge_ind = phi_ind * nslices + theta_ind;
        slice_weights[wedge_ind + s->zoom_props->nr_zoom_cells] +=
          cell_weights[cid];
      }
    }
  }

#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
  /* Ensure the weights are in the correct range. */
  double sum = 0.0;
  for (int i = 0; i < s->zoom_props->nr_zoom_cells + nwedges; i++) {
    sum += slice_weights[i];
  }

  /* Keep the sum of particles across all ranks in the range of IDX_MAX. */
  if (sum > (double)(IDX_MAX - 10000)) {
    double vscale = (double)(IDX_MAX - 10000) / sum;
    for (int k = 0; k < s->nr_cells; k++) slice_weights[k] *= vscale;
  }
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
                         int nslices, float slice_width,
                         double *slice_weights) {

    /* Get some useful constants. */
    const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                              s->zoom_props->cdim[2]};
    int iedge = 0;

    /* Find adjacency arrays for zoom cells. */
    edge_loop(zoom_cdim, 0, s, /*adjncy*/ NULL,
              /*xadj*/ NULL, counts, edges, &iedge, nslices, slice_width,
              slice_weights);
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
 */
static void graph_init_zoom(struct space *s, int periodic, idx_t *weights_e,
                            idx_t *adjncy, int *nadjcny, idx_t *xadj,
                            int *nxadj, int nslices, double slice_width) {

  /* Loop over all cells in the space. */
  *nadjcny = 0;

  /* Get some useful constants. */
  const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                            s->zoom_props->cdim[2]};
  int iedge = 0;

  /* Find adjacency arrays for zoom cells. */
  edge_loop(zoom_cdim, 0, s, adjncy, xadj, /*counts*/ NULL, /*edges*/ NULL,
            &iedge, nslices, slice_width, /*slice_weights*/ NULL);

  /* Set the number of adjacncy entries. */
  *nadjcny = iedge;

  /* If given set METIS xadj. */
  if (xadj != NULL) {
    xadj[s->nr_cells] = iedge;
    *nxadj = s->nr_cells;
  }
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
 */
static void pick_parmetis_zoom(int nodeID, struct space *s, int nregions,
                               double *vertexw, double *edgew, int refine,
                               int adaptive, float itr, int *celllist,
                               int nwedges, int nedges,
                               int nslices, double slice_width) {

  int res;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  /* Total number of cells. */
  int ncells = s->zoom_props->nr_zoom_cells + nwedges;  /* Total number of cells. */
  int ncells = s->nr_cells;

  /* Total number of edges. */
  int nedges;
  if (s->with_zoom_region) {
    nedges = s->zoom_props->nr_edges;
  } else {
    nedges = 26 * s->nr_cells;
  }
  
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
  if (s->with_zoom_region) {
    for (int cid = vtxdist[nodeID]; cid < vtxdist[nodeID + 1]; cid++)
      nr_my_edges += s->cells_top[cid].nr_vertex_edges; 
  } else {
    nr_my_edges = nverts * 26;
  }

  idx_t *xadj = NULL;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
    error("Failed to allocate xadj buffer.");

  idx_t *adjncy = NULL;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * nr_my_edges)) == NULL)
    error("Failed to allocate adjncy array.");

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

    /* Define the cell graph. Keeping the edge weights association. */
    int nadjcny = 0;
    int nxadj = 0;
    graph_init_zoom(s, s->periodic, full_weights_e, full_adjncy, &nadjcny,
                    std_xadj, &nxadj, nslices, slice_width);
    
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
    if (!refine) {

      /* No old partition was given, so we need to construct the existing
       * partition from the cells, if one existed. */
      int nsum = 0;
      for (int i = 0; i < s->nr_cells; i++) {
        celllist[i] = s->cells_top[i].nodeID;
        nsum += celllist[i];
      }

      /* If no previous partition then all nodeIDs will be set to 0. */
      if (nsum == 0) permute = 0;
    }

    if (permute) {
      int *permcelllist = NULL;
      if ((permcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
        error("Failed to allocate perm celllist array");
      permute_regions(newcelllist, celllist, nregions, ncells, permcelllist);

      /* And keep. */
      memcpy(celllist, permcelllist, sizeof(int) * ncells);
      free(permcelllist);

    } else {
      memcpy(celllist, newcelllist, sizeof(int) * ncells);
    }
    free(newcelllist);
  }

  /* And everyone gets a copy. */
  res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD);
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
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights. Need to be in the range of idx_t.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array. Need to be in the range of
 *        idx_t.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 */
void pick_metis_zoom(int nodeID, struct space *s, int nregions,
                     double *vertexw, double *edgew, int *celllist,
                     int nwedges, int nedges,
                     int nslices, double slice_width) {

  /* Total number of cells. */
  int ncells = s->zoom_props->nr_zoom_cells + nwedges;

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* Only one node needs to calculate this. */
  if (nodeID == 0) {

    /* Allocate adjacency and weights arrays . */
    idx_t *xadj;
    if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
      error("Failed to allocate xadj buffer.");
    idx_t *adjncy;
    if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * nedges)) == NULL)
      error("Failed to allocate adjncy array.");
    idx_t *weights_v = NULL;
    if (vertexw != NULL)
      if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate vertex weights array");
    idx_t *weights_e = NULL;
    if (edgew != NULL)
      if ((weights_e = (idx_t *)malloc(nedges * sizeof(idx_t))) == NULL)
        error("Failed to allocate edge weights array");
    idx_t *regionid;
    if ((regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate regionid array");

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          weights_v[k] = vertexw[k];
        } else {
          weights_v[k] = 0;
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
               nslices, slice_width);

#ifdef SWIFT_DEBUG_CHECKS

    message("There are %d edges and %d adjacncies", nedges, nadjcny)
    
    /* Check all adjacencies are set. */
    int failed = 0;
    for (int k = 0; k < nedges; k++) {
      if ((idx_t)adjncy[k] < 0) {
        message("No adjacency at %d", k);
        failed++;
      }
    }
    if (failed > 0) error("%d adjacencies not set", failed);
#endif

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
    idx_t idx_ncells = ncells;
    idx_t idx_nregions = nregions;
    idx_t objval;

    /* Dump graph in METIS format */
    /* dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy, weights_v, */
    /*                NULL, weights_e); */

    if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v, NULL,
                            weights_e, &idx_nregions, NULL, NULL, options,
                            &objval, regionid) != METIS_OK)
      error("Call to METIS_PartGraphKway failed.");

    /* Check that the regionids are ok. */
    for (int k = 0; k < ncells; k++) {
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
  int res = MPI_Bcast(celllist, ncells, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");
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
void split_metis_zoom(struct space *s, int nregions, int *celllist,
                      int nslices, double slice_width) {

  const int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* First do the zoom cells. */
  for (int i = 0; i < nr_zoom_cells; i++)
    s->cells_top[i].nodeID = celllist[i];

  /* Now we need to loop over all the background cells and assign them based
   * on their slice. */
  
  double r, theta, phi;
    
  /* Define variables for selection */
  const int bkg_cell_offset = s->zoom_props->tl_cell_offset;
  const int buffer_cell_offset = s->zoom_props->buffer_cell_offset;

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
        int wedge_ind = phi_ind * nslices + theta_ind + nr_zoom_cells;
        s->cells_top[cid].nodeID = celllist[wedge_ind];
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
        int wedge_ind = phi_ind * nslices + theta_ind + nr_zoom_cells;
        s->cells_top[cid].nodeID = celllist[wedge_ind];
      }
    }
  }

  /* To check or visualise the partition dump all the cells. */
  /*if (engine_rank == 0) dumpCellRanks("metis_partition", s->cells_top,
                                      s->nr_cells);*/
}
#endif
