/* Config parameters. */
#include "../config.h"

#include <float.h>

#include "cell.h"
#include "engine.h"
#include "proxy.h"
#include "space.h"
#include "zoom_region.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Define some values, shouldn't need to change these. */
#define zoom_boost_factor 1.1 // Multiply zoom region by this to give a buffer.
#define neighbour_cell_delta 1 // How many layers of neighbours do we go out?

/**
 * @brief Read parameter file for "ZoomRegion" properties, and initialize the zoom_region struct.
 *
 * @param params Swift parameter structure.
 * @param s The space
 */
void zoom_region_init(struct swift_params *params, struct space *s) {
#ifdef WITH_ZOOM_REGION
  /* Are we running with a zoom region? */
  s->with_zoom_region = parser_get_opt_param_int(params, "ZoomRegion:enable", 0);

  /* If so... */
  if (s->with_zoom_region) {
    /* Zoom region properties are stored in a structure. */
    s->zoom_props = (struct zoom_region_properties *)malloc(
        sizeof(struct zoom_region_properties));
    if (s->zoom_props == NULL)
      error("Error allocating memory for the zoom parameters.");

  }
#endif
}

/**
 * @brief For a given particle location, what TL cell does it belong to?
 *
 * Slightly more complicated in the zoom case, as there are now two embedded TL grids.
 * 
 * First see if the particle is within the zoom bounds, then find its TL cell. 
 *
 * @param cdim Cell dimentions of the TL grid (same for natural TL and zoom grid).
 * @param x, y, z Location of particle.
 * @param i, j, k Location of particle (as integer fraction of box, 0 to 1).
 * @param s The space.
 */
int cell_getid_zoom(const int cdim[3], const double x, const double y,
                    const double z, const struct space *s,
                    const int i, const int j, const int k) {
#ifdef WITH_ZOOM_REGION
  int cell_id;

  if (s->with_zoom_region) {

    /* Properties of the zoom region. */
    const struct zoom_region_properties *zoom_props = s->zoom_props;
    const int zoom_cell_offset = zoom_props->tl_cell_offset;
    const double zoom_region_bounds[6] = {
        zoom_props->region_bounds[0], zoom_props->region_bounds[1],
        zoom_props->region_bounds[2], zoom_props->region_bounds[3],
        zoom_props->region_bounds[4], zoom_props->region_bounds[5]};
    const double ih_x_zoom = zoom_props->iwidth[0];
    const double ih_y_zoom = zoom_props->iwidth[1];
    const double ih_z_zoom = zoom_props->iwidth[2];

    /* Are the passed coordinates within the zoom region? */
    if (x > zoom_region_bounds[0] && x < zoom_region_bounds[1] &&
        y > zoom_region_bounds[2] && y < zoom_region_bounds[3] &&
        z > zoom_region_bounds[4] && z < zoom_region_bounds[5]) {
    
      /* Which zoom TL cell are we in? */
      const int zoom_index =
          cell_getid(cdim, (x - zoom_region_bounds[0]) * ih_x_zoom,
                     (y - zoom_region_bounds[2]) * ih_y_zoom,
                     (z - zoom_region_bounds[4]) * ih_z_zoom);
      cell_id = zoom_cell_offset + zoom_index;
#ifdef SWIFT_DEBUG_CHECKS
      if (zoom_index < 0 || zoom_index >= cdim[0] * cdim[1] * cdim[2])
        error("zoom_index out of range %i (%f %f %f)", cell_id, x, y, z);
#endif
    /* If not then treat it like normal, and find the natural TL cell. */
    } else {
        cell_id = cell_getid(cdim, i, j, k);
    }
  /* Normal non-zoom-region case. */
  } else {
    cell_id = cell_getid(cdim, i, j, k);
  }

  return cell_id;
#else
  return 0;
#endif
}

/**
 * @brief Compute the extent/bounds of the zoom region using the high-res DM particles.
 *
 * The min/max [x,y,z] for each particle is found, and the CoM of these particles is computed.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void construct_zoom_region(struct space *s, int verbose) {
#ifdef WITH_ZOOM_REGION
  double new_zoom_boundary[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
  const size_t nr_gparts = s->nr_gparts;
  double mtot = 0.0;
  double com[3] = {0.0, 0.0, 0.0};

  /* Find the min/max location in each dimension for each mask particle, and their COM. */
  for (size_t k = 0; k < nr_gparts; k++) {
    if (s->gparts[k].type != swift_type_dark_matter) continue;

    if (s->gparts[k].x[0] < new_zoom_boundary[0])
      new_zoom_boundary[0] = s->gparts[k].x[0];
    if (s->gparts[k].x[0] > new_zoom_boundary[1])
      new_zoom_boundary[1] = s->gparts[k].x[0];
    if (s->gparts[k].x[1] < new_zoom_boundary[2])
      new_zoom_boundary[2] = s->gparts[k].x[1];
    if (s->gparts[k].x[1] > new_zoom_boundary[3])
      new_zoom_boundary[3] = s->gparts[k].x[1];
    if (s->gparts[k].x[2] < new_zoom_boundary[4])
      new_zoom_boundary[4] = s->gparts[k].x[2];
    if (s->gparts[k].x[2] > new_zoom_boundary[5])
      new_zoom_boundary[5] = s->gparts[k].x[2];

    mtot += s->gparts[k].mass;
    com[0] += s->gparts[k].x[0] * s->gparts[k].mass;
    com[1] += s->gparts[k].x[1] * s->gparts[k].mass;
    com[2] += s->gparts[k].x[2] * s->gparts[k].mass;
  }

#ifdef WITH_MPI
  /* Share answers amoungst nodes. */

  /* Boundary. */
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[0], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[1], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[2], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[3], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[4], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[5], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  /* CoM. */
  MPI_Allreduce(MPI_IN_PLACE, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /* Finalize CoM calcuation. */
  const double imass = 1.0 / mtot;
  com[0] *= imass;
  com[1] *= imass;
  com[2] *= imass;

  /* Store result. */
  for (int k = 0; k < 3; k++) s->zoom_props->com[k] = com[k];

  s->zoom_props->dim[0] = (new_zoom_boundary[1] - new_zoom_boundary[0]) * zoom_boost_factor;
  s->zoom_props->dim[1] = (new_zoom_boundary[3] - new_zoom_boundary[2]) * zoom_boost_factor;
  s->zoom_props->dim[2] = (new_zoom_boundary[5] - new_zoom_boundary[4]) * zoom_boost_factor;

  if (verbose)
    message("com: [%f %f %f] dim: [%f %f %f]",
          com[0], com[1], com[2], s->zoom_props->dim[0], s->zoom_props->dim[1],
          s->zoom_props->dim[2]);
#endif
}

/**
 * @brief Build the TL cells, with a zoom region.
 *
 * This replaces the loop in space_regrid when running with a zoom region.
 * 
 * Construct an additional set of TL "zoom" cells embedded within the TL cell structure
 * with the dimensions of each cell structure being the same (with differing widths).
 *
 * Therefore the new TL cell structure is 2*cdim**3, with the "natural" TL cells ocupying the 
 * first half of the TL cell list, and the "zoom" TL cells ocupying the second half.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void construct_tl_cells_with_zoom_region(struct space *s, const int *cdim, const float dmin,
        const integertime_t ti_current, const int verbose) {
#ifdef WITH_ZOOM_REGION
  
  /* We are recomputing the boundary of the zoom region. */
  double zoom_region_bounds[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
  const int zoom_cell_offset = cdim[0] * cdim[1] * cdim[2];
  float dmin_zoom = 0.f;

  /* Loop over top level cells twice, second time is for the zoom region. */
  for (int n = 0; n < 2; n++) {
    if (n == 1 && !s->with_zoom_region) continue;

    /* Set the cell location and sizes. */
    for (int i = 0; i < cdim[0]; i++)
      for (int j = 0; j < cdim[1]; j++)
        for (int k = 0; k < cdim[2]; k++) {
          const size_t cid = cell_getid(cdim, i, j, k);

          struct cell *restrict c;
          if (n == 0) {
            /* Natural top level cells. */
            c = &s->cells_top[cid];
            c->loc[0] = i * s->width[0];
            c->loc[1] = j * s->width[1];
            c->loc[2] = k * s->width[2];
            c->width[0] = s->width[0];
            c->width[1] = s->width[1];
            c->width[2] = s->width[2];
            if (s->with_self_gravity)
              c->grav.multipole = &s->multipoles_top[cid];
            c->tl_cell_type = tl_cell;
            c->dmin = dmin;
          } else {
            /* Zoom region top level cells. */
            c = &s->cells_top[cid + zoom_cell_offset];
            c->loc[0] = i * s->zoom_props->width[0] + zoom_region_bounds[0];
            c->loc[1] = j * s->zoom_props->width[1] + zoom_region_bounds[2];
            c->loc[2] = k * s->zoom_props->width[2] + zoom_region_bounds[4];
            c->width[0] = s->zoom_props->width[0];
            c->width[1] = s->zoom_props->width[1];
            c->width[2] = s->zoom_props->width[2];
            if (s->with_self_gravity)
              c->grav.multipole = &s->multipoles_top[cid + zoom_cell_offset];
            c->tl_cell_type = zoom_tl_cell;
            c->dmin = dmin_zoom;
          }
          c->depth = 0;
          c->split = 0;
          c->hydro.count = 0;
          c->grav.count = 0;
          c->stars.count = 0;
          c->sinks.count = 0;
          c->top = c;
          c->super = c;
          c->hydro.super = c;
          c->grav.super = c;
          c->hydro.ti_old_part = ti_current;
          c->grav.ti_old_part = ti_current;
          c->stars.ti_old_part = ti_current;
          c->sinks.ti_old_part = ti_current;
          c->black_holes.ti_old_part = ti_current;
          c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
          c->mpi.tag = -1;
          c->mpi.recv = NULL;
          c->mpi.send = NULL;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
          cell_assign_top_level_cell_index(c, s->cdim, s->dim, s->iwidth);
#endif
          /* Only do what comes next first time round. */
          if (n == 1) continue;
          /* Is this top level cell within the zoom region? */
          if (s->with_zoom_region &&
              (c->loc[0] + c->width[0] >
               s->zoom_props->com[0] - s->zoom_props->dim[0] / 2.) &&
              (c->loc[0] <
               s->zoom_props->com[0] + s->zoom_props->dim[0] / 2.) &&
              (c->loc[1] + c->width[1] >
               s->zoom_props->com[1] - s->zoom_props->dim[1] / 2.) &&
              (c->loc[1] <
               s->zoom_props->com[1] + s->zoom_props->dim[1] / 2.) &&
              (c->loc[2] + c->width[2] >
               s->zoom_props->com[2] - s->zoom_props->dim[2] / 2.) &&
              (c->loc[2] <
               s->zoom_props->com[2] + s->zoom_props->dim[2] / 2.)) {

            /* Tag this top level cell as part of the zoom region. */
            c->tl_cell_type = void_tl_cell;

            /* Update the bounds of the zoom region. */
            if (c->loc[0] < zoom_region_bounds[0])
              zoom_region_bounds[0] = c->loc[0];
            if (c->loc[0] + c->width[0] > zoom_region_bounds[1])
              zoom_region_bounds[1] = c->loc[0] + c->width[0];
            if (c->loc[1] < zoom_region_bounds[2])
              zoom_region_bounds[2] = c->loc[1];
            if (c->loc[1] + c->width[1] > zoom_region_bounds[3])
              zoom_region_bounds[3] = c->loc[1] + c->width[1];
            if (c->loc[2] < zoom_region_bounds[4])
              zoom_region_bounds[4] = c->loc[2];
            if (c->loc[2] + c->width[2] > zoom_region_bounds[5])
              zoom_region_bounds[5] = c->loc[2] + c->width[2];
          }
        }

    /* Compute size of top level zoom cells on first iteration. */
    if (n == 1) continue;
    if (s->with_zoom_region) {
      s->zoom_props->dim[0] = zoom_region_bounds[1] - zoom_region_bounds[0];
      s->zoom_props->dim[1] = zoom_region_bounds[3] - zoom_region_bounds[2];
      s->zoom_props->dim[2] = zoom_region_bounds[5] - zoom_region_bounds[4];

      for (int l = 0; l < 3; l++) {
        s->zoom_props->width[l] = s->zoom_props->dim[l] / cdim[l];
        s->zoom_props->iwidth[l] = 1 / s->zoom_props->width[l];
        s->zoom_props->cdim[l] = cdim[l];
      }

      for (int l = 0; l < 6; l++)
        s->zoom_props->region_bounds[l] = zoom_region_bounds[l];
      s->zoom_props->tl_cell_offset = zoom_cell_offset;
      dmin_zoom = min3(s->zoom_props->width[0], s->zoom_props->width[1],
                       s->zoom_props->width[2]);
    }
  }

  /* Now find what cells neighbour the zoom region. */
  if (s->with_zoom_region) find_neighbouring_cells(s, verbose);

#endif
}

/**
 * @brief Find what TL cells surround the zoom region.
 *
 * When interacting "natural" TL cells and "zoom" TL cells, it helps to know what natural TL
 * cells surround the zoom region. These cells then get tagged as "tl_cell_neighbour".
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_neighbouring_cells(struct space *s, const int verbose) {
#ifdef WITH_ZOOM_REGION
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const int periodic = s->periodic;
  struct cell *cells = s->cells_top;

  const int delta_m = neighbour_cell_delta; // Should compute this, but how?
  const int delta_p = neighbour_cell_delta;

  int neighbour_count = 0;

  /* Loop over each cell in the space to find the neighbouring top level cells
   * surrounding the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Only interested in cells hosting zoom top level cells. */
        if (cells[cid].tl_cell_type != void_tl_cell) continue;

        /* Loop over all its direct neighbours. */
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

              /* Get the cell ID of the neighbour. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              if (cells[cjd].tl_cell_type == tl_cell) {

                /* Record that we've found a neighbour. */
                cells[cjd].tl_cell_type = tl_cell_neighbour;
                neighbour_count++;
              }
            }
          }
        }
      }
    }
  }

  if (verbose)
    message("%i cells neighbouring the zoom region", neighbour_count);
#endif
}

/**
 * @brief Minimum distance between two TL cells with different sizes.
 *
 * @param ci, cj The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
double cell_min_dist2_diff_size(const struct cell *restrict ci,
                                const struct cell *restrict cj,
                                const int periodic, const double dim[3]) {
#ifdef WITH_ZOOM_REGION

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->width[0] == cj->width[0]) error("x cells of same size!");
  if (ci->width[1] == cj->width[1]) error("y cells of same size!");
  if (ci->width[2] == cj->width[2]) error("z cells of same size!");
#endif

  const double cix = ci->loc[0] + ci->width[0]/2.;
  const double ciy = ci->loc[1] + ci->width[1]/2.;
  const double ciz = ci->loc[2] + ci->width[2]/2.;

  const double cjx = cj->loc[0] + cj->width[0]/2.;
  const double cjy = cj->loc[1] + cj->width[1]/2.;
  const double cjz = cj->loc[2] + cj->width[2]/2.;

  const double diag_ci2 = ci->width[0] * ci->width[0] + ci->width[1] * ci->width[1] + ci->width[2] * ci->width[2];
  const double diag_cj2 = cj->width[0] * cj->width[0] + cj->width[1] * cj->width[1] + cj->width[2] * cj->width[2];

  /* Get the distance between the cells */
  double dx = cix - cjx;
  double dy = ciy - cjy;
  double dz = ciz - cjz;

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  /* Minimal distance between any 2 particles in the two cells */
  const double dist2 = r2 - (diag_ci2/2. + diag_cj2/2.);

  return dist2; 
#else
  return 0;
#endif
}

/**
 * @brief Compute pair tasks between natural TL cells and zoom TL cells.
 *
 * Loop over the zoom TL cells, then loop over the "neighbours".
 *
 * Apply standard MM/max_distance check to see if we need a pair gravity
 * task between the two cells.
 */
void engine_make_self_gravity_tasks_mapper_between_toplevels(void *map_data,
                                                             int num_elements,
                                                             void *extra_data) {
#ifdef WITH_ZOOM_REGION
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;
  const int nr_cells = s->nr_cells;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Only compute for zoom TL cells. */
    if (cid < s->zoom_props->tl_cell_offset) continue;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    for (int cjd = 0; cjd < nr_cells; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

      /* Only considering "neighbours". */
      if (cj->tl_cell_type != tl_cell_neighbour) continue;

      /* Avoid duplicates, empty cells and completely foreign pairs */
      if (cj->grav.count == 0 || (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      /* Recover the multipole information */
      const struct gravity_tensors *multi_i = ci->grav.multipole;
      const struct gravity_tensors *multi_j = cj->grav.multipole;

      if (multi_i == NULL && ci->nodeID != nodeID)
        error("Multipole of ci was not exchanged properly via the proxies");
      if (multi_j == NULL && cj->nodeID != nodeID)
        error("Multipole of cj was not exchanged properly via the proxies");

      /* Minimal distance between any pair of particles */
      const double min_radius2 =
          cell_min_dist2_diff_size(ci, cj, periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0 ?*/
      if (periodic && min_radius2 > max_distance2) continue;

      /* Are the cells too close for a MM interaction ? */
      if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0)) {

        /* Ok, we need to add a direct pair calculation */
        scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0, ci,
                          cj);
#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI
        /* Let's cross-check that we had a proxy for that cell */
        if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[cj->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", cj->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == cj) {
              break;
            } 
          if (n == p->nr_cells_in)
            error( 
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cjd);
        } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[ci->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", ci->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == ci) {
              break;
            } 
          if (n == p->nr_cells_in)
            error( 
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cid);
        }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */

      }
    } 
  }
#endif
}

/**
 * @brief Minimum distance between two TL cells.
 *
 * Generic wrapper, don't know if the TL cells are the same size or not at time of calling.
 *
 * @param ci, cj The two TL cells.
 * @param periodic Account for periodicity?
 * @param dim The boxsize.
 */
double cell_min_dist2(const struct cell *restrict ci,
                      const struct cell *restrict cj, const int periodic,
                      const double dim[3]) {
#ifdef WITH_ZOOM_REGION
  double dist2;

  /* Two natural TL cells. */
  if (ci->tl_cell_type <= 2 && cj->tl_cell_type <= 2) {
    dist2 = cell_min_dist2_same_size(ci, cj, periodic, dim);
  /* Two zoom TL cells. */
  } else if (ci->tl_cell_type == zoom_tl_cell &&
             cj->tl_cell_type == zoom_tl_cell) {
    dist2 = cell_min_dist2_same_size(ci, cj, periodic, dim);
  /* A mix of natural and zoom TL cells. */
  } else {
    dist2 = cell_min_dist2_diff_size(ci, cj, periodic, dim);
  }

  return dist2;
#else
  return 0;
#endif
}
