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
  double widths[3] = {0.0, 0.0, 0.0};

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

  /* Assign each axis to array. */
  widths[0] = (new_zoom_boundary[1] - new_zoom_boundary[0]);
  widths[1] = (new_zoom_boundary[3] - new_zoom_boundary[2]);
  widths[2] = (new_zoom_boundary[5] - new_zoom_boundary[4]);

  /* Get the maximum axis length and assign it to the zoom region dimension. */
  double max_width = 0;
  for (int k = 0; k < 3; k++) {
      if (widths[k] > max_width)
          max_width = widths[k];
  }
  s->zoom_props->dim[0] = max_width * zoom_boost_factor;
  s->zoom_props->dim[1] = max_width * zoom_boost_factor;
  s->zoom_props->dim[2] = max_width * zoom_boost_factor;

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
 * Therefore the new TL cell structure is 2*cdim**3, with the "natural" TL cells occupying the
 * first half of the TL cell list, and the "zoom" TL cells ocupying the second half.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void construct_tl_cells_with_zoom_region(struct space *s, const int *cdim, const float dmin,
        const integertime_t ti_current, int verbose) {
#ifdef WITH_ZOOM_REGION
  
  /* We are recomputing the boundary of the zoom region. */
  double zoom_region_bounds[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
  const int zoom_cell_offset = cdim[0] * cdim[1] * cdim[2];
  float dmin_zoom = 0.f;
  double widths[3] = {0.0, 0.0, 0.0};
  double mid_points[3] = {0.0, 0.0, 0.0};

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
          	/* On the second loop we need to assign zoom cell information to natural cells. */
            c = &s->cells_top[cid];
            c->nr_zoom_cells = s->width[0] / s->zoom_props->width[0];

            /* We need to update the cell types after increasing the zoom region size in the first interation */
				    if (c->loc[0] > zoom_region_bounds[0] && c->loc[0] < zoom_region_bounds[1] &&
				        c->loc[1] > zoom_region_bounds[2] && c->loc[1] < zoom_region_bounds[3] &&
				        c->loc[2] > zoom_region_bounds[4] && c->loc[2] < zoom_region_bounds[5]) {
				    	/* Tag this top level cell as part of the zoom region. */
              c->tl_cell_type = void_tl_cell;
				    }

            /* Assign the start and end indices for the zoom cells within this cell */
			      if (c->tl_cell_type == void_tl_cell) {
			        c->start_i = (c->loc[0] - zoom_region_bounds[0]) * s->zoom_props->iwidth[0];
			        c->start_j = (c->loc[1] - zoom_region_bounds[2]) * s->zoom_props->iwidth[1];
			        c->start_k = (c->loc[2] - zoom_region_bounds[4]) * s->zoom_props->iwidth[2];
			        c->end_i = (c->loc[0] - zoom_region_bounds[0] + c->width[0]) * s->zoom_props->iwidth[0];
			        c->end_j = (c->loc[1] - zoom_region_bounds[2] + c->width[1]) * s->zoom_props->iwidth[1];
			        c->end_k = (c->loc[2] - zoom_region_bounds[4] + c->width[2]) * s->zoom_props->iwidth[2];
			      }

            /* Zoom region top level cells. */
            c = &s->cells_top[cid + zoom_cell_offset];
            c->loc[0] = i * s->zoom_props->width[0] + zoom_region_bounds[0];
            c->loc[1] = j * s->zoom_props->width[1] + zoom_region_bounds[2];
            c->loc[2] = k * s->zoom_props->width[2] + zoom_region_bounds[4];
            c->parent_tl_cid = cell_getid(cdim, (int)(c->loc[0] * s->iwidth[0]),
            		(int)(c->loc[1] * s->iwidth[1]), (int)(c->loc[2] * s->iwidth[2]));
            c->width[0] = s->zoom_props->width[0];
            c->width[1] = s->zoom_props->width[1];
            c->width[2] = s->zoom_props->width[2];
            if (s->with_self_gravity)
              c->grav.multipole = &s->multipoles_top[cid + zoom_cell_offset];
            c->tl_cell_type = zoom_tl_cell;
            c->dmin = dmin_zoom;
            c->nr_zoom_cells = s->width[0] / s->zoom_props->width[0];
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

      /* Assign each axis to array. */
      widths[0] = (zoom_region_bounds[1] - zoom_region_bounds[0]);
      widths[1] = (zoom_region_bounds[3] - zoom_region_bounds[2]);
      widths[2] = (zoom_region_bounds[5] - zoom_region_bounds[4]);

      /* Get the maximum axis length of the zoom region. */
      double max_width = 0;
      for (int k = 0; k < 3; k++) {
          if (widths[k] > max_width)
              max_width = widths[k];
      }

      struct cell *restrict cent_c;
      const int cent_cid = cell_getid(cdim, (int)(s->zoom_props->com[0] * s->iwidth[0]),
      		(int)(s->zoom_props->com[1] * s->iwidth[1]), (int)(s->zoom_props->com[2] * s->iwidth[2]));;
      cent_c = &s->cells_top[cent_cid];

      for (int k = 0; k < 3; k++) {
      	mid_points[k] = cent_c->loc[k] + (s->width[k] / 2);
      	zoom_region_bounds[k * 2] = mid_points[k] - (max_width / 2);
      	zoom_region_bounds[(k * 2) + 1] = mid_points[k] + (max_width / 2);
      }

      /* Overwrite zoom region properties. */
      s->zoom_props->dim[0] = max_width;
      s->zoom_props->dim[1] = max_width;
      s->zoom_props->dim[2] = max_width;
      for (int l = 0; l < 3; l++) {
        s->zoom_props->width[l] = s->zoom_props->dim[l] / cdim[l];
        s->zoom_props->iwidth[l] = 1 / s->zoom_props->width[l];
        s->zoom_props->cdim[l] = cdim[l];
      }

      if (verbose) {
        message("zoom_region_cent: [%f %f %f]",
        		mid_points[0], mid_points[1], mid_points[2]);
        message("tl_cell_width: [%f %f %f] zoom_cell_width: [%f %f %f] dim: [%f %f %f]",
	      		s->width[0], s->width[1], s->width[2],
	      		s->zoom_props->width[0], s->zoom_props->width[1], s->zoom_props->width[2],
	      		max_width, max_width, max_width);
        message("nr_tl_cells_in_zoom_region: [%f %f %f] nr_zoom_cells_in_tl_cell: [%f %f %f]",
        		max_width / s->width[0], max_width / s->width[1], max_width / s->width[2],
        		s->width[0] / s->zoom_props->width[0], s->width[1] / s->zoom_props->width[1],
        		s->width[2] / s->zoom_props->width[2]);
      }

      /* Set the number of zoom cells in a natural cell */

      for (int l = 0; l < 6; l++)
        s->zoom_props->region_bounds[l] = zoom_region_bounds[l];
      s->zoom_props->tl_cell_offset = zoom_cell_offset;
      dmin_zoom = min3(s->zoom_props->width[0], s->zoom_props->width[1],
                       s->zoom_props->width[2]);
      s->zoom_props->nr_zoom_cells = s->width[0] / s->zoom_props->width[0];
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

///**
// * @brief Compute pair tasks between natural TL cells and zoom TL cells.
// *
// * Loop over the zoom TL cells, then loop over the "neighbours".
// *
// * Apply standard MM/max_distance check to see if we need a pair gravity
// * task between the two cells.
// */
//void engine_make_self_gravity_tasks_mapper_between_toplevels(void *map_data,
//                                                             int num_elements,
//                                                             void *extra_data) {
//#ifdef WITH_ZOOM_REGION
//  struct engine *e = (struct engine *)extra_data;
//  struct space *s = e->s;
//  struct scheduler *sched = &e->sched;
//  const int nodeID = e->nodeID;
//  const int periodic = s->periodic;
//  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
//  struct cell *cells = s->cells_top;
//  const double max_distance = e->mesh->r_cut_max;
//  const double max_distance2 = max_distance * max_distance;
//  const int nr_cells = s->nr_cells;
//
//  /* Loop through the elements, which are just byte offsets from NULL. */
//  for (int ind = 0; ind < num_elements; ind++) {
//
//    /* Get the cell index. */
//    const int cid = (size_t)(map_data) + ind;
//
//    /* Only compute for zoom TL cells. */
//    if (cid < s->zoom_props->tl_cell_offset) continue;
//
//    /* Get the cell */
//    struct cell *ci = &cells[cid];
//
//    /* Skip cells without gravity particles */
//    if (ci->grav.count == 0) continue;
//
//    for (int cjd = 0; cjd < nr_cells; cjd++) {
//
//      /* Get the cell */
//      struct cell *cj = &cells[cjd];
//
//      /* Only considering "neighbours". */
//      if (cj->tl_cell_type != tl_cell_neighbour) continue;
//
//      /* Avoid duplicates, empty cells and completely foreign pairs */
//      if (cj->grav.count == 0 || (ci->nodeID != nodeID && cj->nodeID != nodeID))
//        continue;
//
//      /* Recover the multipole information */
//      const struct gravity_tensors *multi_i = ci->grav.multipole;
//      const struct gravity_tensors *multi_j = cj->grav.multipole;
//
//      if (multi_i == NULL && ci->nodeID != nodeID)
//        error("Multipole of ci was not exchanged properly via the proxies");
//      if (multi_j == NULL && cj->nodeID != nodeID)
//        error("Multipole of cj was not exchanged properly via the proxies");
//
//      /* Minimal distance between any pair of particles */
//      const double min_radius2 =
//          cell_min_dist2_diff_size(ci, cj, periodic, dim);
//
//      /* Are we beyond the distance where the truncated forces are 0 ?*/
//      if (periodic && min_radius2 > max_distance2) continue;
//
//      /* Are the cells too close for a MM interaction ? */
//      if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
//                                    /*is_tree_walk=*/0)) {
//
//        /* Ok, we need to add a direct pair calculation */
//        scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0, ci,
//                          cj);
//#ifdef SWIFT_DEBUG_CHECKS
//#ifdef WITH_MPI
//        /* Let's cross-check that we had a proxy for that cell */
//        if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {
//
//          /* Find the proxy for this node */
//          const int proxy_id = e->proxy_ind[cj->nodeID];
//          if (proxy_id < 0)
//            error("No proxy exists for that foreign node %d!", cj->nodeID);
//
//          const struct proxy *p = &e->proxies[proxy_id];
//
//          /* Check whether the cell exists in the proxy */
//          int n = 0;
//          for (; n < p->nr_cells_in; n++)
//            if (p->cells_in[n] == cj) {
//              break;
//            }
//          if (n == p->nr_cells_in)
//            error(
//                "Cell %d not found in the proxy but trying to construct "
//                "grav task!",
//                cjd);
//        } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {
//
//          /* Find the proxy for this node */
//          const int proxy_id = e->proxy_ind[ci->nodeID];
//          if (proxy_id < 0)
//            error("No proxy exists for that foreign node %d!", ci->nodeID);
//
//          const struct proxy *p = &e->proxies[proxy_id];
//
//          /* Check whether the cell exists in the proxy */
//          int n = 0;
//          for (; n < p->nr_cells_in; n++)
//            if (p->cells_in[n] == ci) {
//              break;
//            }
//          if (n == p->nr_cells_in)
//            error(
//                "Cell %d not found in the proxy but trying to construct "
//                "grav task!",
//                cid);
//        }
//#endif /* WITH_MPI */
//#endif /* SWIFT_DEBUG_CHECKS */
//
//      }
//    }
//  }
//#endif
//}

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
#ifdef WITH_ZOOM_REGION
/**
 * @brief Create and fill the proxies including the zoom region.
 *
 * This replaces the function in engine_proxy when running with a zoom region.
 *
 * @param e The #engine.
 */
void engine_makeproxies_with_zoom_region(struct engine *e) {

	if (!e->s->with_zoom_region) engine_makeproxies(e);

#ifdef WITH_MPI
    /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the zoom domain */
  const int zoom_cell_offset = s->zoom_props->tl_cell_offset;
  int zoom_cell_flag = 0;

  /* Some info about the domain */
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Distance between centre of the cell and corners for natural top level cells */
  double r_diag2 = cells[0].width[0] * cells[0].width[0] +
  		             cells[0].width[1] * cells[0].width[1] +
  		             cells[0].width[2] * cells[0].width[2];
  double zoom_r_diag2 = cells[zoom_cell_offset].width[0] * cells[zoom_cell_offset].width[0] +
                        cells[zoom_cell_offset].width[1] * cells[zoom_cell_offset].width[1] +
                        cells[zoom_cell_offset].width[2] * cells[zoom_cell_offset].width[2];

  const double r_diag = 0.5 * sqrt(r_diag2);
  const double zoom_r_diag = 0.5 * sqrt(zoom_r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double natural_r_max = 2 * r_diag;
  const double zoom_r_max = 2 * zoom_r_diag;
  double r_max = natural_r_max; /* initially the loop needs the natural distance */
  const double distance = 2. * r_max * theta_crit_inv;
  const double zoom_distance = 2. * zoom_r_max * theta_crit_inv;

  /* Compute how many cells away we need to walk */
  const int delta_cells = (int)(distance / cells[0].dmin) + 1;
  const int natural_delta_cells = delta_cells;
  const int zoom_delta_cells = (int)(zoom_distance / cells[zoom_cell_offset].dmin) + 1;

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;
  int natural_delta_m = natural_delta_cells;
  int natural_delta_p = natural_delta_cells;
  int zoom_delta_m = zoom_delta_cells;
  int zoom_delta_p = zoom_delta_cells;

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

  /* Special case where every zoom cell is in range of every other one */
  if (zoom_delta_cells >= cdim[0] / 2) {
  	if (cdim[0] % 2 == 0) {
  		zoom_delta_m = cdim[0] / 2;
  		zoom_delta_p = cdim[0] / 2 - 1;
  	} else {
  		zoom_delta_m = cdim[0] / 2;
  		zoom_delta_p = cdim[0] / 2;
  	}
  }
  /* Let's be verbose about this choice */
  if (e->verbose) {
    message(
        "Looking for proxies up to %d top-level cells away (delta_m=%d "
				"delta_p=%d)",
				delta_cells, delta_m, delta_p);
    message(
        "Looking for zoom proxies up to %d top-level zoom cells away (delta_m=%d "
				"delta_p=%d)",
				(zoom_delta_m + zoom_delta_p) / 2, zoom_delta_m, zoom_delta_p);
  }

  /* Loop over cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

  	/* Get the cell including any possible zoom offset */
  	int cid_with_offset = cid;

  	zoom_cell_flag = 0;

  	/* Integer indices of the cell in the natural parent */
  	int natural_i = cid / (cdim[1] * cdim[2]);
  	int natural_j = (cid / cdim[2]) % cdim[1];
  	int natural_k = cid % cdim[2];

  	if (cid >= zoom_cell_offset) {
  		zoom_cell_flag = 1;
  		cid_with_offset -= zoom_cell_offset;

  		/* Integer indices of the cell in the natural parent */
  		int natural_tl_cid = cells[cid].parent_tl_cid;
  		natural_i = natural_tl_cid / (cdim[1] * cdim[2]);
  		natural_j = (natural_tl_cid / cdim[2]) % cdim[1];
  		natural_k = natural_tl_cid % cdim[2];

  		/* Overwrite delta_m and delta_p with the number of zoom cells*/
  		delta_m = zoom_delta_m;
  		delta_p = zoom_delta_p;

  		/* Set upper search distance to zoom */
  		r_max = zoom_r_max;

  	}

  	/* Integer indices of the cell in the top-level grid */
  	const int i = cid_with_offset / (cdim[1] * cdim[2]);
  	const int j = (cid_with_offset / cdim[2]) % cdim[1];
  	const int k = cid_with_offset % cdim[2];

  	/* Loop over all its neighbours in range. */
  	for (int ii = -delta_m; ii <= delta_p; ii++) {
  		int iii = i + ii;
  		if ((!periodic || zoom_cell_flag) && (iii < 0 || iii >= cdim[0])) continue;
  		iii = (iii + cdim[0]) % cdim[0];
  		for (int jj = -delta_m; jj <= delta_p; jj++) {
  			int jjj = j + jj;
  			if ((!periodic || zoom_cell_flag) && (jjj < 0 || jjj >= cdim[1])) continue;
  			jjj = (jjj + cdim[1]) % cdim[1];
  			for (int kk = -delta_m; kk <= delta_p; kk++) {
  				int kkk = k + kk;
  				if ((!periodic || zoom_cell_flag) && (kkk < 0 || kkk >= cdim[2])) continue;
  				kkk = (kkk + cdim[2]) % cdim[2];

  				/* Get the cell ID. */
  				int cjd = cell_getid(cdim, iii, jjj, kkk);
  				if (zoom_cell_flag) cjd += zoom_cell_offset;

          /* Early abort  */
          if (cid >= cjd) continue;

          /* Early abort (both same node) */
          if (cells[cid].nodeID == nodeID && cells[cjd].nodeID == nodeID)
            continue;

          /* Early abort (both foreign node) */
          if (cells[cid].nodeID != nodeID && cells[cjd].nodeID != nodeID)
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
              const double min_dist_CoM2 = cell_min_dist2(
                      &cells[cid], &cells[cjd], periodic, dim);

              /* Are we beyond the distance where the truncated forces are 0
               * but not too far such that M2L can be used? */
              if (periodic && !zoom_cell_flag) {

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
          if (cells[cid].nodeID == nodeID && cells[cjd].nodeID != nodeID) {

            /* Do we already have a relationship with this node? */
            int proxy_id = e->proxy_ind[cells[cjd].nodeID];
            if (proxy_id < 0) {
              if (e->nr_proxies == engine_maxproxies)
                error("Maximum number of proxies exceeded.");

              /* Ok, start a new proxy for this pair of nodes */
              proxy_init(&proxies[e->nr_proxies], e->nodeID,
                         cells[cjd].nodeID);

              /* Store the information */
              e->proxy_ind[cells[cjd].nodeID] = e->nr_proxies;
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
            cells[cid].mpi.sendto |= (1ULL << proxy_id);
          }

          /* Same for the symmetric case? */
          if (cells[cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

            /* Do we already have a relationship with this node? */
            int proxy_id = e->proxy_ind[cells[cid].nodeID];
            if (proxy_id < 0) {
              if (e->nr_proxies == engine_maxproxies)
                error("Maximum number of proxies exceeded.");

              /* Ok, start a new proxy for this pair of nodes */
              proxy_init(&proxies[e->nr_proxies], e->nodeID,
                         cells[cid].nodeID);

              /* Store the information */
              e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
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
            cells[cjd].mpi.sendto |= (1ULL << proxy_id);
          }

          /* For natural top level cell neighbours in the zoom
           * region we need to include the nested zoom cells
           * NOTE: we will never have hydro tasks between grid levels */
          if (!zoom_cell_flag && cells[cjd].tl_cell_type == void_tl_cell) {

          	int start_i = cells[cjd].start_i;
          	int start_j = cells[cjd].start_j;
          	int start_k = cells[cjd].start_k;
            int end_i = cells[cjd].end_i + 1;
            int end_j = cells[cjd].end_j + 1;
            int end_k = cells[cjd].end_k + 1;

            for (int iiii = start_i; iiii < end_i; iiii++) {
              for (int jjjj = start_j; jjjj < end_j; jjjj++) {
                for (int kkkk = start_k; kkkk < end_k; kkkk++) {

                  /* Zoom level neighbour */
                  int zoom_cjd = cell_getid(cdim, iiii, jjjj, kkkk) + zoom_cell_offset;

                  /* Early abort (both same node) */
			            if (cells[cid].nodeID == nodeID && cells[zoom_cjd].nodeID == nodeID)
			              continue;

			            /* Early abort (both foreign node) */
			            if (cells[cid].nodeID != nodeID && cells[zoom_cjd].nodeID != nodeID)
			              continue;

									/* Ensure we are still in the natural neighbour */
									if (zoom_cjd != cjd) continue;

                  proxy_type = 0;

                  /* In the gravity case, check distances using the MAC.
                   * NOTE: The direct neighbour check is done on the natural
                   * grid level (i and iii etc), i.e. all zoom cells within
                   * a natural direct neighbour are direct neighbours */
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
                      const double min_dist_CoM2 = cell_min_dist2(
                              &cells[cid], &cells[zoom_cjd], periodic, dim);

                      /* Are we beyond the distance where the truncated forces are 0
                       * but not too far such that M2L can be used? */
                      if (periodic) {

                        if ((min_dist_CoM2 < max_mesh_dist2) &&
                            !(4. * natural_r_max * natural_r_max <
                              theta_crit * theta_crit * min_dist_CoM2))
                          proxy_type |= (int)proxy_cell_type_gravity;

                      } else {

                        if (!(4. * natural_r_max * natural_r_max <
                              theta_crit * theta_crit * min_dist_CoM2)) {
                          proxy_type |= (int)proxy_cell_type_gravity;
                        }
                      }
                    }
                  }

                  /* Abort if not in range at all */
                  if (proxy_type == proxy_cell_type_none) continue;

                  /* Add to proxies? */
                  if (cells[cid].nodeID == nodeID && cells[zoom_cjd].nodeID != nodeID) {

                    /* Do we already have a relationship with this node? */
                    int proxy_id = e->proxy_ind[cells[zoom_cjd].nodeID];
                    if (proxy_id < 0) {
                      if (e->nr_proxies == engine_maxproxies)
                        error("Maximum number of proxies exceeded.");

                      /* Ok, start a new proxy for this pair of nodes */
                      proxy_init(&proxies[e->nr_proxies], e->nodeID,
                                 cells[zoom_cjd].nodeID);

                      /* Store the information */
                      e->proxy_ind[cells[zoom_cjd].nodeID] = e->nr_proxies;
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
                    proxy_addcell_in(&proxies[proxy_id], &cells[zoom_cjd], proxy_type);
                    proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                    /* Store info about where to send the cell */
                    cells[cid].mpi.sendto |= (1ULL << proxy_id);
                  }

                  /* Same for the symmetric case? */
                  if (cells[zoom_cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

                    /* Do we already have a relationship with this node? */
                    int proxy_id = e->proxy_ind[cells[cid].nodeID];
                    if (proxy_id < 0) {
                      if (e->nr_proxies == engine_maxproxies)
                        error("Maximum number of proxies exceeded.");

                      /* Ok, start a new proxy for this pair of nodes */
                      proxy_init(&proxies[e->nr_proxies], e->nodeID,
                                 cells[cid].nodeID);

                      /* Store the information */
                      e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
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
                    proxy_addcell_out(&proxies[proxy_id], &cells[zoom_cjd], proxy_type);

                    /* Store info about where to send the cell */
                    cells[zoom_cjd].mpi.sendto |= (1ULL << proxy_id);
                  }
                }
              }
            }
          }
        }
      }
    }

    /* For the zoom cells we need to find all natural neighbours */
    if (zoom_cell_flag) {

      /* Loop over all its neighbours in range. */
      for (int ii = -natural_delta_m; ii <= natural_delta_p; ii++) {
        int iii = natural_i + ii;
        if ((!periodic) && (iii < 0 || iii >= cdim[0])) continue;
        iii = (iii + cdim[0]) % cdim[0];
        for (int jj = -natural_delta_m; jj <= natural_delta_p; jj++) {
          int jjj = natural_j + jj;
          if ((!periodic) && (jjj < 0 || jjj >= cdim[1])) continue;
          jjj = (jjj + cdim[1]) % cdim[1];
          for (int kk = -natural_delta_m; kk <= natural_delta_p; kk++) {
            int kkk = natural_k + kk;
            if ((!periodic) && (kkk < 0 || kkk >= cdim[2])) continue;
            kkk = (kkk + cdim[2]) % cdim[2];

            /* Get the cell ID. */
            int nat_cjd = cell_getid(cdim, iii, jjj, kkk);

            /* Early abort (both same node) */
            if (cells[cid].nodeID == nodeID && cells[nat_cjd].nodeID == nodeID)
              continue;

            /* Early abort (both foreign node) */
            if (cells[cid].nodeID != nodeID && cells[nat_cjd].nodeID != nodeID)
              continue;

            int proxy_type = 0;

            /* In the hydro case, only care about direct neighbours */
            if (with_hydro) {

              // MATTHIEU: to do: Write a better expression for the
              // non-periodic case.

              /* This is super-ugly but checks for direct neighbours */
              /* with periodic BC */
              if (((abs(natural_i - iii) <= 1 || abs(natural_i - iii - cdim[0]) <= 1 ||
                    abs(natural_i - iii + cdim[0]) <= 1) &&
                   (abs(natural_j - jjj) <= 1 || abs(natural_j - jjj - cdim[1]) <= 1 ||
                    abs(natural_j - jjj + cdim[1]) <= 1) &&
                   (abs(natural_k - kkk) <= 1 || abs(natural_k - kkk - cdim[2]) <= 1 ||
                    abs(natural_k - kkk + cdim[2]) <= 1)))
                proxy_type |= (int)proxy_cell_type_hydro;
            }

            /* In the gravity case, check distances using the MAC. */
            if (with_gravity) {

              /* First just add the direct neighbours. Then look for
                 some further out if the opening angle demands it */

              /* This is super-ugly but checks for direct neighbours */
              /* with periodic BC */
              if (((abs(natural_i - iii) <= 1 || abs(natural_i - iii - cdim[0]) <= 1 ||
                    abs(natural_i - iii + cdim[0]) <= 1) &&
                   (abs(natural_j - jjj) <= 1 || abs(natural_j - jjj - cdim[1]) <= 1 ||
                    abs(natural_j - jjj + cdim[1]) <= 1) &&
                   (abs(natural_k - kkk) <= 1 || abs(natural_k - kkk - cdim[2]) <= 1 ||
                    abs(natural_k - kkk + cdim[2]) <= 1))) {

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
                const double min_dist_CoM2 = cell_min_dist2(
                        &cells[cid], &cells[nat_cjd], periodic, dim);

                /* Are we beyond the distance where the truncated forces are 0
                 * but not too far such that M2L can be used? */
                if (periodic && !zoom_cell_flag) {

                  if ((min_dist_CoM2 < max_mesh_dist2) &&
                      !(4. * natural_r_max * natural_r_max <
                        theta_crit * theta_crit * min_dist_CoM2))
                    proxy_type |= (int)proxy_cell_type_gravity;

                } else {

                  if (!(4. * natural_r_max * natural_r_max <
                        theta_crit * theta_crit * min_dist_CoM2)) {
                    proxy_type |= (int)proxy_cell_type_gravity;
                  }
                }
              }
            }

            /* Abort if not in range at all */
            if (proxy_type == proxy_cell_type_none) continue;

            /* Add to proxies? */
            if (cells[cid].nodeID == nodeID && cells[nat_cjd].nodeID != nodeID) {

              /* Do we already have a relationship with this node? */
              int proxy_id = e->proxy_ind[cells[nat_cjd].nodeID];
              if (proxy_id < 0) {
                if (e->nr_proxies == engine_maxproxies)
                  error("Maximum number of proxies exceeded.");

                /* Ok, start a new proxy for this pair of nodes */
                proxy_init(&proxies[e->nr_proxies], e->nodeID,
                           cells[nat_cjd].nodeID);

                /* Store the information */
                e->proxy_ind[cells[nat_cjd].nodeID] = e->nr_proxies;
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
              proxy_addcell_in(&proxies[proxy_id], &cells[nat_cjd], proxy_type);
              proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

              /* Store info about where to send the cell */
              cells[cid].mpi.sendto |= (1ULL << proxy_id);
            }

            /* Same for the symmetric case? */
            if (cells[nat_cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

              /* Do we already have a relationship with this node? */
              int proxy_id = e->proxy_ind[cells[cid].nodeID];
              if (proxy_id < 0) {
                if (e->nr_proxies == engine_maxproxies)
                  error("Maximum number of proxies exceeded.");

                /* Ok, start a new proxy for this pair of nodes */
                proxy_init(&proxies[e->nr_proxies], e->nodeID,
                           cells[cid].nodeID);

                /* Store the information */
                e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
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
              proxy_addcell_out(&proxies[proxy_id], &cells[nat_cjd], proxy_type);

              /* Store info about where to send the cell */
              cells[nat_cjd].mpi.sendto |= (1ULL << proxy_id);
            }
          }
        }
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

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions.
 *
 * This replaces the function in engine_maketasks when running with a zoom region.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.a
 * @param extra_data The #engine.

 */
void engine_make_self_gravity_tasks_mapper_with_zoom(void *map_data, int num_elements,
																										 void *extra_data) {

	/* Useful local information */
	struct engine *e = (struct engine *)extra_data;
	struct space *s = e->s;
	struct scheduler *sched = &e->sched;
	const int nodeID = e->nodeID;

	/* Handle on the cells and proxies */
	struct cell *cells = s->cells_top;

	/* Some info about the zoom domain */
	const int zoom_cell_offset = s->zoom_props->tl_cell_offset;

	/* Some info about the domain */
	const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
	const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
	int periodic = s->periodic;

	message("periodic=%d zoom_cell_offset=%d", periodic, zoom_cell_offset)

	/* Get some info about the physics */
	const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
	const double max_mesh_dist = e->mesh->r_cut_max;
	const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

	/* Maximal distance from shifted CoM to any corner */
	const double distance = 2. * cells[0].width[0] * theta_crit_inv;
	const double zoom_distance = 2. * cells[zoom_cell_offset].width[0] * theta_crit_inv;

	/* Compute how many cells away we need to walk */
	const int delta_cells = (int)(distance / cells[0].dmin) + 1;
	const int natural_delta_cells = delta_cells;
	const int zoom_delta_cells = (int)(zoom_distance / cells[zoom_cell_offset].dmin) + 1;

	/* Turn this into upper and lower bounds for loops */
	int delta_m = delta_cells;
	int delta_p = delta_cells;
	int natural_delta_m = natural_delta_cells;
	int natural_delta_p = natural_delta_cells;
	int zoom_delta_m = zoom_delta_cells;
	int zoom_delta_p = zoom_delta_cells;

	/* Special case where every cell is in range of every other one */
	if (delta_cells >= cdim[0] / 2) {
		if (cdim[0] % 2 == 0) {
			delta_m = cdim[0] / 2;
			delta_p = cdim[0] / 2 - 1;
			natural_delta_m = cdim[0] / 2;
			natural_delta_p = cdim[0] / 2 - 1;
		} else {
			delta_m = cdim[0] / 2;
			delta_p = cdim[0] / 2;
			natural_delta_m = cdim[0] / 2;
			natural_delta_p = cdim[0] / 2 - 1;
		}
	}
	if (zoom_delta_cells >= cdim[0] / 2) {
		if (cdim[0] % 2 == 0) {
			zoom_delta_m = cdim[0] / 2;
			zoom_delta_p = cdim[0] / 2 - 1;
		} else {
			zoom_delta_m = cdim[0] / 2;
			zoom_delta_p = cdim[0] / 2;
		}
	}

	/* Loop through the elements, which are just byte offsets from NULL. */
	for (int ind = 0; ind < num_elements; ind++) {

		/* Get the cell index. */
		const int cid = (size_t)(map_data) + ind;
		int cid_with_offset = cid;

		/* Get the cell */
		struct cell *ci = &cells[cid];

		/* Integer grid indices of the cell in the natural cells */
		int natural_i = (ci->loc[0] + (s->zoom_props->width[0] / 2)) * s->iwidth[0];
		int natural_j = (ci->loc[1] + (s->zoom_props->width[1] / 2)) * s->iwidth[1];
		int natural_k = (ci->loc[2] + (s->zoom_props->width[2] / 2)) * s->iwidth[2];

		if (cid >= zoom_cell_offset) {

			/* Remove the offset for zoom cell ids */
			cid_with_offset -= zoom_cell_offset;

			/* Overwrite delta_m and delta_p with the number of zoom cells */
			delta_m = zoom_delta_m;
			delta_p = zoom_delta_p;

		} else {
			/* Overwrite delta_m and delta_p with the number of natural cells */
			delta_m = natural_delta_m;
			delta_p = natural_delta_p;
		}

		/* Integer indices of the cell in the top-level grid */
		const int i = cid_with_offset / (cdim[1] * cdim[2]);
		const int j = (cid_with_offset / cdim[2]) % cdim[1];
		const int k = cid_with_offset % cdim[2];

		/* Skip cells without gravity particles */
		if (ci->grav.count == 0) continue;

		/* If the cell is local build a self-interaction */
		if (ci->nodeID == nodeID) {
			scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
												NULL);
		}

		/* Loop over all its neighbours in range. */
		for (int ii = -delta_m; ii <= delta_p; ii++) {
			int iii = i + ii;
			if ((!periodic || cid >= zoom_cell_offset) && (iii < 0 || iii >= cdim[0])) continue;
			iii = (iii + cdim[0]) % cdim[0];
			for (int jj = -delta_m; jj <= delta_p; jj++) {
				int jjj = j + jj;
				if ((!periodic || cid >= zoom_cell_offset) && (jjj < 0 || jjj >= cdim[1])) continue;
				jjj = (jjj + cdim[1]) % cdim[1];
				for (int kk = -delta_m; kk <= delta_p; kk++) {
					int kkk = k + kk;
					if ((!periodic || cid >= zoom_cell_offset) && (kkk < 0 || kkk >= cdim[2])) continue;
					kkk = (kkk + cdim[2]) % cdim[2];

					/* Get the cell ID. */
					int cjd = cell_getid(cdim, iii, jjj, kkk);
					message("cid=%d cid_with_offset=%d cjd=%d cjd_with_offset=%d ii,jj,kk=%d %d %d i + ii,j + jj,k + kk=%d %d %d iii,jjj,kkk=%d %d %d", cid, cid_with_offset, cjd, cjd + zoom_cell_offset, ii, jj, kk, i + ii, j + jj, k + kk, iii, jjj, kkk);
					if (cid >= zoom_cell_offset) {
						cjd += zoom_cell_offset;
					}

					struct cell *cj = &cells[cjd];

					/* Avoid duplicates, void cells, empty cells and completely foreign pairs */
					if (cid >= cjd || cj->grav.count == 0 ||
					(ci->nodeID != nodeID && cj->nodeID != nodeID))
						continue;

					/* Recover the multipole information */
					const struct gravity_tensors *multi_i = ci->grav.multipole;
					const struct gravity_tensors *multi_j = cj->grav.multipole;

					if (multi_i == NULL && ci->nodeID != nodeID)
						error("Multipole of ci was not exchanged properly via the proxies");
					if (multi_j == NULL && cj->nodeID != nodeID)
						error("Multipole of cj was not exchanged properly via the proxies");

					/* Minimal distance between any pair of particles */
					const double min_radius2 = cell_min_dist2_same_size(ci, cj, periodic, dim);

					/* Are we beyond the distance where the truncated forces are 0 ?*/
					if (periodic && (cid < zoom_cell_offset) && min_radius2 > max_mesh_dist2) continue;

					/* Are the cells too close for a MM interaction ? */
					if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
																		/*is_tree_walk=*/0)) {

						/* Ok, we need to add a direct pair calculation */
						scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
															ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
						/* Ensure both cells are in the same level */
						if (!((ci->tl_cell_type <= 2 && cj->tl_cell_type <= 2) ||
						(ci->tl_cell_type == zoom_tl_cell && cj->tl_cell_type == zoom_tl_cell))) {
							error("Cell %d and cell %d are not the same cell type! (ci=%d, cj=%d)",
										cid, cjd, ci->tl_cell_type, cj->tl_cell_type);
						}
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

//					/* For natural top level cell neighbours in the zoom
//           * region we need to include the nested zoom cells */
//					if (cid < zoom_cell_offset && cj->tl_cell_type == void_tl_cell) {
//
//						int start_i = cj->start_i;
//						int start_j = cj->start_j;
//						int start_k = cj->start_k;
//						int end_i = cells[cjd].end_i + 1;
//						int end_j = cells[cjd].end_j + 1;
//						int end_k = cells[cjd].end_k + 1;
//
//						for (int iiii = start_i; iiii < end_i; iiii++) {
//							for (int jjjj = start_j; jjjj < end_j; jjjj++) {
//								for (int kkkk = start_k; kkkk < end_k; kkkk++) {
//
//									/* Zoom level neighbour */
//									const int zoom_cjd = cell_getid(cdim, iiii, jjjj, kkkk) + zoom_cell_offset;
//
//									struct cell *cj_zoom = &cells[zoom_cjd];
//
//									/* Ensure we are still in the natural neighbour */
//									if (cj_zoom->parent_tl_cid != cjd) continue;
//
//									/* Avoid duplicates, empty cells and completely foreign pairs */
//									if (cj_zoom->grav.count == 0 ||
//									(ci->nodeID != nodeID && cj_zoom->nodeID != nodeID))
//										continue;
//
//									/* Recover the multipole information */
//									const struct gravity_tensors *multi_j_zoom = cj_zoom->grav.multipole;
//
//									if (multi_i == NULL && ci->nodeID != nodeID)
//										error("Multipole of ci was not exchanged properly via the proxies");
//									if (multi_j_zoom == NULL && cj_zoom->nodeID != nodeID)
//										error("Multipole of cj was not exchanged properly via the proxies");
//
//									/* Are the cells too close for a MM interaction ? */
//									if (!cell_can_use_pair_mm(ci, cj_zoom, e, s, /*use_rebuild_data=*/1,
//																						/*is_tree_walk=*/0)) {
//
//										/* Ok, we need to add a direct pair calculation */
//										scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
//																			ci, cj_zoom);
//
//#ifdef SWIFT_DEBUG_CHECKS
//#ifdef WITH_MPI
//
//										/* Let's cross-check that we had a proxy for that cell */
//										if (ci->nodeID == nodeID && cj_zoom->nodeID != engine_rank) {
//
//											/* Find the proxy for this node */
//											const int proxy_id = e->proxy_ind[cj_zoom->nodeID];
//											if (proxy_id < 0)
//												error("No proxy exists for that foreign node %d!", cj_zoom->nodeID);
//
//											const struct proxy *p = &e->proxies[proxy_id];
//
//											/* Check whether the cell exists in the proxy */
//											int n = 0;
//											for (; n < p->nr_cells_in; n++)
//												if (p->cells_in[n] == cj_zoom) {
//													break;
//												}
//											if (n == p->nr_cells_in)
//												error(
//														"Cell %d not found in the proxy but trying to construct "
//														"grav task!",
//														zoom_cjd);
//										} else if (cj_zoom->nodeID == nodeID && ci->nodeID != engine_rank) {
//
//											/* Find the proxy for this node */
//											const int proxy_id = e->proxy_ind[ci->nodeID];
//											if (proxy_id < 0)
//												error("No proxy exists for that foreign node %d!", ci->nodeID);
//
//											const struct proxy *p = &e->proxies[proxy_id];
//
//											/* Check whether the cell exists in the proxy */
//											int n = 0;
//											for (; n < p->nr_cells_in; n++)
//												if (p->cells_in[n] == ci) {
//													break;
//												}
//											if (n == p->nr_cells_in)
//												error(
//														"Cell %d not found in the proxy but trying to construct "
//														"grav task!",
//														cid);
//										}
//#endif /* WITH_MPI */
//#endif /* SWIFT_DEBUG_CHECKS */
//
//									}
//								}
//							}
//						}
//					}
				}
			}
		}

		/* For the zoom cells we need to find all natural cell neighbours */
		if (ci->tl_cell_type == zoom_tl_cell) {

			/* Loop over all its neighbours in range. */
			for (int ii = -natural_delta_m; ii <= natural_delta_p; ii++) {
				int nat_iii = natural_i + ii;
				if ((!periodic) && (nat_iii < 0 || nat_iii >= cdim[0])) continue;
				nat_iii = (nat_iii + cdim[0]) % cdim[0];
				for (int jj = -natural_delta_m; jj <= natural_delta_p; jj++) {
					int nat_jjj = natural_j + jj;
					if ((!periodic) && (nat_jjj < 0 || nat_jjj >= cdim[1])) continue;
					nat_jjj = (nat_jjj + cdim[1]) % cdim[1];
					for (int kk = -natural_delta_m; kk <= natural_delta_p; kk++) {
						int nat_kkk = natural_k + kk;
						if ((!periodic) && (nat_kkk < 0 || nat_kkk >= cdim[2])) continue;
						nat_kkk = (nat_kkk + cdim[2]) % cdim[2];

						/* Get the cell ID. */
						int nat_cjd = cell_getid(cdim, nat_iii, nat_jjj, nat_kkk);
						message("cid=%d cid_with_offset=%d nat_cjd=%d ii,jj,kk=%d %d %d natural_i + ii,natural_j + jj,natural_k + kk=%d %d %d nat_iii,nat_jjj,nat_kkk=%d %d %d", cid, cid_with_offset, nat_cjd, ii, jj, kk, natural_i + ii, natural_j + jj, natural_k + kk, nat_iii, nat_jjj, nat_kkk);

						struct cell *nat_cj = &cells[nat_cjd];

						/* empty cells and completely foreign pairs */
						if (nat_cj->grav.count == 0 ||
						(ci->nodeID != nodeID && nat_cj->nodeID != nodeID))
							continue;

						/* Recover the multipole information */
						const struct gravity_tensors *multi_i = ci->grav.multipole;
						const struct gravity_tensors *nat_multi_j = nat_cj->grav.multipole;

						if (multi_i == NULL && ci->nodeID != nodeID)
							error("Multipole of ci was not exchanged properly via the proxies");
						if (nat_multi_j == NULL && nat_cj->nodeID != nodeID)
							error("Multipole of nat_cj was not exchanged properly via the proxies");

						/* Minimal distance between cells */
						const double nat_min_radius2 = cell_min_dist2_diff_size(ci, nat_cj, periodic, dim);

						/* Are we beyond the distance where the truncated forces are 0 ?*/
						if (periodic && nat_min_radius2 > max_mesh_dist2) continue;

						/* Are the cells too close for a MM interaction ? */
						if (!cell_can_use_pair_mm(ci, nat_cj, e, s, /*use_rebuild_data=*/1,
																			/*is_tree_walk=*/0)) {

							/* Ok, we need to add a direct pair calculation */
							scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
																ci, nat_cj);

#ifdef SWIFT_DEBUG_CHECKS
						/* Ensure both cells are in different levels */
						if ((ci->tl_cell_type <= 2 && nat_cj->tl_cell_type <= 2) ||
						(ci->tl_cell_type == zoom_tl_cell && nat_cj->tl_cell_type == zoom_tl_cell)) {
							error("Cell %d and cell %d are the same cell type! (ci=%d, cj=%d)",
										cid, nat_cjd, ci->tl_cell_type, nat_cj->tl_cell_type);
						}
#ifdef WITH_MPI

							/* Let's cross-check that we had a proxy for that cell */
							if (ci->nodeID == nodeID && nat_cj->nodeID != engine_rank) {

								/* Find the proxy for this node */
								const int proxy_id = e->proxy_ind[nat_cj->nodeID];
								if (proxy_id < 0)
									error("No proxy exists for that foreign node %d!", nat_cj->nodeID);

								const struct proxy *p = &e->proxies[proxy_id];

								/* Check whether the cell exists in the proxy */
								int n = 0;
								for (; n < p->nr_cells_in; n++)
									if (p->cells_in[n] == nat_cj) {
										break;
									}
								if (n == p->nr_cells_in)
									error(
											"Cell %d not found in the proxy but trying to construct "
											"grav task!",
											nat_cjd);
							} else if (nat_cj->nodeID == nodeID && ci->nodeID != engine_rank) {

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
			}
		}
	}
}
#endif /* WITH_ZOOM_REGION */


