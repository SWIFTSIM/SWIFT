/* Config parameters. */
#include "../config.h"

#include <float.h>

#include "engine.h"
#include "space.h"
#include "cell.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

void zoom_region_init(struct swift_params *params, struct space *s) {

  /* Initialise the zoom region. */
  s->with_zoom_region =
      parser_get_opt_param_int(params, "ZoomRegion:enable", 0);

  if (s->with_zoom_region) {
    s->zoom_props = (struct zoom_region_properties *)malloc(
        sizeof(struct zoom_region_properties));
    if (s->zoom_props == NULL)
      error("Error allocating memory for the zoom parameters.");

    s->zoom_props->mask_parttype =
        parser_get_param_int(params, "ZoomRegion:mask_parttype");
    s->zoom_props->max_size =
        parser_get_opt_param_double(params, "ZoomRegion:max_size", FLT_MAX);
    s->zoom_props->boost_factor =
        parser_get_opt_param_double(params, "ZoomRegion:boost_factor", 1.0);
  }
}

void construct_zoom_region(struct space *s, int verbose) {

  double new_zoom_boundary[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
  const size_t nr_gparts = s->nr_gparts;
  double mtot = 0.0;
  double com[3] = {0.0, 0.0, 0.0};

  for (size_t k = 0; k < nr_gparts; k++) {
    if (s->gparts[k].type != s->zoom_props->mask_parttype) continue;

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
  /* Make sure all nodes agree on the zoom region. */
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
  const double imass = 1.0 / mtot;
  com[0] *= imass;
  com[1] *= imass;
  com[2] *= imass;

  for (int k = 0; k < 3; k++) s->zoom_props->loc[k] = com[k];

  s->zoom_props->dim[0] = (new_zoom_boundary[1] - new_zoom_boundary[0]) *
                          s->zoom_props->boost_factor;
  s->zoom_props->dim[1] = (new_zoom_boundary[3] - new_zoom_boundary[2]) *
                          s->zoom_props->boost_factor;
  s->zoom_props->dim[2] = (new_zoom_boundary[5] - new_zoom_boundary[4]) *
                          s->zoom_props->boost_factor;

  for (int k = 0; k < 3; k++)
    if (s->zoom_props->dim[k] > s->zoom_props->max_size)
      s->zoom_props->dim[k] = s->zoom_props->max_size;

  if (verbose)
    message("zoom loc: [%f %f %f] zoom dim: [%f %f %f]", s->zoom_props->loc[0],
            s->zoom_props->loc[1], s->zoom_props->loc[2],
            s->zoom_props->dim[0], s->zoom_props->dim[1],
            s->zoom_props->dim[2]);
}

void check_zoom_region(const struct space *s, const int verbose) {
  if (verbose) {
    message("top level zoom cell dimensions [ %i %i %i ].", s->zoom_props->cdim[0],
			s->zoom_props->cdim[1], s->zoom_props->cdim[2]);
    message(
        "zoom region centered on [ %.3f %.3f %.3f ] with dimensions [ %.3f "
        "%.3f %.3f ]",
        s->zoom_props->loc[0], s->zoom_props->loc[1], s->zoom_props->loc[2],
        s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  }
}

int cell_getid_zoom(const int cdim[3], const double x, const double y, const double z,
	const struct zoom_region_properties *zoom_props,
	const int i, const int j, const int k) {

  /* Properties of the zoom region. */
  const int zoom_cell_offset = zoom_props->tl_cell_offset;
  const double zoom_region_bounds[6] = {
      zoom_props->region_bounds[0], zoom_props->region_bounds[1],
      zoom_props->region_bounds[2], zoom_props->region_bounds[3],
      zoom_props->region_bounds[4], zoom_props->region_bounds[5]};
  const double ih_x_zoom = zoom_props->iwidth[0];
  const double ih_y_zoom = zoom_props->iwidth[1];
  const double ih_z_zoom = zoom_props->iwidth[2];
  int cell_id;

  /* Are the passed coordinates within the zoom region? */
  if (x >= zoom_region_bounds[0] && x <= zoom_region_bounds[1] &&
      y >= zoom_region_bounds[2] && y <= zoom_region_bounds[3] &&
      z >= zoom_region_bounds[4] && z <= zoom_region_bounds[5]) {

    const int zoom_index =
        cell_getid(cdim, (x - zoom_region_bounds[0]) * ih_x_zoom,
                   (y - zoom_region_bounds[2]) * ih_y_zoom,
                   (z - zoom_region_bounds[4]) * ih_z_zoom);
#ifdef SWIFT_DEBUG_CHECKS
  if (zoom_index < 0 || zoom_index >= cdim[0]*cdim[1]*cdim[2])
    error("zoom_index out of range %i (%f %f %f)", cell_id, x, y, z);
#endif
    cell_id = zoom_cell_offset + zoom_index;
  } else {
    cell_id = cell_getid(cdim, i, j, k);
#ifdef SWIFT_DEBUG_CHECKS
  if (cell_id < 0 || cell_id >= cdim[0]*cdim[1]*cdim[2])
    error("cell_id out of range %i (i:%i j:%i k:%i) (x:%f y:%f z:%f)",
            cell_id, i, j, k, x, y, z);
#endif
  }

  return cell_id;
}
