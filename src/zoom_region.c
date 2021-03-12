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

  //  s->zoom_props->max_size =
  //      parser_get_opt_param_double(params, "ZoomRegion:max_size", FLT_MAX);
  //  s->zoom_props->boost_factor =
  //      parser_get_opt_param_double(params, "ZoomRegion:boost_factor", 1.0);
  }
#endif
}

/**
 * @brief Compute the extent/bounds of the zoom region.
 *
 * Constructed from the high-res dark matter particles only (swift_type_dark_matter).
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
