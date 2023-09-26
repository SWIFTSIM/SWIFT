//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_DELAUNAY_2D_H
#define SWIFTSIM_SHADOWSWIFT_DELAUNAY_2D_H

#include "./geometry.h"
#include "./triangle.h"

#include <float.h>

/**
 * @brief Delaunay tessellation.
 *
 * The tessellation contains all the triangles that make up the tessellation,
 * their connectivity is stored implicitly within the triangles themselves.
 */
struct delaunay {

  /* Box geometry: used to set up the initial vertices and triangles and to
   * convert input coordinates to integer coordinates. */

  /*! @brief Anchor of the box containing the delaunay tesselation. */
  double anchor[2];

  /*! @brief Side length of the box containing the delaunay tesselation. */
  double side;

  /*! @brief Inverse side length. */
  double inverse_side;

  /*! @brief Vertex positions. */
  double* vertices;

  /*! @brief Vertex positions, rescaled to the range 1-2. Kept here in case we
   *  want to adopt hybrid geometrical checks (floating point checks when safe,
   *  integer checks when there is a risk of numerical error leading to the
   *  wrong result) to speed things up. */
  double* rescaled_vertices;

  /*! @brief Integer vertices. These are the vertex coordinates that are
   *  actually used during the incremental construction. */
  unsigned long int* integer_vertices;

  /*! @brief Vertex-triangle connections. For every vertex in the tessellation,
   *  this array stores the index of a triangle that contains this vertex
   *  (usually one of the last triangles that was constructed with this vertex
   *  as a member vertex). This array is not required for the incremental
   *  construction algorithm itself, but is indispensable for the conversion
   *  from Delaunay tessellation to Voronoi grid, since it links each input
   *  vertex to one of the triangles it is connected to. Without this array, we
   *  would have no efficient way of finding a triangle that contains a given
   *  vertex. */
  int* vertex_triangles;

  /*! @brief Vertex-triangle connection indices. For every vertex-triangle pair
   *  stored in vertex_triangles, this array contains the index of the vertex
   *  within the vertex list of the triangle. This saves us from having to
   *  loop through the triangle vertices during Voronoi grid construction. */
  int* vertex_triangle_index;

  /*! @brief Index of the corresponding particles in their respective cells. */
  int* vertex_part_idx;

  /*! @brief Next available index within the vertex array. Corresponds to the
   *  actual size of the vertex array. */
  int vertex_index;

  /*! @brief Current size of the vertex array in memory. If vertex_size matches
   *  vertex_index, the memory buffer is full and needs to be expanded. */
  int vertex_size;

  /*! @brief Begin index of the normal vertices. This skips the 3 auxiliary
   *  vertices required for the incremental construction algorithm. */
  int vertex_start;

  /*! @brief End index of the normal vertices. This variable is set by calling
   *  delaunay_consolidate() and contains the offset of the ghost vertices
   *  within the vertex array. */
  int vertex_end;

  /*! @brief Triangles that make up the tessellation. */
  struct triangle* triangles;

  /*! @brief Next available index within the triangle array. Corresponds to the
   *  actual size of the triangle array. */
  int triangle_index;

  /*! @brief Current size of the triangle array in memory. If triangle_size
   *  matches triangle_index, the memory buffer is full and needs to be
   *  expanded. */
  int triangle_size;

  /*! @brief Queue of triangles that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  triangles are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New triangles
   *  created when flipping invalid triangles are also added to the queue. */
  int* queue;

  /*! @brief Next available index in the queue. Determines both the actual size
   *  of the queue as the first element that will be popped from the queue. */
  int queue_index;

  /*! @brief Current size of the queue in memory. If queue_size matches
   *  queue_index, the memory buffer is full and needs to be expanded. */
  int queue_size;

  /*! @brief Index of the last triangle that was accessed. Used as initial
   *  guess for the triangle that contains the next vertex that will be added.
   *  If vertices are added in some sensible order (e.g. in Peano-Hilbert curve
   *  order) then this will greatly speed up the algorithm. */
  int last_triangle;

  /*! @brief Geometry variables. Auxiliary variables used by the exact integer
   *  geometry tests that need to be stored in between tests, since allocating
   *  and deallocating them for every test is too expensive. */
  struct geometry2d geometry;

  /*! @brief Cell sids keeping track of which cells contain a specific ghost
   * vertex. */
  int* ghost_cell_sids;

  /*! @brief Current used size of the ghost vertex bookkeeping arrays
   *  (and next valid index in this array). */
  int ghost_index;

  /*! @brief Current size in memory of the ghost vertex bookkeeping
   *  arrays. More memory needs to be allocated if ngb_index reaches this
   *  value. */
  int ghost_size;

  /*! @brief Array of booleans indicating whether or not neighbouring particles
   * have been tried to be added for a given sid. If this is 0 for a given sid,
   * this means that this cell should get the reflective boundary condition
   * applied for that sid. */
  unsigned long int sid_is_inside_face_mask;
};

inline static void delaunay_init_vertex(struct delaunay* restrict d, int v,
                                        double x, double y, int idx) {

  /* store a copy of the vertex coordinates */
  d->vertices[2 * v] = x;
  d->vertices[2 * v + 1] = y;

  /* compute the rescaled coordinates. We do this because floating point values
     in the range [1,2[ all have the same exponent (0), which guarantees that
     their mantissas form a linear sequence */
  double rescaled_x = 1. + (x - d->anchor[0]) * d->inverse_side;
  double rescaled_y = 1. + (y - d->anchor[1]) * d->inverse_side;

  delaunay_assert(rescaled_x >= 1.);
  delaunay_assert(rescaled_x < 2.);
  delaunay_assert(rescaled_y >= 1.);
  delaunay_assert(rescaled_y < 2.);

  /* store a copy of the rescaled coordinates to apply non-exact tests */
  d->rescaled_vertices[2 * v] = rescaled_x;
  d->rescaled_vertices[2 * v + 1] = rescaled_y;

  /* convert the rescaled coordinates to integer coordinates and store these */
  d->integer_vertices[2 * v] = delaunay_double_to_int(rescaled_x);
  d->integer_vertices[2 * v + 1] = delaunay_double_to_int(rescaled_y);

  /* initialise the variables that keep track of the link between vertices and
     triangles.
     We use negative values so that we can later detect missing links. */
  d->vertex_triangles[v] = -1;
  d->vertex_triangle_index[v] = -1;
  d->vertex_part_idx[v] = idx;
}

/**
 * @brief Add a new vertex with the given coordinates.
 *
 * This function first makes sure there is sufficient memory to store the
 * vertex and all its properties. It then converts the vertex coordinates from
 * box coordinates to coordinates in the range [1,2[ and stores the 52-bit
 * mantissas of these coordintes as integer coordinates.
 * The function also initialises all vertex properties to sensible values.
 *
 * @param d Delaunay tessellation.
 * @param x Horizontal coordinate of the vertex.
 * @param y Vertical coordinate of the vertex.
 * @param idx The index of the corresponding SWIFT particle in its cell.
 * @return Index of the new vertex within the vertex array.
 */
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x,
                                      double y, int idx) {

  /* check the size of the vertex arrays against the allocated memory size */
  if (d->vertex_index == d->vertex_size) {
    /* dynamically grow the size of the arrays with a factor 2 */
    d->vertex_size <<= 1;
    d->vertices = (double*)swift_realloc("delaunay", d->vertices,
                                         d->vertex_size * 2 * sizeof(double));
    d->rescaled_vertices = (double*)swift_realloc(
        "delaunay", d->rescaled_vertices, d->vertex_size * 2 * sizeof(double));
    d->integer_vertices = (unsigned long int*)swift_realloc(
        "delaunay", d->integer_vertices,
        d->vertex_size * 2 * sizeof(unsigned long int));
    d->vertex_triangles = (int*)swift_realloc("delaunay", d->vertex_triangles,
                                              d->vertex_size * sizeof(int));
    d->vertex_triangle_index = (int*)swift_realloc(
        "delaunay", d->vertex_triangle_index, d->vertex_size * sizeof(int));
    d->vertex_part_idx = (int*)swift_realloc("delaunay", d->vertex_part_idx,
                                             d->vertex_size * sizeof(int));
  }

  delaunay_init_vertex(d, d->vertex_index, x, y, idx);

  /* return the vertex index and then increase it by 1.
     After this operation, vertex_index will correspond to the size of the
     vertex arrays and is also the index of the next vertex that will be
     created. */
  return d->vertex_index++;
}

/**
 * @brief Claim a new triangle in the triangle array.
 *
 * This function first ensures that the triangle array is still large enough
 * to hold all triangles, and then returns the index of the next available
 * triangle.
 *
 * @param d Delaunay tessellation.
 * @return Index of the next available triangle in the triangle array.
 */
inline static int delaunay_new_triangle(struct delaunay* restrict d) {

  /* check that we still have triangles available */
  if (d->triangle_index == d->triangle_size) {
    /* no: increase the size of the triangle array with a factor 2 and
       reallocate it in memory */
    d->triangle_size <<= 1;
    d->triangles = (struct triangle*)swift_realloc(
        "delaunay", d->triangles, d->triangle_size * sizeof(struct triangle));
  }

  /* return the triangle index and then increase it by 1.
     After this operation, triangle_index will correspond to the size of the
     triangle array and is also the index of the next triangle that will be
     created. */
  return d->triangle_index++;
}

/**
 * @brief Perform an (expensive) check on the tessellation to see that it is
 * still valid.
 *
 * This function will iterate over all triangles (except the 3 dummy triangles)
 * of the tessellation. For each triangle, it checks that all neighbour
 * relations for that triangle are set correctly (i.e. if triangle A is a
 * neighbour of B, B should always also be a neighbour of A; if triangle A
 * thinks it is the ith neighbour in triangle B, it should also be found in the
 * ith neighbour position in triangle B), and that none fo the vertices of
 * neighbouring triangles are within the circumcircle of the triangle, which
 * would violate the Delaunay criterion.
 *
 * Finally, this function also checks the vertex-triangle information by
 * making sure that all triangle indices stored in vertex_triangles and
 * vertex_triangle_index are correct.
 *
 * The function will abort with an error if any problems are found with the
 * tessellation.
 *
 * This function returns immediately if DELAUNAY_CHECKS in inactive. It adds
 * a significant extra runtime cost and should never be used in production runs.
 *
 * @param d Delaunay tessellation (note that this parameter cannot be const as
 * one might suspect, since the geometrical test variables are used to test
 * the Delaunay criterion and they cannot be read-only).
 */
inline static void delaunay_check_tessellation(struct delaunay* restrict d) {

#ifndef DELAUNAY_CHECKS
  /* return immediately if we do not want to run the expensive checks */
  return;
#endif

  /* Loop over all triangles. We skip the dummy triangles, since their
     neighbour relations are (on purpose) incorrect. */
  for (int i = 3; i < d->triangle_index; ++i) {
    int vt1_0 = d->triangles[i].vertices[0];
    int vt1_1 = d->triangles[i].vertices[1];
    int vt1_2 = d->triangles[i].vertices[2];

    const unsigned long* al = &d->integer_vertices[2 * vt1_0];
    const unsigned long* bl = &d->integer_vertices[2 * vt1_1];
    const unsigned long* cl = &d->integer_vertices[2 * vt1_2];

    const double* ad = &d->rescaled_vertices[2 * vt1_0];
    const double* bd = &d->rescaled_vertices[2 * vt1_1];
    const double* cd = &d->rescaled_vertices[2 * vt1_2];

    /* loop over the 3 neighbours of this triangle */
    for (int j = 0; j < 3; ++j) {
      int ngb = d->triangles[i].neighbours[j];
      /* skip potential dummy neighbours, since they have an invalid vertex
         that will break the Delaunay test. */
      if (ngb < 3) continue;
      int i0 = d->triangles[i].index_in_neighbour[j];

      /* check the mutual neighbour relations */
      if (d->triangles[ngb].neighbours[i0] != i) {
        warning("Wrong neighbour!\n");
        warning("Triangle %i: %i %i %i\n", i, vt1_0, vt1_1, vt1_2);
        warning("Neighbours: %i %i %i\n", d->triangles[i].neighbours[0],
                d->triangles[i].neighbours[1], d->triangles[i].neighbours[2]);
        warning("Index in neighbour: %i %i %i\n",
                d->triangles[i].index_in_neighbour[0],
                d->triangles[i].index_in_neighbour[1],
                d->triangles[i].index_in_neighbour[2]);
        warning("Neighbour triangle %i: %i %i %i\n", ngb,
                d->triangles[ngb].vertices[0], d->triangles[ngb].vertices[1],
                d->triangles[ngb].vertices[2]);
        fprintf(
            stderr, "Neighbours: %i %i %i\n", d->triangles[ngb].neighbours[0],
            d->triangles[ngb].neighbours[1], d->triangles[ngb].neighbours[2]);
        error("Index in neighbour: %i %i %i\n",
              d->triangles[ngb].index_in_neighbour[0],
              d->triangles[ngb].index_in_neighbour[1],
              d->triangles[ngb].index_in_neighbour[2]);
      }

      /* check the neighbour index information */
      if (d->triangles[ngb].index_in_neighbour[i0] != j) {
        warning("Wrong neighbour!\n");
        warning("Triangle %i: %i %i %i\n", i, vt1_0, vt1_1, vt1_2);
        warning("Neighbours: %i %i %i\n", d->triangles[i].neighbours[0],
                d->triangles[i].neighbours[1], d->triangles[i].neighbours[2]);
        warning("Index in neighbour: %i %i %i\n",
                d->triangles[i].index_in_neighbour[0],
                d->triangles[i].index_in_neighbour[1],
                d->triangles[i].index_in_neighbour[2]);
        warning("Neighbour triangle %i: %i %i %i\n", ngb,
                d->triangles[ngb].vertices[0], d->triangles[ngb].vertices[1],
                d->triangles[ngb].vertices[2]);
        fprintf(
            stderr, "Neighbours: %i %i %i\n", d->triangles[ngb].neighbours[0],
            d->triangles[ngb].neighbours[1], d->triangles[ngb].neighbours[2]);
        error("Index in neighbour: %i %i %i\n",
              d->triangles[ngb].index_in_neighbour[0],
              d->triangles[ngb].index_in_neighbour[1],
              d->triangles[ngb].index_in_neighbour[2]);
      }

      /* now get the vertex that is not shared between the triangle and this
         neighbour and make sure it is not in the triangle's circumcircle */
      int vt2_0 = d->triangles[ngb].vertices[i0];

      const unsigned long* dl = &d->integer_vertices[2 * vt2_0];
      const double* dd = &d->rescaled_vertices[2 * vt2_0];

      int testi = geometry2d_in_circle_adaptive(&d->geometry, al, bl, cl, dl,
                                                ad, bd, cd, dd);

      if (testi > 0) {
        warning("Wrong triangle!\n");
        warning("Triangle %i: %i (%g %g) %i (%g %g) %i (%g %g)\n", i, vt1_0,
                d->vertices[2 * vt1_0], d->vertices[2 * vt1_0 + 1], vt1_1,
                d->vertices[2 * vt1_1], d->vertices[2 * vt1_1 + 1], vt1_2,
                d->vertices[2 * vt1_2], d->vertices[2 * vt1_2 + 1]);
        warning("Opposite vertex: %i (%g %g)\n", vt2_0, d->vertices[2 * vt2_0],
                d->vertices[2 * vt2_0 + 1]);
        error("Test result: %i\n", testi);
      }
    }
  }

  /* Loop over all vertex_triangle elements and check that they indeed link
     to triangles that contain their respective vertex */
    for (int i = 0; i < d->vertex_index; ++i) {
      int t = d->vertex_triangles[i];
      int vi = d->vertex_triangle_index[i];
      if (d->triangles[t].vertices[vi] != i) {
        error("Wrong vertex-triangle connection!");
      }
    }
}

/**
 * @brief Reset the Delaunay tessellation without reallocating memory.
 *
 * It sets up a large triangle that contains the entire simulation box and
 * additional buffer space to deal with boundary ghost vertices, and 3
 * additional dummy triangles that provide valid neighbours for the 3 sides of
 * this triangle (these dummy triangles themselves have an invalid tip vertex
 * and are therefore simply placeholders).
 *
 * @param d Delaunay tessellation.
 */
inline static void delaunay_reset(struct delaunay* restrict d,
                                  const double* cell_loc,
                                  const double* cell_width, int vertex_size) {

  if (vertex_size == 0) {
    /* Don't bother for empty cells */
    return;
  }

  /* by overwriting the indices, we invalidate all arrays without changing their
     allocated size */
  /* We reserve vertex_size spots for local vertices. */
  d->vertex_index = 0;
  d->triangle_index = 0;
  d->queue_index = 0;
  d->ghost_index = 0;

  /* determine the size of a box large enough to accommodate the entire
     simulation volume and all possible ghost vertices required to deal with
     boundaries. Note that the box we use is quite big, since it is chosen to
     work for even the very degenerate case of a single vertex that could be
     anywhere in the box. Also note that we convert the generally rectangular
     box to a square. */
  /* We add an extra layer of padding of 1 times cell_width around our box to
   * compensate for particle movements (up to 1 cell side before a rebuild is
   * triggered). */
  double box_anchor[2] = {cell_loc[0] - 2. * cell_width[0],
                          cell_loc[1] - 2. * cell_width[1]};
  double box_side = 10. * fmax(cell_width[0], cell_width[1]);

  /* store the anchor and inverse side_length for the conversion from box
  coordinates to rescaled (integer) coordinates */
  d->anchor[0] = box_anchor[0];
  d->anchor[1] = box_anchor[1];
  d->side = box_side;
  /* the 1.e-13 makes sure converted values are in the range [1, 2[ instead of
   * [1,2] (unlike Springel, 2010) */
  d->inverse_side = (1. - DBL_EPSILON) / box_side;

  /* set up the large triangle and the 3 dummies */
  /* mind the orientation: counterclockwise w.r.t. the z-axis. */
  int v0 = delaunay_new_vertex(d, box_anchor[0], box_anchor[1], -1);
  delaunay_log("Creating vertex %i: %g %g", v0, box_anchor[0], box_anchor[1]);
  int v1 = delaunay_new_vertex(d, box_anchor[0] + box_side, box_anchor[1], -1);
  delaunay_log("Creating vertex %i: %g %g", v1, box_anchor[0] + box_side,
               box_anchor[1]);
  int v2 = delaunay_new_vertex(d, box_anchor[0], box_anchor[1] + box_side, -1);
  delaunay_log("Creating vertex %i: %g %g", v2, box_anchor[0],
               box_anchor[1] + box_side);
  d->vertex_start = d->vertex_index;
  d->vertex_end = 0;

  /* we also create 3 dummy triangles with a fake tip (just to create valid
     neighbours for the big triangle */
  int dummy0 = delaunay_new_triangle(d);
  int dummy1 = delaunay_new_triangle(d);
  int dummy2 = delaunay_new_triangle(d);
  int first_triangle = delaunay_new_triangle(d);
  delaunay_log("Creating triangle %i: %i %i %i", dummy0, v1, v0, -1);
  triangle_init(&d->triangles[dummy0], v1, v0, -1);
  triangle_swap_neighbour(&d->triangles[dummy0], 2, first_triangle, 2);
  delaunay_log("Creating triangle %i: %i %i %i", dummy1, v2, v1, -1);
  triangle_init(&d->triangles[dummy1], v2, v1, -1);
  triangle_swap_neighbour(&d->triangles[dummy1], 2, first_triangle, 0);
  delaunay_log("Creating triangle %i: %i %i %i", dummy2, v0, v2, -1);
  triangle_init(&d->triangles[dummy2], v0, v2, -1);
  triangle_swap_neighbour(&d->triangles[dummy2], 2, first_triangle, 1);
  delaunay_log("Creating triangle %i: %i %i %i", first_triangle, v0, v1, v2);
  triangle_init(&d->triangles[first_triangle], v0, v1, v2);
  triangle_swap_neighbour(&d->triangles[first_triangle], 0, dummy1, 2);
  triangle_swap_neighbour(&d->triangles[first_triangle], 1, dummy2, 2);
  triangle_swap_neighbour(&d->triangles[first_triangle], 2, dummy0, 2);

  /* set up the vertex-triangle links for the initial triangle (not for the
     dummies) */
  d->vertex_triangles[v0] = first_triangle;
  d->vertex_triangle_index[v0] = 0;
  d->vertex_triangles[v1] = first_triangle;
  d->vertex_triangle_index[v1] = 1;
  d->vertex_triangles[v2] = first_triangle;
  d->vertex_triangle_index[v2] = 2;

  d->last_triangle = first_triangle;

  /* Initialize the sid mask so that the sid's with a z-component are already 1
   * (are all threated as internal), since this is the 2D version.
   * Also sid=13 does not correspond to a face and is always set to 1 as well.
   * We only set the sid's corresponding to the cardinal directions to 0
   * (only faces perpendicular to one of the axes can be boundary faces). */
  d->sid_is_inside_face_mask = DEFAULT_SID_MASK;

  /* Perform potential log output and sanity checks */
  delaunay_log("Post init or reset check");
  delaunay_check_tessellation(d);
}

/**
 * @brief Initialize the Delaunay tessellation.
 *
 * This function allocates memory for all arrays that make up the tessellation
 * and initializes the variables used for bookkeeping.
 *
 * @param d Delaunay tesselation.
 * @param hs Spatial extents of the simulation box.
 * @param vertex_size Initial size of the vertex array.
 */
inline static struct delaunay* delaunay_malloc(const double* cell_loc,
                                               const double* cell_width,
                                               int vertex_size) {
  /* Don't bother setting up a Delaunay tessellation for empty cells */
  if (vertex_size == 0) {
    return NULL;
  }

  struct delaunay* d = malloc(sizeof(struct delaunay));

  /* allocate memory for the vertex arrays */
  d->vertices =
      (double*)swift_malloc("delaunay", vertex_size * 2 * sizeof(double));
  d->rescaled_vertices =
      (double*)swift_malloc("delaunay", vertex_size * 2 * sizeof(double));
  d->integer_vertices = (unsigned long int*)swift_malloc(
      "delaunay", vertex_size * 2 * sizeof(unsigned long int));
  d->vertex_triangles =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_triangle_index =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_part_idx =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_size = vertex_size;

  /* allocate memory for the triangle array */
  d->triangle_size = 6 * vertex_size;
  d->triangles = (struct triangle*)swift_malloc(
      "delaunay", d->triangle_size * sizeof(struct triangle));

  /* allocate memory for the queue (note that the queue size of 10 was chosen
     arbitrarily, and a proper value should be chosen based on performance
     measurements) */
  d->queue = (int*)swift_malloc("delaunay", 10 * sizeof(int));
  d->queue_size = 10;

  d->ghost_cell_sids =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->ghost_size = vertex_size;

  /* initialise the structure used to perform exact geometrical tests */
  geometry_init(&d->geometry);

  delaunay_reset(d, cell_loc, cell_width, vertex_size);

  return d;
}

/**
 * @brief Free up all memory associated with the Delaunay tessellation.
 *
 * @param d Delaunay tessellation.
 */
inline static void delaunay_destroy(struct delaunay* restrict d) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(d->vertices != NULL);
#endif
  swift_free("delaunay", d->vertices);
  swift_free("delaunay", d->rescaled_vertices);
  swift_free("delaunay", d->integer_vertices);
  swift_free("delaunay", d->vertex_triangles);
  swift_free("delaunay", d->vertex_triangle_index);
  swift_free("delaunay", d->vertex_part_idx);
  swift_free("delaunay", d->triangles);
  swift_free("delaunay", d->queue);
  swift_free("delaunay", d->ghost_cell_sids);
  geometry_destroy(&d->geometry);

  /* Free delaunay struct itself */
  free(d);
}

/**
 * @brief Get the radius of the circumcircle of the given triangle.
 *
 * @param d Delaunay tessellation.
 * @param t Triangle index.
 * @return Radius of the circumcircle of the triangle.
 */
inline static double delaunay_get_radius(const struct delaunay* restrict d,
                                         int t) {
  int v0 = d->triangles[t].vertices[0];
  int v1 = d->triangles[t].vertices[1];
  int v2 = d->triangles[t].vertices[2];

  double v0x = d->vertices[2 * v0];
  double v0y = d->vertices[2 * v0 + 1];
  double v1x = d->vertices[2 * v1];
  double v1y = d->vertices[2 * v1 + 1];
  double v2x = d->vertices[2 * v2];
  double v2y = d->vertices[2 * v2 + 1];

  double ax = v1x - v0x;
  double ay = v1y - v0y;
  double bx = v2x - v0x;
  double by = v2y - v0y;

  double D = 2. * (ax * by - ay * bx);
  double a2 = ax * ax + ay * ay;
  double b2 = bx * bx + by * by;
  double Rx = (by * a2 - ay * b2) / D;
  double Ry = (ax * b2 - bx * a2) / D;

  return sqrt(Rx * Rx + Ry * Ry);
}

/**
 * @brief Consolidate the Delaunay tessellation. This signals the end of the
 * addition of normal vertices. All vertices added after this point are
 * considered to be ghost vertices.
 *
 * @param d Delaunay tessellation.
 */
inline static void delaunay_consolidate(struct delaunay* restrict d) {
  d->vertex_end = d->vertex_index;
  delaunay_check_tessellation(d);
}

/**
 * @brief Computes the delaunay search radii for a list of particles.
 *
 * For a given generator, we must loop over all the triangles connected to this
 * generator and compute the maximal circum-radius.
 * This is fairly straightforward in the 2D case.
 *
 * @param d The #delaunay tesselation
 * @param pid The indices of the particles to compute the search radius for
 * @param count The number of particles.
 * @param r (return) The search radii.
 * */
inline static void delaunay_get_search_radii(struct delaunay* restrict d,
                                             const struct part* restrict parts,
                                             const int* restrict pid, int count,
                                             /*return*/ double* restrict r) {
  /* Loop over the particles */
  for (int i = 0; i < count; i++) {
    int vi = parts[pid[i]].geometry.delaunay_vertex;
    ;
    int t0 = d->vertex_triangles[vi];
    int vi0 = d->vertex_triangle_index[vi];
    int vi0p1 = (vi0 + 1) % 3;
    double search_radius = 2. * delaunay_get_radius(d, t0);
    int t1 = d->triangles[t0].neighbours[vi0p1];
    int vi1 = d->triangles[t0].index_in_neighbour[vi0p1];
    while (t1 != t0) {
      search_radius = max(search_radius, 2. * delaunay_get_radius(d, t1));
      int vi1p2 = (vi1 + 2) % 3;
      vi1 = d->triangles[t1].index_in_neighbour[vi1p2];
      t1 = d->triangles[t1].neighbours[vi1p2];
    }
    r[i] = search_radius;
  }
}

/**
 * @brief Randomly choose a neighbour from the given two options.
 *
 * @param ngb0 First option.
 * @param ngb1 Second option.
 * @return One of the neighbours. Both options have a 50% likelihood of being
 * returned based on random sampling using rand().
 */
inline static int delaunay_choose(int ngb0, int ngb1) {
  return (rand() % 2) ? ngb0 : ngb1;
}

inline static int delaunay_choose_line_segment_intersection(
    int ngb0, int ngb1, const double* restrict face0_0,
    const double* restrict face0_1, const double* restrict v,
    const double* restrict centroid) {
  /* test line segment */
  if (geometry2d_test_line_segment_intersection(face0_0, face0_1, v,
                                                centroid)) {
    return ngb0;
  }
  return ngb1;
}

/**
 * @brief Test if the given vertex is within the given triangle.
 *
 * This function has 3 possible return values:
 *  - 0: the vertex is not inside the triangle. In this case, t_new is set to
 *    the index of a triangle that is closer to the vertex and should be tested
 *    next.
 *  - 1: the vertex is inside the triangle. t_new is set to t.
 *  - 2: the vertex is on one of the edges of the triangle. t_new is set to the
 *    index on the other side of the edge, while ngb_index is set to the index
 *    of t_new within the neighbour list of t.
 *
 * This function uses arbitrarily exact geometric tests and therefore always
 * returns the correct answer if the vertex coordinates are expressed in
 * (rescaled) integer coordinates. The answer could still be different from the
 * answer that would be obtained when expressing the vertex coordinates in the
 * original box coordinates. However, the exactness of the results in rescaled
 * coordinates guarantees that all geometric tests will always be consistent;
 * if a vertex V is (very) close to the edge between triangles A and B, then the
 * tests for A and B will always agree on whether V is in A or B. Without exact
 * tests, the test for A and B could both be negative, causing the algorithm
 * to get stuck in an endless loop that jumps between both triangles.
 *
 * The geometrical tests used by this function require a certain level of
 * geometrical consistency that will make them fail if the vertices of t are
 * colinear or if the vertex v lies on top of one of the vertices of the
 * triangle. This function will detect these cases and will abort with an error
 * message if they occur.
 *
 * @param d Delaunay tessellation.
 * @param v Vertex index.
 * @param t Triangle index.
 * @param t_new Variable in which the next triangle index is stored (see cases
 * above for the exact meaning).
 * @param ngb_index Variable in which the index of t_new within the neighbour
 * list of t is stored, but only in case 2 (v lies on the edge between t and
 * t_new).
 * @return 0, 1 or 2, as explained above.
 */
inline static int delaunay_test_point_inside_triangle(
    struct delaunay* restrict d, int v, int t, int* t_new, int* ngb_index) {

  delaunay_log("Checking if vertex %i is inside triange %i", v, t);

  /* make sure we are not testing a dummy triangle. Since the dummies lie
     outside the triangle that is supposed to contain all vertices (including
     all possible ghosts), this should never happen. */
  delaunay_assert(t > 2);

  /* get the vertex positions
     a is the vertex we want to test
     b, c, d are the vertices of the triangle */

  int vt0 = d->triangles[t].vertices[0];
  int vt1 = d->triangles[t].vertices[1];
  int vt2 = d->triangles[t].vertices[2];
  delaunay_log("Triangle vertices: %i %i %i", vt0, vt1, vt2);

  const unsigned long* al = &d->integer_vertices[2 * v];
  const unsigned long* bl = &d->integer_vertices[2 * vt0];
  const unsigned long* cl = &d->integer_vertices[2 * vt1];
  const unsigned long* dl = &d->integer_vertices[2 * vt2];

  const double* ad = &d->rescaled_vertices[2 * v];
  const double* bd = &d->rescaled_vertices[2 * vt0];
  const double* cd = &d->rescaled_vertices[2 * vt1];
  const double* dd = &d->rescaled_vertices[2 * vt2];

  double centroid[2];
  geometry2d_compute_centroid_triangle(bd[0], bd[1], cd[0], cd[1], dd[0], dd[1],
                                       centroid);

  /* test all 3 edges of the triangle. This code could potentially be sped up
     by including conditional statements that only execute necessary tests.
     However, it is unclear whether such an approach is beneficial because it
     would introduce additional branching of the code. */
  int testi0 = geometry2d_orient_adaptive(&d->geometry, cl, dl, al, cd, dd, ad);
  int testi1 = geometry2d_orient_adaptive(&d->geometry, dl, bl, al, dd, bd, ad);
  int testi2 = geometry2d_orient_adaptive(&d->geometry, bl, cl, al, bd, cd, ad);

  /* to simplify the evaluation of the various scenarios, we combine the test
     results into a single test value that is then checked within a switch
     statement. This approach is overall cleaner than using a bunch of nested
     conditionals. */
  int testsum = (testi0 + 1) << 4 | (testi1 + 1) << 2 | (testi2 + 1);
  switch (testsum) {
    case 1:
    case 2:
      *t_new = delaunay_choose_line_segment_intersection(
          d->triangles[t].neighbours[0], d->triangles[t].neighbours[1], cd, dd,
          ad, centroid);
      return 0;
    case 4:
    case 8:
      *t_new = delaunay_choose_line_segment_intersection(
          d->triangles[t].neighbours[0], d->triangles[t].neighbours[2], cd, dd,
          ad, centroid);
      return 0;
    case 5:
    case 6:
    case 9:
    case 10:
      /* testi0 is negative */
      *t_new = d->triangles[t].neighbours[0];
      return 0;
    case 16:
    case 32:
      *t_new = delaunay_choose_line_segment_intersection(
          d->triangles[t].neighbours[1], d->triangles[t].neighbours[2], bd, dd,
          ad, centroid);
      return 0;
    case 17:
    case 18:
    case 33:
    case 34:
      /* testi1 is negative */
      *t_new = d->triangles[t].neighbours[1];
      return 0;
    case 20:
    case 24:
    case 36:
    case 40:
      /* testi2 is negative */
      *t_new = d->triangles[t].neighbours[2];
      return 0;
    case 26:
      /* testi0 is zero */
      *t_new = d->triangles[t].neighbours[0];
      *ngb_index = 0;
      return 2;
    case 38:
      /* testi1 is zero */
      *t_new = d->triangles[t].neighbours[1];
      *ngb_index = 1;
      return 2;
    case 41:
      /* testi2 is zero */
      *t_new = d->triangles[t].neighbours[2];
      *ngb_index = 2;
      return 2;
    case 42:
      /* all tests returned 1 (testsum is 32 + 8 + 2) */
      *t_new = t;
      return 1;
    default:
      /* we have ended up in a geometrically impossible scenario. This can
         happen if the triangle vertices are colinear, or if one or multiple
         vertices coincide. */
      return -1;
  }
}

/**
 * @brief Add the given triangle to the queue of triangles that need checking.
 *
 * @param d Delaunay tessellation.
 * @param t Triangle index.
 */
inline static void delaunay_enqueue(struct delaunay* restrict d, int t) {

  /* make sure there is sufficient space in the queue */
  if (d->queue_index == d->queue_size) {
    /* there isn't: increase the size of the queue with a factor 2. */
    d->queue_size <<= 1;
    d->queue =
        (int*)swift_realloc("delaunay", d->queue, d->queue_size * sizeof(int));
  }

  delaunay_log("Enqueuing triangle %i and vertex 2", t);

  /* add the triangle to the queue and advance the queue index */
  d->queue[d->queue_index] = t;
  ++d->queue_index;
}

/**
 * @brief Pop the next triangle to check from the end of the queue.
 *
 * If no more triangles are queued, this function returns a negative value.
 *
 * Note that the returned triangle index is effectively removed from the queue
 * and will be overwritten by subsequent calls to delaunay_enqueue().
 *
 * @param d Delaunay tessellation.
 * @return Index of the next triangle to test, or -1 if the queue is empty.
 */
inline static int delaunay_queue_pop(struct delaunay* restrict d) {
  if (d->queue_index > 0) {
    --d->queue_index;
    return d->queue[d->queue_index];
  } else {
    return -1;
  }
}

/**
 * @brief Check the Delaunay criterion for the given triangle.
 *
 * Per convention, we assume that the check is triggered by inserting the
 * final vertex of the given triangle in the tessellation. In this case, only
 * one check is required: we need to test the circumcircle criterion for the
 * top vertex of the neighbouring triangle opposite the new vertex (so the
 * vertex that is not shared between this triangle and that neighbour).
 *
 * This function is also responsible for performing the edge flip that restores
 * Delaunayhood in case the test fails. This boils down to erasing the edge that
 * is shared between this triangle and its neighbour and constructing a new
 * edge between the newly inserted vertex and the opposite vertex in the
 * neighbour. After this operation, the two newly created triangles need to be
 * checked for Delaunayhood themselves and are added to the end of the queue.
 *
 * As for delaunay_test_point_inside_triangle(), this function makes use of
 * arbitrarily exact geometric tests that guarantee consistency of test results
 * within the rescaled (integer) coordinate frame. This is to avoid situations
 * in which newly created triangles would appear to violate the Delaunay
 * criterion and restore a previously flipped edge, leading to an endless loop
 * and a catastrophic deadlock of the algorithm.
 *
 * @param d Delaunay tessellation.
 * @param t Triangle index.
 */
inline static void delaunay_check_triangle(struct delaunay* restrict d, int t) {

  delaunay_log("Checking triangle %i and vertex 2", t);

  int t2 = d->triangles[t].neighbours[2];
  /* check if we have a neighbour that can be checked (dummies are not real and
     should not be tested) */
  if (t2 < 3) {
    delaunay_log("No neighbour to check");
    return;
  }

  /* get the vertex positions
     a, b, c are the vertices of the triangle
     d is the vertex we want to test */
  int vt1_0 = d->triangles[t].vertices[0];
  int vt1_1 = d->triangles[t].vertices[1];
  int vt1_2 = d->triangles[t].vertices[2];

  delaunay_log("Vertices: %i %i %i", vt1_0, vt1_1, vt1_2);

  int i0 = d->triangles[t].index_in_neighbour[2];
  int vt2_0 = d->triangles[t2].vertices[i0];

  delaunay_log("Opposite vertex: %i", vt2_0);

  const unsigned long* al = &d->integer_vertices[2 * vt1_0];
  const unsigned long* bl = &d->integer_vertices[2 * vt1_1];
  const unsigned long* cl = &d->integer_vertices[2 * vt1_2];
  const unsigned long* dl = &d->integer_vertices[2 * vt2_0];

  const double* ad = &d->rescaled_vertices[2 * vt1_0];
  const double* bd = &d->rescaled_vertices[2 * vt1_1];
  const double* cd = &d->rescaled_vertices[2 * vt1_2];
  const double* dd = &d->rescaled_vertices[2 * vt2_0];

  int testi = geometry2d_in_circle_adaptive(&d->geometry, al, bl, cl, dl, ad,
                                            bd, cd, dd);

  /* check if we need to flip the edge */
  if (testi > 0) {
    delaunay_log("Flipping triangle");

    /* flip edge. We first obtain the indices of the edge vertices in the
       neighbouring triangle (assuming they are correctly oriented) */
    int i1 = (i0 + 1) % 3;
    int i2 = (i0 + 2) % 3;
    delaunay_assert(d->triangles[t2].vertices[i1] == vt1_1);
    delaunay_assert(d->triangles[t2].vertices[i2] == vt1_0);

    /* obtain all the neighbouring information for the 4 neighbours of the
       two triangles that share the edge we want to flip */
    int ngb0 = d->triangles[t].neighbours[1];
    int ngbi0 = d->triangles[t].index_in_neighbour[1];

    int ngb1 = d->triangles[t].neighbours[0];
    int ngbi1 = d->triangles[t].index_in_neighbour[0];

    int ngb2 = d->triangles[t2].neighbours[i2];
    int ngbi2 = d->triangles[t2].index_in_neighbour[i2];

    int ngb3 = d->triangles[t2].neighbours[i1];
    int ngbi3 = d->triangles[t2].index_in_neighbour[i1];

    /* now create the 2 new triangles. We reuse the triangle indices of the
       original triangles. It helps enormously to make a sketch of the 2
       triangles in order to see that the neighbour relations below are
       correct. Also note that for future checks to work, we need to stick to
       the convention that the newly inserted vertex (vt1_2) is the last vertex
       of the new triangles. */
    delaunay_log("Creating triangle %i: %i %i %i", t, vt1_0, vt2_0, vt1_2);
    triangle_init(&d->triangles[t], vt1_0, vt2_0, vt1_2);
    triangle_swap_neighbour(&d->triangles[t], 0, t2, 1);
    triangle_swap_neighbour(&d->triangles[t], 1, ngb0, ngbi0);
    triangle_swap_neighbour(&d->triangles[ngb0], ngbi0, t, 1);
    triangle_swap_neighbour(&d->triangles[t], 2, ngb3, ngbi3);
    triangle_swap_neighbour(&d->triangles[ngb3], ngbi3, t, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t2, vt2_0, vt1_1, vt1_2);
    triangle_init(&d->triangles[t2], vt2_0, vt1_1, vt1_2);
    triangle_swap_neighbour(&d->triangles[t2], 0, ngb1, ngbi1);
    triangle_swap_neighbour(&d->triangles[ngb1], ngbi1, t2, 0);
    triangle_swap_neighbour(&d->triangles[t2], 1, t, 0);
    triangle_swap_neighbour(&d->triangles[t2], 2, ngb2, ngbi2);
    triangle_swap_neighbour(&d->triangles[ngb2], ngbi2, t2, 2);

    /* update the vertex-triangle links */
    d->vertex_triangles[vt1_0] = t;
    d->vertex_triangle_index[vt1_0] = 0;
    d->vertex_triangles[vt1_1] = t2;
    d->vertex_triangle_index[vt1_1] = 1;
    d->vertex_triangles[vt1_2] = t2;
    d->vertex_triangle_index[vt1_2] = 2;
    d->vertex_triangles[vt2_0] = t2;
    d->vertex_triangle_index[vt2_0] = 0;

    /* add the new triangles to the queue for checking */
    delaunay_enqueue(d, t);
    delaunay_enqueue(d, t2);

    /* store the index of t2 as starting point for future vertex insertions */
    d->last_triangle = t2;
  }
}

/**
 * @brief Test the Delaunay criterion for all triangles in the queue after a
 * vertex insertion event.
 *
 * @param d Delaunay tessellation.
 */
inline static void delaunay_check_triangles(struct delaunay* restrict d) {
  int t = delaunay_queue_pop(d);
  while (t >= 0) {
    delaunay_check_triangle(d, t);
    t = delaunay_queue_pop(d);
  }
}

/**
 * @brief Add a new vertex to the tessellation.
 *
 * This function first adds the vertex to the internal list, computing its
 * rescaled (integer) coordinates. It then locates the triangle in the current
 * tessellation that contains the new vertex. This triangle is split into 3 new
 * triangles by combining the 3 edges of the original triangle with the new
 * vertex. If the new vertex degenerately lies on the edge between two existing
 * triangles, those two triangles are split into 4 new triangles by connecting
 * the new vertex with the two vertices that are not on that edge.
 *
 * After the creation of 3 or 4 new triangles, all new triangles are added to a
 * queue and are tested against the Delaunay criterion. If required, additional
 * edge flips are executed to restore Delaunayhood. To enable efficient testing,
 * it is important that all new triangles satisfy two conventions: their
 * vertices should always be ordered counterclockwise, and the newly inserted
 * vertex should always be the final vertex of any newly created triangle.
 *
 * @param d Delaunay tessellation.
 * @param v Vertex index.
 * @param x Horizontal coordinate of the new vertex.
 * @param y Vertical coordinate of the new vertex.
 * @return 0 on success, -1 if the vertex already exists.
 */
inline static int delaunay_finalize_vertex(struct delaunay* restrict d, int v,
                                           double x, double y) {

  delaunay_log("Adding vertex with position %g %g", x, y);

  /* find the triangle that contains the new vertex. Starting from an initial
     guess (the last triangle accessed by the previous insertion), we test if
     the vertex is inside or on the edge of the triangle. If the test finds the
     vertex to be outside of the triangle, a better guess is obtained by taking
     a neighbouring triangle of the current triangle that is closer to the
     current guess. */
  int t0, t1, ngb_index;
  t0 = d->last_triangle;
  int flag = delaunay_test_point_inside_triangle(d, v, t0, &t1, &ngb_index);
  if (flag == -1) {
    return -1;
  }
  int count = 0;
  while (flag == 0) {
    t0 = t1;
    flag = delaunay_test_point_inside_triangle(d, v, t0, &t1, &ngb_index);
    ++count;
    delaunay_assert(count < d->triangle_index);
  }
  if (flag == -1) {
    return -1;
  }

  delaunay_log("Found triangle: %i (flag %i)", t0, flag);
  /* check whether we are in the normal (flag = 1) or degenerate (flag = 2)
     case */
  if (flag == 1) {
    /* normal case: split t0 into 3 new triangles */

    /* retrieve the properties of the old triangle */
    int tv0 = d->triangles[t0].vertices[0];
    int tv1 = d->triangles[t0].vertices[1];
    int tv2 = d->triangles[t0].vertices[2];
    int ngb0 = d->triangles[t0].neighbours[0];
    int ngb1 = d->triangles[t0].neighbours[1];
    int ngb2 = d->triangles[t0].neighbours[2];
    int ix_in_ngb0 = d->triangles[t0].index_in_neighbour[0];
    int ix_in_ngb1 = d->triangles[t0].index_in_neighbour[1];
    int ix_in_ngb2 = d->triangles[t0].index_in_neighbour[2];

    /* create 3 new triangles (drawing the triangles helps to check that the
       code below is indeed correct) */
    t1 = delaunay_new_triangle(d);
    int t2 = delaunay_new_triangle(d);
    delaunay_log("Creating triangle %i: %i %i %i", t0, tv0, tv1, v);
    triangle_init(&d->triangles[t0], tv0, tv1, v);
    triangle_swap_neighbour(&d->triangles[t0], 0, t1, 1);
    triangle_swap_neighbour(&d->triangles[t0], 1, t2, 0);
    triangle_swap_neighbour(&d->triangles[t0], 2, ngb2, ix_in_ngb2);
    triangle_swap_neighbour(&d->triangles[ngb2], ix_in_ngb2, t0, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t1, tv1, tv2, v);
    triangle_init(&d->triangles[t1], tv1, tv2, v);
    triangle_swap_neighbour(&d->triangles[t1], 0, t2, 1);
    triangle_swap_neighbour(&d->triangles[t1], 1, t0, 0);
    triangle_swap_neighbour(&d->triangles[t1], 2, ngb0, ix_in_ngb0);
    triangle_swap_neighbour(&d->triangles[ngb0], ix_in_ngb0, t1, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t2, tv2, tv0, v);
    triangle_init(&d->triangles[t2], tv2, tv0, v);
    triangle_swap_neighbour(&d->triangles[t2], 0, t0, 1);
    triangle_swap_neighbour(&d->triangles[t2], 1, t1, 0);
    triangle_swap_neighbour(&d->triangles[t2], 2, ngb1, ix_in_ngb1);
    triangle_swap_neighbour(&d->triangles[ngb1], ix_in_ngb1, t2, 2);

    /* update the vertex-triangle links */
    d->vertex_triangles[tv0] = t0;
    d->vertex_triangle_index[tv0] = 0;
    d->vertex_triangles[tv1] = t1;
    d->vertex_triangle_index[tv1] = 0;
    d->vertex_triangles[tv2] = t2;
    d->vertex_triangle_index[tv2] = 0;
    d->vertex_triangles[v] = t2;
    d->vertex_triangle_index[v] = 2;

    /* add the triangles to the queue for checking */
    delaunay_enqueue(d, t0);
    delaunay_enqueue(d, t1);
    delaunay_enqueue(d, t2);

    /* touch the last triangle; if nothing else happens, this triangle will be
       the initial guess for the next insertion */
    d->last_triangle = t2;
  } else {

    /* degenerate case: t0 and t1 are split into 4 new triangles */

    delaunay_log("Degenerate insertion!");

    delaunay_assert(d->triangles[t0].neighbours[ngb_index] == t1);

    /* the new vertex lies on the edge separating triangles t0 and t1
       this edge lies opposite vertex ngb_index
       we need to split t0 and t1 into 4 new triangles */
    int i0_1 = (ngb_index + 1) % 3;
    int i0_2 = (ngb_index + 2) % 3;

    int vt0_0 = d->triangles[t0].vertices[ngb_index];
    int vt0_1 = d->triangles[t0].vertices[i0_1];
    int vt0_2 = d->triangles[t0].vertices[i0_2];

    int i1_0 = d->triangles[t0].index_in_neighbour[ngb_index];
    int i1_1 = (i1_0 + 1) % 3;
    int i1_2 = (i1_0 + 2) % 3;

    delaunay_assert(d->triangles[t1].neighbours[i1_0] == t0);
    delaunay_assert(d->triangles[t1].vertices[i1_1] == vt0_2);
    delaunay_assert(d->triangles[t1].vertices[i1_2] == vt0_1);

    int vt1_0 = d->triangles[t1].vertices[i1_0];

    delaunay_log("Triangle %i: %i %i %i", t0, vt0_0, vt0_1, vt0_2);
    delaunay_log("Triangle %i: %i %i %i", t1, vt1_0, vt0_2, vt0_1);

    int ngb0_1 = d->triangles[t0].neighbours[i0_1];
    int ngbi0_1 = d->triangles[t0].index_in_neighbour[i0_1];
    int ngb0_2 = d->triangles[t0].neighbours[i0_2];
    int ngbi0_2 = d->triangles[t0].index_in_neighbour[i0_2];

    int ngb1_1 = d->triangles[t1].neighbours[i1_1];
    int ngbi1_1 = d->triangles[t1].index_in_neighbour[i1_1];
    int ngb1_2 = d->triangles[t1].neighbours[i1_2];
    int ngbi1_2 = d->triangles[t1].index_in_neighbour[i1_2];

    int t2 = delaunay_new_triangle(d);
    int t3 = delaunay_new_triangle(d);

    /* as always, a drawing is very helpful to check the code below */
    delaunay_log("Creating triangle %i: %i %i %i", t0, vt0_0, vt0_1, v);
    triangle_init(&d->triangles[t0], vt0_0, vt0_1, v);
    triangle_swap_neighbour(&d->triangles[t0], 0, t2, 1);
    triangle_swap_neighbour(&d->triangles[t0], 1, t1, 0);
    triangle_swap_neighbour(&d->triangles[t0], 2, ngb0_2, ngbi0_2);
    triangle_swap_neighbour(&d->triangles[ngb0_2], ngbi0_2, t0, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t1, vt0_2, vt0_0, v);
    triangle_init(&d->triangles[t1], vt0_2, vt0_0, v);
    triangle_swap_neighbour(&d->triangles[t1], 0, t0, 1);
    triangle_swap_neighbour(&d->triangles[t1], 1, t3, 0);
    triangle_swap_neighbour(&d->triangles[t1], 2, ngb0_1, ngbi0_1);
    triangle_swap_neighbour(&d->triangles[ngb0_1], ngbi0_1, t1, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t2, vt0_1, vt1_0, v);
    triangle_init(&d->triangles[t2], vt0_1, vt1_0, v);
    triangle_swap_neighbour(&d->triangles[t2], 0, t3, 1);
    triangle_swap_neighbour(&d->triangles[t2], 1, t0, 0);
    triangle_swap_neighbour(&d->triangles[t2], 2, ngb1_1, ngbi1_1);
    triangle_swap_neighbour(&d->triangles[ngb1_1], ngbi1_1, t2, 2);

    delaunay_log("Creating triangle %i: %i %i %i", t3, vt1_0, vt0_2, v);
    triangle_init(&d->triangles[t3], vt1_0, vt0_2, v);
    triangle_swap_neighbour(&d->triangles[t3], 0, t1, 1);
    triangle_swap_neighbour(&d->triangles[t3], 1, t2, 0);
    triangle_swap_neighbour(&d->triangles[t3], 2, ngb1_2, ngbi1_2);
    triangle_swap_neighbour(&d->triangles[ngb1_2], ngbi1_2, t3, 2);

    /* update the triangle-vertex connections */
    d->vertex_triangles[vt0_0] = t0;
    d->vertex_triangle_index[vt0_0] = 0;
    d->vertex_triangles[vt0_1] = t0;
    d->vertex_triangle_index[vt0_1] = 1;
    d->vertex_triangles[vt0_2] = t1;
    d->vertex_triangle_index[vt0_2] = 0;
    d->vertex_triangles[vt1_0] = t2;
    d->vertex_triangle_index[vt1_0] = 1;
    d->vertex_triangles[v] = t0;
    d->vertex_triangle_index[v] = 2;

    /* add the new triangles to the queue for checking */
    delaunay_enqueue(d, t0);
    delaunay_enqueue(d, t1);
    delaunay_enqueue(d, t2);
    delaunay_enqueue(d, t3);

    /* touch the last triangle; if nothing else happens, this triangle will be
       the initial guess for the next insertion */
    d->last_triangle = t3;
  }

  /* now process all triangles in the test queue */
  delaunay_check_triangles(d);
  delaunay_log("Post vertex %i check", v);
  /* perform a consistency test if enabled */
  delaunay_check_tessellation(d);
  return 0;
}

/**
 * @brief Add a local (non ghost) vertex at the given index.
 * @param d Delaunay tessellation
 * @param x, y, z Position of vertex
 * @param idx Index to of the corresponding particle
 */
inline static int delaunay_add_vertex(struct delaunay* d, double x, double y,
                                      double z, int idx) {

  int v = delaunay_new_vertex(d, x, y, idx);
  if (delaunay_finalize_vertex(d, v, x, y) == -1) {
    error("Local vertices cannot be added twice!");
  }
  return v;
}

/**
 * @brief Add a new (ghost) vertex.
 * @param d Delaunay tessellation
 * @param x, y, z Position of vertex
 * @param part_idx Index of the corresponding particle in its SWIFT cell.
 * @param ngb_v_idx Index of a close delaunay vertex.
 * particle.
 */
inline static void delaunay_add_ghost_vertex(struct delaunay* d, double x,
                                             double y, double z, int cell_sid,
                                             int part_idx, int ngb_v_idx) {

  /* Update last triangle to be a triangle connected to the neighbouring
   * vertex. */
#ifdef SWIFT_DEBUG_CHECKS
  if (ngb_v_idx < d->vertex_start || ngb_v_idx >= d->vertex_end)
    error("Invalid neighbour index passed to delaunay_add_new_vertex!");
#endif
  d->last_triangle = d->vertex_triangles[ngb_v_idx];

  /* create the new vertex */
  int v = delaunay_new_vertex(d, x, y, part_idx);
  delaunay_log("Created new vertex with index %i for part %i", v, part_idx);

  if (delaunay_finalize_vertex(d, v, x, y) == -1)
    error("Trying to add the same vertex more than once to this delaunay!");

  /* vertex is new: add neighbour information for bookkeeping */
  if (d->ghost_index == d->ghost_size) {
    d->ghost_size <<= 1;
    d->ghost_cell_sids = (int*)swift_realloc("delaunay", d->ghost_cell_sids,
                                             d->ghost_size * sizeof(int));
  }
  delaunay_assert(d->ghost_index == v - d->vertex_end);
  /* Store info */
  d->ghost_cell_sids[d->ghost_index] = cell_sid;
  ++d->ghost_index;
}

/*! @brief stores the actual coordinates of the vertex at given by idx in out */
inline static void delaunay_get_vertex_at(const struct delaunay* d, int idx,
                                          double* out /*ret*/) {
  out[0] = d->vertices[2 * idx];
  out[1] = d->vertices[2 * idx + 1];
}

inline static void delaunay_write_tessellation(
    const struct delaunay* restrict d, FILE* file, size_t* offset) {

  fprintf(file, "#VertexEnd\t%lu\tNeighbourOffset\t%lu\tnTriangles\t%d\n",
          *offset + d->vertex_end, *offset + d->vertex_end,
          d->triangle_index - 3);

  for (int i = 0; i < d->vertex_index; ++i) {
    fprintf(file, "V\t%lu\t%g\t%g\n", *offset + i, d->vertices[2 * i],
            d->vertices[2 * i + 1]);
  }
  for (int i = 3; i < d->triangle_index; ++i) {
    fprintf(file, "T\t%lu\t%lu\t%lu\n", *offset + d->triangles[i].vertices[0],
            *offset + d->triangles[i].vertices[1],
            *offset + d->triangles[i].vertices[2]);
  }
  *offset += d->vertex_index;
}

/**
 * @brief Output the tessellation to a text file with the given name.
 *
 * The tessellation is output as follows:
 *  1. All the vertices are output in order, as individual lines with the
 *     format "V\tindex\tx\ty\n", where x and y are the original indices of
 *     the vertex in box coordinates.
 *  2. All triangles are output as individual lines with the format
 *     "T\tindex\tindex\tindex\n", where the 3 index values correspond to the
 *     indices of the 3 vertices of that triangle. Since all vertices are
 *     output, these indices should match the order written by 1.
 *
 * @param d Delaunay tessellation (read-only).
 * @param file_name Name of the output file.
 */
inline static void delaunay_print_tessellation(
    const struct delaunay* restrict d, const char* file_name) {

  FILE* file = fopen(file_name, "w");

  size_t offset = 0;
  delaunay_write_tessellation(d, file, &offset);

  fclose(file);
}

#endif  // SWIFTSIM_SHADOWSWIFT_DELAUNAY_2D_H
