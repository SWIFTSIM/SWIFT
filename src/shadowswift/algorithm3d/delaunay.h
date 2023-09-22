//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_DELAUNAY_3D_H
#define SWIFTSIM_SHADOWSWIFT_DELAUNAY_3D_H

#include "../queues.h"
#include "./geometry.h"
#include "tetrahedron.h"

struct delaunay {
  /*! Activity flag, useful for debugging. */
  int active;

  /*! @brief Anchor of the simulation volume. */
  double anchor[3];

  /*! @brief Side length of the simulation volume. */
  double side;

  /*! @brief Inverse side length of the simulation volume. */
  double inverse_side;

  /*! @brief Vertex positions, rescaled to the range 1-2. Kept here in case we
   *  want to adopt hybrid geometrical checks (floating point checks when safe,
   *  integer checks when there is a risk of numerical error leading to the
   *  wrong result) to speed things up. */
  double* rescaled_vertices;

  /*! @brief Integer vertex_indices. These are the vertex coordinates that are
   *  actually used during the incremental construction. */
  unsigned long int* integer_vertices;

  /*! @brief Vertex-tetrahedron connections. For every vertex in the
   * tessellation, this array stores the index of a tetrahedron that contains
   * this vertex (usually one of the last tetrahedra that was constructed with
   * this vertex as a member vertex). This array is not required for the
   * incremental construction algorithm itself, but is indispensable for the
   * conversion from Delaunay tessellation to Voronoi grid, since it links each
   * input vertex to one of the tetrahedra it is connected to. Without this
   * array, we would have no efficient way of finding a tetrahedron that
   * contains a given vertex. */
  int* vertex_tetrahedron_links;

  /*! @brief Vertex-tetrahedron connection indices. For every vertex-tetrahedron
   * pair stored in vertex_tetrahedron_links, this array contains the index of
   * the vertex within the vertex list of the tetrahedron. This saves us from
   * having to loop through the tetrahedron vertex_indices during Voronoi grid
   * construction. */
  int* vertex_tetrahedron_index;

  /*! @brief Index of the corresponding particles in their respective cells. */
  int* vertex_part_idx;

  /*! @brief Next available index within the vertex array. Corresponds to the
   *  actual size of the vertex array. */
  int vertex_index;

  /*! @brief Current size of the vertex array in memory. If vertex_size matches
   *  n_vertices, the memory buffer is full and needs to be expanded. */
  int vertex_size;

  /*! @brief Begin index of the normal vertex_indices. This skips the 4
   * auxiliary vertex_indices required for the incremental construction
   * algorithm. */
  int vertex_start;

  /*! @brief End index of the normal vertex_indices. This variable is set by
   * calling delaunay_consolidate() and contains the offset of the ghost
   * vertex_indices within the vertex array. */
  int vertex_end;

  /*! @brief Tetrahedra that make up the tessellation. */
  struct tetrahedron* tetrahedra;

  /*! @brief Next available index within the tetrahedron array. Corresponds to
   * the actual size of the tetrahedron array. */
  int tetrahedra_index;

  /*! @brief Current size of the tetrahedron array in memory. If
   * tetrahedra_size matches tetrahedra_index, the memory buffer is full and
   * needs to be expanded. */
  int tetrahedra_size;

  /*! @brief Index of the last tetrahedron that was created or modified. Used as
   * initial guess for the tetrahedron that contains the next vertex that will
   * be added. If vertex_indices are added in some sensible order (e.g. in
   * Peano-Hilbert curve order) then this will greatly speed up the algorithm.
   */
  int last_tetrahedron;

  /*! @brief Lifo queue of tetrahedra that need checking during the incremental
   *  construction algorithm. After a new vertex has been added, all new
   *  tetrahedra are added to this queue and are tested to see if the
   *  Delaunay criterion (empty circumcircles) still holds. New tetrahedra
   *  created when flipping invalid tetrahedra are also added to the queue. */
  struct int_lifo_queue tetrahedra_to_check;

  /*! @brief Lifo queue of free spots in the tetrahedra array. Sometimes 3
   * tetrahedra can be split into 2 new ones. This leaves a free spot in the
   * array.
   */
  struct int_lifo_queue free_tetrahedron_indices;

  /*! @brief Array of tetrahedra containing the current vertex */
  struct int_lifo_queue tetrahedra_containing_vertex;

  /*! @brief Geometry variables. Auxiliary variables used by the exact integer
   *  geometry3d tests that need to be stored in between tests, since allocating
   *  and deallocating them for every test is too expensive. */
  struct geometry3d geometry;

  /*! @brief ghost particles cell sids keeping track of which cells
   *  contains a specific ghost vertex. */
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

/* forward declarations */
inline static void delaunay_reset(struct delaunay* restrict d,
                                  const double* cell_loc,
                                  const double* cell_width, int vertex_size);
inline static void delaunay_check_tessellation(struct delaunay* d);
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x,
                                      double y, double z, int idx);
inline static int delaunay_finalize_vertex(struct delaunay* restrict d, int v);
inline static int delaunay_new_tetrahedron(struct delaunay* restrict d);
inline static void delaunay_init_tetrahedron(struct delaunay* d, int t, int v0,
                                             int v1, int v2, int v3);
inline static int get_next_tetrahedron_to_check(struct delaunay* restrict d);
inline static int delaunay_find_tetrahedra_containing_vertex(struct delaunay* d,
                                                             int v);
inline static void delaunay_one_to_four_flip(struct delaunay* d, int v, int t);
inline static void delaunay_two_to_six_flip(struct delaunay* d, int v,
                                            const int* t);
inline static void delaunay_n_to_2n_flip(struct delaunay* d, int v,
                                         const int* t, int n);
inline static void delaunay_finalize_tetrahedron(struct delaunay* d, int t);
inline static void delaunay_check_tetrahedra(struct delaunay* d, int v);
inline static int delaunay_check_tetrahedron(struct delaunay* d, int t, int v);
inline static int positive_permutation(int a, int b, int c, int d);
inline static double delaunay_get_radius(struct delaunay* restrict d, int t);
inline static int delaunay_test_orientation(struct delaunay* restrict d, int v0,
                                            int v1, int v2, int v3);
inline static int delaunay_vertex_is_valid(struct delaunay* restrict d, int v);
inline static void delaunay_get_vertex_at(const struct delaunay* d, int idx,
                                          double* out);
inline static int delaunay_choose_2(int a, int b,
                                    const double* restrict centroid,
                                    const double* a0, const double* a1,
                                    const double* a2, const double* b0,
                                    const double* b1, const double* b2,
                                    const double* restrict v);
inline static int delaunay_choose_3(int a, int b, int c,
                                    const double* restrict centroid,
                                    const double* a0, const double* a1,
                                    const double* a2, const double* b0,
                                    const double* b1, const double* b2,
                                    const double* c0, const double* c1,
                                    const double* c2, const double* restrict v);
inline static int delaunay_choose_random_3(int a, int b, int c);

/**
 * @brief Initialize the Delaunay tessellation.
 *
 * This function allocates memory for all arrays that make up the tessellation
 * and initializes the variables used for bookkeeping.
 *
 * It then sets up a large tetrahedron that contains the entire simulation box
 * and additional buffer space to deal with boundary ghost vertex_indices, and 4
 * additional dummy tetrahedron that provide valid neighbours for the 4 sides of
 * this tetrahedron (these dummy tetrahedra themselves have an invalid tip
 * vertex and are therefore simply placeholders).
 *
 * @param d Delaunay tessellation.
 * @param hs Spatial extents of the simulation box.
 * @param vertex_size Initial size of the vertex array.
 * @param tetrahedra_size Initial size of the tetrahedra array.
 */
inline static struct delaunay* delaunay_malloc(const double* cell_loc,
                                               const double* cell_width,
                                               int vertex_size) {

  /* Don't bother setting up a Delaunay tessellation for empty cells */
  if (vertex_size == 0) {
    return NULL;
  }

  struct delaunay* d = malloc(sizeof(struct delaunay));

  d->active = 1;

  /* allocate memory for the vertex arrays */
  d->vertex_size = vertex_size;
  d->rescaled_vertices =
      (double*)swift_malloc("delaunay", vertex_size * 3 * sizeof(double));
  d->integer_vertices = (unsigned long int*)swift_malloc(
      "delaunay", vertex_size * 3 * sizeof(unsigned long int));
  d->vertex_tetrahedron_links =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_tetrahedron_index =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_part_idx =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));

  /* Allocate memory for the tetrahedra array */
  d->tetrahedra_size = 10 * vertex_size;
  d->tetrahedra = (struct tetrahedron*)swift_malloc(
      "delaunay", d->tetrahedra_size * sizeof(struct tetrahedron));

  /* Allocate queues */
  int_lifo_queue_init(&d->tetrahedra_containing_vertex, 10);
  int_lifo_queue_init(&d->tetrahedra_to_check, 10);
  int_lifo_queue_init(&d->free_tetrahedron_indices, 10);

  /* Allocate the array with the cell information of the ghost particles */
  d->ghost_cell_sids =
      (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->ghost_size = vertex_size;

  /* initialise the structure used to perform exact geometrical tests */
  geometry3d_init(&d->geometry);

  /* Initialize the Delaunay tesselation */
  delaunay_reset(d, cell_loc, cell_width, vertex_size);
  return d;
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

  if (d->active != 1) {
    error("Delaunay tessellation corruption!");
  }

  /* reset the vertex and tetrahedra array indices. */
  d->vertex_index = 0;
  d->tetrahedra_index = 0;

  /* Reset all the queues */
  int_lifo_queue_reset(&d->tetrahedra_to_check);
  int_lifo_queue_reset(&d->free_tetrahedron_indices);
  int_lifo_queue_reset(&d->tetrahedra_containing_vertex);

  /* Reset the sid mask.
   * Sid=13 does not correspond to a face and is always set to 1.
   * We only set the sid's corresponding to the cardinal directions to 0
   * (only faces perpendicular to one of the axes can be boundary faces). */
  d->sid_is_inside_face_mask = DEFAULT_SID_MASK;

  /* determine the size of a box large enough to accommodate the entire
   * simulation volume and all possible ghost vertex_indices required to deal
   * with boundaries. Note that we convert the generally rectangular box to a
   * square. */
  /* We add an extra layer of padding of 1 times cell_width around our box to
   * compensate for particle movements (up to 1 cell side before a rebuild is
   * triggered). */
  double box_anchor[3] = {cell_loc[0] - 2 * cell_width[0],
                          cell_loc[1] - 2 * cell_width[1],
                          cell_loc[2] - 2 * cell_width[2]};
  /* Notice we have to take box_side rather large, because we want to fit the
   * cell and all neighbouring cells inside the first tetrahedron. This comes at
   * a loss of precision in the integer arithmetic, though... A better solution
   * would possibly be to start from 5 tetrahedra forming a cube (box_side would
   * have to be 5 in that case). */
  double box_side = fmax(cell_width[0], cell_width[1]);
  box_side = 15. * fmax(box_side, cell_width[2]);
  /* the 1.e-13 makes sure converted values are in the range [1, 2[ instead of
   * [1,2] (unlike Springel, 2010) */
  d->inverse_side = (1. - 1.e-13) / box_side;
  d->side = box_side;

  /* store the anchor and inverse side_length for the conversion from box
     coordinates to rescaled (integer) coordinates */
  d->anchor[0] = box_anchor[0];
  d->anchor[1] = box_anchor[1];
  d->anchor[2] = box_anchor[2];

  /* set up vertex_indices for large initial tetrahedron */
  int v0 =
      delaunay_new_vertex(d, box_anchor[0], box_anchor[1], box_anchor[2], -1);
  int v1 = delaunay_new_vertex(d, box_anchor[0] + box_side, box_anchor[1],
                               box_anchor[2], -1);
  int v2 = delaunay_new_vertex(d, box_anchor[0], box_anchor[1] + box_side,
                               box_anchor[2], -1);
  int v3 = delaunay_new_vertex(d, box_anchor[0], box_anchor[1],
                               box_anchor[2] + box_side, -1);

  /* Initialise the indices indicating where the local vertices start and end.*/
  d->vertex_start = d->vertex_index;
  d->vertex_end = 0;

  /* Reset the ghost index. */
  d->ghost_index = 0;

  /* Create initial large tetrahedron and 4 dummy neighbours */
  int dummy0 = delaunay_new_tetrahedron(d); /* opposite of v0 */
  int dummy1 = delaunay_new_tetrahedron(d); /* opposite of v1 */
  int dummy2 = delaunay_new_tetrahedron(d); /* opposite of v2 */
  int dummy3 = delaunay_new_tetrahedron(d); /* opposite of v3 */
  int first_tetrahedron = delaunay_new_tetrahedron(d);
  delaunay_log(
      "Creating dummy tetrahedron at %i with vertex_indices: %i %i %i %i",
      dummy0, v1, v2, v3, -1);
  tetrahedron_init(&d->tetrahedra[dummy0], v1, v2, v3, -1);
  delaunay_log(
      "Creating dummy tetrahedron at %i with vertex_indices: %i %i %i %i",
      dummy1, v2, v0, v3, -1);
  tetrahedron_init(&d->tetrahedra[dummy1], v2, v0, v3, -1);
  delaunay_log(
      "Creating dummy tetrahedron at %i with vertex_indices: %i %i %i %i",
      dummy2, v3, v0, v1, -1);
  tetrahedron_init(&d->tetrahedra[dummy2], v3, v0, v1, -1);
  delaunay_log(
      "Creating dummy tetrahedron at %i with vertex_indices: %i %i %i %i",
      dummy3, v0, v2, v1, -1);
  tetrahedron_init(&d->tetrahedra[dummy3], v0, v2, v1, -1);
  delaunay_init_tetrahedron(d, first_tetrahedron, v0, v1, v2, v3);
  delaunay_finalize_tetrahedron(d, first_tetrahedron);

  /* Setup neighbour relations */
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy0], 3, first_tetrahedron, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy1], 3, first_tetrahedron, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy2], 3, first_tetrahedron, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[dummy3], 3, first_tetrahedron, 3);
  tetrahedron_swap_neighbours(&d->tetrahedra[first_tetrahedron], dummy0, dummy1,
                              dummy2, dummy3, 3, 3, 3, 3);

  /* Perform sanity checks */
  delaunay_check_tessellation(d);
  delaunay_log("Passed post init or reset check");
}

inline static void delaunay_destroy(struct delaunay* restrict d) {
  if (d->active != 1) {
    error("Delaunay tessellation corruption!");
  }
  d->active = 0;

  swift_free("delaunay", d->rescaled_vertices);
  swift_free("delaunay", d->integer_vertices);
  swift_free("delaunay", d->vertex_tetrahedron_links);
  swift_free("delaunay", d->vertex_tetrahedron_index);
  swift_free("delaunay", d->vertex_part_idx);
  swift_free("delaunay", d->tetrahedra);
  int_lifo_queue_destroy(&d->tetrahedra_to_check);
  int_lifo_queue_destroy(&d->free_tetrahedron_indices);
  int_lifo_queue_destroy(&d->tetrahedra_containing_vertex);
  geometry3d_destroy(&d->geometry);
  swift_free("delaunay", d->ghost_cell_sids);

  d->rescaled_vertices = NULL;
  d->integer_vertices = NULL;
  d->vertex_tetrahedron_links = NULL;
  d->vertex_tetrahedron_index = NULL;
  d->vertex_part_idx = NULL;
  d->tetrahedra = NULL;
  d->ghost_cell_sids = NULL;

  bzero(d->anchor, 3 * sizeof(double));
  d->side = 0.;
  d->inverse_side = 0.;
  d->vertex_size = 0;
  d->vertex_index = -1;
  d->vertex_start = -1;
  d->vertex_end = -1;
  d->tetrahedra_size = 0;
  d->tetrahedra_index = -1;
  d->last_tetrahedron = -1;
  d->ghost_size = 0;
  d->ghost_index = -1;
  d->sid_is_inside_face_mask = 0;

  /* Free delaunay struct itself */
  free(d);
}

inline static int delaunay_new_tetrahedron(struct delaunay* restrict d) {
  /* check whether there is a free spot somewhere in the array */
  if (!int_lifo_queue_is_empty(&d->free_tetrahedron_indices)) {
    return int_lifo_queue_pop(&d->free_tetrahedron_indices);
  }
  /* Else: check that we still have space for tetrahedrons available */
  if (d->tetrahedra_index == d->tetrahedra_size) {
    d->tetrahedra_size <<= 1;
    d->tetrahedra = (struct tetrahedron*)swift_realloc(
        "delaunay", d->tetrahedra,
        d->tetrahedra_size * sizeof(struct tetrahedron));
  }
  /* return and then increase */
  return d->tetrahedra_index++;
}

/**
 * @brief Utility method to initialize a tetrahedron.
 *
 * If Checks are enabled, this will also check the orientation of the
 * tetrahedron.
 *
 * @param d Delaunay tesselation
 * @param t Index to initialize tetrahedron at
 * @param v0, v1, v2, v3 Indices of the vertex_indices of the tetrahedron
 */
inline static void delaunay_init_tetrahedron(struct delaunay* d, int t, int v0,
                                             int v1, int v2, int v3) {
  delaunay_log(
      "Initializing tetrahedron at %i with vertex_indices: %i %i %i %i", t, v0,
      v1, v2, v3);
#ifdef SWIFT_DEBUG_CHECKS
  const int test = delaunay_test_orientation(d, v0, v1, v2, v3);
  if (test >= 0) {
    fprintf(stderr, "Initializing tetrahedron with incorrect orientation!\n");
    fprintf(stderr, "\tTetrahedron: %i\n\tVertices: %i %i %i %i\n", t, v0, v1,
            v2, v3);
    abort();
  }
#endif
  tetrahedron_init(&d->tetrahedra[t], v0, v1, v2, v3);

  /* Update vertex-tetrahedron links */
  d->vertex_tetrahedron_links[v0] = t;
  d->vertex_tetrahedron_index[v0] = 0;
  d->vertex_tetrahedron_links[v1] = t;
  d->vertex_tetrahedron_index[v1] = 1;
  d->vertex_tetrahedron_links[v2] = t;
  d->vertex_tetrahedron_index[v2] = 2;
  d->vertex_tetrahedron_links[v3] = t;
  d->vertex_tetrahedron_index[v3] = 3;

  /* Touch the last initialized tetrahedron, This will be our next guess upon
   * insertion. */
  d->last_tetrahedron = t;
}

/**
 * @brief Add a new vertex with the given coordinates.
 *
 * This function first makes sure there is sufficient memory to store the
 * vertex and all its properties. It then initializes the vertex.
 *
 * @param d Delaunay tessellation.
 * @param x Horizontal coordinate of the vertex.
 * @param y Vertical coordinate of the vertex.
 * @param z Z position of the vertex.
 * @param idx Index of the corresponding part in its SWIFT cell.
 * @return Index of the new vertex within the vertex array.
 */
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x,
                                      double y, double z, int idx) {
  delaunay_log("Adding new vertex at %i with coordinates: %g %g %g",
               d->vertex_index, x, y, z);
  /* check the size of the vertex arrays against the allocated memory size */
  if (d->vertex_index == d->vertex_size) {
    /* dynamically grow the size of the arrays with a factor 2 */
    d->vertex_size <<= 1;
    d->rescaled_vertices = (double*)swift_realloc(
        "delaunay", d->rescaled_vertices, d->vertex_size * 3 * sizeof(double));
    d->integer_vertices = (unsigned long int*)swift_realloc(
        "delaunay", d->integer_vertices,
        d->vertex_size * 3 * sizeof(unsigned long int));
    d->vertex_tetrahedron_links = (int*)swift_realloc(
        "delaunay", d->vertex_tetrahedron_links, d->vertex_size * sizeof(int));
    d->vertex_tetrahedron_index = (int*)swift_realloc(
        "delaunay", d->vertex_tetrahedron_index, d->vertex_size * sizeof(int));
    d->vertex_part_idx = (int*)swift_realloc("delaunay", d->vertex_part_idx,
                                             d->vertex_size * sizeof(int));
  }

  /* compute the rescaled coordinates. We do this because floating point values
     in the range [1,2[ all have the same exponent (0), which guarantees that
     their mantissas form a linear sequence */
  double rescaledX = 1. + (x - d->anchor[0]) * d->inverse_side;
  double rescaledY = 1. + (y - d->anchor[1]) * d->inverse_side;
  double rescaledZ = 1. + (z - d->anchor[2]) * d->inverse_side;

  delaunay_assert(rescaledX >= 1.);
  delaunay_assert(rescaledX < 2.);
  delaunay_assert(rescaledY >= 1.);
  delaunay_assert(rescaledY < 2.);
  delaunay_assert(rescaledZ >= 1.);
  delaunay_assert(rescaledZ < 2.);

  /* store a copy of the rescaled coordinates to apply non-exact tests */
  d->rescaled_vertices[3 * d->vertex_index] = rescaledX;
  d->rescaled_vertices[3 * d->vertex_index + 1] = rescaledY;
  d->rescaled_vertices[3 * d->vertex_index + 2] = rescaledZ;

  /* convert the rescaled coordinates to integer coordinates and store these */
  d->integer_vertices[3 * d->vertex_index] = delaunay_double_to_int(rescaledX);
  d->integer_vertices[3 * d->vertex_index + 1] =
      delaunay_double_to_int(rescaledY);
  d->integer_vertices[3 * d->vertex_index + 2] =
      delaunay_double_to_int(rescaledZ);

  /* Set the particle index */
  d->vertex_part_idx[d->vertex_index] = idx;

  /* Initialise the variables that keep track of the link between vertex_indices
   * and tetrahedra. We use negative values so that we can later detect missing
   * links. */
  d->vertex_tetrahedron_links[d->vertex_index] = -1;
  d->vertex_tetrahedron_index[d->vertex_index] = -1;

  /* return the vertex index and then increase it by 1.
     After this operation, n_vertices will correspond to the size of the
     vertex arrays and is also the index of the next vertex that will be
     created. */
  return d->vertex_index++;
}

/**
 * @brief Add a local (non ghost) vertex at the given index.
 * @param d Delaunay tessellation
 * @param x, y, z Position of vertex
 * @param idx Index to of the corresponding particle
 */
inline static int delaunay_add_vertex(struct delaunay* d, double x, double y,
                                      double z, int idx) {
  delaunay_assert(d->active == 1);
  delaunay_log("Adding vertex at %i with coordinates: %g %g %g", idx, x, y, z);

  int v = delaunay_new_vertex(d, x, y, z, idx);
  if (delaunay_finalize_vertex(d, v) == -1) {
    error("Vertices cannot be added twice!");
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
  delaunay_assert(d->active == 1);

  /* Update last tetrahedron to be a tetrahedron connected to the neighbouring
   * vertex. */
#ifdef SWIFT_DEBUG_CHECKS
  if (ngb_v_idx < d->vertex_start || d->vertex_end <= ngb_v_idx)
    error("Invalid neighbour index passed to delaunay_add_new_vertex!");
#endif
  d->last_tetrahedron = d->vertex_tetrahedron_links[ngb_v_idx];

  /* Create the new vertex */
  int v = delaunay_new_vertex(d, x, y, z, part_idx);

  if (delaunay_finalize_vertex(d, v) == -1)
    error("Trying to add the same vertex more than once to this delaunay!");

  /* Ghost vertex: add neighbour information */
  if (d->ghost_index == d->ghost_size) {
    d->ghost_size <<= 1;
    d->ghost_cell_sids = (int*)swift_realloc("delaunay", d->ghost_cell_sids,
                                             d->ghost_size * sizeof(int));
  }
  /* Store info */
  delaunay_assert(d->ghost_index == v - d->vertex_end);
  d->ghost_cell_sids[d->ghost_index] = cell_sid;
  d->ghost_index++;
}

/**
 * @brief Finalize adding a new vertex to the tessellation.
 *
 * After this function, the Delaunay tesselation is restored to a valid state.
 * This function locates the tetrahedron in the current tessellation that
 * contains the new vertex. Depending on the case (see below) new tetrahedra are
 * added to the tessellation and some are removed.
 *
 * Cases:
 * 1. Point is fully inside a tetrahedron. In this case, this tetrahedron will
 *    be replaced by 4 new tetrahedra.
 * 2. Point is on face between two tetrahedra. In this case, these two will be
 *    replaced by 6 new ones.
 * 3. Point is on edge of N tetrahedra. These N tetrahedra will be replaced by
 *    2N new ones.
 *
 * @param d Delaunay tessellation.
 * @param v Index of new vertex
 * @return The number of tetrahedra initially containing the new vertex
 */
inline static int delaunay_finalize_vertex(struct delaunay* restrict d, int v) {
#ifdef DELAUNAY_DO_ASSERTIONS
  /* Check that the new vertex falls in the bounding box */
  const unsigned long* vl = &d->integer_vertices[3 * v];
  const unsigned long* d0l = &d->integer_vertices[0];
  const unsigned long* d1l = &d->integer_vertices[3];
  const unsigned long* d2l = &d->integer_vertices[6];
  const unsigned long* d3l = &d->integer_vertices[9];

  const double* vd = &d->rescaled_vertices[3 * v];
  const double* d0d = &d->rescaled_vertices[0];
  const double* d1d = &d->rescaled_vertices[3];
  const double* d2d = &d->rescaled_vertices[6];
  const double* d3d = &d->rescaled_vertices[9];

  delaunay_assert(geometry3d_orient_adaptive(&d->geometry, d0l, d1l, d2l, vl,
                                             d0d, d1d, d2d, vd) &&
                  geometry3d_orient_adaptive(&d->geometry, d0l, d2l, d3l, vl,
                                             d0d, d2d, d3d, vd) &&
                  geometry3d_orient_adaptive(&d->geometry, d0l, d3l, d1l, vl,
                                             d0d, d3d, d1d, vd) &&
                  geometry3d_orient_adaptive(&d->geometry, d1l, d3l, d2l, vl,
                                             d1d, d3d, d2d, vd));
#endif

  int number_of_tetrahedra = delaunay_find_tetrahedra_containing_vertex(d, v);

  if (number_of_tetrahedra == -1) {
    /* Vertex already exists! */
    return -1;
  } else if (number_of_tetrahedra == 1) {
    /* normal case: split 'd->tetrahedra_containing_vertex[0]' into 4 new
     * tetrahedra */
    delaunay_log("Vertex %i lies fully inside tetrahedron %i", v,
                 d->tetrahedra_containing_vertex[0]);
    delaunay_one_to_four_flip(d, v, d->tetrahedra_containing_vertex.values[0]);
  } else if (number_of_tetrahedra == 2) {
    /* point on face: replace the 2 tetrahedra with 6 new ones */
    delaunay_log("Vertex %i on the face between tetrahedra %i and %i", v,
                 d->tetrahedra_containing_vertex[0],
                 d->tetrahedra_containing_vertex[1]);
    delaunay_two_to_six_flip(d, v, d->tetrahedra_containing_vertex.values);
  } else if (number_of_tetrahedra > 2) {
    /* point on edge: replace the N tetrahedra with 2N new ones */
    delaunay_log(
        "Vertex %i lies on the edge shared by tetrahedra %i, %i and %i", v,
        d->tetrahedra_containing_vertex[0], d->tetrahedra_containing_vertex[1],
        d->tetrahedra_containing_vertex[number_of_tetrahedra - 1]);
    delaunay_n_to_2n_flip(d, v, d->tetrahedra_containing_vertex.values,
                          number_of_tetrahedra);
  } else {
    fprintf(stderr, "Unknown case of number of tetrahedra: %i!\n",
            number_of_tetrahedra);
    abort();
  }

  /* Now check all tetrahedra in de queue */
  delaunay_check_tetrahedra(d, v);

  /* perform sanity checks if enabled */
  delaunay_check_tessellation(d);
  delaunay_log("Passed checks after inserting vertex %i", v);

  return number_of_tetrahedra;
}

/**
 * @brief Find tetrahedra containing the given vertex
 * The tetrahedra are stored in d->tetrahedra_containing_vertex
 *
 * @param d Delaunay tessellation
 * @param v The vertex
 * @return The number of tetrahedra containing the given vertex.
 */
inline static int delaunay_find_tetrahedra_containing_vertex(
    struct delaunay* restrict d, const int v) {
  /* Before we do anything: reset the index in the array of tetrahedra
   * containing the current vertex */
  int_lifo_queue_reset(&d->tetrahedra_containing_vertex);

  /* Get integer and rescaled vertex coordinates */
  const unsigned long* el = &d->integer_vertices[3 * v];
  const double* ed = &d->rescaled_vertices[3 * v];

  /* Get the last tetrahedron index */
  int tetrahedron_idx = d->last_tetrahedron;

  while (int_lifo_queue_is_empty(&d->tetrahedra_containing_vertex)) {
    delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
    const struct tetrahedron* tetrahedron = &d->tetrahedra[tetrahedron_idx];
    const int v0 = tetrahedron->vertices[0];
    const int v1 = tetrahedron->vertices[1];
    const int v2 = tetrahedron->vertices[2];
    const int v3 = tetrahedron->vertices[3];

    /* Get pointers to the coordinates of the vertices */
    const unsigned long* al = &d->integer_vertices[3 * v0];
    const unsigned long* bl = &d->integer_vertices[3 * v1];
    const unsigned long* cl = &d->integer_vertices[3 * v2];
    const unsigned long* dl = &d->integer_vertices[3 * v3];

    const double* ad = &d->rescaled_vertices[3 * v0];
    const double* bd = &d->rescaled_vertices[3 * v1];
    const double* cd = &d->rescaled_vertices[3 * v2];
    const double* dd = &d->rescaled_vertices[3 * v3];

#ifdef DELAUNAY_CHECKS
    /* made sure the tetrahedron is correctly oriented */
    if (geometry3d_orient_adaptive(&d->geometry, al, bl, cl, dl, ad, bd, cd,
                                   dd) >= 0) {
      fprintf(stderr, "Incorrect orientation for tetrahedron %i!",
              tetrahedron_idx);
      abort();
    }
#endif

    int non_axis_v_idx[4];
    int next_tetrahedron_idx = -1;

    /* Check whether the point is inside or outside all four faces */
    int tests[4];
#ifdef DELAUNAY_3D_RANDOM_SUP_TET
    int test_flags = 0;
#ifdef DELAUNAY_3D_HAND_VEC
    geometry3d_orient_4(&d->geometry, al, bl, cl, dl, el, ad, bd, cd, dd, ed,
                        tests);
#else
    tests[0] = geometry3d_orient_adaptive(&d->geometry, bl, dl, cl, el, bd, dd,
                                          cd, ed);
    tests[1] = geometry3d_orient_adaptive(&d->geometry, al, cl, dl, el, ad, cd,
                                          dd, ed);
    tests[2] = geometry3d_orient_adaptive(&d->geometry, al, dl, bl, el, ad, dd,
                                          bd, ed);
    tests[3] = geometry3d_orient_adaptive(&d->geometry, al, bl, cl, el, ad, bd,
                                          cd, ed);
#endif
    if (tests[0] > 0) {
      test_flags |= 1;
    }

    if (tests[1] > 0) {
      test_flags |= 1 << 1;
    }

    if (tests[2] > 0) {
      test_flags |= 1 << 2;
    }

    if (tests[3] > 0) {
      test_flags |= 1 << 3;
    }
    switch (test_flags) {
      case 1:
        /* Orientation test BDCE > 0 */
        next_tetrahedron_idx = tetrahedron->neighbours[0];
        break;
      case 2:
        /* Orientation test ACDE > 0 */
        next_tetrahedron_idx = tetrahedron->neighbours[1];
        break;
      case 4:
        /* Orientation test ADBE > 0 */
        next_tetrahedron_idx = tetrahedron->neighbours[2];
        break;
      case 8:
        /* Orientation test ABCE > 0 */
        next_tetrahedron_idx = tetrahedron->neighbours[3];
        break;
      case 3:
        /* Orientation test BDCE and ACDE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[0], tetrahedron->neighbours[1],
            tetrahedron->centroid, bd, dd, cd, ad, cd, dd, ed);
        break;
      case 5:
        /* Orientation test BDCE and ADBE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[0], tetrahedron->neighbours[2],
            tetrahedron->centroid, bd, dd, cd, ad, dd, bd, ed);
        break;
      case 6:
        /* Orientation test ACDE and ADBE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[1], tetrahedron->neighbours[2],
            tetrahedron->centroid, ad, cd, dd, ad, dd, bd, ed);
        break;
      case 9:
        /* Orientation test BDCE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[0], tetrahedron->neighbours[3],
            tetrahedron->centroid, bd, dd, cd, ad, bd, cd, ed);
        break;
      case 10:
        /* Orientation test ACDE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[1], tetrahedron->neighbours[3],
            tetrahedron->centroid, ad, cd, dd, ad, bd, cd, ed);
        break;
      case 12:
        /* Orientation test ADBE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_2(
            tetrahedron->neighbours[2], tetrahedron->neighbours[3],
            tetrahedron->centroid, ad, dd, bd, ad, bd, cd, ed);
        break;
      case 7:
        /* Orientation test BDCE and ACDE and ADBE > 0 */
        next_tetrahedron_idx = delaunay_choose_random_3(
            tetrahedron->neighbours[0], tetrahedron->neighbours[1],
            tetrahedron->neighbours[2]);
        break;
      case 11:
        /* Orientation test BDCE and ACDE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_random_3(
            tetrahedron->neighbours[0], tetrahedron->neighbours[1],
            tetrahedron->neighbours[3]);
        break;
      case 13:
        /* Orientation test BDCE and ADBE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_random_3(
            tetrahedron->neighbours[0], tetrahedron->neighbours[2],
            tetrahedron->neighbours[3]);
        break;
      case 14:
        /* Orientation test ACDE and ADBE and ABCE > 0 */
        next_tetrahedron_idx = delaunay_choose_random_3(
            tetrahedron->neighbours[1], tetrahedron->neighbours[2],
            tetrahedron->neighbours[3]);
        break;
      case 15:
        /* All tests > 0, which is impossible (vertex cannot lay outside 4
         * faces). */
        error("Impossible tetrahedron test results!");
        break;
      case 0:
        /* all tests <= 0: the vertex lies in the current triangle or this is a
         * degeneracy which will be resolved below. */
        break;
      default:
        error("Impossible tetrahedron test result: %n", &test_flags);
    }
#else  // DELAUNAY_3D_RANDOM_SUP_TET
    double min_dist = -1.;
    double dist;
    double centroid[3];
    geometry3d_compute_centroid_tetrahedron(ad[0], ad[1], ad[2], bd[0], bd[1],
                                            bd[2], cd[0], cd[1], cd[2], dd[0],
                                            dd[1], dd[2], centroid);

    unsigned long centroid_ul[3];
    geometry3d_compute_centroid_tetrahedron_exact(
        al[0], al[1], al[2], bl[0], bl[1], bl[2], cl[0], cl[1], cl[2], dl[0],
        dl[1], dl[2], centroid_ul);

    struct shadowswift_ray r;
    shadowswift_ray_init(&r, centroid, ed, centroid_ul, el);
#ifdef DELAUNAY_3D_TRIANGLE_INTERSECTIONS
    tests[3] = geometry3d_orient_adaptive(&d->geometry, al, bl, cl, el, ad, bd,
                                          cd, ed);
    if (tests[3] > 0) {
      /* v outside face opposite of v3 */
      if (geometry3d_ray_triangle_intersect(&d->geometry, &r, ad, bd, cd, al,
                                            bl, cl, &dist)) {
        tetrahedron_idx = tetrahedron->neighbours[3];
        delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
        continue;
      }
      if (isnan(min_dist) || dist < min_dist) {
        delaunay_assert(dist > -1e-13);
        min_dist = dist;
        next_tetrahedron_idx = tetrahedron->neighbours[3];
      }
    }
    tests[1] = geometry3d_orient_adaptive(&d->geometry, al, cl, dl, el, ad, cd,
                                          dd, ed);
    if (tests[1] > 0) {
      /* v outside face opposite of v1 */
      if (geometry3d_ray_triangle_intersect(&d->geometry, &r, ad, cd, dd, al,
                                            cl, dl, &dist)) {
        tetrahedron_idx = tetrahedron->neighbours[1];
        delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
        continue;
      }
      if (isnan(min_dist) || dist < min_dist) {
        delaunay_assert(dist > -1e-13);
        min_dist = dist;
        next_tetrahedron_idx = tetrahedron->neighbours[1];
      }
    }
    tests[2] = geometry3d_orient_adaptive(&d->geometry, al, dl, bl, el, ad, dd,
                                          bd, ed);
    if (tests[2] > 0) {
      /* v outside face opposite of v2 */
      if (geometry3d_ray_triangle_intersect(&d->geometry, &r, ad, dd, bd, al,
                                            dl, bl, &dist)) {
        tetrahedron_idx = tetrahedron->neighbours[2];
        delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
        continue;
      }
      if (isnan(min_dist) || dist < min_dist) {
        delaunay_assert(dist > -1e-13);
        min_dist = dist;
        next_tetrahedron_idx = tetrahedron->neighbours[2];
      }
    }
    tests[0] = geometry3d_orient_adaptive(&d->geometry, bl, dl, cl, el, bd, dd,
                                          cd, ed);
    if (tests[0] > 0) {
      /* v outside face opposite of v0 */
      if (geometry3d_ray_triangle_intersect(&d->geometry, &r, bd, dd, cd, bl,
                                            dl, cl, &dist)) {
        tetrahedron_idx = tetrahedron->neighbours[0];
        delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
        continue;
      }
      if (isnan(min_dist) || dist < min_dist) {
        delaunay_assert(dist > -1e-13);
        min_dist = dist;
        next_tetrahedron_idx = tetrahedron->neighbours[0];
      }
    }
#else  // DELAUNAY_3D_TRIANGLE_INTERSECTIONS
    tests[3] = geometry3d_orient_adaptive(&d->geometry, al, bl, cl, el, ad, bd,
                                          cd, ed);
    dist = geometry3d_ray_plane_intersect(&r, ad, bd, cd);
    if (tests[3] > 0 && (min_dist == -1. || dist < min_dist)) {
      delaunay_assert(dist > -1e-13);
      min_dist = dist;
      next_tetrahedron_idx = tetrahedron->neighbours[3];
    }

    tests[1] = geometry3d_orient_adaptive(&d->geometry, al, cl, dl, el, ad, cd,
                                          dd, ed);
    dist = geometry3d_ray_plane_intersect(&r, ad, cd, dd);
    if (tests[1] > 0 && (min_dist == -1. || dist < min_dist)) {
      delaunay_assert(dist > -1e-13);
      min_dist = dist;
      next_tetrahedron_idx = tetrahedron->neighbours[1];
    }

    tests[2] = geometry3d_orient_adaptive(&d->geometry, al, dl, bl, el, ad, dd,
                                          bd, ed);
    dist = geometry3d_ray_plane_intersect(&r, ad, dd, bd);
    if (tests[2] > 0 && (min_dist == -1. || dist < min_dist)) {
      delaunay_assert(dist > -1e-13);
      min_dist = dist;
      next_tetrahedron_idx = tetrahedron->neighbours[2];
    }

    tests[0] = geometry3d_orient_adaptive(&d->geometry, bl, dl, cl, el, bd, dd,
                                          cd, ed);
    if (tests[0] > 0) {
      if (min_dist == -1.) {
        /* Point inside other faces */
#ifdef DELAUNAY_DO_ASSERTIONS
        dist = geometry3d_ray_plane_intersect(&r, bd, dd, cd);
        delaunay_assert(dist > -1e-13);
#endif
        next_tetrahedron_idx = tetrahedron->neighbours[0];
      } else {
        /* Normal case: compare distance */
        dist = geometry3d_ray_plane_intersect(&r, bd, dd, cd);
        if (dist < min_dist) {
          delaunay_assert(dist > -1e-13);
          next_tetrahedron_idx = tetrahedron->neighbours[0];
        }
      }
    }
#endif  // DELAUNAY_3D_TRIANGLE_INTERSECTIONS
#endif  // DELAUNAY_3D_RANDOM_SUP_TET

    if (next_tetrahedron_idx != -1) {
      tetrahedron_idx = next_tetrahedron_idx;
      delaunay_assert(tetrahedron_idx > 3); /* No dummy tetrahedron? */
      continue;
    }

    /* Point inside tetrahedron, check for degenerate cases */
    delaunay_assert(tests[3] <= 0 && tests[1] <= 0 && tests[2] <= 0 &&
                    tests[0] <= 0);

    int n_zero_tests = 0;
    int_lifo_queue_push(&d->tetrahedra_containing_vertex, tetrahedron_idx);
    if (tests[3] == 0) {
      non_axis_v_idx[n_zero_tests] = 3;
      int_lifo_queue_push(&d->tetrahedra_containing_vertex,
                          tetrahedron->neighbours[3]);
      n_zero_tests++;
    }
    if (tests[2] == 0) {
      non_axis_v_idx[n_zero_tests] = 2;
      int_lifo_queue_push(&d->tetrahedra_containing_vertex,
                          tetrahedron->neighbours[2]);
      n_zero_tests++;
    }
    if (tests[1] == 0) {
      non_axis_v_idx[n_zero_tests] = 1;
      int_lifo_queue_push(&d->tetrahedra_containing_vertex,
                          tetrahedron->neighbours[1]);
      n_zero_tests++;
    }
    if (tests[0] == 0) {
      non_axis_v_idx[n_zero_tests] = 0;
      int_lifo_queue_push(&d->tetrahedra_containing_vertex,
                          tetrahedron->neighbours[0]);
      n_zero_tests++;
    }

    if (n_zero_tests > 2) {
      /* Impossible case, the vertex cannot simultaneously lie in this
       * tetrahedron and 3 or more of its direct neighbours, unless it coincides
       * with an already existing vertex */
      return -1;
    }
    if (n_zero_tests > 1) {
      /* Vertex on edge of tetrahedron. This edge can be shared by any number of
       * tetrahedra, of which we already know three. Find the other ones by
       * rotating around this edge. */
      const int non_axis_idx0 = non_axis_v_idx[0];
      const int non_axis_idx1 = non_axis_v_idx[1];
      int axis_idx0 = (non_axis_idx0 + 1) % 4;
      if (axis_idx0 == non_axis_idx1) {
        axis_idx0 = (axis_idx0 + 1) % 4;
      }
      const int axis_idx1 = 6 - axis_idx0 - non_axis_idx0 - non_axis_idx1;
      delaunay_assert(
          axis_idx0 != axis_idx1 && axis_idx0 != non_axis_idx0 &&
          axis_idx0 != non_axis_idx1 && axis_idx1 != non_axis_idx0 &&
          axis_idx1 != non_axis_idx1 && non_axis_idx0 != non_axis_idx1);

      /* a0 and a1 are the vertex_indices shared by all tetrahedra */
      const int a0 = tetrahedron->vertices[axis_idx0];
      const int a1 = tetrahedron->vertices[axis_idx1];

      /* We now walk around the axis and add all tetrahedra to the list of
       * tetrahedra containing v. */
      const int last_t = d->tetrahedra_containing_vertex.values[1];
      int next_t = d->tetrahedra_containing_vertex.values[2];
      int next_vertex = tetrahedron->index_in_neighbour[non_axis_idx1];

      /* We are going to add next_t and last_t back to the array of tetrahedra
       * containing v, but now with all other tetrahedra that also share the
       * edge in between, so first remove them from the queue. */
      d->tetrahedra_containing_vertex.index -= 2;
      while (next_t != last_t) {
        int_lifo_queue_push(&d->tetrahedra_containing_vertex, next_t);
        next_vertex = (next_vertex + 1) % 4;
        if (d->tetrahedra[next_t].vertices[next_vertex] == a0 ||
            d->tetrahedra[next_t].vertices[next_vertex] == a1) {
          next_vertex = (next_vertex + 1) % 4;
        }
        if (d->tetrahedra[next_t].vertices[next_vertex] == a0 ||
            d->tetrahedra[next_t].vertices[next_vertex] == a1) {
          next_vertex = (next_vertex + 1) % 4;
        }
        delaunay_assert(d->tetrahedra[next_t].vertices[next_vertex] != a0 &&
                        d->tetrahedra[next_t].vertices[next_vertex] != a1);

        const int cur_vertex = next_vertex;
        next_vertex = d->tetrahedra[next_t].index_in_neighbour[cur_vertex];
        next_t = d->tetrahedra[next_t].neighbours[cur_vertex];
      }
      /* Don't forget to add back last_t (which was overwritten) */
      int_lifo_queue_push(&d->tetrahedra_containing_vertex, last_t);
    }
  }
  return d->tetrahedra_containing_vertex.index;
}

/**
 * @brief Replace the given tetrahedron with four new ones by inserting the
 * given new vertex.
 *
 * @image html newvoronoicell_one_to_four_flip.png
 *
 * The original tetrahedron is positively oriented, and hence its vertex_indices
 * are ordered as shown in the figure. We construct four new tetrahedra by
 * replacing one of the four original vertex_indices with the new vertex. If we
 * keep the ordering of the vertex_indices, then the new tetrahedra will also be
 * positively oriented. The new neighbour relations can be easily deduced from
 * the figure, and the new index relations follow automatically from the way we
 * construct the new tetrahedra.
 *
 * For clarity, the common faces of the new tetrahedra are marked in green in
 * the figure.
 *
 * The first new tetrahedron replaces the original tetrahedron, while the three
 * extra new tetrahedra are added to the list.
 *
 * @param d Delaunay tessellation
 * @param v New vertex to insert.
 * @param t Tetrahedron to replace.
 */
inline static void delaunay_one_to_four_flip(struct delaunay* d, int v, int t) {
  delaunay_log("Flipping tetrahedron %i to 4 new ones.", t);

  /* Extract necessary information */
  const int vertices[4] = {
      d->tetrahedra[t].vertices[0], d->tetrahedra[t].vertices[1],
      d->tetrahedra[t].vertices[2], d->tetrahedra[t].vertices[3]};
  const int ngbs[4] = {
      d->tetrahedra[t].neighbours[0], d->tetrahedra[t].neighbours[1],
      d->tetrahedra[t].neighbours[2], d->tetrahedra[t].neighbours[3]};
  const int idx_in_ngbs[4] = {d->tetrahedra[t].index_in_neighbour[0],
                              d->tetrahedra[t].index_in_neighbour[1],
                              d->tetrahedra[t].index_in_neighbour[2],
                              d->tetrahedra[t].index_in_neighbour[3]};

  /* Replace t and create 3 new tetrahedra */
  delaunay_init_tetrahedron(d, t, vertices[0], vertices[1], vertices[2], v);
  const int t1 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, t1, vertices[0], vertices[1], v, vertices[3]);
  const int t2 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, t2, vertices[0], v, vertices[2], vertices[3]);
  const int t3 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, t3, v, vertices[1], vertices[2], vertices[3]);

  /* update neighbour relations */
  tetrahedron_swap_neighbours(&d->tetrahedra[t], t3, t2, t1, ngbs[3], 3, 3, 3,
                              idx_in_ngbs[3]);
  tetrahedron_swap_neighbours(&d->tetrahedra[t1], t3, t2, ngbs[2], t, 2, 2,
                              idx_in_ngbs[2], 2);
  tetrahedron_swap_neighbours(&d->tetrahedra[t2], t3, ngbs[1], t1, t, 1,
                              idx_in_ngbs[1], 1, 1);
  tetrahedron_swap_neighbours(&d->tetrahedra[t3], ngbs[0], t2, t1, t,
                              idx_in_ngbs[0], 0, 0, 0);

  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[0]], idx_in_ngbs[0], t3, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[1]], idx_in_ngbs[1], t2, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[2]], idx_in_ngbs[2], t1, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[3]], idx_in_ngbs[3], t, 3);

  /* enqueue all new/updated tetrahedra for delaunay checks */
  int_lifo_queue_push(&d->tetrahedra_to_check, t);
  int_lifo_queue_push(&d->tetrahedra_to_check, t1);
  int_lifo_queue_push(&d->tetrahedra_to_check, t2);
  int_lifo_queue_push(&d->tetrahedra_to_check, t3);
}

/**
 * @brief Replace the given two tetrahedra with six new ones by inserting the
 * given new vertex.
 *
 * @image html newvoronoicell_two_to_six_flip.png
 *
 * The two positively oriented tetrahedra (0123) and (0134) are replaced with
 * six new ones by replacing the common triangle vertex_indices one at a time:
 * (0125), (0523), (5123), (0154), (0534), and (5134). The new neighbour
 * relations can be easily deduced from the figure, while the new neighbour
 * indices follow automatically from the way we set up the tetrahedra.
 *
 * @param d Delaunay tessellation
 * @param v The new vertex
 * @param t Tetrahedra to replace
 */
inline static void delaunay_two_to_six_flip(struct delaunay* d, int v,
                                            const int* t) {
  /* Find the indices of the vertex_indices of the common triangle in both
   * tetrahedra
   */
  int triangle_indices[2][3];
  int num_vertices = 0;
  struct tetrahedron* t0 = &d->tetrahedra[t[0]];
  struct tetrahedron* t1 = &d->tetrahedra[t[1]];
  for (int current_vertex_idx_in_t0 = 0; current_vertex_idx_in_t0 < 4;
       current_vertex_idx_in_t0++) {
    int test_idx = 0;
    while (test_idx < 4 &&
           t0->vertices[current_vertex_idx_in_t0] != t1->vertices[test_idx]) {
      test_idx++;
    }
    if (test_idx < 4) {
      triangle_indices[0][num_vertices] = current_vertex_idx_in_t0;
      triangle_indices[1][num_vertices] = test_idx;
      num_vertices++;
    }
  }
  delaunay_assert(num_vertices == 3);

  /* Get the vertex of the first tetrahedron not shared with the second
   * tetrahedron */
  int top_idx_in_t0 = 6 - triangle_indices[0][0] - triangle_indices[0][1] -
                      triangle_indices[0][2];
  /* Make sure we have a positive permutation of 0123 */
  if (!positive_permutation(triangle_indices[0][0], triangle_indices[0][1],
                            top_idx_in_t0, triangle_indices[0][2])) {
    int tmp = triangle_indices[0][0];
    triangle_indices[0][0] = triangle_indices[0][1];
    triangle_indices[0][1] = tmp;

    tmp = triangle_indices[1][0];
    triangle_indices[1][0] = triangle_indices[1][1];
    triangle_indices[1][1] = tmp;
  }

  /* Set variables in accordance to figure */
  const int v0_0 = triangle_indices[0][0];
  const int v1_0 = triangle_indices[0][1];
  const int v2_0 = top_idx_in_t0;
  const int v3_0 = triangle_indices[0][2];

  const int v0_1 = triangle_indices[1][0];
  const int v1_1 = triangle_indices[1][1];
  const int v3_1 = triangle_indices[1][2];
  const int v4_1 = d->tetrahedra[t[0]].index_in_neighbour[v2_0];

  // now set some variables to the names in the documentation figure
  const int vert[6] = {
      d->tetrahedra[t[0]].vertices[v0_0], d->tetrahedra[t[0]].vertices[v1_0],
      d->tetrahedra[t[0]].vertices[v2_0], d->tetrahedra[t[0]].vertices[v3_0],
      d->tetrahedra[t[1]].vertices[v4_1], v};

  const int ngbs[6] = {d->tetrahedra[t[0]].neighbours[v0_0],
                       d->tetrahedra[t[1]].neighbours[v0_1],
                       d->tetrahedra[t[1]].neighbours[v1_1],
                       d->tetrahedra[t[0]].neighbours[v1_0],
                       d->tetrahedra[t[0]].neighbours[v3_0],
                       d->tetrahedra[t[1]].neighbours[v3_1]};

  const int idx_in_ngbs[6] = {d->tetrahedra[t[0]].index_in_neighbour[v0_0],
                              d->tetrahedra[t[1]].index_in_neighbour[v0_1],
                              d->tetrahedra[t[1]].index_in_neighbour[v1_1],
                              d->tetrahedra[t[0]].index_in_neighbour[v1_0],
                              d->tetrahedra[t[0]].index_in_neighbour[v3_0],
                              d->tetrahedra[t[1]].index_in_neighbour[v3_1]};

  /* Overwrite the two existing tetrahedra and create 4 new ones */
  delaunay_init_tetrahedron(d, t[0], vert[0], vert[1], vert[2], vert[5]);
  delaunay_init_tetrahedron(d, t[1], vert[0], vert[5], vert[2], vert[3]);
  int tn2 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, tn2, vert[5], vert[1], vert[2], vert[3]);
  int tn3 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, tn3, vert[0], vert[1], vert[5], vert[4]);
  int tn4 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, tn4, vert[0], vert[5], vert[3], vert[4]);
  int tn5 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, tn5, vert[5], vert[1], vert[3], vert[4]);

  /* Update neighbour relations */
  tetrahedron_swap_neighbours(&d->tetrahedra[t[0]], tn2, t[1], tn3, ngbs[4], 3,
                              3, 3, idx_in_ngbs[4]);
  tetrahedron_swap_neighbours(&d->tetrahedra[t[1]], tn2, ngbs[3], tn4, t[0], 1,
                              idx_in_ngbs[3], 3, 1);
  tetrahedron_swap_neighbours(&d->tetrahedra[tn2], ngbs[0], t[1], tn5, t[0],
                              idx_in_ngbs[0], 0, 3, 0);
  tetrahedron_swap_neighbours(&d->tetrahedra[tn3], tn5, tn4, ngbs[5], t[0], 2,
                              2, idx_in_ngbs[5], 2);
  tetrahedron_swap_neighbours(&d->tetrahedra[tn4], tn5, ngbs[2], tn3, t[1], 1,
                              idx_in_ngbs[2], 1, 2);
  tetrahedron_swap_neighbours(&d->tetrahedra[tn5], ngbs[1], tn4, tn3, tn2,
                              idx_in_ngbs[1], 0, 0, 2);

  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[0]], idx_in_ngbs[0], tn2, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[1]], idx_in_ngbs[1], tn5, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[2]], idx_in_ngbs[2], tn4, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[3]], idx_in_ngbs[3], t[1], 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[4]], idx_in_ngbs[4], t[0], 3);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[5]], idx_in_ngbs[5], tn3, 2);

  /* Add new/updated tetrahedra to queue for checking */
  int_lifo_queue_push(&d->tetrahedra_to_check, t[0]);
  int_lifo_queue_push(&d->tetrahedra_to_check, t[1]);
  int_lifo_queue_push(&d->tetrahedra_to_check, tn2);
  int_lifo_queue_push(&d->tetrahedra_to_check, tn3);
  int_lifo_queue_push(&d->tetrahedra_to_check, tn4);
  int_lifo_queue_push(&d->tetrahedra_to_check, tn5);
}

/**
 * @brief Replace the given @f$n@f$ tetrahedra with @f$2n@f$ new ones by
 * inserting the given new vertex.
 *
 * @image html newvoronoicell_n_to_2n_flip.png
 *
 * The @f$n@f$ tetrahedra
 * (v0 v@f$(n+1)@f$ v1 v@f$(n)@f$),
 * @f$...@f$,
 * (v@f$(i-1)@f$ v@f$(n+1)@f$ v@f$(i)@f$ v@f$(n)@f$),
 * @f$...@f$,
 * (v@f$(n-1)@f$ v@f$(n+1)@f$ v0 v@f$(n)@f$)
 * are replaced with the @f$2n@f$ tetrahedra
 * (v0 v@f$(n+2)@f$ v1 v@f$(n)@f$),
 * (v0 v@f$(n+1)@f$ v1 v@f$(n+2)@f$),
 * @f$...@f$,
 * (v@f$(i-1)@f$ v@f$(n+2)@f$ v@f$(i)@f$ v@f$(n)@f$),
 * (v@f$(i-1)@f$ v@f$(n+1)@f$ v@f$(i)@f$ v@f$(n+2)@f$),
 * @f$...@f$,
 * (v@f$(n-1)@f$ v@f$(n+2)@f$ v0 v@f$(n)@f$),
 * (v@f$(n-1)@f$ v@f$(n+1)@f$ v0 v@f$(n+2)@f$).
 *
 * The new neighbour relations can easily be deduced from the figure, while the
 * new neighbour indices are set automatically by the way the new tetrahedra are
 * constructed.
 *
 * @param d Delaunay tessellation
 * @param v The new vertex
 * @param t The tetrahedra to replace
 * @param n The number of tetrahedra to replace
 */
inline static void delaunay_n_to_2n_flip(struct delaunay* d, int v,
                                         const int* t, int n) {
  /* find the indices of the common axis vertex_indices in all tetrahedra */
  int axis_idx_in_tj[n][2];
  int tn_min_1_idx_in_t0 = 0;
  int num_axis = 0;
  struct tetrahedron* t0 = &d->tetrahedra[t[0]];
  for (int cur_v_idx_in_t0 = 0; cur_v_idx_in_t0 < 4; cur_v_idx_in_t0++) {
    int current_vertex_idx_in_tj[n];
    current_vertex_idx_in_tj[0] = cur_v_idx_in_t0;
    int current_vertex_is_axis = 1;
    for (int j = 1; j < n; ++j) {
      struct tetrahedron* tj = &d->tetrahedra[t[j]];
      int test_idx = 0;
      while (test_idx < 4 &&
             t0->vertices[cur_v_idx_in_t0] != tj->vertices[test_idx]) {
        test_idx++;
      }
      current_vertex_is_axis &= (test_idx < 4);
      current_vertex_idx_in_tj[j] = test_idx;
    }
    if (current_vertex_is_axis) {
      for (int j = 0; j < n; ++j) {
        axis_idx_in_tj[j][num_axis] = current_vertex_idx_in_tj[j];
      }
      ++num_axis;
    } else if (current_vertex_idx_in_tj[1] < 4) {
      /* Vertex is not an axis vertex, but is present in t1. This means that the
       * tetrahedron opposite of this vertex must be t_(n-1) (see figure). */
      tn_min_1_idx_in_t0 = current_vertex_idx_in_tj[0];
    }
  }
  /* Found both vertex_indices of the edge shared by all tetrahedra? */
  delaunay_assert(num_axis == 2);

  /* now make sure we give the indices the same meaning as in the figure */
  const int t1_idx_in_t0 =
      6 - axis_idx_in_tj[0][0] - axis_idx_in_tj[0][1] - tn_min_1_idx_in_t0;
  if (!positive_permutation(t1_idx_in_t0, axis_idx_in_tj[0][0],
                            tn_min_1_idx_in_t0, axis_idx_in_tj[0][1])) {
    for (int j = 0; j < n; ++j) {
      const int tmp = axis_idx_in_tj[j][0];
      axis_idx_in_tj[j][0] = axis_idx_in_tj[j][1];
      axis_idx_in_tj[j][1] = tmp;
    }
  }

#ifdef DELAUNAY_CHECKS
  /* Check that the new vertex lies in between the axis vertices */
  const double* v_new = &d->rescaled_vertices[3 * v];
  const double* a0 =
      &d->rescaled_vertices[3 * t0->vertices[axis_idx_in_tj[0][0]]];
  const double* a1 =
      &d->rescaled_vertices[3 * t0->vertices[axis_idx_in_tj[0][1]]];
  const double vec1[3] = {v_new[0] - a0[0], v_new[1] - a0[1], v_new[2] - a0[2]};
  const double vec2[3] = {a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]};
  const double norm1 = geometry3d_dot(vec1, vec1);
  const double norm2 = geometry3d_dot(vec2, vec2);
  const double dot12 = geometry3d_dot(vec1, vec2);
  if (fabs(dot12 / sqrt(norm1 * norm2) - 1) > 1e-6) {
    abort();
  }
  if (dot12 < 0 || dot12 > norm2) {
    abort();
  }

  /* Check that this t0 has the right orientation using this notation */
  int orientation = delaunay_test_orientation(
      d, t0->vertices[t1_idx_in_t0], t0->vertices[axis_idx_in_tj[0][0]],
      t0->vertices[tn_min_1_idx_in_t0], t0->vertices[axis_idx_in_tj[0][1]]);
  if (orientation >= 0) {
    fprintf(stderr, "Incorrect orientation with current notation!");
    abort();
  }
#endif

  /* set some variables to the values in the documentation figure */
  int vert[n + 3], ngbs[2 * n], idx_in_ngb[2 * n];
  int tprev_in_tcur = tn_min_1_idx_in_t0;
  for (int j = 0; j < n; ++j) {
    const int tnext_in_tcur =
        6 - tprev_in_tcur - axis_idx_in_tj[j][0] - axis_idx_in_tj[j][1];
    struct tetrahedron* tj = &d->tetrahedra[t[j]];
    vert[j] = tj->vertices[tnext_in_tcur];
    tprev_in_tcur = tj->index_in_neighbour[tnext_in_tcur];
    ngbs[2 * j] = tj->neighbours[axis_idx_in_tj[j][0]];
    ngbs[2 * j + 1] = tj->neighbours[axis_idx_in_tj[j][1]];
    idx_in_ngb[2 * j] = tj->index_in_neighbour[axis_idx_in_tj[j][0]];
    idx_in_ngb[2 * j + 1] = tj->index_in_neighbour[axis_idx_in_tj[j][1]];
  }
  int vn = d->tetrahedra[t[0]].vertices[axis_idx_in_tj[0][1]];
  int vn_plus_1 = d->tetrahedra[t[0]].vertices[axis_idx_in_tj[0][0]];
  int vn_plus_2 = v;

  /* create n new tetrahedra and overwrite the n existing ones */
  int tn[2 * n];
  for (int j = 0; j < n; j++) {
    tn[2 * j] = t[j];
    tn[2 * j + 1] = delaunay_new_tetrahedron(d);
  }
  for (int j = 0; j < n; j++) {
    /* Non-axis vertices */
    int vj = vert[j];
    int vj_plus_1 = vert[(j + 1) % n];

    /* Upper tetrahedron (connected to axis1 = v_n, see figure) */
    int tn0 = tn[2 * j];
    delaunay_init_tetrahedron(d, tn0, vj, vn_plus_2, vj_plus_1, vn);

    /* Lower tetrahedron (connected to axis0 = vn_plus_1, see figure) */
    int tn1 = tn[2 * j + 1];
    delaunay_init_tetrahedron(d, tn1, vj_plus_1, vn_plus_2, vj, vn_plus_1);

    /* Setup neighbour relations */
    /* Upper tetrahedron (see figure) */
    int t_next_upper = tn[(2 * (j + 1)) % (2 * n)];
    int t_prev_upper = tn[(2 * (j - 1) + 2 * n) % (2 * n)];
    int t_ngb_upper = ngbs[2 * j];
    int idx_in_t_ngb_upper = idx_in_ngb[2 * j];
    tetrahedron_swap_neighbours(&d->tetrahedra[tn0], t_next_upper, t_ngb_upper,
                                t_prev_upper, tn1, 2, idx_in_t_ngb_upper, 0, 3);
    tetrahedron_swap_neighbour(&d->tetrahedra[t_ngb_upper], idx_in_t_ngb_upper,
                               tn0, 1);

    /* Lower tetrahedron (see figure) */
    int t_next_lower = tn[(2 * (j + 1) + 1) % (2 * n)];
    int t_prev_lower = tn[(2 * (j - 1) + 1 + 2 * n) % (2 * n)];
    int t_ngb_lower = ngbs[2 * j + 1];
    int idx_in_t_ngb_lower = idx_in_ngb[2 * j + 1];
    tetrahedron_swap_neighbours(&d->tetrahedra[tn1], t_prev_lower, t_ngb_lower,
                                t_next_lower, tn0, 2, idx_in_t_ngb_lower, 0, 3);
    tetrahedron_swap_neighbour(&d->tetrahedra[t_ngb_lower], idx_in_t_ngb_lower,
                               tn1, 1);
  }

#ifdef DELAUNAY_CHECKS
  for (int j = 0; j < 2 * n; j++) {
    /* Check neighbour relations */
    int t_idx = tn[j];
    struct tetrahedron* this_t = &d->tetrahedra[t_idx];
    for (int i = 0; i < 4; i++) {
      struct tetrahedron* t_ngb = &d->tetrahedra[this_t->neighbours[i]];
      int idx_in_t_ngb = this_t->index_in_neighbour[i];
      if (t_ngb->neighbours[idx_in_t_ngb] != t_idx) {
        fprintf(stderr, "Wrong neighbour relation!");
        abort();
      }
    }
  }
#endif

  /* add new/updated tetrahedra to the queue for checking */
  for (int j = 0; j < 2 * n; j++) {
    int_lifo_queue_push(&d->tetrahedra_to_check, tn[j]);
  }
}

/**
 * @brief Replace the given two tetrahedra with three new tetrahedra.
 *
 * @image html newvoronoicell_two_to_three_flip.png
 *
 * The two positively oriented tetrahedra (v0 v1 v2 v3) and (v0 v1 v3 v4) (that
 * share the red face in the figure) are
 * replaced with three new positively oriented tetrahedra (that share the green
 * faces): (v0 v1 v2 v4), (v0 v4 v2 v3), and (v4 v1 v2 v3).
 *
 * Before the flip, t0 has ngb0, ngb3 and ngb4 as neighbours, while t1 has ngb1,
 * ngb2 and ngb5 as neighbours. After the flip, t'0 has ngb4 and ngb5 as
 * neighbours, t'1 has ngb3 and ngb4 as neighbours, and t'2 has ngb0 and ngb1
 * as neighbours.
 *
 * We first figure out the indices of the common triangle vertex_indices v0, v1
 * and v3 in both tetrahedra. Once we know these, it is very straigthforward to
 * match each index to one of these three vertex_indices (requiring that t0 is
 * positively oriented). We can then get the actual vertex_indices, neighbours
 * and neighbour indices and construct the new tetrahedra.
 *
 * @param d Delaunay tessellation
 * @param t0 First tetrahedron
 * @param t1 Second tetrahedron
 * @param top0 Index of the vertex of the first tetrahedron opposite the second
 * tetrahedron.
 * @param top1 Index of the vertex of the second tetrahedron opposite the first
 * tetrahedron.
 */
inline static void delaunay_two_to_three_flip(struct delaunay* restrict d,
                                              int t0, int t1, int top0,
                                              int top1) {
  /* get the indices of the common triangle of the tetrahedra, and make sure we
   * know which index in tetrahedron0 matches which index in tetrahedron1 */
  int triangle[2][3];
  for (int i = 0; i < 3; ++i) {
    triangle[0][i] = (top0 + i + 1) % 4;
    triangle[1][i] = 0;
    while (d->tetrahedra[t0].vertices[triangle[0][i]] !=
           d->tetrahedra[t1].vertices[triangle[1][i]]) {
      ++triangle[1][i];
    }
  }
  /* Make sure that we start from a positively oriented tetrahedron. The weird
   * index ordering is chosen to match the vertex_indices in the documentation
   * figure
   */
  if (!positive_permutation(triangle[0][1], triangle[0][2], top0,
                            triangle[0][0])) {
    int tmp = triangle[0][1];
    triangle[0][1] = triangle[0][2];
    triangle[0][2] = tmp;

    tmp = triangle[1][1];
    triangle[1][1] = triangle[1][2];
    triangle[1][2] = tmp;
  }
  const int v0_0 = triangle[0][1];
  const int v1_0 = triangle[0][2];
  const int v2_0 = top0;
  const int v3_0 = triangle[0][0];

  const int v0_1 = triangle[1][1];
  const int v1_1 = triangle[1][2];
  const int v3_1 = triangle[1][0];
  const int v4_1 = top1;

  /* Set some variables to the names used in the documentation figure */
  int vert[5], ngbs[6], idx_in_ngb[6];
  vert[0] = d->tetrahedra[t0].vertices[v0_0];
  vert[1] = d->tetrahedra[t0].vertices[v1_0];
  vert[2] = d->tetrahedra[t0].vertices[v2_0];
  vert[3] = d->tetrahedra[t0].vertices[v3_0];
  ngbs[0] = d->tetrahedra[t0].neighbours[v0_0];
  ngbs[3] = d->tetrahedra[t0].neighbours[v1_0];
  ngbs[4] = d->tetrahedra[t0].neighbours[v3_0];
  idx_in_ngb[0] = d->tetrahedra[t0].index_in_neighbour[v0_0];
  idx_in_ngb[3] = d->tetrahedra[t0].index_in_neighbour[v1_0];
  idx_in_ngb[4] = d->tetrahedra[t0].index_in_neighbour[v3_0];
  vert[4] = d->tetrahedra[t1].vertices[v4_1];
  ngbs[1] = d->tetrahedra[t1].neighbours[v0_1];
  ngbs[2] = d->tetrahedra[t1].neighbours[v1_1];
  ngbs[5] = d->tetrahedra[t1].neighbours[v3_1];
  idx_in_ngb[1] = d->tetrahedra[t1].index_in_neighbour[v0_1];
  idx_in_ngb[2] = d->tetrahedra[t1].index_in_neighbour[v1_1];
  idx_in_ngb[5] = d->tetrahedra[t1].index_in_neighbour[v3_1];

  /* overwrite t0 and t1 and create a new tetrahedron */
  delaunay_init_tetrahedron(d, t0, vert[0], vert[1], vert[2], vert[4]);
  delaunay_init_tetrahedron(d, t1, vert[0], vert[4], vert[2], vert[3]);
  const int t2 = delaunay_new_tetrahedron(d);
  delaunay_init_tetrahedron(d, t2, vert[4], vert[1], vert[2], vert[3]);

  /* fix neighbour relations */
  tetrahedron_swap_neighbours(&d->tetrahedra[t0], t2, t1, ngbs[5], ngbs[4], 3,
                              3, idx_in_ngb[5], idx_in_ngb[4]);
  tetrahedron_swap_neighbours(&d->tetrahedra[t1], t2, ngbs[3], ngbs[2], t0, 1,
                              idx_in_ngb[3], idx_in_ngb[2], 1);
  tetrahedron_swap_neighbours(&d->tetrahedra[t2], ngbs[0], t1, ngbs[1], t0,
                              idx_in_ngb[0], 0, idx_in_ngb[1], 0);

  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[0]], idx_in_ngb[0], t2, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[1]], idx_in_ngb[1], t2, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[2]], idx_in_ngb[2], t1, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[3]], idx_in_ngb[3], t1, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[4]], idx_in_ngb[4], t0, 3);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[5]], idx_in_ngb[5], t0, 2);

  /* add new/updated tetrahedrons to queue */
  int_lifo_queue_push(&d->tetrahedra_to_check, t0);
  int_lifo_queue_push(&d->tetrahedra_to_check, t1);
  int_lifo_queue_push(&d->tetrahedra_to_check, t2);
}

/**
 * @brief Replace the given four tetrahedra with four new tetrahedra.
 *
 * @image html newvoronoicell_four_to_four_flip.png
 *
 * The four positively oriented tetrahedra (v0 v1 v2 v3), (v0 v1 v3 v4),
 * (v0 v1 v5 v2), and (v0 v5 v1 v4),
 * that share the edge (v0 v1) are replaced by four new positively oriented
 * tetrahedra that share the edge (v3 v5) (dashed line in figure):
 * (v0 v3 v5 v2), (v1 v5 v3 v2), (v0 v5 v3 v4), and (v1 v3 v5 v4).
 *
 * The original red shared faces are hence replaced by the new green shared
 * faces in the figure; the blue shared faces will also be shared by the new
 * tetrahedra.
 *
 * Originally, t0 has ngb0 and ngb3 as neighbours, t1 has ngb1 and ngb2 as
 * neighbours, t2 has ngb4 and ngb7 as neighbours, and t3 has ngb5 and ngb6 as
 * neighbours. After the flip, t'0 has ngb3 and ngb7 as neighbours, t'1 has ngb0
 * and ngb4 as neighbours, t'2 has ngb2 and ngb6 as neighbours, and t'3 has ngb1
 * and ngb5 as neighbours.
 *
 * The tetrahedra should be given to this routine in the expected order: t0
 * should have t1 and t2 as neighbours, while t3 should be a common neighbour of
 * t1 and t2, but not of t3.
 *
 * The first thing we do is figure out how the internal vertex indices of the
 * four tetrahedra map to the names in the figure. We do this by identifying the
 * indices of the common axis vertex_indices v0 and v1 in all tetrahedra, and
 * the index of v2 in t0. Once we know v0, v1 and v2 in t0, we can deduce the
 * index of v3 in in t0, and require that (v0 v1 v2 v3) is positively oriented
 * (if not, we swap v0 and v1 in all tetrahedra so that it is).
 *
 * Once the indices have been mapped, it is straightforward to deduce the
 * actual vertex_indices, neighbours and neighbour indices, and we can simply
 * construct the new tetrahedra based on the figure.
 *
 * @param d Delaunay tessellation
 * @param t0 First tetrahedron.
 * @param t1 Second tetrahedron.
 * @param t2 Third tetrahedron.
 * @param t3 Fourth tetrahedron.
 */
inline static void delaunay_four_to_four_flip(struct delaunay* restrict d,
                                              int t0, int t1, int t2, int t3) {
  /* the four tetrahedra share an axis, find the indices of the axis points
   * in the four tetrahedra */
  int axis[4][4];
  int num_axis = 0;
  for (int i = 0; i < 4; ++i) {
    int idx_in_t0, idx_in_t1, idx_in_t2, idx_in_t3;
    idx_in_t0 = i;
    idx_in_t1 = 0;
    while (idx_in_t1 < 4 && d->tetrahedra[t0].vertices[idx_in_t0] !=
                                d->tetrahedra[t1].vertices[idx_in_t1]) {
      ++idx_in_t1;
    }
    idx_in_t2 = 0;
    while (idx_in_t2 < 4 && d->tetrahedra[t0].vertices[idx_in_t0] !=
                                d->tetrahedra[t2].vertices[idx_in_t2]) {
      ++idx_in_t2;
    }
    idx_in_t3 = 0;
    while (idx_in_t3 < 4 && d->tetrahedra[t0].vertices[idx_in_t0] !=
                                d->tetrahedra[t3].vertices[idx_in_t3]) {
      ++idx_in_t3;
    }
    if (idx_in_t1 < 4 && idx_in_t2 < 4 && idx_in_t3 < 4) {
      axis[0][num_axis] = idx_in_t0;
      axis[1][num_axis] = idx_in_t1;
      axis[2][num_axis] = idx_in_t2;
      axis[3][num_axis] = idx_in_t3;
      ++num_axis;
    } else {
      if (idx_in_t1 < 4) {
        axis[0][3] = idx_in_t0;
      }
    }
  }
  axis[0][2] = 6 - axis[0][0] - axis[0][1] - axis[0][3];

  /* now make sure we give the indices the same meaning as in the figure, i.e.
   * 'v[0][0]' is the index of vertex 0 in the figure in the first tetrahedron
   * (v0v1v2v3), and so on */
  if (!positive_permutation(axis[0][0], axis[0][1], axis[0][2], axis[0][3])) {
    int tmp = axis[0][0];
    axis[0][0] = axis[0][1];
    axis[0][1] = tmp;

    tmp = axis[1][0];
    axis[1][0] = axis[1][1];
    axis[1][1] = tmp;

    tmp = axis[2][0];
    axis[2][0] = axis[2][1];
    axis[2][1] = tmp;

    tmp = axis[3][0];
    axis[3][0] = axis[3][1];
    axis[3][1] = tmp;
  }

  /* t0 = (v0v1v2v3) */
  const int v0_0 = axis[0][0];
  const int v1_0 = axis[0][1];
  const int v2_0 = axis[0][2];
  const int v3_0 = axis[0][3];

  /* t1 = (v0v1v3v4) */
  const int v0_1 = axis[1][0];
  const int v1_1 = axis[1][1];
  const int v4_1 = d->tetrahedra[t0].index_in_neighbour[v2_0];

  /* t2 = (v0v1v5v2) */
  const int v0_2 = axis[2][0];
  const int v1_2 = axis[2][1];
  const int v5_2 = d->tetrahedra[t0].index_in_neighbour[v3_0];

  /* t3 = (v0v5v1v4) */
  const int v0_3 = axis[3][0];
  const int v1_3 = axis[3][1];

  const int vert[6] = {
      d->tetrahedra[t0].vertices[v0_0], d->tetrahedra[t0].vertices[v1_0],
      d->tetrahedra[t0].vertices[v2_0], d->tetrahedra[t0].vertices[v3_0],
      d->tetrahedra[t1].vertices[v4_1], d->tetrahedra[t2].vertices[v5_2]};

  const int ngbs[8] = {
      d->tetrahedra[t0].neighbours[v0_0], d->tetrahedra[t1].neighbours[v0_1],
      d->tetrahedra[t1].neighbours[v1_1], d->tetrahedra[t0].neighbours[v1_0],
      d->tetrahedra[t2].neighbours[v0_2], d->tetrahedra[t3].neighbours[v0_3],
      d->tetrahedra[t3].neighbours[v1_3], d->tetrahedra[t2].neighbours[v1_2]};

  const int idx_in_ngb[8] = {d->tetrahedra[t0].index_in_neighbour[v0_0],
                             d->tetrahedra[t1].index_in_neighbour[v0_1],
                             d->tetrahedra[t1].index_in_neighbour[v1_1],
                             d->tetrahedra[t0].index_in_neighbour[v1_0],
                             d->tetrahedra[t2].index_in_neighbour[v0_2],
                             d->tetrahedra[t3].index_in_neighbour[v0_3],
                             d->tetrahedra[t3].index_in_neighbour[v1_3],
                             d->tetrahedra[t2].index_in_neighbour[v1_2]};

  /* Replace the tetrahedra */
  /* t0 becomes (v0v3v5v2) */
  delaunay_init_tetrahedron(d, t0, vert[0], vert[3], vert[5], vert[2]);
  /* t1 becomes (v1v5v3v2) */
  delaunay_init_tetrahedron(d, t1, vert[1], vert[5], vert[3], vert[2]);
  /* t2 becomes (v0v5v3v4) */
  delaunay_init_tetrahedron(d, t2, vert[0], vert[5], vert[3], vert[4]);
  /* t3 becomes (v1v3v5v4) */
  delaunay_init_tetrahedron(d, t3, vert[1], vert[3], vert[5], vert[4]);

  /* Setup neighbour information */
  tetrahedron_swap_neighbours(&d->tetrahedra[t0], t1, ngbs[7], ngbs[3], t2, 0,
                              idx_in_ngb[7], idx_in_ngb[3], 3);
  tetrahedron_swap_neighbours(&d->tetrahedra[t1], t0, ngbs[0], ngbs[4], t3, 0,
                              idx_in_ngb[0], idx_in_ngb[4], 3);
  tetrahedron_swap_neighbours(&d->tetrahedra[t2], t3, ngbs[2], ngbs[6], t0, 0,
                              idx_in_ngb[2], idx_in_ngb[6], 3);
  tetrahedron_swap_neighbours(&d->tetrahedra[t3], t2, ngbs[5], ngbs[1], t1, 0,
                              idx_in_ngb[5], idx_in_ngb[1], 3);

  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[0]], idx_in_ngb[0], t1, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[1]], idx_in_ngb[1], t3, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[2]], idx_in_ngb[2], t2, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[3]], idx_in_ngb[3], t0, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[4]], idx_in_ngb[4], t1, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[5]], idx_in_ngb[5], t3, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[6]], idx_in_ngb[6], t2, 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[7]], idx_in_ngb[7], t0, 1);

  /* append updated tetrahedra to queue for checking */
  int_lifo_queue_push(&d->tetrahedra_to_check, t0);
  int_lifo_queue_push(&d->tetrahedra_to_check, t1);
  int_lifo_queue_push(&d->tetrahedra_to_check, t2);
  int_lifo_queue_push(&d->tetrahedra_to_check, t3);
}

/**
 * @brief Replace the given three tetrahedra with two new tetrahedra.
 *
 * @image html newvoronoicell_three_to_two_flip.png
 *
 * The three positively oriented tetrahedra (v0 v1 v2 v4), (v0 v4 v2 v3), and
 * (v4 v1 v2 v3) (with the red common faces in the figure) are
 * replaced with two new positively oriented tetrahedra (with a green common
 * face in the figure): (v0 v1 v2 v3) and (v0 v1 v3 v4).
 *
 * Originally, t0 has ngb4 and ngb5 as neighbours, t1 has ngb2 and ngb3 as
 * neighbours, and t2 has ngb0 and ngb1 as neighbours. After the flip, t'0 has
 * ngb0, ngb3 and ngb4 as neighbours, while t'1 has ngb1, ngb2 and ngb5 as
 * neighbours.
 *
 * We first find the indices of the common axis (v2 v4) in all three tetrahedra,
 * plus the index of v0 (the third common vertex of t0 and t1) in t0. Once we
 * have these we also know the index of v1 in t0, and we can find out which of
 * the two axis indices corresponds to v2 and which to v4 by requiring that the
 * four indices are a positively oriented permutation of 0123. Once this is
 * done, it is very straightforward to obtain the other indices in the other
 * tetrahedra. We can then get the actual vertex_indices, neighbours and
 * neighbour indices, and construct the two new tetrahedra.
 *
 * Note that because this flip removes a tetrahedron, it will free up a spot in
 * the tetrahedra vector. Since removing the tetrahedron from that vector would
 * be very expensive (since it requires a reshuffle of all tetrahedra behind it
 * and requires us to update the neighbour relations for all of these
 * tetrahedra), we just leave it in and keep an extra stack of free spots in the
 * tetrahedra array, which can be filled by other flips.
 *
 * @param d Delaunay tessellation
 * @param t0 First tetrahedron.
 * @param t1 Second tetrahedron.
 * @param t2 Third tetrahedron.
 * @param top_t0 The index of the top vertex in t0 (v1).
 * @param top_t1 The index of the top vertex in t1 (v3).
 * @param non_axis_1_in_t0 The index of the second non axis vertex in t0 (v0).
 * @return The index of the freed tetrahedron.
 */
inline static int delaunay_three_to_two_flip(struct delaunay* restrict d,
                                             int t0, int t1, int t2, int top_t0,
                                             int top_t1, int non_axis_in_t0) {
  /* Get the common axis of the three tetrahedra */
  /* First in t0 */
  int axis[3][2];
  axis[0][0] = (top_t0 + 1) % 4;
  if (axis[0][0] == non_axis_in_t0) axis[0][0] = (axis[0][0] + 1) % 4;
  axis[0][1] = 6 - axis[0][0] - top_t0 - non_axis_in_t0;
  if (!positive_permutation(non_axis_in_t0, top_t0, axis[0][0], axis[0][1])) {
    int tmp = axis[0][0];
    axis[0][0] = axis[0][1];
    axis[0][1] = tmp;
  }
  /* Then also in t1 and t2 */
  const int v_ax0 = d->tetrahedra[t0].vertices[axis[0][0]];
  const int v_ax1 = d->tetrahedra[t0].vertices[axis[0][1]];
  for (int i = 0; i < 4; ++i) {
    if (v_ax0 == d->tetrahedra[t1].vertices[i])
      axis[1][0] = i;
    else if (v_ax1 == d->tetrahedra[t1].vertices[i])
      axis[1][1] = i;
    if (v_ax0 == d->tetrahedra[t2].vertices[i])
      axis[2][0] = i;
    else if (v_ax1 == d->tetrahedra[t2].vertices[i])
      axis[2][1] = i;
  }

  /* Switch to a more convenient naming scheme */
  const int v0_0 = non_axis_in_t0;
  const int v1_0 = top_t0;
  const int v2_0 = axis[0][0];
  const int v4_0 = axis[0][1];

  const int v2_1 = axis[1][0];
  const int v3_1 = top_t1;
  const int v4_1 = axis[1][1];

  const int v2_2 = axis[2][0];
  const int v4_2 = axis[2][1];

  /* set some variables to the names used in the documentation figure */
  const int vert[5] = {
      d->tetrahedra[t0].vertices[v0_0], d->tetrahedra[t0].vertices[v1_0],
      d->tetrahedra[t0].vertices[v2_0], d->tetrahedra[t1].vertices[v3_1],
      d->tetrahedra[t0].vertices[v4_0],
  };

  const int ngbs[6] = {
      d->tetrahedra[t2].neighbours[v4_2], d->tetrahedra[t2].neighbours[v2_2],
      d->tetrahedra[t1].neighbours[v2_1], d->tetrahedra[t1].neighbours[v4_1],
      d->tetrahedra[t0].neighbours[v4_0], d->tetrahedra[t0].neighbours[v2_0]};

  const int idx_in_ngb[6] = {d->tetrahedra[t2].index_in_neighbour[v4_2],
                             d->tetrahedra[t2].index_in_neighbour[v2_2],
                             d->tetrahedra[t1].index_in_neighbour[v2_1],
                             d->tetrahedra[t1].index_in_neighbour[v4_1],
                             d->tetrahedra[t0].index_in_neighbour[v4_0],
                             d->tetrahedra[t0].index_in_neighbour[v2_0]};

  /* Overwrite two new tetrahedra and free the third one. */
  delaunay_init_tetrahedron(d, t0, vert[0], vert[1], vert[2], vert[3]);
  delaunay_init_tetrahedron(d, t1, vert[0], vert[1], vert[3], vert[4]);
  delaunay_log("Deactivating tetrahedron %i", t2);
  tetrahedron_deactivate(&d->tetrahedra[t2]);

  /* update neighbour relations */
  tetrahedron_swap_neighbours(&d->tetrahedra[t0], ngbs[0], ngbs[3], t1, ngbs[4],
                              idx_in_ngb[0], idx_in_ngb[3], 3, idx_in_ngb[4]);
  tetrahedron_swap_neighbours(&d->tetrahedra[t1], ngbs[1], ngbs[2], ngbs[5], t0,
                              idx_in_ngb[1], idx_in_ngb[2], idx_in_ngb[5], 2);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[0]], idx_in_ngb[0], t0, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[1]], idx_in_ngb[1], t1, 0);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[2]], idx_in_ngb[2], t1, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[3]], idx_in_ngb[3], t0, 1);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[4]], idx_in_ngb[4], t0, 3);
  tetrahedron_swap_neighbour(&d->tetrahedra[ngbs[5]], idx_in_ngb[5], t1, 2);

  /* add updated tetrahedra to queue */
  int_lifo_queue_push(&d->tetrahedra_to_check, t0);
  int_lifo_queue_push(&d->tetrahedra_to_check, t1);

  /* return invalidated tetrahedron */
  return t2;
}

/**
 * @brief Set the centroid of the tetrahedron
 *
 * @param d Delaunay tessellation
 * @param t The tetrahedron to finalize
 */
inline static void delaunay_finalize_tetrahedron(struct delaunay* d, int t) {
  struct tetrahedron* tet = &d->tetrahedra[t];
  double* v0 = &d->rescaled_vertices[3 * tet->vertices[0]];
  double* v1 = &d->rescaled_vertices[3 * tet->vertices[1]];
  double* v2 = &d->rescaled_vertices[3 * tet->vertices[2]];
  double* v3 = &d->rescaled_vertices[3 * tet->vertices[3]];

  for (int i = 0; i < 3; i++) {
    tet->centroid[i] = 0.25 * (v0[i] + v1[i] + v2[i] + v3[i]);
  }
}

/**
 * @brief Check the Delaunay criterion for tetrahedra in the queue until the
 * queue is empty.
 * @param d Delaunay tessellation
 * @param v The new vertex that might cause invalidation of tetrahedra.
 */
inline static void delaunay_check_tetrahedra(struct delaunay* d, int v) {
  int n_freed = 0;
  int* freed = (int*)malloc(10 * sizeof(int));
  int freed_size = 10;
  int freed_tetrahedron;
  int t = get_next_tetrahedron_to_check(d);
  while (t >= 0) {
    freed_tetrahedron = delaunay_check_tetrahedron(d, t, v);
    /* Did we free a tetrahedron? */
    if (freed_tetrahedron >= 0) {
      if (n_freed >= freed_size) {
        /* Grow array */
        freed_size <<= 1;
        freed = (int*)realloc(freed, freed_size * sizeof(int));
      }
      freed[n_freed] = freed_tetrahedron;
      n_freed++;
    }
    /* Pop next tetrahedron to check */
    t = get_next_tetrahedron_to_check(d);
  }
  /* Enqueue the newly freed tetrahedra indices */
  for (int i = 0; i < n_freed; i++) {
    int_lifo_queue_push(&d->free_tetrahedron_indices, freed[i]);
  }
  free(freed);
}

/**
 * @brief Check if the given tetrahedron satisfies the empty circumsphere
 * criterion that marks it as a Delaunay tetrahedron.
 *
 * If this check fails, this function also performs the necessary flips. All
 * new tetrahedra created by this function are also pushed to the queue for
 * checking.
 *
 * @param d Delaunay tessellation
 * @param t The tetrahedron to check.
 * @param v The new vertex that might cause invalidation of the tetrahedron.
 * @return Index of freed tetrahedron, or negative if no tetrahedra are freed
 */
inline static int delaunay_check_tetrahedron(struct delaunay* d, const int t,
                                             const int v) {
  struct tetrahedron* tetrahedron = &d->tetrahedra[t];
  const int v0 = tetrahedron->vertices[0];
  const int v1 = tetrahedron->vertices[1];
  const int v2 = tetrahedron->vertices[2];
  const int v3 = tetrahedron->vertices[3];

  /* Determine which vertex is the newly added vertex */
  int top;
  if (v == v0) {
    top = 0;
  } else if (v == v1) {
    top = 1;
  } else if (v == v2) {
    top = 2;
  } else if (v == v3) {
    top = 3;
  } else {
    fprintf(stderr,
            "Checking tetrahedron %i which does not contain the last added "
            "vertex %i",
            t, v);
    abort();
  }

  /* Get neighbouring tetrahedron opposite of newly added vertex */
  const int ngb = tetrahedron->neighbours[top];
  const int idx_in_ngb = tetrahedron->index_in_neighbour[top];
  /* Get the vertex in the neighbouring tetrahedron opposite of t */
  const int v4 = d->tetrahedra[ngb].vertices[idx_in_ngb];

  /* check if we have a neighbour that can be checked (dummies are not real and
     should not be tested) */
  if (ngb < 4) {
    delaunay_log("Dummy neighbour! Skipping checks for %i...", t);
    delaunay_assert(v4 == -1);
    return -1;
  }
  delaunay_assert(v4 != -1);

  const unsigned long* al = &d->integer_vertices[3 * v0];
  const unsigned long* bl = &d->integer_vertices[3 * v1];
  const unsigned long* cl = &d->integer_vertices[3 * v2];
  const unsigned long* dl = &d->integer_vertices[3 * v3];
  const unsigned long* el = &d->integer_vertices[3 * v4];

  const double* ad = &d->rescaled_vertices[3 * v0];
  const double* bd = &d->rescaled_vertices[3 * v1];
  const double* cd = &d->rescaled_vertices[3 * v2];
  const double* dd = &d->rescaled_vertices[3 * v3];
  const double* ed = &d->rescaled_vertices[3 * v4];

  const int test = geometry3d_in_sphere_adaptive(&d->geometry, al, bl, cl, dl,
                                                 el, ad, bd, cd, dd, ed);

  if (test < 0) {
    delaunay_log("Tetrahedron %i was invalidated by adding vertex %i", t, v);
    /* Figure out which flip is needed to restore the tetrahedra */
#ifdef DELAUNAY_3D_HAND_VEC
    int tests[4];
    geometry3d_orient_4(&d->geometry, al, bl, cl, dl, el, ad, bd, cd, dd, ed,
                        tests);
    tests[top] = -1;
#else
    int tests[4] = {-1, -1, -1, -1};
    if (top != 0)
      tests[0] = geometry3d_orient_adaptive(&d->geometry, bl, dl, cl, el, bd,
                                            dd, cd, ed);
    if (top != 1)
      tests[1] = geometry3d_orient_adaptive(&d->geometry, al, cl, dl, el, ad,
                                            cd, dd, ed);
    if (top != 2)
      tests[2] = geometry3d_orient_adaptive(&d->geometry, al, dl, bl, el, ad,
                                            dd, bd, ed);
    if (top != 3)
      tests[3] = geometry3d_orient_adaptive(&d->geometry, al, bl, cl, el, ad,
                                            bd, cd, ed);
#endif
    int i;
    for (i = 0; i < 4 && tests[i] < 0; ++i) {
    }
    if (i == 4) {
      /* v4 inside sphere around v0, v1, v2 and v3: need to do a 2 to 3 flip */
      delaunay_log("Performing 2 to 3 flip with %i and %i", t, ngb);
      delaunay_two_to_three_flip(d, t, ngb, top, idx_in_ngb);
    } else if (tests[i] == 0) {
      /* degenerate case: possible 4 to 4 flip needed. The line that connects v
       * and v4 intersects an edge of the triangle formed by the other 3
       * vertex_indices of t. If that edge is shared by exactly 4 tetrahedra in
       * total, the 2 neighbours are involved in the 4 to 4 flip. If it isn't,
       * we cannot solve this situation now, it will be solved later by another
       * flip. */

      /* get the other involved neighbour of t */
      /* the non_axis point is simply the vertex not present in the relevant
       * orientation test */
      const int other_ngb = d->tetrahedra[t].neighbours[i];
      /* get the index of 'new_vertex' in 'other_ngb', as the neighbour
       * opposite that vertex is the other neighbour we need to check */
      int idx_v_in_other_ngb;
      for (idx_v_in_other_ngb = 0;
           idx_v_in_other_ngb < 4 &&
           d->tetrahedra[other_ngb].vertices[idx_v_in_other_ngb] != v;
           idx_v_in_other_ngb++) {
      }
      const int other_ngbs_ngb =
          d->tetrahedra[other_ngb].neighbours[idx_v_in_other_ngb];
      /* check if other_ngbs_ngb is also a neighbour of ngb. */
      int second_idx_in_ngb =
          tetrahedron_is_neighbour(&d->tetrahedra[ngb], other_ngbs_ngb);
      if (second_idx_in_ngb < 4) {
        delaunay_log("Performing 4 to 4 flip between %i, %i, %i and %i!", t,
                     other_ngb, ngb, other_ngbs_ngb);
        delaunay_four_to_four_flip(d, t, other_ngb, ngb, other_ngbs_ngb);
      } else {
        delaunay_log("4 to 4 with %i and %i flip not possible!", t, ngb);
      }
    } else {
      /* check that this is indeed the only case left */
      delaunay_assert(tests[i] > 0);
      /* Outside: possible 3 to 2 flip.
       * The line that connects 'new_vertex' and 'v4' lies outside an edge of
       * the triangle formed by the other 3 vertex_indices of 'tetrahedron'. We
       * need to check if the neighbouring tetrahedron opposite the non-edge
       * point of that triangle is the same for 't' and 'ngb'. If it is, that is
       * the third tetrahedron for the 3 to 2 flip. If it is not, we cannot
       * solve this faulty situation now, but it will be solved by another flip
       * later on */

      /* get the other involved neighbour of t */
      /* the non_axis point is simply the vertex not present in the relevant
       * orientation test */
      const int other_ngb = d->tetrahedra[t].neighbours[i];
      /* check if other_ngb is also a neighbour of ngb */
      const int other_ngb_idx_in_ngb =
          tetrahedron_is_neighbour(&d->tetrahedra[ngb], other_ngb);
      if (other_ngb_idx_in_ngb < 4) {
        delaunay_log("Performing 3 to 2 flip with %i, %i and %i!", t, ngb,
                     other_ngb);

        return delaunay_three_to_two_flip(
            d, t, ngb, other_ngb, top, d->tetrahedra[t].index_in_neighbour[top],
            i);
      } else {
        delaunay_log("3 to 2 with %i and %i flip not possible!", t, ngb);
      }
    }
  } else {
    delaunay_log("Tetrahedron %i is valid!", t)
        delaunay_finalize_tetrahedron(d, t);
  }
  return -1;
}

/**
 * @brief Get the radius of the circumsphere of the given tetrahedron.
 *
 * @param d Delaunay tessellation.
 * @param t Tetrahedron index.
 * @return Radius of the circumsphere of the given tetrahedron.
 */
inline static double delaunay_get_radius(struct delaunay* restrict d, int t) {
  int v0 = d->tetrahedra[t].vertices[0];
  int v1 = d->tetrahedra[t].vertices[1];
  int v2 = d->tetrahedra[t].vertices[2];
  int v3 = d->tetrahedra[t].vertices[3];

  double* v0d = &d->rescaled_vertices[3 * v0];
  double* v1d = &d->rescaled_vertices[3 * v1];
  double* v2d = &d->rescaled_vertices[3 * v2];
  double* v3d = &d->rescaled_vertices[3 * v3];
  unsigned long* v0ul = &d->integer_vertices[3 * v0];
  unsigned long* v1ul = &d->integer_vertices[3 * v1];
  unsigned long* v2ul = &d->integer_vertices[3 * v2];
  unsigned long* v3ul = &d->integer_vertices[3 * v3];

  return geometry3d_compute_circumradius_adaptive(
      &d->geometry, v0d, v1d, v2d, v3d, v0ul, v1ul, v2ul, v3ul, d->side);
}

inline static void delaunay_compute_circumcentres(struct delaunay* d) {
  /* loop over the non-dummy tetrahedra in the Delaunay tessellation and compute
   * the midpoints of their circumspheres. These happen to be the vertices of
   * the Voronoi grid (because they are the points of equal distance to 3
   * generators, while the Voronoi edges are the lines of equal distance to 2
   * generators) */
  for (int i = 4; i < d->tetrahedra_index; i++) {
    struct tetrahedron* t = &d->tetrahedra[i];
    /* if the tetrahedron is inactive or not linked to a non-ghost, non-dummy
     * vertex, corresponding to an active particle, it is not a grid vertex and
     * we can skip it. */
    if (!t->active) {
      t->circumcenter[0] = NAN;
      t->circumcenter[1] = NAN;
      t->circumcenter[2] = NAN;
      continue;
    }
    /* Get the indices of the vertices of the tetrahedron */
    int v0 = t->vertices[0];
    int v1 = t->vertices[1];
    int v2 = t->vertices[2];
    int v3 = t->vertices[3];
    if ((v0 >= d->vertex_end || v0 < d->vertex_start) &&
        (v1 >= d->vertex_end || v1 < d->vertex_start) &&
        (v2 >= d->vertex_end || v2 < d->vertex_start) &&
        (v3 >= d->vertex_end || v3 < d->vertex_start)) {
      t->circumcenter[0] = NAN;
      t->circumcenter[1] = NAN;
      t->circumcenter[2] = NAN;
      continue;
    }
    /* Check that the vertices are valid */
    delaunay_assert(v0 >= 0 && v1 >= 0 && v2 >= 0 && v3 >= 0);

    /* Extract coordinates from the Delaunay vertices (generators) */
    if (v0 < d->vertex_start) {
      /* Dummy vertex!
       * This could mean that a neighbouring cell of this grids cell is empty,
       * or that we did not add all the necessary ghost vertex_indices to the
       * delaunay tesselation. */
      error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    double* v0d = &d->rescaled_vertices[3 * v0];
    unsigned long* v0ul = &d->integer_vertices[3 * v0];

    if (v1 < d->vertex_start) {
      error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    double* v1d = &d->rescaled_vertices[3 * v1];
    unsigned long* v1ul = &d->integer_vertices[3 * v1];

    if (v2 < d->vertex_start) {
      error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    double* v2d = &d->rescaled_vertices[3 * v2];
    unsigned long* v2ul = &d->integer_vertices[3 * v2];

    if (v3 < d->vertex_start) {
      error(
          "Vertex is part of tetrahedron with Dummy vertex! This could mean "
          "that one of the neighbouring cells is empty.");
    }
    double* v3d = &d->rescaled_vertices[3 * v3];
    unsigned long* v3ul = &d->integer_vertices[3 * v3];

    geometry3d_compute_circumcenter_adaptive(
        &d->geometry, v0d, v1d, v2d, v3d, v0ul, v1ul, v2ul, v3ul,
        t->circumcenter, d->side, d->anchor);

#ifdef SWIFT_DEBUG_CHECKS
    const double cx = t->circumcenter[0];
    const double cy = t->circumcenter[1];
    const double cz = t->circumcenter[2];

    double v0r[3], v1r[3], v2r[3], v3r[3];
    delaunay_get_vertex_at(d, v0, v0r);
    delaunay_get_vertex_at(d, v1, v1r);
    delaunay_get_vertex_at(d, v2, v2r);
    delaunay_get_vertex_at(d, v3, v3r);

    const double r0 =
        sqrt((cx - v0r[0]) * (cx - v0r[0]) + (cy - v0r[1]) * (cy - v0r[1]) +
             (cz - v0r[2]) * (cz - v0r[2]));
    const double r1 =
        sqrt((cx - v1r[0]) * (cx - v1r[0]) + (cy - v1r[1]) * (cy - v1r[1]) +
             (cz - v1r[2]) * (cz - v1r[2]));
    const double r2 =
        sqrt((cx - v2r[0]) * (cx - v2r[0]) + (cy - v2r[1]) * (cy - v2r[1]) +
             (cz - v2r[2]) * (cz - v2r[2]));
    const double r3 =
        sqrt((cx - v3r[0]) * (cx - v3r[0]) + (cy - v3r[1]) * (cy - v3r[1]) +
             (cz - v3r[2]) * (cz - v3r[2]));
    delaunay_assert(double_cmp(r0, r1, 1e5) && double_cmp(r0, r2, 1e5) &&
                    double_cmp(r0, r3, 1e5));
#endif
  }
}

/**
 * @brief Computes the delaunay search radii for a list of particles.
 *
 * For a given generator, we must loop over all the tetrahedra connected to this
 * generator and compute the maximal circumradius.
 * This is done by looping around all the delaunay edges connected to the
 * generator. We use a mask to avoid treating a tetrahedron twice.
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
  /* Prepare for calculation */
  /* Queue containing neighbouring generators of the currently active generators
   * and necessary info to proceed (a tetrahedron containing both and the index
   * of the neighbour in the tetrahedron). */
  struct int3_fifo_queue neighbour_info_queue;
  int3_fifo_queue_init(&neighbour_info_queue, 10);
  /* Flags indicating whether a neighbouring generator has already been added
   * to the queue for processing. */
  int* neighbour_flags = (int*)malloc(d->vertex_index * sizeof(int));
  bzero(neighbour_flags, d->vertex_index * sizeof(int));

  /* Flags indicating whether we already tested a tetrahedron */
  int* tetrahedra_flags = (int*)malloc(d->tetrahedra_index * sizeof(int));
  bzero(tetrahedra_flags, d->tetrahedra_index * sizeof(int));
  struct int_fifo_queue tetrahedra_to_check;
  int_fifo_queue_init(&tetrahedra_to_check, 25);

  /* loop over the particles */
  for (int i = 0; i < count; i++) {
    int gen_idx_in_d = parts[pid[i]].geometry.delaunay_vertex;

    /* start */
    neighbour_flags[gen_idx_in_d] = 1;

    /* Get a tetrahedron containing the central generator */
    int t_idx = d->vertex_tetrahedron_links[gen_idx_in_d];
    int gen_idx_in_t = d->vertex_tetrahedron_index[gen_idx_in_d];

    /* Pick another vertex (generator) from this tetrahedron and add it to the
     * queue */
    int other_v_idx_in_t = (gen_idx_in_t + 1) % 4;
    struct tetrahedron* t = &d->tetrahedra[t_idx];
    int other_v_idx_in_d = t->vertices[other_v_idx_in_t];
    /* Add the vertex info to the queue */
    int3 vertex_info = {
        ._0 = t_idx, ._1 = other_v_idx_in_d, ._2 = other_v_idx_in_t};
    int3_fifo_queue_push(&neighbour_info_queue, vertex_info);
    /* update flag of the other vertex indicating that it was added to the
     * queue*/
    neighbour_flags[other_v_idx_in_d] = 1;

    /* Loop over all neighbouring vertices of the current vertex */
    while (!int3_fifo_queue_is_empty(&neighbour_info_queue)) {
      /* Pop the next axis vertex and corresponding tetrahedron from the queue*/
      int3 info = int3_fifo_queue_pop(&neighbour_info_queue);
      int prev_t_idx = info._0;
      int axis_idx_in_d = info._1;
      int axis_idx_in_t = info._2;

      /* Add tetrahedron to queue */
      if (!tetrahedra_flags[prev_t_idx]) {
        int_fifo_queue_push(&tetrahedra_to_check, prev_t_idx);
        tetrahedra_flags[prev_t_idx] = 1;
      }

      struct tetrahedron* prev_t = &d->tetrahedra[prev_t_idx];

      /* Get a non axis vertex from first_t */
      int non_axis_idx_in_prev_t = (axis_idx_in_t + 1) % 4;
      if (prev_t->vertices[non_axis_idx_in_prev_t] == gen_idx_in_d) {
        non_axis_idx_in_prev_t = (non_axis_idx_in_prev_t + 1) % 4;
      }
      int non_axis_idx_in_d = prev_t->vertices[non_axis_idx_in_prev_t];

      if (!neighbour_flags[non_axis_idx_in_d]) {
        /* Add this vertex and tetrahedron to the neighbour queue and update
         * its flag */
        int3 new_info = {._0 = prev_t_idx,
                         ._1 = non_axis_idx_in_d,
                         ._2 = non_axis_idx_in_prev_t};
        int3_fifo_queue_push(&neighbour_info_queue, new_info);
        /* update flag */
        neighbour_flags[non_axis_idx_in_d] = 1;
      }

      /* Get a neighbouring tetrahedron of first_t sharing the axis */
      int cur_t_idx = prev_t->neighbours[non_axis_idx_in_prev_t];
      int prev_t_idx_in_cur_t =
          prev_t->index_in_neighbour[non_axis_idx_in_prev_t];

      /* Loop around the tetrahedra sharing this axis */
      int first_t_idx = prev_t_idx;
      while (cur_t_idx != first_t_idx) {
        /* Add tetrahedron to queue */
        if (!tetrahedra_flags[cur_t_idx]) {
          int_fifo_queue_push(&tetrahedra_to_check, cur_t_idx);
          tetrahedra_flags[cur_t_idx] = 1;
        }

#ifdef SWIFT_DEBUG_CHECKS
        if (axis_idx_in_d != d->tetrahedra[cur_t_idx].vertices[0] &&
            axis_idx_in_d != d->tetrahedra[cur_t_idx].vertices[1] &&
            axis_idx_in_d != d->tetrahedra[cur_t_idx].vertices[2] &&
            axis_idx_in_d != d->tetrahedra[cur_t_idx].vertices[3])
          error("Axis vertex no longer contained in current tet!");
        if (gen_idx_in_d != d->tetrahedra[cur_t_idx].vertices[0] &&
            gen_idx_in_d != d->tetrahedra[cur_t_idx].vertices[1] &&
            gen_idx_in_d != d->tetrahedra[cur_t_idx].vertices[2] &&
            gen_idx_in_d != d->tetrahedra[cur_t_idx].vertices[3])
          error("Central vertex no longer contained in current tet!");
#endif

        /* Update the variables */
        prev_t_idx = cur_t_idx;
        prev_t = &d->tetrahedra[prev_t_idx];
        /* get the next non axis vertex */
        non_axis_idx_in_prev_t = (prev_t_idx_in_cur_t + 1) % 4;
        non_axis_idx_in_d = prev_t->vertices[non_axis_idx_in_prev_t];
        while (non_axis_idx_in_d == axis_idx_in_d ||
               non_axis_idx_in_d == gen_idx_in_d) {
          non_axis_idx_in_prev_t = (non_axis_idx_in_prev_t + 1) % 4;
          non_axis_idx_in_d = prev_t->vertices[non_axis_idx_in_prev_t];
        }
        /* Add it to the queue if necessary */
        if (!neighbour_flags[non_axis_idx_in_d]) {
          /* Add this vertex and tetrahedron to the queue and update its flag */
          int3 new_info = {._0 = prev_t_idx,
                           ._1 = non_axis_idx_in_d,
                           ._2 = non_axis_idx_in_prev_t};
          int3_fifo_queue_push(&neighbour_info_queue, new_info);
          /* update flag */
          neighbour_flags[non_axis_idx_in_d] = 1;
        }
        /* Get the next tetrahedron sharing the same axis */
        cur_t_idx = prev_t->neighbours[non_axis_idx_in_prev_t];
        prev_t_idx_in_cur_t =
            prev_t->index_in_neighbour[non_axis_idx_in_prev_t];
      }
    }

    /* reset neighbour flags for next particle */
    neighbour_flags[gen_idx_in_d] = 0;
    for (int j = 0; j < neighbour_info_queue.end; j++) {
      delaunay_assert(neighbour_info_queue.values[j]._1 < d->vertex_index);
      neighbour_flags[neighbour_info_queue.values[j]._1] = 0;
    }
#ifdef DELAUNAY_CHECKS
    for (int k = 0; k < d->vertex_index; k++) {
      if (neighbour_flags[k]) {
        fprintf(stderr, "Found nonzero flag at end of get_radius!");
        abort();
      }
    }
#endif

    /* Reset neighbour info queue for next particle */
    int3_fifo_queue_reset(&neighbour_info_queue);

    /* Loop over the queued tetrahedra to calculate the search radius */
    double search_radius = 0.;
    while (!int_fifo_queue_is_empty(&tetrahedra_to_check)) {
      int t_to_check = int_fifo_queue_pop(&tetrahedra_to_check);

      /* Reset flag of current tetrahedron */
#ifdef SWIFT_DEBUG_CHECKS
      if (!tetrahedra_flags[t_to_check])
        error("Encountered tetrahedron which should not have to be checked!");
#endif
      tetrahedra_flags[t_to_check] = 0;

      /* Update search radius */
      search_radius =
          fmax(search_radius, 2. * delaunay_get_radius(d, t_to_check));
    }

    /* Reset tetrahedra queue */
    int_fifo_queue_reset(&tetrahedra_to_check);

    /* Save search radius */
    r[i] = search_radius;
  }

  /* Be clean */
  int3_fifo_queue_destroy(&neighbour_info_queue);
  free(neighbour_flags);
  int_fifo_queue_destroy(&tetrahedra_to_check);
  free(tetrahedra_flags);
}

/**
 * @brief Pop the next active tetrahedron to check from the end of the queue.
 *
 * If no more active tetrahedrons are queued, this function returns a negative
 * value.
 * Note that the returned tetrahedron index is effectively removed from the
 * queue and will be overwritten by subsequent calls to
 * delaunay_tetrahedron_enqueue().
 *
 * @param d Delaunay tessellation.
 * @return Index of the next tetrahedron to test, or -1 if the queue is empty.
 */
inline static int get_next_tetrahedron_to_check(struct delaunay* restrict d) {
  int active = 0;
  int t;
  while (!active && !int_lifo_queue_is_empty(&d->tetrahedra_to_check)) {
    t = int_lifo_queue_pop(&d->tetrahedra_to_check);
    active = d->tetrahedra[t].active;
  }
  return active ? t : -1;
}

/**
 * @brief Consolidate the Delaunay tessellation. This signals the end of the
 * addition of local, active vertices. All vertices added after this point are
 * considered to be ghost vertices.
 *
 * This function also performs some consistency checks if enabled.
 *
 * @param d Delaunay tessellation.
 */
inline static void delaunay_consolidate(struct delaunay* restrict d) {
  /* Set the ghost offset */
  d->vertex_end = d->vertex_index;

  /* perform a consistency test if enabled */
  delaunay_check_tessellation(d);
#ifdef DELAUNAY_CHECKS
  /* loop over all vertex_indices to check vertex-tetrahedron links */
  for (int v = 0; v < d->vertex_end; v++) {
    int t_idx = d->vertex_tetrahedron_links[v];
    struct tetrahedron* t = &d->tetrahedra[t_idx];
    int idx_in_t = d->vertex_tetrahedron_index[v];
    if (v != t->vertices[d->vertex_tetrahedron_index[v]]) {
      fprintf(stderr, "Wrong vertex-tetrahedron link!\n");
      fprintf(stderr, "\tVertex %i at index %i in\n", v, idx_in_t);
      fprintf(stderr, "\ttetrahedron %i: %i %i %i %i\n", t_idx, t->vertices[0],
              t->vertices[1], t->vertices[2], t->vertices[3]);
      abort();
    }
  }
#endif
}

inline static void delaunay_write_tessellation(
    const struct delaunay* restrict d, FILE* file, const size_t* offset) {

  for (int i = d->vertex_start; i < d->vertex_index; ++i) {
    double vertex[3] = {
        (d->rescaled_vertices[3 * i] - 1.) * d->side + d->anchor[0],
        (d->rescaled_vertices[3 * i + 1] - 1.) * d->side + d->anchor[1],
        (d->rescaled_vertices[3 * i + 2] - 1.) * d->side + d->anchor[2]};
    fprintf(file, "V\t%lu\t%g\t%g\t%g\n", *offset + i, vertex[0], vertex[1],
            vertex[2]);
  }
  for (int i = 4; i < d->tetrahedra_index; ++i) {
    if (!d->tetrahedra[i].active) {
      continue;
    }
    fprintf(file, "T\t%lu\t%lu\t%lu\t%lu\n",
            *offset + d->tetrahedra[i].vertices[0],
            *offset + d->tetrahedra[i].vertices[1],
            *offset + d->tetrahedra[i].vertices[2],
            *offset + d->tetrahedra[i].vertices[3]);
  }
}

inline static void delaunay_print_tessellation(
    const struct delaunay* restrict d, const char* file_name) {
  FILE* file = fopen(file_name, "w");

  delaunay_write_tessellation(d, file, 0);

  fclose(file);
}

inline static int delaunay_test_orientation(struct delaunay* restrict d, int v0,
                                            int v1, int v2, int v3) {
  return geometry3d_orient_adaptive(
      &d->geometry, &d->integer_vertices[3 * v0], &d->integer_vertices[3 * v1],
      &d->integer_vertices[3 * v2], &d->integer_vertices[3 * v3],
      &d->rescaled_vertices[3 * v0], &d->rescaled_vertices[3 * v1],
      &d->rescaled_vertices[3 * v2], &d->rescaled_vertices[3 * v3]);
}

inline static int delaunay_vertex_is_valid(struct delaunay* restrict d, int v) {
  return (v >= d->vertex_start && v < d->vertex_end) ||
         (v >= d->vertex_end && v < d->ghost_index + d->vertex_end);
}

/*! @brief Store the *actual* coordinates of the vertex at idx in out. */
inline static void delaunay_get_vertex_at(const struct delaunay* d, int idx,
                                          double* out) {
  out[0] = d->side * (d->rescaled_vertices[3 * idx] - 1.) + d->anchor[0];
  out[1] = d->side * (d->rescaled_vertices[3 * idx + 1] - 1.) + d->anchor[1];
  out[2] = d->side * (d->rescaled_vertices[3 * idx + 2] - 1.) + d->anchor[2];
}

/**
 * @brief Perform an (expensive) check on the tessellation to see that it is
 * still valid.
 *
 * This function will iterate over all tetrahedra (except the 3 dummy
 * tetrahedra) of the tessellation. For each tetrahedron, it checks that all
 * neighbour relations for that tetrahedron are set correctly (i.e. if
 * tetrahedron A is a neighbour of B, B should always also be a neighbour of A;
 * if tetrahedron A thinks it is the ith neighbour in tetrahedron B, it should
 * also be found in the ith neighbour position in tetrahedron B), and that none
 * fo the vertices of neighbouring tetrahedra are within the circumcircle of the
 * tetrahedron, which would violate the Delaunay criterion.
 *
 * Finally, this function also checks the vertex-tetrahedron information by
 * making sure that all tetrahedron indices stored in vertex_tetrahedra and
 * vertex_tetrahedron_index are correct.
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
  /* No expensive checks will be performed */
  return;
#endif

  /* loop over all non-dummy tetrahedra */
  for (int t0 = 4; t0 < d->tetrahedra_index; t0++) {
    /* Skip temporary deleted tetrahedra */
    if (!d->tetrahedra[t0].active) {
      continue;
    }
    int vt0_0 = d->tetrahedra[t0].vertices[0];
    int vt0_1 = d->tetrahedra[t0].vertices[1];
    int vt0_2 = d->tetrahedra[t0].vertices[2];
    int vt0_3 = d->tetrahedra[t0].vertices[3];

    int test_orientation =
        delaunay_test_orientation(d, vt0_0, vt0_1, vt0_2, vt0_3);
    if (test_orientation >= 0) {
      fprintf(stderr, "Tetrahedron %i has incorrect orientation!", t0);
      abort();
    }

    /* loop over neighbours */
    for (int i = 0; i < 4; i++) {
      int t_ngb = d->tetrahedra[t0].neighbours[i];
      /* check neighbour relations */
      int idx_in_ngb = d->tetrahedra[t0].index_in_neighbour[i];
      if (!d->tetrahedra[t_ngb].active) {
        fprintf(stderr, "Tetrahedron %i has an inactive neighbour: %i", t0,
                t_ngb);
      }
      if (d->tetrahedra[t_ngb].neighbours[idx_in_ngb] != t0) {
        fprintf(stderr, "Wrong neighbour!\n");
        fprintf(stderr, "Tetrahedron %i: %i %i %i %i\n", t0, vt0_0, vt0_1,
                vt0_2, vt0_3);
        fprintf(
            stderr, "\tNeighbours: %i %i %i %i\n",
            d->tetrahedra[t0].neighbours[0], d->tetrahedra[t0].neighbours[1],
            d->tetrahedra[t0].neighbours[2], d->tetrahedra[t0].neighbours[3]);
        fprintf(stderr, "\tIndex in neighbour: %i %i %i %i\n",
                d->tetrahedra[t0].index_in_neighbour[0],
                d->tetrahedra[t0].index_in_neighbour[1],
                d->tetrahedra[t0].index_in_neighbour[2],
                d->tetrahedra[t0].index_in_neighbour[3]);
        fprintf(
            stderr, "Neighbour tetrahedron %i: %i %i %i %i\n", t_ngb,
            d->tetrahedra[t_ngb].vertices[0], d->tetrahedra[t_ngb].vertices[1],
            d->tetrahedra[t_ngb].vertices[2], d->tetrahedra[t_ngb].vertices[3]);
        fprintf(stderr, "\tNeighbours: %i %i %i %i\n",
                d->tetrahedra[t_ngb].neighbours[0],
                d->tetrahedra[t_ngb].neighbours[1],
                d->tetrahedra[t_ngb].neighbours[2],
                d->tetrahedra[t_ngb].neighbours[3]);
        fprintf(stderr, "\tIndex in neighbour: %i %i %i %i\n",
                d->tetrahedra[t_ngb].index_in_neighbour[0],
                d->tetrahedra[t_ngb].index_in_neighbour[1],
                d->tetrahedra[t_ngb].index_in_neighbour[2],
                d->tetrahedra[t_ngb].index_in_neighbour[3]);
        abort();
      }
      if (t_ngb < 4) {
        /* Don't check delaunayness for dummy neighbour tetrahedra */
        continue;
      }
      /* check in-sphere criterion for delaunayness */
      int vertex_to_check = d->tetrahedra[t_ngb].vertices[idx_in_ngb];

      unsigned long* al = &d->integer_vertices[3 * vt0_0];
      unsigned long* bl = &d->integer_vertices[3 * vt0_1];
      unsigned long* cl = &d->integer_vertices[3 * vt0_2];
      unsigned long* dl = &d->integer_vertices[3 * vt0_3];
      unsigned long* el = &d->integer_vertices[3 * vertex_to_check];

      double* ad = &d->rescaled_vertices[3 * vt0_0];
      double* bd = &d->rescaled_vertices[3 * vt0_1];
      double* cd = &d->rescaled_vertices[3 * vt0_2];
      double* dd = &d->rescaled_vertices[3 * vt0_3];
      double* ed = &d->rescaled_vertices[3 * vertex_to_check];

      int test = geometry3d_in_sphere_adaptive(&d->geometry, al, bl, cl, dl, el,
                                               ad, bd, cd, dd, ed);
      if (test < 0) {
        fprintf(stderr, "Failed in-sphere test, value: %i!\n", test);
        fprintf(stderr, "\tTetrahedron %i: %i %i %i %i\n", t0, vt0_0, vt0_1,
                vt0_2, vt0_3);
        fprintf(stderr, "\tOpposite vertex: %i\n", vertex_to_check);
        abort();
      }
    }
  }
}

/**
 * @brief Check if abcd is a positive permutation of 0123 (meaning that if
 * 0123 are the vertex_indices of a positively ordered tetrahedron, then abcd
 * are also the vertex_indices of a positively ordered tetrahedron).
 *
 * @param a First index.
 * @param b Second index.
 * @param c Third index.
 * @param d Fourth index.
 * @return True if abcd is a positively oriented permutation of 0123.
 */
inline static int positive_permutation(int a, int b, int c, int d) {
  if ((a + 1) % 4 == b) {
    return c % 2 == 0;
  } else if ((a + 2) % 4 == b) {
    return b * c + a * d > b * d + a * c;
  } else {
    return d % 2 == 0;
  }
}

/** @brief Compute the normal to the triangle formed by `v0`, `v1` and `v2`.
 *
 * The will point towards an observer who sees the points in clockwise order.
 */
inline static void compute_normal(const double* restrict v0,
                                  const double* restrict v1,
                                  const double* restrict v2,
                                  double* restrict n) {
  double a[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
  double b[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};
  geometry3d_cross(b, a, n);
}

/**
 * @brief Calculate the cosine of the angle between the normal of the triangle
 * formed by `v0`, `v1` and `v2` and the direction from the centroid to `v`.
 */
inline static double compute_cos_theta(const double* restrict centroid,
                                       const double* restrict v0,
                                       const double* restrict v1,
                                       const double* restrict v2,
                                       const double* restrict v) {
  double n[3];
  compute_normal(v0, v1, v2, n);
  double dv[3] = {v[0] - centroid[0], v[1] - centroid[1], v[2] - centroid[2]};

  double cos_theta = geometry3d_dot(n, dv);
  double norm_2 = geometry3d_dot(n, n) * geometry3d_dot(dv, dv);
  double i_norm = norm_2 > 0. ? 1. / sqrt(norm_2) : 0.;
  return cos_theta * i_norm;
}

/** @brief Choose one of two candidate faces to cross while searching for a
 * tetrahedron containing `v`.
 *
 * The face which faces the new vertex "the most" is chosen. i.e. the face whose
 * normal makes the smallest angle with the direction of the new vertex `v`
 * (from the centroid of the tetrahedron).
 *
 * @param a, b The candidate neighbouring tetrahedra
 * @param centroid The centroid of the current tetrahedron
 * @param a0, a1, a2 The vertices of the face facing neighbour a
 * @param b0, b1, b2 The vertices of the face facing neighbour b
 * @return The index of the chosen neighbour.
 */
inline static int delaunay_choose_2(int a, int b,
                                    const double* restrict centroid,
                                    const double* a0, const double* a1,
                                    const double* a2, const double* b0,
                                    const double* b1, const double* b2,
                                    const double* restrict v) {
  /* Calculate the cosine of the angle between the normal on face a and the
   * vector from its centroid to the new vertex */
  double cos_theta_a = compute_cos_theta(centroid, a0, a1, a2, v);
  double cos_theta_b = compute_cos_theta(centroid, b0, b1, b2, v);

  if (cos_theta_a > cos_theta_b) {
    return a;
  }
  return b;
}

/** @brief Choose one of three candidate faces to cross while searching for a
 * tetrahedron containing `v`.
 *
 * The face which faces the new vertex "the most" is chosen.
 *
 * @return The index of the chosen face/neighbour.
 */
inline static int delaunay_choose_3(
    int a, int b, int c, const double* restrict centroid, const double* a0,
    const double* a1, const double* a2, const double* b0, const double* b1,
    const double* b2, const double* c0, const double* c1, const double* c2,
    const double* restrict v) {
  /* Calculate the cosine of the angle between the normal on face a and the
   * vector from its centroid to the new vertex */
  double cos_theta_a = compute_cos_theta(centroid, a0, a1, a2, v);
  double cos_theta_b = compute_cos_theta(centroid, b0, b1, b2, v);
  double cos_theta_c = compute_cos_theta(centroid, c0, c1, c2, v);
  if (cos_theta_a > cos_theta_b && cos_theta_a > cos_theta_c) {
    return a;
  } else if (cos_theta_b > cos_theta_a && cos_theta_b > cos_theta_c) {
    return b;
  }
  return c;
}

/** @brief Choose one of three candidate faces to cross while searching for a
 * tetrahedron containing `v`.
 *
 * The face will be chosen at random.
 *
 * @return The index of the chosen face/neighbour.
 */
inline static int delaunay_choose_random_3(int a, int b, int c) {
  /* Calculate the cosine of the angle between the normal on face a and the
   * vector from its centroid to the new vertex */
  int r = rand();
  if (r < RAND_MAX / 3) {
    return a;
  } else if (r < 2 * (RAND_MAX / 3)) {
    return b;
  }
  return c;
}

#endif  // SWIFTSIM_SHADOWSWIFT_DELAUNAY_3D_H
