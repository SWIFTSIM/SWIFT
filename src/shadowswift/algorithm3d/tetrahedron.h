//
// Created by yuyttenh on 10/05/22.
//

#ifndef SWIFTSIM_GITLAB_TETRAHEDRON_H
#define SWIFTSIM_GITLAB_TETRAHEDRON_H

typedef uint8_t tetrahedron_flags_t;
enum tetrahedron_flags {
  tetrahedron_flag_none = 0,
  tetrahedron_flag_has_vertex = (1 << 0),
  tetrahedron_flag_valid = (1 << 1),
  tetrahedron_flag_invalid = (1 << 2),
  tetrahedron_flag_validated =
      tetrahedron_flag_invalid | tetrahedron_flag_valid,
};

/**
 * @brief Tetrahedron
 *
 * A tetrahedron connects 4 points in 3D space, has 6 edges and 4 faces.
 *
 * In this struct, we store the indices of the 4 vertex_indices, the indices of
 * the neighbouring tetrahedra (that share a face with this tetrahedron) and the
 * indices of this tetrahedron in the neighbour lists of its neighbours.
 *
 * Conventions:
 * 1. neighbours[i] is the tetrahedron sharing the face oposite of
 * vertex_indices[i].
 * 2. index_in_neighbours follows the same ordering as neighbours
 */
struct tetrahedron {
  /*! @brief Indices in the associated delaunay tesselation of the particles
   * that make up the tetrahedron. */
  int vertices[4];

  /*! @brief Indices of the neighbour tetrahedra. */
  int neighbours[4];

  /*! @brief Index of this tetrahedron in the neighbour lists of its neighbours.
   * */
  int index_in_neighbour[4];

  union {
    /*! @brief The centroid of this tetrahedron. Used during construction. */
    double centroid[3];

    /*! @brief The circumcenter of this tetrahedron. Used for Voronoi
     * construction. */
    double circumcenter[3];
  };

  /*! @brief Indicates whether or not a tetrahedron is active (or has been
   * invalidated) */
  uint8_t active;

  tetrahedron_flags_t _flags;
};

/**
 * @brief Initialize a tetrahedron with the given vertex_indices.
 *
 * Neighbour information is set to nonsensical negative values
 *
 * @param t Tetrahedron
 * @param v0, v1, v2, v3 Vertices
 */
inline static void tetrahedron_init(struct tetrahedron *t, int v0, int v1,
                                    int v2, int v3) {
  t->vertices[0] = v0;
  t->vertices[1] = v1;
  t->vertices[2] = v2;
  t->vertices[3] = v3;

  t->neighbours[0] = -1;
  t->neighbours[1] = -1;
  t->neighbours[2] = -1;
  t->neighbours[3] = -1;

  t->index_in_neighbour[0] = -1;
  t->index_in_neighbour[1] = -1;
  t->index_in_neighbour[2] = -1;
  t->index_in_neighbour[3] = -1;

  t->active = 1;
  t->_flags = tetrahedron_flag_none;
}

/**
 * @brief Deactivates the given tetrahedron (set all values to nonsensical
 * negative values) and sets the active flag to 0
 * @param t Tetrahedron
 */
inline static void tetrahedron_deactivate(struct tetrahedron *restrict t) {
  t->vertices[0] = -1;
  t->vertices[1] = -1;
  t->vertices[2] = -1;
  t->vertices[3] = -1;

  t->neighbours[0] = -1;
  t->neighbours[1] = -1;
  t->neighbours[2] = -1;
  t->neighbours[3] = -1;

  t->index_in_neighbour[0] = -1;
  t->index_in_neighbour[1] = -1;
  t->index_in_neighbour[2] = -1;
  t->index_in_neighbour[3] = -1;

  t->active = 0;
  t->_flags = tetrahedron_flag_none;
}

/**
 * @brief Replace the neighbour at the given index with the given new neighbour.
 *
 * @param t Tetrahedron
 * @param index Index of the neighbour in the list of neighbours.
 * @param neighbour New neighbour.
 * @param index_in_neighbour Index of this tetrahedron in the neighbour list of
 * the new neighbour.
 */
inline static void tetrahedron_swap_neighbour(struct tetrahedron *restrict t,
                                              int index, int neighbour,
                                              int index_in_neighbour) {
  t->neighbours[index] = neighbour;
  t->index_in_neighbour[index] = index_in_neighbour;
}

/**
 * @brief Replace all neighbour relations at once.
 *
 * @param t Tetrahedron.
 * @param n0, n1, n2, n3 New neighbours.
 * @param idx_in_n0, idx_in_n1, idx_in_n2, idx_in_n3 Indices of this tetrahedron
 * the neighbour lists of its new neighbours.
 */
inline static void tetrahedron_swap_neighbours(struct tetrahedron *restrict t,
                                               int n0, int n1, int n2, int n3,
                                               int idx_in_n0, int idx_in_n1,
                                               int idx_in_n2, int idx_in_n3) {
  tetrahedron_swap_neighbour(t, 0, n0, idx_in_n0);
  tetrahedron_swap_neighbour(t, 1, n1, idx_in_n1);
  tetrahedron_swap_neighbour(t, 2, n2, idx_in_n2);
  tetrahedron_swap_neighbour(t, 3, n3, idx_in_n3);
}

/**
 * @brief Find the index of ngb in t
 * @param t Tetrahedron
 * @param ngb Index of potential neighbour in delaunay tesselation
 * @return Index of ngb in t or 4 if ngb is not a neighbour of t
 */
inline static int tetrahedron_is_neighbour(struct tetrahedron *restrict t,
                                           int ngb) {
  int i;
  for (i = 0; i < 4 && t->neighbours[i] != ngb; i++) {
  }
  return i;
}

#endif  // SWIFTSIM_GITLAB_TETRAHEDRON_H
