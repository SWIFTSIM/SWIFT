//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_VORONOI_3D_H
#define SWIFTSIM_SHADOWSWIFT_VORONOI_3D_H

#include "part.h"
#include "shadowswift/algorithm3d/delaunay.h"
#include "shadowswift/algorithm3d/geometry.h"
#include "shadowswift/algorithm3d/tetrahedron.h"
#include "shadowswift/queues.h"

/*! @brief The sid order in which the faces are stored in their array */
static const int face_sid_order[28] = {
    /* Local        */ 13,
    /* Boundary     */ 27,
    /* Cardinal directions (6 options) */
    /* ( 1,  0,  0) */ 22,
    /* (-1,  0,  0) */ 4,
    /* ( 0,  1,  0) */ 16,
    /* ( 0, -1,  0) */ 10,
    /* ( 0,  0,  1) */ 14,
    /* ( 0,  0, -1) */ 12,
    /* Change in 2 dimensions (12 options) */
    /* ( 1,  1,  0) */ 25,
    /* ( 1,  0,  1) */ 23,
    /* ( 1, -1,  0) */ 19,
    /* ( 1,  0, -1) */ 21,
    /* (-1,  1,  0) */ 7,
    /* (-1,  0,  1) */ 5,
    /* (-1, -1,  0) */ 1,
    /* (-1,  0, -1) */ 3,
    /* ( 0,  1,  1) */ 17,
    /* ( 0,  1, -1) */ 15,
    /* ( 0, -1,  1) */ 11,
    /* ( 0, -1, -1) */ 9,
    /* Change in all 3 directions (8 options) */
    /* ( 1,  1,  1) */ 26,
    /* ( 1,  1, -1) */ 24,
    /* ( 1, -1,  1) */ 20,
    /* (-1,  1,  1) */ 8,
    /* ( 1, -1, -1) */ 18,
    /* (-1,  1, -1) */ 6,
    /* (-1, -1,  1) */ 2,
    /* (-1, -1, -1) */ 0,
};

/*! @brief the position of the sid in the corresponding face arrays
 * (the inverse of the above). */
static const int face_sid_index[28] = {27, 14, 26, 15, 3,  13, 25, 12, 23, 19,
                                       5,  18, 7,  0,  6,  17, 4,  16, 24, 10,
                                       22, 11, 2,  9,  21, 8,  20, 1};

/**
 * @brief Voronoi interface.
 *
 * An interface is a connection between two neighbouring Voronoi cells. It is
 * completely defined by the indices of the generators that generate the two
 * neighbouring cells, a surface area and a midpoint position.
 */
struct voronoi_pair {
  /*! idx of the particle on the left of this pair in its respective swift
   * cell. Since the left particle is always local this is also the index of the
   * corresponding cell in this voronoi tesselation. */
  int left_idx;

  /*! idx of the particle on the right of this pair in its respective swift cell
   * if that cell is the same as the cell holding this Voronoi tesselation (i.e.
   * the particle is local) or in the super cell of its respective swift cell if
   * that swift cell is foreign. For local particles, this is also the index of
   * the corresponding cell in this voronoi tesselation. */
  int right_idx;

  /*! Real sid of this pair (boundary faces are stored under sid 27) */
  int sid;

  /*! Surface area of the interface. */
  double surface_area;

  /*! Midpoint of the interface. */
  double midpoint[3];

#ifdef VORONOI_STORE_FACES
  /*! Vertices of the interface. */
  double *vertices;

  /*! Number of vertices of this face. */
  int n_vertices;
#endif
};

struct voronoi {

  /*! @brief Voronoi cell pairs. We store these per (SWIFT) cell, i.e. pairs[0]
   *  contains all pairs that are completely contained within this cell, while
   *  pairs[1] corresponds to pairs crossing the boundary between this cell and
   *  the cell with coordinates that are lower in all coordinate directions (the
   *  cell to the left, front, bottom, sid=0), and so on. */
  struct voronoi_pair *pairs_flat;
  struct voronoi_pair *pairs[28];

  /*! @brief Current number of pairs per cell direction. */
  int pair_count[28];

  /*! @brief Allocated number of pairs */
  int pair_size;

  /*! @brief Total number of occupied pairs (all sids) */
  int pair_index;

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /*! @brief cell pair connections. Concatenation of the indices of the faces of
   * all parts. */
  int *cell_pair_connections;

  int cell_pair_connections_size;

  int cell_pair_connections_index;
#endif

  /*! @brief Flag indicating whether this voronoi struct is active (has memory
   * allocated) */
  int active;

  /*! @brief The absolute minimal surface area of faces in this voronoi
   * tessellation */
  double min_surface_area;

#ifdef VORONOI_STORE_FACES
  /*! @brief The face vertices of all the faces concatenated in one big array.
   * Every face has a pointer to the start of its vertices */
  double *face_vertices;
  int face_vertices_size;
  int face_vertices_index;
#endif
};

/* Forward declarations */
inline static int voronoi_new_face(struct voronoi *restrict v,
                                   const struct delaunay *restrict d,
                                   int left_part_idx_in_d,
                                   int right_part_idx_in_d,
                                   struct part *restrict parts, double area,
                                   double *restrict centroid,
                                   double *restrict vertices, int n_vertices);
inline static void voronoi_check_grid(struct voronoi *v,
                                      const struct delaunay *d,
                                      const struct part *parts);
inline static void voronoi_finalize(struct voronoi *v, const struct delaunay *d,
                                    struct part *parts);
inline static void voronoi_destroy(struct voronoi *restrict v);

inline static struct voronoi *voronoi_malloc(int number_of_cells, double dmin) {

  /* Allocate memory */
  struct voronoi *v = malloc(sizeof(struct voronoi));

  /* Allocate one chunk of memory for the voronoi pairs (faces). */
  v->pair_size = 9 * number_of_cells;
  v->pair_index = 0;
  v->pairs_flat = (struct voronoi_pair *)swift_malloc(
      "voronoi", v->pair_size * sizeof(struct voronoi_pair));
  for (int sid = 0; sid < 28; ++sid) {
    v->pair_count[sid] = 0;
    v->pairs[sid] = NULL;
  }

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /* Allocate memory for the cell_pair connections */
  v->cell_pair_connections_size = 17 * number_of_cells;
  v->cell_pair_connections_index = 0;
  v->cell_pair_connections =
      (int *)swift_malloc("voronoi", v->cell_pair_connections_size *
                                         sizeof(*v->cell_pair_connections));
#endif

#ifdef VORONOI_STORE_FACES
  /* Allocate memory for the face_vertices */
  v->face_vertices_size = 100 * number_of_cells;
  v->face_vertices = (double *)swift_malloc(
      "voronoi", 3 * v->face_vertices_size * sizeof(double));
  v->face_vertices_index = 0;
#endif

  v->min_surface_area = MIN_REL_FACE_SIZE * dmin;
  v->active = 1;

  return v;
}

inline static void voronoi_reset(struct voronoi *restrict v,
                                 int number_of_cells, double dmin) {
  voronoi_assert(v->active);

  /* reset indices for the voronoi pairs (faces). */
  v->pair_index = 0;
  for (int sid = 0; sid < 28; ++sid) {
    v->pair_count[sid] = 0;
    v->pairs[sid] = NULL;
  }

#ifdef VORONOI_STORE_FACES
  /* Reset face_vertices index */
  v->face_vertices_index = 0;
#endif

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /* Reset the cell_pair connections */
  v->cell_pair_connections_index = 0;
#endif

  v->min_surface_area = MIN_REL_FACE_SIZE * dmin;
}

/**
 * @brief Build the Voronoi grid based on the given Delaunay tessellation.
 *
 * This function allocates the memory for the Voronoi grid arrays and creates
 * the grid in linear time by
 *  1. Computing the grid vertices as the midpoints of the circumspheres of the
 *     Delaunay tetrahedra.
 *  2. Looping over all vertices and for each generator looping over all
 *     tetrahedra that link to that vertex. This is done by looping around all
 *     the Delaunay edges connected to that generator. While looping around an
 *     edge, for each tetrahedron, we add the edges of that tetrahedron which
 *     are connected to the current generator to a queue of next edges to loop
 *     around (if we did not already do so).
 *
 * During the second step, the geometrical properties (cell centroid, volume
 * and face midpoint, area) are computed as well.
 *
 * @param v Voronoi grid.
 * @param d Delaunay tessellation (read-only).
 * @param parts The particle array of the local cell.
 */
inline static void voronoi_build(struct voronoi *v, struct delaunay *d,
                                 struct part *parts) {

  /* the number of cells equals the number of non-ghost and non-dummy
     vertex_indices in the Delaunay tessellation */
  voronoi_assert(d->vertex_end > 0);

  /* Compute the circumcentres */
  delaunay_compute_circumcentres(d);

  /* Allocate memory for the neighbour flags and initialize them to 0 (will be
   * freed at the end of this function!) */
  int *neighbour_flags = (int *)malloc(d->vertex_index * sizeof(int));
  for (int i = 0; i < d->vertex_index; i++) {
    neighbour_flags[i] = 0;
  }

  /* Allocate a tetrahedron_vertex_queue (will be freed at the end of this
   * function!) */
  struct int3_fifo_queue neighbour_info_q;
  int3_fifo_queue_init(&neighbour_info_q, 10);

  /* The size of the array used to temporarily store the vertices of the voronoi
   * faces in */
  int face_vertices_size = 10;
  /* Temporary array to store face vertices in (will be freed at the end of this
   * function!) */
  double *face_vertices =
      (double *)malloc(3 * face_vertices_size * sizeof(double));
  int face_vertices_index = 0;

  /* loop over all cell generators, and hence over all non-ghost, non-dummy
     Delaunay vertex_indices */
  for (int gen_idx_in_d = d->vertex_start; gen_idx_in_d < d->vertex_end;
       gen_idx_in_d++) {

    /* Get the corresponding particle idx */
    int p_idx = d->vertex_part_idx[gen_idx_in_d];

    /* First reset the tetrahedron_vertex_queue */
    int3_fifo_queue_reset(&neighbour_info_q);

    /* Set the flag of the central generator so that we never pick it as
     * possible neighbour */
    neighbour_flags[gen_idx_in_d] = 1;

    /* Create a new voronoi cell for this generator */
    struct part *p = &parts[p_idx];
    double volume = 0.;
    double centroid[3] = {0., 0., 0.};
    int nface = 0;
#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
    int pair_connections_offset = v->cell_pair_connections_index;
#else
    int pair_connections_offset = 0;
#endif
    double min_ngb_dist2 = DBL_MAX;
    double min_ngb_dist_pos[3] = {0., 0., 0.};
    double generator_pos[3] = {p->x[0], p->x[1], p->x[2]};

    /* get the generator position, we use it during centroid/volume
       calculations */
    voronoi_assert(gen_idx_in_d < d->vertex_end);
    double ax = p->x[0];
    double ay = p->x[1];
    double az = p->x[2];

    /* Get a tetrahedron containing the central generator */
    int t_idx = d->vertex_tetrahedron_links[gen_idx_in_d];
    int gen_idx_in_t = d->vertex_tetrahedron_index[gen_idx_in_d];

    /* Pick another vertex (generator) from this tetrahedron and add it to the
     * queue */
    int other_v_idx_in_t = (gen_idx_in_t + 1) % 4;
    struct tetrahedron *t = &d->tetrahedra[t_idx];
    int other_v_idx_in_d = t->vertices[other_v_idx_in_t];
    int3 info = {._0 = t_idx, ._1 = other_v_idx_in_d, ._2 = other_v_idx_in_t};
    int3_fifo_queue_push(&neighbour_info_q, info);
    /* update flag of the other vertex */
    neighbour_flags[other_v_idx_in_d] = 1;

    int counter_outer = 0;
    while (!int3_fifo_queue_is_empty(&neighbour_info_q)) {
      if (counter_outer++ >= VORONOI_CONSTRUCTION_MAX_NGB_ITER)
        error(
            "Voronoi cell construction did not finish within allowed number of "
            "iterations! This probably means something went wrong with the "
            "Delaunay construction.");

      /* Initialize the area and centroid of this face */
      double face_area = 0.;
      double face_centroid[3] = {0., 0., 0.};

      /* Pop the next axis vertex and corresponding tetrahedron from the queue
       */
      info = int3_fifo_queue_pop(&neighbour_info_q);
      int first_t_idx = info._0;
      int axis_idx_in_d = info._1;
      int axis_idx_in_t = info._2;
      voronoi_assert(axis_idx_in_d >= d->vertex_start);
      struct tetrahedron *first_t = &d->tetrahedra[first_t_idx];

      /* Get a non axis vertex from first_t */
      int non_axis_idx_in_first_t = (axis_idx_in_t + 1) % 4;
      if (first_t->vertices[non_axis_idx_in_first_t] == gen_idx_in_d) {
        non_axis_idx_in_first_t = (non_axis_idx_in_first_t + 1) % 4;
      }
      int non_axis_idx_in_d = first_t->vertices[non_axis_idx_in_first_t];

      if (!neighbour_flags[non_axis_idx_in_d]) {
        /* Add this vertex and tetrahedron to the queue and update its flag */
        int3 new_info = {._0 = first_t_idx,
                         ._1 = non_axis_idx_in_d,
                         ._2 = non_axis_idx_in_first_t};
        int3_fifo_queue_push(&neighbour_info_q, new_info);
        neighbour_flags[non_axis_idx_in_d] |= 1;
      }

      /* Get a neighbouring tetrahedron of first_t sharing the axis */
      int cur_t_idx = first_t->neighbours[non_axis_idx_in_first_t];
      struct tetrahedron *cur_t = &d->tetrahedra[cur_t_idx];
      int prev_t_idx_in_cur_t =
          first_t->index_in_neighbour[non_axis_idx_in_first_t];

      /* Get a neighbouring tetrahedron of cur_t that is not first_t, sharing
       * the same axis */
      int next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
      while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
             cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
        next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
      }
      int next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];

      /* Get the next non axis vertex and add it to the queue if necessary */
      int next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
      if (!neighbour_flags[next_non_axis_idx_in_d]) {
        int3 new_info = {._0 = cur_t_idx,
                         ._1 = next_non_axis_idx_in_d,
                         ._2 = next_t_idx_in_cur_t};
        int3_fifo_queue_push(&neighbour_info_q, new_info);
        neighbour_flags[next_non_axis_idx_in_d] |= 1;
      }

      /* Get the coordinates of the voronoi vertex of the new face */
      const double *vor_vertex0 = d->tetrahedra[first_t_idx].circumcenter;
#ifdef VORONOI_STORE_FACES
      memcpy(&face_vertices[0], vor_vertex0, 3 * sizeof(*vor_vertex0));
      face_vertices_index = 1;
#endif

      /* Loop around the axis */
      int counter_inner = 0;
      while (next_t_idx != first_t_idx) {
        if (counter_inner++ >= VORONOI_CONSTRUCTION_MAX_NGB_ITER)
          error(
              "Voronoi face construction did not finish within allowed number "
              "of iterations! This probably means something went wrong with "
              "the Delaunay construction.");
#ifdef VORONOI_STORE_FACES
        if (face_vertices_index + 6 > face_vertices_size) {
          face_vertices_size <<= 1;
          face_vertices = (double *)realloc(
              face_vertices, 3 * face_vertices_size * sizeof(double));
        }
#endif
        /* Get the coordinates of the voronoi vertex corresponding to cur_t and
         * next_t */
        const double *vor_vertex1 = d->tetrahedra[cur_t_idx].circumcenter;
        const double *vor_vertex2 = d->tetrahedra[next_t_idx].circumcenter;
#ifdef VORONOI_STORE_FACES
        memcpy(&face_vertices[3 * face_vertices_index], vor_vertex1,
               3 * sizeof(*vor_vertex1));
        memcpy(&face_vertices[3 * face_vertices_index + 3], vor_vertex2,
               3 * sizeof(*vor_vertex2));
        face_vertices_index += 2;
#endif

        /* Update cell volume and tetrahedron_centroid */
        double tetrahedron_centroid[3] = {0., 0., 0.};
        const double V = geometry3d_compute_centroid_volume_tetrahedron(
            ax, ay, az, vor_vertex0[0], vor_vertex0[1], vor_vertex0[2],
            vor_vertex1[0], vor_vertex1[1], vor_vertex1[2], vor_vertex2[0],
            vor_vertex2[1], vor_vertex2[2], tetrahedron_centroid);
        volume += V;
        centroid[0] += V * tetrahedron_centroid[0];
        centroid[1] += V * tetrahedron_centroid[1];
        centroid[2] += V * tetrahedron_centroid[2];
        voronoi_assert(V >= 0.);

        /* Update face area and centroid */
        double triangle_area = geometry3d_compute_area_triangle(
            vor_vertex0[0], vor_vertex0[1], vor_vertex0[2], vor_vertex1[0],
            vor_vertex1[1], vor_vertex1[2], vor_vertex2[0], vor_vertex2[1],
            vor_vertex2[2]);
        double triangle_centroid[3];
        geometry3d_compute_centroid_triangle(
            vor_vertex0[0], vor_vertex0[1], vor_vertex0[2], vor_vertex1[0],
            vor_vertex1[1], vor_vertex1[2], vor_vertex2[0], vor_vertex2[1],
            vor_vertex2[2], triangle_centroid);
        face_area += triangle_area;
        face_centroid[0] += triangle_area * triangle_centroid[0];
        face_centroid[1] += triangle_area * triangle_centroid[1];
        face_centroid[2] += triangle_area * triangle_centroid[2];

        /* Update variables */
        prev_t_idx_in_cur_t = cur_t->index_in_neighbour[next_t_idx_in_cur_t];
        cur_t_idx = next_t_idx;
        cur_t = &d->tetrahedra[cur_t_idx];
        next_t_idx_in_cur_t = (prev_t_idx_in_cur_t + 1) % 4;
        while (cur_t->vertices[next_t_idx_in_cur_t] == gen_idx_in_d ||
               cur_t->vertices[next_t_idx_in_cur_t] == axis_idx_in_d) {
          next_t_idx_in_cur_t = (next_t_idx_in_cur_t + 1) % 4;
        }
        next_t_idx = cur_t->neighbours[next_t_idx_in_cur_t];
        /* Get the next non axis vertex and add it to the queue if necessary */
        next_non_axis_idx_in_d = cur_t->vertices[next_t_idx_in_cur_t];
        if (!neighbour_flags[next_non_axis_idx_in_d]) {
          int3 new_info = {._0 = cur_t_idx,
                           ._1 = next_non_axis_idx_in_d,
                           ._2 = next_t_idx_in_cur_t};
          int3_fifo_queue_push(&neighbour_info_q, new_info);
          neighbour_flags[next_non_axis_idx_in_d] |= 1;
        }
      }
      /* Finalize the face centroid */
      voronoi_assert(face_area >= 0.);
      if (face_area > 0.) {
        face_centroid[0] /= face_area;
        face_centroid[1] /= face_area;
        face_centroid[2] /= face_area;
      }
      if (voronoi_new_face(v, d, gen_idx_in_d, axis_idx_in_d, parts, face_area,
                           face_centroid, face_vertices, face_vertices_index)) {
#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
        /* The face is not degenerate */
        nface++;
#endif
        /* Update the minimal dist to a neighbouring generator */
        double ngb_pos[3];
        delaunay_get_vertex_at(d, axis_idx_in_d, ngb_pos);
        double dx[3] = {ngb_pos[0] - generator_pos[0],
                        ngb_pos[1] - generator_pos[1],
                        ngb_pos[2] - generator_pos[2]};
        double dist = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
        if (dist < min_ngb_dist2) {
          min_ngb_dist2 = dist;
          min_ngb_dist_pos[0] = ngb_pos[0];
          min_ngb_dist_pos[1] = ngb_pos[1];
          min_ngb_dist_pos[2] = ngb_pos[2];
        }
      }
    }

    /* Compute the actual centroid. */
    voronoi_assert(volume > 0.);
    centroid[0] /= volume;
    centroid[1] /= volume;
    centroid[2] /= volume;

    /* Estimate distance from centroid to nearest (to generator) face */
    double face[3] = {0.5 * (min_ngb_dist_pos[0] + p->x[0]),
                      0.5 * (min_ngb_dist_pos[1] + p->x[1]),
                      0.5 * (min_ngb_dist_pos[2] + p->x[2])};
    double dx_gen[3] = {face[0] - p->x[0], face[1] - p->x[1],
                        face[2] - p->x[2]};
    double dx_cen[3] = {face[0] - centroid[0], face[1] - centroid[1],
                        face[2] - centroid[2]};
    double dist = (dx_cen[0] * dx_gen[0] + dx_cen[1] * dx_gen[1] +
                   dx_cen[2] * dx_gen[2]) /
                  (0.5 * sqrt(min_ngb_dist2));

    p->geometry.volume = (float)volume;
    p->geometry.centroid[0] = (float)(centroid[0] - p->x[0]);
    p->geometry.centroid[1] = (float)(centroid[1] - p->x[1]);
    p->geometry.centroid[2] = (float)(centroid[2] - p->x[2]);
    p->geometry.nface = nface;
    p->geometry.pair_connections_offset = pair_connections_offset;
    p->geometry.min_face_dist = (float)dist;

    /* reset flags for all neighbours of this cell */
    neighbour_flags[gen_idx_in_d] = 0;
    for (int j = 0; j < neighbour_info_q.end; j++) {
      voronoi_assert(neighbour_info_q.values[j]._1 < d->vertex_index);
      neighbour_flags[neighbour_info_q.values[j]._1] = 0;
    }
#ifdef VORONOI_CHECKS
    for (int i = 0; i < d->vertex_index; i++) {
      voronoi_assert(neighbour_flags[i] == 0);
    }
#endif
  }
#ifdef VORONOI_STORE_FACES
  free(voronoi_vertices);
#endif
  free(neighbour_flags);
  int3_fifo_queue_destroy(&neighbour_info_q);
  voronoi_finalize(v, d, parts);
  voronoi_check_grid(v, d, parts);
  free(face_vertices);
}

/**
 * @brief Sort the faces according to their sid and set the face sid offsets and
 * recompute the cell face connections if necessary.
 *
 * @param v The voronoi struct
 * @param d The delaunay tesselation from which this voronoi was built
 * @param parts The array of local particles.
 */
inline static void voronoi_finalize(struct voronoi *v, const struct delaunay *d,
                                    struct part *parts) {
  /* Set some trackers */
  int counts[28];
  int offsets[28];
  int offset = 0;
  for (int i = 0; i < 28; i++) {
    int sid = face_sid_order[i];
    counts[sid] = 0;
    offsets[sid] = offset;
    /* Also set the (to be) correct pointers for the sids */
    v->pairs[sid] = &v->pairs_flat[offset];
    offset += v->pair_count[sid];
  }

  /* Loop over the faces and swap until they are in the right order */
  int face_idx = 0;
  int i = 0;
  int sid = 13;
  while (face_idx < v->pair_index) {
    struct voronoi_pair *face = &v->pairs_flat[face_idx];
    int face_sid = face->sid;
    if (face->sid & 1 << 5) {
      face_sid = 27;
      face->sid &= ~(1 << 5);
    }
    if (face_sid == sid) {
      /* Face at correct location: continue */
      counts[sid]++;
      voronoi_assert(counts[sid] <= v->pair_count[sid]);
      face_idx++;
      /* Go to next sid? */
      while (i < 27 && face_idx == offsets[sid] + v->pair_count[sid]) {
        i++;
        sid = face_sid_order[i];
        face_idx = offsets[sid] + counts[sid];
      }
      voronoi_assert(face_idx >= 0 && face_idx <= offset);
    } else {
      /* Swap this face to the correct sid */
      voronoi_assert(face_sid_index[face_sid] > face_sid_index[sid]);
      int idx_dest = offsets[face_sid] + counts[face_sid];
      voronoi_assert(idx_dest < v->pair_index);
      struct voronoi_pair tmp = v->pairs_flat[idx_dest];
      memcpy(&v->pairs_flat[idx_dest], &v->pairs_flat[face_idx],
             sizeof(struct voronoi_pair));
      memcpy(&v->pairs_flat[face_idx], &tmp, sizeof(struct voronoi_pair));
      counts[face_sid]++;
      voronoi_assert(counts[face_sid] <= v->pair_count[face_sid]);
      voronoi_assert(v->pairs_flat[face_idx].sid < 59);
    }
  }

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /* Recompute the cell face links */
#ifdef SWIFT_DEBUG_CHECKS
  int *face_counts_old = (int *)malloc((d->vertex_end - d->vertex_start) *
                                       sizeof(*face_counts_old));
#endif
  for (int j = d->vertex_start; j < d->vertex_end; j++) {
    struct part *p = &parts[d->vertex_part_idx[j]];
#ifdef SWIFT_DEBUG_CHECKS
    face_counts_old[j - d->vertex_start] = p->geometry.nface;
#endif
    p->geometry.nface = 0;
  }

  int idx = 0;
  /* First treat internal faces */
  voronoi_assert(face_sid_order[0] == 13);
  for (; idx < v->pair_count[13]; idx++) {
    const struct voronoi_pair *face = &v->pairs_flat[idx];
    struct part *left = &parts[face->left_idx];
    struct part *right = &parts[face->right_idx];

    /* Left part is always active */
    voronoi_assert(left->geometry.delaunay_vertex >= d->vertex_start &&
                   left->geometry.delaunay_vertex < d->vertex_end);
    v->cell_pair_connections[left->geometry.pair_connections_offset +
                             left->geometry.nface] = idx;
    left->geometry.nface++;

    /* Only update right face links if active */
    if (right->geometry.delaunay_vertex >= 0) {
      v->cell_pair_connections[right->geometry.pair_connections_offset +
                               right->geometry.nface] = idx;
      right->geometry.nface++;
    }
  }
  /* Now do the external faces (only left part needs updating */
  for (; idx < v->pair_index; idx++) {
    const struct voronoi_pair *face = &v->pairs_flat[idx];
    struct part *left = &parts[face->left_idx];

    /* Left part is always active */
    voronoi_assert(left->geometry.delaunay_vertex >= d->vertex_start &&
                   left->geometry.delaunay_vertex < d->vertex_end);
    v->cell_pair_connections[left->geometry.pair_connections_offset +
                             left->geometry.nface] = idx;
    left->geometry.nface++;
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int j = 0; j < 28; j++) {
    if (counts[j] != v->pair_count[j])
      error("Incorrect face count %i for sid %i! Should be %i", counts[j], j,
            v->pair_count[j]);
    int face_count = v->pair_count[j];
    const struct voronoi_pair *faces = v->pairs[j];
    assert(face_count == counts[j]);
    for (int k = 0; k < face_count; k++)
      if (faces->sid != j) error("Face ended up under wrong sid!");
  }
  /* Check whether the face counts match */
  int count_total = 0;
  for (int j = d->vertex_start; j < d->vertex_end; j++) {
    int counts_old = face_counts_old[j - d->vertex_start];
    int p_idx = d->vertex_part_idx[j];
    struct part *p = &parts[p_idx];
    int counts_new = p->geometry.nface;
    if (counts_new != counts_old) error("Face counts do not match!");
    count_total += counts_new;
  }
  if (count_total != v->cell_pair_connections_index)
    error("Total number of cell face connections does not match!");
  free(face_counts_old);
#endif
#endif
}

/**
 * @brief Free up all memory used by the Voronoi grid.
 *
 * @param v Voronoi grid.
 */
inline static void voronoi_destroy(struct voronoi *restrict v) {
  /* Anything allocated? */
  if (!v->active) {
    return;
  }

  if (v->pairs_flat != NULL) swift_free("voronoi", v->pairs_flat);
  v->pair_size = 0;
  v->pair_index = 0;
  for (int i = 0; i < 28; i++) {
    if (v->pairs[i] != NULL) {
      v->pairs[i] = NULL;
      v->pair_count[i] = 0;
    }
  }

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  if (v->cell_pair_connections != NULL)
    swift_free("voronoi", v->cell_pair_connections);
  v->cell_pair_connections_size = 0;
  v->cell_pair_connections_index = 0;
#endif

#ifdef VORONOI_STORE_FACES
  swift_free("voronoi", v->face_vertices);
  v->face_vertices = NULL;
  v->face_vertices_size = 0;
  v->face_vertices_index = 0;
#endif

  v->active = 0;
  v->min_surface_area = -1.;

  free(v);
}

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
inline static void voronoi_add_cell_face_connection(struct voronoi *v,
                                                    int face_idx) {

  if (v->cell_pair_connections_index == v->cell_pair_connections_size) {
    v->cell_pair_connections_size = (3 * v->cell_pair_connections_size) / 2;
    v->cell_pair_connections = (int *)swift_realloc(
        "voronoi", v->cell_pair_connections,
        v->cell_pair_connections_size * sizeof(*v->cell_pair_connections));
  }
  v->cell_pair_connections[v->cell_pair_connections_index++] = face_idx;
}
#endif

/**
 * @brief Add a two particle pair to the grid.
 *
 * This function also adds the correct tuple to the cell_pair_connections queue.
 *
 * The grid connectivity is stored per cell sid. We use the same convention as
 * in `sort_part.h`; i.e.: sid=13 corresponds to particle pairs encountered
 * during a self task (both particles are within the local cell), while sid=0-12
 * and sid 14-26 correspond to particle interactions for which the right
 * neighbour is part of one of the 26 neighbouring cells.  Additionally, any
 * boundary faces are stored under "fictive sid 27".
 *
 * For each pair, we compute and store all the quantities required to compute
 * fluxes between the Voronoi cells: the surface area and midpoint of the
 * interface.
 *
 * @param v Voronoi grid.
 * @param d Delaunay tesselation, dual of this Voronoi grid.
 * @param left_part_idx_in_d The index in the delaunay tesselation of the left
 * vertex of the new pair.
 * @param right_part_idx_in_d The index in the delaunay tesselation of the right
 * vertex of the new pair.
 * @param parts Particle array of local particles.
 * @param area Surface area of the new face.
 * @param centroid Centroid of the new face.
 * @param vertices Corner vertices of the new face.
 * @param n_vertices Number of vertices in the vertices array.
 * @returns 1 if a non-degenerate face was added or found, else 0
 */
inline static int voronoi_new_face(struct voronoi *restrict v,
                                   const struct delaunay *restrict d,
                                   int left_part_idx_in_d,
                                   int right_part_idx_in_d,
                                   struct part *restrict parts, double area,
                                   double *restrict centroid,
                                   double *restrict vertices, int n_vertices) {
  int sid;
  int left_part_idx = d->vertex_part_idx[left_part_idx_in_d];
  int right_part_idx = d->vertex_part_idx[right_part_idx_in_d];

  /* Pair between local active particles? */
  if (right_part_idx_in_d < d->vertex_end) {
    /* Already processed this pair? */
    if (right_part_idx_in_d < left_part_idx_in_d) {
#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
      /* Find the existing pair and add it to the cell_pair_connections. */
      struct part *ngb = &parts[right_part_idx];
      for (int i = 0; i < ngb->geometry.nface; i++) {
        int idx = ngb->geometry.pair_connections_offset + i;
        int face_idx = v->cell_pair_connections[idx];
        voronoi_assert(face_idx < v->pair_index);
        struct voronoi_pair *face = &v->pairs_flat[face_idx];
        if (face->sid == 13 && face->right_idx == left_part_idx) {
          voronoi_add_cell_face_connection(v, face_idx);
          return 1;
        }
      }
      /* If no pair is found, the face must have been degenerate, nothing left
       * to do. */
      return 0;
#else
      return 1;
#endif
    }
    sid = 13;
  } else {
    sid = d->ghost_cell_sids[right_part_idx_in_d - d->vertex_end];
  }

  /* Boundary particle? */
  int actual_sid = sid;
  if (sid & 1 << 5) {
    /* We store all boundary faces under fictive sid 27 */
    sid = 27;
  }

  /* Do we need to extend the pairs array? */
  if (v->pair_index == v->pair_size) {
    v->pair_size = (3 * v->pair_size) / 2;
    v->pairs_flat = (struct voronoi_pair *)swift_realloc(
        "voronoi", v->pairs_flat, v->pair_size * sizeof(struct voronoi_pair));
  }

  /* Grab the next free pair */
  int face_idx = v->pair_index;
  struct voronoi_pair *this_pair = &v->pairs_flat[face_idx];
  this_pair->surface_area = area;
  memcpy(this_pair->midpoint, centroid, 3 * sizeof(this_pair->midpoint[0]));
  /* is the face degenerate? */
  if (this_pair->surface_area < v->min_surface_area) {
    return 0;
  }
  v->pair_index++;

  /* Initialize pair */
  this_pair->left_idx = left_part_idx;
  this_pair->right_idx = right_part_idx;
  this_pair->sid = actual_sid;

#ifdef VORONOI_STORE_FACES
#ifdef SWIFT_DEBUG_CHECKS
  assert(n_vertices > 0);
#endif
  /* Enough space to store new faces? */
  int need_realloc = 0;
  while (v->face_vertices_index + n_vertices >= v->face_vertices_size) {
    v->face_vertices_size <<= 1;
    need_realloc = 1;
  }
  if (need_realloc) {
    v->face_vertices =
        realloc(v->face_vertices, 3 * v->face_vertices_size * sizeof(double));
  }
  this_pair->vertices = &v->face_vertices[3 * v->face_vertices_index];
  this_pair->n_vertices = n_vertices;
  memcpy(this_pair->vertices, vertices, 3 * n_vertices * sizeof(double));
#endif

  voronoi_assert(this_pair->sid < 59);

#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /* Add cell_pair_connection */
  voronoi_add_cell_face_connection(v, face_idx);
#endif

  /* increase index (to signal that we added the face) */
  v->pair_count[sid]++;

  return 1;
}

static inline double voronoi_compute_volume(const struct voronoi *restrict v) {
  double total_volume = 0.;
  /* TODO */
  return total_volume;
}

/**
 * @brief Sanity checks on the grid.
 *
 * Right now, this only checks the total volume of the cells.
 */
inline static void voronoi_check_grid(struct voronoi *v,
                                      const struct delaunay *d,
                                      const struct part *parts) {
#ifdef SWIFT_DEBUG_CHECKS
#ifdef VORONOI_STORE_CELL_FACE_CONNECTIONS
  /* Check cell - face connections */
  int count_total = 0;
  for (int i = d->vertex_start; i < d->vertex_end; i++) {
    int p_idx = d->vertex_part_idx[i];
    int face_connections_offset = parts[p_idx].geometry.pair_connections_offset;
    int nface = parts[p_idx].geometry.nface;
    count_total += nface;
    for (int k = face_connections_offset; k < face_connections_offset + nface;
         k++) {
      const struct voronoi_pair *face =
          &v->pairs_flat[v->cell_pair_connections[k]];
      if (face->sid == 13) {
        if (face->right_idx != p_idx && face->left_idx != p_idx) {
          error("Incorrect cell face link!");
        }
      } else if (face->left_idx != p_idx) {
        error("Incorrect cell face link!");
      }
    }
  }
  if (count_total != v->cell_pair_connections_index)
    error("total face connections count is wrong!");
#endif

  /* Check total volume */
  //  double total_volume = 0.;
  //  for (int j = d->vertex_start; j < d->vertex_end; j++) {
  //    total_volume += parts[d->vertex_part_idx[j]].geometry.volume;
  //  }
  //  fprintf(stderr, "Total volume: %g\n", total_volume);

  /* For each cell check that the total surface area is not bigger than the
   * surface area of a sphere with the same volume */
  double *surface_areas =
      malloc((d->vertex_end - d->vertex_start) * sizeof(double));
  for (int i = 0; i < (d->vertex_end - d->vertex_start); i++)
    surface_areas[i] = 0;

  for (int sid = 0; sid < 28; sid++)
    for (int i = 0; i < v->pair_count[sid]; i++) {
      struct voronoi_pair *pair = &v->pairs[sid][i];
      int left_idx =
          parts[pair->left_idx].geometry.delaunay_vertex - d->vertex_start;
      surface_areas[left_idx] += pair->surface_area;
      if (sid == 13 && parts[pair->right_idx].geometry.delaunay_vertex >= 0) {
        int right_idx =
            parts[pair->right_idx].geometry.delaunay_vertex - d->vertex_start;
        surface_areas[right_idx] += pair->surface_area;
      }
    }

  for (int i = 0; i < (d->vertex_end - d->vertex_start); i++) {
    float volume =
        parts[d->vertex_part_idx[i + d->vertex_start]].geometry.volume;
    double sphere_surface_area =
        4. * M_PI * pow(3. * volume / (4. * M_PI), 2. / 3.);
    voronoi_assert(sphere_surface_area < surface_areas[i]);
  }
  free(surface_areas);

#if defined(VORONOI_STORE_GENERATORS)
  /* Check connectivity */
  for (int i = 0; i < v->number_of_cells; i++) {
    struct voronoi_cell *this_cell = &v->cells[i];
    for (int j = this_cell->pair_connections_offset;
         j < this_cell->pair_connections_offset + this_cell->nface; j++) {
      int2 connection = v->cell_pair_connections.values[j];
      int pair_idx = connection._0;
      int sid = connection._1;
      struct voronoi_pair *pair = &v->pairs[sid][pair_idx];
      assert(i == pair->left_idx || (sid == 13 && i == pair->right_idx));
    }
  }
#endif

#endif
}

/**
 * @brief Write the Voronoi grid information to the given file.
 *
 * The output depends on the configuration. The maximal output contains 3
 * different types of output lines:
 *  - "G\tgx\tgx: x and y position of a single grid generator (optional).
 *  - "C\tcx\tcy\tV\tnface": centroid position, volume and (optionally) number
 *    of faces for a single Voronoi cell.
 *  - "F\tsid\tarea\tcx\tcy\tcz\t(v0x, v0y, v0z)\t...\t(vNx, vNy, vNz)": sid,
 *    area, coordinates of centroid, coordinates of vertices of face (optional).
 *
 * @param v Voronoi grid.
 * @param file File to write to.
 */
inline static void voronoi_write_grid(const struct voronoi *restrict v,
                                      const struct part *parts, int count,
                                      FILE *file, size_t *offset) {
  /* write the faces */
  for (int sid = 0; sid < 28; ++sid) {
    for (int i = 0; i < v->pair_count[sid]; ++i) {
      struct voronoi_pair *pair = &v->pairs[sid][i];
      fprintf(file, "F\t%i\t%g\t%g\t%g\t%g", sid, pair->surface_area,
              pair->midpoint[0], pair->midpoint[1], pair->midpoint[2]);
#ifdef VORONOI_STORE_FACES
      for (int j = 0; j < pair->n_vertices; j++) {
        fprintf(file, "\t(%g, %g, %g)", pair->vertices[3 * j],
                pair->vertices[3 * j + 1], pair->vertices[3 * j + 2]);
      }
#endif
      fprintf(file, "\n");
    }
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_VORONOI_3D_H
