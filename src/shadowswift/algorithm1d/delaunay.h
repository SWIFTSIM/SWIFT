//
// Created by yuyttenh on 9/12/22.
//

#ifndef SWIFTSIM_DELAUNAY_H
#define SWIFTSIM_DELAUNAY_H

struct line {
  /*! @brief Indices of the particles that make up the line. */
  int vertices[2];

  /*! @brief Indices of the neighbour lines. The neighbour at index 0 is the
   *  left neighbour and the neighbour at index 1 is the right neighbour */
  int neighbours[2];
};

inline static void line_set(struct line* l, int v0, int v1, int n0, int n1) {
  l->vertices[0] = v0;
  l->vertices[1] = v1;

  l->neighbours[0] = n0;
  l->neighbours[1] = n1;
}

struct delaunay {
  /*! @brief Anchor of the simulation volume. */
  double anchor;

  /*! @brief Side length of the delaunay tesselation. */
  double side;

  /*! @brief Vertex positions. This array is a copy of the array defined in
   *  main() and we probably want to get rid of it in a SWIFT implementation. */
  double* vertices;

  /*! @brief Flags that keep track whether a local vertex has been added or not
   */
  int* vertex_added;

  /*! @brief Link to the line containing this vertex as its left vertex. */
  int* vertex_line;

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

  /*! @brief lines that make up the tesselation */
  struct line* lines;

  /*! @brief Next available index in the lines array*/
  int line_index;

  /*! @brief Current allocated size of the lines array in memory */
  int line_size;

  /*! @brief Index of the last line that was accessed. Used as initial
   *  guess for the line that contains the next vertex that will be added.
   *  If vertices are added in some sensible order (e.g. in Peano-Hilbert curve
   *  order) then this will greatly speed up the algorithm. */
  int last_line;

  /*! @brief Cell neighbour sids keeping track of which neighbouring cells
   *  contains a specific neighbouring vertex. */
  int* ngb_cell_sids;

  /*! @brief Index of neighbouring particles in their respective cells. */
  int* ngb_part_idx;

  /*! @brief Current used size of the neighbouring vertex bookkeeping arrays
   *  (and next valid index in this array). */
  int ngb_index;

  /*! @brief Current size in memory of the neighbouring vertex bookkeeping
   *  arrays. More memory needs to be allocated if ngb_index reaches this
   *  value. */
  int ngb_size;

  /*! @brief Offset of the neighbouring vertices within the vertex array (so
   *  that neighbouring information for vertex v is stored in
   *  ngb_cell_sids[v-ngb_offset]). */
  int ngb_offset;

  /*! @brief Array of booleans indicating whether or not neighbouring particles
   * have been tried to be added for a given sid. If this is 0 for a given sid,
   * this means that this cell should get the reflective boundary condition
   * applied for that sid. */
  unsigned long int sid_is_inside_face_mask;
};

inline static int delaunay_add_vertex(struct delaunay* restrict d, int v,
                                      double x);
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x);
inline static void delaunay_init_vertex(struct delaunay* restrict d, int v,
                                        double x);
inline static int delaunay_new_line(struct delaunay* restrict d);
inline static void delaunay_check(const struct delaunay* d);

/**
 * @brief Reset the Delaunay tessellation without reallocating memory.
 *
 * It sets up a large line that contains the entire simulation box and
 * additional buffer space to deal with boundary ghost vertices, and 2
 * additional dummy lines that provide valid neighbours for the 2 sides of
 * this line (these dummy lines themselves have an invalid tip vertex
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
  d->vertex_index = vertex_size;
  d->line_index = 0;
  d->ngb_index = 0;

  /* Reset vertex_added flags */
  bzero(d->vertex_added, vertex_size * sizeof(int));

  /* Initialize box large enough to hold all local particles and ghost
   * particles from neighbouring cells */
  d->anchor = cell_loc[0] - cell_width[0];
  d->side = 3. * cell_width[0];

  /* set up the large triangle and the 3 dummies */
  /* mind the orientation: counterclockwise w.r.t. the z-axis. */
  int v0 = delaunay_new_vertex(d, d->anchor);
  delaunay_log("Creating vertex %i: %g %g", v0, box_anchor[0], box_anchor[1]);
  int v1 = delaunay_new_vertex(d, d->anchor + d->side);
  delaunay_log("Creating vertex %i: %g %g", v1, box_anchor[0] + box_side,
               box_anchor[1]);

  /* Now create a big initial line that will contain all possible generators
   * with 2 dummy lines with a fake tip as its neighbours */
  int dummy0 = delaunay_new_line(d);
  int dummy1 = delaunay_new_line(d);
  int line0 = delaunay_new_line(d);
  line_set(&d->lines[dummy0], -1, v0, -1, line0);
  line_set(&d->lines[line0], v0, v1, dummy0, dummy1);
  line_set(&d->lines[dummy1], v1, -1, line0, -1);
  d->vertex_line[v0] = line0;
  d->vertex_line[v1] = dummy1;

  /* Set up other parameters to valid values */
  d->last_line = line0;

  d->vertex_start = 0;
  d->vertex_end = vertex_size;

  d->ngb_offset = d->vertex_index;

  /* Initialize the sid mask so that the sid's with a z-component are already 1
   * (are all treated as internal), since this is the 2D version.
   * Also, sid=13 does not correspond to a face and is always set to 1 as well.
   * We only set the sid's corresponding to the cardinal directions to 0
   * (only faces perpendicular to one of the axes can be boundary faces). */
  d->sid_is_inside_face_mask = DEFAULT_SID_MASK;

  delaunay_check(d);
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
  vertex_size += 2;
  d->vertex_added = (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertices = (double*)swift_malloc("delaunay", vertex_size * sizeof(double));
  d->vertex_line = (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->vertex_size = vertex_size;

  /* allocate memory for the triangle array */
  d->line_size = vertex_size;
  d->lines = (struct line*)swift_malloc("delaunay",
                                        d->line_size * sizeof(struct line));

  d->ngb_cell_sids = (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->ngb_part_idx = (int*)swift_malloc("delaunay", vertex_size * sizeof(int));
  d->ngb_size = vertex_size;

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
  assert(d->active);
  assert(d->vertices != NULL);
#endif
  swift_free("delaunay", d->vertex_added);
  swift_free("delaunay", d->vertices);
  swift_free("delaunay", d->vertex_line);
  swift_free("delaunay", d->lines);
  swift_free("delaunay", d->ngb_cell_sids);
  swift_free("delaunay", d->ngb_part_idx);

  d->vertex_added = NULL;
  d->vertices = NULL;
  d->vertex_line = NULL;
  d->lines = NULL;
  d->ngb_cell_sids = NULL;
  d->ngb_part_idx = NULL;

  d->anchor = 0.f;
  d->side = 0;
  d->vertex_index = -1;
  d->vertex_size = 0;
  d->vertex_start = -1;
  d->vertex_end = -1;
  d->line_index = -1;
  d->line_size = 0;
  d->last_line = -1;
  d->ngb_index = -1;
  d->ngb_size = 0;
  d->ngb_offset = -1;
  d->sid_is_inside_face_mask = 0;

  /* Free delaunay struct itself */
  free(d);
}

/**
 * @brief Get the radius of the circumcircle of the given triangle.
 *
 * @param d Delaunay tessellation.
 * @param l Line index.
 * @return Radius of the circumcircle of the triangle.
 */
inline static double delaunay_get_radius(const struct delaunay* restrict d,
                                         int l) {
  int* vertices = d->lines[l].vertices;
  return 0.5 * (d->vertices[vertices[1]] - d->vertices[vertices[0]]);
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
                                             const int* restrict pid, int count,
                                             /*return*/ double* restrict r) {
  for (int i = 0; i < count; i++) {
    int idx = pid[i];
    int right_l_idx = d->vertex_line[idx];
    double radius = delaunay_get_radius(d, right_l_idx);
    /* Now also get the line connected to this generator on the left */
    int left_l_idx = d->lines[right_l_idx].neighbours[0];
    radius = fmax(radius, delaunay_get_radius(d, left_l_idx));
    r[i] = 2. * radius;
  }
}

inline static void delaunay_add_local_vertex(struct delaunay* restrict d, int v,
                                             double x, double y, double z,
                                             int ngb_idx) {
  /* Vertex already added? */
  v = v + d->vertex_start;
  if (d->vertex_added[v]) return;

  /* Update last triangle to be a triangle connected to the neighbouring
   * vertex. */
  if (ngb_idx >= 0) {
#ifdef SWIFT_DEBUG_CHECKS
    if (d->vertex_end - d->vertex_start <= ngb_idx)
      error("Invalid neighbour index passed to delaunay_add_new_vertex!");
#endif
    if (d->vertex_added[ngb_idx + d->vertex_start])
      d->last_line = d->vertex_line[ngb_idx + d->vertex_start];
  }

  delaunay_init_vertex(d, v, x);
  if (delaunay_add_vertex(d, v, x) == -1) {
    error("Local vertices cannot be added twice!");
  }
  d->vertex_added[v] = 1;
}

inline static void delaunay_add_new_vertex(struct delaunay* d, double x,
                                           double y, double z, int cell_sid,
                                           int part_idx, int ngb_idx,
                                           int is_boundary_particle) {
  /* Update last triangle to be a triangle connected to the neighbouring
   * vertex. */
  if (ngb_idx >= 0) {
#ifdef SWIFT_DEBUG_CHECKS
    if (d->vertex_end - d->vertex_start <= ngb_idx)
      error("Invalid neighbour index passed to delaunay_add_new_vertex!");
#endif
    if (d->vertex_added[ngb_idx]) d->last_line = d->vertex_line[ngb_idx];
  }

  /* create the new vertex */
  int v = delaunay_new_vertex(d, x);
  delaunay_log("Created new vertex with index %i", v);

  int flag = delaunay_add_vertex(d, v, x);
  if (flag == -1) {
    /* vertex already exists. Delete the last vertex. */
    --d->vertex_index;
  } else {
    /* vertex is new: add neighbour information for bookkeeping */
    if (d->ngb_index == d->ngb_size) {
      d->ngb_size <<= 1;
      d->ngb_cell_sids = (int*)swift_realloc("delaunay", d->ngb_cell_sids,
                                             d->ngb_size * sizeof(int));
      d->ngb_part_idx = (int*)swift_realloc("delaunay", d->ngb_part_idx,
                                            d->ngb_size * sizeof(int));
    }
    delaunay_assert(d->ngb_index == v - d->ngb_offset);
    /* Mark as boundary particle? */
    if (is_boundary_particle) cell_sid |= 1 << 5;
    /* Store info */
    d->ngb_cell_sids[d->ngb_index] = cell_sid;
    d->ngb_part_idx[d->ngb_index] = part_idx;
    ++d->ngb_index;
  }
}

/**
 * @brief Add a new vertex to the tessellation.
 *
 * This function finalizes the addition of a new vertex to the tessellation.
 * It locates which line contains the vertex and splits it into two lines.
 * No additional checks are needed in 1D.
 *
 * @param d Delaunay tessellation.
 * @param v Vertex index.
 * @param x Horizontal coordinate of the new vertex.
 * @return 0 on success, -1 if the vertex already exists.
 */
inline static int delaunay_add_vertex(struct delaunay* restrict d, int v,
                                      double x) {
  /* Find the line containing the new vertex */
  int l0_idx = d->last_line;
  while (1) {
    /* Make sure the candidate is valid */
    if (l0_idx < 2) {
      error("Trying to test a dummy line!");
    }
    struct line* candidate = &d->lines[l0_idx];

    if (x < d->vertices[candidate->vertices[0]]) {
      /* The vertex lies strictly on the left of this line */
      l0_idx = candidate->neighbours[0];
      continue;
    } else if (x > d->vertices[candidate->vertices[1]]) {
      /* The vertex lies strictly on the right of this line */
      l0_idx = candidate->neighbours[1];
      continue;
    } else {
      /* This vertex lies inside this line */
      /* Does it already exist? */
      if (d->vertices[candidate->vertices[0]] == x ||
          d->vertices[candidate->vertices[1]] == x)
        return -1;
      /* All is good */
      break;
    }
  }

  /* Split the line into two new lines */
  int new_l_idx = delaunay_new_line(d);
  struct line* l0 = &d->lines[l0_idx];
  struct line* new_l = &d->lines[new_l_idx];
  int v0 = l0->vertices[0];
  int v1 = l0->vertices[1];
  int n0 = l0->neighbours[0];
  int n1 = l0->neighbours[1];
  struct line* l1 = &d->lines[n1];
  line_set(l0, v0, v, n0, new_l_idx);
  line_set(new_l, v, v1, l0_idx, n1);
  /* The right neighbour of the new line must also update its left neighbour */
  line_set(l1, l1->vertices[0], l1->vertices[1], new_l_idx, l1->neighbours[1]);
  /* Also update vertex line links */
  d->vertex_line[v] = new_l_idx;

  delaunay_check(d);

  return 0;
}

/*! @brief stores the actual coordinates of the vertex at given by idx in out */
inline static void delaunay_get_vertex_at(const struct delaunay* d, int idx,
                                          double* out /*ret*/) {
  out[0] = d->vertices[idx];
}

inline static void delaunay_write_tessellation(
    const struct delaunay* restrict d, FILE* file, size_t* offset) {
  fprintf(file, "#VertexEnd\t%lu\tNeighbourOffset\t%lu\tnLines\t%d\n",
          *offset + d->vertex_end, *offset + d->ngb_offset, d->line_index - 2);

  for (int i = 0; i < d->vertex_index; ++i) {
    fprintf(file, "V\t%lu\t%g\n", *offset + i, d->vertices[i]);
  }
  for (int i = 2; i < d->line_index; ++i) {
    fprintf(file, "T\t%lu\t%lu\n", *offset + d->lines[i].vertices[0],
            *offset + d->lines[i].vertices[1]);
  }
  *offset += d->vertex_index;
}

/**
 * @brief Add a new vertex with the given coordinates.
 *
 * This function makes sure there is sufficient memory to store the
 * vertex and all its properties.
 * The function also initialises all vertex properties to sensible values.
 *
 * @param d Delaunay tessellation.
 * @param x Horizontal coordinate of the vertex.
 * @return Index of the new vertex within the vertex array.
 */
inline static int delaunay_new_vertex(struct delaunay* restrict d, double x) {
  /* check the size of the vertex arrays against the allocated memory size */
  if (d->vertex_index == d->vertex_size) {
    /* dynamically grow the size of the arrays with a factor 2 */
    d->vertex_size <<= 1;
    d->vertices = (double*)swift_realloc("delaunay", d->vertices,
                                         d->vertex_size * sizeof(double));
    d->vertex_line = (int*)swift_realloc("delaunay", d->vertex_line,
                                         d->vertex_size * sizeof(int));
  }

  delaunay_init_vertex(d, d->vertex_index, x);

  /* return the vertex index and then increase it by 1.
     After this operation, vertex_index will correspond to the size of the
     vertex arrays and is also the index of the next vertex that will be
     created. */
  return d->vertex_index++;
}

/**
 * @brief Initialize a vertex' properties
 *
 * @param d Delaunay tessellation
 * @param x The coordinate of the vertex
 */
inline static void delaunay_init_vertex(struct delaunay* restrict d, int v,
                                        double x) {
  d->vertices[v] = x;
  d->vertex_line[v] = -1;
}

/**
 * @brief Claim a new line in the line array.
 *
 * This function first ensures that the line array is still large enough
 * to hold all lines, and then returns the index of the next available
 * line.
 *
 * @param d Delaunay tessellation.
 * @return Index of the next available line in the line array.
 */
inline static int delaunay_new_line(struct delaunay* restrict d) {
  /* check that we still have triangles available */
  if (d->line_index == d->line_size) {
    /* no: increase the size of the line array with a factor 2 and
       reallocate it in memory */
    d->line_size <<= 1;
    d->lines = (struct line*)swift_realloc("delaunay", d->lines,
                                           d->line_size * sizeof(struct line));
  }

  /* return the line index and then increase it by 1.
     After this operation, line_index will correspond to the size of the
     line array and is also the index of the next line that will be
     created. */
  return d->line_index++;
}

inline static void delaunay_check(const struct delaunay* d) {
#ifdef DELAUNAY_CHECKS
  /* Check linkages */
  for (int i = 0; i < d->vertex_index; i++) {
    if (i < d->vertex_end && !d->vertex_added[i]) continue;
    if (i != d->lines[d->vertex_line[i]].vertices[0]) {
      error("Vertex-line links mixed up!");
    }
  }
  for (int i = 2; i < d->line_index; i++) {
    const struct line* l = &d->lines[i];
    const struct line* l0 = &d->lines[l->neighbours[0]];
    const struct line* l1 = &d->lines[l->neighbours[1]];
    if (l0->neighbours[1] != i) error("Line-line links mixed up!");
    if (l1->neighbours[0] != i) error("Line-line links mixed up!");
    if (l0->vertices[1] != l->vertices[0]) error("Line-vertex links mixed up!");
    if (l1->vertices[0] != l->vertices[1]) error("Line-vertex links mixed up!");
  }

  /* Check delaunay criterion */
  for (int i = 2; i < d->line_index; i++) {
    const struct line* l = &d->lines[i];
    for (int j = 0; j < d->vertex_index; j++) {
      if (j < d->vertex_end && !d->vertex_added[j]) continue;
      int v0 = l->vertices[0];
      int v1 = l->vertices[1];
      if (j == v0 || j == v1) continue;
      if (d->vertices[v0] < d->vertices[j] && d->vertices[j] < d->vertices[v1])
        error("Delaunay criterion violated!");
    }
  }
#endif
}
#endif  // SWIFTSIM_DELAUNAY_H
