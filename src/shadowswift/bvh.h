//
// Created by yuyttenh on 14/06/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_BVH_H
#define SWIFTSIM_SHADOWSWIFT_BVH_H

#include "hydro.h"

#include <float.h>

/** @file This file contains structs and functions used to construct and utilize
 * a BVH (bounding volume heap) of SWIFT hydro particles.
 *
 * Right now, the bounding volumes are constructed using the particles' search
 * radii, which is only defined in the shadowswift hydro scheme.
 **/

/** @brief A simple axes-aligned bounding box */
typedef struct bbox {
  /*!@brief The lower left corner. */
  double anchor[3];

  /*! @brief The corner diagonally opposite anchor. */
  double opposite[3];
} bbox_t;

/**
 * @brief Tests whether a given positions is contained in a given bounding box.
 *
 * @param bbox The bounding box to test
 * @param x, y, z The coordinate to test
 * @returns 1 when the coordinate is contained in the bounding box, else 0.
 **/
inline static int bbox_contains(const bbox_t *bbox, double x, double y,
                                double z) {
  return bbox->anchor[0] <= x && x < bbox->opposite[0] &&
         bbox->anchor[1] <= y && y < bbox->opposite[1] &&
         bbox->anchor[2] <= z && z < bbox->opposite[2];
}

/**
 * @brief Initialize a new bounding box so that it encloses two given bboxes
 *
 * @param bbox1, bbox2 The bboxes to wrap.
 * @param result (Return) The bbox which will be initialized by this function.
 **/
inline static void bbox_wrap(const bbox_t *bbox1, const bbox_t *bbox2,
                             bbox_t *result) {
  result->anchor[0] = min(bbox1->anchor[0], bbox2->anchor[0]);
  result->anchor[1] = min(bbox1->anchor[1], bbox2->anchor[1]);
  result->anchor[2] = min(bbox1->anchor[2], bbox2->anchor[2]);
  result->opposite[0] = max(bbox1->opposite[0], bbox2->opposite[0]);
  result->opposite[1] = max(bbox1->opposite[1], bbox2->opposite[1]);
  result->opposite[2] = max(bbox1->opposite[2], bbox2->opposite[2]);
}

/**
 * @brief Clip a coordinate to a bbox.
 *
 * The coordinate is clipped in all 3 directions by the extent of the bbox.
 *
 * @param self The bbox used for clipping
 * @param loc The coordinate to clip
 * @param out (Return) The clipped coordinate.
 **/
inline static void bbox_clip(const bbox_t *self, const double loc[3],
                             double out[3]) {
  out[0] = fmin(self->opposite[0], fmax(self->anchor[0], loc[0]));
  out[1] = fmin(self->opposite[1], fmax(self->anchor[1], loc[1]));
  out[2] = fmin(self->opposite[2], fmax(self->anchor[2], loc[2]));
}

/**
 * @brief Compute the square of the shortest possible distance from a given
 * coordinate to a bbox.
 *
 * This is equal to the distance between the given coordinate and that
 * coordinate clipped by the bbox.
 *
 * @param self The bbox
 * @param loc The coordinate
 * @return The squared distance.
 **/
inline static double bbox_distance2(const bbox_t *self, const double loc[3]) {
  double clipped[3];
  bbox_clip(self, loc, clipped);
  clipped[0] -= loc[0];
  clipped[1] -= loc[1];
  clipped[2] -= loc[2];
  return clipped[0] * clipped[0] + clipped[1] * clipped[1] +
         clipped[2] * clipped[2];
}

/**
 * @brief Construct a bbox around some hydro particles.
 *
 * The bbox will bound the
 * search radii of the particles (`part.geometry.search_radius`) and not the
 * smoothing length. This means that for now, this function is only compatible
 * with the shadowswift hydro scheme.
 *
 * @param parts The array containing the particles to wrap.
 * @param pid The indices of the particles we actually want to wrap in the
 * `parts` array.
 * @param count The length of the `pid` array.
 * @return A new bbox that wraps the particles.
 **/
inline static bbox_t bbox_from_parts(const struct part *parts, const int *pid,
                                     int count) {
#ifndef MOVING_MESH
  error("Should not be calling this function!");
#else
  double min_x = DBL_MAX;
  double max_x = -DBL_MAX;
  double min_y = DBL_MAX;
  double max_y = -DBL_MAX;
  double min_z = DBL_MAX;
  double max_z = -DBL_MAX;

  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    min_x = min(min_x, p->x[0] - p->geometry.search_radius);
    max_x = max(max_x, p->x[0] + p->geometry.search_radius);
    min_y = min(min_y, p->x[1] - p->geometry.search_radius);
    max_y = max(max_y, p->x[1] + p->geometry.search_radius);
    min_z = min(min_z, p->x[2] - p->geometry.search_radius);
    max_z = max(max_z, p->x[2] + p->geometry.search_radius);
  }

  bbox_t result = {.anchor = {min_x, min_y, min_z},
                   .opposite = {max_x, max_y, max_z}};
  return result;
#endif
}

/**
 * @brief Simple enum for the cartesian directions.
 **/
enum direction { X_axis, Y_axis, Z_axis };

/*! @brief The number of particles/elements to store in a BVH leaf (meaning
 * that we don't refine a node further once it has this number of elements or
 * less. */
#define BVH_DATA_SIZE 3

/**
 * @brief A node of a BVH.
 **/
struct flat_bvh_node {
  /*! @brief The bounding box of this node and all elements of its children */
  bbox_t bbox;

  union {
    /*! @brief The children of this node */
    struct {
      /*! @brief The left subtree */
      int right;
      /*! @brief The right subtree */
      int left;
    } children;

    /*! @brief The data/elements contained in this node. */
    int data[BVH_DATA_SIZE + 1];
  };

  /*! @brief Whether this node is a leaf (stores `data`) or an internal node
   * (has children). */
  uint8_t is_leaf;
};

/**
 * @brief A BVH (Bounding Volume Heap) struct. Nodes are stored contiguously in
 * memory in a flat array.
 *
 * This is the most efficient BVH implementation and should be used by default.
 **/
struct flat_bvh {
  /*! @brief The nodes of this BVH. */
  struct flat_bvh_node *nodes;

  /*! @brief The allocated size of the `nodes` array */
  int size;

  /*! @brief The actual number of nodes. */
  int count;
};

/**
 * @brief Allocates a BVH.
 *
 * @param size The number of nodes to reserve memory for.
 * @returns A pointer to the newly allocated BVH
 **/
inline static struct flat_bvh *flat_bvh_malloc(int size) {
  struct flat_bvh *bvh = malloc(sizeof(*bvh));
  bvh->nodes = (struct flat_bvh_node *)malloc(size * sizeof(*bvh->nodes));
  bvh->size = size;
  bvh->count = 0;
  return bvh;
}

/**
 * @brief Reset this BVH without freeing any memory. This effectively clears
 * the bvh.
 *
 * @param bvh The BVH to reset.
 **/
inline static void flat_bvh_reset(struct flat_bvh *bvh) { bvh->count = 0; }

/**
 * @brief Deallocate this BVH.
 *
 * @param bvh The BVH to destroy.
 **/
inline static void flat_bvh_destroy(struct flat_bvh *bvh) {
  free(bvh->nodes);
  free(bvh);
}

/**
 * @brief Add a new node to a BVH and return its index. This function ensures
 * that there is enough memory for the new node.
 *
 * @param bvh The BVH to add a node to.
 * @returns The index of the newly added, but uninitialized node.
 **/
inline static int flat_bvh_new_node(struct flat_bvh *bvh) {
  if (bvh->count == bvh->size) {
    bvh->size <<= 1;
    bvh->nodes = (struct flat_bvh_node *)realloc(
        bvh->nodes, bvh->size * sizeof(*bvh->nodes));
  }
  return bvh->count++;
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * @param bvh The BVH under construction
 * @param node_id The index of the root of the current subtree under
 * construction.
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current subtree.
 * @param count The length of the `pid` array.
 **/
void flat_bvh_populate_rec(struct flat_bvh *bvh, int node_id,
                           const struct part *parts, int *pid, int count);

/**
 * @brief Public interface for constructing a BVH from hydro particles.
 *
 * The search radii of the particles are bounded, not the smoothing lengths.
 * This makes this function incompatible with any hydro schemes other than
 * shadowswift.
 *
 * Nodes of the BVH are split at their midpoint in the direction in which their
 * width is largest. This produces a BVH of lesser quality than e.g. using a
 * median split, but is much faster to compute.
 *
 * @param bvh The BVH struct with some memory allocated.
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current subtree.
 * @param count The length of the `pid` array.
 **/
inline static void flat_bvh_populate(struct flat_bvh *bvh, struct part *parts,
                                     int *pid, int count) {
  flat_bvh_reset(bvh);
  int root = flat_bvh_new_node(bvh);
  flat_bvh_populate_rec(bvh, root, parts, pid, count);
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < bvh->count; i++) {
    struct flat_bvh_node *node = &bvh->nodes[i];
    if (node->is_leaf) {
      assert(node->data[BVH_DATA_SIZE] <= BVH_DATA_SIZE);
    } else {
      assert(node->children.left < bvh->count && 0 <= node->children.left);
      assert(node->children.right < bvh->count && 0 <= node->children.right);
    }
  }
#endif
}

/**
 * @brief Finds a particle from this bvh that contains a given candidate
 * position in its search radius.
 *
 * Particles that are further away than a safety radius `r` from the candidate
 * are discarded and not considered a hit.
 *
 * @param bvh The #flat_bvh to search
 * @param node_id The node to start searching from (initially the root).
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @param r2 The square of the safety radius.
 * @returns The index of the "hit", i.e.: the particle that contains the
 * candidate position in it's search radius and is itself contained in the
 * safety radius of the candidate position. If no hit is found, -1 is returned.
 */
int flat_bvh_hit_rec(const struct flat_bvh *bvh, int node_id,
                     struct part *parts, double x, double y, double z,
                     double r2);

/**
 * @brief Finds a particle from this bvh that contains a given candidate
 * position in its search radius. (public interface)
 *
 * Particles that are further away than a safety radius `r` from the candidate
 * are discarded and not considered a hit.
 *
 * @param bvh The #flat_bvh to search
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @param r2 The square of the safety radius.
 * @returns The index of the "hit", i.e.: the particle that contains the
 * candidate position in it's search radius and is itself contained in the
 * safety radius of the candidate position. If no hit is found, -1 is returned.
 */
inline static int flat_bvh_hit(const struct flat_bvh *bvh, struct part *parts,
                               double x, double y, double z, double r2) {
  return flat_bvh_hit_rec(bvh, 0, parts, x, y, z, r2);
}

/**
 * @brief Get the anchor of the bounding box of this BVH (i.e. of the root
 * node).
 *
 * @param bvh The BVH.
 * @param anchor (Return) Array to store the anchor in.
 **/
inline static void flat_bvh_get_anchor(const struct flat_bvh *bvh,
                                       double *anchor) {
  anchor[0] = bvh->nodes[0].bbox.anchor[0];
  anchor[1] = bvh->nodes[0].bbox.anchor[1];
  anchor[2] = bvh->nodes[0].bbox.anchor[2];
}

/**
 * @brief Get the 3 dimensional width of the bounding box of this BVH (i.e. of
 * the root node).
 *
 * @param bvh The BVH.
 * @param anchor (Return) Array to store the 3D width in.
 **/
inline static void flat_bvh_get_width(const struct flat_bvh *bvh,
                                      double *width) {
  struct flat_bvh_node *node = &bvh->nodes[0];
  width[0] = node->bbox.opposite[0] - node->bbox.anchor[0];
  width[1] = node->bbox.opposite[1] - node->bbox.anchor[1];
  width[2] = node->bbox.opposite[2] - node->bbox.anchor[2];
}

/**
 * @brief A BVH (Bounding Volume Heap) struct. Each node is allocated
 * independently.
 *
 * This implementation is not recommended due to poor cache efficiency, but is
 * still available as reference.
 **/
struct BVH {
  /*! @brief The bounding box of this BVH */
  bbox_t bbox;

  /*! @brief The left subtree of this BVH node (NULL for leafs) */
  struct BVH *left;

  /*! @brief The right subtree of this BVH node (NULL for leafs) */
  struct BVH *right;

  /*! @brief The data stored in this BVH node */
  int *data;

  /*! @brief The number of elements stored in the `data` array. */
  int count;
};

/**
 * @brief Free any data associated to this BVH.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The BVH to destroy.
 **/
inline static void bvh_destroy(struct BVH *bvh) {

  if (bvh->left != NULL) {
    bvh_destroy(bvh->left);
  }
  if (bvh->right != NULL) {
    bvh_destroy(bvh->right);
  }

  if (bvh->data != NULL) {
    free(bvh->data);
  }

  free(bvh);
}

/**
 * @brief Computes the ordering of 2 indices using some associated data (in this
 * case: coordinates).
 *
 * @param a Pointer to the index of the first element (will be cast to (int *))
 * @param b Pointer to the index of the second element (will be cast to (int *))
 * @param arg Pointer to the coordinates array (will be cast to (double *))
 * @returns -1 if arg[*a] < arg[*b] 0 if they are equal and 1 otherwise.
 **/
inline static int cmp(const void *a, const void *b, void *arg) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  double *coords = (double *)arg;
  double ad = coords[ai];
  double bd = coords[bi];

  if (ad < bd) {
    return -1;
  } else if (ad > bd) {
    return 1;
  } else {
    return 0;
  }
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * This method splits nodes by their median position along the direction of the
 * maximal width. This produces a well balanced BVH, but is slower to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param coords Array of 3 arrays containing the x, y and z coordinates of the
 * particles to be added respectively.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
void bvh_populate_rec(struct BVH *bvh, const struct part *parts,
                      double **coords, int *restrict pid, int count);

/**
 * @brief Construct a BVH from hydro particles (public interface).
 *
 * This method splits nodes by their median position along the direction of the
 * maximal width. This produces a well balanced BVH, but is slower to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
inline static void bvh_populate(struct BVH *bvh, const struct part *parts,
                                int *restrict pid, int count, int n_parts) {
  double **coords = malloc(3 * sizeof(double *));
  coords[0] = malloc(n_parts * sizeof(double));
  coords[1] = malloc(n_parts * sizeof(double));
  coords[2] = malloc(n_parts * sizeof(double));

  for (int i = 0; i < n_parts; i++) {
    coords[0][i] = parts[i].x[0];
    coords[1][i] = parts[i].x[1];
    coords[2][i] = parts[i].x[2];
  }

  bvh_populate_rec(bvh, parts, coords, pid, count);

  /* be clean */
  free(coords[0]);
  free(coords[1]);
  free(coords[2]);
  free(coords);
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * This method splits nodes by their midpoint, which produces a less balanced
 * BVH, but is faster to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
void bvh_populate_rec_midpoint(struct BVH *bvh, const struct part *parts,
                               int *pid, int count);

/**
 * @brief Construct a BVH from hydro particles (public interface).
 *
 * This method splits nodes by their midpoint, which produces a less balanced
 * BVH, but is faster to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
inline static void bvh_populate_midpoint(struct BVH *bvh,
                                         const struct part *parts,
                                         int *restrict pid, int count) {
  bvh_populate_rec_midpoint(bvh, parts, pid, count);
}

/**
 * @brief Is this BVH a leaf node?
 *
 * @param bvh The BVH.
 * @returns 1 if leaf 0 otherwise.
 **/
inline static int bvh_is_leaf(const struct BVH *bvh) {
#ifdef SWIFT_DEBUG_CHECKS
  if (bvh->left == NULL) {
    assert(bvh->right == NULL);
    assert(bvh->count > 0);
    assert(bvh->data != NULL);
  } else {
    assert(bvh->data == NULL);
    assert(bvh->count == 0);
  }
#endif
  return bvh->left == NULL;
}

/**
 * @brief Finds a particle from this bvh that contains a given candidate
 * position in its search radius. (public interface)
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The #flat_bvh to search
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @returns The index of the "hit", i.e.: the particle that contains the
 * candidate position in it's search radius. If no hit is found, -1 is returned.
 */
inline static int bvh_hit(const struct BVH *bvh, const struct part *parts,
                          double x, double y, double z) {
#ifndef MOVING_MESH
  error("Should not be calling this function!");
#else
  /* Anything to do here? */
  if (!bbox_contains(&bvh->bbox, x, y, z)) return -1;

  /* Are we in a leaf? */
  if (bvh_is_leaf(bvh)) {

    /* Check individual particles */
    for (int i = 0; i < bvh->count; i++) {
      int pid = bvh->data[i];
      const struct part *p = &parts[pid];
      double r2 = p->geometry.search_radius * p->geometry.search_radius;
      double dx[3] = {p->x[0] - x, p->x[1] - y, p->x[2] - z};
      double dx2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (dx2 <= r2) return pid;
    }

    /* no hit */
    return -1;
  }

  /* We are not in a leaf: recurse */
  int pid = bvh_hit(bvh->left, parts, x, y, z);
  if (pid >= 0) return pid;
  pid = bvh_hit(bvh->right, parts, x, y, z);
  if (pid >= 0) return pid;

  /* No hit */
  return -1;
#endif
}

/**
 * Get the width of this BVH.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The BVH.
 * @param width (Return) Array in which the width will be stored.
 **/
inline static void bvh_get_width(const struct BVH *bvh, double *width) {
  width[0] = bvh->bbox.opposite[0] - bvh->bbox.anchor[0];
  width[1] = bvh->bbox.opposite[1] - bvh->bbox.anchor[1];
  width[2] = bvh->bbox.opposite[2] - bvh->bbox.anchor[2];
}

#endif  // SWIFTSIM_SHADOWSWIFT_BVH_H
