//
// Created by yuyttenh on 14/06/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_BVH_H
#define SWIFTSIM_SHADOWSWIFT_BVH_H

#include "hydro.h"

#include <float.h>

typedef struct bbox {
  double anchor[3];
  double opposite[3];
} bbox_t;

inline static int bbox_contains(const bbox_t *bbox, double x, double y,
                                double z) {
  return bbox->anchor[0] <= x && x < bbox->opposite[0] &&
         bbox->anchor[1] <= y && y < bbox->opposite[1] &&
         bbox->anchor[2] <= z && z < bbox->opposite[2];
}

inline static void bbox_wrap(const bbox_t *bbox1, const bbox_t *bbox2, bbox_t *result) {
  result->anchor[0] = min(bbox1->anchor[0], bbox2->anchor[0]);
  result->anchor[1] = min(bbox1->anchor[1], bbox2->anchor[1]);
  result->anchor[2] = min(bbox1->anchor[2], bbox2->anchor[2]);
  result->opposite[0] = max(bbox1->opposite[0], bbox2->opposite[0]);
  result->opposite[1] = max(bbox1->opposite[1], bbox2->opposite[1]);
  result->opposite[2] = max(bbox1->opposite[2], bbox2->opposite[2]);
}

inline static void bbox_clip(const bbox_t *self, const double loc[3],
                             double out[3]) {
  out[0] = fmin(self->opposite[0], fmax(self->anchor[0], loc[0]));
  out[1] = fmin(self->opposite[1], fmax(self->anchor[1], loc[1]));
  out[2] = fmin(self->opposite[2], fmax(self->anchor[2], loc[2]));
}

inline static double bbox_distance2(const bbox_t *self, const double loc[3]) {
  double clipped[3];
  bbox_clip(self, loc, clipped);
  clipped[0] -= loc[0];
  clipped[1] -= loc[1];
  clipped[2] -= loc[2];
  return clipped[0] * clipped[0] + clipped[1] * clipped[1] +
         clipped[2] * clipped[2];
}

inline static bbox_t bbox_from_parts(const struct part *parts, const int *pid,
                                     int count) {
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
}

enum direction { X_axis, Y_axis, Z_axis };

#define BVH_DATA_SIZE 3

struct flat_bvh_node {
  bbox_t bbox;
  union {
    struct {
      int left;
      int right;
    } children;
    int data[BVH_DATA_SIZE + 1];
  };
  uint8_t is_leaf;
};

struct flat_bvh {
  struct flat_bvh_node *nodes;
  int size;
  int count;
};

inline static struct flat_bvh *flat_bvh_malloc(int size) {
  struct flat_bvh *bvh = malloc(sizeof(*bvh));
  bvh->nodes = (struct flat_bvh_node *)malloc(size * sizeof(*bvh->nodes));
  bvh->size = size;
  bvh->count = 0;
  return bvh;
}

inline static void flat_bvh_reset(struct flat_bvh *bvh) {
  bvh->count = 0;
}

inline static void flat_bvh_destroy(struct flat_bvh *bvh) {
  free(bvh->nodes);
  free(bvh);
}

inline static int flat_bvh_new_node(struct flat_bvh *bvh) {
  if (bvh->count == bvh->size) {
    bvh->size <<= 1;
    bvh->nodes = (struct flat_bvh_node *)realloc(
        bvh->nodes, bvh->size * sizeof(*bvh->nodes));
  }
  return bvh->count++;
}

void flat_bvh_populate_rec(struct flat_bvh *bvh, int node_id, const struct part *parts,
                           int *pid, int count);

inline static void flat_bvh_populate(struct flat_bvh *bvh, struct part *parts, int *pid, int count) {
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
 * @brief Finds a particle from this bvh that contains a given position in its
 * search radius (recursive version).
 *
 * Particles that are further away than a safety radius `h` are discarded
 *
 * @param bvh The #flat_bvh to search
 * @param node_id The node to start searching from (initially the root).
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @param r2 The square of the safety radius.
 */
int flat_bvh_hit_rec(const struct flat_bvh *bvh, int node_id,
                     struct part *parts, double x, double y, double z, double r2);

/**
 * @brief Finds a particle from this bvh that contains a given position in its
 * search radius.
 *
 * Particles that are further away than a safety radius `h` are discarded
 *
 * @param bvh The #flat_bvh to search
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @param r2 The square of the safety radius.
 */
inline static int flat_bvh_hit(const struct flat_bvh *bvh, struct part *parts,
                               double x, double y, double z, double r2) {
  return flat_bvh_hit_rec(bvh, 0, parts, x, y, z, r2);
}

inline static void flat_bvh_get_anchor(const struct flat_bvh *bvh, double *anchor) {
  anchor[0] = bvh->nodes[0].bbox.anchor[0];
  anchor[1] = bvh->nodes[0].bbox.anchor[1];
  anchor[2] = bvh->nodes[0].bbox.anchor[2];
}

inline static void flat_bvh_get_width(const struct flat_bvh *bvh, double *width) {
  struct flat_bvh_node *node = &bvh->nodes[0];
  width[0] = node->bbox.opposite[0] - node->bbox.anchor[0];
  width[1] = node->bbox.opposite[1] - node->bbox.anchor[1];
  width[2] = node->bbox.opposite[2] - node->bbox.anchor[2];
}

struct BVH {
  bbox_t bbox;

  struct BVH *left;

  struct BVH *right;

  int *data;

  int count;
};

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

void bvh_populate_rec_midpoint(struct BVH *bvh, const struct part *parts,
                               int *pid, int count);

void bvh_populate_rec(struct BVH *bvh, const struct part *parts,
                      double **coords, int *restrict pid, int count);

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

inline static void bvh_populate_midpoint(struct BVH *bvh,
                                         const struct part *parts,
                                         int *restrict pid, int count,
                                         int n_parts) {
  bvh_populate_rec_midpoint(bvh, parts, pid, count);
}

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

inline static int bvh_hit(const struct BVH *bvh, const struct part *parts,
                          double x, double y, double z) {
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
}

inline static void bvh_get_width(const struct BVH *bvh, double *width) {
  width[0] = bvh->bbox.opposite[0] - bvh->bbox.anchor[0];
  width[1] = bvh->bbox.opposite[1] - bvh->bbox.anchor[1];
  width[2] = bvh->bbox.opposite[2] - bvh->bbox.anchor[2];
}

#endif  // SWIFTSIM_SHADOWSWIFT_BVH_H
