//
// Created by yuyttenh on 14/06/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_BVH_H
#define SWIFTSIM_SHADOWSWIFT_BVH_H

#include "hydro.h"

#include <float.h>

typedef struct BBox {
  double anchor[3];
  double opposite[3];
} BBox;

inline static int bbox_contains(const BBox *bbox, double x, double y,
                                double z) {
  return bbox->anchor[0] <= x && x < bbox->opposite[0] &&
         bbox->anchor[1] <= y && y < bbox->opposite[1] &&
         bbox->anchor[2] <= z && z < bbox->opposite[2];
}

inline static BBox bbox_wrap(const BBox *bbox1, const BBox *bbox2) {
  BBox result = {.anchor = {min(bbox1->anchor[0], bbox2->anchor[0]),
                            min(bbox1->anchor[1], bbox2->anchor[1]),
                            min(bbox1->anchor[2], bbox2->anchor[2])},
                 .opposite = {max(bbox1->opposite[0], bbox2->opposite[0]),
                              max(bbox1->opposite[1], bbox2->opposite[1]),
                              max(bbox1->opposite[2], bbox2->opposite[2])}};
  return result;
}

inline static BBox bbox_from_parts(const struct part *parts, const int *pid,
                                   int count) {
  double min_x = DBL_MAX;
  double max_x = -DBL_MAX;
  double min_y = DBL_MAX;
  double max_y = -DBL_MAX;
  double min_z = DBL_MAX;
  double max_z = -DBL_MAX;

  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    min_x = min(min_x, p->x[0] - p->h);
    max_x = max(max_x, p->x[0] + p->h);
    min_y = min(min_y, p->x[1] - p->h);
    max_y = max(max_y, p->x[1] + p->h);
    min_z = min(min_z, p->x[2] - p->h);
    max_z = max(max_z, p->x[2] + p->h);
  }

  BBox result = {.anchor = {min_x, min_y, min_z},
                 .opposite = {max_x, max_y, max_z}};
  return result;
}

enum direction { X_axis, Y_axis, Z_axis };

struct BVH {
  BBox bbox;

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
      double r2 = p->h * p->h;
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

#endif  // SWIFTSIM_SHADOWSWIFT_BVH_H
