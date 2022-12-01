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

inline static BBox bbox_from_coords(const double *coords[3],
                                    const double *search_radii, const int *idx,
                                    int count) {
  double min_x = DBL_MAX;
  double max_x = -DBL_MAX;
  double min_y = DBL_MAX;
  double max_y = -DBL_MAX;
  double min_z = DBL_MAX;
  double max_z = -DBL_MAX;

  for (int i = 0; i < count; i++) {
    const double x[3] = {coords[0][idx[i]], coords[1][idx[i]],
                         coords[2][idx[i]]};
    const double search_radius = search_radii[idx[i]];
    min_x = min(min_x, x[0] - search_radius);
    max_x = max(max_x, x[0] + search_radius);
    min_y = min(min_y, x[1] - search_radius);
    max_y = max(max_y, x[1] + search_radius);
    min_z = min(min_z, x[2] - search_radius);
    max_z = max(max_z, x[2] + search_radius);
  }

  BBox result = {.anchor = {min_x, min_y, min_z},
                 .opposite = {max_x, max_y, max_z}};
  return result;
}

enum direction { X_axis, Y_axis, Z_axis };

struct BVH {
  /*! Bounding box of all the particles in this BVH. This is the union of
   * the bounding boxes of the subtrees if any or jus the bounding box of the
   * particles in this leaf */
  BBox bbox;

  /*! Left subtree (if any) */
  struct BVH *left;

  /*! Right subtree (if any) */
  struct BVH *right;

  struct {
    /*! Array of particle indices (with respect to the cell for which this bvh
     * was constructed) */
    int *pid;
    /*! Array of particle search radii (as used for construction) */
    double *radius;
  } data;

  /*! The number of particles stored in this BVH and its subtrees */
  int count;
};

void bvh_destroy(struct BVH *bvh);

void bvh_clear(struct BVH *bvh);

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

void bvh_populate_rec(struct BVH *bvh, const double *coords[3],
                      const double *search_radii, const int *pid, int *idx,
                      int count);

inline static void bvh_populate(struct BVH *bvh, const struct part *parts,
                                double *search_radii, const int *pid,
                                int count) {
  double *coords[3] = {malloc(count * sizeof(*coords[0])),
                       malloc(count * sizeof(*coords[1])),
                       malloc(count * sizeof(*coords[2]))};
  int *idx = malloc(count * sizeof(*idx));
  int malloced_search_radii = 0;
  if (search_radii == NULL) {
    search_radii = malloc(count * sizeof(*search_radii));
    malloced_search_radii = 1;
  }

  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    coords[0][i] = p->x[0];
    coords[1][i] = p->x[1];
    coords[2][i] = p->x[2];
    idx[i] = i;
    if (malloced_search_radii) {
      search_radii[i] = p->h;
    }
  }

  bvh_populate_rec(bvh, (const double **)coords, search_radii, pid, idx, count);

  /* be clean */
  free(coords[0]);
  free(coords[1]);
  free(coords[2]);
  free(idx);
  if (malloced_search_radii) free(search_radii);
}

inline static int bvh_is_leaf(const struct BVH *bvh) {
#ifdef SWIFT_DEBUG_CHECKS
  if (bvh->left == NULL) {
    assert(bvh->right == NULL);
    assert(bvh->count > 0);
    assert(bvh->data.pid != NULL);
    assert(bvh->data.radius != NULL);
  } else {
    assert(bvh->data.pid == NULL);
    assert(bvh->data.radius == NULL);
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
      int pid = bvh->data.pid[i];
      const struct part *p = &parts[pid];
      double r2 = bvh->data.radius[i] * bvh->data.radius[i];
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
