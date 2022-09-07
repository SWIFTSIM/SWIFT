#include "bvh.h"

#include <stdlib.h>
#include <string.h>

void bvh_populate_rec(struct BVH *bvh, const struct part *parts,
                      double **coords, int *restrict pid,
                      int count) {
  if (count <= 4) {
    /* Set unused fields of this bvh */
    bvh->left = NULL;
    bvh->right = NULL;

    /* Set used fields for leaf */
    size_t size = count * sizeof(int);
    bvh->data = malloc(size);
    memcpy(bvh->data, pid, size);
    bvh->count = count;
    bvh->bbox = bbox_from_parts(parts, pid, count);
    return;
  }

  /* Determine split_direction (direction with the largest spread) */
  double min_x = DBL_MAX;
  double max_x = -DBL_MAX;
  double min_y = DBL_MAX;
  double max_y = -DBL_MAX;
  double min_z = DBL_MAX;
  double max_z = -DBL_MAX;

  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    min_x = min(min_x, p->x[0]);
    max_x = max(max_x, p->x[0]);
    min_y = min(min_y, p->x[1]);
    max_y = max(max_y, p->x[1]);
    min_z = min(min_z, p->x[2]);
    max_z = max(max_z, p->x[2]);
  }

  enum direction split_direction;
  if (max_x - min_x >= max_y - min_y && max_x - min_x >= max_z - min_z) {
    split_direction = X_axis;
  } else if (max_y - min_y >= max_x - min_x && max_y - min_y >= max_z - min_z) {
    split_direction = Y_axis;
  } else {
    split_direction = Z_axis;
  }

  /* Sort particles along splitting direction and apply median splitting */
  qsort_r(pid, count, sizeof(int), &cmp, coords[split_direction]);
  int median_idx = count / 2 + 1;

  /* Populate the left and right subtree of this bvh */
  bvh->left = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->left, parts, coords, pid, median_idx);
  bvh->right = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->right, parts, coords, &pid[median_idx],
                   count - median_idx);

  /* Set the other fields of this bvh */
  bvh->bbox = bbox_wrap(&bvh->left->bbox, &bvh->right->bbox);
  bvh->data = NULL;
  bvh->count = 0;
}

