#include "bvh.h"

void bvh_destroy(struct BVH *bvh) {

  if (bvh->left != NULL) {
    bvh_destroy(bvh->left);
  }
  if (bvh->right != NULL) {
    bvh_destroy(bvh->right);
  }

  if (bvh->data.pid != NULL) {
    free(bvh->data.pid);
#ifdef SWIFT_DEBUG_CHECKS
    if (bvh->data.radius == NULL)
      error("Inconsistent data for BVH leaf!");
#endif
    free(bvh->data.radius);
  }

  free(bvh);
}

void bvh_populate_rec(struct BVH* bvh, const double* coords[3],
                      const double* search_radii, const int* pid, int* idx,
                      int count) {
  if (count <= 4) {
    /* Set unused fields of this bvh */
    bvh->left = NULL;
    bvh->right = NULL;

    /* Set used fields for leaf */
    bvh->data.pid = malloc(count * sizeof(*bvh->data.pid));
    bvh->data.radius = malloc(count * sizeof(*bvh->data.radius));
    for (int i = 0; i < count; i++) {
      bvh->data.pid[i] = pid[idx[i]];
      bvh->data.radius[i] = search_radii[idx[i]];
    }
    bvh->count = count;
    bvh->bbox = bbox_from_coords(coords, search_radii, idx, count);
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
    const double x[3] = {coords[0][idx[i]], coords[1][idx[i]],
                         coords[2][idx[i]]};
    min_x = min(min_x, x[0]);
    max_x = max(max_x, x[0]);
    min_y = min(min_y, x[1]);
    max_y = max(max_y, x[1]);
    min_z = min(min_z, x[2]);
    max_z = max(max_z, x[2]);
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
  qsort_r(idx, count, sizeof(int), &cmp, (void*)coords[split_direction]);
  int median_idx = (count + 1) / 2;

  /* Populate the left and right subtree of this bvh */
  bvh->left = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->left, coords, search_radii, pid, idx, median_idx);
  bvh->right = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->right, coords, search_radii, pid, &idx[median_idx],
                   count - median_idx);

  /* Set the other fields of this bvh */
  bvh->bbox = bbox_wrap(&bvh->left->bbox, &bvh->right->bbox);
  bvh->data.pid = NULL;
  bvh->data.radius = NULL;
  bvh->count = count;
}
