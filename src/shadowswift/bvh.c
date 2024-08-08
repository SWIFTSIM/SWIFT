#include "bvh.h"

#include "qselect.h"

#include <stdlib.h>

/**
 * @file Function implementations from bvh.h.
 **/

inline static int bvh_split_midpoint(double midpoint,
                                     enum direction split_direction,
                                     const struct part *parts, int *pid,
                                     int count);
inline static int cmp(const void *a, const void *b, void *arg);

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
int flat_bvh_hit_rec(const bvh_t *bvh, int node_id, struct part *parts,
                     double x, double y, double z, double r2) {
  flat_bvh_node_t *node = &bvh->nodes[node_id];

  /* Anything to do here? */
  if (!bbox_contains(&node->bbox, x, y, z)) return -1;

  /* Are we in a leaf? */
  if (node->is_leaf) {

    /* Check individual particles */
    for (int i = 0; i < node->data[BVH_DATA_SIZE]; i++) {
      int pid = node->data[i];
      const struct part *p = &parts[pid];
      double hit_r2 = p->geometry.search_radius * p->geometry.search_radius;
      double dx[3] = {p->x[0] - x, p->x[1] - y, p->x[2] - z};
      double dx2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (dx2 <= r2 && dx2 < hit_r2) return pid;
    }

    /* no hit */
    return -1;
  }

  /* Else, recurse */
#ifdef SHADOWSWIFT_BVH_INSERT_BFO
  /* Check if hit with central part of this node */
  int pid = node->data[BVH_DATA_SIZE];
  const struct part *p = &parts[pid];
  double hit_r2 = p->geometry.search_radius * p->geometry.search_radius;
  double dx[3] = {p->x[0] - x, p->x[1] - y, p->x[2] - z};
  double dx2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  if (dx2 <= r2 && dx2 < hit_r2) return pid;
    /* Else check subtrees */
#endif
  int hit = flat_bvh_hit_rec(bvh, node->children.left, parts, x, y, z, r2);
  if (hit == -1) {
    hit = flat_bvh_hit_rec(bvh, node->children.right, parts, x, y, z, r2);
  }
  return hit;
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
void flat_bvh_populate_rec(bvh_t *bvh, int node_id, const struct part *parts,
                           int *pid, int count) {

#ifdef SWIFT_DEBUG_CHECKS
  if (count < 0) {
    error("Trying to populate BVH with negative particle count!");
  }
#endif

  flat_bvh_node_t *node = &bvh->nodes[node_id];

  /* Construct leaf node? */
  if (count <= BVH_DATA_SIZE) {
    /* Copy data */
    node->data[BVH_DATA_SIZE] = count;
    for (int i = 0; i < count; i++) node->data[i] = pid[i];

    /* Set remaining fields for leaf */
    node->bbox = bbox_from_parts(parts, pid, count);
    node->is_leaf = 1;
    return;
  }

  /* Else, build bvh recursively */

  /* Determine split_direction (direction with the largest spread) */
  double anchor[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  double opposite[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    for (int k = 0; k < 3; k++) {
      anchor[k] = fmin(anchor[k], p->x[k]);
      opposite[k] = fmax(opposite[k], p->x[k]);
    }
  }
  double width[3] = {opposite[0] - anchor[0], opposite[1] - anchor[1],
                     opposite[2] - anchor[2]};
  enum direction split_direction;
  if (width[0] >= width[1] && width[0] >= width[2]) {
    split_direction = X_axis;
  } else if (width[1] >= width[0] && width[1] >= width[2]) {
    split_direction = Y_axis;
  } else {
    split_direction = Z_axis;
  }

#if BVH_SPLITTING_METHOD == BVH_MIDPOINT_SPLIT
  /* Now flip particles around, until they are divided into  two groups along
   * the split direction (those below and above the midpoint) */
  double midpoint = 0.5 * (anchor[split_direction] + opposite[split_direction]);
  int pivot = bvh_split_midpoint(midpoint, split_direction, parts, pid, count);
#ifdef SHADOWSWIFT_BVH_INSERT_BFO
  /* Finally find central particle of this node. Take it from the largest
   * subtree */
  int pid_central;
  if (pivot > (count - pivot)) {
    /* Search for the right-most particle in the left subtree */
    int idx_max = 0;
    float v_max = parts[pid[0]].x[split_direction];
    for (int i = 1; i < pivot; i++) {
      float v = parts[pid[i]].x[split_direction];
      if (v > v_max) {
        idx_max = i;
        v_max = v;
      }
    }
    /* Save it and flip the end of this subtrees array in its place */
    pid_central = pid[idx_max];
    pid[idx_max] = pid[pivot - 1];
    pid[pivot - 1] = pid_central;
  } else {
    /* Search for the left-most particle in the right subtree */
    int idx_min = pivot;
    float v_min = parts[pid[pivot]].x[split_direction];
    for (int i = pivot + 1; i < count; i++) {
      float v = parts[pid[i]].x[split_direction];
      if (v < v_min) {
        idx_min = i;
        v_min = v;
      }
    }
    /* Save it and flip the end of this subtrees array in its place */
    pid_central = pid[idx_min];
    pid[idx_min] = pid[pivot];
    pid[pivot] = pid_central;
    /* Also increase pivot to indicate that the start of the right subtrees
     * array has shifted */
    pivot++;
  }
#endif
#elif BVH_SPLITTING_METHOD == BVH_MEDIAN_SPLIT
  /* The pivot will be the index of the median. */
  int pivot = (count - 1) / 2 + (count - 1) % 2;
  struct {
    enum direction split_direction;
    const struct part *parts;
  } data = {
      .split_direction = split_direction,
      .parts = parts,
  };
  qselect_r(pivot, pid, count, sizeof(int), &cmp, &data);
#ifdef SHADOWSWIFT_BVH_INSERT_BFO
  int pid_central = pid[pivot];
  pivot++;
#endif
#else
#error "Unkown or undefined BVH_SPLITTING_METHOD!"
#endif

  if (pivot >= count)
    error("Failed to partition elements in bvh construction!");

  /* Populate the left and right subtrees */
  int left_id = flat_bvh_new_node(bvh);
#ifdef SHADOWSWIFT_BVH_INSERT_BFO
  flat_bvh_populate_rec(bvh, left_id, parts, pid, pivot - 1);
#else
  flat_bvh_populate_rec(bvh, left_id, parts, pid, pivot);
#endif
  int right_id = flat_bvh_new_node(bvh);
  flat_bvh_populate_rec(bvh, right_id, parts, &pid[pivot], count - pivot);

  /* Nodes might have been moved */
  node = &bvh->nodes[node_id];
  /* Update fields of this node */
  node->children.left = left_id;
  node->children.right = right_id;
  node->is_leaf = 0;
  bbox_wrap(&bvh->nodes[left_id].bbox, &bvh->nodes[right_id].bbox, &node->bbox);
#ifdef SHADOWSWIFT_BVH_INSERT_BFO
  /* Store the central part separately and make sure it is contained in the
   * bbox */
  node->data[BVH_DATA_SIZE] = pid_central;
  const struct part *p_central = &parts[pid_central];
  bbox_t bbox_part;
  for (int i = 0; i < 3; i++) {
    bbox_part.anchor[i] = p_central->x[i] - p_central->geometry.search_radius;
    bbox_part.opposite[i] = p_central->x[i] + p_central->geometry.search_radius;
  }
  bbox_wrap(&node->bbox, &bbox_part, &node->bbox);
#endif
#ifdef SWIFT_DEBUG_CHECKS
  if (node->is_leaf) {
    assert(node->data[BVH_DATA_SIZE] <= BVH_DATA_SIZE);
  } else {
    assert(node->children.left < bvh->count && 0 <= node->children.left);
    assert(node->children.right < bvh->count && 0 <= node->children.right);
  }
#endif
}

/** @brief Reorder the pid array so that particles to the left of #midpoint come
 * first and particles to the right of #midpoint afterwards.
 * @param midpoint The chosen midpoint
 * @param split_direction The direction along which the particles must be
 * ordered.
 * @param parts The particle array
 * @param pid An array containing the indices of the particles that must be
 * reordered.
 * @param count The length of the #pid array.
 * @return The index of the first part-id (in pid) to the right of #midpoint.
 * This is also equal to the number of parts to the left of #midpoint.
 */
inline static int bvh_split_midpoint(double midpoint,
                                     enum direction split_direction,
                                     const struct part *parts, int *pid,
                                     int count) {
  int head = 0, tail = count - 1;
  while (head < tail) {
    /* Find the first part that belongs in the second group */
    while (parts[pid[head]].x[split_direction] <= midpoint && head < tail) {
      head++;
    }
    /* Find the last part that belongs in the first group */
    while (parts[pid[tail]].x[split_direction] > midpoint && head < tail) {
      tail--;
    }
    if (head < tail) {
      /* Swap the elements at head and tail */
      int tmp = pid[head];
      pid[head] = pid[tail];
      pid[tail] = tmp;
      head++;
      tail--;
    }
  }
  return head;
}

/**
 * @brief Computes the ordering of 2 indices using some associated data (in this
 * case: particle coordinates).
 *
 * @param a Pointer to the index of the first element (will be cast to (int *))
 * @param b Pointer to the index of the second element (will be cast to (int *))
 * @param arg Pointer to structure containing the split direction and a const
 * pointer to the particles array (used to access the particle coordinates).
 * @returns -1 if arg[*a] < arg[*b] 0 if they are equal and 1 otherwise.
 **/
inline static int cmp(const void *a, const void *b, void *arg) {
  int ai = *(int *)a;
  int bi = *(int *)b;
  typedef struct {
    enum direction split_direction;
    const struct part *parts;
  } arg_t;
  arg_t *data = (arg_t *)arg;
  double ad = data->parts[ai].x[data->split_direction];
  double bd = data->parts[bi].x[data->split_direction];

  if (ad < bd) {
    return -1;
  } else if (ad > bd) {
    return 1;
  } else {
    return 0;
  }
}
