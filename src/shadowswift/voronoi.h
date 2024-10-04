//
// Created by yuyttenh on 10/4/24.
//

#ifndef SWIFTSIM_SHADOWSWIFT_VORONOI_H
#define SWIFTSIM_SHADOWSWIFT_VORONOI_H

struct voronoi {
  int pair_count[27];
};

struct voronoi_pair {};

__attribute__((always_inline)) INLINE static void voronoi_destroy(struct voronoi* v) {}

#endif  // SWIFTSIM_SHADOWSWIFT_VORONOI_H
