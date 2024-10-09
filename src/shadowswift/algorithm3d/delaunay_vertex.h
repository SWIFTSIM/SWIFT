//
// Created by yuyttenh on 5/7/24.
//

#ifndef SWIFTSIM_DELAUNAY_VERTEX_H
#define SWIFTSIM_DELAUNAY_VERTEX_H

typedef struct delaunay_vertex {
  union {
    double x_f64[3];
    uint64_t x_u64[3];
  };
} delaunay_vertex_t;

#endif  // SWIFTSIM_DELAUNAY_VERTEX_H
