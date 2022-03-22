//
// Created by yuyttenh on 21/03/22.
//

#ifndef SWIFTSIM_CELL_GRID_H
#define SWIFTSIM_CELL_GRID_H

#include "voronoi/voronoi.h"

struct cell_grid {
  struct cell *super;

  int active;

  int unsplittable_flag;

  struct voronoi *voronoi;

  struct delaunay *delaunay;

  struct link *construction;

  struct task *ghost;
};

struct pcell_grid {
  int face_counts[26];

  struct voronoi_pair faces[];
};

#endif  // SWIFTSIM_CELL_GRID_H
