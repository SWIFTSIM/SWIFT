//
// Created by yuyttenh on 21/03/22.
//

#ifndef SWIFTSIM_CELL_GRID_H
#define SWIFTSIM_CELL_GRID_H

#include "shadowswift/voronoi.h"
#include "shadowswift/delaunay.h"
#include "const.h"

enum construction_level {
  above_construction_level,
  on_construction_level,
  below_construction_level,
  uninitialized,
};

struct cell_grid {
  /*! Flag indicating the construction level location. */
  enum construction_level construction_level;

  /*! Pointer to shallowest parent of this cell used in any pair construction
   * task. Can be above the construction level of this cell. We need to drift at
   * this level. */
  struct cell *super;

  /*! Flag that indicates whether this cell is unsplittable or has a directly
   * neighbouring cell on the same level that is unsplittable. */
  int unsplittable_flag;

  /*! Whether this cell can safely be split for grid construction (cell can be
   * split if all of it's progeny is complete). */
  int split;

  /*! Pointer to the voronoi struct of this cell (if any) */
  struct voronoi *voronoi;

  /*! Pointer to the delaunay struct of this cell (if any) */
  struct delaunay *delaunay;

  /*! Indices sorting the particles of this cell according to their hilbert
   * ordering */
  int *hilbert_r_sort;

  /*! Linked list of this cells construction tasks. */
  struct link *construction;

  /*! Pointer to this cells construction ghost task. */
  struct task *ghost;

  /*! Time of last construction */
  integertime_t ti_old;
};

struct pcell_grid {
  int face_counts[26];

  struct voronoi_pair faces[];
};

#endif  // SWIFTSIM_CELL_GRID_H
