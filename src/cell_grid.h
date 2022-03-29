//
// Created by yuyttenh on 21/03/22.
//

#ifndef SWIFTSIM_CELL_GRID_H
#define SWIFTSIM_CELL_GRID_H

#include "shadowswift/voronoi.h"
#include "shadowswift/delaunay.h"

struct cell_grid {
  /*! Pointer to the parent cell of this cell containing the Voronoi grid (if
   * any). */
  struct cell *super;

  /*! Flag that indicates whether this cell is unsplittable or has a directly
   * neighbouring cell on the same level that is unsplittable. */
  int unsplittable_flag;

  /*! Whether this cell can safely be split for grid construction (a cell can be
   * split if all 27 small cubes of dimensions 1/3 the dimensions of this cell
   * making up the cell, contain at least 1 hydro particle). */
  /* TODO send this over MPI (add it in cell_pack and cell_unpack) */
  int split;

  /*! The maximal search radius of any particle in the voronoi tessellation of
   * this cell. */
  double r_max;

  struct voronoi *voronoi;

  struct delaunay *delaunay;

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
