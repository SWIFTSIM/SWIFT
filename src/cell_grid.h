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
  struct cell *construction_level;

  /*! Pointer to shallowest parent of this cell used in any pair construction
   * task. Can be above the construction level of this cell. We need to drift at
   * this level. */
  struct cell *super;

  /*! Flag that indicates whether this cell is unsplittable or has a directly
   * neighbouring cell on the same level that is unsplittable. */
  int unsplittable_flag;

  /*! Whether this cell can safely be split for grid construction (cell can be
   * split if all of it's progeny is complete). */
  /* TODO send this over MPI (add it in cell_pack and cell_unpack) */
  int split;

  /*! Flag indicating whether this cell satisfies the completeness criterion.
   * The voronoi grid of a cell is guaranteed to be completed by only particles
   * of its directly neighbouring cells if when we would split that cell in
   * thirds along each dimension (i.e. in 27 smaller cells), every small cube
   * would contain at least one particle. */
   int complete;

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
