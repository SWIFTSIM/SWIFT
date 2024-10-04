//
// Created by yuyttenh on 21/03/22.
//

#ifndef SWIFTSIM_CELL_GRID_H
#define SWIFTSIM_CELL_GRID_H

#include "const.h"
#include "shadowswift/voronoi.h"

/*! @brief Enum indicating the completeness for the Voronoi mesh of this cell.
 *
 * A cell is considered complete when it and its neighbours on the same level in
 * the AMR have at least one particle in every 1/27th cube of the cell (obtained
 * by dividing cells in three along all axes).
 *
 * The Voronoi grid can safely be constructed on any level where the cell is
 * complete. */
enum grid_completeness {
  grid_invalidated_completeness = 0,
  grid_complete,
  grid_incomplete,
};

struct cell_grid {
  /*! Pointer to parent where the grid is constructed. */
  struct cell *construction_level;

  /*! Pointer to shallowest parent of this cell used in any pair construction
   * task. Can be above the construction level of this cell. We need to drift at
   * this level. */
  struct cell *super;

  /*! Whether this cell is complete (at least one particle in all 27 sub-cells
   * if this cell is divided in thirds along each axis). */
  enum grid_completeness self_completeness;

  /*! Whether this cell is itself complete and has directly neighbouring cell
   * on the same level in all directions which are also complete. */
  int complete;

#ifdef WITH_MPI
  /*! Flags indicating whether we should send the faces for the corresponding
   * SIDs over MPI */
  int send_flags;
#endif

  /*! Pointer to the voronoi struct of this cell (if any) */
  struct voronoi *voronoi;

  /*! Pointer to this cells construction task. */
  struct task *construction;

  /*! Linked list of this cells outgoing construction synchronization tasks
   * (i.e. for cells that need this cell for their construction task) */
  struct link *sync_out;

  /*! Linked list of this cells incoming construction synchronization tasks
   * (i.e. cells needed for this cell's construction task) */
  struct link *sync_in;

  /*! Time of last construction */
  integertime_t ti_old;
};

struct pcell_faces {
  size_t counts[27];

  struct voronoi_pair faces[];
};

/*! @brief Enum used to indicate whether a cell is above, below or on the
 * construction level. Only used in the packed cell representation */
enum grid_construction_level {
  above_construction_level,
  on_construction_level,
  below_construction_level
};

#endif  // SWIFTSIM_CELL_GRID_H
