//
// Created by yuyttenh on 21/03/22.
//

#ifndef SWIFTSIM_CELL_GRID_H
#define SWIFTSIM_CELL_GRID_H

#include "const.h"
#include "shadowswift/bvh.h"
#include "shadowswift/delaunay.h"
#include "shadowswift/voronoi.h"

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

  /*! Whether this cell is complete (at least one particle in all 27 sub-cells
   * if this cell is divided in thirds along each axis). */
  int self_complete;

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

  /*! Pointer to the delaunay struct of this cell (if any) */
  struct delaunay *delaunay;

#ifdef SHADOWSWIFT_BVH
  /*! Pointer to the bvh struct of this cell (if any) */
  struct BVH *bvh;
#endif

  /*! Indices sorting the particles of this cell according to their hilbert
   * ordering */
  int *hilbert_r_sort;

  /*! Linked list of this cells construction tasks. */
  struct link *construction;

#ifdef SHADOWSWIFT_BVH
  /*! Pointer to this cells BVH construction task. */
  struct task *build_bvh;
#endif

  /*! Pointer to this cells construction ghost task. */
  struct task *ghost;

  /*! Time of last construction */
  integertime_t ti_old;

  /*! @brief Array of booleans indicating whether or not neighbouring particles
   * have been tried to be added for a given sid. If this is 0 for a given sid,
   * this means that this cell should get the reflective boundary condition
   * applied for that sid. */
  unsigned long int sid_is_inside_face_mask;
};

/*! @brief The default sid mask for marking faces as inside */
#ifdef HYDRO_DIMENSION_2D
/*                         111_101_111_101_111_101_111_101_111 */
#define DEFAULT_SID_MASK 0b111101111101111101111101111ul;
#elif defined(HYDRO_DIMENSION_3D)
/*                         111_101_111_101_010_101_111_101_111 */
#define DEFAULT_SID_MASK 0b111101111101010101111101111ul;
#endif

struct pcell_faces {
  size_t counts[27];

  struct voronoi_pair faces[];
};

void cell_grid_construct_local_voronoi(struct cell *restrict c,
                                       const struct engine *e, const int *pid,
                                       const int count, struct voronoi *vortess,
                                       struct part *parts_out);

void cell_grid_doself_apoptosis(struct cell *c, struct engine *e, int p_idx);

void cell_grid_dopair_apoptosis(struct cell *restrict ci,
                                struct cell *restrict cj, struct engine *e,
                                int p_idx);

#endif  // SWIFTSIM_CELL_GRID_H
