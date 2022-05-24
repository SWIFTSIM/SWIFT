//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_VORONOI_H
#define SWIFTSIM_SHADOWSWIFT_VORONOI_H

/*! @brief Store the edges of faces (so that the actual Voronoi grid can be
 *  reconstructed). */
#ifdef SHADOWSWIFT_OUTPUT_GRIDS
#define VORONOI_STORE_FACES
#endif

/*! @brief Store cell generators. */
//#define VORONOI_STORE_GENERATORS

/*! @brief Activate runtime assertions. */
//#define VORONOI_DO_ASSERTIONS

/*! @brief Activate extra checks */
//#define VORONOI_CHECKS

/*! @brief The minimal relative face size in 1D of voronoi faces */
#define MIN_REL_FACE_SIZE 1e-7


/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#if defined(VORONOI_DO_ASSERTIONS) || defined(VORONOI_CHECKS)
#define voronoi_assert(condition)                                     \
  if (!(condition)) {                                                 \
    fprintf(stderr, "%s:%s():%i: Condition failed: " #condition "\n", \
            __FILE__, __FUNCTION__, __LINE__);                        \
    abort();                                                          \
  }
#else
#define voronoi_assert(condition)
#endif

/**
 * @brief Voronoi cell.
 *
 * A cell stores geometrical information about a Voronoi cell: its volume and
 * the location of its centroid (for compatibility reasons, these must always be
 * 3D vectors (also for the 2D voronoi algorithm).
 */
struct voronoi_cell {
  /*! Cell volume. */
  double volume;

  /*! Cell centroid. */
  double centroid[3];

#ifdef VORONOI_STORE_GENERATORS
  /*! Position of the cell generator. */
  double x[3];
#endif

  /*! Number of faces of this cell. */
  int nface;

  /*! cell_pair_connections offset */
  int pair_connections_offset;
};

#ifdef HYDRO_DIMENSION_2D
#include "algorithm2d/voronoi.h"
#elif defined(HYDRO_DIMENSION_3D)
#include "algorithm3d/voronoi.h"
#else
#error "Only 2D and 3D schemes are supported by ShadowSWIFT!"
#endif

#endif  // SWIFTSIM_SHADOWSWIFT_VORONOI_H
