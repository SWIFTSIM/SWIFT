/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *                             Yolan Uyttenhove (Yolan.Uyttenhove@UGent.be)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFTSIM_SHADOWSWIFT_VORONOI_H
#define SWIFTSIM_SHADOWSWIFT_VORONOI_H

#include "error.h"

/*! @brief Store the edges of faces (so that the actual Voronoi grid can be
 *  reconstructed). */
#ifdef SHADOWSWIFT_OUTPUT_GRIDS
#define VORONOI_STORE_FACES
#endif

/*! @brief Whether to store cell-face connections. */
// #define VORONOI_STORE_CELL_FACE_CONNECTIONS

#ifdef SWIFT_DEBUG_CHECKS
/*! @brief Activate runtime assertions. */
#define VORONOI_DO_ASSERTIONS

/*! @brief Activate extra (expensive) checks */
// #define VORONOI_CHECKS
#endif

/*! @brief The minimal relative face size in 1D of voronoi faces */
#define MIN_REL_FACE_SIZE 0.0

/*! @brief The maximum number of vertices on a voronoi face (3D) */
#define VORONOI_CONSTRUCTION_MAX_FACES 1000

/*! @brief The maximum number of vertices on a voronoi face (3D) */
#define VORONOI_CONSTRUCTION_MAX_FACE_VERTICES 1000

/**
 *@brief Evaluate the given condition and abort if it evaluates to true.
 *
 * This macro is similar to the standard assert() macro.
 * This macro is only defined when VORONOI_DO_ASSERTIONS is active.
 */
#if defined(VORONOI_DO_ASSERTIONS) || defined(VORONOI_CHECKS)
#define voronoi_assert(condition)                                     \
  if (!(condition)) {                                                 \
    error("%s:%s():%i: Condition failed: " #condition "\n", __FILE__, \
          __FUNCTION__, __LINE__);                                    \
  }
#else
#define voronoi_assert(condition)
#endif

#ifdef MOVING_MESH
#ifdef HYDRO_DIMENSION_1D
#include "algorithm1d/voronoi.h"
#elif defined(HYDRO_DIMENSION_2D)
#include "algorithm2d/voronoi.h"
#elif defined(HYDRO_DIMENSION_3D)
#include "algorithm3d/voronoi.h"
#else
#error "Unknown hydro dimensionality!"
#endif
#else
struct voronoi_pair {
  double midpoint[3];
  double surface_area;
  int left_idx;
  int right_idx;
  int sid;
};

struct voronoi {
  struct voronoi_pair* pairs[28];
  int pair_count[28];
};

static inline void voronoi_destroy(struct voronoi* restrict v) {}
#endif

inline static struct voronoi_pair* voronoi_get_local_faces(struct voronoi* v,
                                                           int* count) {
  *count = v->pair_count[13];
  return v->pairs[13];
}

inline static struct voronoi_pair* voronoi_get_sid_faces(struct voronoi* v,
                                                         int sid, int* count) {
  *count = v->pair_count[sid];
  return v->pairs[sid];
}

inline static struct voronoi_pair* voronoi_get_boundary_faces(struct voronoi* v,
                                                              int* count) {
  *count = v->pair_count[27];
  return v->pairs[27];
}

#endif  // SWIFTSIM_SHADOWSWIFT_VORONOI_H
