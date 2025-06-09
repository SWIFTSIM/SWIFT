#ifndef SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H
#define SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H

#include "const.h"

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void fvpm_reset_centroids(
    struct part* restrict p) {

  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param p Particle.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void fvpm_normalise_centroid(
    struct part* restrict p, const float wcount) {

  const float norm = kernel_norm / wcount;
  p->geometry.centroid[0] *= norm;
  p->geometry.centroid[1] *= norm;
  p->geometry.centroid[2] *= norm;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param p Particle (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_left(
    struct part* restrict p, const float* dx, const float w) {

  p->geometry.centroid[0] -= dx[0] * w;
  p->geometry.centroid[1] -= dx[1] * w;
  p->geometry.centroid[2] -= dx[2] * w;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param p Particle (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_right(
    struct part* restrict p, const float* dx, const float w) {

  p->geometry.centroid[0] += dx[0] * w;
  p->geometry.centroid[1] += dx[1] * w;
  p->geometry.centroid[2] += dx[2] * w;
}

#endif /* SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H */
