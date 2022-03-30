//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION "No piecewise slope limiter"
#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION "No cell wide slope limiter"

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part *p) {}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj, float r) {}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part *p) {}

/**
 * @brief Slope limit face
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, const float *xij_i,
    const float *xij_j, float r) {}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
