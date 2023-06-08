//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H

#include <float.h>

#ifdef SHADOWSWIFT_SLOPE_LIMITER_PER_FACE

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION \
  "GIZMO piecewise slope limiter (Hopkins 2015)"
#include "hydro_slope_limiters_face.h"

#else

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION "No piecewise slope limiter"

/**
 * @brief Empty for no piece wise slope limiters.
 *
 * @param Wi Hydrodynamic variables of particle i.
 * @param Wj Hydrodynamic variables of particle j.
 * @param dWi Difference between the hydrodynamic variables of particle i at the
 * position of particle i and at the interface position.
 * @param dWj Difference between the hydrodynamic variables of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, const float *xij_i,
    const float *xij_j, float r) {}

#endif

#ifdef SHADOWSWIFT_SLOPE_LIMITER_CELL_WIDE

#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION \
  "Cell wide slope limiter (Springel 2010)"

#include "hydro_slope_limiters_cell.h"

#elif defined(SHADOWSWIFT_SLOPE_LIMITER_MESHLESS)

#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION \
  "Meshless slope limiter (Gizmo 2015)"

#include "hydro_slope_limiters_meshless.h"

#else

#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION "No cell wide slope limiter"

/**
 * @brief Empty for no cell wide slope limiters.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param dx vector pointing from pi to the centroid of the face between pi and
 * pj.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj,
                               const float *dx) {}

/**
 * @brief Empty for no cell wide slope limiters.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part *p) {}

#endif

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limiter_prepare(
    struct part *p) {
  p->limiter.rho[0] = p->rho;
  p->limiter.rho[1] = p->rho;
  p->limiter.v[0][0] = p->v[0];
  p->limiter.v[0][1] = p->v[0];
  p->limiter.v[1][0] = p->v[1];
  p->limiter.v[1][1] = p->v[1];
  p->limiter.v[2][0] = p->v[2];
  p->limiter.v[2][1] = p->v[2];
  p->limiter.P[0] = p->P;
  p->limiter.P[1] = p->P;
  p->limiter.A[0] = p->A;
  p->limiter.A[1] = p->A;

#ifdef SHADOWSWIFT_SLOPE_LIMITER_MESHLESS
  p->limiter.r_max = 0.f;
#endif

#ifdef SHADOWSWIFT_SLOPE_LIMITER_CELL_WIDE
  p->limiter.extrapolations.rho[0] = 0.f;
  p->limiter.extrapolations.rho[1] = 0.f;
  p->limiter.extrapolations.v[0][0] = 0.f;
  p->limiter.extrapolations.v[0][1] = 0.f;
  p->limiter.extrapolations.v[1][0] = 0.f;
  p->limiter.extrapolations.v[1][1] = 0.f;
  p->limiter.extrapolations.v[2][0] = 0.f;
  p->limiter.extrapolations.v[2][1] = 0.f;
  p->limiter.extrapolations.P[0] = 0.f;
  p->limiter.extrapolations.P[1] = 0.f;
  p->limiter.extrapolations.A[0] = 0.f;
  p->limiter.extrapolations.A[1] = 0.f;
#endif
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
