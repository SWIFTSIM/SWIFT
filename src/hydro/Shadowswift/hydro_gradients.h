//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H

#ifdef SHADOWSWIFT_GRADIENTS

#ifndef SHADOWSWIFT_GRADIENTS_WLS

/* Green-Gauss gradient estimates (default) */
#include "hydro_gradients_shadowswift.h"

#else

/* Weighted least square gradient estimates */
#include "hydro_gradients_wls.h"

#endif

#else

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/**
 * @brief Update the gradient estimation for a particle using a given
 * neighbouring particle.
 *
 * @param pi Particle we are updating.
 * @param pj Particle we are using to update pi.
 */
__attribute__((always_inline)) INLINE void hydro_gradients_collect(
    struct part* restrict pi, const struct part* restrict pj,
    const double* restrict centroid, const double* restrict dx, double r,
    float surface_area) {}

/**
 * @brief Extrapolate all the quantities over the given distance.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_extrapolate(
    const struct part* p, const float* dx, float* dW) {
  dW[0] = 0.f;
  dW[1] = 0.f;
  dW[2] = 0.f;
  dW[3] = 0.f;
  dW[4] = 0.f;
}

/**
 * @brief Gradients reconstruction. Empty for no gradients.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, const float* dx,
    float r, const float* xij_i, float dt, float* Wi, float* Wj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {}
#endif

/**
 * @brief Initialize gradient variables.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part* p) {
  p->gradients.rho[0] = 0.0f;
  p->gradients.rho[1] = 0.0f;
  p->gradients.rho[2] = 0.0f;

  p->gradients.v[0][0] = 0.0f;
  p->gradients.v[0][1] = 0.0f;
  p->gradients.v[0][2] = 0.0f;

  p->gradients.v[1][0] = 0.0f;
  p->gradients.v[1][1] = 0.0f;
  p->gradients.v[1][2] = 0.0f;

  p->gradients.v[2][0] = 0.0f;
  p->gradients.v[2][1] = 0.0f;
  p->gradients.v[2][2] = 0.0f;

  p->gradients.P[0] = 0.0f;
  p->gradients.P[1] = 0.0f;
  p->gradients.P[2] = 0.0f;

#ifdef SHADOWSWIFT_GRADIENTS_WLS
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->gradients.matrix_wls[i][j] = 0.f;
    }
  }
#endif
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
