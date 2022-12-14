//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H

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

#ifdef SHADOWSWIFT_GRADIENTS

__attribute__((always_inline)) INLINE static void
hydro_gradients_apply_extrapolation(float* W, const float* restrict dW) {
  if (-dW[0] > W[0]) {
    message(
        "Gradient extrapolation would lead to unphysical Density!"
        "Skipping extrapolation!");
  } else if (-dW[4] > W[4]) {
    message(
        "Gradient extrapolation would lead to unphysical Pressure!"
        "Skipping extrapolation!");
  } else if (-dW[5] > W[5]) {
    message(
        "Gradient extrapolation would lead to unphysical Entropy!"
        "Skipping extrapolation!");
  } else {
    W[0] += dW[0];
    W[1] += dW[1];
    W[2] += dW[2];
    W[3] += dW[3];
    W[4] += dW[4];
    W[5] += dW[5];
  }
}

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


#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
