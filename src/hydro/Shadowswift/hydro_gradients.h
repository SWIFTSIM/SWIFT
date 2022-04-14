//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/* Forward declarations */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, const float *xij_i,
    const float *xij_j, float r);

/**
 * @brief Initialize gradient variables
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
}

/**
 * @brief Add the gradient estimate for a single quantity due to a particle pair
 * to the total gradient for that quantity
 *
 * This corresponds to one term of equation (21) in Springel (2010).
 *
 * @param qL Value of the quantity on the left.
 * @param qR Value of the quantity on the right.
 * @param cLR Vector pointing from the midpoint of the particle pair to the
 * geometrical centroid of the face in between the particles.
 * @param xLR Vector pointing from the right particle to the left particle.
 * @param rLR Distance between two particles.
 * @param A Surface area of the face in between the particles.
 * @param grad Current value of the gradient for the quantity (is updated).
 */
__attribute__((always_inline)) INLINE void hydro_gradients_single_quantity(
    float qL, float qR, const double* cLR, const double* xLR, double rLR,
    float A, float* grad) {

  grad[0] += A * ((qR - qL) * cLR[0] / rLR - 0.5f * (qL + qR) * xLR[0] / rLR);
  grad[1] += A * ((qR - qL) * cLR[1] / rLR - 0.5f * (qL + qR) * xLR[1] / rLR);
  grad[2] += A * ((qR - qL) * cLR[2] / rLR - 0.5f * (qL + qR) * xLR[2] / rLR);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {
  float volume = p->geometry.volume;

  p->gradients.rho[0] /= volume;
  p->gradients.rho[1] /= volume;
  p->gradients.rho[2] /= volume;

  p->gradients.v[0][0] /= volume;
  p->gradients.v[0][1] /= volume;
  p->gradients.v[0][2] /= volume;
  p->gradients.v[1][0] /= volume;
  p->gradients.v[1][1] /= volume;
  p->gradients.v[1][2] /= volume;
  p->gradients.v[2][0] /= volume;
  p->gradients.v[2][1] /= volume;
  p->gradients.v[2][2] /= volume;

  p->gradients.P[0] /= volume;
  p->gradients.P[1] /= volume;
  p->gradients.P[2] /= volume;
}

/**
 * @brief Gradients time extrapolation (makes scheme second order in time).
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_extrapolate_in_time(const struct part* p, const float* W,
                                    float dt, float* dW) {
  const float div_v =
      p->gradients.v[0][0] + p->gradients.v[1][1] + p->gradients.v[2][2];

  dW[0] = -0.5f * dt *
           (W[0] * div_v + W[1] * p->gradients.rho[0] +
            W[2] * p->gradients.rho[1] + W[3] * p->gradients.rho[2]);

  if (W[0] != 0.0f) {
    const float rho_inv = 1.f / W[0];
    dW[1] = -0.5f * dt * (W[1] * div_v + rho_inv * p->gradients.P[0]);
    dW[2] = -0.5f * dt * (W[2] * div_v + rho_inv * p->gradients.P[1]);
    dW[3] = -0.5f * dt * (W[3] * div_v + rho_inv * p->gradients.P[2]);
  } else {
    dW[1] = 0.0f;
    dW[2] = 0.0f;
    dW[3] = 0.0f;
  }
  dW[4] = -0.5f * dt *
           (hydro_gamma * W[4] * div_v + W[1] * p->gradients.P[0] +
            W[2] * p->gradients.P[1] + W[3] * p->gradients.P[2]);

  /* Sanity check */
  if (W[0] + dW[0] < 0) {
    dW[0] = 0.f;
  }
  if (W[4] + dW[4] < 0) {
    dW[4] = 0.f;
  }
}

/**
 * @brief Extrapolate the given gradient over the given distance.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float
hydro_gradients_extrapolate_single_quantity(const float* gradient,
                                            const float* dx) {

  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Extrapolate all the quantities over the given distance.
 *
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_extrapolate(
    const struct part* p, const float* dx, float* dW) {

  float drho[3], dvx[3], dvy[3], dvz[3], dP[3];
  hydro_part_get_gradients(p, drho, dvx, dvy, dvz, dP);

  dW[0] += hydro_gradients_extrapolate_single_quantity(drho, dx);
  dW[1] += hydro_gradients_extrapolate_single_quantity(dvx, dx);
  dW[2] += hydro_gradients_extrapolate_single_quantity(dvy, dx);
  dW[3] += hydro_gradients_extrapolate_single_quantity(dvz, dx);
  dW[4] += hydro_gradients_extrapolate_single_quantity(dP, dx);
}

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, const float* dx,
    float r, const float* xij_i, float dt, float* Wi, float* Wj) {

  /* perform gradient reconstruction in space and time */
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8) */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];
  float drho_j[3], dvx_j[3], dvy_j[3], dvz_j[3], dP_j[3];
  hydro_part_get_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);
  hydro_part_get_gradients(pj, drho_j, dvx_j, dvy_j, dvz_j, dP_j);

  float dWi[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
  hydro_gradients_extrapolate_in_time(pi, Wi, dt, dWi);
  hydro_gradients_extrapolate(pi, xij_i, dWi);
  float dWj[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
  hydro_gradients_extrapolate_in_time(pj, Wj, dt, dWj);
  hydro_gradients_extrapolate(pj, xij_j, dWj);

  /* Apply the slope limiter at this interface */
  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  Wi[0] += dWi[0];
  Wi[1] += dWi[1];
  Wi[2] += dWi[2];
  Wi[3] += dWi[3];
  Wi[4] += dWi[4];

  Wj[0] += dWj[0];
  Wj[1] += dWj[1];
  Wj[2] += dWj[2];
  Wj[3] += dWj[3];
  Wj[4] += dWj[4];
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
