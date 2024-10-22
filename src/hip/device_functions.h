#ifndef DEVICE_FUNCTIONS_H
#define DEVICE_FUNCTIONS_H
#include "../../config.h"

/* Local headers. */
//#include "../dimension.h"
//#include "../error.h"
//#include "../inline.h"
//#include "../minmax.h"
//#include "../vector.h"

// Is this even necessary? Probably not as our code will operate differently
#define num_cuda_threads 128
#define hydro_dimension 3.f

/// Here we define stuff from kernel_hydro.h when using cubic_spline_kernel.
/// Will worry about sorting 'if statements for different kernels later////
/* First some powers of gamma = H/h */
#define kernel_gamma ((float)(1.825742))
#define kernel_gamma_inv ((float)(1. / kernel_gamma))
#define kernel_gamma2 ((float)(kernel_gamma * kernel_gamma))
#define kernel_ivals 2
#define kernel_degree 3 /*!< Degree of the polynomial */
#define kernel_gamma_dim ((float)(kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_dim_plus_one                                              \
  ((float)(kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma))
#define kernel_gamma_inv_dim                                                   \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma)))
#define kernel_gamma_inv_dim_plus_one                                          \
  ((float)(1. / (kernel_gamma * kernel_gamma * kernel_gamma * kernel_gamma)))
#define kernel_ivals_f ((float)kernel_ivals) /*!< Number of branches */
#define kernel_constant ((float)(16. * M_1_PI))
/*! Cosmology default beta=3.0.
 * Alpha can be set in the parameter file.
 * Beta is defined as in e.g. Price (2010) Eqn (103) */
#define const_viscosity_beta 3.0f
#ifdef WITH_CUDA
extern "C" {
#endif
/**
 * @brief Returns the argument to the power given by the dimension plus one
 *
 * Computes \f$x^{d+1}\f$.
 */
__device__ float d_pow_dimension_plus_one(float x) {

#if defined(HYDRO_DIMENSION_3D)

  const float x2 = x * x;
  return x2 * x2;

#elif defined(HYDRO_DIMENSION_2D)

  return x * x * x;

#elif defined(HYDRO_DIMENSION_1D)

  return x * x;

#else

  error("The dimension is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power three adiabatic index minus five over
 * two.
 *
 * Computes \f$x^{(3\gamma - 5)/2}\f$.
 *
 * @param x Argument
 */
__device__ float d_pow_three_gamma_minus_five_over_two(float x) {
#if defined(HYDRO_GAMMA_5_3)

  return 1.f; /* x^(0) */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, -0.4f); /* x^(-2/5) */

#elif defined(HYDRO_GAMMA_4_3)

  return 1.f / sqrtf(x); /* x^(-1/2) */

#elif defined(HYDRO_GAMMA_2_1)

  return sqrtf(x); /* x^(1/2) */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Computes the kernel function and its derivative.
 *
 * The kernel function needs to be mutliplied by \f$h^{-d}\f$ and the gradient
 * by \f$h^{-(d+1)}\f$, where \f$d\f$ is the dimensionality of the problem.
 *
 * Returns 0 if \f$u > \gamma = H/h\f$.
 *
 * @param u The ratio of the distance to the smoothing length \f$u = x/h\f$.
 * @param W (return) The value of the kernel function \f$W(x,h)\f$.
 * @param dW_dx (return) The norm of the gradient of \f$|\nabla W(x,h)|\f$.
 */
__device__ void d_kernel_deval(float u, float *restrict W,
                               float *restrict dW_dx) {

  /* Go to the range [0,1[ from [0,H[ */
  const float x = u * kernel_gamma_inv;

  /* Pick the correct branch of the kernel */
  const int temp = (int)(x * kernel_ivals_f);
  const int ind = temp > kernel_ivals ? kernel_ivals : temp;
  static const float kernel_coeffs[(kernel_degree + 1) * (kernel_ivals + 1)] = {
      3.f,  -3.f, 0.f,  0.5f, /* 0 < u < 0.5 */
      -1.f, 3.f,  -3.f, 1.f,  /* 0.5 < u < 1 */
      0.f,  0.f,  0.f,  0.f}; /* 1 < u */
  const float *const coeffs = &kernel_coeffs[ind * (kernel_degree + 1)];
  /* First two terms of the polynomial ... */
  float w = coeffs[0] * x + coeffs[1];
  float dw_dx = coeffs[0];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_degree; k++) {
    dw_dx = dw_dx * x + w;
    w = x * w + coeffs[k];
  }

  w = max(w, 0.f);
  dw_dx = min(dw_dx, 0.f);

  /* Return everything */
  *W = w * kernel_constant * kernel_gamma_inv_dim;
  *dW_dx = dw_dx * kernel_constant * kernel_gamma_inv_dim_plus_one;
}

#ifdef WITH_CUDA
}
#endif

#endif // DEVICE_FUNCTIONS_H
