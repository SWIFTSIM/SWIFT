#include <assert.h>
#include <math.h>
#define INLINE inline

/* Mock dependencies */

const float Y0 = 200e6;

const int test_mat_id = 0;

float method_yield_weakening_thermal_xi(void) {
  return 1.2f;
}

float material_T_melt(const int mat_id) {
  return 933.f;
}

float material_rho_0(const int mat_id) {
  return 2700.f;
}

float method_yield_weakening_density_mult_param(void) {
  return 0.85f;
}

float method_yield_weakening_density_pow_param(void) {
  return 4.f;
}

/* Include functions to test */
#include "../../src/strength/strength_yield_stress_weakening.h"

static void test_thermal_weakening(
    const float *temperatures,
    const float *expected,
    int n,
    float frac_tol)
{
  for (int i = 0; i < n; i++) {

    const float T = temperatures[i];

    /* Reset Y */
    float Y = Y0;

    /* Apply temperature weakening to Y */
    yield_weakening_apply_temperature_to_yield_stress(&Y, test_mat_id, T);

    /* Check result */
    float scale = fmaxf(fabsf(expected[i]), Y0);
    assert(fabsf(Y - expected[i]) <= frac_tol * scale);
  }
}

static void test_density_weakening(
    const float *densities,
    const float *expected,
    int n,
    float frac_tol)
{
  for (int i = 0; i < n; i++) {
    const float rho = densities[i];
    
    /* Reset Y */
    float Y = Y0;

    /* Apply density weakening to Y */
    yield_weakening_apply_density_to_yield_stress(&Y, test_mat_id, rho);

    /* Check result */
    float scale = fmaxf(fabsf(expected[i]), Y0);
    assert(fabsf(Y - expected[i]) <= frac_tol * scale);
  }
}