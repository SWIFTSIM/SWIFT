#include <assert.h>
#include <math.h>
#define INLINE inline

/* Enable yield stress weakening methods */
#define STRENGTH_YIELD_STRESS_WEAKENING_THERMAL
#define STRENGTH_YIELD_STRESS_WEAKENING_DENSITY

/* Mock dependencies */
const float Y0 = 200e6;

const int dummy_mat_id = 0;

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

/* Dummy temperature calculation */
static float gas_temperature_from_internal_energy(const float density, const float u,
                                     const int mat_id) {

  /* Return the value of the specific internal energy */
  return u;
}

/* Include functions to test */
#include "../../src/strength/strength_yield_stress_weakening.h"

static void test_thermal_weakening_enabled(void) {
  /* Set constants */
  const float xi     = method_yield_weakening_thermal_xi();
  const float T_melt = material_T_melt(dummy_mat_id);

   /* Test temperatures in each regime */
  const float temperatures[4] = {
      0.f,             /* zero temperature edge case */
      0.5f * T_melt,   /* below T_melt */
      T_melt,          /* at T_melt */
      2.0f * T_melt    /* above T_melt */
  };

  /* Expected results for each temperature */
  const float expected[4] = {
      Y0,                                                  /* zero temperature edge case */
      Y0 * tanhf(xi * (T_melt / temperatures[1] - 1.f)),   /* below T_melt */
      0.f,                                                 /* at T_melt */
      0.f                                                  /* above T_melt */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(temperatures) / sizeof(temperatures[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  for (int i = 0; i < n; i++) {

    /* Set internal energies to have same values as temperatures. */
    const float u = temperatures[i];
    const float density = 1.f;

    /* Reset Y */
    float Y = Y0;

    /* Apply temperature weakening to Y */
    yield_weakening_apply_temperature_to_yield_stress(&Y, dummy_mat_id, density, u);

    /* Check result */
    float scale = fmaxf(fabsf(expected[i]), Y0);
    assert(fabsf(Y - expected[i]) <= frac_tol * scale);
  }
}

static void test_density_weakening_enabled(void) {
  /* Set constants */
  const float rho_0    = material_rho_0(dummy_mat_id);
  const float a        = method_yield_weakening_density_mult_param();
  const float b        = method_yield_weakening_density_pow_param();
  const float rho_weak = a * rho_0;

  /* Test densities in each regime */
  const float densities[4] = {
      0.f,               /* zero density edge case */
      0.5f * rho_weak,   /* below rho_weak */
      rho_weak,          /* at rho_weak */
      2.0f * rho_weak    /* above rho_weak */
  };

  /* Expected results for each density */
  const float expected[4] = {
      0.f,                                     /* zero density edge case */
      Y0 * powf(densities[1] / rho_weak, b),   /* below rho_weak */
      Y0,                                      /* at rho_weak */
      Y0                                       /* above rho_weak */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(densities) / sizeof(densities[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  for (int i = 0; i < n; i++) {
    const float rho = densities[i];
    
    /* Reset Y */
    float Y = Y0;

    /* Apply density weakening to Y */
    yield_weakening_apply_density_to_yield_stress(&Y, dummy_mat_id, rho);

    /* Check result */
    float scale = fmaxf(fabsf(expected[i]), Y0);
    assert(fabsf(Y - expected[i]) <= frac_tol * scale);
  }
}

int main(void) {
  test_thermal_weakening_enabled();
  test_density_weakening_enabled();
  return 0;
}