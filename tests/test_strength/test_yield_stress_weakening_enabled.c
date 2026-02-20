#include <assert.h>
#include <math.h>
#define INLINE inline

/* Enable yield stress weakening methods */
#define STRENGTH_YIELD_STRESS_WEAKENING_THERMAL
#define STRENGTH_YIELD_STRESS_WEAKENING_DENSITY

#include "test_yield_stress_weakening.c"

static void test_thermal_weakening_enabled(void)
{
  /* Set constants */
  const float xi     = method_yield_weakening_thermal_xi();
  const float T_melt = material_T_melt(test_mat_id);

   /* Test temperatures in each regime */
  const float temperatures[4] = {
      0,               /* zero temperature edge case */
      0.5f * T_melt,   /* below T_melt */
      T_melt,          /* at T_melt */
      2.0f * T_melt    /* above T_melt */
  };

  /* Expected results for each temperature */
  const float expected[4] = {
      Y0,                                                  /* zero temperature edge case */
      Y0 * tanhf(xi * (T_melt / temperatures[1] - 1.f)),   /* below T_melt */
      Y0 * tanhf(xi * (T_melt / temperatures[2] - 1.f)),   /* at T_melt */
      Y0 * tanhf(xi * (T_melt / temperatures[3] - 1.f))    /* above T_melt */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(temperatures) / sizeof(temperatures[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  test_thermal_weakening(temperatures, expected, n, frac_tol);
}

static void test_density_weakening_enabled(void)
{
  /* Set constants */
  const float rho_0    = material_rho_0(test_mat_id);
  const float a        = method_yield_weakening_density_mult_param();
  const float b        = method_yield_weakening_density_pow_param();
  const float rho_weak = a * rho_0;

  /* Test densities in each regime */
  const float densities[4] = {
      0,                 /* zero density edge case */
      0.5f * rho_weak,   /* below rho_weak */
      rho_weak,          /* at rho_weak */
      2.0f * rho_weak    /* above rho_weak */
  };

  /* Expected results for each density */
  const float expected[4] = {
      0,                                       /* zero density edge case */
      Y0 * powf(densities[1] / rho_weak, b),   /* below rho_weak */
      Y0,                                      /* at rho_weak */
      Y0                                       /* above rho_weak */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(densities) / sizeof(densities[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  test_density_weakening(densities, expected, n, frac_tol);
}

int main(void)
{
  test_thermal_weakening_enabled();
  test_density_weakening_enabled();
  return 0;
}