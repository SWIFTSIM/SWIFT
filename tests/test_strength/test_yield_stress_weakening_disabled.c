#include <assert.h>
#include <math.h>
#define INLINE inline

/* DO NOT Enable yield stress weakening methods. i.e. don't do */
/* #define STRENGTH_YIELD_STRESS_WEAKENING_THERMAL */
/* #define STRENGTH_YIELD_STRESS_WEAKENING_DENSITY */

#include "test_yield_stress_weakening.c"

static void test_thermal_weakening_disabled(void)
{
  /* Set constants */
  const float T_melt = material_T_melt(test_mat_id);

   /* Test temperatures in each regime */
  const float temperatures[4] = {
      0,                /* zero temperature edge case */
      0.5f * T_melt,   /* below T_melt */
      T_melt,          /* at T_melt */
      2.0f * T_melt    /* above T_melt */
  };

  /* Expected results for each temperature */
  const float expected[4] = {
      Y0,   /* zero temperature edge case */
      Y0,   /* below T_melt */
      Y0,   /* at T_melt */
      Y0,   /* above T_melt */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(temperatures) / sizeof(temperatures[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  test_thermal_weakening(temperatures, expected, n, frac_tol);
}

static void test_density_weakening_disabled(void)
{
  /* Set constants */
  const float rho_0    = material_rho_0(test_mat_id);
  const float a        = method_yield_weakening_density_mult_param();
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
      Y0,   /* zero density edge case */
      Y0,   /* below rho_weak */
      Y0,   /* at rho_weak */
      Y0    /* above rho_weak */
  };

  /* Number of test cases and fractional tolerance for comparisons */
  const int n = sizeof(densities) / sizeof(densities[0]);
  const float frac_tol = 1e-6f;

  /* Run test */
  test_density_weakening(densities, expected, n, frac_tol);
}

int main(void)
{
  test_thermal_weakening_disabled();
  test_density_weakening_disabled();
  return 0;
}