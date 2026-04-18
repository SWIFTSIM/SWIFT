#include <assert.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy constants. */
#define mat_phase_solid  1
#define mat_phase_fluid  0

/* Dummy material parameters. */
static float material_Y_0(const int mat_id)  { return 200e6f; }
static float material_Y_M(const int mat_id)  { return 600e6f; }
static float material_mu_i(const int mat_id) { return 2.0f; }
static float material_mu_d(const int mat_id) { return 0.8f; }

/* Dummy weakening. */
static int density_weakening_calls = 0;
static int temperature_weakening_calls = 0;
static void yield_weakening_apply_density_to_yield_stress(
    float *yield_stress, const int mat_id, const float density) {
  density_weakening_calls++;
}

static void yield_weakening_apply_temperature_to_yield_stress(
    float *yield_stress, const int mat_id, const float density, const float u) {
  temperature_weakening_calls++;
}

/* Dummy sym_matrix. */
struct sym_matrix {
  union {
    struct { float elements[6]; };
    struct { float xx, yy, zz, xy, xz, yz; };
  };
};

/* Dummy J_2 calculation. */
static float __attribute__((unused)) strength_compute_deviatoric_sym_matrix_J_2(
    struct sym_matrix M) {
  return 0.5f * (M.xx*M.xx + M.yy*M.yy + M.zz*M.zz)
       + M.xy*M.xy + M.xz*M.xz + M.yz*M.yz;
}

#define INLINE inline

#define STRENGTH_YIELD_STRESS_COLLINS
#include "../../src/strength/strength_yield_stress.h"

/* Test that damage does not affect the deviatoric stress tensor. */
static void test_damaged_deviatoric_stress_tensor_unchanged(void) {
  struct sym_matrix S;
  S.xx =  3.f; S.yy = -1.f; S.zz = -2.f;
  S.xy =  1.f; S.xz =  0.5f; S.yz = -0.5f;

  const float tol = 1e-6f;

  /* Tensor must be unchanged for any damage value. */
  struct sym_matrix S_d0   = yield_compute_damaged_deviatoric_stress_tensor(S, 0.f);
  struct sym_matrix S_d05  = yield_compute_damaged_deviatoric_stress_tensor(S, 0.5f);
  struct sym_matrix S_d1   = yield_compute_damaged_deviatoric_stress_tensor(S, 1.f);

  assert(fabsf(S_d0.xx -  S.xx) <= tol);
  assert(fabsf(S_d05.xx - S.xx) <= tol);
  assert(fabsf(S_d1.xx -  S.xx) <= tol);

  assert(fabsf(S_d0.yy -  S.yy) <= tol);
  assert(fabsf(S_d05.yy - S.yy) <= tol);
  assert(fabsf(S_d1.yy -  S.yy) <= tol);

  assert(fabsf(S_d0.zz -  S.zz) <= tol);
  assert(fabsf(S_d05.zz - S.zz) <= tol);
  assert(fabsf(S_d1.zz -  S.zz) <= tol);

  assert(fabsf(S_d0.xy -  S.xy) <= tol);
  assert(fabsf(S_d05.xy - S.xy) <= tol);
  assert(fabsf(S_d1.xy -  S.xy) <= tol);

  assert(fabsf(S_d0.xz -  S.xz) <= tol);
  assert(fabsf(S_d05.xz - S.xz) <= tol);
  assert(fabsf(S_d1.xz -  S.xz) <= tol);

  assert(fabsf(S_d0.yz -  S.yz) <= tol);
  assert(fabsf(S_d05.yz - S.yz) <= tol);
  assert(fabsf(S_d1.yz -  S.yz) <= tol);
}


/* Test that damaged yield stress linearly interpolates between intact and
 * damaged values. */
static void test_damaged_yield_stress_interpolation(void) {
  const float Y_intact  = 1e6f;
  const float Y_damaged = 2e5f;
  const float tol = 1e-6f;

  /* At damage=0 should return fully intact. */
  assert(fabsf(
    yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.f) -
    Y_intact) <= tol);

  /* At damage=1 should return fully damaged. */
  assert(fabsf(
    yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 1.f) -
    Y_damaged) <= tol);

  /* At damage=0.5 should return midpoint. */
  const float expected_mid = 0.5f * Y_intact + 0.5f * Y_damaged;
  assert(fabsf(
    yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.5f) -
    expected_mid) <= tol);
}


/* Test yield stress of fully intact material:
 * - Returns Y_0 at zero or negative pressure.
 * - Increases with positive pressure, saturating towards Y_M.
 * - Returns 0 for non-solid phase. */
static void test_yield_stress_fully_intact(void) {
  const int mat_id = 0;
  const float Y_0  = material_Y_0(mat_id);
  const float Y_M  = material_Y_M(mat_id);
  const float mu_i = material_mu_i(mat_id);
  const float tol  = 1e-3f;

  /* Non-solid returns zero. */
  assert(fabsf(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_fluid, 1e5f)) <= tol);

  /* Zero pressure returns Y_0. */
  assert(fabsf(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, 0.f) -
    Y_0) <= tol);

  /* Positive pressure: verify against analytical formula. */
  const float pressure = 1e5f;
  const float expected = Y_0 + mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
  assert(fabsf(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, pressure) -
    expected) <= tol);

  /* Yield stress must increase with pressure. */
  const float Y_low  = yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, 1e4f);
  const float Y_high = yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, 1e6f);
  assert(Y_high > Y_low);

  /* Yield stress must stay below Y_M. */
  const float Y_very_high = yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, 1e12f);
  assert(Y_very_high < Y_M);
}


/* Test yield stress of fully damaged material:
 * - Returns 0 for non-positive pressure.
 * - Returns min(mu_d, Y_intact) for positive pressure.
 * - Returns 0 for non-solid phase. */
static void test_yield_stress_fully_damaged(void) {
  const int mat_id = 0;
  const float mu_d = material_mu_d(mat_id);
  const float tol  = 1e-6f;

  /* Non-solid returns zero. */
  assert(fabsf(
    yield_compute_yield_stress_fully_damaged(mat_id, mat_phase_fluid, 1e5f)) <= tol);

  /* Zero or negative pressure returns zero. */
  assert(fabsf(
    yield_compute_yield_stress_fully_damaged(mat_id, mat_phase_solid, 0.f)) <= tol);
  assert(fabsf(
    yield_compute_yield_stress_fully_damaged(mat_id, mat_phase_solid, -1e5f)) <= tol);

  /* Positive pressure: result must be <= mu_d and <= Y_intact. */
  const float pressure = 1e5f;
  const float Y_damaged = yield_compute_yield_stress_fully_damaged(mat_id, mat_phase_solid, pressure);
  const float Y_intact  = yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, pressure);
  const float expected = fminf(mu_d * pressure, Y_intact);
  assert(fabsf(Y_damaged - expected) <= tol);
}

/* Test full yield stress computation:
 * - Returns 0 for non-solid phase.
 * - Reduces to fully intact for damage=0.
 * - Reduces to fully damaged for damage=1.
 * - Interpolates correctly for intermediate damage. */
static void test_yield_compute_yield_stress(void) {
  const int mat_id = 0;
  const float pressure = 1e5f;
  const float density = 1.f;
  const float u = 1e6f;
  const float tol = 1e-6f;

  /* Non-solid gives 0 */
  assert(fabsf(
    yield_compute_yield_stress(mat_id, mat_phase_fluid,
                               density, pressure, u, 0.5f)) <= tol);

  /* Check weakening calls */
  assert(density_weakening_calls == 0);
  assert(temperature_weakening_calls == 0);

  /* Compute reference intact and damaged values */
  const float Y_intact =
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, pressure);

  const float Y_damaged =
    yield_compute_yield_stress_fully_damaged(mat_id, mat_phase_solid, pressure);

  /* damage = 0 gives fully intact */
  const float Y0 =
    yield_compute_yield_stress(mat_id, mat_phase_solid,
                              density, pressure, u, 0.f);

  assert(fabsf(Y0 - Y_intact) <= tol);

  /* Check weakening calls */
  assert(density_weakening_calls == 1);
  assert(temperature_weakening_calls == 1);

  /* damage = 1 gives fully damaged */
  const float Y1 =
    yield_compute_yield_stress(mat_id, mat_phase_solid,
                              density, pressure, u, 1.f);

  assert(fabsf(Y1 - Y_damaged) <= tol);

  /* Check weakening calls */
  assert(density_weakening_calls == 2);
  assert(temperature_weakening_calls == 2);

  /* Intermediate damage gives linear interpolation */
  const float damage = 0.3f;
  const float expected =
    (1.f - damage) * Y_intact + damage * Y_damaged;

  const float Yd =
    yield_compute_yield_stress(mat_id, mat_phase_solid,
                              density, pressure, u, damage);

  assert(fabsf(Yd - expected) <= tol);

  /* Check weakening calls */
  assert(density_weakening_calls == 3);
  assert(temperature_weakening_calls == 3);
}

/* Test yield_apply_yield_stress_to_sym_matrix. */
static void test_apply_yield_stress(void) {
  const float tol = 1e-6f;

  /* Case 1: stress exceeds yield => f < 1, M should be scaled down. */
  struct sym_matrix S;
  S.xx = 4.f; S.yy = -1.f; S.zz = -3.f;
  S.xy = 0.f; S.xz =  0.f; S.yz =  0.f;

  /* J_2 = 0.5*(16 + 1 + 9) = 13. */
  const float J_2 = 13.f;
  const float yield_stress = 1.f;
  const float expected_f = fminf(yield_stress / sqrtf(J_2), 1.f);

  struct sym_matrix M = S;
  yield_apply_yield_stress_to_sym_matrix(&M, S, 0.f, 0.f, yield_stress);

  assert(fabsf(M.xx - expected_f * S.xx) <= tol);
  assert(fabsf(M.yy - expected_f * S.yy) <= tol);
  assert(fabsf(M.zz - expected_f * S.zz) <= tol);
  assert(fabsf(M.xy - expected_f * S.xy) <= tol);
  assert(fabsf(M.xz - expected_f * S.xz) <= tol);
  assert(fabsf(M.yz - expected_f * S.yz) <= tol);

  /* Case 2: stress within yield => f = 1, M unchanged. */
  struct sym_matrix S2;
  S2.xx = 0.001f; S2.yy = -0.0005f; S2.zz = -0.0005f;
  S2.xy = 0.f; S2.xz = 0.f; S2.yz = 0.f;

  struct sym_matrix M2 = S2;
  yield_apply_yield_stress_to_sym_matrix(&M2, S2, 0.f, 0.f, yield_stress);

  assert(fabsf(M2.xx - S2.xx) <= tol);
  assert(fabsf(M2.yy - S2.yy) <= tol);
  assert(fabsf(M2.zz - S2.zz) <= tol);
  assert(fabsf(M2.xy - S2.xy) <= tol);
  assert(fabsf(M2.xz - S2.xz) <= tol);
  assert(fabsf(M2.yz - S2.yz) <= tol);

  /* Case3: trivial case. */
  struct sym_matrix S0 = {0};
  yield_apply_yield_stress_to_sym_matrix(&S0, S0, 0.f, 0.f, yield_stress);
  assert(fabsf(S0.xx) <= tol);
  assert(fabsf(S0.yy) <= tol);
  assert(fabsf(S0.zz) <= tol);
  assert(fabsf(S0.xy) <= tol);
  assert(fabsf(S0.xz) <= tol);
  assert(fabsf(S0.yz) <= tol);

  /* Check that sqrt(J2_new) <= Y */
  const float J2_new = strength_compute_deviatoric_sym_matrix_J_2(M);
  assert(sqrtf(J2_new) - yield_stress <= tol * yield_stress);
}


int main(void) {
  test_damaged_deviatoric_stress_tensor_unchanged();
  test_damaged_yield_stress_interpolation();
  test_yield_stress_fully_intact();
  test_yield_stress_fully_damaged();
  test_yield_compute_yield_stress();
  test_apply_yield_stress();

  return 0;
}