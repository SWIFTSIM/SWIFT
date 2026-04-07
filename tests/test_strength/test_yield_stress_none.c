#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy constants. */
#define mat_phase_solid  1
#define mat_phase_fluid  0

/* Dummy material parameters. */
static float material_Y_0(const int mat_id) { return FLT_MAX; }

/* Dummy no weakening. */
static void __attribute__((unused)) yield_weakening_apply_density_to_yield_stress(
    float *yield_stress, const int mat_id, const float density) {}

static void __attribute__((unused)) yield_weakening_apply_temperature_to_yield_stress(
    float *yield_stress, const int mat_id, const float density, const float temperature) {}

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

/* Don't define a yield stress scheme */
#include "../../src/strength/strength_yield_stress.h"

static INLINE int within_tol(float a, float b, float rel_tol)
{
  const float magnitude = fmaxf(fabsf(a), fabsf(b));
  return fabsf(a - b) <= rel_tol * magnitude;
}


/* Test that yield_compute_damaged_deviatoric_stress_tensor scales all elements
 * by (1 - damage). */
static void test_damaged_deviatoric_stress_tensor(void)
{
  struct sym_matrix S;
  S.xx =  3.f; S.yy = -1.f; S.zz = -2.f;
  S.xy =  1.f; S.xz =  0.5f; S.yz = -0.5f;

  const float damage = 0.4f;
  const float scale = 1.f - damage;

  struct sym_matrix S_damaged = yield_compute_damaged_deviatoric_stress_tensor(S, damage);

  const float tol = 1e-6f;
  assert(within_tol(S_damaged.xx, scale * S.xx, tol));
  assert(within_tol(S_damaged.yy, scale * S.yy, tol));
  assert(within_tol(S_damaged.zz, scale * S.zz, tol));
  assert(within_tol(S_damaged.xy, scale * S.xy, tol));
  assert(within_tol(S_damaged.xz, scale * S.xz, tol));
  assert(within_tol(S_damaged.yz, scale * S.yz, tol));

  /* Test 0 damage case. */
  struct sym_matrix S_d0 = yield_compute_damaged_deviatoric_stress_tensor(S, 0.f);
  assert(within_tol(S_d0.xx, S.xx, tol));
  assert(within_tol(S_d0.yy, S.yy, tol));
  assert(within_tol(S_d0.zz, S.zz, tol));
  assert(within_tol(S_d0.xy, S.xy, tol));
  assert(within_tol(S_d0.xz, S.xz, tol));
  assert(within_tol(S_d0.yz, S.yz, tol));
}


/* Test that fully damaged material (damage=1) gives a zero deviatoric stress
 * tensor, i.e. acts as a fluid. */
static void test_damaged_deviatoric_stress_tensor_fully_damaged(void)
{
  struct sym_matrix S;
  S.xx =  3.f; S.yy = -1.f; S.zz = -2.f;
  S.xy =  1.f; S.xz =  0.5f; S.yz = -0.5f;

  struct sym_matrix S_damaged = yield_compute_damaged_deviatoric_stress_tensor(S, 1.f);

  const float tol = 1e-6f;
  assert(within_tol(S_damaged.xx, 0.f, tol));
  assert(within_tol(S_damaged.yy, 0.f, tol));
  assert(within_tol(S_damaged.zz, 0.f, tol));
  assert(within_tol(S_damaged.xy, 0.f, tol));
  assert(within_tol(S_damaged.xz, 0.f, tol));
  assert(within_tol(S_damaged.yz, 0.f, tol));
}

/* Damage does not affect yield stress */
static void test_damaged_yield_stress_ignores_damage(void)
{
  const float Y_intact  = 1e6f;
  const float Y_damaged = 0.f;

  const float tol = 1e-6f;
  assert(within_tol(
      yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.f),
      Y_intact, tol));

  assert(within_tol(
      yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.5f),
      Y_intact, tol));

  assert(within_tol(
      yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 1.f),
      Y_intact, tol));
}

/* Fully intact yield stress is FLT_MAX */
static void test_yield_stress_fully_intact(void)
{
  const int mat_id = 0;
  const float Y_0 = material_Y_0(mat_id);
  const float pressure = 1e5f;
  const float tol = 1e-6f;

  assert(within_tol(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, pressure),
    Y_0, tol));

  assert(within_tol(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_fluid, pressure),
    0.f, tol));

  /* Negative pressure still FLT_MAX */
  assert(within_tol(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, -1e6f),
    Y_0, tol));
}


/* Combined yield stress */
static void test_yield_compute_yield_stress(void)
{
  const int mat_id = 0;
  const float pressure = 1e5f;
  const float density = 1.f;
  const float temperature = 300.f;
  const float damage = 0.6f;
  const float tol = 1e-6f;

  const float Y0 = FLT_MAX;

  /* Solid: FLT_MAX */
  const float Y_solid =
    yield_compute_yield_stress(mat_id, mat_phase_solid,
                               density, pressure, temperature, damage);

  assert(within_tol(Y_solid, Y0, tol));

  /* Fluid: 0 */
  const float Y_fluid =
    yield_compute_yield_stress(mat_id, mat_phase_fluid,
                               density, pressure, temperature, damage);

  assert(within_tol(Y_fluid, 0.f, tol));
}


/* Apply yield stress should do nothing */
static void test_apply_yield_stress(void)
{
  const float tol = 1e-6f;

  struct sym_matrix S;
  S.xx = 4.f; S.yy = -1.f; S.zz = -3.f;
  S.xy = 0.2f; S.xz = -0.1f; S.yz = 0.3f;

  struct sym_matrix M = S;

  yield_apply_yield_stress_to_sym_matrix(
      &M, S, 1.f, 1.f, FLT_MAX);

  /* Must remain unchanged */
  assert(within_tol(M.xx, S.xx, tol));
  assert(within_tol(M.yy, S.yy, tol));
  assert(within_tol(M.zz, S.zz, tol));
  assert(within_tol(M.xy, S.xy, tol));
  assert(within_tol(M.xz, S.xz, tol));
  assert(within_tol(M.yz, S.yz, tol));

  /* Zero tensor remains zero */
  struct sym_matrix S0 = {0};
  yield_apply_yield_stress_to_sym_matrix(
      &S0, S0, 0.f, 0.f, FLT_MAX);

  assert(within_tol(S0.xx, 0.f, tol));
  assert(within_tol(S0.yy, 0.f, tol));
  assert(within_tol(S0.zz, 0.f, tol));
  assert(within_tol(S0.xy, 0.f, tol));
  assert(within_tol(S0.xz, 0.f, tol));
  assert(within_tol(S0.yz, 0.f, tol));
}


int main(void)
{
  test_damaged_deviatoric_stress_tensor();
  test_damaged_deviatoric_stress_tensor_fully_damaged();
  test_damaged_yield_stress_ignores_damage();
  test_yield_stress_fully_intact();
  test_yield_compute_yield_stress();
  test_apply_yield_stress();

  return 0;
}