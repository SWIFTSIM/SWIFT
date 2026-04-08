#include <assert.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy constants. */
#define mat_phase_solid  1
#define mat_phase_fluid  0

/* Dummy material parameters. */
static float material_Y_0(const int mat_id) { return 200e6f; }

/* Dummy weakening. */
static int density_weakening_calls = 0;
static int temperature_weakening_calls = 0;
static void yield_weakening_apply_density_to_yield_stress(
    float *yield_stress, const int mat_id, const float density)
{
  density_weakening_calls++;
}

static void yield_weakening_apply_temperature_to_yield_stress(
    float *yield_stress, const int mat_id, const float density, const float u)
{
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

#define STRENGTH_YIELD_STRESS_BENZ_ASPHAUG
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


/* Test that damage does not affect the yield stress: damaged yield stress
 * equals fully intact yield stress regardless of damage value. */
static void test_damaged_yield_stress_ignores_damage(void)
{
  const float Y_intact  = 1e6f;
  const float Y_damaged = 0.f;

  const float tol = 1e-6f;
  assert(within_tol(yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.f),  Y_intact, tol));
  assert(within_tol(yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 0.5f), Y_intact, tol));
  assert(within_tol(yield_compute_damaged_yield_stress(Y_intact, Y_damaged, 1.f),  Y_intact, tol));
}


/* Test that yield stress for solid equals Y_0 and for non-solid equals zero. */
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

  /* Negative pressure should still return Y_0 for solid. */
  assert(within_tol(
    yield_compute_yield_stress_fully_intact(mat_id, mat_phase_solid, -1e6f),
    Y_0, tol));
}

/* Test that yield stress for solid equals Y_0 and for non-solid equals zero. */
static void test_yield_compute_yield_stress(void)
{
  const int mat_id = 0;
  const float pressure = 1e5f;
  const float density = 1.f;
  const float u = 1e6f;
  const float damage = 0.6f;
  const float tol = 1e-6f;

  const float Y0 = material_Y_0(mat_id);

  /* Solid should return intact yield stress regardless of damage. */
  const float Y_solid =
    yield_compute_yield_stress(mat_id, mat_phase_solid,
                               pressure, density, u, damage);

  assert(within_tol(Y_solid, Y0, tol));

  /* Check weakening calls */
  assert(density_weakening_calls == 1);
  assert(temperature_weakening_calls == 1);

  /* Fluid should return zero. */
  const float Y_fluid =
    yield_compute_yield_stress(mat_id, mat_phase_fluid,
                               pressure, density, u, damage);

  assert(within_tol(Y_fluid, 0.f, tol));

  /* Check weakening calls */
  assert(density_weakening_calls == 1);
  assert(temperature_weakening_calls == 1);
}

/* Test yield_apply_yield_stress_to_sym_matrix:
 * - When J_2 is below the yield limit, f = Y^2 / (3*J_2) < 1 and M is scaled down.
 * - When J_2 is above the yield limit, f = 1 and M is unchanged. */
static void test_apply_yield_stress(void)
{
  const float tol = 1e-6f;

  /* Case 1: stress exceeds yield => f < 1, M should be scaled down.
   * Use a purely deviatoric tensor with known J_2. */
  struct sym_matrix S;
  S.xx = 4.f; S.yy = -1.f; S.zz = -3.f;
  S.xy = 0.f; S.xz =  0.f; S.yz =  0.f;

  /* J_2 = 0.5*(16 + 1 + 9) = 13. */
  const float J_2 = 0.5f * (S.xx*S.xx + S.yy*S.yy + S.zz*S.zz);
  const float yield_stress = 1.f;
  const float expected_f = fminf((yield_stress * yield_stress) / (3.f * J_2), 1.f);

  struct sym_matrix M = S;
  yield_apply_yield_stress_to_sym_matrix(&M, S, 0.f, 0.f, yield_stress);

  assert(within_tol(M.xx, expected_f * S.xx, tol));
  assert(within_tol(M.yy, expected_f * S.yy, tol));
  assert(within_tol(M.zz, expected_f * S.zz, tol));
  assert(within_tol(M.xy, expected_f * S.xy, tol));
  assert(within_tol(M.xz, expected_f * S.xz, tol));
  assert(within_tol(M.yz, expected_f * S.yz, tol));

  /* Case 2: stress within yield => f = 1, M unchanged. */
  struct sym_matrix S2;
  S2.xx = 0.001f; S2.yy = -0.0005f; S2.zz = -0.0005f;
  S2.xy = 0.f; S2.xz = 0.f; S2.yz = 0.f;

  struct sym_matrix M2 = S2;
  yield_apply_yield_stress_to_sym_matrix(&M2, S2, 0.f, 0.f, yield_stress);

  assert(within_tol(M2.xx, S2.xx, tol));
  assert(within_tol(M2.yy, S2.yy, tol));
  assert(within_tol(M2.zz, S2.zz, tol));
  assert(within_tol(M2.xy, S2.xy, tol));
  assert(within_tol(M2.xz, S2.xz, tol));
  assert(within_tol(M2.yz, S2.yz, tol));

  /* Case3: trivial case. */
  struct sym_matrix S0 = {0};
  yield_apply_yield_stress_to_sym_matrix(&S0, S0, 0.f, 0.f, yield_stress);
  assert(within_tol(S0.xx, 0.f, tol));
  assert(within_tol(S0.yy, 0.f, tol));
  assert(within_tol(S0.zz, 0.f, tol));
  assert(within_tol(S0.xy, 0.f, tol));
  assert(within_tol(S0.xz, 0.f, tol));
  assert(within_tol(S0.yz, 0.f, tol));

  /* Check that 3 * J2_new <= Y^2 */
  const float J2_new = strength_compute_deviatoric_sym_matrix_J_2(M);
  assert(3.f * J2_new - yield_stress * yield_stress <= tol * yield_stress * yield_stress);
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