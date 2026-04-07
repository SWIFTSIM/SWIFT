#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy structs. */
struct p_strength_data {
  float damage;
  float tensile_damage;
  float shear_damage;
  float dD_dt;
};

struct xp_strength_data {
  float damage_full;
  float tensile_damage_full;
  float shear_damage_full;
};

struct part {
  struct p_strength_data strength_data;
};

struct xpart {
  struct xp_strength_data strength_data;
};

/* Dummy sym_matrix */
struct sym_matrix {
  union {
    struct { float elements[6]; };
    struct { float xx, yy, zz, xy, xz, yz; };
  };
};

/* Include "none" damage functions */
#include "../../src/strength/damage/damage_none.h"

static int within_tol(float a, float b, float tol) {
  const float scale = fmaxf(fabsf(a), fabsf(b));
  return fabsf(a - b) <= tol * scale;
}

/* Getters always return 0 */
static void test_getters(void) {
  struct part p = {0};
  struct xpart xp = {0};

  assert(within_tol(strength_get_damage(&p), 0.f, 1e-6f));
  assert(within_tol(strength_get_damage_full(&xp), 0.f, 1e-6f));
}

/* Setters do nothing, fields remain unmodified */
static void test_setters(void) {
  struct part p = {0};
  struct xpart xp = {0};

  p.strength_data.damage = 0.5f;
  xp.strength_data.damage_full = 0.7f;

  strength_set_damage(&p, 1.f);
  strength_set_damage_full(&xp, 1.f);

  assert(within_tol(p.strength_data.damage, 0.5f, 1e-6f));
  assert(within_tol(xp.strength_data.damage_full, 0.7f, 1e-6f));
}

/* Damage timestep limiter does nothing */
static void test_timestep_damage(void) {
  struct part p = {0};
  float dt = 1.f;
  strength_compute_timestep_damage(&dt, &p);
  assert(within_tol(dt, 1.f, 1e-6f));
}

/* Stress tensor computation does nothing */
static void test_stress_tensor(void) {
  struct sym_matrix damaged_deviatoric_stress_tensor;
  damaged_deviatoric_stress_tensor.xx =  1.f; damaged_deviatoric_stress_tensor.yy = 2.f; damaged_deviatoric_stress_tensor.zz = 3.f;
  damaged_deviatoric_stress_tensor.xy =  4.f; damaged_deviatoric_stress_tensor.xz = 5.f; damaged_deviatoric_stress_tensor.yz = 6.f;
  struct sym_matrix stress_tensor = {0};

  const float pressure = 10.f;
  const float damage = 0.5f;

  damage_compute_stress_tensor(&stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  /* Should remain unchanged */
  assert(within_tol(stress_tensor.xx, 0.f, 1e-6f));
  assert(within_tol(stress_tensor.yy, 0.f, 1e-6f));
  assert(within_tol(stress_tensor.zz, 0.f, 1e-6f));
  assert(within_tol(stress_tensor.xy, 0.f, 1e-6f));
  assert(within_tol(stress_tensor.xz, 0.f, 1e-6f));
  assert(within_tol(stress_tensor.yz, 0.f, 1e-6f));
}

/* Evolution functions do nothing. */
static void test_evolve_functions(void) {
  struct part p = {0};
  struct xpart xp = {0};
  struct sym_matrix stress = {0};
  float damage = 0.3f, tensile = 0.2f, shear = 0.1f;

  damage_evolve(&damage, &tensile, &shear, &p, stress, 0, 1.f, 1.f, 1.f, 1.f);
  damage_predict_evolve(&p, stress, 0, 1.f, 1.f, 1.f, 1.f);
  damage_kick_evolve(&p, &xp, stress, 0, 1.f, 1.f, 1.f, 1.f);

  /* All values remain unmodified */
  assert(within_tol(damage, 0.3f, 1e-6f));
  assert(within_tol(tensile, 0.2f, 1e-6f));
  assert(within_tol(shear, 0.1f, 1e-6f));
}

/* dD/dt computation does nothing */
static void test_compute_dD_dt(void) {
  struct part p = {0};
  struct sym_matrix stress = {0};
  p.strength_data.dD_dt = 5.f;

  damage_compute_dD_dt(&p, stress, 0, 1.f, 1.f, 1.f);

  assert(within_tol(p.strength_data.dD_dt, 5.f, 1e-6f));
}

/* Initialisation does nothing */
static void test_first_init(void) {
  struct part p = {0};
  struct xpart xp = {0};

  /* Fill with non-zero values */
  p.strength_data.damage = 0.5f;
  p.strength_data.tensile_damage = 0.6f;
  p.strength_data.shear_damage = 0.7f;
  xp.strength_data.damage_full = 0.8f;
  xp.strength_data.tensile_damage_full = 0.9f;
  xp.strength_data.shear_damage_full = 1.0f;

  strength_first_init_part_damage(&p, &xp);

  /* Values should remain unchanged */
  assert(within_tol(p.strength_data.damage, 0.5f, 1e-6f));
  assert(within_tol(p.strength_data.tensile_damage, 0.6f, 1e-6f));
  assert(within_tol(p.strength_data.shear_damage, 0.7f, 1e-6f));
  assert(within_tol(xp.strength_data.damage_full, 0.8f, 1e-6f));
  assert(within_tol(xp.strength_data.tensile_damage_full, 0.9f, 1e-6f));
  assert(within_tol(xp.strength_data.shear_damage_full, 1.0f, 1e-6f));
}

int main(void) {
  test_getters();
  test_setters();
  test_timestep_damage();
  test_stress_tensor();
  test_evolve_functions();
  test_compute_dD_dt();
  test_first_init();

  return 0;
}