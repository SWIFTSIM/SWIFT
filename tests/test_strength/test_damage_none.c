#include <assert.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy structs. */
struct p_strength_data {
  float damage;
  float tensile_damage;
  float shear_damage;
  float damage_accumulation_timescale;
  int number_of_flaws;
  float activation_thresholds[40]; //### hardcoded length
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

/* Getters always return 0 */
static void test_getters(void) {
  struct part p = {0};
  struct xpart xp = {0};

  const float tol = 1e-6f;
  assert(fabsf(strength_get_damage(&p)) <= tol);
  assert(fabsf(strength_get_damage_full(&xp)) <= tol);
}

/* Setters do nothing, fields remain unmodified */
static void test_setters(void) {
  struct part p = {0};
  struct xpart xp = {0};

  p.strength_data.damage = 0.5f;
  xp.strength_data.damage_full = 0.7f;

  strength_set_damage(&p, 1.f);
  strength_set_damage_full(&xp, 1.f);

  const float tol = 1e-6f;
  assert(fabsf(p.strength_data.damage - 0.5f) <= tol);
  assert(fabsf(xp.strength_data.damage_full - 0.7f) <= tol);
}

/* Damage timestep limiter does nothing */
static void test_timestep_damage(void) {
  struct part p = {0};
  float dt = 1.f;
  strength_compute_timestep_damage(&dt, &p);

  const float tol = 1e-6f;
  assert(fabsf(dt - 1.f) <= tol);
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
  const float tol = 1e-6f;
  assert(fabsf(stress_tensor.xx) <= tol);
  assert(fabsf(stress_tensor.yy) <= tol);
  assert(fabsf(stress_tensor.zz) <= tol);
  assert(fabsf(stress_tensor.xy) <= tol);
  assert(fabsf(stress_tensor.xz) <= tol);
  assert(fabsf(stress_tensor.yz) <= tol);
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
  const float tol = 1e-6f;
  assert(fabsf(damage - 0.3f) <= tol);
  assert(fabsf(tensile - 0.2f) <= tol);
  assert(fabsf(shear - 0.1f) <= tol);
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
  const float tol = 1e-6f;
  assert(fabsf(p.strength_data.damage - 0.5f) <= tol);
  assert(fabsf(p.strength_data.tensile_damage - 0.6f) <= tol);
  assert(fabsf(p.strength_data.shear_damage - 0.7f) <= tol);
  assert(fabsf(xp.strength_data.damage_full - 0.8f) <= tol);
  assert(fabsf(xp.strength_data.tensile_damage_full - 0.9f) <= tol);
  assert(fabsf(xp.strength_data.shear_damage_full - 1.0f) <= tol);
}

int main(void) {
  test_getters();
  test_setters();
  test_timestep_damage();
  test_stress_tensor();
  test_evolve_functions();
  test_first_init();

  return 0;
}