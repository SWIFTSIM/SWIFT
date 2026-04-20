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

/* Dummy sym_matrix. */
struct sym_matrix {
  union {
    struct { float elements[6]; };
    struct { float xx, yy, zz, xy, xz, yz; };
  };
};

/* Dummy dependencies. */

/* Tensile / shear getters + setters */
static INLINE float damage_get_tensile_damage(const struct part *p) {
  return p->strength_data.tensile_damage;
}
static INLINE float damage_get_shear_damage(const struct part *p) {
  return p->strength_data.shear_damage;
}
static INLINE float damage_get_tensile_damage_full(const struct xpart *xp) {
  return xp->strength_data.tensile_damage_full;
}
static INLINE float damage_get_shear_damage_full(const struct xpart *xp) {
  return xp->strength_data.shear_damage_full;
}

static INLINE void damage_set_tensile_damage(struct part *p, float damage) {
  p->strength_data.tensile_damage = damage;
}
static INLINE void damage_set_shear_damage(struct part *p, float damage) {
  p->strength_data.shear_damage = damage;
}
static INLINE void damage_set_tensile_damage_full(struct xpart *xp, float damage) {
  xp->strength_data.tensile_damage_full = damage;
}
static INLINE void damage_set_shear_damage_full(struct xpart *xp, float damage) {
  xp->strength_data.shear_damage_full = damage;
}

/* Dummy evolution models. */
static INLINE void damage_tensile_evolve(
    float *tensile_damage, struct part *p, struct sym_matrix stress_tensor,
    int mat_id, float mass, float density, float damage, float dt) {
  *tensile_damage += fminf(0.1f * dt, 1.f);
}

static INLINE void damage_shear_evolve(
    float *shear_damage, struct part *p,
    int mat_id, float density, float u) {
  *shear_damage = fminf(*shear_damage + 0.05f, 1.f);
}

static INLINE void damage_tensile_compute_cbrtD_dt(
    float *cbrtD_dt, int *number_of_activated_flaws,
    int number_of_flaws, const float *activation_thresholds,
    struct sym_matrix stress_tensor,
    int mat_id, float mass, float density, float damage) {
    *cbrtD_dt = 2.f;
    *number_of_activated_flaws = 1; // dummy positive
}

/*  Include core damage scheme. */
#include "../../src/strength/damage/damage_core.h"

/* Basic getter and setter consistency */
static void test_damage_get_set(void) {
  struct part p = {0};
  struct xpart xp = {0};

  strength_set_damage(&p, 0.3f);
  strength_set_damage_full(&xp, 0.7f);

  const float tol = 1e-6f;
  assert(fabsf(strength_get_damage(&p) - 0.3f) <= tol);
  assert(fabsf(strength_get_damage_full(&xp) - 0.7f) <= tol);
}


/* Damage timestep limiter behaves as expected */
static void test_damage_timestep_limiter(void) {
  struct part p = {0};
  p.strength_data.damage_accumulation_timescale = 0.1f;

  float dt = 1.f;
  float dt_before = dt;

  strength_compute_timestep_damage(&dt, &p);

  /* dt should be reduced. */
  assert(dt < dt_before);
}


/* If already small enough, dt should remain unchanged */
static void test_damage_timestep_no_change(void) {
  struct part p = {0};
  p.strength_data.damage_accumulation_timescale = 1000.f;

  float dt = 1.f;
  float dt_before = dt;

  strength_compute_timestep_damage(&dt, &p);

  const float tol = 1e-6f;
  assert(fabsf(dt - dt_before) <= tol);
}


/* Stress tensor with positive pressure */
static void test_damage_stress_tensor_positive_pressure(void) {
  struct sym_matrix damaged_deviatoric_stress_tensor;
  damaged_deviatoric_stress_tensor.xx =  3.f; damaged_deviatoric_stress_tensor.yy = -1.f; damaged_deviatoric_stress_tensor.zz = -2.f;
  damaged_deviatoric_stress_tensor.xy =  1.f; damaged_deviatoric_stress_tensor.xz =  0.5f; damaged_deviatoric_stress_tensor.yz = -0.5f;
  struct sym_matrix stress_tensor;

  const float pressure = 10.f;
  const float damage = 0.5f;

  damage_compute_stress_tensor(&stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  /* Subtract pressure */
  const float tol = 1e-6f;
  assert(fabsf(stress_tensor.xx - (damaged_deviatoric_stress_tensor.xx - pressure)) <= tol);
  assert(fabsf(stress_tensor.yy - (damaged_deviatoric_stress_tensor.yy - pressure)) <= tol);
  assert(fabsf(stress_tensor.zz - (damaged_deviatoric_stress_tensor.zz - pressure)) <= tol);
}


/* Stress tensor with negative pressure weakened by damage */
static void test_damage_stress_tensor_negative_pressure(void) {
  struct sym_matrix damaged_deviatoric_stress_tensor;
  damaged_deviatoric_stress_tensor.xx =  3.f; damaged_deviatoric_stress_tensor.yy = -1.f; damaged_deviatoric_stress_tensor.zz = -2.f;
  damaged_deviatoric_stress_tensor.xy =  1.f; damaged_deviatoric_stress_tensor.xz =  0.5f; damaged_deviatoric_stress_tensor.yz = -0.5f;
  struct sym_matrix stress_tensor;

  const float pressure = -10.f;
  const float damage = 0.5f;

  damage_compute_stress_tensor(&stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);

  const float expected_subtract = (1.f - damage) * pressure;

  const float tol = 1e-6f;
  assert(fabsf(stress_tensor.xx - (damaged_deviatoric_stress_tensor.xx - expected_subtract)) <= tol);
  assert(fabsf(stress_tensor.yy - (damaged_deviatoric_stress_tensor.yy - expected_subtract)) <= tol);
  assert(fabsf(stress_tensor.zz - (damaged_deviatoric_stress_tensor.zz - expected_subtract)) <= tol);
}


/* Damage evolution combines tensile + shear and clamps to 1 */
static void test_damage_evolve_combination(void) {
  struct part p = {0};

  float damage = 0.f;
  float tensile = 0.f;
  float shear = 0.f;

  struct sym_matrix stress_tensor = {0};

  damage_evolve(&damage, &tensile, &shear,
                &p, stress_tensor, 0, 1.f, 1.f, 1.f, 1.f);

  /* tensile = 0.1, shear = 0.05 gives total = 0.15 */
  const float tol = 1e-6f;
  assert(fabsf(damage - 0.15f) <= tol);
}


/* Clamp at damage = 1 */
static void test_damage_evolve_clamp(void) {
  struct part p = {0};

  float damage = 0.9f;
  float tensile = 0.9f;
  float shear = 0.5f;

  struct sym_matrix stress_tensor = {0};

  damage_evolve(&damage, &tensile, &shear,
                &p, stress_tensor, 0, 1.f, 1.f, 1.f, 1.f);

  const float tol = 1e-6f;
  assert(fabsf(damage - 1.f) <= tol);
}


/* Initialisation sets everything to zero */
static void test_first_init(void) {
  struct part p = {0};
  struct xpart xp = {0};

  /* Fill with non-zero values to check overwrite */
  p.strength_data.damage = 0.5f;
  p.strength_data.tensile_damage = 0.6f;
  p.strength_data.shear_damage = 0.7f;

  xp.strength_data.damage_full = 0.8f;
  xp.strength_data.tensile_damage_full = 0.9f;
  xp.strength_data.shear_damage_full = 1.0f;

  strength_first_init_part_damage(&p, &xp);

  const float tol = 1e-6f;
  assert(fabsf(p.strength_data.damage - 0.f) <= tol);
  assert(fabsf(p.strength_data.tensile_damage - 0.f) <= tol);
  assert(fabsf(p.strength_data.shear_damage - 0.f) <= tol);
  assert(fabsf(xp.strength_data.damage_full - 0.f) <= tol);
  assert(fabsf(xp.strength_data.tensile_damage_full - 0.f) <= tol);
  assert(fabsf(xp.strength_data.shear_damage_full - 0.f) <= tol);
}

int main(void) {
  test_damage_get_set();
  test_damage_timestep_limiter();
  test_damage_timestep_no_change();
  test_damage_stress_tensor_positive_pressure();
  test_damage_stress_tensor_negative_pressure();
  test_damage_evolve_combination();
  test_damage_evolve_clamp();
  test_first_init();

  return 0;
}