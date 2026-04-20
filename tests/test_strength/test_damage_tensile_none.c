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
    float activation_thresholds[40];
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

/* Include the "none" tensile damage implementation */
#define STRENGTH_DAMAGE
#define STRENGTH_DAMAGE_TENSILE_NONE
#include "../../src/strength/strength_damage.h"

/* Test getters and setters always return zero */
static void test_getters_setters(void) {
    struct part p = {0};
    struct xpart xp = {0};
    const float tol = 1e-6f;

    assert(fabsf(damage_get_tensile_damage(&p)) <= tol);
    assert(fabsf(damage_get_tensile_damage_full(&xp)) <= tol);

    damage_set_tensile_damage(&p, 0.8f);
    damage_set_tensile_damage_full(&xp, 0.9f);

    assert(fabsf(damage_get_tensile_damage(&p)) <= tol);
    assert(fabsf(damage_get_tensile_damage_full(&xp)) <= tol);
}

/* Test compute_cbrtD_dt always returns zero */
static void test_compute_cbrtD_dt(void) {
    struct part p = {0};
    p.strength_data.number_of_flaws = 3;
    p.strength_data.activation_thresholds[0] = 0.1f;
    p.strength_data.activation_thresholds[1] = 0.2f;
    p.strength_data.activation_thresholds[2] = 0.3f;

    struct sym_matrix stress = {0};
    stress.xx = 1.f; stress.yy = -0.5f; stress.zz = 0.2f;

    float cbrtD_dt = 1.f;
    int activated_flaws = -1;

    damage_tensile_compute_cbrtD_dt(&cbrtD_dt, &activated_flaws, p.strength_data.number_of_flaws,
                                    p.strength_data.activation_thresholds, stress, 0, 1.f, 1.f, 0.f);

                                    
    const float tol = 1e-6f;
    assert(fabsf(cbrtD_dt) <= tol);
    assert(activated_flaws == 0);
}

/* Test apply timestep does not change damage */
static void test_apply_timestep(void) {
    float tensile_damage = 0.f;
    const int activated_flaws = 2;
    const int number_of_flaws = 4;
    const float dt = 0.1f;
    const float cbrtD_dt = 1.f;

    damage_tensile_apply_timestep_to_tensile_damage(&tensile_damage, cbrtD_dt, activated_flaws, number_of_flaws, dt);                                    
    const float tol = 1e-6f;
    assert(fabsf(tensile_damage) <= tol);
}

/* Test evolve function does not change damage */
static void test_evolve_function(void) {
    struct part p = {0};
    struct sym_matrix stress = {0};
    stress.xx = 0.25f; stress.yy = 0.1f; stress.zz = 0.f;

    float tensile_damage = 0.f;

    damage_tensile_evolve(&tensile_damage, &p, stress, 0, 1.f, 1.f, 0.f, 0.1f);
                                    
    const float tol = 1e-6f;
    assert(fabsf(tensile_damage) <= tol);
}

int main(void) {
    test_getters_setters();
    test_compute_cbrtD_dt();
    test_apply_timestep();
    test_evolve_function();

    return 0;
}