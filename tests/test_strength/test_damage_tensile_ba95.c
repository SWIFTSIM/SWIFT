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

/* Dummy material properties that give a Young's modulus of 1 */
static float material_shear_mod(int mat_id) { return 0.375f; }
static float material_bulk_mod(int mat_id)  { return 1.f; }

/* Dummy eigenvalue function (simplified: returns diagonal elements as principal stresses) */
static void sym_matrix_compute_eigenvalues(float out[3], struct sym_matrix stress) {
    out[0] = stress.xx; out[1] = stress.yy; out[2] = stress.zz;
}

/* Include the BA95 tensile damage implementation */
#define STRENGTH_DAMAGE
#define STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG
#include "../../src/strength/strength_damage.h"

static int within_tol(float a, float b, float tol) {
    const float scale = fmaxf(fabsf(a), fabsf(b));
    return fabsf(a - b) <= tol * scale;
}


/* Test getters and setters */
static void test_getters_setters(void) {
    struct part p = {0};
    struct xpart xp = {0};

    p.strength_data.tensile_damage = 0.3f;
    xp.strength_data.tensile_damage_full = 0.5f;

    assert(within_tol(damage_get_tensile_damage(&p), 0.3f, 1e-6f));
    assert(within_tol(damage_get_tensile_damage_full(&xp), 0.5f, 1e-6f));

    damage_set_tensile_damage(&p, 0.8f);
    damage_set_tensile_damage_full(&xp, 0.9f);

    assert(within_tol(p.strength_data.tensile_damage, 0.8f, 1e-6f));
    assert(within_tol(xp.strength_data.tensile_damage_full, 0.9f, 1e-6f));
}

/* Test compute_cbrtD_dt with zero flaws */
static void test_compute_cbrtD_dt_zero_flaws(void) {
    struct part p = {0};
    struct sym_matrix stress = {0};
    stress.xx = 1.f; stress.yy = 2.f; stress.zz = 3.f; // positive stress so not in tension

    float cbrtD_dt = 1.f;
    int activated_flaws = -1;

    damage_tensile_compute_cbrtD_dt(&cbrtD_dt, &activated_flaws, 0, p.strength_data.activation_thresholds,
                                    stress, 0, 1.f, 1.f, 0.f);

    /* No flaws so no damage accumulation */
    assert(within_tol(cbrtD_dt, 0.f, 1e-6f));
    assert(activated_flaws == 0);
}

/* Test compute_cbrtD_dt with flaws but compressive stress (max principal ≤ 0) */
static void test_compute_cbrtD_dt_compression_stress(void) {
    struct part p = {0};
    p.strength_data.number_of_flaws = 3;
    p.strength_data.activation_thresholds[0] = 0.1f;
    p.strength_data.activation_thresholds[1] = 0.2f;
    p.strength_data.activation_thresholds[2] = 0.3f;

    struct sym_matrix stress = {0};
    stress.xx = -1.f; stress.yy = -0.5f; stress.zz = -0.2f; // all compressive

    float cbrtD_dt = 1.f;
    int activated_flaws = -1;

    damage_tensile_compute_cbrtD_dt(&cbrtD_dt, &activated_flaws, p.strength_data.number_of_flaws,
                                    p.strength_data.activation_thresholds, stress, 0, 1.f, 1.f, 0.f);

    /* Compressive stress so no damage accumulation even though flaws exist */
    assert(within_tol(cbrtD_dt, 0.f, 1e-6f));
    assert(activated_flaws == 0);
}

/* Test compute_cbrtD_dt with flaws and max principal stress in tension */
static void test_compute_cbrtD_dt_active_flaws(void) {
    struct part p = {0};
    p.strength_data.number_of_flaws = 3;
    p.strength_data.activation_thresholds[0] = 0.1f;
    p.strength_data.activation_thresholds[1] = 0.2f;
    p.strength_data.activation_thresholds[2] = 0.5f;

    struct sym_matrix stress = {0};
    stress.xx = 0.3f; stress.yy = 0.1f; stress.zz = 0.f;

    float cbrtD_dt = 0.f;
    int activated_flaws = 0;

    damage_tensile_compute_cbrtD_dt(&cbrtD_dt, &activated_flaws, p.strength_data.number_of_flaws,
                                    p.strength_data.activation_thresholds, stress, 0, 1.f, 1.f, 0.f);

    assert(activated_flaws == 2); // thresholds 0.1 and 0.2 are exceeded
    assert(cbrtD_dt > 0.f);
}

/* Test apply timestep to tensile damage */
static void test_apply_timestep(void) {
    float tensile_damage = 0.f;
    const int activated_flaws = 8;
    const int number_of_flaws = 1000;
    const float dt = 0.1f;
    const float cbrtD_dt = 1.f;

    /* First timestep should increase but remain below max */
    damage_tensile_apply_timestep_to_tensile_damage(&tensile_damage, cbrtD_dt, activated_flaws, number_of_flaws, dt);
    float max_damage = (float)activated_flaws / (float)number_of_flaws; // = 0.008. See B&A99 for correction of cbrt placement here
    assert(within_tol(tensile_damage, 0.001f, 1e-6f));     // 0.1 * 1.f = 0.1 increase in cbrt(damage), so damage increases by 0.1^3 = 0.001
    assert(tensile_damage < max_damage); // still below max

    /* Second timestep cumulative is exactly max value */
    damage_tensile_apply_timestep_to_tensile_damage(&tensile_damage, cbrtD_dt, activated_flaws, number_of_flaws, dt);
    assert(within_tol(tensile_damage, max_damage, 1e-6f));

    /* Third timestep cumulative damage exceeds max, so should be capped to it */
    damage_tensile_apply_timestep_to_tensile_damage(&tensile_damage, cbrtD_dt, activated_flaws, number_of_flaws, dt);
    assert(within_tol(tensile_damage, max_damage, 1e-6f)); // ensure capped at maximum allowed
}

/* Test fully damaged particle */
static void test_full_damage(void) {
    float tensile_damage = 1.f;
    const int activated_flaws = 1;
    const int number_of_flaws = 2;
    const float dt = 0.1f;
    const float cbrtD_dt = 1.f;

    damage_tensile_apply_timestep_to_tensile_damage(&tensile_damage, cbrtD_dt, activated_flaws, number_of_flaws, dt);
    float max_damage = 1.f;  // fully damaged
    assert(within_tol(tensile_damage, max_damage, 1e-6f));
}

/* Test evolve function */
static void test_evolve_function(void) {
    struct part p = {0};
    p.strength_data.number_of_flaws = 3;
    p.strength_data.activation_thresholds[0] = 0.1f;
    p.strength_data.activation_thresholds[1] = 0.2f;
    p.strength_data.activation_thresholds[2] = 0.3f;

    struct sym_matrix stress = {0};
    stress.xx = 0.25f; stress.yy = 0.1f; stress.zz = 0.f;

    float tensile_damage = 0.f;

    damage_tensile_evolve(&tensile_damage, &p, stress, 0, 1.f, 1.f, 0.f, 0.1f);

    assert(tensile_damage > 0.f);
    assert(tensile_damage <= 1.f);
}

int main(void) {
    test_getters_setters();
    test_compute_cbrtD_dt_zero_flaws();
    test_compute_cbrtD_dt_compression_stress();
    test_compute_cbrtD_dt_active_flaws();
    test_full_damage();
    test_apply_timestep();
    test_evolve_function();

    return 0;
}


