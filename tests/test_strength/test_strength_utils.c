#include <assert.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy sym_matrix. */
struct sym_matrix {
  union {
    struct { float elements[6]; };
    struct { float xx, yy, zz, xy, xz, yz; };
  };
};

#include "../../src/strength/strength_utilities.h"

/* Test J2 invariant for a deviatoric tensor. */
static void test_J2_deviatoric(void) {
  struct sym_matrix M = {0};

  /* Arbitrary components with zero trace */
  M.xx =  4.f;
  M.yy = -1.f;
  M.zz = -3.f;
  M.xy =  2.f;
  M.xz = -1.f;
  M.yz =  0.5f;

  /* Manual computation of J2 */
  const float expected =
      0.5f * (M.xx*M.xx + M.yy*M.yy + M.zz*M.zz)
    + (M.xy*M.xy + M.xz*M.xz + M.yz*M.yz);

  const float J2 = strength_compute_deviatoric_sym_matrix_J_2(M);

  const float tol = 1e-6f;
  assert(fabsf(J2 - expected) <= tol);
}


/* Test J2 invariant for hydrostatic tensor equals zero. */
static void test_J2_hydrostatic(void) {
  struct sym_matrix M = {0};

  const float p = 10.f;
  M.xx = p;
  M.yy = p;
  M.zz = p;

  const float J2 = strength_compute_sym_matrix_J_2(M);

  const float tol = 1e-6f ;
  assert(fabsf(J2) <= tol);
}


/* Test consistency between J2 from full tensor and J2 from deviatoric. */
static void test_J2_consistency(void) {
  struct sym_matrix M = {0};

  /* Arbitrary symmetric tensor with non-zero trace */
  M.xx =  4.f;
  M.yy = -1.f;
  M.zz =  2.f;
  M.xy =  1.f;
  M.xz = -0.5f;
  M.yz =  0.25f;

  /* Compute J2 from full tensor */
  const float J2_full = strength_compute_sym_matrix_J_2(M);

  /* Manually construct deviatoric tensor from M */
  const float trace = M.xx + M.yy + M.zz;
  struct sym_matrix M_dev = M;
  M_dev.xx -= trace / 3.f;
  M_dev.yy -= trace / 3.f;
  M_dev.zz -= trace / 3.f;

  /* Compute deviatoric J2 */
  const float J2_dev = strength_compute_deviatoric_sym_matrix_J_2(M_dev);

  const float tol = 1e-6f;
  assert(fabsf(J2_full - J2_dev) <= tol);
}


/* Test strain-rate tensor. */
static void test_strain_rate_tensor(void) {
  float dv[3][3] = {
    {1.f, 2.f, 4.f},
    {4.f, 2.f, 6.f},
    {2.f, 8.f, 3.f}
  };

  float strain_rate_tensor[3][3];
  strength_compute_strain_rate_tensor(strain_rate_tensor, dv);

  const float tol = 1e-6f;

  /* Diagonal elements passed through unchanged */
  assert(fabsf(strain_rate_tensor[0][0] - dv[0][0]) <= tol);
  assert(fabsf(strain_rate_tensor[1][1] - dv[1][1]) <= tol);
  assert(fabsf(strain_rate_tensor[2][2] - dv[2][2]) <= tol);

  /* Off-diagonal symmetry */
  assert(fabsf(strain_rate_tensor[0][1] - strain_rate_tensor[1][0]) <= tol);
  assert(fabsf(strain_rate_tensor[0][2] - strain_rate_tensor[2][0]) <= tol);
  assert(fabsf(strain_rate_tensor[1][2] - strain_rate_tensor[2][1]) <= tol);

  /* Correct off-diagonal values */
  assert(fabsf(strain_rate_tensor[0][1] - 0.5f * (dv[0][1] + dv[1][0])) <= tol);
  assert(fabsf(strain_rate_tensor[0][2] - 0.5f * (dv[0][2] + dv[2][0])) <= tol);
  assert(fabsf(strain_rate_tensor[1][2] - 0.5f * (dv[1][2] + dv[2][1])) <= tol);
}


/* Test rotation-rate tensor. */
static void test_rotation_rate_tensor(void) {
  float dv[3][3] = {
    {1.f, 2.f, 4.f},
    {4.f, 2.f, 6.f},
    {2.f, 8.f, 3.f}
  };

  float rotation_rate_tensor[3][3];
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, dv);

  const float tol = 1e-6f;

  /* Diagonal should be zero */
  assert(fabsf(rotation_rate_tensor[0][0]) <= tol);
  assert(fabsf(rotation_rate_tensor[1][1]) <= tol);
  assert(fabsf(rotation_rate_tensor[2][2]) <= tol);

  /* Antisymmetry */
  assert(fabsf(rotation_rate_tensor[0][1] + rotation_rate_tensor[1][0]) <= tol);
  assert(fabsf(rotation_rate_tensor[0][2] + rotation_rate_tensor[2][0]) <= tol);
  assert(fabsf(rotation_rate_tensor[1][2] + rotation_rate_tensor[2][1]) <= tol);

  /* Correct off-diagonal values */
  assert(fabsf(rotation_rate_tensor[0][1] - 0.5f * (dv[1][0] - dv[0][1])) <= tol);
  assert(fabsf(rotation_rate_tensor[0][2] - 0.5f * (dv[2][0] - dv[0][2])) <= tol);
  assert(fabsf(rotation_rate_tensor[1][2] - 0.5f * (dv[2][1] - dv[1][2])) <= tol);
}


/* Test rotation-rate tensor vanishes for symmetric dv. */
static void test_rotation_rate_tensor_symmetric_input(void) {
  float dv[3][3] = {
    {1.f, 2.f, 3.f},
    {2.f, 4.f, 5.f},
    {3.f, 5.f, 6.f}
  };

  float rotation_rate_tensor[3][3];
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, dv);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(rotation_rate_tensor[i][j]) <= tol);
    }
  }
}


/* Test rotation term M*R - R*M for identity M gives zero. */
static void test_rotation_term_identity_M(void) {
  float R[3][3] = {
    { 0.f, -1.f,  0.f},
    { 1.f,  0.f,  0.f},
    { 0.f,  0.f,  0.f}
  };
  float M[3][3] = {0};
  float rotation_term[3][3];

  for (int i = 0; i < 3; i++) {
    M[i][i] = 1.f;
  }

  strength_compute_rotation_term(rotation_term, R, M);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(rotation_term[i][j]) <= tol);
    }
  }
}

/* Test that rigid-body rotation generates no new stress. */
static void test_rigid_body_rotation_no_stress_generation(void) {
    /* Initial stress set to zero everywhere */
    float S[3][3] = {
        {0.f, 0.f, 0.f},
        {0.f, 0.f, 0.f},
        {0.f, 0.f, 0.f}
    };

    /* Rigid rotation around z-axis */
    const float omega = 1.0f;
    float dv[3][3] = {
        { 0.f, -omega, 0.f },
        { omega,  0.f, 0.f },
        { 0.f,    0.f, 0.f }
    };

    const float mu = 100.f;
    float strain_rate[3][3];
    float rotation_rate[3][3];
    float rotation_term[3][3];
    float dS_dt[3][3];

    strength_compute_strain_rate_tensor(strain_rate, dv);
    strength_compute_rotation_rate_tensor(rotation_rate, dv);
    strength_compute_rotation_term(rotation_term, rotation_rate, S);

    /* Calculate dS/dt */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dS_dt[i][j] = 2.f * mu * strain_rate[i][j] + rotation_term[i][j];
        }
    }

    /* Verify that rigid rotation does not generate stress */
    const float tol = 1e-6f;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            assert(fabsf(dS_dt[i][j]) <= tol);
        }
    }
}

/* Test stress rotation to detect sign errors */
static void test_stress_rotation_45deg(void) {
    /* Initial stress */
    float S[3][3] = {
        {1.f, 0.f, 0.f},
        {0.f, 0.f, 0.f},
        {0.f, 0.f, 0.f}
    };

    /* Rigid rotation around z-axis */
    const float omega = 0.5f;
    float dv[3][3] = {
        { 0.f, -omega, 0.f },
        { omega,  0.f, 0.f },
        { 0.f,    0.f, 0.f }
    };

    /* Multi-step integration */
    const int nsteps = 1000000;
    const float dt_total = (M_PI / 4.f) / omega;
    const float dt_step = dt_total / nsteps;

    float rotation_rate[3][3];
    float rotation_term[3][3];
    float S_rot[3][3];
    memcpy(S_rot, S, sizeof(S));
    for (int step = 0; step < nsteps; step++) {
        /* Compute tensors */
        strength_compute_rotation_rate_tensor(rotation_rate, dv);
        strength_compute_rotation_term(rotation_term, rotation_rate, S_rot);

        /* Update stress */
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                S_rot[i][j] += rotation_term[i][j] * dt_step;
            }
        }
    }

    /* Check final values after pi/4 rotation */
    const float tol = 1e-3f;
    assert(fabsf(S_rot[0][0] - 0.5f) <= tol);
    assert(fabsf(S_rot[1][1] - 0.5f) <= tol);
    assert(fabsf(S_rot[0][1] - 0.5f) <= tol);
    assert(fabsf(S_rot[1][0] - 0.5f) <= tol);
}

/* Test that rotating stress by 2*pi returns to the original. */
static void test_stress_rotation_full_circle(void) {
  float S[3][3] = {
    {2.f, 1.f, 0.f},
    {1.f, 0.f, 0.f},
    {0.f, 0.f, 0.5f}
  };
  float S_init[3][3];
  memcpy(S_init, S, sizeof(S));

  const float omega = 0.5f;
  float dv[3][3] = {
    { 0.f, -omega, 0.f},
    { omega,  0.f, 0.f},
    { 0.f,    0.f, 0.f}
  };

  const int nsteps = 1000000;
  const float dt_total = (2.f * (float)M_PI) / omega;
  const float dt_step = dt_total / nsteps;

  float rotation_rate[3][3];
  float rotation_term[3][3];

  strength_compute_rotation_rate_tensor(rotation_rate, dv);

  for (int step = 0; step < nsteps; step++) {
    strength_compute_rotation_term(rotation_term, rotation_rate, S);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        S[i][j] += rotation_term[i][j] * dt_step;
      }
    }
  }

  const float tol = 1e-3f;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(S[i][j] - S_init[i][j]) <= tol);
    }
  }
}


int main(void) {
  test_J2_deviatoric();
  test_J2_hydrostatic();
  test_J2_consistency();
  test_strain_rate_tensor();
  test_rotation_rate_tensor();
  test_rotation_rate_tensor_symmetric_input();
  test_rotation_term_identity_M();
  test_rigid_body_rotation_no_stress_generation();
  test_stress_rotation_45deg();
  test_stress_rotation_full_circle();

  return 0;
}