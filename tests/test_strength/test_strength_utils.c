#include <math.h>
#include <string.h>

#define INLINE inline

#include "../../src/strength/strength_utilities.h"

/* Test J2 invariant for a deviatoric tensor. */
static void test_J2_deviatoric(void)
{
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

  const float tol = 1e-6f * expected;

  assert(fabsf(J2 - expected) <= tol);
}


/* Test J2 invariant for hydrostatic tensor equals zero. */
static void test_J2_hydrostatic(void)
{
  struct sym_matrix M = {0};

  const float p = 10.f;
  M.xx = p;
  M.yy = p;
  M.zz = p;

  const float J2 = strength_compute_sym_matrix_J_2(M);

  const float tol = 1e-6f * p * p;

  assert(fabsf(J2) <= tol);
}


/* Test consistency between J2 from full tensor and J2 from deviatoric. */
static void test_J2_consistency(void)
{
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

  const float tol = 1e-6f * fabsf(J2_full);

  assert(fabsf(J2_full - J2_dev) <= tol);
}


/* Test strain-rate tensor. */
static void test_strain_rate_tensor(void)
{
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
static void test_rotation_rate_tensor(void)
{
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
static void test_rotation_rate_tensor_symmetric_input(void)
{
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
static void test_rotation_term_identity_M(void)
{
  float R[3][3] = {
    { 0.f, -1.f,  0.f},
    { 1.f,  0.f,  0.f},
    { 0.f,  0.f,  0.f}
  };
  float M[3][3] = {0};
  float rot_term[3][3];

  for (int i = 0; i < 3; i++) {
    M[i][i] = 1.f;
  }

  strength_compute_rotation_term(rot_term, R, M);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(rot_term[i][j]) <= tol);
    }
  }
}


/* Test rotation term M*R - R*M for an arbitrary case. */
static void test_rotation_term_arbitrary(void)
{
  /* Arbitrary antisymmetric R */
  float R[3][3] = {
    { 0.f, -2.f,  1.f},
    { 2.f,  0.f, -3.f},
    {-1.f,  3.f,  0.f}
  };

  /* Arbitrary M */
  float M[3][3] = {
      {2.f,  1.f, -0.5f},
      {3.f,  5.f,  0.8f},
      {0.1f, -1.f,  3.f}
  };

  float rot_term[3][3];
  strength_compute_rotation_term(rot_term, R, M);

  /* Manual calculation of rotation term */
  float MR[3][3] = {{0.f}};
  float RM[3][3] = {{0.f}};
  float expected[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        MR[i][j] += M[i][k] * R[k][j];
        RM[i][j] += R[i][k] * M[k][j];
      }
      expected[i][j] = MR[i][j] - RM[i][j];
    }
  }

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(rot_term[i][j] - expected[i][j]) <= tol);
    }
  }
}


/* Test strain tensor evolution. */
static void test_evolve_strain_tensor(void)
{
  float strain_tensor[3][3] = {
    {1.f,  0.5f, 0.2f},
    {0.5f, 2.f,  0.3f},
    {0.2f, 0.3f, 1.5f}
  };
  float strain_rate[3][3] = {
    {0.1f, 0.4f, 0.2f},
    {0.4f, 0.2f, 0.1f},
    {0.2f, 0.1f, 0.3f}
  };
  float rot_term[3][3] = {
    { 0.f, -0.1f,  0.05f},
    { 0.1f, 0.f,  -0.02f},
    {-0.05f, 0.02f, 0.f}
  };

  const float dt = 0.1f;

  /* Manually compute expected result */
  float expected[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      expected[i][j] = strain_tensor[i][j] + (strain_rate[i][j] + rot_term[i][j]) * dt;
    }
  }

  strength_evolve_strain_tensor(strain_tensor, strain_rate, rot_term, dt);

  const float tol = 1e-6f;

  /* Diagonal elements */
  assert(fabsf(strain_tensor[0][0] - expected[0][0]) <= tol);
  assert(fabsf(strain_tensor[1][1] - expected[1][1]) <= tol);
  assert(fabsf(strain_tensor[2][2] - expected[2][2]) <= tol);

  /* Symmetry of result */
  assert(fabsf(strain_tensor[0][1] - strain_tensor[1][0]) <= tol);
  assert(fabsf(strain_tensor[0][2] - strain_tensor[2][0]) <= tol);
  assert(fabsf(strain_tensor[1][2] - strain_tensor[2][1]) <= tol);
}



/* Test rotation tensor evolution preserves orthogonality. */
static void test_rotation_tensor_orthogonality(void)
{
  float R[3][3] = {0};
  float rotation_rate[3][3] = {0};

  for (int i = 0; i < 3; i++) {
    R[i][i] = 1.f;
  }

  rotation_rate[0][1] = -1.5f;
  rotation_rate[1][0] =  1.5f;

  strength_evolve_rotation_tensor(R, rotation_rate, 0.3f);

  /* Check R * R^T = I */
  float RRT[3][3] = {{0.f}};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        RRT[i][j] += R[i][k] * R[j][k];
      }
    }
  }

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const float expected = (i == j) ? 1.f : 0.f;
      assert(fabsf(RRT[i][j] - expected) <= tol);
    }
  }
}


/* Test rotation tensor evolution of an arbitrary tensor. */
static void test_rotation_of_arbitrary_tensor(void)
{
  float R[3][3] = {0};
  float rotation_rate[3][3] = {0};

  for (int i = 0; i < 3; i++) {
    R[i][i] = 1.f;
  }

  const float omega = 2.f;
  const float dt = 0.4f;
  const float theta = omega * dt;

  rotation_rate[0][1] = -omega;
  rotation_rate[1][0] =  omega;

  strength_evolve_rotation_tensor(R, rotation_rate, dt);

  float M[3][3] = {
    { 2.f,   1.f,  0.5f},
    { 1.f,   3.f, -0.2f},
    { 0.5f, -0.2f, 1.f }
  };

  /* Compute M' = R * M * R^T using evolved R */
  float RM[3][3] = {{0.f}};
  float Mprime[3][3] = {{0.f}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        RM[i][j] += R[i][k] * M[k][j];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        Mprime[i][j] += RM[i][k] * R[j][k];
      }
    }
  }

  /* Compute exact M' using analytic rotation matrix */
  float R_exact[3][3] = {
    { cosf(theta), -sinf(theta), 0.f},
    { sinf(theta),  cosf(theta), 0.f},
    { 0.f,          0.f,         1.f}
  };

  float RM_exact[3][3] = {{0.f}};
  float Mprime_exact[3][3] = {{0.f}};

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        RM_exact[i][j] += R_exact[i][k] * M[k][j];
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        Mprime_exact[i][j] += RM_exact[i][k] * R_exact[j][k];
      }
    }
  }

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(Mprime[i][j] - Mprime_exact[i][j]) <= tol);
    }
  }
}


/* Test rotation tensor evolution does nothing for negligibly small theta. */
static void test_rotation_evolution_zero_theta(void)
{
  float R[3][3] = {0};
  float rotation_rate[3][3] = {0};

  for (int i = 0; i < 3; i++) {
    R[i][i] = 1.f;
  }

  rotation_rate[0][1] = -1.f;
  rotation_rate[1][0] =  1.f;

  /* dt so small that early return should be triggered */
  const float dt = 1e-10f;

  float R_before[3][3];
  memcpy(R_before, R, sizeof(R));

  strength_evolve_rotation_tensor(R, rotation_rate, dt);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(R[i][j] - R_before[i][j]) <= tol);
    }
  }
}


int main(void)
{
  test_J2_deviatoric();
  test_J2_hydrostatic();
  test_J2_consistency();
  test_strain_rate_tensor();
  test_rotation_rate_tensor();
  test_rotation_rate_tensor_symmetric_input();
  test_rotation_term_identity_M();
  test_rotation_term_arbitrary();
  test_evolve_strain_tensor();
  test_rotation_tensor_orthogonality();
  test_rotation_of_arbitrary_tensor();
  test_rotation_evolution_zero_theta();

  return 0;
}