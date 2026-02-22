#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy particle. */
struct strength_data {
  float principal_stress_eigen[3];
};

struct part {
  float h;
  struct strength_data strength_data;
};

#define STRENGTH_ARTIFICIAL_STRESS_BASIS_INDP
#include "../../src/strength/strength_artificial_stress.h"
#include "../../src/symmetric_matrix.h"


/* Initialise a 3x3 matrix to zero. */
static void zero_tensor(float T[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      T[i][j] = 0.f;
    }
  }
}

/* Copy a 3x3 matrix. */
static void copy_tensor(float dst[3][3], const float src[3][3]) {
  memcpy(dst, src, 3 * 3 * sizeof(float));
}


/* Test no artificial stress when particle separation is large. */
static void test_no_stress_large_separation(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  /* All positive principal stresses so stress would be applied if separation
   * were small. */
  pi.strength_data.principal_stress_eigen[0] = 1.f;
  pi.strength_data.principal_stress_eigen[1] = 2.f;
  pi.strength_data.principal_stress_eigen[2] = 3.f;
  pj.strength_data.principal_stress_eigen[0] = 1.f;
  pj.strength_data.principal_stress_eigen[1] = 2.f;
  pj.strength_data.principal_stress_eigen[2] = 3.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3], pairwise_stress_tensor_j_before[3][3];

  /* Arbitrary stress tensors. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = (float)(i + j + 1);
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);
  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  /* Large separation: r/h large compared with eta_crit. */
  const float r = 10.f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(pairwise_stress_tensor_i[i][j] - pairwise_stress_tensor_i_before[i][j]) <= tol);
      assert(fabsf(pairwise_stress_tensor_j[i][j] - pairwise_stress_tensor_j_before[i][j]) <= tol);
    }
  }
}


/* Test no artificial stress is applied when all principal stresses are non-positive. */
static void test_no_stress_no_pos_eigen(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  /* All non-positive principal stresses. */
  pi.strength_data.principal_stress_eigen[0] = -3.f;
  pi.strength_data.principal_stress_eigen[1] = -1.f;
  pi.strength_data.principal_stress_eigen[2] =  0.f;
  pj.strength_data.principal_stress_eigen[0] = -5.f;
  pj.strength_data.principal_stress_eigen[1] = -2.f;
  pj.strength_data.principal_stress_eigen[2] =  0.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3], pairwise_stress_tensor_j_before[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = (float)(i + j + 1);
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);
  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  /* Close separation so factor would be non-zero. */
  const float r = 0.f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(pairwise_stress_tensor_i[i][j] - pairwise_stress_tensor_i_before[i][j]) <= tol);
      assert(fabsf(pairwise_stress_tensor_j[i][j] - pairwise_stress_tensor_j_before[i][j]) <= tol);
    }
  }
}

/* Test diagonal correction properties: off-diagonals are unchanged, all three
 * diagonal corrections are equal, and correction equals the maximum principal stress at r=0. */
static void test_correction(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3];

  zero_tensor(pairwise_stress_tensor_i);
  zero_tensor(pairwise_stress_tensor_j);
  pairwise_stress_tensor_i[0][0] = 4.f;
  pairwise_stress_tensor_i[1][1] = 2.f;
  pairwise_stress_tensor_i[2][2] = 1.f;
  pairwise_stress_tensor_i[0][1] = pairwise_stress_tensor_i[1][0] = 5.f;
  pairwise_stress_tensor_i[0][2] = pairwise_stress_tensor_i[2][0] = 3.f;
  pairwise_stress_tensor_i[1][2] = pairwise_stress_tensor_i[2][1] = 2.f;

  /* Compute eigenvalues of the tensor we just constructed. */
  struct sym_matrix original;
  get_sym_matrix_from_matrix(&original, pairwise_stress_tensor_i);
  float eigenvalues_before[3];
  sym_matrix_compute_eigenvalues(eigenvalues_before, original);

  /* Find max eigenvalue. */
  float max_eigen = eigenvalues_before[0];
  if (eigenvalues_before[1] > max_eigen) { max_eigen = eigenvalues_before[1]; }
  if (eigenvalues_before[2] > max_eigen) { max_eigen = eigenvalues_before[2]; }

  pi.strength_data.principal_stress_eigen[0] = eigenvalues_before[0];
  pi.strength_data.principal_stress_eigen[1] = eigenvalues_before[1];
  pi.strength_data.principal_stress_eigen[2] = eigenvalues_before[2];

  pj.strength_data.principal_stress_eigen[0] = -1.f;
  pj.strength_data.principal_stress_eigen[1] = -1.f;
  pj.strength_data.principal_stress_eigen[2] = -1.f;

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);

  const float tol = 1e-6f;
  const float r = 0.f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
      pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  /* Off-diagonals must be unchanged. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (i != j) {
        assert(fabsf(pairwise_stress_tensor_i[i][j] -
                     pairwise_stress_tensor_i_before[i][j]) <= tol);
      }
    }
  }

  /* All diagonal corrections must equal max_eigen. */
  for (int i = 0; i < 3; i++) {
    const float corr = pairwise_stress_tensor_i_before[i][i] -
                       pairwise_stress_tensor_i[i][i];
    assert(fabsf(corr - max_eigen) <= tol);
  }

  /* Compute eigenvalues of corrected tensor and verify all are <= 0. */
  struct sym_matrix corrected;
  get_sym_matrix_from_matrix(&corrected, pairwise_stress_tensor_i);
  float eigenvalues_after[3];
  sym_matrix_compute_eigenvalues(eigenvalues_after, corrected);

  for (int i = 0; i < 3; i++) {
    assert(eigenvalues_after[i] <= 0.f);
  }
}

/* Test that particles i and j are treated independently. */
static void test_particles_treated_independently(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  /* Only pi has tensile principal stress. */
  pi.strength_data.principal_stress_eigen[0] = 5.f;
  pi.strength_data.principal_stress_eigen[1] = 1.f;
  pi.strength_data.principal_stress_eigen[2] = 2.f;

  pj.strength_data.principal_stress_eigen[0] = -3.f;
  pj.strength_data.principal_stress_eigen[1] = -1.f;
  pj.strength_data.principal_stress_eigen[2] = -2.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_j_before[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = 1.f;
    }
  }

  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  const float r = 0.f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  const float tol = 1e-6f;

  /* pairwise_stress_tensor_j should be completely unchanged. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(pairwise_stress_tensor_j[i][j] - pairwise_stress_tensor_j_before[i][j]) <= tol);
    }
  }

  /* pairwise_stress_tensor_i diagonals should have been reduced. */
  for (int i = 0; i < 3; i++) {
    assert(pairwise_stress_tensor_i[i][i] < 1.f);
  }
}

int main(void)
{
  test_no_stress_large_separation();
  test_no_stress_no_pos_eigen();
  test_correction();
  test_particles_treated_independently();

  return 0;
}