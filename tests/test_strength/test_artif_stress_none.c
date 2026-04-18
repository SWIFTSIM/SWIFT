#include <assert.h>
#include <math.h>
#include <string.h>

#define INLINE inline

/* Dummy particle. */
struct part {
  float h;
};

/* Don't define any artificial stress method */
#include "../../src/strength/strength_artificial_stress.h"

/* Copy a 3x3 matrix. */
static void copy_tensor(float dst[3][3], const float src[3][3]) {
  memcpy(dst, src, 3 * 3 * sizeof(float));
}

/* No artificial stress. */
static void test_artificial_stress_none(void) {
  struct part pi = {0}, pj = {0};

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3], pairwise_stress_tensor_j_before[3][3];

  /* Fill tensors with arbitrary non-zero values */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = (float)(i + j + 1); // Arbitrary non-zero values to check they are unchanged.
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);
  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  /* Test a range of r values */
  const float r_values[4] = {0.f, 0.3f, 1.f, 2.f};
  const float tol = 1e-6f;

  for (int k = 0; k < 4; k++) {

    /* Reset tensors each iteration */
    copy_tensor(pairwise_stress_tensor_i, pairwise_stress_tensor_i_before);
    copy_tensor(pairwise_stress_tensor_j, pairwise_stress_tensor_j_before);

    const float r = r_values[k];

    artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
        pairwise_stress_tensor_i,
        pairwise_stress_tensor_j,
        &pi, &pj, r);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        assert(fabsf(pairwise_stress_tensor_i[i][j] - pairwise_stress_tensor_i_before[i][j]) <= tol);
        assert(fabsf(pairwise_stress_tensor_j[i][j] - pairwise_stress_tensor_j_before[i][j]) <= tol);
      }
    }
  }
}


int main(void) {
  test_artificial_stress_none();

  return 0;
}