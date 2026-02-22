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

/* Dummy method parameters. */
static float method_artif_stress_n(void)       { return 4.f; }
static float method_artif_stress_epsilon(void) { return 0.2f; }

/* Dummy kernel. */
static void kernel_eval(const float q, float *const W) {
  if (q < 1.f) {
    *W = 1.f - q;
  } else {
    *W = 0.f;
  }
}

#define STRENGTH_ARTIFICIAL_STRESS_MON2000
#include "../../src/strength/strength_artificial_stress.h"


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

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3], pairwise_stress_tensor_j_before[3][3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = (float)(i + j + 1);
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);
  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  /* r >= h means W(r/h) = 0, so factor = 0. */
  const float r = 2.f;

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

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3], pairwise_stress_tensor_j_before[3][3];

  /* All non-positive elements. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = -(float)(i + j + 1);
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);
  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

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


/* Test that positive off-diagonal elements are also modified, unlike the
 * basis-independent method. */
static void test_off_diagonals_modified_when_positive(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3];

  /* All positive elements including off-diagonals. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] = pairwise_stress_tensor_j[i][j] = 5.f;
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);

  const float r = 0.3f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  /* All elements should have been reduced. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(pairwise_stress_tensor_i[i][j] < pairwise_stress_tensor_i_before[i][j]);
    }
  }
}


/* Test that negative elements are left unchanged while positive ones in the
 * same tensor are reduced. */
static void test_negative_elements_unchanged(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3];

  /* Mix of positive and negative elements. */
  pairwise_stress_tensor_i[0][0] =  3.f; pairwise_stress_tensor_i[0][1] = -1.f; pairwise_stress_tensor_i[0][2] =  2.f;
  pairwise_stress_tensor_i[1][0] = -2.f; pairwise_stress_tensor_i[1][1] =  4.f; pairwise_stress_tensor_i[1][2] = -0.5f;
  pairwise_stress_tensor_i[2][0] =  1.f; pairwise_stress_tensor_i[2][1] = -3.f; pairwise_stress_tensor_i[2][2] =  2.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_j[i][j] = 0.f;
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);

  const float r = 0.3f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (pairwise_stress_tensor_i_before[i][j] <= 0.f) {
        /* Negative elements must be unchanged. */
        assert(fabsf(pairwise_stress_tensor_i[i][j] - pairwise_stress_tensor_i_before[i][j]) <= tol);
      } else {
        /* Positive elements must have been reduced. */
        assert(pairwise_stress_tensor_i[i][j] < pairwise_stress_tensor_i_before[i][j]);
      }
    }
  }
}


/* Test that the reduction is proportional to each element itself, i.e. the
 * same fractional reduction is applied to all positive elements. */
static void test_fractional_reduction_uniform(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_i_before[3][3];

  /* All positive with distinct values. */
  pairwise_stress_tensor_i[0][0] = 1.f; pairwise_stress_tensor_i[0][1] = 2.f; pairwise_stress_tensor_i[0][2] = 3.f;
  pairwise_stress_tensor_i[1][0] = 4.f; pairwise_stress_tensor_i[1][1] = 5.f; pairwise_stress_tensor_i[1][2] = 6.f;
  pairwise_stress_tensor_i[2][0] = 7.f; pairwise_stress_tensor_i[2][1] = 8.f; pairwise_stress_tensor_i[2][2] = 9.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_j[i][j] = 0.f;
    }
  }

  copy_tensor(pairwise_stress_tensor_i_before, pairwise_stress_tensor_i);

  const float r = 0.3f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  /* Check all elements have the same ratio pairwise_stress_tensor_i / pairwise_stress_tensor_i_before. */
  const float ratio_ref = pairwise_stress_tensor_i[0][0] / pairwise_stress_tensor_i_before[0][0];
  const float tol = 1e-6f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const float ratio = pairwise_stress_tensor_i[i][j] / pairwise_stress_tensor_i_before[i][j];
      assert(fabsf(ratio - ratio_ref) <= tol);
    }
  }
}


/* Test that particles i and j are treated independently. */
static void test_particles_treated_independently(void)
{
  struct part pi = {0}, pj = {0};
  pi.h = 1.f;
  pj.h = 1.f;

  float pairwise_stress_tensor_i[3][3], pairwise_stress_tensor_j[3][3];
  float pairwise_stress_tensor_j_before[3][3];

  /* pairwise_stress_tensor_i all positive, pairwise_stress_tensor_j all negative. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pairwise_stress_tensor_i[i][j] =  2.f;
      pairwise_stress_tensor_j[i][j] = -2.f;
    }
  }

  copy_tensor(pairwise_stress_tensor_j_before, pairwise_stress_tensor_j);

  const float r = 0.3f;

  artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, &pi, &pj, r);

  const float tol = 1e-6f;

  /* pairwise_stress_tensor_j must be completely unchanged. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(fabsf(pairwise_stress_tensor_j[i][j] - pairwise_stress_tensor_j_before[i][j]) <= tol);
    }
  }

  /* pairwise_stress_tensor_i must have been reduced. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      assert(pairwise_stress_tensor_i[i][j] < 2.f);
    }
  }
}

// ### Add a test_separation_factor when the handling of the kernel factor is finalised

int main(void)
{
  test_no_stress_large_separation();
  test_no_stress_no_pos_eigen();
  test_off_diagonals_modified_when_positive();
  test_negative_elements_unchanged();
  test_fractional_reduction_uniform();
  test_particles_treated_independently();

  return 0;
}