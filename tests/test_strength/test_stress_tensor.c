#include <math.h>
#include <string.h>
#include <assert.h>

#define INLINE inline

/* Dummy constants. */
#define mat_phase_solid 1
#define mat_phase_fluid 0

/* Dummy sym_matrix. */
struct sym_matrix {
  union {
    struct { float elements[6]; };
    struct { float xx, yy, zz, xy, xz, yz; };
  };
};

/* Dummy particle structs. */
struct strength_data {
  struct sym_matrix deviatoric_stress_tensor;
  struct sym_matrix stress_tensor;
  struct sym_matrix dS_dt;
};

struct part {
  int mat_id;
  int phase;
  struct strength_data strength_data;
};

struct xp_strength_data {
  struct sym_matrix deviatoric_stress_tensor_full;
};

struct xpart {
  struct xp_strength_data strength_data;
};

/* Dummy hydro props */
struct hydro_props {
  float CFL_condition;
};

static float material_shear_mod(int mat_id) { return 2.f; }

static float norm_sym_matrix(const struct sym_matrix *M) {
  return sqrtf(M->xx*M->xx + M->yy*M->yy + M->zz*M->zz +
               M->xy*M->xy + M->xz*M->xz + M->yz*M->yz);
}

static void zero_sym_matrix(struct sym_matrix *M) {
  memset(M, 0, sizeof(*M));
}

static void get_matrix_from_sym_matrix(float out[3][3], const struct sym_matrix *M) {
  out[0][0]=M->xx;
  out[1][1]=M->yy;
  out[2][2]=M->zz;
  out[0][1]=out[1][0]=M->xy;
  out[0][2]=out[2][0]=M->xz;
  out[1][2]=out[2][1]=M->yz;
}

static void get_sym_matrix_from_matrix(struct sym_matrix *M, float in[3][3]) {
  M->xx=in[0][0];
  M->yy=in[1][1];
  M->zz=in[2][2];
  M->xy=in[0][1];
  M->xz=in[0][2];
  M->yz=in[1][2];
}

/* Physics hooks */
static float strength_get_damage(struct part *p) { return 0.5f; }

static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(struct sym_matrix S, float damage) {
  const float f = 1.f - damage;
  S.xx *= f; S.yy *= f; S.zz *= f;
  S.xy *= f; S.xz *= f; S.yz *= f;
  return S;
}

static void damage_compute_stress_tensor(struct sym_matrix *S,
                                         struct sym_matrix Sd,
                                         float pressure,
                                         float damage) {
  S->xx = Sd.xx - pressure;
  S->yy = Sd.yy - pressure;
  S->zz = Sd.zz - pressure;
  S->xy = Sd.xy;
  S->xz = Sd.xz;
  S->yz = Sd.yz;
}

/* Artificial stress */
static int artif_calls = 0;

static void artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
    float A[3][3], float B[3][3],
    const struct part *pi, const struct part *pj, float r)
{
  artif_calls++;
}

/* Core tested kernels */
static void strength_compute_strain_rate_tensor(float out[3][3],
                                                const float dv[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      out[i][j] = 0.5f * (dv[i][j] + dv[j][i]);
    }
  }
}

static void strength_compute_rotation_rate_tensor(float out[3][3],
                                                  const float dv[3][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      out[i][j] = 0.5f * (dv[i][j] - dv[j][i]);
    }
  }
}

static void strength_compute_rotation_term(float out[3][3],
                                           float R[3][3],
                                           float S[3][3]) {
  out[0][0] =  0.8f; out[0][1] = -0.1f; out[0][2] =  0.5f;
  out[1][0] = -0.1f; out[1][1] =  0.2f; out[1][2] = -0.3f;
  out[2][0] =  0.5f; out[2][1] = -0.3f; out[2][2] = -0.4f;
}

#include "../../src/strength/strength_stress_tensor.h"

static int within_tol(float a, float b, float tol) {
  float scale = fmaxf(fabsf(a), fabsf(b));
  if (scale < 1e-12f) scale = 1.f;
  return fabsf(a - b) <= tol * scale;
}


static void test_compute_timestep_reduction(void)
{
  struct part p = {0};
  struct hydro_props hydro_props = {0};
  hydro_props.CFL_condition = 0.1f;

  p.strength_data.deviatoric_stress_tensor.xx = 2.f;
  p.strength_data.dS_dt.xx = 4.f;

  float dt = 1.f;

  strength_compute_timestep_stress_tensor(&dt, &p, &hydro_props);

  assert(dt < 1.f);
}

static void test_compute_timestep_zero_rate(void)
{
  struct part p = {0};
  struct hydro_props hydro_props = {0};
  hydro_props.CFL_condition = 0.1f;

  float dt = 1.f;

  strength_compute_timestep_stress_tensor(&dt, &p, &hydro_props);

  assert(within_tol(dt, 1.f, 1e-6f));
}

static void test_wave_speed(void)
{
  struct part p = {.phase = mat_phase_solid};

  float wave_speed = 1.f;
  strength_compute_max_wave_speed_stress_tensor(&wave_speed, &p, 2.f, 1.f);

  assert(wave_speed > 2.f);

  p.phase = mat_phase_fluid;

  float wave_speed_2 = 1.f;
  strength_compute_max_wave_speed_stress_tensor(&wave_speed_2, &p, 2.f, 1.f);

  assert(within_tol(wave_speed_2, 1.f, 1e-6f));
}

static void test_compute_stress_tensor(void)
{
  struct part p = {0};
  p.strength_data.deviatoric_stress_tensor.xx = 2.f;

  strength_compute_stress_tensor(&p, 1.f);

  assert(p.strength_data.stress_tensor.xx < 2.f);
}

static void test_pairwise_stress_condition(void)
{
  struct part pi = {0}, pj = {0};
  float A[3][3] = {{0}}, B[3][3] = {{0}};

  pi.phase = mat_phase_solid;
  pj.phase = mat_phase_solid;

  artif_calls = 0;
  strength_set_pairwise_stress_tensors(A, B, &pi, &pj, 0.1f);
  assert(artif_calls == 1);

  pi.phase = mat_phase_solid;
  pj.phase = mat_phase_fluid;

  artif_calls = 0;
  strength_set_pairwise_stress_tensors(A, B, &pi, &pj, 0.1f);
  assert(artif_calls == 0);
}

static void test_compute_dS_dt_pure_trace(void)
{
  struct part p = {0};
  p.mat_id = 0;

  /* dv = I  =>  strain_rate = I,  trace = 3
   * deviatoric strain_rate = I - I = 0
   * 2*mu * 0 + rotation_term = rotation_term */
  float dv[3][3] = {
    {1.f,0.f,0.f},
    {0.f,1.f,0.f},
    {0.f,0.f,1.f}
  };

  stress_tensor_compute_dS_dt(&p, dv);

  const float tol = 1e-6f;

  assert(within_tol(p.strength_data.dS_dt.xx,  0.8f, tol));
  assert(within_tol(p.strength_data.dS_dt.yy,  0.2f, tol));
  assert(within_tol(p.strength_data.dS_dt.zz, -0.4f, tol));
  assert(within_tol(p.strength_data.dS_dt.xy, -0.1f, tol));
  assert(within_tol(p.strength_data.dS_dt.xz,  0.5f, tol));
  assert(within_tol(p.strength_data.dS_dt.yz, -0.3f, tol));
}

/* Test general dv gives a traceless dS/dt. */
static void test_compute_dS_dt_traceless(void)
{
  struct part p = {0};
  p.mat_id = 0;

  
  float dv[3][3] = {
    {1.f,2.f,3.f},
    {0.f,-1.f,4.f},
    {5.f,6.f,2.f}
  };

  stress_tensor_compute_dS_dt(&p, dv);

  const float mu = 2.f;
  const float tol = 1e-5f;

  /* strain_rate diagonal and off-diagonal */
  const float exx = 1.f,  eyy = -1.f, ezz = 2.f;
  const float exy = 1.f,  exz = 4.f,  eyz = 5.f;
  const float trace = exx + eyy + ezz;  /* = 2 */

  assert(within_tol(p.strength_data.dS_dt.xx, 2.f*mu*(exx - trace/3.f) +  0.8f, tol));
  assert(within_tol(p.strength_data.dS_dt.yy, 2.f*mu*(eyy - trace/3.f) +  0.2f, tol));
  assert(within_tol(p.strength_data.dS_dt.zz, 2.f*mu*(ezz - trace/3.f) + -0.4f, tol));
  assert(within_tol(p.strength_data.dS_dt.xy, 2.f*mu*exy               + -0.1f, tol));
  assert(within_tol(p.strength_data.dS_dt.xz, 2.f*mu*exz               +  0.5f, tol));
  assert(within_tol(p.strength_data.dS_dt.yz, 2.f*mu*eyz               + -0.3f, tol));
}

static void test_compute_dS_dt_pure_shear(void)
{
  struct part p = {0};
  p.mat_id = 0;

  float dv_xy = 0.5f;
  float mu = material_shear_mod(p.mat_id);

  float dv[3][3] = {
    {0.f,dv_xy,0.f},
    {0.f,0.f,0.f},
    {0.f,0.f,0.f}
  };

  stress_tensor_compute_dS_dt(&p, dv);

  const float strain_rate_xy = 0.5f * dv_xy;
  const float tol = 1e-6f;

  assert(within_tol(p.strength_data.dS_dt.xx,  0.f                   +  0.8f, tol));
  assert(within_tol(p.strength_data.dS_dt.yy,  0.f                   +  0.2f, tol));
  assert(within_tol(p.strength_data.dS_dt.zz,  0.f                   + -0.4f, tol));
  assert(within_tol(p.strength_data.dS_dt.xy,  2.f*mu*strain_rate_xy + -0.1f, tol));
  assert(within_tol(p.strength_data.dS_dt.xz,  0.f                   +  0.5f, tol));
  assert(within_tol(p.strength_data.dS_dt.yz,  0.f                   + -0.3f, tol));
}

static void test_compute_dS_dt_rotation_only(void)
{
  struct part p = {0};
  p.mat_id = 0;

  float dv[3][3] = {
    {0.f,1.f,0.f},
    {-1.f,0.f,0.f},
    {0.f,0.f,0.f}
  };

  stress_tensor_compute_dS_dt(&p, dv);

  const float tol = 1e-6f;

  assert(within_tol(p.strength_data.dS_dt.xx,  0.8f, tol));
  assert(within_tol(p.strength_data.dS_dt.yy,  0.2f, tol));
  assert(within_tol(p.strength_data.dS_dt.zz, -0.4f, tol));
  assert(within_tol(p.strength_data.dS_dt.xy, -0.1f, tol));
  assert(within_tol(p.strength_data.dS_dt.xz,  0.5f, tol));
  assert(within_tol(p.strength_data.dS_dt.yz, -0.3f, tol));
}



static void test_evolve_solid(void)
{
  struct part p = {0};
  p.strength_data.dS_dt.xx = 1.f;
  p.strength_data.dS_dt.yy = 2.f;
  p.strength_data.dS_dt.zz = 3.f;
  p.strength_data.dS_dt.xy = 4.f;
  p.strength_data.dS_dt.xz = 5.f;
  p.strength_data.dS_dt.yz = 6.f;
  struct sym_matrix S = {0};
  S.xx = 10.f;
  S.yy = 20.f;
  S.zz = 30.f;
  S.xy = 40.f;
  S.xz = 50.f;
  S.yz = 60.f;


  stress_tensor_evolve_deviatoric_stress_tensor(&S, &p, mat_phase_solid, 0.1f);

  assert(within_tol(S.xx, 10.1f, 1e-6f));
  assert(within_tol(S.yy, 20.2f, 1e-6f));
  assert(within_tol(S.zz, 30.3f, 1e-6f));
  assert(within_tol(S.xy, 40.4f, 1e-6f));
  assert(within_tol(S.xz, 50.5f, 1e-6f));
  assert(within_tol(S.yz, 60.6f, 1e-6f));
}

static void test_evolve_fluid(void)
{
  struct part p = {0};
  p.strength_data.dS_dt.xx = 1.f;
  p.strength_data.dS_dt.yy = 2.f;
  p.strength_data.dS_dt.zz = 3.f;
  p.strength_data.dS_dt.xy = 4.f;
  p.strength_data.dS_dt.xz = 5.f;
  p.strength_data.dS_dt.yz = 6.f;
  struct sym_matrix S = {0};
  S.xx = 10.f;
  S.yy = 20.f;
  S.zz = 30.f;
  S.xy = 40.f;
  S.xz = 50.f;
  S.yz = 60.f;

  stress_tensor_evolve_deviatoric_stress_tensor(&S, &p, mat_phase_fluid, 0.1f);

  assert(within_tol(S.xx, 0.f, 1e-6f));
  assert(within_tol(S.yy, 0.f, 1e-6f));
  assert(within_tol(S.zz, 0.f, 1e-6f));
  assert(within_tol(S.xy, 0.f, 1e-6f));
  assert(within_tol(S.xz, 0.f, 1e-6f));
  assert(within_tol(S.yz, 0.f, 1e-6f));
}

static void test_first_init(void)
{
  struct part p = {0};
  struct xpart xp = {0};

  strength_first_init_part_stress_tensor(&p, &xp);

  assert(within_tol(p.strength_data.deviatoric_stress_tensor.xx, 0.f, 1e-6f));
  assert(within_tol(p.strength_data.deviatoric_stress_tensor.yy, 0.f, 1e-6f));
  assert(within_tol(p.strength_data.deviatoric_stress_tensor.zz, 0.f, 1e-6f));
  assert(within_tol(p.strength_data.deviatoric_stress_tensor.xy, 0.f, 1e-6f));
  assert(within_tol(p.strength_data.deviatoric_stress_tensor.xz, 0.f, 1e-6f));
  assert(within_tol(p.strength_data.deviatoric_stress_tensor.yz, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.xx, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.yy, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.zz, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.xy, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.xz, 0.f, 1e-6f));
  assert(within_tol(xp.strength_data.deviatoric_stress_tensor_full.yz, 0.f, 1e-6f));
}

int main(void)
{
  test_compute_timestep_reduction();
  test_compute_timestep_zero_rate();
  test_wave_speed();
  test_compute_stress_tensor();
  test_pairwise_stress_condition();
  test_compute_dS_dt_pure_trace();
  test_compute_dS_dt_traceless();
  test_compute_dS_dt_pure_shear();
  test_compute_dS_dt_rotation_only();
  test_evolve_solid();
  test_evolve_fluid();
  test_first_init();

  return 0;
}