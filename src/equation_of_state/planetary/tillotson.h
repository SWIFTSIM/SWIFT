/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
 *               2018   Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_TILLOTSON_EQUATION_OF_STATE_H
#define SWIFT_TILLOTSON_EQUATION_OF_STATE_H

/**
 * @file equation_of_state/planetary/tillotson.h
 *
 * Contains the Tillotson EOS functions for
 * equation_of_state/planetary/equation_of_state.h
 *
 */

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "adiabatic_index.h"
#include "common_io.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "sesame.h"
#include "units.h"

// Tillotson parameters
struct Til_params {
  float rho_0, a, b, A, B, u_0, u_iv, u_cv, alpha, beta;
  float eta_min, eta_zero, P_min, C_V;
  float *A1_u_cold;
  float rho_cold_min, rho_cold_max;
  float rho_min, rho_max;
  enum eos_planetary_material_id mat_id;
};

// Parameter values for each material (SI units)
INLINE static void set_Til_iron(struct Til_params *mat,
                                enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 7800.0f;
  mat->a = 0.5f;
  mat->b = 1.5f;
  mat->A = 1.28e11f;
  mat->B = 1.05e11f;
  mat->u_0 = 9.5e6f;
  mat->u_iv = 2.4e6f;
  mat->u_cv = 8.67e6f;
  mat->alpha = 5.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.0f;
  mat->eta_zero = 0.0f;
  mat->P_min = -FLT_MAX;//0.01f;
  mat->C_V = 449.0f;
  mat->rho_cold_min = 100.0f;
  mat->rho_cold_max = 1.0e5f;
  mat->rho_min = 1.0f;
  mat->rho_max = 1.0e5f;
}
INLINE static void set_Til_granite(struct Til_params *mat,
                                   enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 2680.0f;
  mat->a = 0.5f;
  mat->b = 1.3f;
  mat->A = 1.8e10f;
  mat->B = 1.8e10f;
  mat->u_0 = 1.6e7f;
  mat->u_iv = 3.5e6f;
  mat->u_cv = 1.8e7f;
  mat->alpha = 5.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.0f;
  mat->eta_zero = 0.0f;
  mat->P_min = -FLT_MAX;//0.01f;
  mat->C_V = 790.0f;
  mat->rho_cold_min = 100.0f;
  mat->rho_cold_max = 1.0e5f;
  mat->rho_min = 1.0f;
  mat->rho_max = 1.0e5f;
}
INLINE static void set_Til_basalt(struct Til_params *mat,
                                  enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 2700.0f;
  mat->a = 0.5f;
  mat->b = 1.5f;
  mat->A = 2.67e10f;
  mat->B = 2.67e10f;
  mat->u_0 = 4.87e8f;
  mat->u_iv = 4.72e6f;
  mat->u_cv = 1.82e7f;
  mat->alpha = 5.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.0f;
  mat->eta_zero = 0.0f;
  mat->P_min = -FLT_MAX;//0.01f;
  mat->C_V = 790.0f;
  mat->rho_cold_min = 100.0f;
  mat->rho_cold_max = 1.0e5f;
  mat->rho_min = 1.0f;
  mat->rho_max = 1.0e5f;
}
INLINE static void set_Til_water(struct Til_params *mat,
                                 enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 998.0f;
  mat->a = 0.7f;
  mat->b = 0.15f;
  mat->A = 2.18e9f;
  mat->B = 1.325e10f;
  mat->u_0 = 7.0e6f;
  mat->u_iv = 4.19e5f;
  mat->u_cv = 2.69e6f;
  mat->alpha = 10.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.925f;
  mat->eta_zero = 0.875f;
  mat->P_min = -FLT_MAX;//0.01f;
  mat->C_V = 4186.0f;
  mat->rho_cold_min = 100.0f;
  mat->rho_cold_max = 1.0e5f;
  mat->rho_min = 1.0f;
  mat->rho_max = 1.0e5f;
}
INLINE static void set_Til_ice(struct Til_params *mat,
                               enum eos_planetary_material_id mat_id) {
  mat->mat_id = mat_id;
  mat->rho_0 = 1293.0f;
  mat->a = 0.3f;
  mat->b = 0.1f;
  mat->A = 1.07e10f;
  mat->B = 6.5e10f;
  mat->u_0 = 1.0e7f;
  mat->u_iv = 7.73e5f;
  mat->u_cv = 3.04e6f;
  mat->alpha = 10.0f;
  mat->beta = 5.0f;
  mat->eta_min = 0.925f;
  mat->eta_zero = 0.875f;
  mat->P_min = -FLT_MAX;//0.0f;
  mat->C_V = 2093.0f;
  mat->rho_cold_min = 100.0f;
  mat->rho_cold_max = 1.0e5f;
  mat->rho_min = 1.0f;
  mat->rho_max = 1.0e5f;
}

/*
    Read the parameters from a file.

    File contents
    -------------
    # header (5 lines)
    rho_0 (kg/m3)  a (-)  b (-)  A (Pa)  B (Pa)
    u_0 (J)  u_iv (J)  u_cv (J)  alpha (-)  beta (-)
    eta_min (-)  eta_zero (-)  P_min (Pa)  C_V (J kg^-1 K^-1)
    rho_cold_min (kg/m3)  rho_cold_max (kg/m3)  rho_min (kg/m3)  rho_max (kg/m3)
*/
INLINE static void set_Til_custom(struct Til_params *mat,
                                  enum eos_planetary_material_id mat_id,
                                  char *param_file) {
  mat->mat_id = mat_id;

  // Load table contents from file
  FILE *f = fopen(param_file, "r");
  if (f == NULL)
    error("Failed to open the Tillotson EoS file '%s'", param_file);

  // Skip header lines
  skip_lines(f, 5);

  // Read parameters (SI)
  int c;
  c = fscanf(f, "%f %f %f %f %f", &mat->rho_0, &mat->a, &mat->b, &mat->A,
             &mat->B);
  if (c != 5) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f %f", &mat->u_0, &mat->u_iv, &mat->u_cv,
             &mat->alpha, &mat->beta);
  if (c != 5) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f", &mat->eta_min, &mat->eta_zero, &mat->P_min,
             &mat->C_V);
  if (c != 4) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f", &mat->rho_cold_min, &mat->rho_cold_max,
             &mat->rho_min, &mat->rho_max);
  if (c != 4) error("Failed to read the Tillotson EoS file %s", param_file);
}

// Convert to internal units
INLINE static void convert_units_Til(struct Til_params *mat,
                                     const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  int N = 10000;

  // SI to cgs
  mat->rho_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mat->A *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  mat->B *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  mat->u_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_iv *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_cv *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->P_min *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

  for (int i = 0; i < N; i++) {
    mat->A1_u_cold[i] *=
        units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }

  mat->C_V *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin

  mat->rho_cold_min *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mat->rho_cold_max *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mat->rho_min *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mat->rho_max *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);

  // cgs to internal
  mat->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->A /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->B /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  mat->u_0 /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_iv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->u_cv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  mat->P_min /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

  for (int i = 0; i < N; i++) {
    mat->A1_u_cold[i] /=
        units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }

  mat->C_V /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin

  mat->rho_cold_min /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->rho_cold_max /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->rho_min /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->rho_max /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
}

// gas_internal_energy_from_entropy
INLINE static float Til_internal_energy_from_entropy(
    float density, float entropy, const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_pressure_from_entropy
INLINE static float Til_pressure_from_entropy(float density, float entropy,
                                              const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float Til_entropy_from_pressure(float density, float pressure,
                                              const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float Til_soundspeed_from_entropy(float density, float entropy,
                                                const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float Til_entropy_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float Til_pressure_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  const float eta = density / mat->rho_0;
  const float eta_sq = eta * eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (mat->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w;
  float P_c, P_e, P;

  // Condensed or cold
  P_c =
      (mat->a + mat->b * w_inv) * density * u + mat->A * mu + mat->B * mu * mu;
  if (eta < mat->eta_zero) {
    P_c = 0.f;
  } else if (eta < mat->eta_min) {
    P_c *= (eta - mat->eta_zero) / (mat->eta_min - mat->eta_zero);
  }
  // Expanded and hot
  P_e = mat->a * density * u +
        (mat->b * density * u * w_inv + mat->A * mu * expf(-mat->beta * nu)) *
            expf(-mat->alpha * nu * nu);

  // Condensed or cold state
  if ((1.f < eta) || (u < mat->u_iv)) {
    P = P_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (mat->u_cv < u)) {
    P = P_e;
  }
  // Hybrid state
  else {
    P = ((u - mat->u_iv) * P_e + (mat->u_cv - u) * P_c) /
        (mat->u_cv - mat->u_iv);
  }

  // Minimum pressure
  if (P < mat->P_min) {
    P = mat->P_min;
  }
    
  // For aluminium. This should be P_min: P_min = -Y0
        // This should be made P_min in mat properties
    // not sure whether to fix or have weakening
    /*
    float Y_0 = 200e6;
  if (P < -Y_0) {
    P = -Y_0;
  }
    */
    
    // With rho weakening
    float rho0 = 2700.f;
  float rho_weak = 0.85f * rho0;
    float Y_min = -200e6;
    if (density < rho_weak){
        Y_min *= powf(density / rho_weak, 4.f);
    } 
    
    if (P < Y_min) {
    P = Y_min;
  }
  

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float Til_internal_energy_from_pressure(
    float density, float P, const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float Til_soundspeed_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  const float rho_0_inv = 1.f / mat->rho_0;
  const float eta = density * rho_0_inv;
  const float rho_inv = 1.f / density;
  const float eta_sq = eta * eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (mat->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w;
  const float w_inv_sq = w_inv * w_inv;
  const float exp_beta = expf(-mat->beta * nu);
  const float exp_alpha = expf(-mat->alpha * nu * nu);
  float P_c, P_e, c_sq_c, c_sq_e, c_sq;

  // Condensed or cold
  P_c =
      (mat->a + mat->b * w_inv) * density * u + mat->A * mu + mat->B * mu * mu;
  if (eta < mat->eta_zero) {
    P_c = 0.f;
  } else if (eta < mat->eta_min) {
    P_c *= (eta - mat->eta_zero) / (mat->eta_min - mat->eta_zero);
  }
  c_sq_c = P_c * rho_inv * (1.f + mat->a + mat->b * w_inv) +
           mat->b * (w - 1.f) * w_inv_sq * (2.f * u - P_c * rho_inv) +
           rho_inv * (mat->A + mat->B * (eta_sq - 1.f));

  // Expanded and hot
  P_e = mat->a * density * u +
        (mat->b * density * u * w_inv + mat->A * mu * exp_beta) * exp_alpha;

  c_sq_e =
      P_e * rho_inv * (1.f + mat->a + mat->b * w_inv * exp_alpha) +
      (mat->b * density * u * w_inv_sq / eta_sq *
           (rho_inv / mat->u_0 * (2.f * u - P_e * rho_inv) +
            2.f * mat->alpha * nu * w * rho_0_inv) +
       mat->A * rho_0_inv *
           (1.f + mu / eta_sq * (mat->beta + 2.f * mat->alpha * nu - eta)) *
           exp_beta) *
          exp_alpha;

  // Condensed or cold state
  if ((1.f < eta) || (u < mat->u_iv)) {
    c_sq = c_sq_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (mat->u_cv < u)) {
    c_sq = c_sq_e;
  }
  // Hybrid state
  else {
    c_sq = ((u - mat->u_iv) * c_sq_e + (mat->u_cv - u) * c_sq_c) /
           (mat->u_cv - mat->u_iv);
  }

  c_sq = fmaxf(c_sq, mat->A * rho_0_inv);

  return sqrtf(c_sq);
}

// gas_soundspeed_from_pressure
INLINE static float Til_soundspeed_from_pressure(float density, float P,
                                                 const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// Compute u cold
INLINE static float compute_u_cold(float density, struct Til_params *mat,
                                   enum eos_planetary_material_id mat_id) {
  float rho_0, x, u_cold, drho;
  int N = 10000;

  rho_0 = mat->rho_0;
  drho = (density - rho_0) / N;
  x = rho_0;
  u_cold = 1e-9;

  for (int i = 0; i < N; i++) {
    x += drho;
    u_cold +=
        Til_pressure_from_internal_energy(x, u_cold, mat) * drho / (x * x);
  }

  return u_cold;

}

// Compute A1_u_cold
INLINE static void set_Til_u_cold(struct Til_params *mat,
                                  enum eos_planetary_material_id mat_id) {

  int N = 10000;
  float rho_min = 100.f;
  float rho_max = 100000.f;
  float rho, drho;

  // Allocate table memory
  mat->A1_u_cold = (float *)malloc(N * sizeof(float));

  rho = rho_min;
  drho = (rho_max - rho_min) / (N - 1);

  for (int i = 0; i < N; i++) {
    mat->A1_u_cold[i] = compute_u_cold(rho, mat, mat_id);
    rho += drho;
  }
}

// Compute u cold fast from precomputed values
INLINE static float compute_fast_u_cold(float density,
                                        const struct Til_params *mat) {

  int N = 10000;
  float rho_min = mat->rho_cold_min;
  float rho_max = mat->rho_cold_max;
  float drho, u_cold;
  int a, b;

  drho = (rho_max - rho_min) / (N - 1);

  a = (int)((density - rho_min) / drho);
  b = a + 1;

  if (a >= 0 && a < (N - 1)) {
    u_cold = mat->A1_u_cold[a];
    u_cold += ((mat->A1_u_cold[b] - mat->A1_u_cold[a]) / drho) *
              (density - rho_min - a * drho);
  } else if (density < rho_min) {
    u_cold = mat->A1_u_cold[0];
  } else {
    u_cold = mat->A1_u_cold[N - 1];
    u_cold += ((mat->A1_u_cold[N - 1] - mat->A1_u_cold[N - 2]) / drho) *
              (density - rho_max);
  }
  return u_cold;
}

// gas_temperature_from_internal_energy
INLINE static float Til_temperature_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

  float u_cold, T;

  u_cold = compute_fast_u_cold(density, mat);

  T = (u - u_cold) / (mat->C_V);
  if (T < 0.f) {
    T = 0.f;
  }

  return T;
}

// gas_pressure_from_density_and_temperature
INLINE static float Til_pressure_from_temperature(
    float density, float T, const struct Til_params *mat) {

  float u, P;

  u = compute_fast_u_cold(density, mat) + mat->C_V * T;
  P = Til_pressure_from_internal_energy(density, u, mat);

  return P;
}

// gas_density_from_pressure_and_temperature
INLINE static float Til_density_from_pressure_and_temperature(
    float P, float T, const struct Til_params *mat) {

  float rho_min = mat->rho_min;
  float rho_max = mat->rho_max;
  float rho_mid = (rho_min + rho_max) / 2.f;
  float P_min, P_mid, P_max;
  float P_des;
  float tolerance = 0.001 * rho_min;
  int counter = 0;
  int max_counter = 200;
  float f0, f2;

  // Check for P == 0 or T == 0
  if (P <= mat->P_min) {
    P_des = mat->P_min;
  } else {
    P_des = P;
  }

  P_min = Til_pressure_from_temperature(rho_min, T, mat);
  P_mid = Til_pressure_from_temperature(rho_mid, T, mat);
  P_max = Til_pressure_from_temperature(rho_max, T, mat);

  // quick fix?
  if (P_des < P_min) {
    P_des = P_min;
  }

  if (P_des >= P_min && P_des <= P_max) {
    while ((rho_max - rho_min) > tolerance && counter < max_counter) {

      P_min = Til_pressure_from_temperature(rho_min, T, mat);
      P_mid = Til_pressure_from_temperature(rho_mid, T, mat);
      P_max = Til_pressure_from_temperature(rho_max, T, mat);

      f0 = P_des - P_min;
      f2 = P_des - P_mid;

      if ((f0 * f2) > 0) {
        rho_min = rho_mid;
      } else {
        rho_max = rho_mid;
      }

      rho_mid = (rho_min + rho_max) / 2.f;
      counter += 1;
    }
  } else {
    error("Error in Til_density_from_pressure_and_temperature");
    return 0.f;
  }
  return rho_mid;
}

// gas_density_from_pressure_and_internal_energy
INLINE static float Til_density_from_pressure_and_internal_energy(
    float P, float u, float rho_ref, float rho_sph,
    const struct Til_params *mat) {

  if (P <= mat->P_min || u == 0) {
    return rho_sph;
  }

  // These are needed in root finding iteration
  float eta_iter;
  float eta_iter_sq;
  float mu_iter;
  float nu_iter;
  float w_iter;
  float w_iter_inv;
  float exp1;
  float exp2;

  // Derivatives
  float dw_inv_drho_iter;
  float dmu_drho_iter;
  float dmu_sq_drho_iter;
  float dexp1_drho_iter;
  float dexp2_drho_iter;

  // We start search on the same curve as rho_ref, since this is most likely
  // curve to find rho on
  float eta_ref = rho_ref / mat->rho_0;

  /*
  There are 5 possible curves:
    1: cold_min
    2: cold
    3: hybrid_min
    4: hybrid
    5: hot

  These curves cover different eta ranges within three different regions of u:
  u REGION 1 (u < u_iv):
      eta < eta_min:           cold_min
      eta_min < eta:           cold

 u REGION 2 (u_iv < u < u_cv):
      eta < eta_min:           hybrid_min
      eta_min < eta < 1:       hybrid
      1 < eta:                 cold

  u REGION 3 (u_cv < u):
      eta < 1:                 hot
      1 < eta:                 cold

  NOTE: for a lot of EoS, eta_min = 0, so search this region last if given the
 option to save time for most EoS
  */

  // Based on our u region, what possible curves can rho be on? Ordered based on
  // order we search for roots . Numbers correspond to curves in order given
  // above. 0 is a dummy which breaks the loop.
  // int possible_curves[3];
  int possible_curves[3];
  // u REGION 1
  if (u <= mat->u_iv) {
    if (eta_ref <= mat->eta_min) {
      possible_curves[0] = 1;
      possible_curves[1] = 2;
      possible_curves[2] = 0;
    } else {
      possible_curves[0] = 2;
      possible_curves[1] = 1;
      possible_curves[2] = 0;
    }
    // u REGION 2
  } else if (u <= mat->u_cv) {
    if (eta_ref <= mat->eta_min) {
      possible_curves[0] = 3;
      possible_curves[1] = 4;
      possible_curves[2] = 2;
    } else if (eta_ref <= 1) {
      possible_curves[0] = 4;
      possible_curves[1] = 2;
      possible_curves[2] = 3;
    } else {
      possible_curves[0] = 2;
      possible_curves[1] = 4;
      possible_curves[2] = 3;
    }

    // u REGION 1
  } else {
    if (eta_ref <= 1) {
      possible_curves[0] = 5;
      possible_curves[1] = 2;
      possible_curves[2] = 0;
    } else {
      possible_curves[0] = 2;
      possible_curves[1] = 5;
      possible_curves[2] = 0;
    }
  }
  // Newton-Raphson
  int max_iter = 10;
  float tol = 1e-5;

  // loops over possible curves
  for (int i = 0; i < 3; i++) {

    int curve = possible_curves[i];

    // if there are only two possible curves, break when we get to three and
    // haven't found a root
    if (curve == 0) {
      break;
    }

    // Start iteration at reference value
    float rho_iter = rho_ref;

    // Constrain our initial guess to be on the curve we're currently looking at
    // in the first loop, this is already satisfied.
    if (i > 0) {

      // u REGION 1
      if (u <= mat->u_iv) {
        if (curve == 1) {
          if (rho_iter > mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
        // u REGION 2
      } else if (u <= mat->u_cv) {
        if (curve == 3) {
          if (rho_iter > mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else if (curve == 4) {
          if (rho_iter < mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
          if (rho_iter > mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }

        // u REGION 3
      } else {
        if (curve == 5) {
          if (rho_iter > mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
      }
    }

    // Set this to an arbitrary number so we definitely dont't think we converge
    // straigt away
    float last_rho_iter = -1e5;

    // Now iterate
    for (int j = 0; j < max_iter; j++) {

      eta_iter = rho_iter / mat->rho_0;
      eta_iter_sq = eta_iter * eta_iter;
      mu_iter = eta_iter - 1.0f;
      nu_iter = 1.0f / eta_iter - 1.0f;
      w_iter = u / (mat->u_0 * eta_iter_sq) + 1.0f;
      w_iter_inv = 1.0f / w_iter;
      exp1 = expf(-mat->beta * nu_iter);
      exp2 = expf(-mat->alpha * nu_iter * nu_iter);

      // Derivatives
      dw_inv_drho_iter =
          (2.f * mat->u_0 * u * eta_iter / mat->rho_0) /
          ((u + mat->u_0 * eta_iter_sq) * (u + mat->u_0 * eta_iter_sq));
      dmu_drho_iter = 1.f / mat->rho_0;
      dmu_sq_drho_iter =
          2.f * rho_iter / (mat->rho_0 * mat->rho_0) - 2.f / mat->rho_0;
      dexp1_drho_iter = mat->beta * mat->rho_0 * exp1 / (rho_iter * rho_iter);
      dexp2_drho_iter = 2.f * mat->alpha * mat->rho_0 *
                        (mat->rho_0 - rho_iter) * exp2 /
                        (rho_iter * rho_iter * rho_iter);

      // Use P_fraction to determine whether we've converged on a root
      float P_fraction;

      // Newton-Raphson
      float P_c_iter, dP_c_drho_iter, P_h_iter, dP_h_drho_iter;

      // if "cold" or "hybrid"
      if (curve == 2 || curve == 4) {
        P_c_iter = (mat->a + mat->b * w_iter_inv) * rho_iter * u +
                   mat->A * mu_iter + mat->B * mu_iter * mu_iter - P;
        dP_c_drho_iter = (mat->a + mat->b * w_iter_inv) * u +
                         mat->b * u * rho_iter * dw_inv_drho_iter +
                         mat->A * dmu_drho_iter + mat->B * dmu_sq_drho_iter;
        P_fraction = P_c_iter / P;

        // if curve is cold then we've got everything we need
        if (curve == 2) {
          rho_iter -= P_c_iter / dP_c_drho_iter;
          // Don't use these:
          P_h_iter = 0.f;
          dP_h_drho_iter = 0.f;
        }
        // if "cold_min" or "hybrid_min"
        // Can only have one version of the cold curve, therefore either use the
        // min version or the normal version hence "else if"
      } else if (curve == 1 || curve == 3) {
        P_c_iter = ((mat->a + mat->b * w_iter_inv) * rho_iter * u +
                    mat->A * mu_iter + mat->B * mu_iter * mu_iter) *
                       (eta_iter - mat->eta_zero) /
                       (mat->eta_min - mat->eta_zero) -
                   P;
        dP_c_drho_iter =
            ((mat->a + mat->b * w_iter_inv) * u +
             mat->b * u * rho_iter * dw_inv_drho_iter + mat->A * dmu_drho_iter +
             mat->B * dmu_sq_drho_iter) *
                (eta_iter - mat->eta_zero) / (mat->eta_min - mat->eta_zero) +
            ((mat->a + mat->b * w_iter_inv) * rho_iter * u + mat->A * mu_iter +
             mat->B * mu_iter * mu_iter) *
                (1 / (mat->rho_0 * (mat->eta_min - mat->eta_zero)));
        P_fraction = P_c_iter / P;

        // if curve is cold_min then we've got everything we need
        if (curve == 1) {
          rho_iter -= P_c_iter / dP_c_drho_iter;
          // Don't use these:
          P_c_iter = 0.f;
          dP_c_drho_iter = 0.f;
        }
      }

      // if "hybrid_min" or "hybrid" or "hot"
      if (curve == 3 || curve == 4 || curve == 5) {
        P_h_iter =
            mat->a * rho_iter * u +
            (mat->b * rho_iter * u * w_iter_inv + mat->A * mu_iter * exp1) *
                exp2 -
            P;
        dP_h_drho_iter =
            mat->a * u +
            (mat->b * u * w_iter_inv +
             mat->b * u * rho_iter * dw_inv_drho_iter +
             mat->A * mu_iter * dexp1_drho_iter +
             mat->A * exp1 * dmu_drho_iter) *
                exp2 +
            (mat->b * rho_iter * u * w_iter_inv + mat->A * mu_iter * exp1) *
                dexp2_drho_iter;
        P_fraction = P_h_iter / P;

        // if curve is hot then we've got everything we need
        if (curve == 5) {
          rho_iter -= P_h_iter / dP_h_drho_iter;
        }
      }

      // If we are on a hybrid or hybrid_min curve, we combie hot and cold
      // curves
      if (curve == 3 || curve == 4) {
        float P_hybrid_iter =
            ((u - mat->u_iv) * P_h_iter + (mat->u_cv - u) * P_c_iter) /
            (mat->u_cv - mat->u_iv);
        float dP_hybrid_drho_iter = ((u - mat->u_iv) * dP_h_drho_iter +
                                     (mat->u_cv - u) * dP_c_drho_iter) /
                                    (mat->u_cv - mat->u_iv);
        rho_iter -= P_hybrid_iter / dP_hybrid_drho_iter;
        P_fraction = P_hybrid_iter / P;
      }

      // Now we have to constrain the new rho_iter to the curve we're on
      // u REGION 1
      if (u <= mat->u_iv) {
        if (curve == 1) {
          if (rho_iter > mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
        // u REGION 2
      } else if (u <= mat->u_cv) {
        if (curve == 3) {
          if (rho_iter > mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
        } else if (curve == 4) {
          if (rho_iter < mat->eta_min * mat->rho_0) {
            rho_iter = mat->eta_min * mat->rho_0;
          }
          if (rho_iter > mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }

        // u REGION 3
      } else {
        if (curve == 5) {
          if (rho_iter > mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < mat->rho_0) {
            rho_iter = mat->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
      }

      // Either we've converged ...
      if (fabs(P_fraction) < tol) {
        return rho_iter;

        // ... or we're stuck at the boundary ...
      } else if (rho_iter == last_rho_iter) {
        break;
      }

      // ... or we loop again
      last_rho_iter = rho_iter;
    }
  }
  return rho_sph;
}

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */
