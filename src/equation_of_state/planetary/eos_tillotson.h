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
#include "eos_setup.h"
#include "equation_of_state.h"
#include "inline.h"
#include "physical_constants.h"
#include "units.h"
#include "utilities.h"

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
INLINE static void set_Til_iron(struct Til_params *Til,
                                enum eos_planetary_material_id mat_id) {
  Til->mat_id = mat_id;
  Til->rho_0 = 7800.0f;
  Til->a = 0.5f;
  Til->b = 1.5f;
  Til->A = 1.28e11f;
  Til->B = 1.05e11f;
  Til->u_0 = 9.5e6f;
  Til->u_iv = 2.4e6f;
  Til->u_cv = 8.67e6f;
  Til->alpha = 5.0f;
  Til->beta = 5.0f;
  Til->eta_min = 0.925f;
  Til->eta_zero = 0.875f;
  Til->C_V = 449.0f;
  Til->rho_cold_min = 100.0f;
  Til->rho_cold_max = 1.0e5f;
  Til->rho_min = 1.0f;
  Til->rho_max = 1.0e5f;
#ifdef MATERIAL_STRENGTH
  Til->P_min = -FLT_MIN;
#else
  Til->P_min = 0.01f;
#endif /* MATERIAL_STRENGTH */
}
INLINE static void set_Til_granite(struct Til_params *Til,
                                   enum eos_planetary_material_id mat_id) {
  Til->mat_id = mat_id;
  Til->rho_0 = 2680.0f;
  Til->a = 0.5f;
  Til->b = 1.3f;
  Til->A = 1.8e10f;
  Til->B = 1.8e10f;
  Til->u_0 = 1.6e7f;
  Til->u_iv = 3.5e6f;
  Til->u_cv = 1.8e7f;
  Til->alpha = 5.0f;
  Til->beta = 5.0f;
  Til->eta_min = 0.925f;
  Til->eta_zero = 0.875f;
  Til->C_V = 790.0f;
  Til->rho_cold_min = 100.0f;
  Til->rho_cold_max = 1.0e5f;
  Til->rho_min = 1.0f;
  Til->rho_max = 1.0e5f;
#ifdef MATERIAL_STRENGTH
  Til->P_min = -FLT_MIN;
#else
  Til->P_min = 0.01f;
#endif /* MATERIAL_STRENGTH */
}
INLINE static void set_Til_basalt(struct Til_params *Til,
                                  enum eos_planetary_material_id mat_id) {
  Til->mat_id = mat_id;
  Til->rho_0 = 2700.0f;
  Til->a = 0.5f;
  Til->b = 1.5f;
  Til->A = 2.67e10f;
  Til->B = 2.67e10f;
  Til->u_0 = 4.87e8f;
  Til->u_iv = 4.72e6f;
  Til->u_cv = 1.82e7f;
  Til->alpha = 5.0f;
  Til->beta = 5.0f;
  Til->eta_min = 0.925f;
  Til->eta_zero = 0.875f;
  Til->C_V = 790.0f;
  Til->rho_cold_min = 100.0f;
  Til->rho_cold_max = 1.0e5f;
  Til->rho_min = 1.0f;
  Til->rho_max = 1.0e5f;
#ifdef MATERIAL_STRENGTH
  Til->P_min = -FLT_MIN;
#else
  Til->P_min = 0.01f;
#endif /* MATERIAL_STRENGTH */
}
INLINE static void set_Til_water(struct Til_params *Til,
                                 enum eos_planetary_material_id mat_id) {
  Til->mat_id = mat_id;
  Til->rho_0 = 998.0f;
  Til->a = 0.7f;
  Til->b = 0.15f;
  Til->A = 2.18e9f;
  Til->B = 1.325e10f;
  Til->u_0 = 7.0e6f;
  Til->u_iv = 4.19e5f;
  Til->u_cv = 2.69e6f;
  Til->alpha = 10.0f;
  Til->beta = 5.0f;
  Til->eta_min = 0.925f;
  Til->eta_zero = 0.875f;
  Til->C_V = 4186.0f;
  Til->rho_cold_min = 100.0f;
  Til->rho_cold_max = 1.0e5f;
  Til->rho_min = 1.0f;
  Til->rho_max = 1.0e5f;
#ifdef MATERIAL_STRENGTH
  Til->P_min = -FLT_MIN;
#else
  Til->P_min = 0.01f;
#endif /* MATERIAL_STRENGTH */
}
INLINE static void set_Til_ice(struct Til_params *Til,
                               enum eos_planetary_material_id mat_id) {
  Til->mat_id = mat_id;
  Til->rho_0 = 1293.0f;
  Til->a = 0.3f;
  Til->b = 0.1f;
  Til->A = 1.07e10f;
  Til->B = 6.5e10f;
  Til->u_0 = 1.0e7f;
  Til->u_iv = 7.73e5f;
  Til->u_cv = 3.04e6f;
  Til->alpha = 10.0f;
  Til->beta = 5.0f;
  Til->eta_min = 0.925f;
  Til->eta_zero = 0.875f;
  Til->C_V = 2093.0f;
  Til->rho_cold_min = 100.0f;
  Til->rho_cold_max = 1.0e5f;
  Til->rho_min = 1.0f;
  Til->rho_max = 1.0e5f;
#ifdef MATERIAL_STRENGTH
  Til->P_min = -FLT_MIN;
#else
  Til->P_min = 0.0f;
#endif /* MATERIAL_STRENGTH */
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
INLINE static void set_Til_custom(struct Til_params *Til,
                                  enum eos_planetary_material_id mat_id,
                                  char *param_file) {
  Til->mat_id = mat_id;

  // Load table contents from file
  FILE *f = fopen(param_file, "r");
  if (f == NULL)
    error("Failed to open the Tillotson EoS file '%s'", param_file);

  // Skip header lines
  skip_lines(f, 5);

  // Read parameters (SI)
  int c;
  c = fscanf(f, "%f %f %f %f %f", &Til->rho_0, &Til->a, &Til->b, &Til->A,
             &Til->B);
  if (c != 5) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f %f", &Til->u_0, &Til->u_iv, &Til->u_cv,
             &Til->alpha, &Til->beta);
  if (c != 5) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f", &Til->eta_min, &Til->eta_zero, &Til->P_min,
             &Til->C_V);
  if (c != 4) error("Failed to read the Tillotson EoS file %s", param_file);

  c = fscanf(f, "%f %f %f %f", &Til->rho_cold_min, &Til->rho_cold_max,
             &Til->rho_min, &Til->rho_max);
  if (c != 4) error("Failed to read the Tillotson EoS file %s", param_file);
}

// Convert to internal units
INLINE static void convert_units_Til(struct Til_params *Til,
                                     const struct unit_system *us) {

  struct unit_system si;
  units_init_si(&si);

  int N = 10000;

  // SI to cgs
  Til->rho_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  Til->A *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  Til->B *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);
  Til->u_0 *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->u_iv *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->u_cv *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->P_min *= units_cgs_conversion_factor(&si, UNIT_CONV_PRESSURE);

  for (int i = 0; i < N; i++) {
    Til->A1_u_cold[i] *=
        units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }

  Til->C_V *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin

  Til->rho_cold_min *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  Til->rho_cold_max *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  Til->rho_min *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  Til->rho_max *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);

  // cgs to internal
  Til->rho_0 /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  Til->A /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  Til->B /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);
  Til->u_0 /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->u_iv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->u_cv /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  Til->P_min /= units_cgs_conversion_factor(us, UNIT_CONV_PRESSURE);

  for (int i = 0; i < N; i++) {
    Til->A1_u_cold[i] /=
        units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }

  Til->C_V /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin

  Til->rho_cold_min /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  Til->rho_cold_max /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  Til->rho_min /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  Til->rho_max /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
}

// gas_internal_energy_from_entropy
INLINE static float Til_internal_energy_from_entropy(
    float density, float entropy, const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_pressure_from_entropy
INLINE static float Til_pressure_from_entropy(float density, float entropy,
                                              const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_pressure
INLINE static float Til_entropy_from_pressure(float density, float pressure,
                                              const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_entropy
INLINE static float Til_soundspeed_from_entropy(float density, float entropy,
                                                const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_entropy_from_internal_energy
INLINE static float Til_entropy_from_internal_energy(
    float density, float u, const struct Til_params *Til) {

  return 0.f;
}

// gas_pressure_from_internal_energy
INLINE static float Til_pressure_from_internal_energy(
    float density, float u, const struct Til_params *Til) {

  const float eta = density / Til->rho_0;
  const float eta_sq = eta * eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (Til->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w;
  float P_c, P_e, P;

  // Condensed or cold
  P_c =
      (Til->a + Til->b * w_inv) * density * u + Til->A * mu + Til->B * mu * mu;
  if (eta < Til->eta_zero) {
    P_c = 0.f;
  } else if (eta < Til->eta_min) {
    P_c *= (eta - Til->eta_zero) / (Til->eta_min - Til->eta_zero);
  }
  // Expanded and hot
  P_e = Til->a * density * u +
        (Til->b * density * u * w_inv + Til->A * mu * expf(-Til->beta * nu)) *
            expf(-Til->alpha * nu * nu);

  // Condensed or cold state
  if ((1.f < eta) || (u < Til->u_iv)) {
    P = P_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (Til->u_cv < u)) {
    P = P_e;
  }
  // Hybrid state
  else {
    P = ((u - Til->u_iv) * P_e + (Til->u_cv - u) * P_c) /
        (Til->u_cv - Til->u_iv);
  }

  // Minimum pressure
  if (P < Til->P_min) {
    P = Til->P_min;
  }

  return P;
}

// gas_internal_energy_from_pressure
INLINE static float Til_internal_energy_from_pressure(
    float density, float P, const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// gas_soundspeed_from_internal_energy
INLINE static float Til_soundspeed_from_internal_energy(
    float density, float u, const struct Til_params *Til) {

  const float rho_0_inv = 1.f / Til->rho_0;
  const float eta = density * rho_0_inv;
  const float rho_inv = 1.f / density;
  const float eta_sq = eta * eta;
  const float mu = eta - 1.f;
  const float nu = 1.f / eta - 1.f;
  const float w = u / (Til->u_0 * eta_sq) + 1.f;
  const float w_inv = 1.f / w;
  const float w_inv_sq = w_inv * w_inv;
  const float exp_beta = expf(-Til->beta * nu);
  const float exp_alpha = expf(-Til->alpha * nu * nu);
  float P_c, P_e, c_sq_c, c_sq_e, c_sq;

  // Condensed or cold
  P_c =
      (Til->a + Til->b * w_inv) * density * u + Til->A * mu + Til->B * mu * mu;
  if (eta < Til->eta_zero) {
    P_c = 0.f;
  } else if (eta < Til->eta_min) {
    P_c *= (eta - Til->eta_zero) / (Til->eta_min - Til->eta_zero);
  }
  c_sq_c = P_c * rho_inv * (1.f + Til->a + Til->b * w_inv) +
           Til->b * (w - 1.f) * w_inv_sq * (2.f * u - P_c * rho_inv) +
           rho_inv * (Til->A + Til->B * (eta_sq - 1.f));

  // Expanded and hot
  P_e = Til->a * density * u +
        (Til->b * density * u * w_inv + Til->A * mu * exp_beta) * exp_alpha;

  c_sq_e =
      P_e * rho_inv * (1.f + Til->a + Til->b * w_inv * exp_alpha) +
      (Til->b * density * u * w_inv_sq / eta_sq *
           (rho_inv / Til->u_0 * (2.f * u - P_e * rho_inv) +
            2.f * Til->alpha * nu * w * rho_0_inv) +
       Til->A * rho_0_inv *
           (1.f + mu / eta_sq * (Til->beta + 2.f * Til->alpha * nu - eta)) *
           exp_beta) *
          exp_alpha;

  // Condensed or cold state
  if ((1.f < eta) || (u < Til->u_iv)) {
    c_sq = c_sq_c;
  }
  // Expanded and hot state
  else if ((eta < 1.f) && (Til->u_cv < u)) {
    c_sq = c_sq_e;
  }
  // Hybrid state
  else {
    c_sq = ((u - Til->u_iv) * c_sq_e + (Til->u_cv - u) * c_sq_c) /
           (Til->u_cv - Til->u_iv);
  }

  c_sq = fmaxf(c_sq, Til->A * rho_0_inv);

  return sqrtf(c_sq);
}

// gas_soundspeed_from_pressure
INLINE static float Til_soundspeed_from_pressure(float density, float P,
                                                 const struct Til_params *Til) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

// Compute u cold
INLINE static float compute_u_cold(float density, struct Til_params *Til,
                                   enum eos_planetary_material_id mat_id) {
  float rho_0, x, u_cold, drho;
  int N = 10000;

  rho_0 = Til->rho_0;
  drho = (density - rho_0) / N;
  x = rho_0;
  u_cold = 1e-9;

  for (int i = 0; i < N; i++) {
    x += drho;
    u_cold +=
        Til_pressure_from_internal_energy(x, u_cold, Til) * drho / (x * x);
  }

  return u_cold;
}

// Compute A1_u_cold
INLINE static void set_Til_u_cold(struct Til_params *Til,
                                  enum eos_planetary_material_id mat_id) {

  int N = 10000;
  float rho_min = 100.f;
  float rho_max = 100000.f;
  float rho, drho;

  // Allocate table memory
  Til->A1_u_cold = (float *)malloc(N * sizeof(float));

  rho = rho_min;
  drho = (rho_max - rho_min) / (N - 1);

  for (int i = 0; i < N; i++) {
    Til->A1_u_cold[i] = compute_u_cold(rho, Til, mat_id);
    rho += drho;
  }
}

// Compute u cold fast from precomputed values
INLINE static float compute_fast_u_cold(float density,
                                        const struct Til_params *Til) {

  int N = 10000;
  float rho_min = Til->rho_cold_min;
  float rho_max = Til->rho_cold_max;
  float drho, u_cold;
  int a, b;

  drho = (rho_max - rho_min) / (N - 1);

  a = (int)((density - rho_min) / drho);
  b = a + 1;

  if (a >= 0 && a < (N - 1)) {
    u_cold = Til->A1_u_cold[a];
    u_cold += ((Til->A1_u_cold[b] - Til->A1_u_cold[a]) / drho) *
              (density - rho_min - a * drho);
  } else if (density < rho_min) {
    u_cold = Til->A1_u_cold[0];
  } else {
    u_cold = Til->A1_u_cold[N - 1];
    u_cold += ((Til->A1_u_cold[N - 1] - Til->A1_u_cold[N - 2]) / drho) *
              (density - rho_max);
  }
  return u_cold;
}

// gas_temperature_from_internal_energy
INLINE static float Til_temperature_from_internal_energy(
    float density, float u, const struct Til_params *Til) {

  float u_cold, T;

  u_cold = compute_fast_u_cold(density, Til);

  T = (u - u_cold) / (Til->C_V);
  if (T < 0.f) {
    T = 0.f;
  }

  return T;
}

// gas_pressure_from_density_and_temperature
INLINE static float Til_pressure_from_temperature(
    float density, float T, const struct Til_params *Til) {

  float u, P;

  u = compute_fast_u_cold(density, Til) + Til->C_V * T;
  P = Til_pressure_from_internal_energy(density, u, Til);

  return P;
}

// gas_density_from_pressure_and_temperature
INLINE static float Til_density_from_pressure_and_temperature(
    float P, float T, const struct Til_params *Til) {

  float rho_min = Til->rho_min;
  float rho_max = Til->rho_max;
  float rho_mid = (rho_min + rho_max) / 2.f;
  float P_min, P_mid, P_max;
  float P_des;
  float tolerance = 0.001 * rho_min;
  int counter = 0;
  int max_counter = 200;
  float f0, f2;

  // Check for P == 0 or T == 0
  if (P <= Til->P_min) {
    P_des = Til->P_min;
  } else {
    P_des = P;
  }

  P_min = Til_pressure_from_temperature(rho_min, T, Til);
  P_mid = Til_pressure_from_temperature(rho_mid, T, Til);
  P_max = Til_pressure_from_temperature(rho_max, T, Til);

  // quick fix?
  if (P_des < P_min) {
    P_des = P_min;
  }

  if (P_des >= P_min && P_des <= P_max) {
    while ((rho_max - rho_min) > tolerance && counter < max_counter) {

      P_min = Til_pressure_from_temperature(rho_min, T, Til);
      P_mid = Til_pressure_from_temperature(rho_mid, T, Til);
      P_max = Til_pressure_from_temperature(rho_max, T, Til);

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
    const struct Til_params *Til) {

  if (P <= Til->P_min || u == 0) {
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
  float eta_ref = rho_ref / Til->rho_0;

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
  if (u <= Til->u_iv) {
    if (eta_ref <= Til->eta_min) {
      possible_curves[0] = 1;
      possible_curves[1] = 2;
      possible_curves[2] = 0;
    } else {
      possible_curves[0] = 2;
      possible_curves[1] = 1;
      possible_curves[2] = 0;
    }
    // u REGION 2
  } else if (u <= Til->u_cv) {
    if (eta_ref <= Til->eta_min) {
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
      if (u <= Til->u_iv) {
        if (curve == 1) {
          if (rho_iter > Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
        // u REGION 2
      } else if (u <= Til->u_cv) {
        if (curve == 3) {
          if (rho_iter > Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else if (curve == 4) {
          if (rho_iter < Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
          if (rho_iter > Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }

        // u REGION 3
      } else {
        if (curve == 5) {
          if (rho_iter > Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->rho_0) {
            rho_iter = Til->rho_0;
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

      eta_iter = rho_iter / Til->rho_0;
      eta_iter_sq = eta_iter * eta_iter;
      mu_iter = eta_iter - 1.0f;
      nu_iter = 1.0f / eta_iter - 1.0f;
      w_iter = u / (Til->u_0 * eta_iter_sq) + 1.0f;
      w_iter_inv = 1.0f / w_iter;
      exp1 = expf(-Til->beta * nu_iter);
      exp2 = expf(-Til->alpha * nu_iter * nu_iter);

      // Derivatives
      dw_inv_drho_iter =
          (2.f * Til->u_0 * u * eta_iter / Til->rho_0) /
          ((u + Til->u_0 * eta_iter_sq) * (u + Til->u_0 * eta_iter_sq));
      dmu_drho_iter = 1.f / Til->rho_0;
      dmu_sq_drho_iter =
          2.f * rho_iter / (Til->rho_0 * Til->rho_0) - 2.f / Til->rho_0;
      dexp1_drho_iter = Til->beta * Til->rho_0 * exp1 / (rho_iter * rho_iter);
      dexp2_drho_iter = 2.f * Til->alpha * Til->rho_0 *
                        (Til->rho_0 - rho_iter) * exp2 /
                        (rho_iter * rho_iter * rho_iter);

      // Use P_fraction to determine whether we've converged on a root
      float P_fraction;

      // Newton-Raphson
      float P_c_iter, dP_c_drho_iter, P_h_iter, dP_h_drho_iter;

      // if "cold" or "hybrid"
      if (curve == 2 || curve == 4) {
        P_c_iter = (Til->a + Til->b * w_iter_inv) * rho_iter * u +
                   Til->A * mu_iter + Til->B * mu_iter * mu_iter - P;
        dP_c_drho_iter = (Til->a + Til->b * w_iter_inv) * u +
                         Til->b * u * rho_iter * dw_inv_drho_iter +
                         Til->A * dmu_drho_iter + Til->B * dmu_sq_drho_iter;
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
        P_c_iter = ((Til->a + Til->b * w_iter_inv) * rho_iter * u +
                    Til->A * mu_iter + Til->B * mu_iter * mu_iter) *
                       (eta_iter - Til->eta_zero) /
                       (Til->eta_min - Til->eta_zero) -
                   P;
        dP_c_drho_iter =
            ((Til->a + Til->b * w_iter_inv) * u +
             Til->b * u * rho_iter * dw_inv_drho_iter + Til->A * dmu_drho_iter +
             Til->B * dmu_sq_drho_iter) *
                (eta_iter - Til->eta_zero) / (Til->eta_min - Til->eta_zero) +
            ((Til->a + Til->b * w_iter_inv) * rho_iter * u + Til->A * mu_iter +
             Til->B * mu_iter * mu_iter) *
                (1 / (Til->rho_0 * (Til->eta_min - Til->eta_zero)));
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
            Til->a * rho_iter * u +
            (Til->b * rho_iter * u * w_iter_inv + Til->A * mu_iter * exp1) *
                exp2 -
            P;
        dP_h_drho_iter =
            Til->a * u +
            (Til->b * u * w_iter_inv +
             Til->b * u * rho_iter * dw_inv_drho_iter +
             Til->A * mu_iter * dexp1_drho_iter +
             Til->A * exp1 * dmu_drho_iter) *
                exp2 +
            (Til->b * rho_iter * u * w_iter_inv + Til->A * mu_iter * exp1) *
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
            ((u - Til->u_iv) * P_h_iter + (Til->u_cv - u) * P_c_iter) /
            (Til->u_cv - Til->u_iv);
        float dP_hybrid_drho_iter = ((u - Til->u_iv) * dP_h_drho_iter +
                                     (Til->u_cv - u) * dP_c_drho_iter) /
                                    (Til->u_cv - Til->u_iv);
        rho_iter -= P_hybrid_iter / dP_hybrid_drho_iter;
        P_fraction = P_hybrid_iter / P;
      }

      // Now we have to constrain the new rho_iter to the curve we're on
      // u REGION 1
      if (u <= Til->u_iv) {
        if (curve == 1) {
          if (rho_iter > Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }
        // u REGION 2
      } else if (u <= Til->u_cv) {
        if (curve == 3) {
          if (rho_iter > Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
        } else if (curve == 4) {
          if (rho_iter < Til->eta_min * Til->rho_0) {
            rho_iter = Til->eta_min * Til->rho_0;
          }
          if (rho_iter > Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else {
          error("Error in Til_density_from_pressure_and_internal_energy");
        }

        // u REGION 3
      } else {
        if (curve == 5) {
          if (rho_iter > Til->rho_0) {
            rho_iter = Til->rho_0;
          }
        } else if (curve == 2) {
          if (rho_iter < Til->rho_0) {
            rho_iter = Til->rho_0;
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

// material_phase_state_from_internal_energy
INLINE static float Til_phase_state_from_internal_energy(
    float density, float u, const struct mat_params *Til,
    const struct Til_params *Til_eos) {

#ifdef MATERIAL_STRENGTH
  switch (Til->phase_state) {
    case mat_phase_state_fluid:
      return mat_phase_state_fluid;

    case mat_phase_state_solid:
      return mat_phase_state_solid;

    case mat_phase_state_variable: {
      const float T = Til_temperature_from_internal_energy(density, u, Til_eos);

      if (T > Til->T_melt) {
        return mat_phase_state_fluid;
      } else {
        return mat_phase_state_solid;
      }
    }

    default:
      return mat_phase_state_fluid;
  }
#else
  return mat_phase_state_fluid;
#endif /* MATERIAL_STRENGTH */
}

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */
