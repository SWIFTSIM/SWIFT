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
#include "units.h"

// Tillotson parameters
struct Til_params {
  float rho_0, a, b, A, B, u_0, u_iv, u_cv, alpha, beta, eta_min, eta_zero,
      P_min;
  float *A1_u_cold;
  float CV;
  float rho_min_A1_u_cold, rho_max_A1_u_cold;
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
  mat->P_min = 0.01f;
  mat->CV = 449.f;
  mat->rho_min_A1_u_cold = 100.f;
  mat->rho_max_A1_u_cold = 100000.f;
  mat->rho_min = 1.f;
  mat->rho_max = 100000.f;
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
  mat->P_min = 0.01f;
  mat->CV = 790.f;
  mat->rho_min_A1_u_cold = 100.f;
  mat->rho_max_A1_u_cold = 100000.f;
  mat->rho_min = 1.f;
  mat->rho_max = 100000.f;
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
  mat->P_min = 0.01f;
  mat->CV = 790.f;
  mat->rho_min_A1_u_cold = 100.f;
  mat->rho_max_A1_u_cold = 100000.f;
  mat->rho_min = 1.f;
  mat->rho_max = 100000.f;
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
  mat->P_min = 0.01f;
  mat->CV = 4186.f;
  mat->rho_min_A1_u_cold = 100.f;
  mat->rho_max_A1_u_cold = 100000.f;
  mat->rho_min = 1.f;
  mat->rho_max = 100000.f;
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
  mat->P_min = 0.0f;
  mat->CV = 2093.f;
  mat->rho_min_A1_u_cold = 100.f;
  mat->rho_max_A1_u_cold = 100000.f;
  mat->rho_min = 1.f;
  mat->rho_max = 100000.f;
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
    mat->A1_u_cold[i] *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }
  
  mat->CV *= units_cgs_conversion_factor(&si, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin
  
  mat->rho_min_A1_u_cold *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
  mat->rho_max_A1_u_cold *= units_cgs_conversion_factor(&si, UNIT_CONV_DENSITY);
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
    mat->A1_u_cold[i] /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  }
  
  mat->CV /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY_PER_UNIT_MASS);
  // Entropy units don't work? using internal kelvin
  
  mat->rho_min_A1_u_cold /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
  mat->rho_max_A1_u_cold /= units_cgs_conversion_factor(us, UNIT_CONV_DENSITY);
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
INLINE static float compute_u_cold(float density,
                                 struct Til_params *mat,
                                 enum eos_planetary_material_id mat_id) {
    float rho_0, x, u_cold, drho;
    int N = 10000;
    
    rho_0 = mat->rho_0;
    drho = (density - rho_0) / N;
    x = rho_0;
    u_cold = 1e-9;
    
    for (int i = 0; i < N; i++) {
        x += drho;
        u_cold += Til_pressure_from_internal_energy(x, u_cold, mat) * drho / (x*x);
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
    float rho_min = mat->rho_min_A1_u_cold;
    float rho_max = mat->rho_max_A1_u_cold;
    float drho, u_cold;
    int a, b;
    
    drho = (rho_max - rho_min) / (N - 1);

    a = (int)((density - rho_min) / drho);
    b = a + 1;
    
    if (a >= 0 && a < (N - 1)){
        u_cold = mat->A1_u_cold[a];
        u_cold += ((mat->A1_u_cold[b] - mat->A1_u_cold[a]) / drho) * (
            density - rho_min - a * drho
        );
    } else if (density < rho_min){
        u_cold = mat->A1_u_cold[0];
    } else {
        u_cold = mat->A1_u_cold[N - 1];
        u_cold += (
            (mat->A1_u_cold[N - 1] - mat->A1_u_cold[N - 2]) / drho
        ) * (density - rho_max);
    }
    return u_cold;

}

// gas_temperature_from_internal_energy
INLINE static float Til_temperature_from_internal_energy(
    float density, float u, const struct Til_params *mat) {

    float u_cold, T;
    
    u_cold = compute_fast_u_cold(density, mat);
    
    T = (u - u_cold)/(mat->CV);
    if (T < 0.f){
      T = 0.f;
    }

    return T;
    
}

// gas_pressure_from_density_and_temperature
INLINE static float Til_pressure_from_temperature(
    float density, float T, const struct Til_params *mat) {

  float u, P;
  
  u = compute_fast_u_cold(density, mat) + mat->CV * T;
  P = Til_pressure_from_internal_energy(density, u, mat);
  
  return P;
}

// gas_density_from_pressure_and_temperature
INLINE static float Til_density_from_pressure_and_temperature(
    float P, float T, const struct Til_params *mat) {

    float rho_min = mat->rho_min;
    float rho_max = mat->rho_max;
    float rho_mid = (rho_min + rho_max)/2.f;
    float P_min, P_mid, P_max;
    float P_des;
    float tolerance = 0.001*rho_min;
    int counter = 0;
    int max_counter = 200;
    float f0, f2;
    
    // Check for P == 0 or T == 0
    if (P <= mat->P_min){
        P_des = mat->P_min;
    } else {
        P_des = P;
    }
    
    P_min = Til_pressure_from_temperature(rho_min, T, mat);
    P_mid = Til_pressure_from_temperature(rho_mid, T, mat);
    P_max = Til_pressure_from_temperature(rho_max, T, mat);
    
    // quick fix?
    if (P_des < P_min){
        P_des = P_min;
    }
    
    if (P_des >= P_min && P_des <= P_max){
        while ((rho_max - rho_min) > tolerance && counter < max_counter){
        
            P_min = Til_pressure_from_temperature(rho_min, T, mat);
            P_mid = Til_pressure_from_temperature(rho_mid, T, mat);
            P_max = Til_pressure_from_temperature(rho_max, T, mat);
            
            f0 = P_des - P_min;
            f2 = P_des - P_mid;

            if ((f0 * f2) > 0){
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
    float P, float u,  float rho_ref, const struct Til_params *mat) {

  error("This EOS function is not yet implemented!");

  return 0.f;
}

#endif /* SWIFT_TILLOTSON_EQUATION_OF_STATE_H */
