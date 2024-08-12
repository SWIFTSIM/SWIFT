/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_YIELD_H
#define SWIFT_PLANETARY_YIELD_H
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/strength_yield.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_damage.h"
#include "strength_utilities.h"

/**
 * @brief Adjust the yield stress depending on the density
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
adjust_yield_stress_by_temperature(struct part *restrict p, const float Y,
                                   const float density) {

  float yield_stress = Y;

  #if defined(STRENGTH_YIELD_THERMAL_SOFTENING)
    const float temperature =
        gas_temperature_from_internal_energy(density, p->u, p->mat_id);
    //  This is from Emsenhuber+2017
    //  See Senft+Stewart2007 for why this comes here and not just for intact
    float xi = method_yield_thermal_soft_xi();
    float T_melt = material_T_melt(p->mat_id);
    yield_stress *= tanhf(xi * (T_melt / temperature - 1.f));
  #endif

  return yield_stress;
}

/**
 * @brief Adjust the yield stress depending on the density
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
adjust_yield_stress_by_density(struct part *restrict p, const float Y,
                               const float density) {

  float yield_stress = Y;

  #if defined(STRENGTH_YIELD_DENSITY_SOFTENING)
    const float rho_0 = material_rho_0(p->mat_id);
    const float a = method_yield_density_soft_mult_param();
    const float b = method_yield_density_soft_pow_param();

    float rho_weak = a * rho_0;
    if (density < rho_weak) {
      yield_stress *= powf(density / rho_weak, b);
    }
  #endif

  return yield_stress;
}

/**
 * @brief Compute the intact yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_intact(
    struct part *restrict p, const float pressure) {

  float yield_stress_intact = 0.f;

  if (p->phase_state != mat_phase_state_fluid) {
    #if defined(STRENGTH_YIELD_BENZ_ASPHAUG)
      // Constant yield stress
      yield_stress_intact = material_Y_0(p->mat_id);

    #elif defined(STRENGTH_YIELD_COLLINS)
      const float mu_i = material_mu_i(p->mat_id);
      const float Y_0 = material_Y_0(p->mat_id);
      const float Y_M = material_Y_M(p->mat_id);

      // Should be able to decrease if negative pressures until yield_stress=0?
      // miluphcuda does this so maybe not wrong?
      if (pressure > 0.f) {
        yield_stress_intact = Y_0;
        if (Y_M != Y_0) {
          yield_stress_intact +=
              mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
        }
      }
    #endif
  }

  return yield_stress_intact;
}

/**
 * @brief Compute the damaged yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_damaged(
    struct part *restrict p, const float pressure,
    const float yield_stress_intact) {

  float yield_stress_damaged = 0.f;

  if (p->phase_state != mat_phase_state_fluid) {
    #if defined(STRENGTH_YIELD_COLLINS)
      const float mu_d = material_mu_d(p->mat_id);

      if (pressure > 0.f) {
        yield_stress_damaged = mu_d * pressure;
      }

      // Maybe yield_stress_damaged also needs have a max value indep of
      // yield_stress_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018

      yield_stress_damaged = min(yield_stress_damaged, yield_stress_intact);
    #endif
  }

  return yield_stress_damaged;
}

/**
 * @brief Calculates the yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress(
    struct part *restrict p, const float density, const float pressure) {

  float yield_stress = 0.f;

  if (p->phase_state != mat_phase_state_fluid) {
    #if defined(STRENGTH_YIELD_BENZ_ASPHAUG)
      // Constant yield stress
      yield_stress = material_Y_0(p->mat_id);

      yield_stress =
          adjust_yield_stress_by_density(p, yield_stress_intact, density);

      yield_stress =
          adjust_yield_stress_by_temperature(p, yield_stress, density);

    #elif defined(STRENGTH_YIELD_COLLINS)
      // ...
      float yield_stress_intact = compute_yield_stress_intact(p, pressure);
      float yield_stress_damaged =
          compute_yield_stress_damaged(p, pressure, yield_stress_intact);

      yield_stress =
          adjust_yield_stress_by_density(p, yield_stress_intact, density);

      // ...
      adjust_yield_stress_by_damage(p, yield_stress_intact, yield_stress_damaged);

      yield_stress =
          adjust_yield_stress_by_temperature(p, yield_stress, density);
    #endif
  }

  return yield_stress;
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
adjust_deviatoric_stress_tensor_by_yield_stress(
    struct part *restrict p, struct sym_matrix *deviatoric_stress_tensor,
    const float density, const float pressure) {

  #if defined(STRENGTH_YIELD_BENZ_ASPHAUG)
    // Yield stress
    float yield_stress = compute_yield_stress(p, density, pressure);

    float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);

    // ...
    float f = min((yield_stress * yield_stress) / (3.f * J_2), 1.f);

    //## should have some dt dependence?
    for (int i = 0; i < 6; i++) deviatoric_stress_tensor->elements[i] *= f;

  #elif defined(STRENGTH_YIELD_COLLINS)
    // Yield stress
    float yield_stress = compute_yield_stress(p, density, pressure);

    float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);

    // ...
    float f = min(yield_stress / sqrtf(J_2), 1.f);

    // ## should have some dt dependence?
    for (int i = 0; i < 6; i++) deviatoric_stress_tensor->elements[i] *= f;
  #endif
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_YIELD_H */
