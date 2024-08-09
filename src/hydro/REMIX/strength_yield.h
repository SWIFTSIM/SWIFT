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
adjust_yield_stress_by_density(struct part *restrict p, const float Y,
                               const float density) {

  float yield_stress = Y;

  // #if defined(STRENGTH_YIELD_###)
  //  aluminium hard-coded for now. See Luther et al. 2022 appendix. this should
  //  be 0.85*ref density.

  // ### Note that although rho_0 is an eos parameter in Til, we also want
  // ###  it to be a material parameter since it's not in most EoS
  const float rho_0 = material_rho_0(p->mat_id);
  // ### These are constants associated with the method that do no depend on
  // materials
  // ### Should these be treated differently than e.g. rho_0?
  const float a = material_yield_density_soft_mult_param(p->mat_id);
  const float b = material_yield_density_soft_pow_param(p->mat_id);

  float rho_weak = a * rho_0;
  if (density < rho_weak) {
    yield_stress *= powf(density / rho_weak, b);
  }
  // #elif defined(STRENGTH_YIELD_###)
  // #endif

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

  // #if defined(STRENGTH_YIELD_###)
  const float mu_i = material_mu_i(p->mat_id);
  const float Y_0 = material_Y_0(p->mat_id);
  const float Y_M = material_Y_M(p->mat_id);

  // Should be able to decrease if negative pressures until yield_stress=0?
  // miluphcuda does this so maybe not wrong?
  if (pressure > 0.f) {
    yield_stress_intact = Y_0;
    if (Y_M != Y_0) {
      yield_stress_intact += mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
    }
  }
  // #elif defined(STRENGTH_YIELD_###)
  // #endif

  return yield_stress_intact;
}

/**
 * @brief Compute the damaged yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_damaged(
    struct part *restrict p, const float pressure, const float yield_stress_intact) {

  float yield_stress_damaged = 0.f;

  // #if defined(STRENGTH_YIELD_###)
  const float mu_d = material_mu_d(p->mat_id);

  if (pressure > 0.f) {
    yield_stress_damaged = mu_d * pressure;
  }

  // Maybe yield_stress_damaged also needs have a max value indep of
  // yield_stress_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018

  yield_stress_damaged = min(yield_stress_damaged, yield_stress_intact);
  // #elif defined(STRENGTH_YIELD_###)
  // #endif

  return yield_stress_damaged;
}

/**
 * @brief Adjust the yield stress depending on the density
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
adjust_yield_stress_by_temperature(struct part *restrict p, const float Y,
                                   const float temperature) {

  float yield_stress = Y;

  // #if defined(STRENGTH_YIELD_###)
  //  This is from Emsenhuber+2017
  //  See Senft+Stewart2007 for why this comes here and not just for intact
  float xi = material_yield_thermal_soft_xi(p->mat_id);
  float T_melt = material_T_melt(p->mat_id);
  yield_stress *= tanhf(xi * (T_melt / temperature - 1.f));
  // #endif

  return yield_stress;
}

/**
 * @brief Calculates the yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress(
    struct part *restrict p, const float density, const float pressure,
    const float temperature) {

  float yield_stress = 0.f;

  // #if defined(STRENGTH_YIELD_###)
  if (p->phase_state != mat_phase_state_fluid) {
    // Jutzi 2015 notation
    // From Collins (2004)

    float yield_stress_intact = compute_yield_stress_intact(p, pressure);
    float yield_stress_damaged = compute_yield_stress_damaged(p, pressure, yield_stress_intact);

    yield_stress = adjust_yield_stress_by_density(p, yield_stress_intact, density);

    // ...
    yield_stress = (1.f - p->damage) * yield_stress_intact +
                   p->damage * yield_stress_damaged;

    yield_stress = adjust_yield_stress_by_temperature(p, yield_stress, temperature);
  }
  // #elif defined(STRENGTH_YIELD_###)
  // #endif

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
    const float density, const float pressure, const float temperature) {

    // Yield stress
    float yield_stress = compute_yield_stress(p, density, pressure, temperature);
  /*
//#if defined(STRENGTH_YIELD_###)
    float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);

    // ...
    float f = min((yield_stress * yield_stress) / (3.f * J_2), 1.f);

    //## should have some dt dependence?
    for (int i = 0; i < 6; i++) deviatoric_stress_tensor->elements[i] *= f;
  */
  // #elif defined(STRENGTH_YIELD_###)
  float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);

  // ...
  float f = min(yield_stress / sqrtf(J_2), 1.f);

  // ## should have some dt dependence?
  for (int i = 0; i < 6; i++) {
    deviatoric_stress_tensor->elements[i] *= f;
  }

  // #endif
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_YIELD_H */
