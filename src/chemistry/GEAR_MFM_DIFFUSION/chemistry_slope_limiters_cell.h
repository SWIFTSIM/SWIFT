/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_CELL_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_CELL_H

#include "hydro.h"
#include "chemistry_getters.h"

/**
 * @file src/rt/GEAR/rt_slope_limiters_cell.h
 * @brief File containing routines concerning the cell slope
 * limiter for the GEAR RT scheme. (= fist slope limiting step
 * that limits gradients such that they don't predict new extrema
 * at neighbour praticle's positions )
 *
 * */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_cell_init(
    struct part* p) {

  for (int i=0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.limiter[i].metal_density[0] = FLT_MAX;
    p->chemistry_data.limiter[i].metal_density[1] = -FLT_MAX;
  }

  p->chemistry_data.limiter_maxr = -FLT_MAX;
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * grdient loop.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_cell_collect(struct part* pi, struct part* pj, float r, int g) {

  struct chemistry_part_data* chdi = &pi->chemistry_data;

  /* Basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  chdi->limiter[g].metal_density[0] = min(chemistry_part_get_metal_density(pi, g),
                                           chdi->limiter[g].metal_density[0]);
  chdi->limiter[g].metal_density[1] = max(chemistry_part_get_metal_density(pj, g),
                                           chdi->limiter[g].metal_density[1]);

  pi->chemistry_data.limiter_maxr = max(r, pi->chemistry_data.limiter_maxr);
}


/**
 * @brief Slope-limit the given quantity. Result will be written directly
 * to double gradient[3].
 *
 * TODO: Experiment with the value of beta. For now, take stability over
 * diffusivity.
 *
 * @param gradient the gradient of the quantity
 * @param maxr maximal distance to any neighbour of the particle
 * @param value the current value of the quantity
 * @param valmin the minimal value amongst all neighbours of the quantity
 * @param valmax the maximal value amongst all neighbours of the quantity
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_quantity(
    double gradient[3], const float maxr, const double value, const double valmin,
    const double valmax, const float condition_number) {

  const double beta_min = 1.0;
  const double beta_max = 2.0;

  double gradtrue = sqrtf(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                         gradient[2] * gradient[2]);
  if (gradtrue != 0.0) {
    gradtrue *= maxr;
    const double gradtrue_inv = 1.0 / gradtrue;
    const double gradmax = valmax - value;
    const double gradmin = value - valmin;

    const double beta_2 = min(1.0, const_gizmo_max_condition_number/condition_number);
    const double beta = max(beta_min, beta_max*beta_2);
    /* const double beta = 1.0;  */
    const double min_temp =
        min(gradmax * gradtrue_inv, gradmin * gradtrue_inv) * beta;
    const double alpha = min(1.0, min_temp);
    gradient[0] *= alpha;
    gradient[1] *= alpha;
    gradient[2] *= alpha;
  }
}

/**
 * @brief Slope limit cell gradients.
 * This is done in chemistry_gradients_finalise().
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_cell(
    struct part* p) {

  const float maxr = p->chemistry_data.limiter_maxr;
  struct chemistry_part_data* chd = &p->chemistry_data;

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chemistry_slope_limit_quantity(/*gradient=*/ chd->gradients[i].nabla_otimes_q,
				   /*maxr=    */ maxr,
				   /*value=   */ chemistry_part_get_metal_density(p, i),
				   /*valmin=  */ chd->limiter[i].metal_density[0],
				   /*valmax=  */ chd->limiter[i].metal_density[1],
				   /*condition_number*/ chd->geometry.condition_number);
  }
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_CELL_H */
