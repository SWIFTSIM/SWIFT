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

#include "chemistry_getters.h"
#include "hydro.h"

#define GIZMO_SLOPE_LIMITER_BETA_MIN 1.0
#define GIZMO_SLOPE_LIMITER_BETA_MAX 2.0

/**
 * @file src/chemistry/GEAR_MFM_FIDDUSION/chemistry_slope_limiters_cell.h
 * @brief File containing routines concerning the cell slope
 * limiter for the GEAR MFM diffusion scheme. (= fist slope limiting step
 * that limits gradients such that they don't predict new extrema
 * at neighbour praticle's positions )
 *
 * */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_cell_init(struct part* p) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.limiter.metal_density[i][0] = FLT_MAX;
    p->chemistry_data.limiter.metal_density[i][1] = -FLT_MAX;
  }

  p->chemistry_data.limiter.v[0][0] = FLT_MAX;
  p->chemistry_data.limiter.v[0][1] = -FLT_MAX;
  p->chemistry_data.limiter.v[1][0] = FLT_MAX;
  p->chemistry_data.limiter.v[1][1] = -FLT_MAX;
  p->chemistry_data.limiter.v[2][0] = FLT_MAX;
  p->chemistry_data.limiter.v[2][1] = -FLT_MAX;

  p->chemistry_data.limiter.maxr = -FLT_MAX;
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * grdient loop.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_cell_collect(struct part* pi, struct part* pj, float r) {

  struct chemistry_part_data* chdi = &pi->chemistry_data;

  /* Basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chdi->limiter.metal_density[i][0] =
        min(chemistry_part_get_metal_density(pi, i),
            chdi->limiter.metal_density[i][0]);
    chdi->limiter.metal_density[i][1] =
        max(chemistry_part_get_metal_density(pj, i),
            chdi->limiter.metal_density[i][1]);
  }

  pi->chemistry_data.limiter.v[0][0] =
      min(pj->v[0], pi->chemistry_data.limiter.v[0][0]);
  pi->chemistry_data.limiter.v[0][1] =
      max(pj->v[0], pi->chemistry_data.limiter.v[0][1]);
  pi->chemistry_data.limiter.v[1][0] =
      min(pj->v[1], pi->chemistry_data.limiter.v[1][0]);
  pi->chemistry_data.limiter.v[1][1] =
      max(pj->v[1], pi->chemistry_data.limiter.v[1][1]);
  pi->chemistry_data.limiter.v[2][0] =
      min(pj->v[2], pi->chemistry_data.limiter.v[2][0]);
  pi->chemistry_data.limiter.v[2][1] =
      max(pj->v[2], pi->chemistry_data.limiter.v[2][1]);

  pi->chemistry_data.limiter.maxr = max(r, pi->chemistry_data.limiter.maxr);
}

/**
 * @brief Slope-limit the given quantity. Result will be written directly
 * to double gradient[3].
 *
 * TODO: Experiment with the value of beta. For now, take stability over
 * diffusivity with beta=1.
 *
 * @param gradient the gradient of the quantity
 * @param maxr maximal distance to any neighbour of the particle
 * @param value the current value of the quantity
 * @param valmin the minimal value amongst all neighbours of the quantity
 * @param valmax the maximal value amongst all neighbours of the quantity
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_quantity(double gradient[3], const float maxr,
                               const double value, const double valmin,
                               const double valmax,
                               const float condition_number) {

  double gradtrue =
      sqrtf(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
            gradient[2] * gradient[2]);
  if (gradtrue != 0.0) {
    gradtrue *= maxr;
    const double gradtrue_inv = 1.0 / gradtrue;
    const double gradmax = valmax - value;
    const double gradmin = value - valmin;
    const double beta_2 =
        min(1.0, const_gizmo_max_condition_number / condition_number);
    const double beta = max(GIZMO_SLOPE_LIMITER_BETA_MIN,
                            GIZMO_SLOPE_LIMITER_BETA_MAX * beta_2);
    /* const double beta = 1.0; /\* Choose stability *\/ */
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

  struct chemistry_part_data* chd = &p->chemistry_data;
  const float N_cond = chd->geometry.condition_number;
  const float maxr = chd->limiter.maxr;
  const float vxlim[2] = {chd->limiter.v[0][0], chd->limiter.v[0][1]};
  const float vylim[2] = {chd->limiter.v[1][0], chd->limiter.v[1][1]};
  const float vzlim[2] = {chd->limiter.v[2][0], chd->limiter.v[2][1]};

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chemistry_slope_limit_quantity(
        /*gradient=*/chd->gradients.nabla_otimes_q[i],
        /*maxr=    */ maxr,
        /*value=   */ chemistry_part_get_metal_density(p, i),
        /*valmin=  */ chd->limiter.metal_density[i][0],
        /*valmax=  */ chd->limiter.metal_density[i][1],
        /*condition_number*/ N_cond);
  }

  /* Use doubles sice chemistry_slope_limit_quantity() accepts double arrays. */
  double gradvx[3], gradvy[3], gradvz[3];

  /* Get the velocity gradients and cast them as double */
  gradvx[0] = p->chemistry_data.gradients.v[0][0];
  gradvx[1] = p->chemistry_data.gradients.v[0][1];
  gradvx[2] = p->chemistry_data.gradients.v[0][2];
  gradvy[0] = p->chemistry_data.gradients.v[1][0];
  gradvy[1] = p->chemistry_data.gradients.v[1][1];
  gradvy[2] = p->chemistry_data.gradients.v[1][2];
  gradvz[0] = p->chemistry_data.gradients.v[2][0];
  gradvz[1] = p->chemistry_data.gradients.v[2][1];
  gradvz[2] = p->chemistry_data.gradients.v[2][2];

  /* Slope limit the velocity gradient */
  chemistry_slope_limit_quantity(gradvx, maxr, p->v[0], vxlim[0], vxlim[1],
                                 N_cond);
  chemistry_slope_limit_quantity(gradvy, maxr, p->v[1], vylim[0], vylim[1],
                                 N_cond);
  chemistry_slope_limit_quantity(gradvz, maxr, p->v[2], vzlim[0], vzlim[1],
                                 N_cond);

  /* Set the velocity gradient values */
  p->chemistry_data.gradients.v[0][0] = gradvx[0];
  p->chemistry_data.gradients.v[0][1] = gradvx[1];
  p->chemistry_data.gradients.v[0][2] = gradvx[2];
  p->chemistry_data.gradients.v[1][0] = gradvy[0];
  p->chemistry_data.gradients.v[1][1] = gradvy[1];
  p->chemistry_data.gradients.v[1][2] = gradvy[2];
  p->chemistry_data.gradients.v[2][0] = gradvz[0];
  p->chemistry_data.gradients.v[2][1] = gradvz[1];
  p->chemistry_data.gradients.v[2][2] = gradvz[2];
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_CELL_H */
