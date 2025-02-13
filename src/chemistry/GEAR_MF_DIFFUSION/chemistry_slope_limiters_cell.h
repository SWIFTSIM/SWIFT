/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SLOPE_LIMITERS_CELL_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SLOPE_LIMITERS_CELL_H

#include "chemistry_getters.h"
#include "hydro.h"

#define GIZMO_SLOPE_LIMITER_BETA_MIN 1.0
#define GIZMO_SLOPE_LIMITER_BETA_MAX 2.0

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/chemistry_slope_limiters_cell.h
 * @brief File containing routines concerning the cell slope limiter for the
 * GEAR MF diffusion scheme. (= fist slope limiting step that limits gradients
 * such that they don't predict new extrema at neighbour particle's positions)
 * */

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_cell_init(struct part* p) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    p->chemistry_data.limiter.Z[i][0] = FLT_MAX;
    p->chemistry_data.limiter.Z[i][1] = -FLT_MAX;
  }

  p->chemistry_data.limiter.rho[0] = FLT_MAX;
  p->chemistry_data.limiter.rho[1] = -FLT_MAX;

  p->chemistry_data.limiter.v[0][0] = FLT_MAX;
  p->chemistry_data.limiter.v[0][1] = -FLT_MAX;
  p->chemistry_data.limiter.v[1][0] = FLT_MAX;
  p->chemistry_data.limiter.v[1][1] = -FLT_MAX;
  p->chemistry_data.limiter.v[2][0] = FLT_MAX;
  p->chemistry_data.limiter.v[2][1] = -FLT_MAX;

  p->chemistry_data.limiter.v_tilde[0][0] = FLT_MAX;
  p->chemistry_data.limiter.v_tilde[0][1] = -FLT_MAX;
  p->chemistry_data.limiter.v_tilde[1][0] = FLT_MAX;
  p->chemistry_data.limiter.v_tilde[1][1] = -FLT_MAX;
  p->chemistry_data.limiter.v_tilde[2][0] = FLT_MAX;
  p->chemistry_data.limiter.v_tilde[2][1] = -FLT_MAX;

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

  struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;

  /* Basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chi->limiter.Z[i][0] =
      min(chemistry_get_metal_mass_fraction(pj, i), chi->limiter.Z[i][0]);
    chi->limiter.Z[i][1] =
      max(chemistry_get_metal_mass_fraction(pj, i), chi->limiter.Z[i][1]);

    /* Ensures metalicity never exceeds 1 */
    chi->limiter.Z[i][1] =
      min(1.0, chi->limiter.Z[i][1]);
  }

  chi->limiter.v[0][0] = min(pj->v[0], chi->limiter.v[0][0]);
  chi->limiter.v[0][1] = max(pj->v[0], chi->limiter.v[0][1]);
  chi->limiter.v[1][0] = min(pj->v[1], chi->limiter.v[1][0]);
  chi->limiter.v[1][1] = max(pj->v[1], chi->limiter.v[1][1]);
  chi->limiter.v[2][0] = min(pj->v[2], chi->limiter.v[2][0]);
  chi->limiter.v[2][1] = max(pj->v[2], chi->limiter.v[2][1]);

  chi->limiter.rho[0] = min(pj->rho, chi->limiter.rho[0]);
  chi->limiter.rho[1] = max(pj->rho, chi->limiter.rho[1]);

  chi->limiter.v_tilde[0][0] = min(chj->filtered.rho_v[0] / chj->filtered.rho,
                                   chi->limiter.v_tilde[0][0]);
  chi->limiter.v_tilde[0][1] = max(chj->filtered.rho_v[0] / chj->filtered.rho,
                                   chi->limiter.v_tilde[0][1]);
  chi->limiter.v_tilde[1][0] = min(chj->filtered.rho_v[1] / chj->filtered.rho,
                                   chi->limiter.v_tilde[1][0]);
  chi->limiter.v_tilde[1][1] = max(chj->filtered.rho_v[1] / chj->filtered.rho,
                                   chi->limiter.v_tilde[1][1]);
  chi->limiter.v_tilde[2][0] = min(chj->filtered.rho_v[2] / chj->filtered.rho,
                                   chi->limiter.v_tilde[2][0]);
  chi->limiter.v_tilde[2][1] = max(chj->filtered.rho_v[2] / chj->filtered.rho,
                                   chi->limiter.v_tilde[2][1]);

  chi->limiter.maxr = max(r, chi->limiter.maxr);
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
 * @param condition_number Condition number of the particle geometric matrix E
 * @param pos_preserve Flag to use a positivity-preserving slope limiter for
 * quantities that should not become negative (e.g. densities, metallicities).
 */
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_quantity(double gradient[3], const float maxr,
                               const double value, const double valmin,
                               const double valmax,
                               const float condition_number,
                               const int pos_preserve) {

  double gradtrue = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                         gradient[2] * gradient[2]);

  /* If valmin == valmax, this is ill-defined. Do not slope-limit.
     This case can happen e.g. in the MetalDiffusionOnePeak example where we
     only have one particule with non zero metallicity. The slope limiter then
     computes a high value for alpha (1e15). */
  if (gradtrue != 0.0) {
    gradtrue *= maxr;
    const double gradtrue_inv = 1.0 / gradtrue;
    const double gradmax = valmax - value;
    const double gradmin = value - valmin;

    /* Stability factor based on condition number (better than 1.0) */
    const double beta_2 =
        min(1.0, const_gizmo_max_condition_number / condition_number);
    const double beta = max(GIZMO_SLOPE_LIMITER_BETA_MIN,
                            GIZMO_SLOPE_LIMITER_BETA_MAX * beta_2);

    /* Base slope limiter */
    const double min_temp =
        min(gradmax * gradtrue_inv, gradmin * gradtrue_inv) * beta;
    double alpha = min(1.0, min_temp);

    /* Positivity-preserving mechanism */
    if (pos_preserve) {
      const double shoot_tol = 0.0;  // Allowable overshoot tolerance
      const double overshoot_corr = gradmin + shoot_tol * gradmax;
      const double f_corr_overshoot =
          (overshoot_corr < gradmax) ? overshoot_corr : gradmax;

      /* Compute fmin for positivity preservation */
      const double min_value =
          (value + valmin) * 0.5;  // Halfway between central and min value
      const double positive_definite_min =
          (value > 0.0) ? value * FLT_MIN : 0.0;
      const double fmin =
          (f_corr_overshoot < min_value) ? f_corr_overshoot : min_value;
      const double final_fmin =
          (fmin > positive_definite_min) ? fmin : positive_definite_min;

      /* Compute cfac for positivity preservation
	 Note: gradtrue has been multiplied by maxr already. Hence, we don't need
	 to divide (value - final_fmin) by maxr. */
      const double cfac_pos_preserve = (value - final_fmin) / gradtrue;
      alpha = (cfac_pos_preserve < alpha) ? cfac_pos_preserve : alpha;
    }

    /* Apply the limiting factor to the gradient */
    if (alpha < 1.0) {
      gradient[0] *= alpha;
      gradient[1] *= alpha;
      gradient[2] *= alpha;
    }
  } /* gradtrue != 0 */
}

/**
 * @brief Slope limit cell gradients.
 * This is done in chemistry_gradients_finalise().
 *
 * @param p Particle.
 * @param cd The global properties of the chemistry scheme.
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_cell(
    struct part* p, const struct chemistry_global_data* cd) {

  struct chemistry_part_data* chd = &p->chemistry_data;
  const float N_cond = chd->geometry_condition_number;
  const float maxr = chd->limiter.maxr;
  const float vxlim[2] = {chd->limiter.v[0][0], chd->limiter.v[0][1]};
  const float vylim[2] = {chd->limiter.v[1][0], chd->limiter.v[1][1]};
  const float vzlim[2] = {chd->limiter.v[2][0], chd->limiter.v[2][1]};
  const float rholim[2] = {chd->limiter.rho[0], chd->limiter.rho[1]};
  const float vx_tilde_lim[2] = {chd->limiter.v_tilde[0][0],
                                 chd->limiter.v_tilde[0][1]};
  const float vy_tilde_lim[2] = {chd->limiter.v_tilde[1][0],
                                 chd->limiter.v_tilde[1][1]};
  const float vz_tilde_lim[2] = {chd->limiter.v_tilde[2][0],
                                 chd->limiter.v_tilde[2][1]};
  const float v_tilde[3] = {chd->filtered.rho_v[0] / chd->filtered.rho,
                            chd->filtered.rho_v[1] / chd->filtered.rho,
                            chd->filtered.rho_v[2] / chd->filtered.rho};

#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
  /* The positivity preserving better reduces artificial diffusion */
  const int mode = 1;
#else
  /* For parabolic diffusion, it's better to use the regular cell slope
     limiter. The Riemann solver takes care of artifical diffusion. */
  const int mode = 0;
#endif
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    chemistry_slope_limit_quantity(
        /*gradient=*/ chd->gradients.Z[i],
        /*maxr=    */ maxr,
        /*value=   */ chemistry_get_metal_mass_fraction(p, i),
        /*valmin=  */ chd->limiter.Z[i][0],
        /*valmax=  */ chd->limiter.Z[i][1],
        /*condition_number*/ N_cond,
        /*pos_preserve*/ mode);
  }

  /* Use doubles sice chemistry_slope_limit_quantity() accepts double arrays. */
  double gradrho[3], gradvx[3], gradvy[3], gradvz[3], gradvx_tilde[3], gradvy_tilde[3],
      gradvz_tilde[3];

  /* Get the velocity gradients and cast them as double */
  gradvx[0] = chd->gradients.v[0][0];
  gradvx[1] = chd->gradients.v[0][1];
  gradvx[2] = chd->gradients.v[0][2];
  gradvy[0] = chd->gradients.v[1][0];
  gradvy[1] = chd->gradients.v[1][1];
  gradvy[2] = chd->gradients.v[1][2];
  gradvz[0] = chd->gradients.v[2][0];
  gradvz[1] = chd->gradients.v[2][1];
  gradvz[2] = chd->gradients.v[2][2];

  gradrho[0] = chd->gradients.rho[0];
  gradrho[1] = chd->gradients.rho[1];
  gradrho[2] = chd->gradients.rho[2];

  gradvx_tilde[0] = chd->filtered.grad_v_tilde[0][0];
  gradvx_tilde[1] = chd->filtered.grad_v_tilde[0][1];
  gradvx_tilde[2] = chd->filtered.grad_v_tilde[0][2];
  gradvy_tilde[0] = chd->filtered.grad_v_tilde[1][0];
  gradvy_tilde[1] = chd->filtered.grad_v_tilde[1][1];
  gradvy_tilde[2] = chd->filtered.grad_v_tilde[1][2];
  gradvz_tilde[0] = chd->filtered.grad_v_tilde[2][0];
  gradvz_tilde[1] = chd->filtered.grad_v_tilde[2][1];
  gradvz_tilde[2] = chd->filtered.grad_v_tilde[2][2];

  /* Slope limit the velocity gradients */
  chemistry_slope_limit_quantity(gradvx, maxr, p->v[0], vxlim[0], vxlim[1],
                                 N_cond, 0);
  chemistry_slope_limit_quantity(gradvy, maxr, p->v[1], vylim[0], vylim[1],
                                 N_cond, 0);
  chemistry_slope_limit_quantity(gradvz, maxr, p->v[2], vzlim[0], vzlim[1],
                                 N_cond, 0);
  chemistry_slope_limit_quantity(gradvx_tilde, maxr, v_tilde[0],
                                 vx_tilde_lim[0], vx_tilde_lim[1], N_cond, 0);
  chemistry_slope_limit_quantity(gradvy_tilde, maxr, v_tilde[1],
                                 vy_tilde_lim[0], vy_tilde_lim[1], N_cond, 0);
  chemistry_slope_limit_quantity(gradvz_tilde, maxr, v_tilde[2],
                                 vz_tilde_lim[0], vz_tilde_lim[1], N_cond, 0);

  /* Slope limit density gradient */
  chemistry_slope_limit_quantity(gradrho, maxr, p->rho, rholim[0], rholim[1], N_cond, 1);

  /* Set the velocity gradient values */
  chd->gradients.v[0][0] = gradvx[0];
  chd->gradients.v[0][1] = gradvx[1];
  chd->gradients.v[0][2] = gradvx[2];
  chd->gradients.v[1][0] = gradvy[0];
  chd->gradients.v[1][1] = gradvy[1];
  chd->gradients.v[1][2] = gradvy[2];
  chd->gradients.v[2][0] = gradvz[0];
  chd->gradients.v[2][1] = gradvz[1];
  chd->gradients.v[2][2] = gradvz[2];

  chd->filtered.grad_v_tilde[0][0] = gradvx_tilde[0];
  chd->filtered.grad_v_tilde[0][1] = gradvx_tilde[1];
  chd->filtered.grad_v_tilde[0][2] = gradvx_tilde[2];
  chd->filtered.grad_v_tilde[1][0] = gradvy_tilde[0];
  chd->filtered.grad_v_tilde[1][1] = gradvy_tilde[1];
  chd->filtered.grad_v_tilde[1][2] = gradvy_tilde[2];
  chd->filtered.grad_v_tilde[2][0] = gradvz_tilde[0];
  chd->filtered.grad_v_tilde[2][1] = gradvz_tilde[1];
  chd->filtered.grad_v_tilde[2][2] = gradvz_tilde[2];
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_SLOPE_LIMITERS_CELL_H */
