/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_VISCOSITY_H
#define SWIFT_PLANETARY_HYDRO_VISCOSITY_H

/**
 * @file Planetary/hydro_viscosity.h
 * @brief Utilities for hydro viscosity calculations under various options.
 */

#include "const.h"
#include "hydro_misc_utils.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Prepares extra viscosity parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_viscosity(struct part *restrict p) {

#ifdef PLANETARY_QUAD_VISC
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      p->Dinv[i][j] = 0.f;
      p->E[i][j] = 0.f;
    }
  }
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Extra density interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_viscosity(struct part *restrict pi,
                                          struct part *restrict pj,
                                          const float dx[3], const float wi,
                                          const float wj, const float wi_dx,
                                          const float wj_dx) {

#ifdef PLANETARY_QUAD_VISC
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float hi_inv = 1.f / pi->h;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;

  const float hj_inv = 1.f / pj->h;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float wj_dr = hjd_inv * wj_dx;

  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Inverse of D matrix (eq 20 in Rosswog 2020) */
      pi->Dinv[i][j] += pj->mass * dx[i] * dx[j] * wi_dr * r_inv;
      pj->Dinv[i][j] += pi->mass * dx[i] * dx[j] * wj_dr * r_inv;

      /* E matrix (second part of eq 19 in Rosswog 2020) */
      pi->E[i][j] += pj->mass * (pi->v[i] - pj->v[i]) * dx[j] * wi_dr * r_inv;
      pj->E[i][j] += pi->mass * (pi->v[i] - pj->v[i]) * dx[j] * wj_dr * r_inv;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  // ## This can be changed if swift matrix inversion function works
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Dinv[2][2] = 1.f;
  pj->Dinv[2][2] = 1.f;
#endif

#endif
}

/**
 * @brief Extra density interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_viscosity(struct part *restrict pi,
                                                 const struct part *restrict pj,
                                                 const float dx[3],
                                                 const float wi,
                                                 const float wi_dx) {

#ifdef PLANETARY_QUAD_VISC
  const float r = sqrtf(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float h_inv = 1.f / pi->h;
  const float hid_inv = pow_dimension_plus_one(h_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;

  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Inverse of D matrix (eq 20 in Rosswog 2020) */
      pi->Dinv[i][j] += pj->mass * dx[i] * dx[j] * wi_dr * r_inv;

      /* E matrix (second part of eq 19 in Rosswog 2020) */
      pi->E[i][j] += pj->mass * (pi->v[i] - pj->v[i]) * dx[j] * wi_dr * r_inv;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  // ## This can be changed if swift matrix inversion function works
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Dinv[2][2] = 1.f;
#endif

#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Finishes extra viscosity parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_viscosity(struct part *restrict p) {

#ifdef PLANETARY_QUAD_VISC
  // ## Can we do this with inbuilt function?

  int i, j, k;

  /* In this section we carry out matrix inversion to find D and calculate
   * dv_aux */

  float determinant = 0.f;
  /* We normalise the Dinv matrix to the mean of its 9 elements to stop us
   * hitting float precision limits during matrix inversion process */
  float mean_Dinv = (p->Dinv[0][0] + p->Dinv[0][1] + p->Dinv[0][2] +
                     p->Dinv[1][0] + p->Dinv[1][1] + p->Dinv[1][2] +
                     p->Dinv[2][0] + p->Dinv[2][1] + p->Dinv[2][2]) /
                    9.f;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Normalise Dinv to mean of its values */
      p->Dinv[i][j] = p->Dinv[i][j] / mean_Dinv;

      /* Aux dv (eq 19 in Rosswog 2020) */
      p->dv_aux[i][j] = 0.f;
    }
  }

  for (i = 0; i < 3; i++) {
    /* Matrix Dinv det */
    determinant +=
        (p->Dinv[0][i] * (p->Dinv[1][(i + 1) % 3] * p->Dinv[2][(i + 2) % 3] -
                          p->Dinv[1][(i + 2) % 3] * p->Dinv[2][(i + 1) % 3]));
  }

  float D[3][3];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Find D from inverse of Dinv */
      D[i][j] = ((p->Dinv[(i + 1) % 3][(j + 1) % 3] *
                  p->Dinv[(i + 2) % 3][(j + 2) % 3]) -
                 (p->Dinv[(i + 1) % 3][(j + 2) % 3] *
                  p->Dinv[(i + 2) % 3][(j + 1) % 3])) /
                (determinant * mean_Dinv);
      if (isnan(D[i][j]) || isinf(D[i][j])) {
        D[i][j] = 0.f;
      }

      for (k = 0; k < 3; ++k) {
        /* Calculate dv_aux (eq 19 in Rosswog 2020) */
        p->dv_aux[i][k] += D[i][j] * p->E[k][j];
      }
    }
  }
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Prepares extra viscosity parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_viscosity(struct part *restrict p) {

#ifdef PLANETARY_QUAD_VISC
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->dv_no_C[i][j] = 0.f;

      for (int k = 0; k < 3; ++k) {
        p->ddv_no_C[i][j][k] = 0.f;
      }
    }
  }

  p->N_grad = 0.f;
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Extra gradient interaction between two particles
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_viscosity(struct part *restrict pi,
                                           struct part *restrict pj,
                                           const float dx[3], const float wi,
                                           const float wj, const float wi_dx,
                                           const float wj_dx) {

#ifdef PLANETARY_QUAD_VISC
  float volume_i = pi->mass / pi->rho;
  float volume_j = pj->mass / pj->rho;

  planetary_smoothing_correction_tweak_volume(&volume_i, pi);
  planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  /* Set velocity derivative elements */
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      /* Gradients from eq 18 in Rosswog 2020 (without C multiplied) */
      pi->dv_no_C[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wi * volume_j;
      pj->dv_no_C[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wj * volume_i;

      for (int k = 0; k < 3; ++k) {
        /* Gradients from eq 18 in Rosswog 2020 (without C multiplied). Note
         * that we now use dv_aux to get second derivative*/
        pi->ddv_no_C[i][j][k] +=
            (pi->dv_aux[i][j] - pj->dv_aux[i][j]) * dx[k] * wi * volume_j;
        pj->ddv_no_C[i][j][k] +=
            (pi->dv_aux[i][j] - pj->dv_aux[i][j]) * dx[k] * wj * volume_i;
      }
    }
  }

  /* Number of neighbours. Needed for eta_crit factor in slope limiter */
  pi->N_grad += 1.f;
  pj->N_grad += 1.f;
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Extra gradient interaction between two particles (non-symmetric)
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_viscosity(
    struct part *restrict pi, const struct part *restrict pj, const float dx[3],
    const float wi, const float wi_dx) {

#ifdef PLANETARY_QUAD_VISC

  float volume_j = pj->mass / pj->rho;

  planetary_smoothing_correction_tweak_volume(&volume_j, pj);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      /* Gradients from eq 18 in Rosswog 2020 (without C multiplied)*/
      pi->dv_no_C[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wi * volume_j;

      for (int k = 0; k < 3; ++k) {
        /* Gradients from eq 18 in Rosswog 2020 (without C multiplied). Note
         * that we now use dv_aux to get second derivative*/
        pi->ddv_no_C[i][j][k] +=
            (pi->dv_aux[i][j] - pj->dv_aux[i][j]) * dx[k] * wi * volume_j;
      }
    }
  }

  /* Number of neighbours. Needed for eta_crit factor in slope limiter */
  pi->N_grad += 1.f;
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Finishes extra viscosity parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_viscosity(struct part *restrict p) {

#ifdef PLANETARY_QUAD_VISC
  int i, j, k, l;

  /* In this section we:
      1) calculate dv (eq 18 in Rosswog 2020);
      2) calculate ddv (eq 18 in Rosswog 2020 but using the dv_aux instead of
         v to get second derivative).
  */

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      p->dv[i][j] = 0.f;
      for (k = 0; k < 3; ++k) {
        p->ddv[i][j][k] = 0.f;
      }
    }
  }

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!p->is_h_max) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; ++k) {
          /* calculate dv (eq 18 in Rosswog 2020) */
          p->dv[i][k] += p->C[i][j] * p->dv_no_C[k][j];
          for (l = 0; l < 3; ++l) {
            /* calculate ddv (eq 18 in Rosswog 2020) */
            p->ddv[i][l][k] += p->C[i][j] * p->ddv_no_C[l][k][j];
          }
        }
      }
    }
  }
#endif /* PLANETARY_QUAD_VISC */
}

/**
 * @brief Returns particle viscous pressures
 */
__attribute__((always_inline)) INLINE static void hydro_set_Qi_Qj(
    float *Qi, float *Qj, const struct part *restrict pi,
    const struct part *restrict pj, const float dx[3], const float a,
    const float H) {

#ifdef PLANETARY_QUAD_VISC
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);

  const float hi_inv = 1.0f / pi->h;
  const float hj_inv = 1.0f / pj->h;

  /* Quadratically reconstructed velocities at the halfway point between
   * particles */
  float vtilde_i[3], vtilde_j[3];

  /* Some parameters for artificial visc. Taken from Rosswog 2020 */
  float alpha = 1.f;
  float beta = 2.f;
  float epsilon = 0.1;

  /* Square of eta (eq 16 in Rosswog 2020) */
  float eta_i_2 = r2 * hi_inv * hi_inv;
  float eta_j_2 = r2 * hj_inv * hj_inv;

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!pi->is_h_max && !pj->is_h_max) {

    /* eq 23 in Rosswog 2020 */
    float eta_ab = min(r * hi_inv, r * hj_inv);

    /* A numerators and denominators (eq 22 in Rosswog 2020) */
    float A_i_v = 0.f;
    float A_j_v = 0.f;

    /* Terms in square brackets in Rosswog 2020 eq 17 */
    float v_quad_i[3] = {0};
    float v_quad_j[3] = {0};

    /* eq 23 in Rosswog 2020 set to constant */
    float eta_crit = 0.5f;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        /* Get the A numerators and denominators (eq 22 in Rosswog 2020). dv
         * is from eq 18 */
        A_i_v += pi->dv[i][j] * dx[i] * dx[j];
        A_j_v += pj->dv[i][j] * dx[i] * dx[j];

        /* Terms in square brackets in Rosswog 2020 eq 17. Add in FIRST
         * derivative terms */
        v_quad_i[j] -= 0.5 * pi->dv[i][j] * dx[i];
        v_quad_j[j] += 0.5 * pj->dv[i][j] * dx[i];

        for (int k = 0; k < 3; ++k) {
          /* Terms in square brackets in Rosswog 2020 eq 17. Add in SECOND
           * derivative terms */
          v_quad_i[j] += 0.125 * pi->ddv[i][j][k] * dx[i] * dx[k];
          v_quad_j[j] += 0.125 * pj->ddv[i][j][k] * dx[i] * dx[k];
        }
      }
    }

    /* Slope limiter (eq 21 in Rosswog 2020) */
    float phi_i_v =
        min(1.f, 4 * A_i_v / A_j_v / (1 + A_i_v / A_j_v) / (1 + A_i_v / A_j_v));
    phi_i_v = max(0.f, phi_i_v);

    float phi_j_v =
        min(1.f, 4 * A_j_v / A_i_v / (1 + A_j_v / A_i_v) / (1 + A_j_v / A_i_v));
    phi_j_v = max(0.f, phi_j_v);

    if (eta_ab < eta_crit) {
      phi_i_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
      phi_j_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
    }
    /* These are here to catch division by 0. In this case phi tends to 0 anyway
     */
    if (isnan(phi_i_v) || isinf(phi_i_v)) {
      phi_i_v = 0.f;
    }
    if (isnan(phi_j_v) || isinf(phi_j_v)) {
      phi_j_v = 0.f;
    }

    for (int i = 0; i < 3; ++i) {
      /* Assemble the reconstructed velocity (eq 17 in Rosswog 2020) */
      vtilde_i[i] = pi->v[i] + phi_i_v * v_quad_i[i];
      vtilde_j[i] = pj->v[i] + phi_j_v * v_quad_j[i];
    }

  } else {
    for (int i = 0; i < 3; ++i) {
      /* If h=h_max don't reconstruct velocity */
      vtilde_i[i] = pi->v[i];
      vtilde_j[i] = pj->v[i];
    }
  }

  /* Finally assemble eq 15 in Rosswog 2020 */
  float mu_i = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hi_inv / (eta_i_2 + epsilon * epsilon));
  float mu_j = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hj_inv / (eta_j_2 + epsilon * epsilon));

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Get viscous pressure terms (eq 14 in Rosswog 2020) */
  *Qi = pi->rho * (-alpha * ci * mu_i + beta * mu_i * mu_i);
  *Qj = pj->rho * (-alpha * cj * mu_j + beta * mu_j * mu_j);

#else /* !PLANETARY_QUAD_VISC */

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;

  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij;
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Density factors for GDF or standard equations */
  float rho_factor_i, rho_factor_j;
  hydro_set_rho_factors(&rho_factor_i, &rho_factor_j, pi, pj);

  /* Artificial viscosity terms, as pressure in Rosswog2020 framework, S2.2.1 */
  *Qi = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rho_factor_i /
        (pi->rho + pj->rho);
  *Qj = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rho_factor_j /
        (pi->rho + pj->rho);

#endif /* PLANETARY_QUAD_VISC */
}

#endif /* SWIFT_PLANETARY_HYDRO_VISCOSITY_H */
