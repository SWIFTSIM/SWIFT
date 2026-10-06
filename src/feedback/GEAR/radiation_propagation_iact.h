/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H
#define SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_propagation_iact.h
 * @brief Gas-gas density-, gradient- and force-loop hooks for the hyperbolic
 * M1-relaxation propagation of the per-band u and specific_flux fields.
 *
 * `div(F)` is in the force loop because it needs both sides of every pair
 * credited, which only the force loop's type-2 dispatch guarantees. The
 * stored flux is the reduced flux `Ft = F_true/c_hyp`; each pair term uses the
 * receiver's own `c_hyp`. Inputs are comoving, accumulators physical (factor
 * `1/a`). All loops read `feedback_data.rho_prev`, not `p->rho`, which is a
 * partial sum during the density loop. See
 * theory/GEAR/Radiation/02_fuv_isrf.tex.
 */

#include "dimension.h"
#include "kernel_hydro.h"
#include "radiation.h"

#include <math.h>

/**
 * @brief Band contribution to particle i's `div(F)` accumulator, and the
 * mirrored, mass-weighted, opposite-sign contribution to j's.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr Particle i's kernel-gradient term, h_i^-(dim+1) * dW/dq.
 * @param wj_dr Particle j's kernel-gradient term, h_j^-(dim+1) * dW/dq.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param F_i Particle i's tracked reduced flux `Ft = F_true/c_hyp` (this band).
 * @param F_j Particle j's tracked reduced flux.
 * @param c_i Particle i's #feedback_part_data.c_hyp, the receiver-side speed
 * that turns `div(Ft)` into the accumulated `div(F)`.
 * @param c_j Particle j's #feedback_part_data.c_hyp, same role.
 * @param a_factor_comoving_to_physical `1/a`.
 * @param div_F_i (return, accumulated) Particle i's div(F) accumulator.
 * @param div_F_j (return, accumulated) Particle j's div(F) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_divergence_accumulate_band(const float dx[3], float r_inv,
                                     float wi_dr, float wj_dr, float mi,
                                     float mj, float rho_i, float rho_j,
                                     const float F_i[3], const float F_j[3],
                                     float c_i, float c_j,
                                     float a_factor_comoving_to_physical,
                                     float *div_F_i, float *div_F_j) {
  const float Fi_dot_dx = F_i[0] * dx[0] + F_i[1] * dx[1] + F_i[2] * dx[2];
  const float Fj_dot_dx = F_j[0] * dx[0] + F_j[1] * dx[1] + F_j[2] * dx[2];

  const float Phi_ij =
      (Fi_dot_dx / rho_i * wi_dr * r_inv + Fj_dot_dx / rho_j * wj_dr * r_inv) *
      a_factor_comoving_to_physical;

  *div_F_i += c_i * mj * Phi_ij;
  *div_F_j += -c_j * mi * Phi_ij;
}

/**
 * @brief Band contribution to the kernel-mean `|rho_prev*u_prev|` accumulator,
 * the field scale used by #radiation_update_dissipation_alpha_band.
 *
 * Built from the `u_prev` and `rho_prev` snapshots. It takes no
 * comoving-to-physical factor: its consumer divides it into `rho_prev*u`,
 * which has the same `a^3` weight.
 *
 * @param wi Particle i's kernel value, W(r/h_i)*h_i^-dim.
 * @param wj Particle j's kernel value, W(r/h_j)*h_j^-dim.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param u_i_prev Particle i's snapshotted specific field (this band).
 * @param u_j_prev Particle j's snapshotted specific field (this band).
 * @param ngb_mean_abs_u_V_i (return, accumulated) Particle i's accumulator.
 * @param ngb_mean_abs_u_V_j (return, accumulated) Particle j's accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_reference_accumulate_band(float wi, float wj, float mi,
                                                float mj, float rho_i,
                                                float rho_j, float u_i_prev,
                                                float u_j_prev,
                                                float *ngb_mean_abs_u_V_i,
                                                float *ngb_mean_abs_u_V_j) {

  *ngb_mean_abs_u_V_i += (mj / rho_j) * wi * fabsf(rho_j * u_j_prev);
  *ngb_mean_abs_u_V_j += (mi / rho_i) * wj * fabsf(rho_i * u_i_prev);
}

/**
 * @brief Band contribution to particle i's negativity-triggered
 * artificial-dissipation source term, and the mirrored contribution to j's.
 *
 * The coefficient is `alpha_ij = max(trigger_i, trigger_j, floor_i,
 * floor_j)`; each side uses its own `c_hyp`. The live `u` read here is `u^n`.
 * See theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Artificial dissipation".
 *
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param c_i Particle i's #feedback_part_data.c_hyp.
 * @param c_j Particle j's #feedback_part_data.c_hyp.
 * @param alpha_trigger_i Particle i's
 * #feedback_isrf_operator_data.dissipation_alpha_trigger (this band).
 * @param alpha_trigger_j Particle j's, same field.
 * @param alpha_floor_i Particle i's
 * #feedback_isrf_operator_data.dissipation_alpha_floor (this band).
 * @param alpha_floor_j Particle j's, same field.
 * @param u_i Particle i's live specific field `u^n` (this band).
 * @param u_j Particle j's live specific field `u^n` (this band).
 * @param a_factor_comoving_to_physical `1/a`.
 * @param dissipation_u_i (return, accumulated) Particle i's accumulator.
 * @param dissipation_u_j (return, accumulated) Particle j's accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_force_accumulate_band(
    float wi_dr, float wj_dr, float mi, float mj, float rho_i, float rho_j,
    float c_i, float c_j, float alpha_trigger_i, float alpha_trigger_j,
    float alpha_floor_i, float alpha_floor_j, float u_i, float u_j,
    float a_factor_comoving_to_physical, float *dissipation_u_i,
    float *dissipation_u_j) {

  const float d_ij = rho_i * u_i - rho_j * u_j;

  /* Not nested: max() is a statement expression, and nesting trips
   * -Wshadow. */
  const float alpha_trigger_ij = max(alpha_trigger_i, alpha_trigger_j);
  const float alpha_floor_ij = max(alpha_floor_i, alpha_floor_j);
  const float alpha_ij = max(alpha_trigger_ij, alpha_floor_ij);

  const float Wbar_ij = 0.5f * (wi_dr + wj_dr);

  const float shape_ij =
      d_ij * Wbar_ij / (rho_i * rho_j) * a_factor_comoving_to_physical;

  *dissipation_u_i += mj * (alpha_ij * c_i * shape_ij);
  *dissipation_u_j += -mi * (alpha_ij * c_j * shape_ij);
}

/**
 * @brief M1 closure coefficients for one particle, one band, from its own `(u,
 * F, c_M)`; `c_M` normalises the flux (1 for the reduced flux).
 *
 * `f = min(1, |F|/(c_M*u))` for `u > 0`, else 0;
 * `chi(f) = (3+4f^2)/(5+2*sqrt(4-3f^2))`. See
 * theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Moment equations and closure".
 *
 * `F2`, `|F|`, `c_M*u` and `f` are in double: `F.F` underflows in float32 for
 * `|F| < 1.1e-19` and would make a faint beam isotropic.
 *
 * @param u This band's specific field `u^n` for this particle.
 * @param F This particle's tracked reduced flux (this band).
 * @param c_M The speed that normalises `F` against `u`: 1 for the tracked
 * reduced flux.
 * @param n (return) The flux direction `F/|F|`, zero at `F = 0`.
 * @param iso_coeff (return) `(1-chi)/2`.
 * @param aniso_coeff (return) `(3chi-1)/2`.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_m1_closure_coefficients_band(float u, const float F[3], float c_M,
                                           float n[3], float *iso_coeff,
                                           float *aniso_coeff) {

  const double F2 = (double)F[0] * (double)F[0] + (double)F[1] * (double)F[1] +
                    (double)F[2] * (double)F[2];
  const double F_inv = (F2 > 0.) ? 1. / sqrt(F2) : 0.;
  const double Fmag = F2 * F_inv; /* sqrt(F2), no second sqrt call */
  n[0] = (float)(F[0] * F_inv);
  n[1] = (float)(F[1] * F_inv);
  n[2] = (float)(F[2] * F_inv);

  const double denom = (double)c_M * (double)u;
  const float f = (denom > 0.) ? (float)min(Fmag / denom, 1.) : 0.f;

  const float sq = 4.f - 3.f * f * f;
  const float chi = (3.f + 4.f * f * f) / (5.f + 2.f * sqrtf(sq));

  *iso_coeff = 0.5f * (1.f - chi);
  *aniso_coeff = 0.5f * (3.f * chi - 1.f);
}

/**
 * @brief Assemble the M1 closure tensor
 * `D = iso_coeff I + aniso_coeff (n dyadic n)`.
 *
 * @param n The flux direction.
 * @param iso_coeff `(1-chi)/2`.
 * @param aniso_coeff `(3chi-1)/2`.
 * @param D (return) The 3x3 closure tensor.
 */
__attribute__((always_inline)) INLINE static void
radiation_build_m1_closure_tensor(const float n[3], float iso_coeff,
                                  float aniso_coeff, float D[3][3]) {

  for (int a = 0; a < 3; a++) {
    D[a][0] = aniso_coeff * n[a] * n[0];
    D[a][1] = aniso_coeff * n[a] * n[1];
    D[a][2] = aniso_coeff * n[a] * n[2];
    D[a][a] += iso_coeff;
  }
}

/**
 * @brief M1 closure tensor `D(f)` for one particle, one band, computed from
 * scratch. See #radiation_get_m1_closure_coefficients_band.
 *
 * @param u This band's specific field `u^n` for this particle.
 * @param F This particle's tracked reduced flux (this band).
 * @param c_M The speed that normalises `F` against `u`: 1 for the tracked
 * reduced flux.
 * @param D (return) The 3x3 closure tensor.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_m1_closure_tensor_band(float u, const float F[3], float c_M,
                                     float D[3][3]) {

  float n[3], iso_coeff, aniso_coeff;
  radiation_get_m1_closure_coefficients_band(u, F, c_M, n, &iso_coeff,
                                             &aniso_coeff);
  radiation_build_m1_closure_tensor(n, iso_coeff, aniso_coeff, D);
}

/**
 * @brief Cache every band's M1 closure tensor on the particle.
 *
 * Must run after the last write of `u` and `specific_flux` before a gradient
 * loop reads the particle: at drift time and at first init (the initial pass
 * has no drift).
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
radiation_cache_m1_closure_part(struct part *p) {

  struct feedback_part_data *fd = &p->feedback_data;
  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    const int m = radiation_isrf_operator_owner[o];
    const struct feedback_isrf_moment_data *moment = &fd->isrf_moment[m];
    struct feedback_isrf_operator_data *op = &fd->isrf_operator[o];
    radiation_get_m1_closure_tensor_band(moment->u, moment->specific_flux, 1.f,
                                         op->m1_closure_D);
  }
}

/**
 * @brief Band contribution to both particles' `grad(u)` accumulators: the
 * anisotropic M1 pressure-tensor divergence `1/rho * div(D(f)*rho*u)`.
 *
 * Each particle uses its own kernel derivative and closure tensor, with no
 * grad-h factor. Paired with #radiation_divergence_accumulate_band this is
 * exactly skew-adjoint. See theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Spatial
 * operators".
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param u_i Particle i's specific field `u^n` (this band).
 * @param u_j Particle j's specific field `u^n` (this band).
 * @param D_i Particle i's own M1 closure tensor (this band), cached by
 * #radiation_cache_m1_closure_part.
 * @param D_j Particle j's own M1 closure tensor (this band), same source.
 * @param a_factor_comoving_to_physical `1/a`, folded into `fac_i` and `fac_j`.
 * @param grad_u_i (return, accumulated) Particle i's grad(u) accumulator.
 * @param grad_u_j (return, accumulated) Particle j's grad(u) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_gradient_accumulate_band(const float dx[3], float r_inv, float wi_dr,
                                   float wj_dr, float mi, float mj, float rho_i,
                                   float rho_j, float u_i, float u_j,
                                   const float D_i[3][3], const float D_j[3][3],
                                   float a_factor_comoving_to_physical,
                                   float grad_u_i[3], float grad_u_j[3]) {

  const float rho_i_inv = 1.f / rho_i;
  const float rho_j_inv = 1.f / rho_j;

  float temp_i[3], temp_j[3];
  for (int k = 0; k < 3; k++) {
    const float Di_dot_dx =
        D_i[k][0] * dx[0] + D_i[k][1] * dx[1] + D_i[k][2] * dx[2];
    const float Dj_dot_dx =
        D_j[k][0] * dx[0] + D_j[k][1] * dx[1] + D_j[k][2] * dx[2];
    temp_i[k] = Di_dot_dx * rho_i * u_i * r_inv;
    temp_j[k] = Dj_dot_dx * rho_j * u_j * r_inv;
  }

  const float fac_i =
      mj * rho_i_inv * rho_i_inv * wi_dr * a_factor_comoving_to_physical;
  const float fac_j =
      mi * rho_j_inv * rho_j_inv * wj_dr * a_factor_comoving_to_physical;

  for (int k = 0; k < 3; k++) {
    grad_u_i[k] += -(temp_i[k] - temp_j[k]) * fac_i;
    grad_u_j[k] += -(temp_i[k] - temp_j[k]) * fac_j;
  }
}

/**
 * @brief Density-loop propagation interaction between two particles
 * (symmetric): both particles' trigger reference accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor (unused, kept for interface conformance).
 * @param H Current Hubble parameter (unused, same reason).
 * @param us Unit system (unused, same reason).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_propagation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wj;
  kernel_eval(r * hi_inv, &wi);
  kernel_eval(r * hj_inv, &wj);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Slowest time bin in each particle's kernel. */
  fdi->max_ngb_time_bin = max(fdi->max_ngb_time_bin, pj->time_bin);
  fdj->max_ngb_time_bin = max(fdj->max_ngb_time_bin, pi->time_bin);

  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    const int m = radiation_isrf_operator_owner[o];
    const struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    const struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];
    struct feedback_isrf_operator_data *op_j = &fdj->isrf_operator[o];

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, moment_i->u_prev, moment_j->u_prev,
        &op_i->ngb_mean_abs_u_V, &op_j->ngb_mean_abs_u_V);
  }
}

/**
 * @brief Density-loop interaction between two particles (non-symmetric): only
 * particle i's trigger reference accumulator is updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor (unused, see #runner_iact_isrf_propagation).
 * @param H Current Hubble parameter (unused, same reason).
 * @param us Unit system (unused, same reason).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_propagation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H,
                                    const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wj;
  kernel_eval(r * hi_inv, &wi);
  kernel_eval(r * hj_inv, &wj);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  fdi->max_ngb_time_bin = max(fdi->max_ngb_time_bin, pj->time_bin);

  for (int o = 0; o < ISRF_OPERATOR_COUNT; o++) {
    const int m = radiation_isrf_operator_owner[o];
    const struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    const struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_ngb_mean_abs_u_V = 0.f;

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, moment_i->u_prev, moment_j->u_prev,
        &op_i->ngb_mean_abs_u_V, &unused_ngb_mean_abs_u_V);
  }
}

/**
 * @brief `grad(u)` interaction between two particles (symmetric): both
 * particles' accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const float a_factor_comoving_to_physical = 1.f / a;

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    const enum radiation_isrf_operator o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];
    const struct feedback_isrf_operator_data *op_j = &fdj->isrf_operator[o];

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, moment_i->u, moment_j->u,
        op_i->m1_closure_D, op_j->m1_closure_D, a_factor_comoving_to_physical,
        moment_i->grad_u, moment_j->grad_u);
  }
}

/**
 * @brief `grad(u)` interaction between two particles (non-symmetric):
 * only particle i's accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const float a_factor_comoving_to_physical = 1.f / a;

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    const struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    const enum radiation_isrf_operator o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];
    const struct feedback_isrf_operator_data *op_j = &fdj->isrf_operator[o];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_grad_u[3] = {0.f, 0.f, 0.f};

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, moment_i->u, moment_j->u,
        op_i->m1_closure_D, op_j->m1_closure_D, a_factor_comoving_to_physical,
        moment_i->grad_u, unused_grad_u);
  }
}

/**
 * @brief One band of a pair whose members have unequal steps: the finer member
 * may take its rate, the coarser one gets its share as pending amounts.
 *
 * The coarser share is integrated over the finer step and stored divided by
 * c_hyp, so the pair is booked once, from the finer timeline. The amount
 * carries no relaxation factor, so the finer member must have a = 0.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param moment_i Particle i's moment (this band).
 * @param moment_j Particle j's moment (this band).
 * @param op_i Particle i's operator of this band.
 * @param op_j Particle j's operator of this band.
 * @param c_i Particle i's #feedback_part_data.c_hyp.
 * @param c_j Particle j's #feedback_part_data.c_hyp.
 * @param a_factor_comoving_to_physical `1/a`.
 * @param H Current Hubble parameter.
 * @param i_is_fine 1 if particle i has the shorter step.
 * @param fine_rate 1 to add the finer member's rate here.
 * @param dt_fine The finer member's step.
 */
__attribute__((always_inline)) INLINE static void radiation_cross_bin_pair_band(
    const float dx[3], float r_inv, float wi_dr, float wj_dr, float mi,
    float mj, float rho_i, float rho_j,
    struct feedback_isrf_moment_data *moment_i,
    struct feedback_isrf_moment_data *moment_j,
    const struct feedback_isrf_operator_data *op_i,
    const struct feedback_isrf_operator_data *op_j, float c_i, float c_j,
    float a_factor_comoving_to_physical, float H, int i_is_fine, int fine_rate,
    float dt_fine) {

  const struct feedback_isrf_operator_data *op_fine = i_is_fine ? op_i : op_j;
  if (op_fine->kappa != 0.f || H != 0.f)
    error(
        "The cross-bin pending deposit needs a = 0 (no dust, no expansion): "
        "the finer member has kappa %e, H %e.",
        op_fine->kappa, H);

  /* The coarser side is formed with c_hyp = 1, which makes it c-free. */
  const float cf_i = i_is_fine ? c_i : 1.f;
  const float cf_j = i_is_fine ? 1.f : c_j;
  float div[2] = {0.f, 0.f};
  float diss[2] = {0.f, 0.f};
  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, moment_i->specific_flux,
      moment_j->specific_flux, cf_i, cf_j, a_factor_comoving_to_physical,
      &div[0], &div[1]);
  radiation_dissipation_force_accumulate_band(
      wi_dr, wj_dr, mi, mj, rho_i, rho_j, cf_i, cf_j,
      op_i->dissipation_alpha_trigger, op_j->dissipation_alpha_trigger,
      op_i->dissipation_alpha_floor, op_j->dissipation_alpha_floor, moment_i->u,
      moment_j->u, a_factor_comoving_to_physical, &diss[0], &diss[1]);

  const int f = i_is_fine ? 0 : 1;
  struct feedback_isrf_moment_data *fine = i_is_fine ? moment_i : moment_j;
  struct feedback_isrf_moment_data *coarse = i_is_fine ? moment_j : moment_i;
  if (fine_rate) {
    fine->div_specific_flux += div[f];
    fine->dissipation_u += diss[f];
  }
  /* Gains to u over the finer step: -div and +diss. */
  coarse->pending_transport_u += -div[1 - f] * dt_fine;
  coarse->pending_dissipation_u += diss[1 - f] * dt_fine;
}

/**
 * @brief All bands of a pair whose members have unequal steps, see
 * #radiation_cross_bin_pair_band. Out of line, so the same-step path keeps
 * its code generation.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param fdi Particle i's feedback data.
 * @param fdj Particle j's feedback data.
 * @param a_factor_comoving_to_physical `1/a`.
 * @param H Current Hubble parameter.
 * @param i_is_fine 1 if particle i has the shorter step.
 * @param fine_rate 1 to add the finer member's rate here.
 * @param dt_fine The finer member's step.
 */
__attribute__((noinline)) static void radiation_cross_bin_pair(
    const float dx[3], float r_inv, float wi_dr, float wj_dr, float mi,
    float mj, struct feedback_part_data *fdi, struct feedback_part_data *fdj,
    float a_factor_comoving_to_physical, float H, int i_is_fine, int fine_rate,
    float dt_fine) {

  const float c_i = fdi->c_hyp;
  const float c_j = fdj->c_hyp;
  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    const enum radiation_isrf_operator o = radiation_isrf_moment_to_operator[m];
    radiation_cross_bin_pair_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, fdi->rho_prev, fdj->rho_prev,
        &fdi->isrf_moment[m], &fdj->isrf_moment[m], &fdi->isrf_operator[o],
        &fdj->isrf_operator[o], c_i, c_j, a_factor_comoving_to_physical, H,
        i_is_fine, fine_rate, dt_fine);
  }
}

/**
 * @brief Force-loop interaction between two particles (symmetric): both
 * particles' `div(F)` and dissipation terms are updated.
 *
 * Runs after the extra ghost has relaxed `specific_flux` and set `alpha`.
 * With unequal steps, see #radiation_cross_bin_pair.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_dissipation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);
  /* Read outside the band loop: the band writes may alias c_hyp. */
  const float c_i = fdi->c_hyp;
  const float c_j = fdj->c_hyp;

  const float a_factor_comoving_to_physical = 1.f / a;

  const float s_i = fdi->dt_active;
  const float s_j = fdj->dt_active;
  if (s_i != s_j && s_i != 0.f && s_j != 0.f) {
#ifdef SWIFT_DEBUG_CHECKS
    if (s_i < 0.f || s_j < 0.f)
      error("Symmetric ISRF pair with an inactive member (%e, %e).", s_i, s_j);
#endif
    radiation_cross_bin_pair(dx, r_inv, wi_dr, wj_dr, mi, mj, fdi, fdj,
                             a_factor_comoving_to_physical, H, s_i < s_j,
                             /*fine_rate=*/1, min(s_i, s_j));
    return;
  }

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    const enum radiation_isrf_operator o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];
    const struct feedback_isrf_operator_data *op_j = &fdj->isrf_operator[o];

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, moment_i->specific_flux,
        moment_j->specific_flux, c_i, c_j, a_factor_comoving_to_physical,
        &moment_i->div_specific_flux, &moment_j->div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, c_i, c_j,
        op_i->dissipation_alpha_trigger, op_j->dissipation_alpha_trigger,
        op_i->dissipation_alpha_floor, op_j->dissipation_alpha_floor,
        moment_i->u, moment_j->u, a_factor_comoving_to_physical,
        &moment_i->dissipation_u, &moment_j->dissipation_u);
  }
}

/**
 * @brief Force-loop interaction between two particles (non-symmetric): only
 * particle i's `div(F)` and dissipation terms are updated.
 *
 * With unequal steps, an inactive j receives its pending share, and an i
 * coarser than an active j books its own share as pending.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (only its pending amounts may be updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_dissipation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    struct part *restrict pj, const float a,
                                    const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);
  const float c_i = fdi->c_hyp;
  const float c_j = fdj->c_hyp;

  const float a_factor_comoving_to_physical = 1.f / a;

  /* An active coarser j books this pair in its own visit. */
  const float s_i = fdi->dt_active;
  const float s_j = fdj->dt_active;
  if (s_i != s_j && s_i != 0.f && s_j != 0.f && s_j < s_i) {
#ifdef SWIFT_DEBUG_CHECKS
    if (s_i < 0.f)
      error("Inactive particle %lld visits its neighbours.", pi->id);
    if (s_j < 0.f && !(fdj->dt_prev > fdi->dt_prev))
      error("Inactive particle %lld is not coarser than %lld.", pj->id, pi->id);
#endif
    /* s_j < 0: j inactive, i finer. 0 < s_j < s_i: i coarser. */
    const int i_is_fine = s_j < 0.f;
    radiation_cross_bin_pair(dx, r_inv, wi_dr, wj_dr, mi, mj, fdi,
                             &pj->feedback_data, a_factor_comoving_to_physical,
                             H, i_is_fine,
                             /*fine_rate=*/i_is_fine, i_is_fine ? s_i : s_j);
    return;
  }

  for (int m = 0; m < ISRF_MOMENT_COUNT; m++) {
    struct feedback_isrf_moment_data *moment_i = &fdi->isrf_moment[m];
    const struct feedback_isrf_moment_data *moment_j = &fdj->isrf_moment[m];
    const enum radiation_isrf_operator o = radiation_isrf_moment_to_operator[m];
    const struct feedback_isrf_operator_data *op_i = &fdi->isrf_operator[o];
    const struct feedback_isrf_operator_data *op_j = &fdj->isrf_operator[o];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_div_specific_flux = 0.f;
    float unused_dissipation_u = 0.f;

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, moment_i->specific_flux,
        moment_j->specific_flux, c_i, c_j, a_factor_comoving_to_physical,
        &moment_i->div_specific_flux, &unused_div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, c_i, c_j,
        op_i->dissipation_alpha_trigger, op_j->dissipation_alpha_trigger,
        op_i->dissipation_alpha_floor, op_j->dissipation_alpha_floor,
        moment_i->u, moment_j->u, a_factor_comoving_to_physical,
        &moment_i->dissipation_u, &unused_dissipation_u);
  }
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
