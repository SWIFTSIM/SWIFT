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
 * @brief Gas-gas density-loop, gradient-loop and force-loop hooks for the
 * hyperbolic M1-relaxation propagation of the per-band u and specific_flux
 * fields.
 *
 * `div(F)` sits in the force loop rather than a type-1 loop (density,
 * gradient), which reaches particle i only for r < H_i: it needs both sides of
 * every pair credited, and only the force loop's type-2 dispatch fires both
 * sides whenever either kernel reaches. Derivations are in
 * theory/GEAR/Radiation/02_fuv_isrf.tex, secs. "Spatial operators" and "The
 * consistent variable-speed operators".
 *
 * Nothing writes `u` between the drift snapshot and the end-force ghost, so
 * the live `u` the gradient and dissipation loops read is `u^n`.
 *
 * The inputs are comoving and the accumulators physical. Each spatial operator
 * carries one net inverse length, so each closes with
 * `a_factor_comoving_to_physical = 1/a`.
 *
 * All three read `feedback_data.rho_prev`, not `p->rho`: the density loop runs
 * while `p->rho` is still a partial sum, and the skew-adjoint pairing of the
 * divergence and the gradient needs both built from the same densities.
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
 * @param F_i Particle i's tracked flux (this band), or the reduced flux
 * `Ft = F_true/c_hyp` under #isrf_c_hyp_consistent_variable_c.
 * @param F_j Particle j's tracked flux, same convention.
 * @param c_i Particle i's #feedback_part_data.c_hyp, used only under
 * #isrf_c_hyp_consistent_variable_c.
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
  /* Both schemes share the coefficient Phi_ij and differ only in a trailing
   * scalar. Reassociation is disabled so they stay bit-for-bit identical at
   * uniform c_hyp (tests/testRadiationISRFForceDispatchConservation.c). This
   * holds under clang only: GCC has no block-scoped equivalent. */
  {
#if defined(__clang__)
#pragma clang fp reassociate(off) contract(off) reciprocal(off)
#endif
    const float Fi_dot_dx = F_i[0] * dx[0] + F_i[1] * dx[1] + F_i[2] * dx[2];
    const float Fj_dot_dx = F_j[0] * dx[0] + F_j[1] * dx[1] + F_j[2] * dx[2];

    const float Phi_ij = (Fi_dot_dx / rho_i * wi_dr * r_inv +
                          Fj_dot_dx / rho_j * wj_dr * r_inv) *
                         a_factor_comoving_to_physical;

    if (isrf_c_hyp_consistent_variable_c) {
      /* F_i and F_j already hold the reduced flux, so Phi_ij is already
       * div(Ft)'s shared coefficient and each side takes its OWN (receiver)
       * c_hyp. The pair then conserves sum m_i X_i / c_hyp_i, not
       * sum m_i X_i. */
      *div_F_i += c_i * mj * Phi_ij;
      *div_F_j += -c_j * mi * Phi_ij;
      return;
    }

    *div_F_i += mj * Phi_ij;
    *div_F_j += -mi * Phi_ij;
  }
}

/**
 * @brief Band contribution to each particle's kernel-mean `|rho_prev*u_prev|`
 * reference accumulator, the local field scale the negativity trigger divides
 * an undershoot by (#radiation_update_dissipation_alpha_band).
 *
 * Built from the `u_prev` and `rho_prev` snapshots, so the value does not
 * drift across a particle's h-iterations. Takes no comoving-to-physical
 * conversion: its only consumer divides it into `rho_prev*u`, which carries
 * the same `a^3` weight, so converting here would introduce a bias.
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
 * `v_sig,ij = alpha_ij * min(c_hyp_i, c_hyp_j)` with
 * `alpha_ij = max(trigger_i, trigger_j, floor_i, floor_j)`; see
 * theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Artificial dissipation". Runs
 * after the extra ghost has set this step's `alpha`, on the live `u`, which is
 * still `u^n` there. No mutual-reach gate is needed: the force loop fires both
 * sides whenever either kernel reaches.
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
 * #feedback_isrf_band_data.dissipation_alpha_trigger (this band).
 * @param alpha_trigger_j Particle j's, same field.
 * @param alpha_floor_i Particle i's
 * #feedback_isrf_band_data.dissipation_alpha_floor (this band).
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

  /* Same fp-reassociation hazard and clang-only guard as
   * #radiation_divergence_accumulate_band's pragma. */
  {
#if defined(__clang__)
#pragma clang fp reassociate(off) contract(off) reciprocal(off)
#endif
    const float d_ij = rho_i * u_i - rho_j * u_j;

    /* Split out rather than nested: max() expands to a statement
     * expression with its own locals, which -Wshadow rejects when
     * nested. */
    const float alpha_trigger_ij = max(alpha_trigger_i, alpha_trigger_j);
    const float alpha_floor_ij = max(alpha_floor_i, alpha_floor_j);
    const float alpha_ij = max(alpha_trigger_ij, alpha_floor_ij);

    const float Wbar_ij = 0.5f * (wi_dr + wj_dr);

    const float shape_ij =
        d_ij * Wbar_ij / (rho_i * rho_j) * a_factor_comoving_to_physical;

    if (isrf_c_hyp_consistent_variable_c) {
      /* Two receiver-side speeds, not one shared minimum. The pair then
       * conserves sum m_i X_i / c_hyp_i, not sum m_i X_i. */
      *dissipation_u_i += mj * (alpha_ij * c_i * shape_ij);
      *dissipation_u_j += -mi * (alpha_ij * c_j * shape_ij);
      return;
    }

    const float Psi_ij = alpha_ij * min(c_i, c_j) * shape_ij;

    *dissipation_u_i += mj * Psi_ij;
    *dissipation_u_j += -mi * Psi_ij;
  }
}

/**
 * @brief M1 closure coefficients for one particle, one band, built from its
 * own `(u, F, c_M)`. `c_M` is #feedback_part_data.c_hyp, reinterpreted as the
 * fastest M1 characteristic (`f=1`), not a new field.
 *
 * `f = min(1, |F|/(c_M*u))` for `u > 0`, `f = 0` otherwise;
 * `chi(f) = (3+4f^2)/(5+2*sqrt(4-3f^2))`;
 * `D(f) = (1-chi)/2 I + (3chi-1)/2 (n dyadic n)`, `n = F/|F|`. See
 * theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Moment equations and closure".
 *
 * `F2`, `|F|`, `c_M*u` and `f` are formed in double because in float32 `F.F`
 * underflows to zero once `|F| < sqrt(FLT_MIN) ~ 1.1e-19` (internal units),
 * which would turn a faint beam into an isotropic closure. The zero guards
 * keep `F = 0` at `n = 0`, `f = 0` rather than a NaN.
 *
 * @param u This band's specific field `u^n` for this particle.
 * @param F This particle's tracked flux (this band).
 * @param c_M This particle's own #feedback_part_data.c_hyp.
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
 * @param F This particle's tracked flux (this band).
 * @param c_M This particle's own #feedback_part_data.c_hyp.
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
 * Must run after the last write of `u`, `specific_flux` and `c_hyp` preceding
 * a gradient loop that reads the particle. Three call sites cover this: the
 * drift-time reset, the density ghost once this step's `c_hyp` is known, and
 * first init, since the initial ti = 0 pass reaches the gradient loop without
 * a drift.
 *
 * Under #isrf_c_hyp_consistent_variable_c the stored flux is already reduced
 * by `c_hyp`, so the closure is built with `c_M = 1` and yields the same `f`.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
radiation_cache_m1_closure_part(struct part *p) {

  struct feedback_part_data *fd = &p->feedback_data;
  const float c_M = isrf_c_hyp_consistent_variable_c ? 1.f : fd->c_hyp;
  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *band = &fd->isrf_band[b];
    radiation_get_m1_closure_tensor_band(band->u, band->specific_flux, c_M,
                                         band->m1_closure_D);
  }
}

/**
 * @brief Band contribution to both particles' `grad(u)` accumulators: the
 * anisotropic M1 pressure-tensor divergence `1/rho * div(D(f)*rho*u)`.
 *
 * Each particle uses its own kernel derivative and own closure tensor, with no
 * shared average and no grad-h `forcef` factor. That is the deliberate
 * complement of #radiation_divergence_accumulate_band's shared-coefficient
 * construction, and pairing the two is what makes them exactly skew-adjoint.
 * See theory/GEAR/Radiation/02_fuv_isrf.tex sec. "Spatial operators".
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

  /* Slowest clock in each particle's own kernel, read by
   * radiation_end_density_propagation once the h-iteration converges. */
  fdi->max_ngb_time_bin = max(fdi->max_ngb_time_bin, pj->time_bin);
  fdj->max_ngb_time_bin = max(fdj->max_ngb_time_bin, pi->time_bin);

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, bi->u_prev, bj->u_prev,
        &bi->ngb_mean_abs_u_V, &bj->ngb_mean_abs_u_V);
  }
}

/**
 * @brief Density-loop propagation interaction between two particles
 * (non-symmetric): only particle i's trigger reference accumulator is
 * updated.
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

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_ngb_mean_abs_u_V = 0.f;

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, bi->u_prev, bj->u_prev,
        &bi->ngb_mean_abs_u_V, &unused_ngb_mean_abs_u_V);
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

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->u, bj->u,
        bi->m1_closure_D, bj->m1_closure_D, a_factor_comoving_to_physical,
        bi->grad_u, bj->grad_u);
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

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_grad_u[3] = {0.f, 0.f, 0.f};

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->u, bj->u,
        bi->m1_closure_D, bj->m1_closure_D, a_factor_comoving_to_physical,
        bi->grad_u, unused_grad_u);
  }
}

/**
 * @brief Force-loop propagation interaction between two particles
 * (symmetric): both particles' `div(F)` and dissipation accumulators are
 * updated.
 *
 * Runs after the extra ghost has relaxed this step's `specific_flux` and set
 * the dissipation `alpha`. The dispatch fires both sides of a pair whenever
 * either kernel reaches, which keeps both mirrored pairs whole at h_i != h_j.
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
  /* Read outside the band loop: the band writes may alias c_hyp for the
   * compiler. Both speeds are kept, not just their minimum, because the
   * consistent-variable-c scheme needs each side separately. */
  const float c_i = fdi->c_hyp;
  const float c_j = fdj->c_hyp;

  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->specific_flux,
        bj->specific_flux, c_i, c_j, a_factor_comoving_to_physical,
        &bi->div_specific_flux, &bj->div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, c_i, c_j,
        bi->dissipation_alpha_trigger, bj->dissipation_alpha_trigger,
        bi->dissipation_alpha_floor, bj->dissipation_alpha_floor, bi->u, bj->u,
        a_factor_comoving_to_physical, &bi->dissipation_u, &bj->dissipation_u);
  }
}

/**
 * @brief Force-loop propagation interaction between two particles
 * (non-symmetric): only particle i's `div(F)` and dissipation accumulators
 * are updated.
 *
 * Reached once per side, so a pair dispatched on both sides still receives
 * both mirrored pairs in full; see #runner_iact_isrf_dissipation.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulator not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_dissipation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H) {

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
  /* Read outside the band loop, see #runner_iact_isrf_dissipation. */
  const float c_i = fdi->c_hyp;
  const float c_j = fdj->c_hyp;

  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j is not written here, so discard its side of the pair. */
    float unused_div_specific_flux = 0.f;
    float unused_dissipation_u = 0.f;

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->specific_flux,
        bj->specific_flux, c_i, c_j, a_factor_comoving_to_physical,
        &bi->div_specific_flux, &unused_div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, c_i, c_j,
        bi->dissipation_alpha_trigger, bj->dissipation_alpha_trigger,
        bi->dissipation_alpha_floor, bj->dissipation_alpha_floor, bi->u, bj->u,
        a_factor_comoving_to_physical, &bi->dissipation_u,
        &unused_dissipation_u);
  }
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
