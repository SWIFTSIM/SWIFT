#ifndef SWIFT_FVPM_GEOMETRY_GIZMO_H
#define SWIFT_FVPM_GEOMETRY_GIZMO_H

#include "const.h"
#include "part.h"

#include <config.h>

/**
 * @file Gizmo/fvpm_geometry.h
 * @brief Functions related to the Gizmo FVPM geometry struct collection,
 * in particular the collection of the data required for the matrix needed
 * for gradients.
 * This was moved here so we can cleanly couple GEAR-RT on top of SPH
 * hydrodynamics while avoiding code replication.
 */

#if defined(RT_GEAR) && defined(GIZMO_MFM_SPH)
/* Some functions clash here. MFM resets and does some geometry centroid
 * stuff, while GEAR-RT, which uses MFV, doesn't. So we'd need to split the
 * functions for RT and for hydro use.
 * However, it is very unlikely we'll ever actually use that combination,
 * so leaving it as-is for now. */
#error "Combining GIZMO MFM and GEAR-RT not implemented yet."
#endif

#if defined(GIZMO_MFV_SPH) || defined(RT_GEAR)
#include "./MFV/fvpm_geometry.h"
#elif defined(GIZMO_MFM_SPH)
#include "./MFM/fvpm_geometry.h"
#endif

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
fvpm_part_geometry_well_behaved(const struct part *restrict p) {

  return p->geometry.wcorr > const_gizmo_min_wcorr;
}

/**
 * @brief Collect the data needed for the matrix construction.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_geometry_and_matrix(struct part *restrict pi, const float wi,
                                    const float dx[3]) {
  /* these are eqns. (1) and (2) in the Gizmo theory summary */
  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;
}

/* Role-aware: the shared accumulator carries no left/right information, so
   the caller supplies the sign, exactly as fvpm_update_centroid_left/_right
   do. dx = x_i - x_j, and m_i accumulates d_ij = -dx. */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_first_moment_left(struct part *restrict p, const float dx[3],
                                  const float w) {
  p->geometry.first_moment[0] -= dx[0] * w;
  p->geometry.first_moment[1] -= dx[1] * w;
  p->geometry.first_moment[2] -= dx[2] * w;
}

__attribute__((always_inline)) INLINE static void
fvpm_accumulate_first_moment_right(struct part *restrict p, const float dx[3],
                                   const float w) {
  p->geometry.first_moment[0] += dx[0] * w;
  p->geometry.first_moment[1] += dx[1] * w;
  p->geometry.first_moment[2] += dx[2] * w;
}

__attribute__((always_inline)) INLINE static void fvpm_geometry_init(
    struct part *restrict p) {

  p->geometry.volume = 0.0f;
  p->geometry.matrix_E[0][0] = 0.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 0.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 0.0f;

  p->geometry.first_moment[0] = 0.0f;
  p->geometry.first_moment[1] = 0.0f;
  p->geometry.first_moment[2] = 0.0f;

  /* reset the centroid variables used for the velocity correction in MFV */
  fvpm_reset_centroids(p);

  if (p->time_bin == 0) {
    p->geometry.is_problematic = 0;
  }

  if (p->time_bin == 0) {
    /* This only works for the idealised examples... */
    p->geometry.area = 0.0;
    p->geometry.area_sum[0] = 0.0;
    p->geometry.area_sum[1] = 0.0;
    p->geometry.area_sum[2] = 0.0;

    p->geometry.area_sum_plus[0] = 0.0;
    p->geometry.area_sum_plus[1] = 0.0;
    p->geometry.area_sum_plus[2] = 0.0;

    p->geometry.area_sum_minus[0] = 0.0;
    p->geometry.area_sum_minus[1] = 0.0;
    p->geometry.area_sum_minus[2] = 0.0;

    p->geometry.area_sum1[0] = 0.0;
    p->geometry.area_sum1[1] = 0.0;
    p->geometry.area_sum1[2] = 0.0;

    p->geometry.area_sum2[0] = 0.0;
    p->geometry.area_sum2[1] = 0.0;
    p->geometry.area_sum2[2] = 0.0;
  }
}

/**
 * @brief Sets the geometry fields to sensible values when #part has 0 ngbs.
 *
 * @param p the particle to work on
 */
__attribute__((always_inline)) INLINE static void
fvpm_geometry_part_has_no_neighbours(struct part *restrict p) {

  /* Re-set problematic values */
  p->geometry.volume = 1.0f;
  p->geometry.matrix_E[0][0] = 1.0f;
  p->geometry.matrix_E[0][1] = 0.0f;
  p->geometry.matrix_E[0][2] = 0.0f;
  p->geometry.matrix_E[1][0] = 0.0f;
  p->geometry.matrix_E[1][1] = 1.0f;
  p->geometry.matrix_E[1][2] = 0.0f;
  p->geometry.matrix_E[2][0] = 0.0f;
  p->geometry.matrix_E[2][1] = 0.0f;
  p->geometry.matrix_E[2][2] = 1.0f;
  p->geometry.condition_number = 1.f;

  /* No neighbours means no centroid to centre on: fall back to the same
     identity matrix_E gets, so any consumer reading matrix_E_centred still
     gets a sane, finite matrix. */
  p->geometry.first_moment[0] = 0.0f;
  p->geometry.first_moment[1] = 0.0f;
  p->geometry.first_moment[2] = 0.0f;
  p->geometry.centroid_offset[0] = 0.0f;
  p->geometry.centroid_offset[1] = 0.0f;
  p->geometry.centroid_offset[2] = 0.0f;
  p->geometry.matrix_E_centred[0][0] = 1.0f;
  p->geometry.matrix_E_centred[0][1] = 0.0f;
  p->geometry.matrix_E_centred[0][2] = 0.0f;
  p->geometry.matrix_E_centred[1][0] = 0.0f;
  p->geometry.matrix_E_centred[1][1] = 1.0f;
  p->geometry.matrix_E_centred[1][2] = 0.0f;
  p->geometry.matrix_E_centred[2][0] = 0.0f;
  p->geometry.matrix_E_centred[2][1] = 0.0f;
  p->geometry.matrix_E_centred[2][2] = 1.0f;
  p->geometry.centring_margin = 0.0f;
  p->geometry.centring_disabled = 1;

  /* reset the centroid variables used for the velocity correction in MFV */
  fvpm_reset_centroids(p);

  /* TODO: To be defined */
  p->geometry.is_problematic = 0;  
  p->geometry.area = 0.0;
  p->geometry.area_sum[0] = 0.0;
  p->geometry.area_sum[1] = 0.0;
  p->geometry.area_sum[2] = 0.0;
  p->geometry.area_sum_plus[0] = 0.0;
  p->geometry.area_sum_plus[1] = 0.0;
  p->geometry.area_sum_plus[2] = 0.0;
  p->geometry.area_sum_minus[0] = 0.0;
  p->geometry.area_sum_minus[1] = 0.0;
  p->geometry.area_sum_minus[2] = 0.0;

  p->geometry.area_sum1[0] = 0.0;
  p->geometry.area_sum1[1] = 0.0;
  p->geometry.area_sum1[2] = 0.0;

  p->geometry.area_sum2[0] = 0.0;
  p->geometry.area_sum2[1] = 0.0;
  p->geometry.area_sum2[2] = 0.0;
}

/**
 * @brief Finish the computation of the matrix.
 *
 * @param p the particle to work on
 * @param ihdim 1/h^{dim}
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_volume_and_matrix(struct part *restrict p, const float ihdim) {

  /* Step 1: capture omega' = Sum_{j!=i} w_ij, raw, before geometry.volume is
     overwritten below. This is candidate A's kernel-sum denominator. */
  const float omega_prime = p->geometry.volume;

  /* Step 2 (unchanged): Final operation on the geometry. We multiply with the
   * smoothing kernel normalization ih3 and calculate the volume. */
  const float volume_inv = ihdim * (p->geometry.volume + kernel_root);
  const float volume = 1.0f / volume_inv;
  p->geometry.volume = volume;

  /* Step 3: compute the kernel-centred offset c and the centred matrix E^c
     from the raw matrix_E and first_moment, before the ihdim scaling below.
     centred is left 0 (no centring) if omega' is too small relative to
     kernel_root: a denormal omega' gives an enormous |c| and a downdate that
     annihilates E in float. */
  int centred = 0;
  float c[3] = {0.f, 0.f, 0.f};
  float Ec[3][3] = {{0.f}};

  /* -1 marks "SPD margin never evaluated" (omega' floor or step 8 tripped first); any other value is the computed q, disabled or not. */
  p->geometry.centring_margin = -1.0f;

  if (omega_prime > const_fvpm_min_omega_prime * kernel_root) {
    const float omega_inv = 1.0f / omega_prime;
    for (int k = 0; k < 3; k++) c[k] = p->geometry.first_moment[k] * omega_inv;
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        Ec[k][l] = (p->geometry.matrix_E[k][l] - omega_prime * c[k] * c[l]) * ihdim;
    centred = 1;
  }

  /* Step 4 (unchanged): we multiply with the smoothing kernel normalization */
  p->geometry.matrix_E[0][0] *= ihdim;
  p->geometry.matrix_E[0][1] *= ihdim;
  p->geometry.matrix_E[0][2] *= ihdim;
  p->geometry.matrix_E[1][0] *= ihdim;
  p->geometry.matrix_E[1][1] *= ihdim;
  p->geometry.matrix_E[1][2] *= ihdim;
  p->geometry.matrix_E[2][0] *= ihdim;
  p->geometry.matrix_E[2][1] *= ihdim;
  p->geometry.matrix_E[2][2] *= ihdim;

  /* Step 5 (unchanged position): normalise the centroids for MFV. */
  fvpm_normalise_centroid(p, p->density.wcount);

  /* Step 6 (unchanged): Check the condition number to see if we have a stable
     geometry. */
  const float condition_number_E =
      p->geometry.matrix_E[0][0] * p->geometry.matrix_E[0][0] +
      p->geometry.matrix_E[0][1] * p->geometry.matrix_E[0][1] +
      p->geometry.matrix_E[0][2] * p->geometry.matrix_E[0][2] +
      p->geometry.matrix_E[1][0] * p->geometry.matrix_E[1][0] +
      p->geometry.matrix_E[1][1] * p->geometry.matrix_E[1][1] +
      p->geometry.matrix_E[1][2] * p->geometry.matrix_E[1][2] +
      p->geometry.matrix_E[2][0] * p->geometry.matrix_E[2][0] +
      p->geometry.matrix_E[2][1] * p->geometry.matrix_E[2][1] +
      p->geometry.matrix_E[2][2] * p->geometry.matrix_E[2][2];

  p->geometry.condition_number = 0.0f;
  const int inv_status = invert_dimension_by_dimension_matrix(p->geometry.matrix_E);
  if (inv_status != 0) {
    /* something went wrong in the inversion; force bad condition number */
    p->geometry.condition_number = const_gizmo_max_condition_number + 1.0f;
  } else {
    const float condition_number_Einv =
        p->geometry.matrix_E[0][0] * p->geometry.matrix_E[0][0] +
        p->geometry.matrix_E[0][1] * p->geometry.matrix_E[0][1] +
        p->geometry.matrix_E[0][2] * p->geometry.matrix_E[0][2] +
        p->geometry.matrix_E[1][0] * p->geometry.matrix_E[1][0] +
        p->geometry.matrix_E[1][1] * p->geometry.matrix_E[1][1] +
        p->geometry.matrix_E[1][2] * p->geometry.matrix_E[1][2] +
        p->geometry.matrix_E[2][0] * p->geometry.matrix_E[2][0] +
        p->geometry.matrix_E[2][1] * p->geometry.matrix_E[2][1] +
        p->geometry.matrix_E[2][2] * p->geometry.matrix_E[2][2];

    p->geometry.condition_number =
        hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);
  }

  /* Step 8: a failed production inversion leaves matrix_E holding partially
     eliminated garbage, which q (step 9) would otherwise read as B. */
  if (inv_status != 0) centred = 0;

  /* Step 9: the SPD margin. Candidate A's guarantee that E^c is SPD whenever
     E is holds in exact arithmetic; in float, near-coincident neighbours can
     flip a face's orientation, so the Loewner condition is checked
     explicitly rather than trusted. */
  if (centred) {
    /* B = matrix_E after step 7 = E_hat^{-1}, with E_hat = ihdim.E_raw.
       E_hat^c = ihdim.(E_raw - omega' c(x)c) = E_hat - (ihdim.omega') c(x)c, so
       the Loewner condition for E_hat^c > 0 is exactly
       (ihdim.omega').c^T E_hat^{-1} c < 1. */
    const float (*B)[3] = (const float (*)[3])p->geometry.matrix_E;
    float Bc[3];
    for (int k = 0; k < 3; k++)
      Bc[k] = B[k][0] * c[0] + B[k][1] * c[1] + B[k][2] * c[2];
    const float q = ihdim * omega_prime * (c[0] * Bc[0] + c[1] * Bc[1] + c[2] * Bc[2]);
    p->geometry.centring_margin = q;

    /* Negated form, so NaN and inf fail closed. */
    if (!(q < 1.0f - const_fvpm_centring_spd_margin)) centred = 0;
  }

  /* Step 10: under separate storage there are two inversions; a failed
     centred inversion leaves Ec holding the same kind of garbage. */
  if (centred) {
    if (invert_dimension_by_dimension_matrix(Ec) != 0) centred = 0;
  }

  /* Step 11: store, with a defined fallback so consumers need only one code
     path. With centroid_offset = 0 and matrix_E_centred = matrix_E, a
     centring-disabled particle's centred face is bit-identical to its
     uncentred face. */
  if (centred) {
    for (int k = 0; k < 3; k++) {
      p->geometry.centroid_offset[k] = c[k];
      for (int l = 0; l < 3; l++) p->geometry.matrix_E_centred[k][l] = Ec[k][l];
    }
    p->geometry.centring_disabled = 0;
  } else {
    for (int k = 0; k < 3; k++) {
      p->geometry.centroid_offset[k] = 0.0f;
      for (int l = 0; l < 3; l++)
        p->geometry.matrix_E_centred[k][l] = p->geometry.matrix_E[k][l];
    }
    p->geometry.centring_disabled = 1;
  }

  /* Step 12 (unchanged): reads only condition_number, computed above from the
     untouched matrix_E, so wcorr and the fallback decision are unaffected by
     centring. */
  if (p->geometry.condition_number > const_gizmo_max_condition_number &&
      p->geometry.wcorr > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
          const_gizmo_max_condition_number, p->geometry.condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            p->geometry.condition_number, const_gizmo_max_condition_number,
            p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    p->geometry.wcorr = const_gizmo_w_correction_factor * p->geometry.wcorr;
  }
}

/**
 * @brief Compute the face area between i and j.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param Bi Matrix B for particle i.
 * @param Bj Matrix B for particle j.
 * @param ci Kernel centroid offset for particle i (0 for no centring).
 * @param cj Kernel centroid offset for particle j (0 for no centring).
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param (return) A The face area between i and j.
 * @param (return) A1 The particle-i-only term of A (row-sum term).
 * @param (return) A2 The particle-j-only term of A (column-sum term).
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_face_area_vector(const struct part *restrict pi, const struct part *restrict pj,
				 float Bi[3][3], float Bj[3][3],
				 const float ci[3], const float cj[3], const float r2,
				 const float dx[3], const float hi,
				 const float hj, float A[3], float A1[3],
				 float A2[3]) {

  /* Get some useful quantities */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;
  const float Vi = pi->geometry.volume;
  const float Vj = pj->geometry.volume;

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute kernel of pj. */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute (square of) area */
  /* eqn. (7) */
  float Xi = Vi;
  float Xj = Vj;
  if (fvpm_part_geometry_well_behaved(pi) &&
      fvpm_part_geometry_well_behaved(pj)) {
    /* in principle, we use Vi and Vj as weights for the left and right
     * contributions to the generalized surface vector.
     * However, if Vi and Vj are very different (because they have very
     * different smoothing lengths), then the expressions below are more
     * stable. */
#ifdef GIZMO_VOLUME_CORRECTION
    if (fabsf(Vi - Vj) / min(Vi, Vj) > 1.5f * hydro_dimension) {
      Xi = (Vi * hj + Vj * hi) / (hi + hj);
      Xj = Xi;
    }
#endif
    const float dxi[3] = {dx[0] + ci[0], dx[1] + ci[1], dx[2] + ci[2]};
    const float dxj[3] = {dx[0] - cj[0], dx[1] - cj[1], dx[2] - cj[2]};
    for (int k = 0; k < 3; k++) {
      /* minus sign: dx = pi->x - pj->x = -d_ij */
      A1[k] = -Xi * (Bi[k][0] * dxi[0] + Bi[k][1] * dxi[1] + Bi[k][2] * dxi[2]) *
		 wi * hi_inv_dim;
      A2[k] = -Xj * (Bj[k][0] * dxj[0] + Bj[k][1] * dxj[1] + Bj[k][2] * dxj[2]) *
		 wj * hj_inv_dim;
      A[k] = A1[k] + A2[k];
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    const float hidp1 = pow_dimension_plus_one(hi_inv);
    const float hjdp1 = pow_dimension_plus_one(hj_inv);
    const float Anorm_i = -hidp1 * Vi * Vi * wi_dx * r_inv;
    const float Anorm_j = -hjdp1 * Vj * Vj * wj_dx * r_inv;
    for (int k = 0; k < 3; k++) {
      A1[k] = -Anorm_i * dx[k];
      A2[k] = -Anorm_j * dx[k];
      A[k] = A1[k] + A2[k];
    }
  }
}

/**
 * @brief Accumulate the face area vector and norm.
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param interaction_mode 0 if non-symmetric interaction, 1 if symmetric.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_total_face_area_vector_and_norm(struct part *restrict pi,
						struct part *restrict pj, const float r2,
						const float dx[3],
						const float hi, const float hj,
						const int interaction_mode) {
  
  /* Initialize local variables. This is the diagnostic face: it reads
     matrix_E_centred/centroid_offset, so it is centred as of candidate A.
     Production faces (runner_iact_fluxes_common, rt_iact.h) are untouched
     by Stage 1a and stay on matrix_E with zero offsets. */
  float Bi[3][3];
  float Bj[3][3];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E_centred[k][l];
      Bj[k][l] = pj->geometry.matrix_E_centred[k][l];
    }
  }
  const float ci[3] = {pi->geometry.centroid_offset[0],
                       pi->geometry.centroid_offset[1],
                       pi->geometry.centroid_offset[2]};
  const float cj[3] = {pj->geometry.centroid_offset[0],
                       pj->geometry.centroid_offset[1],
                       pj->geometry.centroid_offset[2]};

  /* Compute (square of) area */
  float A[3] = {0.0, 0.0, 0.0};
  float A1[3] = {0.0, 0.0, 0.0};
  float A2[3] = {0.0, 0.0, 0.0};
  fvpm_compute_face_area_vector(pi, pj, Bi, Bj, ci, cj, r2, dx, hi, hj, A, A1, A2);
  const float Anorm2 = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];
  const float Anorm = sqrtf(Anorm2);

#ifdef SWIFT_DEBUG_CHECKS
  {
    /* Swap EVERY role: pi<->pj, Bi<->Bj, ci<->cj, hi<->hj, and NEGATE dx.
       r2 is unchanged. Outputs go to separate locals so the accumulation
       below is untouched. Bit-exact equality is achievable (see the plan's
       item 9a derivation), so a tolerance would hide a real sign error. */
    const float dx_swapped[3] = {-dx[0], -dx[1], -dx[2]};
    float As[3], A1s[3], A2s[3];
    fvpm_compute_face_area_vector(pj, pi, Bj, Bi, cj, ci, r2, dx_swapped,
                                  hj, hi, As, A1s, A2s);
    for (int k = 0; k < 3; k++) {
      if (As[k] != -A[k])
        error("Face antisymmetry violated: A_ji[%d] = %.9e, -A_ij[%d] = %.9e "
              "(pi=%lld, pj=%lld)", k, As[k], k, -A[k], pi->id, pj->id);
    }
  }
#endif

  if (pi->time_bin == 0 && pj->time_bin == 0) {
    /* All particles are logged */    
    /* message("[%lld %lld]", pi->id, pj->id); */

    /* Update the total face area */
    pi->geometry.area += Anorm;
    if (interaction_mode == 1) {
      pj->geometry.area += Anorm;
    }

    /* Update the face area vectorial sum */
    pi->geometry.area_sum[0] += A[0];
    pi->geometry.area_sum[1] += A[1];
    pi->geometry.area_sum[2] += A[2];
    if (interaction_mode == 1) {
      /* We add a minus sign since the faces are antisymmetric */
      pj->geometry.area_sum[0] -= A[0];
      pj->geometry.area_sum[1] -= A[1];
      pj->geometry.area_sum[2] -= A[2];
    }

    /* Split of area_sum: A1 uses only particle i's own data, A2 depends on neighbour j. */
    pi->geometry.area_sum1[0] += A1[0];
    pi->geometry.area_sum1[1] += A1[1];
    pi->geometry.area_sum1[2] += A1[2];
    pi->geometry.area_sum2[0] += A2[0];
    pi->geometry.area_sum2[1] += A2[1];
    pi->geometry.area_sum2[2] += A2[2];
    if (interaction_mode == 1) {
      /* For A_ji the two terms swap roles and flip sign (verified against A_ji = -A_ij). */
      pj->geometry.area_sum1[0] -= A2[0];
      pj->geometry.area_sum1[1] -= A2[1];
      pj->geometry.area_sum1[2] -= A2[2];
      pj->geometry.area_sum2[0] -= A1[0];
      pj->geometry.area_sum2[1] -= A1[1];
      pj->geometry.area_sum2[2] -= A1[2];
    }
  } else {
    if (pi->geometry.is_problematic == 2) {
      const float Si[3] = {pi->geometry.area_sum[0], pi->geometry.area_sum[1],
                           pi->geometry.area_sum[2]};
      float Si_times_Aij = Si[0] * A[0] + Si[1] * A[1] + Si[2] * A[2];

      if (Si_times_Aij > 0.0) {
	pi->geometry.area_sum_plus[0] += Si_times_Aij*A[0];
	pi->geometry.area_sum_plus[1] += Si_times_Aij*A[1];
	pi->geometry.area_sum_plus[2] += Si_times_Aij*A[2];
      } else {
	pi->geometry.area_sum_minus[0] += Si_times_Aij*A[0];
	pi->geometry.area_sum_minus[1] += Si_times_Aij*A[1];
        pi->geometry.area_sum_minus[2] += Si_times_Aij*A[2];
      }

      message(
	      "[%lld %lld, i] Debug: Si_times_ai = %e, Sum_+ = (%e %e %e), Sum_- = (%e %e %e)",
	      pi->id, pj->id, Si_times_Aij, pi->geometry.area_sum_plus[0], pi->geometry.area_sum_plus[1],
	      pi->geometry.area_sum_plus[2], pi->geometry.area_sum_minus[0],
	      pi->geometry.area_sum_minus[1], pi->geometry.area_sum_minus[2]);
    }

    if (pj->geometry.is_problematic == 2) {
      const float Sj[3] = {pj->geometry.area_sum[0], pj->geometry.area_sum[1],
                           pj->geometry.area_sum[2]};

      /* The minus is because A = Aij and we need Aji = - Aij */      
      float Sj_times_Aij = - Sj[0] * A[0] - Sj[1] * A[1] - Sj[2] * A[2];

      if (Sj_times_Aij > 0.0) {
	/* The minus is because A = Aij and we need Aji = - Aij */        
	pj->geometry.area_sum_plus[0] -= Sj_times_Aij*A[0];
	pj->geometry.area_sum_plus[1] -= Sj_times_Aij*A[1];
	pj->geometry.area_sum_plus[2] -= Sj_times_Aij*A[2];
      } else {
	/* The minus is because A = Aij and we need Aji = - Aij */        
	pj->geometry.area_sum_minus[0] -= Sj_times_Aij*A[0];
	pj->geometry.area_sum_minus[1] -= Sj_times_Aij*A[1];
	pj->geometry.area_sum_minus[2] -= Sj_times_Aij*A[2];
      }
      
      message(
	      "[%lld %lld, j] Debug: Sj_times_ai = %e, Sum_+ = (%e %e %e), Sum_- = (%e %e %e)",
	      pi->id, pj->id, Sj_times_Aij, pj->geometry.area_sum_plus[0], pj->geometry.area_sum_plus[1],
	      pj->geometry.area_sum_plus[2], pj->geometry.area_sum_minus[0],
	      pj->geometry.area_sum_minus[1], pj->geometry.area_sum_minus[2]);
    }
    
  }    
}

/**
 * @brief Check that the total face area is close to 0.0.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
fvpm_check_total_face_area_vector_sum(struct part *p) {

  const float threshold = 1e-2;
  const float area_threshold = p->geometry.area*threshold;

  if ((fabsf(p->geometry.area_sum[0]) > area_threshold ||
       fabsf(p->geometry.area_sum[1]) > area_threshold ||
       fabsf(p->geometry.area_sum[2]) > area_threshold)) {
    warning(
        "[%lld] Sum A_ij strongly deviating from 0! A_tot = %e, Sum_j A_ij = ( "
        "%e %e %e ), Sum_j T1_ij = ( %e %e %e ), Sum_j T2_ij = ( %e %e %e ).",
        p->id, p->geometry.area, p->geometry.area_sum[0],
        p->geometry.area_sum[1], p->geometry.area_sum[2],
        p->geometry.area_sum1[0], p->geometry.area_sum1[1],
        p->geometry.area_sum1[2], p->geometry.area_sum2[0],
        p->geometry.area_sum2[1], p->geometry.area_sum2[2]);

    /* 0 = no problem, 2 = has just been detected skip this timestep's force, 1
       = it's safe to correct */
    if (p->geometry.is_problematic == 0) {
      p->geometry.is_problematic = 2;
    } else if (p->geometry.is_problematic == 2) {
      p->geometry.is_problematic = 1;
    } else {
      p->geometry.is_problematic = 1;      
    }

    /* If the particle is not problematic, flag it for the next-timestep. If is
       problematic, unflag it to avoid redoing the computations in the next time
       step. */
    if (p->geometry.is_problematic == 1) {

      const float beta_i[3] = {
        -p->geometry.area_sum_minus[0] / p->geometry.area_sum_plus[0],
        -p->geometry.area_sum_minus[1] / p->geometry.area_sum_plus[1],
        -p->geometry.area_sum_minus[2] / p->geometry.area_sum_plus[2],        
      };

      message(
          "[%lld] Debug: Sum_+ = (%e %e %e), Sum_- = (%e %e %e), alpha = (%e "
          "%e %e)",
          p->id, p->geometry.area_sum_plus[0], p->geometry.area_sum_plus[1],
          p->geometry.area_sum_plus[2], p->geometry.area_sum_minus[0],
          p->geometry.area_sum_minus[1], p->geometry.area_sum_minus[2],
          beta_i[0], beta_i[1], beta_i[2]);
      
    }
  }
}

#endif /* SWIFT_FVPM_GEOMETRY_GIZMO_H */
