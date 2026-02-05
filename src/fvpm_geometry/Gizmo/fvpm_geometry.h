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

  /* reset the centroid variables used for the velocity correction in MFV */
  fvpm_reset_centroids(p);

  p->geometry.area = 0.0;
  p->geometry.area_sum[0] = 0.0;
  p->geometry.area_sum[1] = 0.0;
  p->geometry.area_sum[2] = 0.0;
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

  /* reset the centroid variables used for the velocity correction in MFV */
  fvpm_reset_centroids(p);

  /* TODO: To be defined */
  p->geometry.area = 0.0;
  p->geometry.area_sum[0] = 0.0;
  p->geometry.area_sum[1] = 0.0;
  p->geometry.area_sum[2] = 0.0;
}

/**
 * @brief Finish the computation of the matrix.
 *
 * @param p the particle to work on
 * @param ihdim 1/h^{dim}
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_volume_and_matrix(struct part *restrict p, const float ihdim) {

  /* Final operation on the geometry. */
  /* we multiply with the smoothing kernel normalization ih3 and calculate the
   * volume */
  const float volume_inv = ihdim * (p->geometry.volume + kernel_root);
  const float volume = 1.0f / volume_inv;
  p->geometry.volume = volume;

  /* we multiply with the smoothing kernel normalization */
  p->geometry.matrix_E[0][0] *= ihdim;
  p->geometry.matrix_E[0][1] *= ihdim;
  p->geometry.matrix_E[0][2] *= ihdim;
  p->geometry.matrix_E[1][0] *= ihdim;
  p->geometry.matrix_E[1][1] *= ihdim;
  p->geometry.matrix_E[1][2] *= ihdim;
  p->geometry.matrix_E[2][0] *= ihdim;
  p->geometry.matrix_E[2][1] *= ihdim;
  p->geometry.matrix_E[2][2] *= ihdim;

  /* normalise the centroids for MFV */
  fvpm_normalise_centroid(p, p->density.wcount);

  /* Check the condition number to see if we have a stable geometry. */
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
  if (invert_dimension_by_dimension_matrix(p->geometry.matrix_E) != 0) {
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
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param (return) A The face area between i and j.
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_face_area_vector(const struct part *restrict pi, const struct part *restrict pj,
				 float Bi[3][3], float Bj[3][3], const float r2,
				 const float dx[3], const float hi,
				 const float hj, float A[3]) {

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
    for (int k = 0; k < 3; k++) {
      /* we add a minus sign since dx is pi->x - pj->x */
      A[k] = -Xi * (Bi[k][0] * dx[0] + Bi[k][1] * dx[1] + Bi[k][2] * dx[2]) *
		 wi * hi_inv_dim -
	     Xj * (Bj[k][0] * dx[0] + Bj[k][1] * dx[1] + Bj[k][2] * dx[2]) *
		 wj * hj_inv_dim;
    }
  } else {
    /* ill condition gradient matrix: revert to SPH face area */
    const float hidp1 = pow_dimension_plus_one(hi_inv);
    const float hjdp1 = pow_dimension_plus_one(hj_inv);
    const float Anorm =
	-(hidp1 * Vi * Vi * wi_dx + hjdp1 * Vj * Vj * wj_dx) * r_inv;
    A[0] = -Anorm * dx[0];
    A[1] = -Anorm * dx[1];
    A[2] = -Anorm * dx[2];
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

  /* Initialize local variables */
  float Bi[3][3];
  float Bj[3][3];
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }

  /* Compute (square of) area */
  float A[3] = {0.0, 0.0, 0.0};
  fvpm_compute_face_area_vector(pi, pj, Bi, Bj, r2, dx, hi, hj, A);
  const float Anorm2 = A[0] * A[0] + A[1] * A[1] + A[2] * A[2];
  const float Anorm = sqrtf(Anorm2);

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
    pj->geometry.area_sum[0] -= A[1];
    pj->geometry.area_sum[1] -= A[1];
    pj->geometry.area_sum[2] -= A[2];
  }
}

/**
 * @brief Check that the total face area is close to 0.0.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void
fvpm_check_total_face_area_vector_sum(const struct part *p) {

  const float threshold = 1e-2;

  if (p->geometry.area_sum[0] > threshold ||
      p->geometry.area_sum[1] > threshold ||
      p->geometry.area_sum[2] > threshold) {
    warning("[%lld] Sum A_ij strongly deviating from 0! A_tot = %e, Sum_j A_ij = ( %e %e %e ).", p->id,
	    p->geometry.area, p->geometry.area_sum[0], p->geometry.area_sum[1],
	    p->geometry.area_sum[2]);
  }
}

#endif /* SWIFT_FVPM_GEOMETRY_GIZMO_H */
