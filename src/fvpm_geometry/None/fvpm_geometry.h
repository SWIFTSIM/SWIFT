#ifndef SWIFT_FVPM_GEOMETRY_NONE_H
#define SWIFT_FVPM_GEOMETRY_NONE_H

/**
 * @file None/fvpm_geometry.h
 * @brief Functions related to the Gizmo FVPM geometry struct collection,
 * in particular the collection of the data required for the matrix needed
 * for gradients. Empty definitions for when FVPM are unused.
 */

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void fvpm_reset_centroids(
    struct part *restrict p) {}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param p Particle.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void fvpm_normalise_centroid(
    struct part *restrict p, const float wcount) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param p Particle (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_left(
    struct part *restrict p, const float *dx, const float w) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param p Particle (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_right(
    struct part *restrict p, const float *dx, const float w) {}

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
fvpm_part_geometry_well_behaved(const struct part *restrict p) {
  return 0;
}

/**
 * @brief Collect the data needed for the matrix construction.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_geometry_and_matrix(struct part *restrict pi, const float wi,
                                    const float dx[3]) {}

__attribute__((always_inline)) INLINE static void fvpm_geometry_init(
    struct part *restrict p) {}

/**
 * @brief Sets the geometry fields to sensible values when #part has 0 ngbs.
 *
 * @param p the particle to work on
 */
__attribute__((always_inline)) INLINE static void
fvpm_geometry_part_has_no_neighbours(struct part *restrict p) {}

/**
 * @brief Finish the computation of the matrix.
 *
 * @param p the particle to work on
 * @param ihdim 1/h^{dim}
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_volume_and_matrix(struct part *restrict p, const float ihdim) {}

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
__attribute__((always_inline)) INLINE static void fvpm_compute_face_area_vector(
    const struct part *pi, const struct part *pj, float Bi[3][3],
    float Bj[3][3], const float r2, const float dx[3], const float hi,
    const float hj, float A[3]) {}

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
 * @param (return) A The face area between i and j.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_total_face_area_vector_and_norm(struct part *pi,
                                                struct part *pj, const float r2,
                                                const float dx[3],
                                                const float hi, const float hj,
                                                const int interaction_mode) {}

#endif /* SWIFT_FVPM_GEOMETRY_NONE_H */
