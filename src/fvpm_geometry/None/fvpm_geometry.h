#ifndef SWIFT_FVPM_GEOMETRY_NONE_H
#define SWIFT_FVPM_GEOMETRY_NONE_H

/**
 * @file fvpm_geometry.h
 * @brief Functions related to the Gizmo FVPM geometry struct collection,
 * in particular the collection of the data required for the matrix needed
 * for gradients.
 * This was moved here so we can cleanly couple GEAR-RT on top of SPH
 * hydrodynamics while avoiding code replication.
 */

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
fvpm_part_geometry_well_behaved(const struct part* restrict p) {
  return 0;
}


/**
 * @brief Collect the data needed for the matrix construction.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_geometry_and_matrix(struct part* restrict pi,
                                              const float wi,
                                              const float dx[3]) {}

__attribute__((always_inline)) INLINE static void fvpm_geometry_init(
    struct part* restrict p) {}


/**
 * @brief Finish the computation of the matrix.
 *
 * @param p the particle to work on
 * @param ihdim 1/h^{dim}
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_volume_and_matrix(struct part* restrict p, const float ihdim) {}


#endif /* SWIFT_FVPM_GEOMETRY_NONE_H */
