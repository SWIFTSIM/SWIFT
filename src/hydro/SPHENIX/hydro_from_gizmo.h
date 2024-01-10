#ifndef SWIFT_SPHENIX_HYDRO_FROM_GIZMO_H
#define SWIFT_SPHENIX_HYDRO_FROM_GIZMO_H

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
hydro_part_geometry_well_behaved(const struct part* restrict p) {

  return p->geometry.wcorr > const_gizmo_min_wcorr;
}



#endif
