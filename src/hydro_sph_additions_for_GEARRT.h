#ifndef HYDRO_SPH_ADDITIONS_FOR_GEARRT_H
#define HYDRO_SPH_ADDITIONS_FOR_GEARRT_H



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


/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_reset_centroids(struct part* restrict p) {

  p->geometry.centroid[0] = 0.0f;
  p->geometry.centroid[1] = 0.0f;
  p->geometry.centroid[2] = 0.0f;
}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param p Particle.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_normalise_centroid(struct part* restrict p,
                                    const float wcount) {

  const float norm = kernel_norm / wcount;
  p->geometry.centroid[0] *= norm;
  p->geometry.centroid[1] *= norm;
  p->geometry.centroid[2] *= norm;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param p Particle (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_update_centroid_left(struct part* restrict p, const float* dx,
                                      const float w) {

  p->geometry.centroid[0] -= dx[0] * w;
  p->geometry.centroid[1] -= dx[1] * w;
  p->geometry.centroid[2] -= dx[2] * w;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param p Particle (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void
hydro_velocities_update_centroid_right(struct part* restrict p, const float* dx,
                                       const float w) {

  p->geometry.centroid[0] += dx[0] * w;
  p->geometry.centroid[1] += dx[1] * w;
  p->geometry.centroid[2] += dx[2] * w;
}



// TODO: Maybe move this into rt_additions.h
__attribute__((always_inline)) INLINE static void gearrt_density_accumulate_geometry_and_matrix(struct part* restrict pi, const float wi, const float dx[3]) {
#ifdef RT_GEAR
  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  hydro_velocities_update_centroid_left(pi, dx, wi);
#endif
  }

__attribute__((always_inline)) INLINE static void gearrt_geometry_init(struct part* restrict p){
#ifdef RT_GEAR

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
  hydro_velocities_reset_centroids(p);
#endif
  }


__attribute__((always_inline)) INLINE static void gearrt_compute_volume_and_matrix(struct part* restrict p, const float ihdim){

#ifdef RT_GEAR

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
  hydro_velocities_normalise_centroid(p, p->density.wcount);

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

  float condition_number = 0.0f;
  if (invert_dimension_by_dimension_matrix(p->geometry.matrix_E) != 0) {
    /* something went wrong in the inversion; force bad condition number */
    condition_number = const_gizmo_max_condition_number + 1.0f;
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

    condition_number =
        hydro_dimension_inv * sqrtf(condition_number_E * condition_number_Einv);
  }

  if (condition_number > const_gizmo_max_condition_number &&
      p->geometry.wcorr > const_gizmo_min_wcorr) {
#ifdef GIZMO_PATHOLOGICAL_ERROR
    error("Condition number larger than %g (%g)!",
            const_gizmo_max_condition_number, condition_number);
#endif
#ifdef GIZMO_PATHOLOGICAL_WARNING
    message("Condition number too large: %g (> %g, p->id: %llu)!",
            condition_number, const_gizmo_max_condition_number, p->id);
#endif
    /* add a correction to the number of neighbours for this particle */
    p->geometry.wcorr = const_gizmo_w_correction_factor * p->geometry.wcorr;
  }
#endif

}



#endif
