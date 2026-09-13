#ifndef SWIFT_FVPM_GEOMETRY_STRUCT_GIZMO_H
#define SWIFT_FVPM_GEOMETRY_STRUCT_GIZMO_H

/**
 * @file Gizmo/fvpm_geometry_struct.h
 * @brief Struct related to the Gizmo FVPM geometry particle data collection,
 * in particular the collection of the data required for the matrix needed
 * for gradients.
 * This was moved here so we can cleanly couple GEAR-RT on top of SPH
 * hydrodynamics while avoiding code replication.
 */

/* Geometrical quantities used for hydro. */
struct fvpm_geometry_struct {

  /* Volume of the particle. */
  float volume;

  /* Geometrical shear matrix used to calculate second order accurate
     gradients */
  float matrix_E[3][3];

  /* Centroid of the "cell". */
  float centroid[3];

  /* Correction factor for wcount. */
  float wcorr;

  /*! Condition number of matrix_E (eq C1) */
  float condition_number;

  /*! Total particle area */
  float area;

  /*! Sum of the face area vectors */
  float area_sum[3];

  float area_sum_plus[3];
  float area_sum_minus[3];

  /*! Row-sum term of area_sum: -V_i (B_i . dx) w_i h_i^-dim, uses only particle i's own data */
  float area_sum1[3];

  /*! Column-sum term of area_sum: -V_j (B_j . dx) w_j h_j^-dim, depends on neighbour j's data */
  float area_sum2[3];

  /*! Raw first moment m_i = Sum_{j!=i} d_ij W_ij, accumulated in the density
      loop. Not normalised, no h^-nu factor. Also the source of MFV's
      geometry.centroid (see fvpm_normalise_centroid). */
  float first_moment[3];

  /*! Kernel centroid offset c_i = m_i / omega'_i, omega' excluding the self term. */
  float centroid_offset[3];

  /*! B^c_i = (E_i - omega'_i c_i (x) c_i)^{-1}, kernel-normalised. SEPARATE from
      matrix_E, which is never modified by candidate A. */
  float matrix_E_centred[3][3];

  /*! SPD margin q = ihdim.omega'.c^T B c. Centring is applied only if q < 1-tau. */
  float centring_margin;

  char is_problematic;

  /*! 1 if centring was not applied for this particle this step. */
  char centring_disabled;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Item 9b: role-aware Sum_j psitilde^c_j(x_i), accumulated in the
      gradient loop and reset in hydro_gradients_init (once per gradient
      loop, not once per h-iteration). Must equal 0 up to round-off. */
  float psi_c_sum[3];

  /*! Item 9b: Sum_j |psitilde^c_j(x_i)|, the normalisation for the check. */
  float psi_c_abs_sum[3];
#endif

#if defined(RT_GEAR) && defined(FVPM_RT_FACE_CLOSURE_DIAGNOSTIC)
  /* Separate from area/area_sum1/2 above: RT_GEAR can run on gizmo-mfv
     hydro too, so those fields (hydro gradient loop) and these (RT
     transport loop) must not share storage. Behind a flag (off by
     default) because the check this feeds is an unrated warning() that
     would flood production runs at every h-discontinuity otherwise. */

  /*! Total RT face area, summed over the transport/flux loop. */
  float rt_area;

  /*! Sum of the RT face area vectors over the transport/flux loop. */
  float rt_area_sum[3];

  /*! Row-sum term of rt_area_sum (uses only particle i's own data). */
  float rt_area_sum1[3];

  /*! Column-sum term of rt_area_sum (depends on neighbour j's data). */
  float rt_area_sum2[3];
#endif
};

#endif /* SWIFT_FVPM_GEOMETRY_STRUCT_GIZMO_H */
