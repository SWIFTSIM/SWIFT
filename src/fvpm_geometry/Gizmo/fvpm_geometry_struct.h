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

  char is_problematic;

#if defined(RT_GEAR) && defined(FVPM_RT_FACE_CLOSURE_DIAGNOSTIC)
  /* Separate from area/area_sum1/2 above: RT_GEAR can run on gizmo-mfv
     hydro too, so those fields (hydro gradient loop) and these (RT
     transport loop) must not share storage. Behind a flag (off by
     default, like FVPM_ATTEMPT1_FACE_RESCALING) because the check this
     feeds is an unrated warning() that would flood production runs at
     every h-discontinuity otherwise. */

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
