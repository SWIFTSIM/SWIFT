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
};

#endif /* SWIFT_FVPM_GEOMETRY_STRUCT_GIZMO_H */
