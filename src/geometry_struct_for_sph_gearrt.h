#ifndef GEOMETRY_STRUCT_FOR_GEARRT_H
#define GEOMETRY_STRUCT_FOR_GEARRT_H

/* Geometrical quantities used for hydro. */
struct geometry_struct_for_rt {

#ifdef RT_GEAR
  /* Volume of the particle. */
  float volume;

  /* Geometrical shear matrix used to calculate second order accurate
     gradients */
  float matrix_E[3][3];

  /* Centroid of the "cell". */
  float centroid[3];

  /* Correction factor for wcount. */
  float wcorr;
#endif
};

#endif
