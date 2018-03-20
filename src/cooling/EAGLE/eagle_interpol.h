#ifndef SWIFT_EAGLE_INTERPOL_H
#define SWIFT_EAGLE_INTERPOL_H

#include "cooling_read_table.h"

/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE void get_index_1d(float *table, int ntable, double x, int *i, float *dx) {
  float dxm1;
  const float EPS = 1.e-4;

  dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

  if ((float)x <= table[0] + EPS) {
    *i = 0;
    *dx = 0;
  } else if ((float)x >= table[ntable - 1] - EPS) {
    *i = ntable - 2;
    *dx = 1;
  } else {
    *i = (int)floor(((float)x - table[0]) * dxm1);
    *dx = ((float)x - table[*i]) * dxm1;
  }
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_1d(float *table, int i, float dx) {
  float result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_1d_dbl(double *table, int i, float dx) {
  double result;

  result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_2d(float *table, int i, int j, float dx, float dy, int size) {
  float result;
  int index[4];
  
  index[0] = row_major_index_2d(i,j,size);
  index[1] = row_major_index_2d(i,j+1,size);
  index[2] = row_major_index_2d(i+1,j,size);
  index[3] = row_major_index_2d(i+1,j+1,size);
  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE double interpol_2d_dbl(double *table, int i, int j, double dx, double dy, int size) {
  double result;
  int index[4];
  
  index[0] = row_major_index_2d(i,j,size);
  index[1] = row_major_index_2d(i,j+1,size);
  index[2] = row_major_index_2d(i+1,j,size);
  index[3] = row_major_index_2d(i+1,j+1,size);
  result = (1 - dx) * (1 - dy) * table[index[0]] + (1 - dx) * dy * table[index[1]] +
           dx * (1 - dy) * table[index[2]] + dx * dy * table[index[3]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_3d(float *table, int i, int j, int k, float dx, float dy,
                  float dz, int size) {
  float result;
  int index[8];

  index[0] = row_major_index_3d(i,j,k,size);
  index[1] = row_major_index_3d(i,j,k+1,size);
  index[2] = row_major_index_3d(i,j+1,k,size);
  index[3] = row_major_index_3d(i,j+1,k+1,size);
  index[4] = row_major_index_3d(i+1,j,k,size);
  index[5] = row_major_index_3d(i+1,j,k+1,size);
  index[6] = row_major_index_3d(i+1,j+1,k,size);
  index[7] = row_major_index_3d(i+1,j+1,k+1,size);
  result = (1 - dx) * (1 - dy) * (1 - dz) * table[index[0]] +
           (1 - dx) * (1 - dy) * dz * table[index[1]] +
           (1 - dx) * dy * (1 - dz) * table[index[2]] +
           (1 - dx) * dy * dz * table[index[3]] +
           dx * (1 - dy) * (1 - dz) * table[index[4]] +
           dx * (1 - dy) * dz * table[index[5]] +
           dx * dy * (1 - dz) * table[index[6]] +
           dx * dy * dz * table[index[7]];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

__attribute__((always_inline)) INLINE float interpol_4d(float *table, int i, int j, int k, int l, float dx, 
			float dy, float dz, float dw, int size) {
  float result;
  int index[16];

  index[0]  = row_major_index_4d(i,j,k,l,size);
  index[1]  = row_major_index_4d(i,j,k,l+1,size);
  index[2]  = row_major_index_4d(i,j,k+1,l,size);
  index[3]  = row_major_index_4d(i,j,k+1,l+1,size);
  index[4]  = row_major_index_4d(i,j+1,k,l,size);
  index[5]  = row_major_index_4d(i,j+1,k,l+1,size);
  index[6]  = row_major_index_4d(i,j+1,k+1,l,size);
  index[7]  = row_major_index_4d(i,j+1,k+1,l+1,size);
  index[8]  = row_major_index_4d(i+1,j,k,l,size);
  index[9]  = row_major_index_4d(i+1,j,k,l+1,size);
  index[10] = row_major_index_4d(i+1,j,k+1,l,size);
  index[11] = row_major_index_4d(i+1,j,k+1,l+1,size);
  index[12] = row_major_index_4d(i+1,j+1,k,l,size);
  index[13] = row_major_index_4d(i+1,j+1,k,l+1,size);
  index[14] = row_major_index_4d(i+1,j+1,k+1,l,size);
  index[15] = row_major_index_4d(i+1,j+1,k+1,l+1,size);

  result = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table[index[0]] +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table[index[1]] +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table[index[2]] +
           (1 - dx) * (1 - dy) * dz * dw * table[index[3]] +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table[index[4]] +
           (1 - dx) * dy * (1 - dz) * dw * table[index[5]] +
           (1 - dx) * dy * dz * (1 - dw) * table[index[6]] +
           (1 - dx) * dy * dz * dw * table[index[7]] +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table[index[8]] +
           dx * (1 - dy) * (1 - dz) * dw * table[index[9]] +
           dx * (1 - dy) * dz * (1 - dw) * table[index[10]] +
           dx * (1 - dy) * dz * dw * table[index[11]] +
           dx * dy * (1 - dz) * (1 - dw) * table[index[12]] +
           dx * dy * (1 - dz) * dw * table[index[13]] +
           dx * dy * dz * (1 - dw) * table[index[14]] +
           dx * dy * dz * dw * table[index[15]];

  return result;
}

#endif
