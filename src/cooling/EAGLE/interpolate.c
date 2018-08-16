/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "interpolate.h"
#include "chemistry.h"
#include "cooling_struct.h"
#include "eagle_cool_tables.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

static int get_redshift_index_first_call = 0;
static int get_redshift_index_previous = -1;
  
const float EPS = 1.e-4;

/**
 * @brief Returns the 1d index of element with 2d indices i,j
 * from a flattened 2d array in row major order
 *
 * @param i, j Indices of element of interest
 * @param nx, ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_2d(int i, int j,
                                                             int nx, int ny) {
  return i * ny + j;
}

/**
 * @brief Returns the 1d index of element with 3d indices i,j,k
 * from a flattened 3d array in row major order
 *
 * @param i, j, k Indices of element of interest
 * @param nx, ny, nz Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_3d(int i, int j,
                                                             int k, int nx,
                                                             int ny, int nz) {
  return i * ny * nz + j * nz + k;
}

/**
 * @brief Returns the 1d index of element with 4d indices i,j,k,l
 * from a flattened 4d array in row major order
 *
 * @param i, j, k, l Indices of element of interest
 * @param nx, ny, nz, nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_4d(int i, int j,
                                                             int k, int l,
                                                             int nx, int ny,
                                                             int nz, int nw) {
  return i * ny * nz * nw + j * nz * nw + k * nw + l;
}

/*
 * @brief This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 *
 * @param table Pointer to table of values
 * @param ntable Size of the table
 * @param x Value whose position within the array we are interested in
 * @param i Pointer to the index whose corresponding table value
 * is the greatest value in the table less than x
 * @param dx Pointer to offset of x within table cell
 */
__attribute__((always_inline)) INLINE void get_index_1d(float *table,
                                                        int ntable, double x,
                                                        int *i, float *dx) {

  float dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

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
 * @brief Returns the position i of a value z in the redshift table
 * and computes the displacement dz needed for the interpolation.
 * Since the redshift table is not evenly spaced, compare z with each
 * table value in decreasing order starting with the previous redshift index
 *
 * @param z Redshift whose position within the redshift array we are interested
 * in
 * @param z_index Pointer to the index whose corresponding redshift
 * is the greatest value in the redshift table less than x
 * @param dz Pointer to offset of z within redshift cell
 * @param cooling Pointer to cooling structure containing redshift table
 */
__attribute__((always_inline)) INLINE void get_redshift_index(
    float z, int *z_index, float *dz,
    const struct cooling_function_data *restrict cooling) {
  int i, iz;

  if (get_redshift_index_first_call == 0) {
    get_redshift_index_first_call = 1;
    get_redshift_index_previous = cooling->N_Redshifts - 2;

    /* this routine assumes cooling_redshifts table is in increasing order. Test
     * this. */
    for (i = 0; i < cooling->N_Redshifts - 2; i++)
      if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
        error("[get_redshift_index]: table should be in increasing order\n");
      }
  }

  /* before the earliest redshift or before hydrogen reionization, flag for
   * collisional cooling */
  if (z > cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts;
    *dz = 0.0;
  }
  /* from reionization use the cooling tables */
  else if (z > cooling->Redshifts[cooling->N_Redshifts - 1] &&
           z <= cooling->reionisation_redshift) {
    *z_index = cooling->N_Redshifts + 1;
    *dz = 0.0;
  }
  /* at the end, just use the last value */
  else if (z <= cooling->Redshifts[0]) {
    *z_index = 0;
    *dz = 0.0;
  } else {
    /* start at the previous index and search */
    for (iz = get_redshift_index_previous; iz >= 0; iz--) {
      if (z > cooling->Redshifts[iz]) {
        *dz = (z - cooling->Redshifts[iz]) /
              (cooling->Redshifts[iz + 1] - cooling->Redshifts[iz]);

        get_redshift_index_previous = *z_index = iz;

        break;
      }
    }
  }
}

/*
 * @brief Performs 1d interpolation
 *
 * @param table Pointer to 1d table of values
 * @param i Index of cell we are interpolating
 * @param dx Offset within cell
 */
__attribute__((always_inline)) INLINE float interpol_1d(float *table, int i,
                                                        float dx) {

  return (1 - dx) * table[i] + dx * table[i + 1];
}

/*
 * @brief Performs 1d interpolation (double precision)
 *
 * @param table Pointer to 1d table of values
 * @param i Index of cell we are interpolating
 * @param dx Offset within cell
 */
__attribute__((always_inline)) INLINE double interpol_1d_dbl(double *table,
                                                             int i, float dx) {

  return (1 - dx) * table[i] + dx * table[i + 1];
}

/*
 * @brief Performs 2d interpolation
 *
 * @param table Pointer to flattened 2d table of values
 * @param i,j Indices of cell we are interpolating
 * @param dx,dy Offset within cell
 * @param nx,ny Table dimensions
 * @param upper Pointer to value set to the table value at
 * the when dy = 1 (used for calculating derivatives)
 * @param lower Pointer to value set to the table value at
 * the when dy = 0 (used for calculating derivatives)
 */
__attribute__((always_inline)) INLINE float interpol_2d(float *table, int i,
                                                        int j, float dx,
                                                        float dy, int nx,
                                                        int ny, double *upper,
							double *lower) {
  float result;
  int index[4];

  index[0] = row_major_index_2d(i, j, nx, ny);
  index[1] = row_major_index_2d(i, j + 1, nx, ny);
  index[2] = row_major_index_2d(i + 1, j, nx, ny);
  index[3] = row_major_index_2d(i + 1, j + 1, nx, ny);

  result = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];
  dy = 1.0;
  *upper = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];
  dy = 0.0;
  *lower = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];

  return result;
}

/*
 * @brief Performs 2d interpolation (double precision)
 *
 * @param table Pointer to flattened 2d table of values
 * @param i,j Indices of cell we are interpolating
 * @param dx,dy Offset within cell
 * @param nx,ny Table dimensions
 * @param upper Pointer to value set to the table value at
 * the when dy = 1 (used for calculating derivatives)
 * @param lower Pointer to value set to the table value at
 * the when dy = 0 (used for calculating derivatives)
 */
__attribute__((always_inline)) INLINE double interpol_2d_dbl(
    double *table, int i, int j, double dx, double dy, int nx, int ny,
    double *upper, double *lower) {
  double result;
  int index[4];

  index[0] = row_major_index_2d(i, j, nx, ny);
  index[1] = row_major_index_2d(i, j + 1, nx, ny);
  index[2] = row_major_index_2d(i + 1, j, nx, ny);
  index[3] = row_major_index_2d(i + 1, j + 1, nx, ny);

  result = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];
  dy = 1.0;
  *upper = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];
  dy = 0.0;
  *lower = (1 - dx) * (1 - dy) * table[index[0]] +
           (1 - dx) * dy * table[index[1]] + 
	   dx * (1 - dy) * table[index[2]] +
           dx * dy * table[index[3]];

  return result;
}

/*
 * @brief Performs 3d interpolation
 *
 * @param table Pointer to flattened 3d table of values
 * @param i,j,k Indices of cell we are interpolating
 * @param dx,dy,dz Offset within cell
 * @param nx,ny,nz Table dimensions
 * @param upper Pointer to value set to the table value at
 * the when dz = 1 (used for calculating derivatives)
 * @param lower Pointer to value set to the table value at
 * the when dz = 0 (used for calculating derivatives)
 */
__attribute__((always_inline)) INLINE float interpol_3d(float *table, int i,
                                                        int j, int k, float dx,
                                                        float dy, float dz,
                                                        int nx, int ny,
                                                        int nz, double *upper, 
							double *lower) {
  float result;
  int index[8];

  index[0] = row_major_index_3d(i, j, k, nx, ny, nz);
  index[1] = row_major_index_3d(i, j, k + 1, nx, ny, nz);
  index[2] = row_major_index_3d(i, j + 1, k, nx, ny, nz);
  index[3] = row_major_index_3d(i, j + 1, k + 1, nx, ny, nz);
  index[4] = row_major_index_3d(i + 1, j, k, nx, ny, nz);
  index[5] = row_major_index_3d(i + 1, j, k + 1, nx, ny, nz);
  index[6] = row_major_index_3d(i + 1, j + 1, k, nx, ny, nz);
  index[7] = row_major_index_3d(i + 1, j + 1, k + 1, nx, ny, nz);
  
  float table0  = table[index[0]];
  float table1  = table[index[1]];
  float table2  = table[index[2]];
  float table3  = table[index[3]];
  float table4  = table[index[4]];
  float table5  = table[index[5]];
  float table6  = table[index[6]];
  float table7  = table[index[7]];

  result = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 +
           (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 +
           dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 +
           dx * dy * dz * table7;
  dz = 1.0;
  *upper = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 +
           (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 +
           dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 +
           dx * dy * dz * table7;
  dz = 0.0;
  *lower = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 +
           (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 +
           dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 +
           dx * dy * dz * table7;

  return result;
}

/*
 * @brief Performs 4d interpolation
 *
 * @param table Pointer to flattened 4d table of values
 * @param i,j,k,l Indices of cell we are interpolating
 * @param dx,dy,dz,dw Offset within cell
 * @param nx,ny,nz,nw Table dimensions
 * @param upper Pointer to value set to the table value at
 * the when dw = 1 when used for interpolating table 
 * depending on metal species, dz = 1 otherwise 
 * (used for calculating derivatives)
 * @param upper Pointer to value set to the table value at
 * the when dw = 0 when used for interpolating table 
 * depending on metal species, dz = 0 otherwise 
 * (used for calculating derivatives)
 */
__attribute__((always_inline)) INLINE float interpol_4d(
    float *table, int i, int j, int k, int l, float dx, float dy, float dz,
    float dw, int nx, int ny, int nz, int nw, double *upper, double *lower) {
  float result;
  int index[16];

  index[0] = row_major_index_4d(i, j, k, l, nx, ny, nz, nw);
  index[1] = row_major_index_4d(i, j, k, l + 1, nx, ny, nz, nw);
  index[2] = row_major_index_4d(i, j, k + 1, l, nx, ny, nz, nw);
  index[3] = row_major_index_4d(i, j, k + 1, l + 1, nx, ny, nz, nw);
  index[4] = row_major_index_4d(i, j + 1, k, l, nx, ny, nz, nw);
  index[5] = row_major_index_4d(i, j + 1, k, l + 1, nx, ny, nz, nw);
  index[6] = row_major_index_4d(i, j + 1, k + 1, l, nx, ny, nz, nw);
  index[7] = row_major_index_4d(i, j + 1, k + 1, l + 1, nx, ny, nz, nw);
  index[8] = row_major_index_4d(i + 1, j, k, l, nx, ny, nz, nw);
  index[9] = row_major_index_4d(i + 1, j, k, l + 1, nx, ny, nz, nw);
  index[10] = row_major_index_4d(i + 1, j, k + 1, l, nx, ny, nz, nw);
  index[11] = row_major_index_4d(i + 1, j, k + 1, l + 1, nx, ny, nz, nw);
  index[12] = row_major_index_4d(i + 1, j + 1, k, l, nx, ny, nz, nw);
  index[13] = row_major_index_4d(i + 1, j + 1, k, l + 1, nx, ny, nz, nw);
  index[14] = row_major_index_4d(i + 1, j + 1, k + 1, l, nx, ny, nz, nw);
  index[15] = row_major_index_4d(i + 1, j + 1, k + 1, l + 1, nx, ny, nz, nw);

  float table0  = table[index[0]];
  float table1  = table[index[1]];
  float table2  = table[index[2]];
  float table3  = table[index[3]];
  float table4  = table[index[4]];
  float table5  = table[index[5]];
  float table6  = table[index[6]];
  float table7  = table[index[7]];
  float table8  = table[index[8]];
  float table9  = table[index[9]];
  float table10 = table[index[10]];
  float table11 = table[index[11]];
  float table12 = table[index[12]];
  float table13 = table[index[13]];
  float table14 = table[index[14]];
  float table15 = table[index[15]];

  result = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table0 +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table1 +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table2 +
           (1 - dx) * (1 - dy) * dz * dw * table3 +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table4 +
           (1 - dx) * dy * (1 - dz) * dw * table5 +
           (1 - dx) * dy * dz * (1 - dw) * table6 +
           (1 - dx) * dy * dz * dw * table7 +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table8 +
           dx * (1 - dy) * (1 - dz) * dw * table9 +
           dx * (1 - dy) * dz * (1 - dw) * table10 +
           dx * (1 - dy) * dz * dw * table11 +
           dx * dy * (1 - dz) * (1 - dw) * table12 +
           dx * dy * (1 - dz) * dw * table13 +
           dx * dy * dz * (1 - dw) * table14 +
           dx * dy * dz * dw * table15;
  if (nw == 9) {
    // interpolating metal species
    dz = 1.0;
  } else {
    dw = 1.0;
  }
  *upper = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table0 +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table1 +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table2 +
           (1 - dx) * (1 - dy) * dz * dw * table3 +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table4 +
           (1 - dx) * dy * (1 - dz) * dw * table5 +
           (1 - dx) * dy * dz * (1 - dw) * table6 +
           (1 - dx) * dy * dz * dw * table7 +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table8 +
           dx * (1 - dy) * (1 - dz) * dw * table9 +
           dx * (1 - dy) * dz * (1 - dw) * table10 +
           dx * (1 - dy) * dz * dw * table11 +
           dx * dy * (1 - dz) * (1 - dw) * table12 +
           dx * dy * (1 - dz) * dw * table13 +
           dx * dy * dz * (1 - dw) * table14 +
           dx * dy * dz * dw * table15;
  if (nw == 9) {
    // interpolating metal species
    dz = 0.0;
  } else {
    dw = 0.0;
  }
  *lower = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table0 +
           (1 - dx) * (1 - dy) * (1 - dz) * dw * table1 +
           (1 - dx) * (1 - dy) * dz * (1 - dw) * table2 +
           (1 - dx) * (1 - dy) * dz * dw * table3 +
           (1 - dx) * dy * (1 - dz) * (1 - dw) * table4 +
           (1 - dx) * dy * (1 - dz) * dw * table5 +
           (1 - dx) * dy * dz * (1 - dw) * table6 +
           (1 - dx) * dy * dz * dw * table7 +
           dx * (1 - dy) * (1 - dz) * (1 - dw) * table8 +
           dx * (1 - dy) * (1 - dz) * dw * table9 +
           dx * (1 - dy) * dz * (1 - dw) * table10 +
           dx * (1 - dy) * dz * dw * table11 +
           dx * dy * (1 - dz) * (1 - dw) * table12 +
           dx * dy * (1 - dz) * dw * table13 +
           dx * dy * dz * (1 - dw) * table14 +
           dx * dy * dz * dw * table15;

  return result;
}

/*
 * @brief Interpolates 2d EAGLE table over one of the dimensions,
 * producing 1d table. May be useful if need to use bisection
 * method often with lots of iterations.
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i Bin index of value interpolated over
 * @param d_x Offset of value interpolated over
 * @param n_x Size of dimension interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void construct_1d_table_from_2d(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int array_size, double *result_table, float *ub, float *lb) {

  int index_low = 0, index_high = array_size;
  float d_high, d_low;
  
  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  int index[2];
  for (int i = index_low; i < index_high; i++) {
    index[0] = row_major_index_2d(x_i, i, n_x,
                                  array_size);
    index[1] = row_major_index_2d(x_i + 1, i, n_x,
                                  array_size);

    result_table[i] = (1 - d_x) * table[index[0]] +
                      d_x * table[index[1]];
  }
}

/*
 * @brief Interpolates 3d EAGLE table over two of the dimensions,
 * producing 1d table
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i, y_i Bin indices of values interpolated over
 * @param d_x, d_y Offsets of values interpolated over
 * @param n_x, n_y Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void construct_1d_table_from_3d(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int y_i, float d_y, int n_y, int array_size, double *result_table, float *ub, float *lb) {

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  int index[4];
  for (int i = index_low; i < index_high; i++) {
    index[0] = row_major_index_3d(x_i, y_i, i, n_x,
                                  n_y, array_size);
    index[1] = row_major_index_3d(x_i, y_i + 1, i, n_x,
                                  n_y, array_size);
    index[2] = row_major_index_3d(x_i + 1, y_i, i, n_x,
                                  n_y, array_size);
    index[3] = row_major_index_3d(x_i + 1, y_i + 1, i, n_x,
                                  n_y, array_size);

    result_table[i] = (1 - d_x) * (1 - d_y) * table[index[0]] +
                      (1 - d_x) * d_y * table[index[1]] +
                      d_x * (1 - d_y) * table[index[2]] +
                      d_x * d_y * table[index[3]];
  }
}

/*
 * @brief Interpolates 3d EAGLE table over two of the dimensions,
 * where one of the dimensions consists of entries for each of the
 * metals, producing 1d table.
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i Bin indices of values interpolated over
 * @param d_x Offsets of values interpolated over
 * @param n_x Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param abundance_ratio Array of ratio of particle metal abundances
 * to solar
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void
construct_1d_print_table_from_3d_elements(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int array_size, double *result_table, float *abundance_ratio, float *ub, float *lb) {
  int index[2];

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  for (int j = 0; j < cooling->N_Elements; j++) {
    for (int i = index_low; i < index_high; i++) {
      if (j == 0) result_table[i] = 0.0;
      index[0] = row_major_index_3d(x_i, i, j,
                                    n_x, array_size, cooling->N_Elements);
      index[1] = row_major_index_3d(x_i + 1, i, j,
                                    n_x, array_size, cooling->N_Elements);

      result_table[i + array_size*j] += ((1 - d_x) * table[index[0]] +
                         d_x * table[index[1]])*abundance_ratio[j+2];
    }
  }
}

/*
 * @brief Interpolates 3d EAGLE table over two of the dimensions,
 * where one of the dimensions consists of entries for each of the
 * metals, producing 1d table.
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i Bin indices of values interpolated over
 * @param d_x Offsets of values interpolated over
 * @param n_x Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param abundance_ratio Array of ratio of particle metal abundances
 * to solar
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void
construct_1d_table_from_3d_elements(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int array_size, double *result_table, float *abundance_ratio, float *ub, float *lb) {
  int index[2];

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  for (int j = 0; j < cooling->N_Elements; j++) {
    for (int i = index_low; i < index_high; i++) {
      if (j == 0) result_table[i] = 0.0;
      index[0] = row_major_index_3d(x_i, i, j,
                                    n_x, array_size, cooling->N_Elements);
      index[1] = row_major_index_3d(x_i + 1, i, j,
                                    n_x, array_size, cooling->N_Elements);

      result_table[i] += ((1 - d_x) * table[index[0]] +
                         d_x * table[index[1]])*abundance_ratio[j+2];
    }
  }
}


/*
 * @brief Interpolates 4d EAGLE table over three of the dimensions,
 * producing 1d table
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i, y_i, w_i Bin indices of values interpolated over
 * @param d_x, d_y, d_w Offsets of values interpolated over
 * @param n_x, n_y, n_w Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void construct_1d_table_from_4d(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x, int n_x,
    int y_i, float d_y, int n_y, int w_i, float d_w, int n_w, int array_size, double *result_table, float *ub, float *lb) {

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  int index[8];
  for (int i = index_low; i < index_high; i++) {
    index[0] =
        row_major_index_4d(x_i, y_i, w_i, i, n_x,
                           n_y, n_w, array_size);
    index[1] =
        row_major_index_4d(x_i, y_i, w_i + 1, i, n_x,
                           n_y, n_w, array_size);
    index[2] =
        row_major_index_4d(x_i, y_i + 1, w_i, i, n_x,
                           n_y, n_w, array_size);
    index[3] =
        row_major_index_4d(x_i, y_i + 1, w_i + 1, i, n_x,
                           n_y, n_w, array_size);
    index[4] =
        row_major_index_4d(x_i + 1, y_i, w_i, i, n_x,
                           n_y, n_w, array_size);
    index[5] =
        row_major_index_4d(x_i + 1, y_i, w_i + 1, i, n_x,
                           n_y, n_w, array_size);
    index[6] =
        row_major_index_4d(x_i + 1, y_i + 1, w_i, i, n_x,
                           n_y, n_w, array_size);
    index[7] = row_major_index_4d(x_i + 1, y_i + 1, w_i + 1, i,
                                  n_x, n_y,
                                  n_w, array_size);

    result_table[i] = (1 - d_x) * (1 - d_y) * (1 - d_w) * table[index[0]] +
                      (1 - d_x) * (1 - d_y) * d_w * table[index[1]] +
                      (1 - d_x) * d_y * (1 - d_w) * table[index[2]] +
                      (1 - d_x) * d_y * d_w * table[index[3]] +
                      d_x * (1 - d_y) * (1 - d_w) * table[index[4]] +
                      d_x * (1 - d_y) * d_w * table[index[5]] +
                      d_x * d_y * (1 - d_w) * table[index[6]] +
                      d_x * d_y * d_w * table[index[7]];
  }
}

/*
 * @brief Interpolates 4d EAGLE table over two of the dimensions,
 * producing 1d table. Used for obtaining contributions to cooling
 * rate from each of the metals
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i, y_i Bin indices of values interpolated over
 * @param d_x, d_y Offsets of values interpolated over
 * @param n_x, n_y Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param abundance_ratio Array of ratio of particle metal abundances
 * to solar
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void
construct_1d_print_table_from_4d_elements(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int y_i, float d_y, int n_y, int array_size, double *result_table, float *abundance_ratio, float *ub, float *lb) {

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  int index[4];

  for (int j = 0; j < cooling->N_Elements; j++) {
    for (int i = index_low; i < index_high; i++) {
      if (j == 0) result_table[i] = 0.0;
      index[0] = row_major_index_4d(x_i, y_i, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[1] = row_major_index_4d(x_i, y_i+1, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[2] = row_major_index_4d(x_i+1, y_i, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[3] = row_major_index_4d(x_i+1, y_i+1, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);

      result_table[i+array_size*j] = ((1 - d_x) * (1 - d_y) * table[index[0]] +
                         (1 - d_x) * d_y * table[index[1]] +
                         d_x * (1 - d_y) * table[index[2]] +
                         d_x * d_y * table[index[3]]) * abundance_ratio[j+2];
    }
  }
}

/*
 * @brief Interpolates 4d EAGLE table over three of the dimensions,
 * where one of the dimensions consists of entries for each of the
 * metals, producing 1d table.
 *
 * @param p Particle structure
 * @param cooling Cooling data structure
 * @param cosmo Cosmology structure
 * @param internal_const Physical constants structure
 * @param table Pointer to table we are interpolating
 * @param x_i, y_i Bin indices of values interpolated over
 * @param d_x, d_y Offsets of values interpolated over
 * @param n_x, n_y Sizes of dimensions interpolated over
 * @param array_size Size of output array
 * @param result_table Pointer to 1d interpolated table
 * @param abundance_ratio Array of ratio of particle metal abundances
 * to solar
 * @param ub Upper bound in temperature up to which to interpolate
 * @param lb Lower bound in temperature above which to interpolate
 */
__attribute__((always_inline)) INLINE void
construct_1d_table_from_4d_elements(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo,
    const struct phys_const *internal_const, float *table, int x_i, float d_x,
    int n_x, int y_i, float d_y, int n_y, int array_size, double *result_table, float *abundance_ratio, float *ub, float *lb) {
  int index[4];

  int index_low, index_high;
  float d_high, d_low;

  // only compute the part of the table between the bounds
  if (array_size == cooling->N_Temp) {
    get_index_1d(cooling->Temp, cooling->N_Temp,(*ub)/eagle_log_10,&index_high,&d_high);
    get_index_1d(cooling->Temp, cooling->N_Temp,(*lb)/eagle_log_10,&index_low,&d_low);
    if (index_high < array_size-1) {
      index_high += 2;
    } else {
      index_high = array_size;
    }
    if (index_low > 0) index_low--;
  }

  for (int j = 0; j < cooling->N_Elements; j++) {
    for (int i = index_low; i < index_high; i++) {
      if (j == 0) result_table[i] = 0.0;
      index[0] = row_major_index_4d(x_i, y_i, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[1] = row_major_index_4d(x_i, y_i+1, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[2] = row_major_index_4d(x_i+1, y_i, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);
      index[3] = row_major_index_4d(x_i+1, y_i+1, i, j, n_x,
                                    n_y, array_size, cooling->N_Elements);

      result_table[i] += ((1 - d_x) * (1 - d_y) * table[index[0]] +
                         (1 - d_x) * d_y * table[index[1]] +
                         d_x * (1 - d_y) * table[index[2]] +
                         d_x * d_y * table[index[3]]) * abundance_ratio[j+2];
    }
  }
}

/*
 * @brief Interpolates internal energy from temperature based 
 * on preinterpolated 1d table
 *
 * @param temp Temperature
 * @param temperature_table Table of internal energy values
 * @param cooling Cooling data structure
 */
__attribute__((always_inline)) INLINE double
eagle_convert_temp_to_u_1d_table(double temp, float *temperature_table,
                                 const struct cooling_function_data *restrict
                                     cooling) {

  int temp_i;
  float d_temp, u;

  get_index_1d(temperature_table, cooling->N_Temp, log10(temp), &temp_i,
               &d_temp);

  u = pow(10.0, interpol_1d(cooling->Therm, temp_i, d_temp));

  return u;
}

/*
 * @brief Interpolates temperature from internal energy based on table and
 * calculates the size of the internal energy cell for the specified
 * internal energy.
 *
 * @param log_10_u Log base 10 of internal energy
 * @param delta_u Pointer to size of internal energy cell
 * @param z_i Redshift index
 * @param n_h_i Hydrogen number density index
 * @param He_i Helium fraction index
 * @param d_z Redshift offset
 * @param d_n_h Hydrogen number density offset
 * @param d_He Helium fraction offset
 * @param cooling Cooling data structure
 * @param cosmo Cosmology data structure
 */
__attribute__((always_inline)) INLINE double
eagle_convert_u_to_temp(
    double log_10_u, float *delta_u,
    int z_i, int n_h_i, int He_i, 
    float d_z, float d_n_h, float d_He,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo) {

  int u_i;
  float d_u, logT;
  double upper, lower;

  get_index_1d(cooling->Therm, cooling->N_Temp, log_10_u, &u_i, &d_u);

  if (cosmo->z > cooling->reionisation_redshift){
  logT = interpol_3d(cooling->table.photodissociation_cooling.temperature,n_h_i,He_i,u_i,d_n_h,d_He,d_u,
    cooling->N_nH,cooling->N_He,cooling->N_Temp,&upper,&lower);
  } else if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
  logT = interpol_3d(cooling->table.no_compton_cooling.temperature,n_h_i,He_i,u_i,d_n_h,d_He,d_u,
    cooling->N_nH,cooling->N_He,cooling->N_Temp,&upper,&lower);
  } else {
  logT = interpol_4d(cooling->table.element_cooling.temperature, z_i,n_h_i,He_i,u_i,d_z,d_n_h,d_He,d_u,
    cooling->N_Redshifts,cooling->N_nH,cooling->N_He,cooling->N_Temp,&upper,&lower);
  }

  *delta_u =
      pow(10.0, cooling->Therm[u_i + 1]) - pow(10.0, cooling->Therm[u_i]);

  return logT;
}


/*
 * @brief Interpolates temperature from internal energy based on 
 * 1d preinterpolated table for the given particle. Also
 * calculates the size of the internal energy cell for the specified
 * internal energy
 *
 * @param log_10_u Log base 10 of internal energy
 * @param delta_u Pointer to size of internal energy cell
 * @param temperature_table 1d array of temperatures 
 * @param cooling Cooling data structure
 */
__attribute__((always_inline)) INLINE double
eagle_convert_u_to_temp_1d_table(
    double log_10_u, float *delta_u, double *temperature_table,
    const struct cooling_function_data *restrict cooling) {

  int u_i;
  float d_u, logT;

  get_index_1d(cooling->Therm, cooling->N_Temp, log_10_u, &u_i, &d_u);

  logT = interpol_1d_dbl(temperature_table, u_i, d_u);

  *delta_u =
      pow(10.0, cooling->Therm[u_i + 1]) - pow(10.0, cooling->Therm[u_i]);

  return logT;
}

/**
 * @brief construct all 1d tables
 *
 * @param z_index Redshift index
 * @param d_z Redshift offset
 * @param n_h_i Hydrogen number density index
 * @param d_n_h Hydrogen number density offset
 * @param He_i Helium fraction index
 * @param d_He Helium fraction offset
 * @param phys_const Physical constants data structure
 * @param cosmo Cosmology data structure
 * @param cooling Cooling data structure
 * @param p Particle data structure
 * @param abundance_ratio Array of ratio of particle metal abundances
 * to solar
 * @param H_plus_He_heat_table Pointer to table of values
 * which will contain contribution to cooling from 
 * hydrogen and helium
 * @param H_plus_He_electron_abundance_table Pointer to table of values 
 * which will contain contribution hydrogen and helium 
 * electron abundances
 * @param element_cooling_table Pointer to table of values which
 * whill contain contribution to heating from all of
 * the metals
 * @param element_cooling_print_table Pointer to table of values
 * which will contain contribution to heating from each of 
 * metals (has size n_elements*n_temp)
 * @param element_electron_abundance_table Pointer to table
 * of values which will contain metal electron abundances
 * @param temp_table Pointer to table of values which will 
 * contain log_10(temperature) for conversion of internal energy to
 * temperature
 * @param ub Upper bound in log(internal energy) up to which to produce tables
 * @param lb Lower bound in log(internal energy) above which to produce tables
 */
__attribute__((always_inline)) INLINE void construct_1d_tables(
		int z_index, float dz, int n_h_i, float d_n_h,
		int He_i, float d_He,
                const struct phys_const *restrict phys_const,
                const struct cosmology *restrict cosmo,
                const struct cooling_function_data *restrict cooling,
                const struct part *restrict p,
		float *abundance_ratio,
		double *H_plus_He_heat_table,
		double *H_plus_He_electron_abundance_table,
		double *element_cooling_table,
		double *element_cooling_print_table,
		double *element_electron_abundance_table,
		double *temp_table,
		float *ub, float *lb) {

  if (cosmo->z > cooling->reionisation_redshift) {
    // Photodissociation table
    construct_1d_table_from_3d(p, cooling, cosmo, phys_const,
    		cooling->table.photodissociation_cooling.temperature,
    		n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb);
    construct_1d_table_from_3d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.photodissociation_cooling.H_plus_He_heating, n_h_i, d_n_h, cooling->N_nH,  
    		He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb);
    construct_1d_table_from_3d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.photodissociation_cooling.H_plus_He_electron_abundance,
    		n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He,
    		cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
    construct_1d_table_from_3d_elements(
    		p, cooling, cosmo, phys_const,
    		cooling->table.photodissociation_cooling.metal_heating, 
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_table_from_2d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.photodissociation_cooling.electron_abundance,
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  } else if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    // High redshift table
    construct_1d_table_from_3d(p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.temperature,
    		n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb);
    construct_1d_table_from_3d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.H_plus_He_heating, n_h_i, d_n_h, cooling->N_nH,  
    		He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb);
    construct_1d_table_from_3d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.H_plus_He_electron_abundance,
    		n_h_i, d_n_h, cooling->N_nH, He_i, d_He, cooling->N_He,
    		cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
    construct_1d_table_from_3d_elements(
    		p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.metal_heating, 
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_print_table_from_3d_elements(
    		p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.metal_heating, 
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_print_table, abundance_ratio, ub, lb);
    construct_1d_table_from_2d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.no_compton_cooling.electron_abundance,
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  } else {
    // Normal tables 
    construct_1d_table_from_4d(p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.temperature,
    		z_index, dz, cooling->N_Redshifts, n_h_i, d_n_h,
    		cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, temp_table, ub, lb);
    construct_1d_table_from_4d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.H_plus_He_heating, z_index, dz, cooling->N_Redshifts, n_h_i, d_n_h,
    		cooling->N_nH, He_i, d_He, cooling->N_He, cooling->N_Temp, H_plus_He_heat_table, ub, lb);
    construct_1d_table_from_4d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.H_plus_He_electron_abundance, z_index,
    		dz, cooling->N_Redshifts, n_h_i, d_n_h, cooling->N_nH, He_i, d_He,
    		cooling->N_He, cooling->N_Temp, H_plus_He_electron_abundance_table, ub, lb);
    construct_1d_table_from_4d_elements(
    		p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.metal_heating, z_index, dz, cooling->N_Redshifts, 
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_table, abundance_ratio, ub, lb);
    construct_1d_print_table_from_4d_elements(
    		p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.metal_heating, z_index, dz, cooling->N_Redshifts, 
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_cooling_print_table, abundance_ratio, ub, lb);
    construct_1d_table_from_3d(
    		p, cooling, cosmo, phys_const,
    		cooling->table.element_cooling.electron_abundance, z_index, dz, cooling->N_Redshifts,
    		n_h_i, d_n_h, cooling->N_nH, cooling->N_Temp, element_electron_abundance_table, ub, lb);
  }		
}


