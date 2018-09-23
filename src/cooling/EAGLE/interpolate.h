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
#ifndef SWIFT_INTERPOL_EAGLE_H
#define SWIFT_INTERPOL_EAGLE_H

/**
 * @file src/cooling/EAGLE/interpolate.h
 * @brief Interpolation functions for EAGLE tables
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
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

// Used for identifying first call to get_redshift_index and previous redshift
// index
// static int z_index_initialised = 0;
// static int previous_z_index = -1;

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
 * be evenly spaced, and in increasing order.
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
  // Used for deciding whether to use table endpoints when finding index of a
  // variable
  const float table_bound_tolerance = 1.e-4;

  float dxm1 = (float)(ntable - 1) / (table[ntable - 1] - table[0]);

  if ((float)x <= table[0] + table_bound_tolerance) {
    // x less than first table element, use first element
    *i = 0;
    *dx = 0;
  } else if ((float)x >= table[ntable - 1] - table_bound_tolerance) {
    // x greater than last table element, use last element
    *i = ntable - 2;
    *dx = 1;
  } else {
    // x inside table range, interpolate
    *i = (int)floor(((float)x - table[0]) * dxm1);
#ifdef SWIFT_DEBUG_CHECKS
    if (*i > ntable || *i < 0)
      error(
          "trying to get index for value outside table range. Table size: %d, "
          "calculated index: %d, value: %.5e, table[0]: %.5e, grid size: %.5e",
          ntable, *i, x, table[0], dxm1);
#endif
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
    struct cooling_function_data *restrict cooling) {
  int i, iz;

  if (cooling->z_index_initialised == 0) {
    cooling->z_index_initialised = 1;
    cooling->previous_z_index = cooling->N_Redshifts - 2;

    /* this routine assumes cooling_redshifts table is in increasing order. Test
     * this. */
    for (i = 0; i < cooling->N_Redshifts - 2; i++)
      if (cooling->Redshifts[i + 1] < cooling->Redshifts[i]) {
        error("table should be in increasing order\n");
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
    for (iz = cooling->previous_z_index; iz >= 0; iz--) {
      if (z > cooling->Redshifts[iz]) {
        *dz = (z - cooling->Redshifts[iz]) /
              (cooling->Redshifts[iz + 1] - cooling->Redshifts[iz]);

        cooling->previous_z_index = *z_index = iz;

        break;
      }
    }
  }
}

/*
 * @brief Performs 2d interpolation of table based on indices
 * and offsets. Also computes values of table directly above and
 * below in second dimension by setting dy = 0,1 (used for
 * computing derivatives for Newton's method)
 *
 * @param table Pointer to flattened 2d table of values
 * @param i,j Indices of cell we are interpolating
 * @param dx,dy Offset within cell
 * @param nx,ny Table dimensions
 * @param upper Pointer to value set to the table value
 * when dy = 1 (used for calculating derivatives)
 * @param lower Pointer to value set to the table value
 * when dy = 0 (used for calculating derivatives)
 */
__attribute__((always_inline)) INLINE float interpolate_2d(
    const float *table, int i, int j, float dx, float dy, int nx, int ny,
    double *upper, double *lower) {
  float result;

  const float table0 = table[row_major_index_2d(i, j, nx, ny)];
  const float table1 = table[row_major_index_2d(i, j + 1, nx, ny)];
  const float table2 = table[row_major_index_2d(i + 1, j, nx, ny)];
  const float table3 = table[row_major_index_2d(i + 1, j + 1, nx, ny)];

  result = (1 - dx) * (1 - dy) * table0 + (1 - dx) * dy * table1 +
           dx * (1 - dy) * table2 + dx * dy * table3;
  dy = 1.0;
  *upper = (1 - dx) * (1 - dy) * table0 + (1 - dx) * dy * table1 +
           dx * (1 - dy) * table2 + dx * dy * table3;
  dy = 0.0;
  *lower = (1 - dx) * (1 - dy) * table0 + (1 - dx) * dy * table1 +
           dx * (1 - dy) * table2 + dx * dy * table3;

  return result;
}

/*
 * @brief Performs 3d interpolation of table based on indices
 * and offsets. Also computes values of table directly above and
 * below in third dimension by setting dz = 0,1 (used for
 * computing derivatives for Newton's method)
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
__attribute__((always_inline)) INLINE float interpolate_3d(
    const float *table, int i, int j, int k, float dx, float dy, float dz,
    int nx, int ny, int nz, double *upper, double *lower) {
  float result;

  const float table0 = table[row_major_index_3d(i, j, k, nx, ny, nz)];
  const float table1 = table[row_major_index_3d(i, j, k + 1, nx, ny, nz)];
  const float table2 = table[row_major_index_3d(i, j + 1, k, nx, ny, nz)];
  const float table3 = table[row_major_index_3d(i, j + 1, k + 1, nx, ny, nz)];
  const float table4 = table[row_major_index_3d(i + 1, j, k, nx, ny, nz)];
  const float table5 = table[row_major_index_3d(i + 1, j, k + 1, nx, ny, nz)];
  const float table6 = table[row_major_index_3d(i + 1, j + 1, k, nx, ny, nz)];
  const float table7 =
      table[row_major_index_3d(i + 1, j + 1, k + 1, nx, ny, nz)];

  result = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 + (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 + dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 + dx * dy * dz * table7;
  dz = 1.0;
  *upper = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 + (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 + dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 + dx * dy * dz * table7;
  dz = 0.0;
  *lower = (1 - dx) * (1 - dy) * (1 - dz) * table0 +
           (1 - dx) * (1 - dy) * dz * table1 +
           (1 - dx) * dy * (1 - dz) * table2 + (1 - dx) * dy * dz * table3 +
           dx * (1 - dy) * (1 - dz) * table4 + dx * (1 - dy) * dz * table5 +
           dx * dy * (1 - dz) * table6 + dx * dy * dz * table7;

  return result;
}

/*
 * @brief Performs 4d interpolation of table based on indices
 * and offsets. Also computes values of table directly above and
 * below in third or fourth dimension by setting dz,dw = 0,1 (which
 * dimension this happens for depends on whether this function is
 * used to interpolate one of the tables depending on metal species,
 * and is identified by cooling_is_metal flag. These values are used
 * for computing derivatives for Newton's method)
 *
 * @param table Pointer to flattened 4d table of values
 * @param i,j,k,l Indices of cell we are interpolating
 * @param dx,dy,dz,dw Offset within cell
 * @param nx,ny,nz,nw Table dimensions
 * @param upper Pointer to value set to the table value at
 * the when dw = 1 when used for interpolating table
 * depending on metal species, dz = 1 otherwise
 * (used for calculating derivatives)
 * @param lower Pointer to value set to the table value at
 * the when dw = 0 when used for interpolating table
 * depending on metal species, dz = 0 otherwise
 * (used for calculating derivatives)
 * @param cooling_is_metal If using to interpolate cooling
 * table for one of the metals set to 1, else 0.
 */
__attribute__((always_inline)) INLINE float interpolate_4d(
    const float *table, int i, int j, int k, int l, float dx, float dy,
    float dz, float dw, int nx, int ny, int nz, int nw, double *upper,
    double *lower, const int cooling_is_metal) {
  float result;

  const float table0 = table[row_major_index_4d(i, j, k, l, nx, ny, nz, nw)];
  const float table1 =
      table[row_major_index_4d(i, j, k, l + 1, nx, ny, nz, nw)];
  const float table2 =
      table[row_major_index_4d(i, j, k + 1, l, nx, ny, nz, nw)];
  const float table3 =
      table[row_major_index_4d(i, j, k + 1, l + 1, nx, ny, nz, nw)];
  const float table4 =
      table[row_major_index_4d(i, j + 1, k, l, nx, ny, nz, nw)];
  const float table5 =
      table[row_major_index_4d(i, j + 1, k, l + 1, nx, ny, nz, nw)];
  const float table6 =
      table[row_major_index_4d(i, j + 1, k + 1, l, nx, ny, nz, nw)];
  const float table7 =
      table[row_major_index_4d(i, j + 1, k + 1, l + 1, nx, ny, nz, nw)];
  const float table8 =
      table[row_major_index_4d(i + 1, j, k, l, nx, ny, nz, nw)];
  const float table9 =
      table[row_major_index_4d(i + 1, j, k, l + 1, nx, ny, nz, nw)];
  const float table10 =
      table[row_major_index_4d(i + 1, j, k + 1, l, nx, ny, nz, nw)];
  const float table11 =
      table[row_major_index_4d(i + 1, j, k + 1, l + 1, nx, ny, nz, nw)];
  const float table12 =
      table[row_major_index_4d(i + 1, j + 1, k, l, nx, ny, nz, nw)];
  const float table13 =
      table[row_major_index_4d(i + 1, j + 1, k, l + 1, nx, ny, nz, nw)];
  const float table14 =
      table[row_major_index_4d(i + 1, j + 1, k + 1, l, nx, ny, nz, nw)];
  const float table15 =
      table[row_major_index_4d(i + 1, j + 1, k + 1, l + 1, nx, ny, nz, nw)];

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
           dx * dy * dz * (1 - dw) * table14 + dx * dy * dz * dw * table15;
  if (cooling_is_metal == 1) {
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
           dx * dy * dz * (1 - dw) * table14 + dx * dy * dz * dw * table15;
  if (cooling_is_metal == 1) {
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
           dx * dy * dz * (1 - dw) * table14 + dx * dy * dz * dw * table15;

  return result;
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
__attribute__((always_inline)) INLINE double eagle_convert_u_to_temp(
    double log_10_u, float *delta_u, int n_h_i, int He_i, float d_n_h,
    float d_He, const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo) {

  int u_i;
  float d_u, logT;
  double upper, lower;

  get_index_1d(cooling->Therm, cooling->N_Temp, log_10_u, &u_i, &d_u);

  if (cosmo->z > cooling->reionisation_redshift) {
    logT = interpolate_3d(cooling->table.temperature, n_h_i, He_i, u_i, d_n_h,
                          d_He, d_u, cooling->N_nH, cooling->N_He,
                          cooling->N_Temp, &upper, &lower);
  } else if (cosmo->z > cooling->Redshifts[cooling->N_Redshifts - 1]) {
    logT = interpolate_3d(cooling->table.temperature, n_h_i, He_i, u_i, d_n_h,
                          d_He, d_u, cooling->N_nH, cooling->N_He,
                          cooling->N_Temp, &upper, &lower);
  } else {
    logT = interpolate_4d(cooling->table.temperature, 0, n_h_i, He_i, u_i,
                          cooling->dz, d_n_h, d_He, d_u, 2, cooling->N_nH,
                          cooling->N_He, cooling->N_Temp, &upper, &lower, 0);
  }

  *delta_u =
      exp(cooling->Therm[u_i + 1] * M_LN10) - exp(cooling->Therm[u_i] * M_LN10);

  return logT;
}

#endif
