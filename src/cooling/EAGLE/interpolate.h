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

/**
 * @brief Returns the 1d index of element with 2d indices x,y
 * from a flattened 2d array in row major order
 *
 * @param x, y Indices of element of interest
 * @param Nx, Ny Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_2d(int x, int y,
                                                             int Nx, int Ny) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
#endif
  return x * Ny + y;
}

/**
 * @brief Returns the 1d index of element with 3d indices x,y,z
 * from a flattened 3d array in row major order
 *
 * @param x, y, z Indices of element of interest
 * @param Nx, Ny, Nz Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_3d(int x, int y,
                                                             int z, int Nx,
                                                             int Ny, int Nz) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
#endif
  return x * Ny * Nz + y * Nz + z;
}

/**
 * @brief Returns the 1d index of element with 4d indices x,y,z,w
 * from a flattened 4d array in row major order
 *
 * @param x, y, z, w Indices of element of interest
 * @param Nx, Ny, Nz, Nw Sizes of array dimensions
 */
__attribute__((always_inline)) INLINE int row_major_index_4d(int x, int y,
                                                             int z, int w,
                                                             int Nx, int Ny,
                                                             int Nz, int Nw) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(x < Nx);
  assert(y < Ny);
  assert(z < Nz);
  assert(w < Nw);
#endif
  return x * Ny * Nz * Nw + y * Nz * Nw + z * Nw + w;
}

/**
 * @brief Finds the index of a value in a table and compute delta to nearest
 * element.
 *
 * This function assumes the table is monotonically increasing with a constant
 * difference between adjacent values.
 *
 * The returned difference is expressed in units of the table separation. This
 * means dx = (x - table[i]) / (table[i+1] - table[i]). It is always between
 * 0 and 1.
 *
 * @param table The table to search in.
 * @param size The number of elements in the table.
 * @param x The value to search for.
 * @param i (return) The index in the table of the element.
 * @param *dx (return) The difference between x and table[i]
 */
__attribute__((always_inline)) INLINE void get_index_1d(
    const float *restrict table, const int size, const float x, int *i,
    float *restrict dx) {

  const float delta = (size - 1) / (table[size - 1] - table[0]);

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(float, table, SWIFT_STRUCT_ALIGNMENT);

  if (x < table[0]) {
    /* We are below the first element */
    *i = 0;
    *dx = 0.f;
  } else if (x < table[size - 1]) {
    /* Normal case */
    *i = (x - table[0]) * delta;

#ifdef SWIFT_DEBUG_CHECKS
    if (*i > size || *i < 0) {
      error(
          "trying to get index for value outside table range. Table size: %d, "
          "calculated index: %d, value: %.5e, table[0]: %.5e, grid size: %.5e",
          size, *i, x, table[0], delta);
    }
#endif

    *dx = (x - table[*i]) * delta;
  } else {
    /* We are after the last element */
    *i = size - 2;
    *dx = 1.f;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (*dx < -0.001f || *dx > 1.001f) error("Invalid distance found dx=%e", *dx);
#endif
}

/**
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
 * @param cooling #cooling_function_data structure containing redshift table
 */
__attribute__((always_inline)) INLINE void get_redshift_index(
    float z, int *z_index, float *dz,
    struct cooling_function_data *restrict cooling) {
  int iz;

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
 */
__attribute__((always_inline)) INLINE float interpolate_2d(const float *table,
                                                           int i, int j,
                                                           float dx, float dy,
                                                           int nx, int ny) {
  const float table0 = table[row_major_index_2d(i, j, nx, ny)];
  const float table1 = table[row_major_index_2d(i, j + 1, nx, ny)];
  const float table2 = table[row_major_index_2d(i + 1, j, nx, ny)];
  const float table3 = table[row_major_index_2d(i + 1, j + 1, nx, ny)];

  return (1 - dx) * (1 - dy) * table0 + (1 - dx) * dy * table1 +
         dx * (1 - dy) * table2 + dx * dy * table3;
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
 */
__attribute__((always_inline)) INLINE float interpolate_3d(const float *table,
                                                           int i, int j, int k,
                                                           float dx, float dy,
                                                           float dz, int nx,
                                                           int ny, int nz) {
  const float table0 = table[row_major_index_3d(i, j, k, nx, ny, nz)];
  const float table1 = table[row_major_index_3d(i, j, k + 1, nx, ny, nz)];
  const float table2 = table[row_major_index_3d(i, j + 1, k, nx, ny, nz)];
  const float table3 = table[row_major_index_3d(i, j + 1, k + 1, nx, ny, nz)];
  const float table4 = table[row_major_index_3d(i + 1, j, k, nx, ny, nz)];
  const float table5 = table[row_major_index_3d(i + 1, j, k + 1, nx, ny, nz)];
  const float table6 = table[row_major_index_3d(i + 1, j + 1, k, nx, ny, nz)];
  const float table7 =
      table[row_major_index_3d(i + 1, j + 1, k + 1, nx, ny, nz)];

  return (1 - dx) * (1 - dy) * (1 - dz) * table0 +
         (1 - dx) * (1 - dy) * dz * table1 + (1 - dx) * dy * (1 - dz) * table2 +
         (1 - dx) * dy * dz * table3 + dx * (1 - dy) * (1 - dz) * table4 +
         dx * (1 - dy) * dz * table5 + dx * dy * (1 - dz) * table6 +
         dx * dy * dz * table7;
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
 */
__attribute__((always_inline)) INLINE float interpolate_4d(
    const float *table, int i, int j, int k, int l, float dx, float dy,
    float dz, float dw, int nx, int ny, int nz, int nw) {
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

  return (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table0 +
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
         dx * dy * (1 - dz) * dw * table13 + dx * dy * dz * (1 - dw) * table14 +
         dx * dy * dz * dw * table15;
}

#endif
