/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_FEEDBACK_INTERPOLATE_H
#define SWIFT_EAGLE_FEEDBACK_INTERPOLATE_H

#include "error.h"

/**
 * @brief Returns the 1d index of element with 2d indices i,j
 * from a flattened 2d array in row major order
 *
 * @param i, j Indices of element of interest
 * @param Nx, Ny Sizes of array dimensions
 */
__attribute__((always_inline)) static INLINE int row_major_index_2d(
    const int i, const int j, const int Nx, const int Ny) {
#ifdef SWIFT_DEBUG_CHECKS
  assert(i < Nx);
  assert(j < Ny);
#endif

  return i * Ny + j;
}

/**
 * @brief Returns the 1d index of element with 3d indices i,j,k
 * from a flattened 3d array in row major order
 *
 * @param i, j, k Indices of element of interest
 * @param Nx, Ny, Nz Sizes of array dimensions
 */
__attribute__((always_inline)) static INLINE int row_major_index_3d(
    const int i, const int j, const int k, const int Nx, const int Ny,
    const int Nz) {

#ifdef SWIFT_DEBUG_CHECKS
  assert(i < Nx);
  assert(j < Ny);
  assert(k < Nz);
#endif

  return i * Ny * Nz + j * Nz + k;
}

/**
 * @brief linear interpolation of 1d table at bin i with offset dx
 *
 * @param table array to interpolate
 * @param i index of cell to interpolate
 * @param dx offset within cell to interpolate
 */
__attribute__((always_inline)) static INLINE double interpolate_1d(
    const double* table, const int i, const float dx) {

  const float tx = 1.f - dx;

  return tx * table[i] + dx * table[i + 1];
}

/**
 * @brief linear interpolation of 2d table at bin i,j with offset dx, dy
 *
 * @param table array to interpolate
 * @param i row index of cell to interpolate
 * @param j column index of cell to interpolate
 * @param dx row offset within cell to interpolate
 * @param dy column offset within cell to interpolate
 */
__attribute__((always_inline)) static INLINE double interpolate_2d(
    double** table, const int i, const int j, const float dx, const float dy) {
  const float tx = 1.f - dx;
  const float ty = 1.f - dy;

  double result = tx * ty * table[i][j];
  result += tx * dy * table[i][j + 1];
  result += dx * ty * table[i + 1][j];
  result += dx * dy * table[i + 1][j + 1];

  return result;
}

/**
 * @brief linear interpolation of non-uniformly spaced 1d array, array_y, whose
 * positions are specified in array_x. The function takes an input value in the
 * range of array_x and returns a value interpolated from array_y with the same
 * offset in the corresponding bin.
 *
 * @param array_x array of values indicating positions of the array to be
 * interpolated
 * @param array_y array to interpolate
 * @param size length of array_x and array_y
 * @param x value within range of array_x indicating bin and offset within
 * array_y to interpolate
 */
static INLINE double interpolate_1D_non_uniform(const double* array_x,
                                                const double* array_y,
                                                const int size,
                                                const double x) {
#ifdef SWIFT_DEBUG_CHECKS

  /* Check that x within range of array_x */
  if (x < array_x[0])
    error("interpolating value less than array min. value %.5e array min %.5e",
          x, array_x[0]);
  if (x > array_x[size - 1])
    error(
        "interpolating value greater than array max. value %.5e array max %.5e",
        x, array_x[size - 1]);
#endif

  /* Find bin index and offset of x within array_x */
  int index = 0;
  while (array_x[index] <= x) index++;

  const double offset =
      (array_x[index] - x) / (array_x[index] - array_x[index - 1]);

  /* Interpolate array_y */
  return offset * array_y[index - 1] + (1. - offset) * array_y[index];
}

#endif /* SWIFT_EAGLE_FEEDBACK_INTERPOLATE_H */
